args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 0) {
    output_prefix <- args[[1]]
    output_prefix <- file.path(Sys.getenv('HOME'), 'outputs', 'gcmp_symportal_2021', paste(args[[1]], gsub(':', '-', gsub(' ', '_', Sys.time())), sep = '_'))
} else {
    output_prefix <- file.path(Sys.getenv('HOME'), 'outputs', 'gcmp_symportal_2021', gsub(':', '-', gsub(' ', '_', Sys.time())))
}

dir.create(output_prefix)
outfile <- file(file.path(output_prefix, 'runlog.log'), open = 'wt')
sink(outfile, type = 'output', split = TRUE)
#outfileerr <- file(file.path(output_prefix, 'runlogerr.log'), open = 'wt')
#sink(outfileerr, type = 'message')

library(ape)
library(phytools)
library(phangorn)
library(parallel)
library(rstan)
library(RColorBrewer)
library(paleotree)
library(ggplot2)
library(cmdstanr)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

model_dir <- file.path(Sys.getenv('HOME'), 'scripts/stan_models/cophylogeny/')
source(file.path(model_dir, 'gcmp_symportal_2021', 'statistics', 'lcGLM_functions.r'))

input_prefix <- file.path(Sys.getenv('HOME'), 'data/gcmp_symportal_2021/')
preprocess_prefix <- file.path(Sys.getenv('HOME'), 'outputs/gcmp_symportal_2021/intermediate/')
include_path <- file.path(Sys.getenv('HOME'), 'scripts/stan_models/utility/')
model_name <- 'cophou_mix2'
engine <- c('sampling','advi')[2]
opencl <- FALSE

## Stan options
init_r <- 30
NCores <- 2 #NTrees
NChains <- 1 ## this is per tree; since I'm doing a large number of trees in parallel i'll just do one chain for each
NIterations <- 400 #2^(12 - 1) ## will probably need >10,000? maybe start with 2, check convergence, double it, check, double, check, double, etc.?
max_treedepth <- 15 ## a warning will tell you if this needs to be increased
adapt_delta <- 0.8 ## increase this if you get 'divergences' - even one means your model fit sucks!
adapt_kappa <- 0.75
adapt_t0 <- 10
adapt_gamma <- 0.05
thin <- 1 ## NIterations / thin number of Monte Carlo samples from the fit
##

sampling_commands <- list(sampling = paste(paste0('./', model_name),
                                           paste0('data file=', file.path(output_prefix, 'data.json')),
                                           #paste0('init=', file.path(output_prefix, 'inits.json')),
                                           'init=0.01',
                                           'output',
                                           paste0('file=', file.path(output_prefix, 'samples_sampling.csv')),
                                           paste0('refresh=', 1),
                                           'method=sample algorithm=hmc',
                                           'stepsize=1',
                                           'engine=nuts',
                                           'max_depth=10',
                                           'adapt t0=10',
                                           'delta=0.8',
                                           'kappa=0.75',
                                           'init_buffer=0',
                                           'window=2',
                                           'term_buffer=50',
                                           'num_warmup=200',
                                           'num_samples=200',
                                           ('opencl platform=0 device=0')[opencl],
                                           sep=' '),
                          advi = paste(paste0('./', model_name),
                                       paste0('data file=', file.path(output_prefix, 'data.json')),
                                       #paste0('init=', file.path(output_prefix, 'inits.json')),
                                       'init=0.01',
                                       'output',
                                       paste0('file=', file.path(output_prefix, 'samples_advi.csv')),
                                       paste0('refresh=', 100),
                                       'method=variational algorithm=meanfield',
                                       'grad_samples=1',
                                       'elbo_samples=1',
                                       'iter=3000',
                                       'eta=0.1',
                                       'adapt engaged=0',
                                       'tol_rel_obj=0.001',
                                       'eval_elbo=1',
                                       'output_samples=200',
                                       ('opencl platform=0 device=0')[opencl],
                                       sep=' '))

microbeTreePath <- file.path(preprocess_prefix, 'symbio_phylo/profile_WU_fastmeBal_correlated_chronos.tree') #symportal profile tree
hostTreePath <- file.path(preprocess_prefix, 'host_phylo/combined_trees.newick') #set of Bayesian draws of host species phylogeny
mapfilePath <- file.path(input_prefix, 'occurence/GCMP_symbiodinium_map_r29.txt') #mapping file
fulltablePath <- file.path(input_prefix, 'occurence/species_table.txt') #symporal species table
modelPath <- file.path(model_dir, paste0(model_name, '.stan')) #stan model

## filtration options
minCountSamp <- 5 # minimum sequencing depth for a sample to be included
minSamps <- 1 # minimum number of samples that a sequence variant is present in at the above threshold for it to be included
##

## model options
NTrees <- 2 ## number of random trees to sample and to fit the model to
sampleFactors <- list(location           = 'reef_name', #c('ocean', 'ocean_area', 'reef_name'),
                      date               = 'concatenated_date',
                      tissue_compartment = 'tissue_compartment',
                      sequencing_depth   = 'log_sequencing_depth_scaled')

hostFactors <- list(functional_group     = 'functional_group_sensu_darling',
                    zoox_in_propagules   = 'Symbiodinium.sp..in.propagules',
                    oz_disease_mean      = 'oz_disease_mean')
##

## define the set of genera that we think our unidentified fungid samples could belong to
possibleFungidGenera <- c('Fungia_', 'Danafungia_', 'Cycloseris_', 'Pleuractis_')
## would be good here and for samples ID'd to genus to have some kind of weights for each possible species. i.e. it's really not possible that the acropora samples are A. cervicornis or palmata because their range is wrong, so even though I don't know which Acropora species exactly, I do have some info about the relative probabilities for many of them. similarly, some species might just be super rare in general and should likely be downweighted

sampleTipKey <- 'host_scientific_name'

filterfunction <- function(dfin) {
    levels(dfin[,sampleTipKey]) <- gsub(' ', '_', levels(dfin[,sampleTipKey]))
    undupsamps <- unique(dfin$physical_sample_name)[sapply(unique(dfin$physical_sample_name),
                                                           function(x) sum(x == dfin$physical_sample_name) == 1)]
    df1 <- droplevels(dfin[dfin$physical_sample_name %in% undupsamps,])
    df2 <- droplevels(df1[(df1$tissue_compartment == 'W' |
                           df1$tissue_compartment == 'T' |
                           df1$tissue_compartment == 'S' |
                           df1$tissue_compartment == 'M') &
                           !grepl('Unknown|Missing',
                                  df1[,sampleTipKey],
                                  ignore.case = T),])
    return(df2)
}

contrastfunction <- function(dfin) {
    df2 <- dfin
    df2$functional_group_sensu_darling <- as.factor(df2$functional_group_sensu_darling)
    df2$Symbiodinium.sp..in.propagules <- as.factor(df2$Symbiodinium.sp..in.propagules)
    df2$oz_disease_mean <- as.numeric(df2$oz_disease_mean)
    df2$ocean <- as.factor(df2$ocean)
    contrasts(df2$ocean) <- 'contr.sum'
    df2$ocean_area <- as.factor(df2$ocean_area)
    contrasts(df2$ocean_area) <- 'contr.sum'
    df2$host_scientific_name <- as.factor(df2$host_scientific_name)
    contrasts(df2$host_scientific_name) <- 'contr.sum'
    df2$tissue_compartment <- as.factor(df2$tissue_compartment)
    if('W' %in% levels(df2$tissue_compartment)) {
        levels(df2$tissue_compartment)[levels(df2$tissue_compartment) == 'W'] <- NA
    }
    contrasts(df2$tissue_compartment) <- 'contr.sum'
    df2$reef_name <- as.factor(df2$reef_name)
    contrasts(df2$reef_name) <- 'contr.sum'
    df2$concatenated_date <- as.factor(df2$concatenated_date)
    contrasts(df2$concatenated_date) <- 'contr.sum'
    df2$colony_name <- as.factor(df2$colony_name)
    contrasts(df2$colony_name) <- 'contr.sum'
    levels(df2[,sampleTipKey])[levels(df2[,sampleTipKey]) == 'Homophyllia_hillae'] <- "Homophyllia_bowerbanki"
    levels(df2[,sampleTipKey])[levels(df2[,sampleTipKey]) == 'Pocillopora_eydouxi'] <- "Pocillopora_grandis"
    return(df2)
}

##import data
fulltable <- t(read.table(fulltablePath,
                          header = T,
                          sep = '\t',
                          comment.char = '',
                          row.names = 1,
                          check.names = F))

newtable <- fulltable[-nrow(fulltable),] # last line in symportal file is not part of the table
mode(newtable) <- 'numeric'
##

## import phylogenies
hostTree <- read.tree(hostTreePath)
hostTree <- .compressTipLabel(hostTree)
microbeTree <- read.tree(microbeTreePath)

## import mapping file
map <- read.table(mapfilePath,
                  header = T,
                  sep = '\t',
                  comment.char = '',
                  check.names = F)
rownames(map) <- map[,'#SampleID']
newmap <- filterfunction(map)
newmap[,sampleTipKey] <- as.factor(newmap[,sampleTipKey])
##

####filter out species that i don't currently have a way to incorporate
## identify unique host species in the data, and replace spaces with underscores
levels(newmap[,sampleTipKey]) <- gsub(' ', '_', levels(newmap[,sampleTipKey]))
initial.species <- levels(newmap[,sampleTipKey])
##

## identify the Scleractinian species in the data that do not exist in the template tree
initial.species.missing <- initial.species[!initial.species %in% grep(paste(c(attr(hostTree, "TipLabel"),
                                                                              'Fungid',
                                                                              'not_applicable'),
                                                                            collapse = '|'),
                                                                      initial.species,
                                                                      ignore.case = T,
                                                                      value = T)]
initialGeneraOfUnknowns <- sapply(initial.species.missing, function(x) strsplit(x, '_')[[1]][[1]])
##

treeGenera <- c(unique(sapply(attr(hostTree, "TipLabel"), function(x) strsplit(x, '_')[[1]][[1]])), 'Fungid')

for (i in initial.species.missing) {
    if (!initialGeneraOfUnknowns[i] %in% treeGenera) {
        newmap <- droplevels(newmap[newmap[,sampleTipKey] != i,])
        cat(paste0('filtering out ', i))
    }
}
####


## merge data and mapping file
idx <- rownames(newtable)[rownames(newtable) %in% rownames(newmap) & rowSums(newtable) >= minCountSamp]
abundance <- newtable[idx, colnames(newtable) %in% microbeTree$tip.label]
newermaptemp <- droplevels(newmap[idx,])
##

## define contrasts
newermap <- contrastfunction(newermaptemp)
##

## convert data to presence/absence
ybinary <- apply(abundance, 2, function(x) x > 0)
mode(ybinary) <- 'numeric'
##

## filter microbes that aren't in the minimum number of samples
present <- ybinary[, colSums(ybinary) >= minSamps]
microbeTol <- exp(mean(log(microbeTree$edge.length)))
multiTree <- force.ultrametric(di2multi(microbeTree, tol=microbeTol))
finalMicrobeTree <- reorder(drop.tip(multiTree, multiTree$tip.label[!multiTree$tip.label %in% colnames(present)]), 'pruningwise')
finalMicrobeTree$edge.length <- finalMicrobeTree$edge.length / exp(mean(log(finalMicrobeTree$edge.length)))
present <- present[, finalMicrobeTree$tip.label]
NS <- nrow(present)
##

## generate some summary numbers regarding microbes
microbeTips <- colnames(present)
NT_m <- length(microbeTips)
NI_m <- finalMicrobeTree$Nnode
NN_m <- NI_m + NT_m
sa_m <- rbind(finalMicrobeTree$edge,c(0,NT_m+1))[NN_m:1,]
time_m <- c(finalMicrobeTree$edge.length[1:length(finalMicrobeTree$tip.label)],1,finalMicrobeTree$edge.length[(length(finalMicrobeTree$tip.label)+1):length(finalMicrobeTree$edge.length)])
NDesc_m <- sapply((NT_m+1):NN_m, function(x) sum(x == sa_m[,1]))
microbeTreeDetails <- getTreeDetails(finalMicrobeTree)
##


## prepare sequencing depth input
newermap$sequencing_depth <- rowSums(abundance[,microbeTips])
newermap$log_sequencing_depth <- log(newermap$sequencing_depth)
newermap$log_sequencing_depth_scaled <- scale(newermap$log_sequencing_depth)
##

## extract all the possible species that 'fungid' could refer to
possibleFungidSpecs <- grep(paste0(paste(possibleFungidGenera, collapse = '|'), '_'), attr(hostTree, "TipLabel"), value = T)
##

## identify unique host species in the data, and replace spaces with underscores
study.species <- gsub(' ', '_', levels(newermap[,sampleTipKey]))
##

## identify the Scleractinian species in the data that do not exist in the template tree
study.species.missing <- study.species[!study.species %in% grep(paste(c(attr(hostTree, "TipLabel"),
                                                                        'Fungid',
                                                                        'not_applicable'),
                                                                      collapse = '|'),
                                                                study.species,
                                                                ignore.case = T,
                                                                value = T)]
generaOfUnknowns <- sapply(study.species.missing, function(x) strsplit(x, '_')[[1]][[1]])
##

## generate multiple random samples of the map and trees
sampleMap <- list()
hostTreeSampleNumbers <- sample(1:length(hostTree), NTrees)
hostTreesSampled <- list()
hostTreeDetails <- list()
time_h <- list()
sa_h <- list()
NDesc_h <- list()
idx_host <- list()
for(i in 1:NTrees) {
    sampleMap[[i]] <- newermap
    fungidSps <- grep('Fungid', levels(sampleMap[[i]][,sampleTipKey]))
    ##assign unidentified Fungids to a random member of the group independently for each tree
    levels(sampleMap[[i]][,sampleTipKey])[fungidSps] <- sample(possibleFungidSpecs[!possibleFungidSpecs %in% levels(sampleMap[[i]][,sampleTipKey])],
                                                               length(fungidSps))
    for (j in unique(generaOfUnknowns)) {
        possibleGenera <- attr(hostTree, "TipLabel")[!attr(hostTree, "TipLabel") %in% levels(sampleMap[[i]][,sampleTipKey])]
        if(!any(grepl(j, possibleGenera))) {
            possibleGenera <- attr(hostTree, "TipLabel")
        }
        ## assign unidentified species to a random member of their genus independently for each tree
        levels(sampleMap[[i]][,sampleTipKey])[levels(sampleMap[[i]][,sampleTipKey]) %in% grep(paste0(j, '_'),
                                                                                              study.species.missing,
                                                                                              value = T)] <- sample(grep(paste0(j, '_'),
                                                                                                                         possibleGenera,
                                                                                                                         value = T),
                                                                                                                    sum(generaOfUnknowns == j))
    }

    # replace species tips with colony tips could do stochastic binary trees if doing many samples, but for now just do average
    hostTreesSampled[[i]] <- hostTree[[hostTreeSampleNumbers[[i]]]]
    for(species in levels(sampleMap[[i]][,sampleTipKey])) {
        temp <- hostTreesSampled[[i]]
        idx <- temp$edge[,2] == which(temp$tip.label==species)
        newlen <- temp$edge.length[idx] / 2
        temp$edge.length[idx] <- newlen
        temp <- add.tips(temp, unique(as.character(sampleMap[[i]][sampleMap[[i]][,sampleTipKey] == species, 'colony_name'])), which(temp$tip.label==species), newlen)
        hostTreesSampled[[i]] <- drop.tip(temp, species)
    }

    #filter the tree to only contain the samples
    hostTreesSampled[[i]] <- reorder(drop.tip(hostTreesSampled[[i]],
                                      hostTreesSampled[[i]]$tip.label[!hostTreesSampled[[i]]$tip.label %in% levels(sampleMap[[i]][,'colony_name'])]), 'pruningwise')


    #get some tree stats for later use
    hostTreeDetails[[i]] <- getTreeDetails(hostTreesSampled[[i]])

    hostTreesSampled[[i]]$edge.length <- hostTreesSampled[[i]]$edge.length / exp(mean(log(hostTreesSampled[[i]]$edge.length)))

    hostTreeDetails[[i]] <- getTreeDetails(hostTreesSampled[[i]])

    sa_h[[i]] <- rbind(hostTreesSampled[[i]]$edge,c(0,length(hostTreesSampled[[i]]$tip.label)+1))[(nrow(hostTreesSampled[[i]]$edge)+1):1,]
    time_h[[i]] <- c(hostTreesSampled[[i]]$edge.length[1:length(hostTreesSampled[[i]]$tip.label)],1,hostTreesSampled[[i]]$edge.length[(length(hostTreesSampled[[i]]$tip.label)+1):length(hostTreesSampled[[i]]$edge.length)])
    NDesc_h[[i]] <- sapply((length(hostTreesSampled[[i]]$tip.label)+1):(length(hostTreesSampled[[i]]$edge.length)+1), function(x) sum(x == sa_h[[i]][,1]))
    idx_host[[i]] <- sapply(sampleMap[[i]][,'colony_name'], function(x) which(x == hostTreesSampled[[i]]$tip.label))
}
##

## create ancestry matrices for each host tree
NT_h <- length(hostTreesSampled[[1]]$tip.label)
NI_h <- hostTreesSampled[[1]]$Nnode
NN_h <- NT_h + NI_h
##

## prepare data for the model matrix
idx_s <- 0
modelMat <- matrix(NA,nrow=NS,ncol=0)
for(fact in unlist(sampleFactors)) {
    if(is.factor(newermap[,fact])) {
        n <- length(levels(newermap[,fact]))
    } else {
        n <- 1
    }
    idx_s <- c(idx_s, rep(max(idx_s) + 1, n))
    modelMat <- cbind(modelMat, model.matrix(as.formula(paste0('~0+',fact)), model.frame(newermap, na.action = 'na.pass')))
}
modelMat[is.na(modelMat)] <- 0
idx_s <- idx_s[-1]
##

##
map_tips <- data.frame(newermap[which(newermap$colony_name==hostTreesSampled[[1]]$tip.label[[1]])[[1]],])
for(tip in hostTreesSampled[[1]]$tip.label[-1]) {
    map_tips <- rbind(map_tips, data.frame(newermap[which(newermap$colony_name==tip)[[1]],]))
}
levels(map_tips$functional_group_sensu_darling)[levels(map_tips$functional_group_sensu_darling) %in% c('Missing: Not collected','Unknown')] <- NA
levels(map_tips$Symbiodinium.sp..in.propagules)[levels(map_tips$Symbiodinium.sp..in.propagules) %in% c('Missing: Not collected','Unknown')] <- NA
idx_theta <- 0
thetaModel <- matrix(NA,nrow=NT_h,ncol=0)
for(fact in unlist(hostFactors)) {
    if(is.factor(map_tips[,fact])) {
        n <- length(levels(map_tips[,fact]))
    } else {
        n <- 1
    }
    idx_theta <- c(idx_theta, rep(max(idx_theta) + 1, n))
    thetaModel <- cbind(thetaModel, model.matrix(as.formula(paste0('~0+',fact)), model.frame(map_tips, na.action = 'na.pass')))
}
##

##
NSB_s <- length(sampleFactors)
NB_s <- ncol(modelMat)
##


## collect data to feed to stan
standat <- list()
for (i in 1:NTrees) {
    standat[[i]] <- list(NS           = NS,
                         NI_m         = NI_m,
                         NT_m         = NT_m,
                         NB_s         = NB_s,
                         NSB_s        = NSB_s,
                         idx_s        = idx_s,
                         present      = present,
                         modelMat     = modelMat,
                         global_scale = 0.1 * sd(present),
                         time_h       = time_h[[i]],
                         time_m       = time_m,
                         NI_h         = NI_h,
                         NT_h         = NT_h,
                         self_h       = sa_h[[i]][,2],
                         ancestor_h   = sa_h[[i]][,1],
                         self_m       = sa_m[,2],
                         ancestor_m   = sa_m[,1],
                         idx_host     = idx_host[[i]],
                         NDesc_h      = NDesc_h[[i]],
                         NDesc_m      = NDesc_m,
                         X_s          = modelMat,
                         NB_m         = 0,
                         NB_h         = 0,
                         X_m_raw      = matrix(0, nrow=0, ncol=NT_m),
                         X_h_raw      = matrix(0, nrow=0, ncol=NT_h),
                         NM_X_m       = 0,
                         NM_X_h       = 0,
                         idx_X_m_m_x  = vector(),
                         idx_X_m_m_y  = vector(),
                         idx_X_h_m_x  = vector(),
                         idx_X_h_m_y  = vector())
}

save.image(file.path(output_prefix, 'setup.RData'))

write_stan_json(standat[[1]], file.path(output_prefix, 'data.json'))
#write_stan_json(init, file.path(output_prefix, 'inits.json'))

setwd(cmdstan_path())
system(paste0(c('make ', 'make STAN_OPENCL=true ')[opencl+1], 'STANCFLAGS="--include-paths=', include_path, '" ', file.path(model_dir, model_name)))

setwd(model_dir)
print(sampling_commands[[engine]])
print(date())
system(sampling_commands[[engine]])

stan.fit <- read_cmdstan_csv(file.path(output_prefix, paste0('samples_',engine,'.csv')),
                             format = 'draws_array')

save.image(file.path(output_prefix, paste0('results_',engine,'.RData')))

if(length(stan.fit) == 2) {
    stan.fit.draws <- stan.fit[[2]]
} else {
    stan.fit.draws <- stan.fit$post_warmup_draws
}

sumfunc <- median

colorpal <- colorRampPalette(brewer.pal(9, 'Blues'))
plotcolors <- c('white', colorpal(100), 'black')

mycols <- rev(brewer.pal(11, 'RdYlBu'))
mycols[[6]] <- 'white'
colorpal <- colorRampPalette(mycols)
plotcolorsVar <- c(colorpal(101))

plotM <- force.ultrametric(multi2di(ladderize(finalMicrobeTree)), 'extend')
plotH <- multi2di(ladderize(hostTreesSampled[[1]],right=FALSE))

inheritance <- stan.fit.draws[,,grep('^inheritance\\[.*',dimnames(stan.fit.draws)[[3]]), drop=FALSE]
inheritance <- apply(inheritance,3,sumfunc)
inheritance <- array(inheritance,dim=c(NI_h,NI_m))

inheritance_ia <- inheritance[sa_h[[1]][sapply(1:NT_h, \(x) sa_h[[1]][x == sa_h[[1]][,2],1]),1] - NT_h, sa_m[sapply(1:NT_m, \(x) sa_m[x == sa_m[,2],1]),1] - NT_m]
rownames(inheritance_ia) <- hostTreesSampled[[1]]$tip.label
colnames(inheritance_ia) <- finalMicrobeTree$tip.label

heatmap(t(logit(inheritance_ia)),
        Rowv   = as.dendrogram(as.hclust(plotM)),
        Colv   = as.dendrogram(as.hclust(plotH)),
        col    = plotcolorsVar,
        cexCol = 0.2,
        cexRow = 0.1,
        scale  = 'none')


alpha_h <- stan.fit.draws[,,grep('alpha_h\\[.*',dimnames(stan.fit.draws)[[3]]), drop=FALSE]
alpha_h <- apply(alpha_h,3,sumfunc)
alpha_h <- array(alpha_h,dim=c(NN_h,NN_m))

alpha_h_filt <- alpha_h[1:NT_h,1:NT_m]
rownames(alpha_h_filt) <- hostTreesSampled[[1]]$tip.label
colnames(alpha_h_filt) <- finalMicrobeTree$tip.label

heatmap(t(log(alpha_h_filt)),
        Rowv   = as.dendrogram(as.hclust(ladderize(force.ultrametric(multi2di(finalMicrobeTree))))),
        Colv   = as.dendrogram(as.hclust(multi2di(ladderize(hostTreesSampled[[1]],right=FALSE)))),
        col    = plotcolorsVar,
        cexCol = 0.2,
        cexRow = 0.1,
        scale  = 'none')

alpha_m <- stan.fit.draws[,,grep('alpha_m\\[.*',dimnames(stan.fit.draws)[[3]]), drop=FALSE]
alpha_m <- apply(alpha_m,3,sumfunc)
alpha_m <- array(alpha_m,dim=c(NN_h,NN_m))

alpha_m_filt <- alpha_m[1:NT_h,1:NT_m]
rownames(alpha_m_filt) <- hostTreesSampled[[1]]$tip.label
colnames(alpha_m_filt) <- finalMicrobeTree$tip.label

heatmap(t(log(alpha_m_filt)),
        Rowv   = as.dendrogram(as.hclust(ladderize(force.ultrametric(multi2di(finalMicrobeTree))))),
        Colv   = as.dendrogram(as.hclust(multi2di(ladderize(hostTreesSampled[[1]],right=FALSE)))),
        col    = plotcolorsVar,
        cexCol = 0.2,
        cexRow = 0.1,
        scale  = 'none')


present_anc_logit <- stan.fit.draws[,,grep('present_anc_logit\\[.*',dimnames(stan.fit.draws)[[3]]), drop=FALSE]
present_anc_logit <- apply(present_anc_logit,3,sumfunc)
present_anc_logit <- array(present_anc_logit,dim=c(NN_h,NN_m))

present_anc_logit_filt <- present_anc_logit[1:NT_h,1:NT_m]
rownames(present_anc_logit_filt) <- hostTreesSampled[[1]]$tip.label
colnames(present_anc_logit_filt) <- finalMicrobeTree$tip.label

heatmap(t(present_anc_logit_filt),
        Rowv   = as.dendrogram(as.hclust(ladderize(force.ultrametric(multi2di(finalMicrobeTree))))),
        Colv   = as.dendrogram(as.hclust(multi2di(ladderize(hostTreesSampled[[1]],right=FALSE)))),
        col    = plotcolorsVar,
        cexCol = 0.2,
        cexRow = 0.1,
        scale  = 'none')

sigma_host <- stan.fit.draws[,,grep('sigma_host\\[.*',dimnames(stan.fit.draws)[[3]]), drop=FALSE]
sigma_host <- apply(sigma_host,3,sumfunc)
sigma_host <- array(sigma_host,dim=c(NN_h,NN_m))

sigma_host_filt <- sigma_host[1:NT_h,1:NT_m]
rownames(sigma_host_filt) <- hostTreesSampled[[1]]$tip.label
colnames(sigma_host_filt) <- finalMicrobeTree$tip.label

heatmap(t(log(sigma_host_filt)),
        Rowv   = as.dendrogram(as.hclust(ladderize(force.ultrametric(multi2di(finalMicrobeTree))))),
        Colv   = as.dendrogram(as.hclust(multi2di(ladderize(hostTreesSampled[[1]],right=FALSE)))),
        col    = plotcolorsVar,
        cexCol = 0.2,
        cexRow = 0.1,
        scale  = 'none')

sigma_host[sa_h[[1]][1,2],sa_m[1,2]]


sigma_drift_h <- stan.fit.draws[,,grep('sigma_drift_h\\[.*',dimnames(stan.fit.draws)[[3]]), drop=FALSE]
sigma_drift_h <- apply(sigma_drift_h,3,sumfunc)
sigma_drift_h <- array(sigma_drift_h,dim=c(NN_h,NN_m))

sigma_drift_h_filt <- sigma_drift_h[1:NT_h,1:NT_m]
rownames(sigma_drift_h_filt) <- hostTreesSampled[[1]]$tip.label
colnames(sigma_drift_h_filt) <- finalMicrobeTree$tip.label

heatmap(t(log(sigma_drift_h_filt)),
        Rowv   = as.dendrogram(as.hclust(ladderize(force.ultrametric(multi2di(finalMicrobeTree))))),
        Colv   = as.dendrogram(as.hclust(multi2di(ladderize(hostTreesSampled[[1]],right=FALSE)))),
        col    = plotcolorsVar,
        cexCol = 0.2,
        cexRow = 0.1,
        scale  = 'none')

sigma_drift_m <- stan.fit.draws[,,grep('sigma_drift_m\\[.*',dimnames(stan.fit.draws)[[3]]), drop=FALSE]
sigma_drift_m <- apply(sigma_drift_m,3,sumfunc)
sigma_drift_m <- array(sigma_drift_m,dim=c(NN_h,NN_m))

sigma_drift_m_filt <- sigma_drift_m[1:NT_h,1:NT_m]
rownames(sigma_drift_m_filt) <- hostTreesSampled[[1]]$tip.label
colnames(sigma_drift_m_filt) <- finalMicrobeTree$tip.label

heatmap(t(log(sigma_drift_m_filt)),
        Rowv   = as.dendrogram(as.hclust(ladderize(force.ultrametric(multi2di(finalMicrobeTree))))),
        Colv   = as.dendrogram(as.hclust(multi2di(ladderize(hostTreesSampled[[1]],right=FALSE)))),
        col    = plotcolorsVar,
        cexCol = 0.2,
        cexRow = 0.1,
        scale  = 'none')

present_s <- t(sapply(hostTreesSampled[[1]]$tip.label, function(colony) apply(present[rownames(newermap)[newermap$colony_name == colony],, drop=FALSE], 2, mean)))
heatmap(t(present_s),
        Rowv   = as.dendrogram(as.hclust(ladderize(force.ultrametric(multi2di(finalMicrobeTree))))),
        Colv   = as.dendrogram(as.hclust(multi2di(ladderize(hostTreesSampled[[1]],right=FALSE)))),
        col    = plotcolors,
        cexCol = 0.2,
        cexRow = 0.1,
        scale  = 'none')


sigma_s <- stan.fit.draws[,,grep('sigma_s\\[.*',dimnames(stan.fit.draws)[[3]]), drop=FALSE]
sigma_s <- apply(sigma_s,3,sumfunc)
names(sigma_s) <- names(sampleFactors)

sigma_sigma_host <- sumfunc(stan.fit.draws[,,'sigma_sigma_host' == dimnames(stan.fit.draws)[[3]], drop=FALSE])
sigma_sigma_drift_h <- sumfunc(stan.fit.draws[,,'sigma_sigma_drift_h' == dimnames(stan.fit.draws)[[3]], drop=FALSE])
sigma_sigma_drift_m <- sumfunc(stan.fit.draws[,,'sigma_sigma_drift_m' == dimnames(stan.fit.draws)[[3]], drop=FALSE])
sigma_inheritance <- sumfunc(stan.fit.draws[,,'sigma_inheritance' == dimnames(stan.fit.draws)[[3]], drop=FALSE])
uninherited <- sumfunc(stan.fit.draws[,,'uninherited' == dimnames(stan.fit.draws)[[3]], drop=FALSE])

multiTree <- di2multi(finalMicrobeTree, tol=exp(mean(log(finalMicrobeTree$edge.length)) - sd(log(finalMicrobeTree$edge.length))))

#NMCSamples <- NIterations / thin
#warmup <- NMCSamples / 2
##

#sm <- stan_model(modelPath)

## run the model!
#runStanModel()

## re-fit the model but shuffle the samples (to see if there are any biases generated by the sampling design)
#runStanModel(shuffleSamples = T)

## re-fit the model but ignore all data (sampling from prior to see if there are any biases generated by the model itself)
#runStanModel(noData = T)

## re-fit the model but shuffle the data at the observation level (to see if there are any biases generated by the sampling design)
#runStanModel(shuffleData = T)

## fin
