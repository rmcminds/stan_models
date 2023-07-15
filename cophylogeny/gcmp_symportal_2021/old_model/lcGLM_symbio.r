args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 0) {
    out_prefix <- args[[1]]
    outdir <- file.path('output', paste(out_prefix, gsub(':', '-', gsub(' ', '_', Sys.time())), sep = '_'))
} else {
    outdir <- file.path('output', gsub(':', '-', gsub(' ', '_', Sys.time())))
}

dir.create(outdir)
outfile <- file(file.path(outdir, 'runlog.log'), open = 'wt')
sink(outfile, type = 'output', split = T)
outfileerr <- file(file.path(outdir, 'runlogerr.log'), open = 'wt')
sink(outfileerr, type = 'message')

library(ape)
library(phytools)
library(phangorn)
library(reshape2)
library(parallel)
library(rstan)
library(RColorBrewer)
library(paleotree)
library(ggplot2)
library(geosphere) ##not sure this is being used? would be for an implementation of geographic covariance
library(geiger) ##?
library(shinystan) ##?
library(nlme) ##?
library(picante) ##?
library(Matrix)##?

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source(file.path('gcmp_stan_symportal', 'model', 'lcGLM_functions.r'))

microbeTreePath <- 'output/symbio_phylo/profile_WU_fastmeBal_correlated_chronos.tree' #symportal profile tree
hostTreePath <- 'output/host_phylo/combined_trees.newick' #set of Bayesian draws of host species phylogeny
mapfilePath <- 'raw_data/occurence/GCMP_symbiodinium_map_r29.txt' #mapping file
fulltablePath <- 'raw_data/occurence/species_table.txt' #symporal species table
modelPath <- 'gcmp_stan_symportal/model/logistic_cophylogenetic_GLM_varVar.stan' #stan model
seed <- 123
timeLimit <- 30 * 24 * 60 * 60

## filtration options
minCountSamp <- 5 # minimum sequencing depth for a sample to be included
minSamps <- 1 # minimum number of samples that a sequence variant is present in at the above threshold for it to be included
##

## model options
aveStDPriorExpect <- 1.0
aveStDMetaPriorExpect <- 1.0
NTrees <- 2 ## number of random trees to sample and to fit the model to
groupedFactors <- list(location           = c('ocean', 'ocean_area', 'reef_name'),
                       date               = 'concatenated_date',
                       colony             = 'colony_name',
                       tissue_compartment = 'tissue_compartment',
                       sequencing_depth   = 'log_sequencing_depth_scaled')
##

## Stan options
init_r <- 30
NCores <- 2 #NTrees
NChains <- 1 ## this is per tree; since I'm doing a large number of trees in parallel i'll just do one chain for each
NIterations <- 2^(12 - 1) ## will probably need >10,000? maybe start with 2, check convergence, double it, check, double, check, double, etc.?
max_treedepth <- 15 ## a warning will tell you if this needs to be increased
adapt_delta <- 0.8 ## increase this if you get 'divergences' - even one means your model fit sucks!
adapt_kappa <- 0.75
adapt_t0 <- 10
adapt_gamma <- 0.05
thin <- 2^(2 - 1) ## NIterations / thin number of Monte Carlo samples from the fit
##

## define the set of genera that we think our unidentified fungid samples could belong to
possibleFungidGenera <- c('Fungia_', 'Danafungia_', 'Cycloseris_', 'Pleuractis_')
## would be good here and for samples ID'd to genus to have some kind of weights for each possible species. i.e. it's really not possible that the acropora samples are A. cervicornis or palmata because their range is wrong, so even though I don't know which Acropora species exactly, I do have some info about the relative probabilities for many of them. similarly, some species might just be super rare in general and should likely be downweighted

sampleTipKey <- 'host_scientific_name'

filterfunction <- function(dfin) {
    levels(dfin[,sampleTipKey]) <- gsub(' ', '_', levels(dfin[,sampleTipKey]))
    undupsamps <- levels(dfin$physical_sample_name)[sapply(levels(dfin$physical_sample_name),
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
    contrasts(df2$ocean) <- 'contr.sum'
    contrasts(df2$ocean_area) <- 'contr.sum'
    contrasts(df2$host_scientific_name) <- 'contr.sum'
    if('W' %in% levels(df2$tissue_compartment)) {
        levels(df2$tissue_compartment)[levels(df2$tissue_compartment) == 'W'] <- NA
    }
    contrasts(df2$tissue_compartment) <- 'contr.sum'
    contrasts(df2$reef_name) <- 'contr.sum'
    contrasts(df2$concatenated_date) <- 'contr.sum'
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
##

####filter out species that i don't currently have a way to incorporate
## identify unique host species in the data, and replace spaces with underscores
initial.species <- gsub(' ', '_', levels(newmap[,sampleTipKey]))
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
y.old <- newtable[idx, colnames(newtable) %in% microbeTree$tip.label]
newermaptemp <- droplevels(newmap[idx,])
##

## define contrasts
newermap <- contrastfunction(newermaptemp)
##

## convert data to presence/absence
ybinary <- apply(y.old, 2, function(x) x > 0)
mode(ybinary) <- 'numeric'
##

## filter microbes that aren't in the minimum number of samples
y <- ybinary[, colSums(ybinary) >= minSamps]
finalMicrobeTree <- drop.tip(microbeTree, microbeTree$tip.label[!microbeTree$tip.label %in% colnames(y)])
y <- y[, finalMicrobeTree$tip.label]
##

## generate some summary numbers regarding microbes
microbeTips <- colnames(y)
NMicrobeTips <- length(microbeTips)
NIntMicrobeNodes <- finalMicrobeTree$Nnode
NMicrobeNodes <- NMicrobeTips + NIntMicrobeNodes - 1
microbeTreeDetails <- getTreeDetails(finalMicrobeTree)
##

## create ancestry matrix for microbes
microbeAncestors <- createAncestryMat(NMicrobeNodes,
                                      finalMicrobeTree,
                                      NMicrobeTips,
                                      microbeTips)
colnames(microbeAncestors)[1:NMicrobeTips] <- rownames(microbeAncestors)[1:NMicrobeTips] <- paste0('t', microbeTips)
##



## melt the data into long format to feed to model, and generate summary numbers about samples
senddat <- melt(y, varnames = c('sample', 'tip'), value.name = 'present')
sampleNames <- as.numeric(factor(senddat[,1]))
microbeTipNames <- as.numeric(factor(senddat[,2], levels = microbeTips))
present <- senddat[,3]
NObs <- nrow(senddat)
NSamples <- length(unique(sampleNames))
##

## make the mapping file match the set and order of samples in the data
newermap <- newermap[levels(factor(senddat[,1])),]
##

## prepare sequencing depth input
newermap$sequencing_depth <- rowSums(y.old[,microbeTips])
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

    #filter the tree only contain the sampled (or assigned) species
    hostTreesSampled[[i]] <- ladderize(drop.tip(hostTree[[hostTreeSampleNumbers[[i]]]],
                                                hostTree[[hostTreeSampleNumbers[[i]]]]$tip.label[!hostTree[[hostTreeSampleNumbers[[i]]]]$tip.label %in% levels(sampleMap[[i]][,sampleTipKey])]))
    
    #get some tree stats for later use
    hostTreeDetails[[i]] <- getTreeDetails(hostTreesSampled[[i]])
    
    hostTreesSampled[[i]]$edge.length <- hostTreesSampled[[i]]$edge.length / hostTreeDetails[[i]]$maxNHs
    
    hostTreeDetails[[i]] <- getTreeDetails(hostTreesSampled[[i]])
}
##

## create ancestry matrices for each host tree
NHostTips <- length(hostTreesSampled[[1]]$tip.label)
NIntHostNodes <- hostTreesSampled[[1]]$Nnode
NHostNodes <- NIntHostNodes + NHostTips - 1
hostAncestors <- list()
hostAncestorsExpanded <- list()
for (i in 1:NTrees) {
	hostAncestors[[i]] <- createAncestryMat(NHostNodes,
                                            hostTreesSampled[[i]],
                                            NHostTips,
                                            hostTreesSampled[[i]]$tip.label)
    hostAncestorsExpanded[[i]] <- hostAncestors[[i]][as.character(sampleMap[[i]][,sampleTipKey]),]
    rownames(hostAncestorsExpanded[[i]]) <- rownames(sampleMap[[i]])
}
##

## prepare data for the model matrix
modelform <- as.formula(paste0('~',paste(unlist(groupedFactors), collapse = ' + ')))
NFactors <- length(groupedFactors)
NSubPerFactor <- sapply(groupedFactors, length)
NSubfactors <- sum(NSubPerFactor)
modelMat <- model.matrix(modelform, model.frame(newermap, na.action = 'na.pass'))
modelMat[is.na(modelMat)] <- 0
sumconts <- names(attr(modelMat, "contrasts")[attr(modelMat, "contrasts") == 'contr.sum'])
##

## create matrix relating each 'effect' (categorical and numeric) to the 'factor' that it belongs to
stDAdjust <- rep(1, NSubfactors)
baseLevelMat <- NULL
subfactLevelMat <- matrix(NA, nrow = ncol(modelMat) - 1, ncol = NSubfactors)
colnames(subfactLevelMat) <- unlist(groupedFactors)
NSumTo0 <- 0
for(j in 1:NSubfactors) {
    newColumn <- as.numeric(attr(modelMat, 'assign')[-1] == j)
    if(colnames(subfactLevelMat)[[j]] %in% names(attr(modelMat, "contrasts"))) {
        if(attr(modelMat, "contrasts")[[colnames(subfactLevelMat)[[j]]]] == 'contr.sum') {
            ## if the contrast is a sum-to-zero ('effects') contrast, adjust the scale in preparation for making symmetrical marginal priors
            stDAdjust[[j]] <- 1 / sqrt(1 - 1 / (sum(newColumn) + 1))
            baseLevelMat <- rbind(baseLevelMat, -newColumn)
            subfactLevelMat[,j] <- newColumn * stDAdjust[[j]]
            NSumTo0 <- NSumTo0 + 1
        } else {
            subfactLevelMat[,j] <- newColumn
        }
    } else {
        subfactLevelMat[,j] <- newColumn
    }
}

##
NEffects <- ncol(modelMat) - 1
##

## rename factors that have 'sum contrasts' because by default they get arbitrary names (careful with interpretation of interactions... probably better to add them as separate 'main effects' produced by concatenation)
for(j in sumconts) {
    searchTerms <- paste0('^', j, 1:(nlevels(newermap[,j]) - 1), '$')
    replacementTerms <- paste0(j, levels(newermap[,j])[-nlevels(newermap[,j])])
    for(k in 1:length(searchTerms)) {
        colnames(modelMat) <- sub(searchTerms[[k]], replacementTerms[[k]], colnames(modelMat))
    }
}
rownames(subfactLevelMat) <- colnames(modelMat)[2:ncol(modelMat)]
##

## collect data to feed to stan
standat <- list()
for (i in 1:NTrees) {
    standat[[i]] <- list(NSamples                       = NSamples,
                         NObs                           = NObs,
                         NMicrobeNodes                  = NMicrobeNodes,
                         NMicrobeTips                   = NMicrobeTips,
                         NFactors                       = NFactors,
                         NSubPerFactor                  = NSubPerFactor,
                         NEffects                       = NEffects,
                         present                        = present,
                         sampleNames                    = sampleNames,
                         microbeTipNames                = microbeTipNames,
                         subfactLevelMat                = subfactLevelMat,
                         modelMat                       = cbind(modelMat, hostAncestorsExpanded[[i]]),
                         NSumTo0                        = NSumTo0,
                         baseLevelMat                   = baseLevelMat,
                         microbeAncestorsT              = t(microbeAncestors),
                         microbeTipAncestorsT           = t(cbind(1, microbeAncestors[1:NMicrobeTips, ])),
                         hostAncestors                  = hostAncestors[[i]],
                         hostTipAncestors               = hostAncestors[[i]][1:NHostTips, ],
                         hostEdges                      = hostTreeDetails[[i]]$edgeLengths,
                         microbeEdges                   = microbeTreeDetails$edgeLengths,
                         NHostNodes                     = NHostNodes,
                         NHostTips                      = NHostTips,
                         aveStDPriorExpect              = aveStDPriorExpect,
                         aveStDMetaPriorExpect          = aveStDMetaPriorExpect)
}

NMCSamples <- NIterations / thin
warmup <- NMCSamples / 2
##

sm <- stan_model(modelPath)

## run the model!
runStanModel()

## re-fit the model but shuffle the samples (to see if there are any biases generated by the sampling design)
runStanModel(shuffleSamples = T)

## re-fit the model but ignore all data (sampling from prior to see if there are any biases generated by the model itself)
runStanModel(noData = T)

## re-fit the model but shuffle the data at the observation level (to see if there are any biases generated by the sampling design)
runStanModel(shuffleData = T)

## fin
