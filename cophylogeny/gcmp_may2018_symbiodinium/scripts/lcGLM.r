library(ape)
library(phangorn)
library(geosphere)
library(rstan)
library(geiger)
library(phytools)
library(shinystan)
library(nlme)
library(reshape2)
library(paleotree)
library(parallel)
library(ggplot2)
library(RColorBrewer)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
source(file.path('scripts', 'lcGLM_preparation_functions.r'))
source(file.path('scripts', 'lcGLM_summary_functions.r'))


microbeTreePath <- 'raw_data/species_rooted_mapped_labels.tree' #ASTRAL microbial phylogeny
hostTreePath <- 'raw_data/combined_trees.newick' #set of Bayesian draws of host species phylogeny
mapfilePath <- 'raw_data/GCMP_symbiodinium_map2.txt' #mapping file
fulltablePath <- 'raw_data/species_table.txt' #250 bp deblur otu table output
modelPath <- 'scripts/logistic_cophylogenetic_GLM_varVar.stan' #stan model
seed <- 123
timeLimit <- 30 * 24 * 60 * 60

outdir <- file.path('output', gsub(':', '-', gsub(' ', '_', Sys.time())))

## filtration options
minCountSamp <- 5 # minimum sequencing depth for a sample to be included
minSamps <- 1 # minimum number of samples that a sequence variant is present in at the above threshold for it to be included
##

## model options
aveStDPriorExpect <- 1.0
aveStDMetaPriorExpect <- 0.1
hostOUAlphaPriorExpect <- 1.0
microbeOUAlphaPriorExpect <- 1.0
globalScale <- 50
NTrees <- 10 ## number of random trees to sample and to fit the model to
NHostSplits <- 15 ## desired number of nodes per timeBin
NMicrobeSplits <- 15 ## desired number of nodes per timeBin
modelform <- ~ ocean + ocean_area + reef_name + concatenated_date + colony_name + tissue_compartment + log_sequencing_depth_scaled
##

## Stan options
init_r <- 2
NCores <- NTrees
NChains <- 1 ## this is per tree; since I'm doing a large number of trees in parallel i'll just do one chain for each
NIterations <- 2500 ## will probably need >10,000? maybe start with 2, check convergence, double it, check, double, check, double, etc.?
max_treedepth <- 10 ## a warning will tell you if this needs to be increased
adapt_delta <- 0.8 ## increase this if you get 'divergences' - even one means your model fit sucks!
minMCSamples <- 2000 ## approximate number of Monte Carlo samples to save from the fit
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
    df2 <- droplevels(df1[(df1$tissue_compartment == 'T' |
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
    contrasts(df2$tissue_compartment) <- 'contr.sum'
    contrasts(df2$reef_name) <- 'contr.sum'
    contrasts(df2$concatenated_date) <- 'contr.sum'
    contrasts(df2$colony_name) <- 'contr.sum'
    levels(df2[,sampleTipKey])[levels(df2[,sampleTipKey]) == 'Homophyllia_hillae'] <- "Homophyllia_bowerbanki"
    levels(df2[,sampleTipKey])[levels(df2[,sampleTipKey]) == 'Pocillopora_eydouxi'] <- "Pocillopora_grandis"
    return(df2)
}

## import phylogenies
hostTree <- read.tree(hostTreePath)
hostTree <- .compressTipLabel(hostTree)
microbeTree <- read.tree(microbeTreePath)
##

## root the tree if it's unrooted
if(is.rooted(microbeTree)) {
    microbeTree.root <- reorder(microbeTree, order='pruningwise')
} else {
    microbeTree.root <- reorder(midpoint.root(microbeTree), order='pruningwise')
}
##

## add edge lengths to microbial tree if they're missing, and make tips contemporary if indicated
if(is.null(microbeTree.root$edge.length)) {
    microbeTree.root$edge.length <- rep(1,length(microbeTree.root$edge))
    microbeTree.root <- chronos(microbeTree.root, model = "discrete", control = chronos.control(nb.rate.cat = 1))
    class(microbeTree.root) <- 'phylo'
} else {
    microbeTree.root <- chronos(microbeTree.root)
    class(microbeTree.root) <- 'phylo'
}
##

##import data
fulltable <- t(read.table(fulltablePath, header=T, sep='\t', comment.char='', row.names=1, check.names=F))
newtable <- fulltable[-nrow(fulltable),] # last line in symportal file is not part of the table
mode(newtable) <- 'numeric'
##

## import mapping file
map <- read.table(mapfilePath,header=T,sep='\t',comment.char='',check.names=F)
rownames(map) <- map[,'#SampleID']
newmap <- filterfunction(map)
##

## merge data and mapping file
idx <- rownames(newtable)[rownames(newtable) %in% rownames(newmap) & rowSums(newtable) >= minCountSamp]
y.old <- newtable[idx, colnames(newtable) %in% microbeTree.root$tip.label]
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
finalMicrobeTree <- drop.tip(microbeTree.root, microbeTree.root$tip.label[!microbeTree.root$tip.label %in% colnames(y)])
y <- y[, finalMicrobeTree$tip.label]
##

## generate some summary numbers regarding microbes
microbeTips <- colnames(y)
NMicrobeTips <- length(microbeTips)
NIntMicrobeNodes <- finalMicrobeTree$Nnode
NMicrobeNodes <- NMicrobeTips + NIntMicrobeNodes - 1
microbeTreeDetails <- getTreeDetails(finalMicrobeTree)
microbeEdges <- finalMicrobeTree$edge.length[microbeTreeDetails$edgeOrder]
NMicrobeTimeBinsEst <- ceiling(NIntMicrobeNodes / NMicrobeSplits)
##

## create ancestry matrix for microbes
microbeAncestors <- createAncestryMat(NMicrobeNodes,
                                      finalMicrobeTree,
                                      NMicrobeTips,
                                      microbeTips)
colnames(microbeAncestors)[1:NMicrobeTips] <- rownames(microbeAncestors)[1:NMicrobeTips] <- paste0('t', microbeTips)
##

#divide total evolutionary time into chunks that contain approximately equal numbers of splits
microbeCutPoints <- findCutPoints(finalMicrobeTree,
                                  microbeTreeDetails$maxNHs,
                                  NMicrobeSplits,
                                  NMicrobeTimeBinsEst)

microbeTimeBins <- createTimeBins(microbeCutPoints,
                                  microbeTreeDetails$maxNHs,
                                  microbeTreeDetails$NHs,
                                  finalMicrobeTree,
                                  microbeTreeDetails$edgeOrder)

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

## prepare data for the model matrix
allfactors <- attr(terms.formula(modelform), "term.labels")
NFactors <- length(allfactors)
allfactorder <- sapply(allfactors, function(x) sum(gregexpr(':', x, fixed = TRUE)[[1]] > 0))
modelMat <- model.matrix(modelform, model.frame(newermap, na.action = NULL))
modelMat[is.na(modelMat)] <- 0
sumconts <- names(attr(modelMat, "contrasts")[attr(modelMat, "contrasts") == 'contr.sum'])
##

## create matrix relating each 'effect' (categorical and numeric) to the 'factor' that it belongs to
factLevelMat <- sapply(1:length(allfactors), function (j) {
    as.numeric(attr(modelMat, 'assign')[-1] == j)
})
colnames(factLevelMat) <- c(allfactors)

## ditch the intercept in this matrix because the 'global' intercept and 'main effects' of host and microbe are all estimated separately.
modelMat <- modelMat[, -1]
##

##
NEffects <- ncol(modelMat)
##

## rename factors that have 'sum contrasts' because by default they get arbitrary names (careful with interpretation of interactions... probably better to add them as separate 'main effects' produced by concatenation)
for(j in sumconts) {
    searchTerms <- paste0('^', j, 1:(nlevels(newermap[,j]) - 1), '$')
    replacementTerms <- paste0(j, levels(newermap[,j])[-nlevels(newermap[,j])])
    for(k in 1:length(searchTerms)) {
        colnames(modelMat) <- sub(searchTerms[[k]], replacementTerms[[k]], colnames(modelMat))
    }
}
rownames(factLevelMat) <- colnames(modelMat)
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

### starting here, generate multiple random samples of the map and trees (later to be summarized to define the time Bins and also to be separately used for replicate runs of the model)
NHostTimeBinsEst <- ceiling(length(levels(newermap[,sampleTipKey])) / NHostSplits)
sampleMap <- list()
hostTreesSampled <- list()
hostTreeDetails <- list()
hostCutPoints <- list()
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
    hostTreesSampled[[i]] <- ladderize(drop.tip(hostTree[[i]], hostTree[[i]]$tip.label[!hostTree[[i]]$tip.label %in% levels(sampleMap[[i]][,sampleTipKey])]))
    
    #get some tree stats for later use
    hostTreeDetails[[i]] <- getTreeDetails(hostTreesSampled[[i]])

	#divide total evolutionary time into chunks that contain approximately equal numbers of splits
    hostCutPoints[[i]] <- findCutPoints(hostTreesSampled[[i]],
                                        hostTreeDetails[[i]]$maxNHs,
                                        NHostSplits,
                                        NHostTimeBinsEst)
}
##

## average cut points among all the trees so the bins of time are consistent across replicates
meanHostBoundaries <- apply(sapply(hostCutPoints, function(x) x$boundaries), 1, mean)
meanHostBoundariesRounded <- round(meanHostBoundaries, 1)
##

## create matrices that relate the portion of each host branch that belongs to each time bin
hostTimeBins <- list()
for(i in 1:NTrees) {
    hostTimeBins[[i]] <- createTimeBins(hostCutPoints[[i]],
                                        hostTreeDetails[[i]]$maxNHs,
                                        hostTreeDetails[[i]]$NHs,
                                        hostTreesSampled[[i]],
                                        hostTreeDetails[[i]]$edgeOrder)
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


## collect data to feed to stan
standat <- list()
for (i in 1:NTrees) {
    standat[[i]] <- list(NSamples                       = NSamples,
                         NObs                           = NObs,
                         NMicrobeNodes                  = NMicrobeNodes,
                         NMicrobeTips                   = NMicrobeTips,
                         NFactors                       = NFactors,
                         NEffects                       = NEffects,
                         present                        = present,
                         sampleNames                    = sampleNames,
                         microbeTipNames                = microbeTipNames,
                         factLevelMat                   = factLevelMat,
                         modelMat                       = cbind(cbind(1, modelMat), hostAncestorsExpanded[[i]]),
                         microbeAncestorsT              = t(microbeAncestors),
                         microbeTipAncestorsT           = t(cbind(1, microbeAncestors[1:NMicrobeTips, ])),
                         hostAncestors                  = hostAncestors[[i]],
                         hostTipAncestors               = hostAncestors[[i]][1:NHostTips, ],
                         hostNodeHeights                = hostTreeDetails[[i]]$NHRel,
                         microbeNodeHeights             = microbeTreeDetails$NHRel,
                         microbeEdgeToBin               = t(microbeTimeBins$edgeToBin),
                         NMicrobeTimeBins               = microbeCutPoints$NTimeBins,
                         hostEdgeToBin                  = hostTimeBins[[i]]$edgeToBin,
                         NHostNodes                     = NHostNodes,
                         NHostTimeBins                  = hostCutPoints[[i]]$NTimeBins,
                         aveStDPriorExpect              = aveStDPriorExpect,
                         aveStDMetaPriorExpect          = aveStDMetaPriorExpect,
                         hostOUAlphaPriorExpect         = hostOUAlphaPriorExpect,
                         microbeOUAlphaPriorExpect      = microbeOUAlphaPriorExpect,
                         microbeEdges                   = microbeEdges,
                         globalScale                    = globalScale)
}

thin = max(1, floor(NIterations / minMCSamples))
NMCSamples <- NIterations / thin
warmup <- floor(NMCSamples / 2)
##

## save the inputs
dir.create(outdir, recursive=T)
save.image(file.path(outdir, 'setup.RData'))
save(microbeAncestors, file = file.path(outdir, 'microbeAncestors.RData'))
save(meanHostBoundaries, file = file.path(outdir, 'meanHostBoundaries.RData'))
save(hostAncestors, file = file.path(outdir, 'hostAncestors.RData'))
save(hostAncestorsExpanded, file = file.path(outdir, 'hostAncestorsExpanded.RData'))
save(hostTreesSampled, file = file.path(outdir, 'hostTreesSampled.RData'))
save(hostTreeDetails, file = file.path(outdir, 'hostTreeDetails.RData'))
##

cat('\nFitting model\n')
cat(paste0(as.character(Sys.time()),'\n'))

## run the model!
fit <- mclapply(1:NTrees,
    function(i) {
        setTimeLimit(timeLimit)
        tryCatch({
            stan(file     = modelPath,
            data     = standat[[i]],
            control  = list(adapt_delta   = adapt_delta,
                            max_treedepth = max_treedepth),
            iter     = NIterations,
            thin     = thin,
            chains   = NChains,
            seed     = seed,
            chain_id = (NChains * (i - 1) + (1:NChains)),
            pars     = c('rawMicrobeNodeEffects'),
            include  = FALSE,
            init_r   = init_r)
        }, error = function(e) NA)
    }, mc.preschedule = F,
       mc.cores       = NCores)

cat('\nSaving results\n')
cat(paste0(as.character(Sys.time()),'\n'))

save(fit, file = file.path(outdir,'fit.RData'))
##

## summarize the model fit
summarizeLcGLM()

## fin
