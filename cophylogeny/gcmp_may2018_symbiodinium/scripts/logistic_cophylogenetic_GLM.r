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
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


microbeTreePath <- 'raw_data/species_rooted_mapped_labels.tree' #ASTRAL microbial phylogeny
hostTreePath <- 'raw_data/combined_trees.newick' #set of Bayesian draws of host species phylogeny
mapfilePath <- 'raw_data/GCMP_symbiodinium_map2.txt' #mapping file
fulltablePath <- 'raw_data/species_table.txt' #250 bp deblur otu table output
modelPath <- 'scripts/logistic_cophylogenetic_GLM_varVar.stan' #stan model
seed <- 123

outdir <- file.path('output',gsub(':','-',gsub(' ', '_', Sys.time())))


## filtration options
minCountSamp <- 5 # minimum sequencing depth for a sample to be included
minSamps <- 1 # minimum number of samples that a sequence variant is present in at the above threshold for it to be included
##

## model options
aveStDPriorExpect <- 1.0
aveStDMetaPriorExpect <- 1.0
globalScale <- 50
NTrees <- 2 ## number of random trees to sample and to fit the model to
NSplits <- 15 ## desired number of nodes per host timeBin
ultrametricizeMicrobeTree <- TRUE
##

## Stan options
NChains <- 1 ## this is per tree; since I'm doing a large number of trees in parallel i'll just do one chain for each
NIterations <- 2500 ## will probably need >10,000? maybe start with 2, check convergence, double it, check, double, check, double, etc.?
max_treedepth <- 10 ## a warning will tell you if this needs to be increased
adapt_delta <- 0.8 ## increase this if you get 'divergences' - even one means your model fit sucks!
minMCSamples <- 2000 ## approximate number of Monte Carlo samples to save from the fit
##

## define the set of genera that we think our unidentified fungid samples could belong to
possibleFungidGenera <- c('Fungia_','Danafungia_','Cycloseris_','Pleuractis_')
## would be good here and for samples ID'd to genus to have some kind of weights for each possible species. i.e. it's really not possible that the acropora samples are A. cervicornis or palmata because their range is wrong, so even though I don't know which Acropora species exactly, I do have some info about the relative probabilities for many of them. similarly, some species might just be super rare in general and should likely be downweighted

sampleTipKey <- 'host_scientific_name'

filterfunction <- function(dfin) {
    levels(dfin[,sampleTipKey]) <- gsub(' ','_',levels(dfin[,sampleTipKey]))
    undupsamps <- levels(dfin$physical_sample_name)[sapply(levels(dfin$physical_sample_name), function(x) sum(x==dfin$physical_sample_name) == 1)]
    df1 <- droplevels(dfin[dfin$physical_sample_name %in% undupsamps,])
    df2 <- droplevels(df1[df1$tissue_compartment=='T' | df1$tissue_compartment=='S' | df1$tissue_compartment=='M' & df1[,sampleTipKey] != 'Unknown',])
    df2$tissue_compartment <- relevel(df2$tissue_compartment, 'T')
    contrasts(df2$ocean) <- 'contr.sum'
    contrasts(df2$complex_robust) <- 'contr.sum'
    contrasts(df2$host_scientific_name) <- 'contr.sum'
    contrasts(df2$tissue_compartment) <- 'contr.sum'
    contrasts(df2$reef_name) <- 'contr.sum'
    contrasts(df2$colony_name) <- 'contr.sum'
    contrasts(df2$host_genus) <- 'contr.sum'
	contrasts(df2$host_clade_sensu_fukami) <- 'contr.sum'
	df2$longitude <- as.numeric(as.character(df2$longitude))
    df2$latitude <- as.numeric(as.character(df2$latitude))
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
    if(ultrametricizeMicrobeTree) {
        microbeTree.root <- chronos(microbeTree.root, model = "discrete", control = chronos.control(nb.rate.cat = 1))
        class(microbeTree.root) <- 'phylo'
    }
} else {
    if(ultrametricizeMicrobeTree) {
        microbeTree.root <- chronos(microbeTree.root)
        class(microbeTree.root) <- 'phylo'
    }
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
y.old <- newtable[idx,colnames(newtable) %in% microbeTree.root$tip.label]
newermap <- newmap[idx,]
##

## convert data to presence/absence
ybinary <- apply(y.old, 2, function(x) x > 0)
mode(ybinary) <- 'numeric'
##

## filter microbes that aren't in the minimum number of samples
yb <- ybinary[, colSums(ybinary) >= minSamps]
microbeTree.root.Y <- drop.tip(microbeTree.root, microbeTree.root$tip.label[!microbeTree.root$tip.label %in% colnames(yb)])
yb <- yb[, microbeTree.root.Y$tip.label]
##

## generate some summary numbers regarding microbes
microbeTips <- colnames(yb)
NMicrobeTips <- length(microbeTips)
NIntMicrobeNodes <- microbeTree.root.Y$Nnode
NMicrobeNodes <- NMicrobeTips + NIntMicrobeNodes - 1
microbeEdgeOrder <- order(microbeTree.root.Y$edge[,2])
microbeEdges <- microbeTree.root.Y$edge.length[microbeEdgeOrder]
##

## melt the data into long format to feed to model, and generate summary numbers about samples
senddat <- melt(yb, varnames=c('sample','tip'), value.name='present')
sampleNames <- as.numeric(factor(senddat[,1]))
microbeTipNames <- as.numeric(factor(senddat[,2], levels=microbeTips))
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

## create ancestry matrix for microbes
microbeAncestors <- matrix(0, NMicrobeNodes + 1, NMicrobeNodes + 1)
for(node in 1:(NMicrobeNodes + 1)) {
    microbeAncestors[node, ] <- as.numeric(1:(NMicrobeNodes + 1) %in% c(Ancestors(microbeTree.root.Y, node), node))
}
colnames(microbeAncestors) <- rownames(microbeAncestors) <- paste0('i',1:(NMicrobeNodes + 1))
colnames(microbeAncestors)[1:NMicrobeTips] <- rownames(microbeAncestors)[1:NMicrobeTips] <- paste0('t',colnames(yb))
microbeAncestors <- microbeAncestors[-(NMicrobeTips + 1), -(NMicrobeTips + 1)]
##

## prepare data for the model matrix
modelform <- ~ ocean + reef_name + colony_name + tissue_compartment + log_sequencing_depth_scaled
allfactors <- attr(terms.formula(modelform), "term.labels")
NFactors <- length(allfactors)
allfactorder <- sapply(allfactors, function(x) sum(gregexpr(':', x, fixed=TRUE)[[1]] > 0))
modelMat <- model.matrix(modelform, model.frame(newermap,na.action=NULL))
modelMat[is.na(modelMat)] <- 0
##

## rename factors that have 'sum contrasts' because by default they get arbitrary names
sumconts <- names(attr(modelMat, "contrasts")[attr(modelMat, "contrasts")=='contr.sum'])
for(j in sumconts) {
    colnames(modelMat)[grep(j,colnames(modelMat))] <- paste0(j,levels(newermap[,j])[-nlevels(newermap[,j])]) ##this will not work if there are 'interaction effects' specified in the model formula!!!
}
##

## ditch the intercept in this matrix because the 'global' intercept and 'main effects' of host and microbe are all estimated separately.
modelMat <- modelMat[,-1]
##

##
NEffects <- ncol(modelMat)
##

## create matrix relating each 'effect' (categorical and numeric) to the 'factor' that it belongs to
factLevelMat <- matrix(NA, NEffects, NFactors)
colnames(factLevelMat) <- c(allfactors)
rownames(factLevelMat) <- colnames(modelMat)
remainder <- colnames(modelMat)
for(fact in names(sort(allfactorder,decreasing=T))) {
    matches <- remainder
    for(y in strsplit(fact,':')[[1]]) {
        matches <- grep(y,matches,value=T)
    }
    factLevelMat[,fact] <- as.numeric(colnames(modelMat) %in% matches)
    remainder <- remainder[!remainder %in% matches]
}
##

## identify unique host species in the data, and replace spaces with underscores
study.species <- gsub(' ', '_', levels(newermap[,sampleTipKey]))
##

## identify the host species in the data that do not exist in the template tree
study.species.missing <- study.species[!study.species %in% c(attr(hostTree, "TipLabel"),'Fungid_sp','not_applicable')]
generaOfUnknowns <- sapply(study.species.missing, function(x) strsplit(x, '_')[[1]][[1]])
##

### starting here, generate multiple random samples of the map and trees (later to be summarized to define the time Bins and also to be separately used for replicate runs of the model)
NTimeBins <- ceiling(length(levels(newermap[,sampleTipKey])) / NSplits)
sampleMap <- list()
hostTreesSampled <- list()
splittimes <- list()
boundaries <- matrix(NA,nrow=NTrees, ncol=NTimeBins-1)
for(i in 1:NTrees) {
    sampleMap[[i]] <- newermap
    levels(sampleMap[[i]][,sampleTipKey])[levels(sampleMap[[i]][,sampleTipKey]) == 'Fungid_sp'] <- sample(grep(paste0(paste(possibleFungidGenera,collapse='|'),'_'), attr(hostTree, "TipLabel"), value=T), 1) ##assign unidentified Fungids to a random member of the group independently for each tree
    for (j in unique(generaOfUnknowns)) {
        possibleGenera <- attr(hostTree, "TipLabel")[!attr(hostTree, "TipLabel") %in% levels(sampleMap[[i]][,sampleTipKey])]
        levels(sampleMap[[i]][,sampleTipKey])[levels(sampleMap[[i]][,sampleTipKey]) %in% grep(paste0(j,'_'),study.species.missing,value=T)] <- sample(grep(paste0(j,'_'), possibleGenera, value=T), sum(generaOfUnknowns == j)) ## assign unidentified species to a random member of their genus independently for each tree
    }

	#filter the tree only contain the sampled (or assigned) species
	hostTreesSampled[[i]] <- ladderize(drop.tip(hostTree[[i]],hostTree[[i]]$tip.label[!hostTree[[i]]$tip.label %in% levels(sampleMap[[i]][,sampleTipKey])]))

	#divide total evolutionary time into chunks that contain approximately equal numbers of splits
	lttHostTreesSampled <- ltt(hostTreesSampled[[i]],log.lineages=F,plot=F)
	temp <- max(nodeHeights(hostTreesSampled[[i]])) - lttHostTreesSampled$times[-length(lttHostTreesSampled$times)]
	splittimes[[i]] <- split(temp, ceiling(seq_along(temp)/NSplits))
    if(length(splittimes[[i]][[NTimeBins]]) < NSplits/2) {
        splittimes[[i]][[NTimeBins - 1]] <- c(splittimes[[i]][[NTimeBins - 1]], splittimes[[i]][[NTimeBins]])
        splittimes[[i]][[NTimeBins]] <- NULL
    }
    
    #cut points in each phylogeny that would result in approximately equal numbers of splits per bin of time
    boundaries[i,] <- sapply(2:NTimeBins, function(x) max(splittimes[[i]][[x]]))

}
##

## average cut points among all the trees so the bins of time are consistent across replicates
meanBoundaries <- apply(boundaries,2,mean)
meanBoundariesRounded <- round(meanBoundaries,1)
##

## create matrices that relate the portion of each host branch that belongs to each time bin
timeBinSizes <- list()
relativeTimeBinSizes <- list()
edgeToBin <- list()
hostEdgeOrder <- list()
nh <- list()
maxNH <- list()
nhRel <- list()
for(i in 1:NTrees) {
    nh[[i]] <- nodeHeights(hostTreesSampled[[i]])
    maxNH[[i]] <- max(nh[[i]])
    timeBinSizes[[i]] <- sapply(2:(length(meanBoundaries)+2),function(x) c(maxNH[[i]],meanBoundaries,0)[x-1] - c(maxNH[[i]],meanBoundaries,0)[x]) #size in millions of years of each bin (consistent across replicates /except/ for the oldest bin
    relativeTimeBinSizes[[i]] <- timeBinSizes[[i]] / sum(timeBinSizes[[i]]) #proportion of total tree height belonging to each bin (differs for each replicate due to the varying sizes of the oldest bin)
    nd <- maxNH[[i]] - nh[[i]]
    edgeToBin[[i]] <- matrix(NA,ncol=NTimeBins,nrow=nrow(hostTreesSampled[[i]]$edge)) #create a matrix for each replicate that describes the proportion of each time bin that each branch exists for
    for(j in 1:NTimeBins) {
        if(j == 1) {
            allin <- which(nd[,2] >= meanBoundaries[j])
			allout <- which(nd[,1] <= meanBoundaries[j])
            cedge <- which((nd[,1] > meanBoundaries[j]) & (nd[,2] < meanBoundaries[j]))
            edgeToBin[[i]][cedge,j] <- nd[cedge,1] - meanBoundaries[j]
        } else if(j == NTimeBins) {
            allin <- which(nd[,1] <= meanBoundaries[j-1])
            allout <- which(nd[,2] >= meanBoundaries[j-1])
            cedge <- which((nd[,1] > meanBoundaries[j-1]) & (nd[,2] < meanBoundaries[j-1]))
            edgeToBin[[i]][cedge,j] <- meanBoundaries[j-1] - nd[cedge,2]
        } else {
            allin <- which((nd[,1] <= meanBoundaries[j-1]) & (nd[,2] >= meanBoundaries[j]))
            allout <- which((nd[,1] <= meanBoundaries[j]) | (nd[,2] >= meanBoundaries[j-1]))
            cedge1 <- which((nd[,1] < meanBoundaries[j-1]) & (nd[,1] > meanBoundaries[j]) & (nd[,2] < meanBoundaries[j]))
            edgeToBin[[i]][cedge1,j] <- nd[cedge1,1] - meanBoundaries[j]
            cedge2 <- which((nd[,1] > meanBoundaries[j-1]) & (nd[,2] < meanBoundaries[j-1]) & (nd[,2] > meanBoundaries[j]))
            edgeToBin[[i]][cedge2,j] <- meanBoundaries[j-1] - nd[cedge2,2]
            cedge3 <- which((nd[,1] > meanBoundaries[j-1]) & (nd[,2] < meanBoundaries[j]))
            edgeToBin[[i]][cedge3,j] <- meanBoundaries[j-1] - meanBoundaries[j]
        }
        edgeToBin[[i]][allin,j] <- hostTreesSampled[[i]]$edge.length[allin]
        edgeToBin[[i]][allout,j] <- 0

		#edgeToBin[[i]][,j] <- edgeToBin[[i]][,j] / timeBinSizes[[i]][j]
    }
    edgeToBin[[i]] <- edgeToBin[[i]] / maxNH[[i]] #this replaces the commented line above

    hostEdgeOrder[[i]] <- order(hostTreesSampled[[i]]$edge[,2])
    rownames(edgeToBin[[i]]) <- hostTreesSampled[[i]]$edge[,2]
    edgeToBin[[i]] <- edgeToBin[[i]][hostEdgeOrder[[i]],]
    
    nhRel[[i]] <- nh[[i]] / maxNH[[i]]
    rownames(nhRel[[i]]) <- hostTreesSampled[[i]]$edge[,2]
    nhRel[[i]] <- nhRel[[i]][hostEdgeOrder[[i]],]
}
##

## create ancestry matrices for each host tree
NHostTips <- length(hostTreesSampled[[1]]$tip.label)
NIntHostNodes <- hostTreesSampled[[1]]$Nnode
NHostNodes <- NIntHostNodes + NHostTips - 1
hostAncestors <- list()
hostAncestorsExpanded <- list()
for (i in 1:NTrees) {
	hostAncestors[[i]] <- matrix(NA, NHostNodes+1, NHostNodes+1)
	for(node in 1:(NHostNodes+1)) {
	    hostAncestors[[i]][node, ] <- as.numeric(1:(NHostNodes+1) %in% c(Ancestors(hostTreesSampled[[i]], node), node))
	}
    colnames(hostAncestors[[i]]) <- rownames(hostAncestors[[i]]) <- paste0('i',1:(NHostNodes+1))
    colnames(hostAncestors[[i]])[1:NHostTips] <- rownames(hostAncestors[[i]])[1:NHostTips] <- hostTreesSampled[[i]]$tip.label
    hostAncestors[[i]] <- hostAncestors[[i]][-(NHostTips + 1), -(NHostTips + 1)]
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
                         hostNodeHeights                = nhRel[[i]],
                         edgeToBin                      = edgeToBin[[i]],
                         NHostNodes                     = NHostNodes,
                         NTimeBins                      = NTimeBins,
                         aveStDPriorExpect              = aveStDPriorExpect,
                         aveStDMetaPriorExpect          = aveStDMetaPriorExpect,
                         microbeEdges                   = microbeEdges,
                         globalScale                    = globalScale)
}

thin = max(1, floor(NIterations / minMCSamples))
NMCSamples <- NIterations / thin
warmup <- floor(NMCSamples / 2)
##

## save the inputs
dir.create(outdir, recursive=T)
save.image(file.path(outdir,'setup.RData'))
save(microbeAncestors,file=file.path(outdir,'microbeAncestors.RData'))
save(meanBoundaries,file=file.path(outdir,'meanBoundaries.RData'))
save(hostAncestors,file=file.path(outdir,'hostAncestors.RData'))
save(hostAncestorsExpanded,file=file.path(outdir,'hostAncestorsExpanded.RData'))
save(hostTreesSampled,file=file.path(outdir,'hostTreesSampled.RData'))
save(timeBinSizes,file=file.path(outdir,'timeBinSizes.RData'))
##

## run the model!
fit <- mclapply(1:NTrees, function(i) stan(file     = modelPath,
                                           data     = standat[[i]],
                                           control  = list(adapt_delta   = adapt_delta,
                                                           max_treedepth = max_treedepth),
                                           iter     = NIterations,
                                           thin     = thin,
                                           chains   = NChains,
                                           seed     = seed,
                                           chain_id = (NChains * (i - 1) + (1:NChains)),
                                           pars     = c('rawMicrobeNodeEffects'),
                                           include  = FALSE))

save(fit, file=file.path(outdir,'fit.RData'))
##

## summarize the results separately for each sampled host tree
for(i in 1:NTrees) {
    
    check_hmc_diagnostics(fit[[i]])

    currplotdir <- file.path(outdir,paste0('tree_',i),'plots')
    currtabledir <- file.path(outdir,paste0('tree_',i),'tables')
    currdatadir <- file.path(outdir,paste0('tree_',i),'data')

    dir.create(currplotdir, recursive=T)
    dir.create(currtabledir, recursive=T)
    dir.create(currdatadir, recursive=T)

    ## plot the sampled tree with the time bins marked
    pdf(file=file.path(currplotdir,'sampledTree.pdf'), width=25, height=15)
    plot(hostTreesSampled[[i]],cex=0.75)
    for(age in max(nodeHeights(hostTreesSampled[[i]])) - meanBoundaries) {
        lines(x = c(age, age), y = c(1, length(hostTreesSampled[[i]]$tip.label)), lwd=1)
    }
    graphics.off()
    ##

    ## variance partitioning
    stDProps <- extract(fit[[i]], pars='stDProps')[[1]]
    colnames(stDProps) <- c(paste0('ADiv.',colnames(factLevelMat)), paste0('Specificity.',colnames(factLevelMat)), 'ADiv.host', 'host.specificity', 'microbe.prevalence')

    pdf(file=file.path(currplotdir,'scalesboxes.pdf'), width=25, height=15)
    boxplot(stDProps, cex.axis=0.5, las=2)
    graphics.off()

    save(stDProps,file=file.path(currdatadir,'stDProps.RData'))
    ##
    
    ## proportion of variance explained by each time bin
    metaVarProps <- extract(fit[[i]], pars='metaVarProps')[[1]]
    colnames(metaVarProps) <- c('prevalence','adiv','specificty')
    pdf(file=file.path(currplotdir,'metaVarProps_boxes.pdf'), width=25, height=15)
    boxplot(metaVarProps, cex.axis=0.5, las=2)
    graphics.off()
    save(metaVarProps,file=file.path(currdatadir,'metaVarProps.RData'))
    ##
    
    ## compare rates of host evolution in each time bin
    relativeEvolRates <- exp(extract(fit[[i]], pars='logRelativeEvolRates')[[1]])
    colnames(relativeEvolRates) <- c(paste0('before ',meanBoundariesRounded[1],' mya'), paste0(meanBoundariesRounded[1],' - ',meanBoundariesRounded[2],' mya'), paste0(meanBoundariesRounded[2],' - ',meanBoundariesRounded[3],' mya'), paste0(meanBoundariesRounded[3],' - ',meanBoundariesRounded[4],' mya'), paste0(meanBoundariesRounded[4],' - present'))

    pdf(file=file.path(currplotdir,'evolRatesRelToMRCA.pdf'), width=25, height=15)
    boxplot(relativeEvolRates, xlab='Time Period', ylab='Rate of Evolution Relative to MRCA')
    graphics.off()
    
    pdf(file=file.path(currplotdir,'logEvolRatesRelToMRCA.pdf'), width=25, height=15)
    boxplot(log(relativeEvolRates), xlab='Time Period', ylab='Log Rate of Evolution Relative to MRCA')
    graphics.off()

	save(relativeEvolRates,file=file.path(currdatadir,'relativeEvolRates.RData'))
    
    ##
    
    ## summarize effects
    currsubtabledir <- file.path(currtabledir, 'nodeEffects')
    dir.create(currsubtabledir, recursive=T)

    scaledMicrobeNodeEffects <- array(extract(fit[[i]], pars='scaledMicrobeNodeEffects', permuted=F, inc_warmup=T),
                                      dim=c(NMCSamples,
                                            NChains,
                                            NEffects + NHostNodes + 1,
                                            NMicrobeNodes + 1),
                                      dimnames=list(sample  = NULL,
                                                    chain   = NULL,
                                                    effect  = c('microbePrevalence', colnames(modelMat)[1:NEffects], colnames(hostAncestors[[i]])),
                                                    taxnode = c('alphaDiversity', colnames(microbeAncestors))))
                                                    
    save(scaledMicrobeNodeEffects, file = file.path(currdatadir, 'scaledMicrobeNodeEffects.RData'))
                                                    
    baseLevelEffects <- array(NA,
                              dim=c(NMCSamples,
                                    NChains,
                                    length(sumconts),
                                    NMicrobeNodes + 1),
                              dimnames=list(sample  = NULL,
                                            chain   = NULL,
                                            effect  = sumconts,
                                            taxnode = c('alphaDiversity', colnames(microbeAncestors))))
    for(j in 1:NMCSamples) {
        for(k in 1:NChains) {
            for(m in sumconts) {
                baseLevelEffects[j,k,m,] <- -colSums(scaledMicrobeNodeEffects[j,k,rownames(factLevelMat)[factLevelMat[,m]==1],])
            }
        }
    }
    
    save(baseLevelEffects, file = file.path(currdatadir, 'baseLevelEffects.RData'))

    for(l in 1:(NEffects + NHostNodes + 1)) {
        yeah <- monitor(array(scaledMicrobeNodeEffects[,,l,],
                              dim = c(NMCSamples, NChains, NMicrobeNodes + 1)),
                        warmup = warmup,
                        probs = c(0.05, 0.95))
        rownames(yeah) <- c('alphaDiversity', rownames(microbeAncestors))
        cat('\t', file = file.path(currsubtabledir, paste0(dimnames(scaledMicrobeNodeEffects)[[3]][l], '.txt')))
        write.table(yeah, file = file.path(currsubtabledir, paste0(dimnames(scaledMicrobeNodeEffects)[[3]][l], '.txt')), sep='\t', quote=F,append=T)
    }
    
    for(m in sumconts) {
        yeah <- monitor(array(baseLevelEffects[,,m,],
                              dim = c(NMCSamples, NChains, NMicrobeNodes + 1)),
                        warmup = warmup,
                        probs = c(0.05, 0.95))
        rownames(yeah) <- c('alphaDiversity', rownames(microbeAncestors))
        cat('\t', file = file.path(currsubtabledir, paste0(m, levels(newermap[,m])[nlevels(newermap[,m])], '.txt')))
        write.table(yeah, file = file.path(currsubtabledir, paste0(m, levels(newermap[,m])[nlevels(newermap[,m])], '.txt')), sep='\t', quote=F,append=T)
    }
    ##
    
    ## see if any pairs of clades have higher variance among their descendants (maybe suggesting codiversification)
    currsubtabledir <- file.path(currtabledir, 'codivEffects')
    dir.create(currsubtabledir, recursive=T)
    
    sums <- summary(fit[[i]], pars='phyloLogVarMultRaw', probs=c(0.05,0.95), use_cache = F)

    sums3d <- array(NA, dim=c(NHostNodes - NHostTips, NMicrobeNodes - NMicrobeTips, ncol(sums$summary)))
    for(effect in 1:(NHostNodes - NHostTips)) {
        sums3d[effect,,] <- sums$summary[(effect-1) * (NMicrobeNodes - NMicrobeTips) + (1:(NMicrobeNodes - NMicrobeTips)),]
    }
    dimnames(sums3d) <- list(colnames(hostAncestors[[i]])[(NHostTips + 1):NHostNodes], rownames(microbeAncestors)[(NMicrobeTips + 1):NMicrobeNodes], colnames(sums$summary))
    factorfilenames <- colnames(hostAncestors[[i]])[(NHostTips + 1):NHostNodes]

    for(effect in 1:(NHostNodes - NHostTips)) {
        cat('\t', file=file.path(currsubtabledir, paste0(factorfilenames[[effect]],'.txt')))
        write.table(sums3d[effect,,], file=file.path(currsubtabledir, paste0(factorfilenames[[effect]],'.txt')), sep='\t', quote=F,append=T)
    }
    ##
    
    ## summarize the mean branch lengths of the hosts
    sums <- summary(fit[[i]], pars='hostScales', probs=c(0.05,0.95), use_cache = F)
    newEdges <- sums$summary[,'mean']^2
    hostTreesSampled[[i]]$edge.length <- newEdges[order(hostEdgeOrder[[i]])]
    pdf(file=file.path(currplotdir,'hostTreeWEstimatedEdgeLengths.pdf'), width=25, height=15)
    plot(hostTreesSampled[[i]], cex=0.5)
    graphics.off()
    ##
}
##

## summarize results for parameters that can be interpretted across all sampled host trees
currplotdir <- file.path(outdir,'alltrees','plots')
currtabledir <- file.path(outdir,'alltrees','tables')
currdatadir <- file.path(outdir,'alltrees','data')

dir.create(currplotdir, recursive=T)
dir.create(currtabledir, recursive=T)
dir.create(currdatadir, recursive=T)

allfit <- sflist2stanfit(fit)

stDProps <- extract(allfit, pars='stDProps')[[1]]
colnames(stDProps) <- c(paste0('ADiv.',colnames(factLevelMat)), paste0('Specificity.',colnames(factLevelMat)), 'ADiv.host', 'host.specificity', 'microbe.prevalence')

pdf(file=file.path(currplotdir,'scalesboxes.pdf'), width=25, height=15)
boxplot(stDProps, cex.axis=0.5, las=2)
graphics.off()

save(stDProps,file=file.path(currdatadir,'stDProps.RData'))

relativeEvolRates <- exp(extract(allfit, pars='logRelativeEvolRates')[[1]])
colnames(relativeEvolRates) <- c(paste0('before ',meanBoundariesRounded[1],' mya'), paste0(meanBoundariesRounded[1],' - ',meanBoundariesRounded[2],' mya'), paste0(meanBoundariesRounded[2],' - ',meanBoundariesRounded[3],' mya'), paste0(meanBoundariesRounded[3],' - ',meanBoundariesRounded[4],' mya'), paste0(meanBoundariesRounded[4],' - present'))

pdf(file=file.path(currplotdir,'evolRatesRelToMRCA.pdf'), width=25, height=15)
boxplot(relativeEvolRates, xlab='Time Period', ylab='Rate of Evolution Relative to MRCA')
graphics.off()

pdf(file=file.path(currplotdir,'logEvolRatesRelToMRCA.pdf'), width=25, height=15)
boxplot(log(relativeEvolRates), xlab='Time Period', ylab='Log Rate of Evolution Relative to MRCA')
graphics.off()

save(relativeEvolRates,file=file.path(currdatadir,'relativeEvolRates.RData'))

## summarize the mean branch lengths of the microbes
sums <- summary(allfit, pars='microbeScales', probs=c(0.05,0.95), use_cache = F)
newEdges <- sums$summary[,'mean']^2
microbeTree.root.Y$edge.length <- newEdges[order(microbeEdgeOrder)]
pdf(file=file.path(currplotdir,'microbeTreeWEstimatedEdgeLengths.pdf'), width=25, height=15)
plot(microbeTree.root.Y, cex=0.5)
graphics.off()
##

## summarize effects
currsubtabledir <- file.path(currtabledir, 'nodeEffects')
dir.create(currsubtabledir, recursive=T)

scaledMicrobeNodeEffects <- array(extract(allfit, pars='scaledMicrobeNodeEffects', permuted=F, inc_warmup=T),
                                  dim=c(NMCSamples,
                                        NChains * NTrees,
                                        NEffects + NHostNodes + 1,
                                        NMicrobeNodes + 1),
                                  dimnames=list(sample  = NULL,
                                                chain   = NULL,
                                                effect  = c('microbePrevalence', colnames(modelMat)[1:NEffects], colnames(hostAncestors[[i]])),
                                                taxnode = c('alphaDiversity', colnames(microbeAncestors))))
                                                
save(scaledMicrobeNodeEffects, file = file.path(currdatadir, 'scaledMicrobeNodeEffects.RData'))
                                                
baseLevelEffects <- array(NA,
                          dim=c(NMCSamples,
                                NChains * NTrees,
                                length(sumconts),
                                NMicrobeNodes + 1),
                          dimnames=list(sample  = NULL,
                                        chain   = NULL,
                                        effect  = sumconts,
                                        taxnode = c('alphaDiversity', colnames(microbeAncestors))))
for(j in 1:NMCSamples) {
    for(k in 1:(NChains * NTrees)) {
        for(m in sumconts) {
            baseLevelEffects[j,k,m,] <- -colSums(scaledMicrobeNodeEffects[j,k,rownames(factLevelMat)[factLevelMat[,m]==1],])
        }
    }
}

save(baseLevelEffects, file = file.path(currdatadir, 'baseLevelEffects.RData'))

for(l in 1:(NEffects + NHostNodes + 1)) {
    yeah <- monitor(array(scaledMicrobeNodeEffects[,,l,],
                          dim = c(NMCSamples, NChains * NTrees, NMicrobeNodes + 1)),
                    warmup = warmup,
                    probs = c(0.05, 0.95))
    rownames(yeah) <- c('alphaDiversity', rownames(microbeAncestors))
    cat('\t', file = file.path(currsubtabledir, paste0(dimnames(scaledMicrobeNodeEffects)[[3]][l], '.txt')))
    write.table(yeah, file = file.path(currsubtabledir, paste0(dimnames(scaledMicrobeNodeEffects)[[3]][l], '.txt')), sep='\t', quote=F,append=T)
}

for(m in sumconts) {
    yeah <- monitor(array(baseLevelEffects[,,m,],
                          dim = c(NMCSamples, NChains * NTrees, NMicrobeNodes + 1)),
                    warmup = warmup,
                    probs = c(0.05, 0.95))
    rownames(yeah) <- c('alphaDiversity', rownames(microbeAncestors))
    cat('\t', file = file.path(currsubtabledir, paste0(m, levels(newermap[,m])[nlevels(newermap[,m])], '.txt')))
    write.table(yeah, file = file.path(currsubtabledir, paste0(m, levels(newermap[,m])[nlevels(newermap[,m])], '.txt')), sep='\t', quote=F,append=T)
}
##

allfit <- c(fit, allfit)

## make some diagnostic plots
for(i in 1:(NTrees + 1)) {
    currdiagnosticdir <- file.path(outdir, 'diagnostics', if(i <= NTrees) {paste0('tree_', i)} else {'allfit'})
    dir.create(currdiagnosticdir, recursive=T)
    
    pars <- c('stDProps', 'hostMetaVarProps', 'timeBinMetaVar', 'metaVarProps', 'aveStD')
    
    pdf(file=file.path(currdiagnosticdir,'pairs_grabBag.pdf'), width=50, height=50)
    pairs(allfit[[i]], pars=pars)
    graphics.off()
    
    pdf(file=file.path(currdiagnosticdir,'pairs_grabBagLog.pdf'), width=50, height=50)
    pairs(allfit[[i]], pars=pars, log=T)
    graphics.off()
    
    pars <- paste0('phyloLogVarMultPrev[',1:20,']')
    
    pdf(file=file.path(currdiagnosticdir,'pairs_phyloLogVarMultPrev.pdf'), width=50, height=50)
    pairs(allfit[[i]], pars=pars, log=T)
    graphics.off()
    
}
##

## fin
