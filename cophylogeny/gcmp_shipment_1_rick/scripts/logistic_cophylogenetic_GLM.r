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



microbeTreePath <- 'raw_data/gg_constrained_fastttree.tre' #ML microbial phylogeny
hostTreePath <- 'raw_data/combined_trees.newick' #set of Bayesian draws of host species phylogeny
mapfilePath <- 'raw_data/GCMP_EMP_map_r28.txt' #mapping file
fulltablePath <- 'raw_data/reference-hit.txt' #250 bp deblur otu table output
taxAssignmentPath <- 'raw_data/reference-hit.seqs_tax_assignments.txt' #greegenes taxonomy
modelPath <- 'scripts/logistic_cophylogenetic_GLM_varVar.stan' #stan model
seed <- 123
timeLimit <- 30 * 24 * 60 * 60

outdir <- file.path('output',gsub(':', '-', gsub(' ', '_', Sys.time())))

## filtration options
minCountSamp <- 100 # minimum sequencing depth for a sample to be included
minPercent <- 0 # minimum percent of a sample composed of a sequence variant for it to pass the below filter
minSamps <- 1 # minimum number of samples that a sequence variant is present in at the above threshold for it to be included
##

## model options
aveStDPriorExpect <- 1.0
aveStDMetaPriorExpect <- 1.0
hostOUAlphaPriorExpect <- 1.0
microbeOUAlphaPriorExpect <- 1.0
globalScale <- 50
NTrees <- 1 ## number of random trees to sample and to fit the model to
NSplits <- 30 ## desired number of nodes per host timeBin
##

## Stan options
init_r <- 2
NCores <- 1
NChains <- 1 ## this is per tree; since I'm doing a large number of trees in parallel i'll just do one chain for each
NIterations <- 1500 ## will probably need >10,000? maybe start with 2, check convergence, double it, check, double, check, double, etc.?
max_treedepth <- 10 ## a warning will tell you if this needs to be increased
adapt_delta <- 0.8 ## increase this if you get 'divergences' - even one means your model fit sucks!
minMCSamples <- 2000 ## approximate number of Monte Carlo samples to save from the fit
##

## define the set of genera that we think our unidentified fungid samples could belong to
possibleFungidGenera <- c('Fungia_','Danafungia_','Cycloseris_','Pleuractis_')
## would be good here and for samples ID'd to genus to have some kind of weights for each possible species. i.e. it's really not possible that the acropora samples are A. cervicornis or palmata because their range is wrong, so even though I don't know which Acropora species exactly, I do have some info about the relative probabilities for many of them. similarly, some species might just be super rare in general and should likely be downweighted

taxaToExclude <- c('Aiptasia_sp','Entacmaea_quadricolor','Heteractis_aurora','Macrorhynchia_philippina','Mnemiopsis_sp','Sarcophyton_sp','Sinularia_polydactyla','Stichodactyla__gigantea','Xenia_umbellata') # don't currently have a way to incorporate these into the phylogeny...

sampleTipKey <- 'host_scientific_name'

filterfunction <- function(dfin) {
    levels(dfin[,sampleTipKey]) <- gsub(' ','_',levels(dfin[,sampleTipKey]))
    df2 <- droplevels(dfin[(dfin$tissue_compartment=='T' | dfin$tissue_compartment=='S' | dfin$tissue_compartment=='M') & !grepl(paste(c('Unknown|Missing',taxaToExclude),collapse='|'),dfin[,sampleTipKey],ignore.case=T) ,])
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
    levels(df2[,sampleTipKey])[levels(df2[,sampleTipKey]) == 'Pseudosiderastrea_tayami'] <- "Pseudosiderastrea_tayamai"
    return(df2)
}

subsetdetails <- list(taxlevel = 'Order', taxname = 'o__Rickettsiales')

## import host phylogenies
hostTree <- read.tree(hostTreePath)
hostTree <- .compressTipLabel(hostTree)
##

## import mapping file
map <- read.table(mapfilePath,header=T,sep='\t',comment.char='',check.names=F)
rownames(map) <- map[,'#SampleID']
newmap <- filterfunction(map)
##

## import data; merge data and mapping file
fulltable <- t(read.table(fulltablePath, header=T, sep='\t', skip=1, comment.char='', row.names=1, check.names=F))
idx <- rownames(fulltable)[rownames(fulltable) %in% rownames(newmap) & rowSums(fulltable) >= minCountSamp]
y.old <- fulltable[idx, ]
newermaptemp <- droplevels(newmap[idx,])
##

## define contrasts
newermap <- contrastfunction(newermaptemp)
##

y.old.filt <- t(apply(y.old,1,function(x) {
    temp <- x/sum(x)
    return(temp > minPercent & !is.na(temp))
}))
y.old.binary <- apply(y.old,2,function(x) x > 0)
mode(y.old.binary) <- 'numeric'

y.binary.filtered <- y.old.binary[,colSums(y.old.filt) >= minSamps]


taxdat <- read.table(taxAssignmentPath,sep='\t',stringsAsFactors=F, row.names=1)
x <- strsplit(taxdat[,1],'; ')
most <- max(sapply(x,length))
parsedtax <- lapply(x,function(x) {length(x) <- most; return(x)})
tax <- do.call('rbind',parsedtax)
rownames(tax) <- rownames(taxdat)
colnames(tax) <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')


## import microbe tree
microbeTree <- read.tree(microbeTreePath)
##

## filter the dataset to include only the monophyletic group of Rickettsiales plus some sampled outgroups
rickAnc <- getMRCA(microbeTree, rownames(tax[tax[,'Family'] == 'f__Anaplasmataceae' | tax[,'Family'] == 'f__Rickettsiaceae' & !is.na(tax[,'Family']),]))
rickTips <- microbeTree$tip.label[Descendants(microbeTree, rickAnc)[[1]]]

otherRick <- rownames(tax[tax[,'Order']=='o__Rickettsiales' & !is.na(tax[,'Order']),])
myOtherRick <- colnames(y.binary.filtered)[colnames(y.binary.filtered) %in% otherRick]
myOtherRickSampled <- sample(myOtherRick[!myOtherRick %in% rickTips], ceiling(length(rickTips) / 20))

myOthersSampled <- sample(colnames(y.binary.filtered)[!colnames(y.binary.filtered) %in% c(rickTips, myOtherRickSampled)], ceiling(length(rickTips) / 20))

y.binary.filtered.subgroup <- y.binary.filtered[,colnames(y.binary.filtered) %in% c(rickTips, myOtherRickSampled,  myOthersSampled)]

microbeTree.Y <- drop.tip(microbeTree, microbeTree$tip.label[!microbeTree$tip.label %in% colnames(y.binary.filtered.subgroup)])
##


## root the tree if it's unrooted
if(is.rooted(microbeTree.Y)) {
    microbeTree.Y.root <- reorder(microbeTree.Y, order='pruningwise')
} else {
    microbeTree.Y.root <- reorder(midpoint.root(microbeTree.Y), order='pruningwise')
}
##

## add edge lengths to microbial tree if they're missing, and make tips contemporary
if(is.null(microbeTree.Y.root$edge.length)) {
    microbeTree.Y.root$edge.length <- rep(1, length(microbeTree.Y.root$edge))
    microbeTree.Y.root <- chronos(microbeTree.Y.root, control = chronos.control(tol = 1e-12)) #default tolerance of 1e-8 seems to lead to underflow and edge lengths of 0.
    class(microbeTree.Y.root) <- 'phylo'
} else {
    microbeTree.Y.root <- chronos(microbeTree.Y.root, control = chronos.control(tol = 1e-12))
    class(microbeTree.Y.root) <- 'phylo'
}
##


## summarize the putative taxa to be estimated
microbes <- factor(colnames(y.binary.filtered.subgroup),levels=microbeTree.Y.root$tip.label)
microbeNames <- levels(microbes)
NMicrobeTips <- length(microbeNames)
##

## sort the OTU table so its entries match the tree's tip labels
y <- y.binary.filtered.subgroup[,microbeNames]
##

microbeTips <- colnames(y)
NMicrobeTips <- length(microbeTips)


senddat <- melt(y, varnames = c('sample','tip'), value.name = 'present')
sampleNames <- as.numeric(factor(senddat[,1]))
microbeTipNames <- as.numeric(factor(senddat[,2], levels = microbeTips))
present <- senddat[,3]
NObs <- nrow(senddat)
NSamples <- length(unique(sampleNames))

newermap <- newermap[levels(factor(senddat[,1])),]


NIntMicrobeNodes <- microbeTree.Y.root$Nnode
NMicrobeNodes <- NMicrobeTips + NIntMicrobeNodes - 1

microbeEdgeOrder <- order(microbeTree.Y.root$edge[,2])
microbeEdges <- microbeTree.Y.root$edge.length[microbeEdgeOrder]
NMicrobeTimeBins <- ceiling(NIntMicrobeNodes / NSplits)
microbeNHs <- nodeHeights(microbeTree.Y.root)
maxMicrobeNHs <- max(microbeNHs)
##

#divide total evolutionary time into chunks that contain approximately equal numbers of splits
lttMicrobeTree <- ltt(microbeTree.Y.root, log.lineages = F, plot = F)
temp <- maxMicrobeNHs - lttMicrobeTree$times[-length(lttMicrobeTree$times)]
microbeSplitTimes <- split(temp, ceiling(seq_along(temp) / NSplits))
if(length(microbeSplitTimes[[NMicrobeTimeBins]]) < NSplits/2) {
    microbeSplitTimes[[NMicrobeTimeBins - 1]] <- c(microbeSplitTimes[[NMicrobeTimeBins - 1]],
                                                   microbeSplitTimes[[NMicrobeTimeBins]])
    microbeSplitTimes[[NMicrobeTimeBins]] <- NULL
    NMicrobeTimeBins <- NMicrobeTimeBins - 1
}

#cut points in each phylogeny that would result in approximately equal numbers of splits per bin of time
microbeBoundaries <- sapply(2:NMicrobeTimeBins, function(x) max(microbeSplitTimes[[x]]))
microbeBoundariesRounded <- round(microbeBoundaries, 1)
 microbeTimeBinSizes <- sapply(2:(length(microbeBoundaries)+2), function(x) c(maxMicrobeNHs, microbeBoundaries, 0)[x-1] - c(maxMicrobeNHs, microbeBoundaries, 0)[x]) #size of each bin
relativeMicrobeTimeBinSizes <- microbeTimeBinSizes / sum(microbeTimeBinSizes) #proportion of total tree height belonging to each bin (differs for each replicate due to the varying sizes of the oldest bin)
ndMicrobe <- maxMicrobeNHs - microbeNHs
microbeEdgeToBin <- matrix(NA, ncol = NMicrobeTimeBins, nrow = nrow(microbeTree.Y.root$edge)) #create a matrix for each replicate that describes the proportion of each time bin that each branch exists for
for(j in 1:NMicrobeTimeBins) {
    if(j == 1) {
        allin <- which(ndMicrobe[,2] >= microbeBoundaries[j])
        allout <- which(ndMicrobe[,1] <= microbeBoundaries[j])
        cedge <- which((ndMicrobe[,1] > microbeBoundaries[j]) & (ndMicrobe[,2] < microbeBoundaries[j]))
        microbeEdgeToBin[cedge,j] <- ndMicrobe[cedge,1] - microbeBoundaries[j]
    } else if(j == NMicrobeTimeBins) {
        allin <- which(ndMicrobe[,1] <= microbeBoundaries[j-1])
        allout <- which(ndMicrobe[,2] >= microbeBoundaries[j-1])
        cedge <- which((ndMicrobe[,1] > microbeBoundaries[j-1]) & (ndMicrobe[,2] < microbeBoundaries[j-1]))
        microbeEdgeToBin[cedge,j] <- microbeBoundaries[j-1] - ndMicrobe[cedge,2]
    } else {
        allin <- which((ndMicrobe[,1] <= microbeBoundaries[j-1]) & (ndMicrobe[,2] >= microbeBoundaries[j]))
        allout <- which((ndMicrobe[,1] <= microbeBoundaries[j]) | (ndMicrobe[,2] >= microbeBoundaries[j-1]))
        cedge1 <- which((ndMicrobe[,1] <= microbeBoundaries[j-1]) & (ndMicrobe[,1] > microbeBoundaries[j]) & (ndMicrobe[,2] < microbeBoundaries[j]))
        microbeEdgeToBin[cedge1,j] <- ndMicrobe[cedge1,1] - microbeBoundaries[j]
        cedge2 <- which((ndMicrobe[,1] > microbeBoundaries[j-1]) & (ndMicrobe[,2] < microbeBoundaries[j-1]) & (ndMicrobe[,2] >= microbeBoundaries[j]))
        microbeEdgeToBin[cedge2,j] <- microbeBoundaries[j-1] - ndMicrobe[cedge2,2]
        cedge3 <- which((ndMicrobe[,1] > microbeBoundaries[j-1]) & (ndMicrobe[,2] < microbeBoundaries[j]))
        microbeEdgeToBin[cedge3,j] <- microbeBoundaries[j-1] - microbeBoundaries[j]
    }
    microbeEdgeToBin[allin,j] <- microbeTree.Y.root$edge.length[allin]
    microbeEdgeToBin[allout,j] <- 0
 }
microbeEdgeToBin <- microbeEdgeToBin / maxMicrobeNHs
rownames(microbeEdgeToBin) <- microbeTree.Y.root$edge[,2]
microbeEdgeToBin <- microbeEdgeToBin[microbeEdgeOrder,]

allnodes <- unique(microbeTree.Y.root$edge[,1])
NSamplesRaw <- nrow(y)
microbeAncestors <- matrix(0, NMicrobeNodes + 1, NMicrobeNodes + 1)
for(node in 1:(NMicrobeNodes + 1)) {
    microbeAncestors[, node] <- as.numeric(1:(NMicrobeNodes + 1) %in% c(Ancestors(microbeTree.Y.root, node), node))
}
##

colnames(microbeAncestors) <- rownames(microbeAncestors) <- paste0('i',1:(NMicrobeNodes + 1))
colnames(microbeAncestors)[1:NMicrobeTips] <- rownames(microbeAncestors)[1:NMicrobeTips] <- paste0('t',colnames(y))

microbeAncestors <- microbeAncestors[-(NMicrobeTips + 1), -(NMicrobeTips + 1)]
microbeNHRel <- microbeNHs / maxMicrobeNHs
rownames(microbeNHRel) <- microbeTree.Y.root$edge[,2]
microbeNHRel <- microbeNHRel[microbeEdgeOrder,]
##

newermap$sequencing_depth <- rowSums(y.old[,microbeTips])
newermap$log_sequencing_depth <- log(newermap$sequencing_depth)
newermap$log_sequencing_depth_scaled <- scale(newermap$log_sequencing_depth)



modelform <- ~ ocean + ocean_area + reef_name + concatenated_date + colony_name + tissue_compartment + log_sequencing_depth_scaled
allfactors <- attr(terms.formula(modelform), "term.labels")
NFactors <- length(allfactors)
allfactorder <- sapply(allfactors, function(x) sum(gregexpr(':', x, fixed = TRUE)[[1]] > 0))

modelMat <- model.matrix(modelform, model.frame(newermap,na.action = NULL))
modelMat[is.na(modelMat)] <- 0
sumconts <- names(attr(modelMat, "contrasts")[attr(modelMat, "contrasts") == 'contr.sum'])

## create matrix relating each 'effect' (categorical and numeric) to the 'factor' that it belongs to
factLevelMat <- sapply(1:length(allfactors), function (j) {
    as.numeric(attr(modelMat, 'assign')[-1] == j)
})
colnames(factLevelMat) <- c(allfactors)

modelMat <- modelMat[,-1]
NEffects <- ncol(modelMat)

## rename factors that have 'sum contrasts' because by default they get arbitrary names (careful with interpretation of interactions... probably better to add them as separate 'main effects' produced by concatenation)
for(j in sumconts) {
    searchTerms <- paste0('^', j, 1:(nlevels(newermap[,j]) - 1), '$')
    replacementTerms <- paste0(j, levels(newermap[,j])[-nlevels(newermap[,j])])
    for(k in 1:length(searchTerms)) {
        colnames(modelMat) <- sub(searchTerms[[k]], replacementTerms[[k]], colnames(modelMat))
    }
}
rownames(factLevelMat) <- colnames(modelMat)

NHostTimeBins <- ceiling(length(levels(newermap[,sampleTipKey])) / NSplits)

## extract all the possible species that 'fungid' could refer to
possibleFungidSpecs <- grep(paste0(paste(possibleFungidGenera,collapse='|'),'_'), attr(hostTree, "TipLabel"), value=T)

## identify unique Scleractinian species in the data, and replace spaces with underscores
study.species <- gsub(' ', '_', levels(newermap[,sampleTipKey]))

## identify the Scleractinian species in the data that do not exist in the template tree
study.species.missing <- study.species[!study.species %in% grep(paste(c(attr(hostTree, "TipLabel"),'Fungid','not_applicable'),collapse='|'), study.species, ignore.case=T, value=T)]
generaOfUnknowns <- sapply(study.species.missing, function(x) strsplit(x, '_')[[1]][[1]])


### starting here, generate multiple random samples of the map and trees (later to be summarized to define the time Bins and also to be separately used for replicate runs of the model)
sampleMap <- list()
hostTreesSampled <- list()
hostSplitTimes <- list()
hostBoundaries <- matrix(NA, nrow = NTrees, ncol = NHostTimeBins - 1)
for(i in 1:NTrees) {
    sampleMap[[i]] <- newermap
    fungidSps <- grep('Fungid',levels(sampleMap[[i]][,sampleTipKey]))
    levels(sampleMap[[i]][,sampleTipKey])[fungidSps] <- sample(possibleFungidSpecs[!possibleFungidSpecs %in% levels(sampleMap[[i]][,sampleTipKey])], length(fungidSps)) ##assign unidentified Fungids to a random member of the group independently for each tree
    for (j in unique(generaOfUnknowns)) {
        possibleGenera <- attr(hostTree, "TipLabel")[!attr(hostTree, "TipLabel") %in% levels(sampleMap[[i]][,sampleTipKey])]
        if(!any(grepl(j,possibleGenera))) {
            possibleGenera <- attr(hostTree, "TipLabel")
        }
        levels(sampleMap[[i]][,sampleTipKey])[levels(sampleMap[[i]][,sampleTipKey]) %in% grep(paste0(j,'_'),study.species.missing,value=T)] <- sample(grep(paste0(j,'_'), possibleGenera, value=T), sum(generaOfUnknowns == j)) ## assign unidentified species to a random member of their genus independently for each tree
    }

    #filter the tree only contain the sampled (or assigned) species
    hostTreesSampled[[i]] <- ladderize(drop.tip(hostTree[[i]],hostTree[[i]]$tip.label[!hostTree[[i]]$tip.label %in% levels(sampleMap[[i]][,sampleTipKey])]))
    
    #divide total evolutionary time into chunks that contain approximately equal numbers of splits
    lttHostTreesSampled <- ltt(hostTreesSampled[[i]],log.lineages=F,plot=F)
    temp <- max(nodeHeights(hostTreesSampled[[i]])) - lttHostTreesSampled$times[-length(lttHostTreesSampled$times)]
    hostSplitTimes[[i]] <- split(temp, ceiling(seq_along(temp)/NSplits))
    if(length(hostSplitTimes[[i]][[NHostTimeBins]]) < NSplits/2) {
        hostSplitTimes[[i]][[NHostTimeBins - 1]] <- c(hostSplitTimes[[i]][[NHostTimeBins - 1]], hostSplitTimes[[i]][[NHostTimeBins]])
        hostSplitTimes[[i]][[NHostTimeBins]] <- NULL
        NHostTimeBins <- NHostTimeBins - 1
        hostBoundaries <- hostBoundaries <- matrix(NA,nrow=NTrees, ncol=NHostTimeBins-1)
    }
    
    #cut points in each phylogeny that would result in approximately equal numbers of splits per bin of time
    hostBoundaries[i,] <- sapply(2:NHostTimeBins, function(x) max(hostSplitTimes[[i]][[x]]))
    
}


#average cut points among all the trees so the bins of time are consistent across replicates
meanHostBoundaries <- apply(hostBoundaries,2,mean)
meanHostBoundariesRounded <- round(meanHostBoundaries,1)
##

NHostTips <- length(hostTreesSampled[[1]]$tip.label)
NIntHostNodes <- hostTreesSampled[[1]]$Nnode
NHostNodes <- NIntHostNodes + NHostTips - 1


hostTimeBinSizes <- list()
relativeHostTimeBinSizes <- list()
hostEdgeToBin <- list()
hostEdgeOrder <- list()
nh <- list()
maxNH <- list()
nhRel <- list()
for(i in 1:NTrees) {
    nh[[i]] <- nodeHeights(hostTreesSampled[[i]])
    maxNH[[i]] <- max(nh[[i]])
    
    #size in millions of years of each bin (consistent across replicates /except/ for the oldest bin
    hostTimeBinSizes[[i]] <- sapply(2:(length(meanHostBoundaries) + 2),
                                function(x) {c(maxNH[[i]], meanHostBoundaries, 0)[x - 1]
                                             - c(maxNH[[i]], meanHostBoundaries, 0)[x]})
    relativeHostTimeBinSizes[[i]] <- hostTimeBinSizes[[i]] / sum(hostTimeBinSizes[[i]]) #proportion of total tree height belonging to each bin (differs for each replicate due to the varying sizes of the oldest bin)
    nd <- maxNH[[i]] - nh[[i]]
    hostEdgeToBin[[i]] <- matrix(NA,ncol=NHostTimeBins,nrow=nrow(hostTreesSampled[[i]]$edge)) #create a matrix for each replicate that describes the proportion of each time bin that each branch exists for
    for(j in 1:NHostTimeBins) {
        if(j == 1) {
            allin <- which(nd[,2] >= meanHostBoundaries[j])
            allout <- which(nd[,1] <= meanHostBoundaries[j])
            cedge <- which((nd[,1] > meanHostBoundaries[j]) & (nd[,2] < meanHostBoundaries[j]))
            hostEdgeToBin[[i]][cedge,j] <- nd[cedge,1] - meanHostBoundaries[j]
        } else if(j == NHostTimeBins) {
            allin <- which(nd[,1] <= meanHostBoundaries[j-1])
            allout <- which(nd[,2] >= meanHostBoundaries[j-1])
            cedge <- which((nd[,1] > meanHostBoundaries[j-1]) & (nd[,2] < meanHostBoundaries[j-1]))
            hostEdgeToBin[[i]][cedge,j] <- meanHostBoundaries[j-1] - nd[cedge,2]
        } else {
            allin <- which((nd[,1] <= meanHostBoundaries[j-1]) & (nd[,2] >= meanHostBoundaries[j]))
            allout <- which((nd[,1] <= meanHostBoundaries[j]) | (nd[,2] >= meanHostBoundaries[j-1]))
            cedge1 <- which((nd[,1] <= meanHostBoundaries[j-1]) & (nd[,1] > meanHostBoundaries[j]) & (nd[,2] < meanHostBoundaries[j]))
            hostEdgeToBin[[i]][cedge1,j] <- nd[cedge1,1] - meanHostBoundaries[j]
            cedge2 <- which((nd[,1] > meanHostBoundaries[j-1]) & (nd[,2] < meanHostBoundaries[j-1]) & (nd[,2] >= meanHostBoundaries[j]))
            hostEdgeToBin[[i]][cedge2,j] <- meanHostBoundaries[j-1] - nd[cedge2,2]
            cedge3 <- which((nd[,1] > meanHostBoundaries[j-1]) & (nd[,2] < meanHostBoundaries[j]))
            hostEdgeToBin[[i]][cedge3,j] <- meanHostBoundaries[j-1] - meanHostBoundaries[j]
        }
        hostEdgeToBin[[i]][allin,j] <- hostTreesSampled[[i]]$edge.length[allin]
        hostEdgeToBin[[i]][allout,j] <- 0
        
    }
    hostEdgeToBin[[i]] <- hostEdgeToBin[[i]] / maxNH[[i]]
    
    hostEdgeOrder[[i]] <- order(hostTreesSampled[[i]]$edge[,2])
    rownames(hostEdgeToBin[[i]]) <- hostTreesSampled[[i]]$edge[,2]
    hostEdgeToBin[[i]] <- hostEdgeToBin[[i]][hostEdgeOrder[[i]],]
    
    nhRel[[i]] <- nh[[i]] / maxNH[[i]]
    rownames(nhRel[[i]]) <- hostTreesSampled[[i]]$edge[,2]
    nhRel[[i]] <- nhRel[[i]][hostEdgeOrder[[i]],]
}
##

## create ancestry matrices for each host tree
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
                         microbeNodeHeights             = microbeNHRel,
                         microbeEdgeToBin               = t(microbeEdgeToBin),
                         hostEdgeToBin                  = hostEdgeToBin[[i]],
                         NHostNodes                     = NHostNodes,
                         NHostTimeBins                  = NHostTimeBins,
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
dir.create(outdir, recursive = T)
save.image(file.path(outdir,'setup.RData'))
save(microbeAncestors, file = file.path(outdir, 'microbeAncestors.RData'))
save(meanHostBoundaries, file = file.path(outdir, 'meanHostBoundaries.RData'))
save(hostAncestors, file = file.path(outdir, 'hostAncestors.RData'))
save(hostAncestorsExpanded, file = file.path(outdir, 'hostAncestorsExpanded.RData'))
save(hostTreesSampled, file = file.path(outdir, 'hostTreesSampled.RData'))
save(hostTimeBinSizes, file = file.path(outdir, 'hostTimeBinSizes.RData'))
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

fitModes <- sapply(fit, function(x) x@mode)

## create parental ancestry matrix for microbes (only sums effect of a given node and its direct parent, not all ancestors
microbeParents <- matrix(0, NMicrobeNodes + 1, NMicrobeNodes + 1)
for(node in 1:(NMicrobeNodes + 1)) {
    microbeParents[node, ] <- as.numeric(1:(NMicrobeNodes + 1) %in% c(Ancestors(microbeTree.Y.root, node, 'parent'), node))
}
colnames(microbeParents) <- rownames(microbeParents) <- paste0('i',1:(NMicrobeNodes + 1))
colnames(microbeParents)[1:NMicrobeTips] <- rownames(microbeParents)[1:NMicrobeTips] <- paste0('t', colnames(y))
microbeParentsT <- t(cbind(1, microbeParents[-(NMicrobeTips + 1), -(NMicrobeTips + 1)]))

hostParents <- list()
microbeAncestorsT <- t(cbind(1, microbeAncestors))

colorpal <- colorRampPalette(brewer.pal(9,'Blues'))
plotcolors <- c('white', colorpal(100),'black')

## summarize the results separately for each sampled host tree
for(i in 1:NTrees) {
    
    if (fitModes[[i]] != 0) next
    
    check_hmc_diagnostics(fit[[i]])

    currplotdir <- file.path(outdir, paste0('tree_',i), 'plots')
    currtabledir <- file.path(outdir, paste0('tree_',i), 'tables')
    currdatadir <- file.path(outdir, paste0('tree_',i), 'data')

    dir.create(currplotdir, recursive = T)
    dir.create(currtabledir, recursive = T)
    dir.create(currdatadir, recursive = T)

    ## plot the sampled tree with the time bins marked
    pdf(file = file.path(currplotdir,'sampledHostTree.pdf'), width = 25, height = 15)
    plot(hostTreesSampled[[i]],cex = 0.75)
    for(age in max(nodeHeights(hostTreesSampled[[i]])) - meanHostBoundaries) {
        lines(x = c(age, age), y = c(1, length(hostTreesSampled[[i]]$tip.label)), lwd = 1)
    }
    graphics.off()
    ##
    
    ## plot the sampled tree with the time bins marked
    pdf(file = file.path(currplotdir,'sampledMicrobeTree.pdf'), width = 25, height = 15)
    plot(microbeTree.Y.root,cex = 0.75)
    for(age in max(nodeHeights(microbeTree.Y.root)) - microbeBoundaries) {
        lines(x = c(age, age), y = c(1, length(microbeTree.Y.root$tip.label)), lwd = 1)
    }
    graphics.off()
    ##

    ## variance partitioning
    stDProps <- extract(fit[[i]], pars = 'stDProps')[[1]]
    colnames(stDProps) <- c(paste0('ADiv.',colnames(factLevelMat)),
                            paste0('Specificity.', colnames(factLevelMat)),
                            'ADiv.host',
                            'host.specificity',
                            'microbe.prevalence')

    pdf(file = file.path(currplotdir,'scalesboxes.pdf'), width = 25, height = 15)
    boxplot(stDProps, cex.axis = 0.5, las = 2)
    graphics.off()

    save(stDProps, file = file.path(currdatadir, 'stDProps.RData'))
    ##
    
    ## proportion of variance explained by each time bin
    metaVarProps <- extract(fit[[i]], pars = 'metaVarProps')[[1]]
    colnames(metaVarProps) <- c('prevalence', 'adiv', 'specificty')
    pdf(file = file.path(currplotdir, 'metaVarProps_boxes.pdf'), width = 25, height = 15)
    boxplot(metaVarProps, cex.axis = 0.5, las = 2)
    graphics.off()
    save(metaVarProps,file = file.path(currdatadir, 'metaVarProps.RData'))
    ##
    
    ## compare rates of host evolution in each time bin
    relativeHostEvolRates <- exp(extract(fit[[i]], pars = 'logRelativeHostEvolRates')[[1]])
    colnames(relativeHostEvolRates) <- sapply(1:(length(meanHostBoundariesRounded) + 1),
                                              function(x) paste0(c('before', meanHostBoundariesRounded)[[x]],
                                                                 ' - ',
                                                                 c(meanHostBoundariesRounded, 'present')[[x]],
                                                                 ' mya'))

    pdf(file = file.path(currplotdir, 'evolRatesHost.pdf'), width = 25, height = 15)
    boxplot(relativeHostEvolRates, xlab = 'Time Period', ylab = 'Rate of Evolution Relative to Mean')
    graphics.off()
    
    pdf(file = file.path(currplotdir, 'logEvolRatesHost.pdf'), width = 25, height = 15)
    boxplot(log(relativeHostEvolRates), xlab = 'Time Period', ylab = 'Log Rate of Evolution Relative to Mean')
    graphics.off()

	save(relativeHostEvolRates, file = file.path(currdatadir, 'relativeHostEvolRates.RData'))
    ##
    
    ## compare rates of microbe evolution in each time bin
    relativeMicrobeEvolRates <- exp(extract(fit[[i]], pars = 'logRelativeMicrobeEvolRates')[[1]])
    
    pdf(file = file.path(currplotdir, 'evolRatesMicrobe.pdf'), width = 25, height = 15)
    boxplot(relativeMicrobeEvolRates, xlab = 'Time Period', ylab = 'Rate of Evolution Relative to Mean')
    graphics.off()
    
    pdf(file = file.path(currplotdir, 'logEvolRatesMicrobe.pdf'), width = 25, height = 15)
    boxplot(log(relativeMicrobeEvolRates), xlab = 'Time Period', ylab = 'Log Rate of Evolution Relative to Mean')
    graphics.off()
    
    save(relativeMicrobeEvolRates, file = file.path(currdatadir, 'relativeMicrobeEvolRates.RData'))
    ##
    
    ## ornstein-uhlenbeck parameters
    hostOUAlpha <- extract(fit[[i]], pars = 'hostOUAlpha')[[1]]
    microbeOUAlpha <- extract(fit[[i]], pars = 'microbeOUAlpha')[[1]]
    OUAlphas <- cbind(hostOUAlpha, microbeOUAlpha)
    colnames(OUAlphas) <- c('host', 'microbe')
    
    pdf(file = file.path(currplotdir, 'OUAlphas.pdf'), width = 25, height = 15)
    boxplot(OUAlphas, xlab = 'Host or microbe', ylab = 'alpha')
    graphics.off()
    
    save(OUAlphas, file = file.path(currdatadir, 'OUAlphas.RData'))
    ##
    
    ## extract effects from model
    scaledMicrobeNodeEffects <- array(extract(fit[[i]],
                                              pars       = 'scaledMicrobeNodeEffects',
                                              permuted   = F,
                                              inc_warmup = T),
                                      dim = c(NMCSamples,
                                              NChains,
                                              NEffects + NHostNodes + 1,
                                              NMicrobeNodes + 1),
                                      dimnames = list(sample  = NULL,
                                                      chain   = NULL,
                                                      effect  = c('microbePrevalence',
                                                                  colnames(modelMat)[1:NEffects],
                                                                  paste0('host', colnames(hostAncestors[[i]]))),
                                                      taxnode = c('alphaDiversity', colnames(microbeAncestors))))
                                                    
    save(scaledMicrobeNodeEffects, file = file.path(currdatadir, 'scaledMicrobeNodeEffects.RData'))
    ##
    
    ## calculate base-level effects (negative sum of all others in category)
    baseLevelEffects <- array(NA,
                              dim = c(NMCSamples,
                                      NChains,
                                      length(sumconts),
                                      NMicrobeNodes + 1),
                              dimnames = list(sample  = NULL,
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
    ##

    ## summarize posterior distributions of effects
    allRes <- NULL
    for(j in 1:(NEffects + NHostNodes + 1)) {
        temp <- monitor(array(scaledMicrobeNodeEffects[,,j,],
                              dim = c(NMCSamples, NChains, NMicrobeNodes + 1)),
                        warmup = warmup,
                        probs = c(0.05, 0.95))
        temp <- cbind(effect = dimnames(scaledMicrobeNodeEffects)[[3]][j], temp)
        rownames(temp) <- c('alphaDiversity', rownames(microbeAncestors))
        allRes <- rbind(allRes, temp)
    }
    ##
    
    ## summarize posterior distibutions of base-level effects
    for(m in sumconts) {
        temp <- monitor(array(baseLevelEffects[,,m,],
                              dim = c(NMCSamples, NChains, NMicrobeNodes + 1)),
                        warmup = warmup,
                        probs = c(0.05, 0.95))
        temp <- cbind(effect = paste0(m, levels(newermap[,m])[nlevels(newermap[,m])]), temp)
        rownames(temp) <- c('alphaDiversity', rownames(microbeAncestors))
        allRes <- rbind(allRes, temp)
    }
    ##
    
    ## write posterior summaries of effects to file
    cat('microbeNode\t', file = file.path(currtabledir, 'nodeEffects.txt'))
    write.table(allRes,
                file   = file.path(currtabledir, 'nodeEffects.txt'),
                sep    = '\t',
                quote  = F,
                append = T)
    ##
    
    ## extract variance modifier terms from the fit model
    phyloLogVarMultPrev <- array(extract(fit[[i]],
                                         pars = 'phyloLogVarMultPrev',
                                         permuted = F,
                                         inc_warmup = T),
                                 dim = c(NMCSamples,
                                         NChains,
                                         NMicrobeNodes),
                                 dimnames = list(sample = NULL,
                                                 chain = NULL,
                                                 taxnode = colnames(microbeAncestors)))
    save(phyloLogVarMultPrev, file = file.path(currdatadir, 'phyloLogVarMultPrev.RData'))
    
    phyloLogVarMultADiv <- array(extract(fit[[i]],
                                         pars = 'phyloLogVarMultADiv',
                                         permuted = F,
                                         inc_warmup = T),
                                 dim = c(NMCSamples,
                                         NChains,
                                         NHostNodes),
                                 dimnames = list(sample = NULL,
                                                 chain = NULL,
                                                 taxnode = colnames(hostAncestors[[i]])))
    save(phyloLogVarMultADiv, file = file.path(currdatadir, 'phyloLogVarMultADiv.RData'))
    
    phyloLogVarMultRaw <- array(extract(fit[[i]],
                                        pars = 'phyloLogVarMultRaw',
                                        permuted = F,
                                        inc_warmup = T),
                                dim = c(NMCSamples,
                                        NChains,
                                        NHostNodes,
                                        NMicrobeNodes),
                                dimnames = list(sample = NULL,
                                                chain = NULL,
                                                hostnode = colnames(hostAncestors[[i]]),
                                                microbenode = colnames(microbeAncestors)))
    save(phyloLogVarMultRaw, file = file.path(currdatadir, 'phyloLogVarMultRaw.RData'))
    
    metaScales <- array(extract(fit[[i]],
                                pars = 'metaScales',
                                permuted = F,
                                inc_warmup = T),
                        dim = c(NMCSamples,
                                NChains,
                                3))
    save(metaScales, file = file.path(currdatadir, 'metaScales.RData'))
    
    matMult <- array(NA,
                     dim = c(NMCSamples,
                             NChains,
                             NHostNodes,
                             NMicrobeNodes),
                     dimnames = list(sample = NULL,
                                     chain = NULL,
                                     hostnode = colnames(hostAncestors[[i]]),
                                     microbenode = colnames(microbeAncestors)))
    ##
    
    ## see if any clades have higher variance among their descendants (maybe suggesting codiversification)
    allRes <- NULL
    for(j in 1:NHostNodes) {
        temp <- monitor(array(phyloLogVarMultRaw[,,j,],
                              dim = c(NMCSamples, NChains, NMicrobeNodes)),
                        warmup = warmup,
                        probs = c(0.05, 0.95))
        temp <- cbind(hostNode = colnames(hostAncestors[[i]])[[j]], temp)
        rownames(temp) <- rownames(microbeAncestors)
        allRes <- rbind(allRes, temp)
    }
    
    cat('microbeNode\t', file = file.path(currtabledir, 'individualPhyloVarianceEffects.txt'))
    write.table(allRes,
                file   = file.path(currtabledir, 'individualPhyloVarianceEffects.txt'),
                sep    = '\t',
                quote  = F,
                append = T)
    ##
    
    ## since the model has many degrees of freedom, it's possible that certain regions of the tree have 'significant' effects but that those effects are 'implemented' via different, nearby nodes in each iteration. We can sum the effects of all ancestral nodes to see if a given node is consistently influenced by all of itself plus its ancestral terms, even if all of those terms are individually inconsistent
    for(j in 1:NMCSamples) {
        for(k in 1:NChains) {
            matMult[j,k,,] <- cbind(1, hostAncestors[[i]]) %*% rbind(c(0, phyloLogVarMultPrev[j,k,] * metaScales[j,k,1]), cbind(phyloLogVarMultADiv[j,k,] * metaScales[j,k,2], phyloLogVarMultRaw[j,k,,] * metaScales[j,k,3])) %*% microbeAncestorsT
        }
    }
    
    allRes <- NULL
    for(j in 1:NHostNodes) {
        temp <- monitor(array(matMult[,,j,],
                              dim = c(NMCSamples, NChains, NMicrobeNodes)),
                        warmup = warmup,
                        probs = c(0.05, 0.95))
        temp <- cbind(hostNode = colnames(hostAncestors[[i]])[[j]], temp)
        rownames(temp) <- rownames(microbeAncestors)
        allRes <- rbind(allRes, temp)
    }
    
    cat('microbeNode\t', file = file.path(currtabledir, 'allSummedPhyloVarianceEffects.txt'))
    write.table(allRes,
                file   = file.path(currtabledir, 'allSummedPhyloVarianceEffects.txt'),
                sep    = '\t',
                quote  = F,
                append = T)
    ##
    
    ## we can also sum the effects of just a node and its parent, as a bit of a middle-ground approach to see whether a given node is significantly different from its grandparent instead of its parent or instead of the overall mean
    hostParents[[i]] <- matrix(NA, NHostNodes+1, NHostNodes+1)
    for(node in 1:(NHostNodes+1)) {
        hostParents[[i]][node, ] <- as.numeric(1:(NHostNodes+1) %in% c(Ancestors(hostTreesSampled[[i]], node, 'parent'), node))
    }
    colnames(hostParents[[i]]) <- rownames(hostParents[[i]]) <- paste0('i',1:(NHostNodes+1))
    colnames(hostParents[[i]])[1:NHostTips] <- rownames(hostParents[[i]])[1:NHostTips] <- hostTreesSampled[[i]]$tip.label
    hostParents[[i]] <- cbind(1, hostParents[[i]][-(NHostTips + 1), -(NHostTips + 1)])
    
    for(j in 1:NMCSamples) {
        for(k in 1:NChains) {
            matMult[j,k,,] <- hostParents[[i]] %*% rbind(c(0, phyloLogVarMultPrev[j,k,] * metaScales[j,k,1]), cbind(phyloLogVarMultADiv[j,k,] * metaScales[j,k,2], phyloLogVarMultRaw[j,k,,] * metaScales[j,k,3])) %*% microbeParentsT
        }
    }
    
    allRes <- NULL
    for(j in 1:NHostNodes) {
        temp <- monitor(array(matMult[,,j,],
                              dim = c(NMCSamples, NChains, NMicrobeNodes)),
                        warmup = warmup,
                        probs = c(0.05, 0.95))
        temp <- cbind(hostNode = colnames(hostAncestors[[i]])[[j]], temp)
        rownames(temp) <- rownames(microbeAncestors)
        allRes <- rbind(allRes, temp)
    }
    
    cat('microbeNode\t', file = file.path(currtabledir, 'parentalSummedPhyloVarianceEffects.txt'))
    write.table(allRes,
                file   = file.path(currtabledir, 'parentalSummedPhyloVarianceEffects.txt'),
                sep    = '\t',
                quote  = F,
                append = T)
    ##
    
    ## see if any host clades have higher variance among their descendants
    sums <- summary(fit[[i]], pars = 'phyloLogVarMultADiv', probs = c(0.05,0.95), use_cache = F)
    rownames(sums$summary) <- colnames(hostAncestors[[i]])
    
    cat('hostNode\t', file = file.path(currtabledir, 'phylogeneticADivEffects.txt'))
    write.table(sums$summary,
                file   = file.path(currtabledir, 'phylogeneticADivEffects.txt'),
                sep    = '\t',
                quote  = F,
                append = T)
    ##
    
    ## see if any microbe clades have higher variance among their descendants
    sums <- summary(fit[[i]], pars = 'phyloLogVarMultPrev', probs = c(0.05,0.95), use_cache = F)
    rownames(sums$summary) <- rownames(microbeAncestors)
    
    cat('microbeNode\t', file = file.path(currtabledir, 'phylogeneticPrevalenceEffects.txt'))
    write.table(sums$summary,
                file   = file.path(currtabledir, 'phylogeneticPrevalenceEffects.txt'),
                sep    = '\t',
                quote  = F,
                append = T)
    ##
    
    ## summarize the mean branch lengths of the microbes
    sums <- summary(fit[[i]], pars = 'microbeScales', probs = c(0.05,0.95), use_cache = F)
    newEdges <- sums$summary[,'mean']^2
    microbeTree.Y.root.newEdges <- microbeTree.Y.root
    microbeTree.Y.root.newEdges$edge.length <- newEdges[order(microbeEdgeOrder)]
    pdf(file = file.path(currplotdir, 'microbeTreeWEstimatedEdgeLengths.pdf'), width = 25, height = 15)
    plot(microbeTree.Y.root.newEdges, cex = 0.5)
    graphics.off()
    ##

    ## summarize the mean branch lengths of the hosts
    sums <- summary(fit[[i]], pars = 'hostScales', probs = c(0.05,0.95), use_cache = F)
    newEdges <- sums$summary[,'mean']^2
    hostTreesSampled.newEdges <- hostTreesSampled[[i]]
    hostTreesSampled.newEdges$edge.length <- newEdges[order(hostEdgeOrder[[i]])]
    pdf(file = file.path(currplotdir,'hostTreeWEstimatedEdgeLengths.pdf'), width = 25, height = 15)
    plot(hostTreesSampled.newEdges, cex = 0.5)
    graphics.off()
    ##
    
    ## plot heatmap of cophylogenetic patterns
    plotmicrobetree <- ladderize(multi2di(microbeTree.Y.root))
    hclmicrobetree <- as.hclust(plotmicrobetree)
    # for each sample in the data, assign its mitotype to a vector
    hostvect <- sampleMap[[i]]$host_scientific_name
    # name the vector with the sample IDs
    names(hostvect) <- rownames(sampleMap[[i]])
    # sort the samples in the vector
    temp <- sampleMap[[i]][order(sampleMap[[i]]$concatenated_date),]
    temp <- temp[order(temp$reef_name),]
    temp <- temp[order(temp$host_scientific_name),]
    for(comp in c('T', 'S', 'M')) {
        temp2 <- temp[temp$tissue_compartment == comp,]
        hostvectTemp <- hostvect[rownames(temp2)]
        # expand the tips into polytomies containing a tip for each sample
        hosttree <- expandTaxonTree(hostTreesSampled[[i]], hostvectTemp, keepBrLen = T)
        hosttree <- drop.tip(hosttree, hosttree$tip.label[!hosttree$tip.label %in% names(hostvectTemp)])
        # convert polytomies into randomly-split, binary subtrees with 0-length branches, and ladderize the whole tree
        hosttree.dichotomous <- ladderize(multi2di(hosttree, random = F), right = F)
        # convert the phylogenetic tree into an hclust object
        hclhosttree <- as.hclust(force.ultrametric(hosttree.dichotomous))
        plotFilt <- as.matrix(t(y)[plotmicrobetree$tip.label, hosttree.dichotomous$tip.label])
        pdf(file   = file.path(currplotdir,
                               paste0('cophylogeny_heatmap_',
                                      comp,
                                      '.pdf')),
            width  = 10,
            height = 10)
        heatmap(plotFilt,
                Rowv   = as.dendrogram(hclmicrobetree),
                Colv   = as.dendrogram(hclhosttree),
                col    = plotcolors,
                cexCol = 0.2,
                cexRow = 0.1,
                scale  = 'none')
        graphics.off()
    }
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

NSuccessTrees <- sum(fitModes == 0)

if (NSuccessTrees > 1) {
    
    allfit <- sflist2stanfit(fit[fitModes == 0])

    stDProps <- extract(allfit, pars = 'stDProps')[[1]]
    colnames(stDProps) <- c(paste0('ADiv.', colnames(factLevelMat)),
                            paste0('Specificity.', colnames(factLevelMat)),
                            'ADiv.host',
                            'host.specificity',
                            'microbe.prevalence')

    pdf(file = file.path(currplotdir, 'scalesboxes.pdf'), width = 25, height = 15)
    boxplot(stDProps, cex.axis = 0.5, las = 2)
    graphics.off()

    save(stDProps, file = file.path(currdatadir, 'stDProps.RData'))

    relativeHostEvolRates <- exp(extract(allfit, pars='logRelativeHostEvolRates')[[1]])
    colnames(relativeHostEvolRates) <- sapply(1:(length(meanHostBoundariesRounded) + 1),
                                              function(x) paste0(c('before', meanHostBoundariesRounded)[[x]],
                                                                 ' - ',
                                                                 c(meanHostBoundariesRounded, 'present')[[x]],
                                                                 ' mya'))

    pdf(file = file.path(currplotdir, 'evolRatesHost.pdf'), width = 25, height = 15)
    boxplot(relativeHostEvolRates, xlab = 'Time Period', ylab = 'Rate of Evolution Relative to Mean')
    graphics.off()

    pdf(file = file.path(currplotdir, 'logEvolRatesHost.pdf'), width = 25, height = 15)
    boxplot(log(relativeHostEvolRates), xlab = 'Time Period', ylab = 'Log Rate of Evolution Relative to Mean')
    graphics.off()

    save(relativeHostEvolRates, file = file.path(currdatadir, 'relativeHostEvolRates.RData'))

    ## compare rates of microbe evolution in each time bin
    relativeMicrobeEvolRates <- exp(extract(allfit, pars = 'logRelativeMicrobeEvolRates')[[1]])

    pdf(file = file.path(currplotdir, 'evolRatesMicrobe.pdf'), width = 25, height = 15)
    boxplot(relativeMicrobeEvolRates, xlab = 'Time Period', ylab = 'Rate of Evolution Relative to Mean')
    graphics.off()

    pdf(file = file.path(currplotdir, 'logEvolRatesMicrobe.pdf'), width = 25, height = 15)
    boxplot(log(relativeMicrobeEvolRates), xlab = 'Time Period', ylab = 'Log Rate of Evolution Relative to Mean')
    graphics.off()

    save(relativeMicrobeEvolRates, file = file.path(currdatadir, 'relativeMicrobeEvolRates.RData'))
    ##

    ## ornstein-uhlenbeck parameters
    hostOUAlpha <- extract(allfit, pars = 'hostOUAlpha')[[1]]
    microbeOUAlpha <- extract(allfit, pars = 'microbeOUAlpha')[[1]]
    OUAlphas <- cbind(hostOUAlpha, microbeOUAlpha)
    colnames(OUAlphas) <- c('host', 'microbe')

    pdf(file = file.path(currplotdir,'OUAlphas.pdf'), width = 25, height = 15)
    boxplot(OUAlphas, xlab = 'Host or microbe', ylab = 'alpha')
    graphics.off()

    save(OUAlphas, file = file.path(currdatadir, 'OUAlphas.RData'))
    ##

    ## summarize the mean branch lengths of the microbes
    sums <- summary(allfit, pars = 'microbeScales', probs = c(0.05,0.95), use_cache = F)
    newEdges <- sums$summary[,'mean']^2
    microbeTree.Y.root.newEdges <- microbeTree.Y.root
    microbeTree.Y.root.newEdges$edge.length <- newEdges[order(microbeEdgeOrder)]
    pdf(file = file.path(currplotdir, 'microbeTreeWEstimatedEdgeLengths.pdf'), width = 25, height = 15)
    plot(microbeTree.Y.root.newEdges, cex = 0.5)
    graphics.off()
    ##

    ## extract effects from model
    scaledMicrobeNodeEffects <- array(extract(allfit, pars='scaledMicrobeNodeEffects', permuted=F, inc_warmup=T),
                                      dim = c(NMCSamples,
                                              NChains * NTrees,
                                              NEffects + NHostNodes + 1,
                                              NMicrobeNodes + 1),
                                      dimnames = list(sample  = NULL,
                                                      chain   = NULL,
                                                      effect  = c('microbePrevalence',
                                                                  colnames(modelMat)[1:NEffects],
                                                                  paste0('host', colnames(hostAncestors[[i]]))),
                                                      taxnode = c('alphaDiversity', colnames(microbeAncestors))))
                                                    
    save(scaledMicrobeNodeEffects, file = file.path(currdatadir, 'scaledMicrobeNodeEffects.RData'))
    ##

    ## calculate base-level effects (negative sum of all others in category)
    baseLevelEffects <- array(NA,
                              dim = c(NMCSamples,
                                      NChains * NTrees,
                                      length(sumconts),
                                      NMicrobeNodes + 1),
                              dimnames = list(sample  = NULL,
                                              chain   = NULL,
                                              effect  = sumconts,
                                              taxnode = c('alphaDiversity', colnames(microbeAncestors))))
    for(j in 1:NMCSamples) {
        for(k in 1:(NChains * NTrees)) {
            for(m in sumconts) {
                baseLevelEffects[j,k,m,] <- -colSums(scaledMicrobeNodeEffects[j, k, rownames(factLevelMat)[factLevelMat[,m] == 1],])
            }
        }
    }

    save(baseLevelEffects, file = file.path(currdatadir, 'baseLevelEffects.RData'))
    ##

    ## summarize posterior distributions of effects
    allRes <- NULL
    for(j in 1:(NEffects + NHostNodes + 1)) {
        temp <- monitor(array(scaledMicrobeNodeEffects[,,j,],
                              dim = c(NMCSamples, NChains * NTrees, NMicrobeNodes + 1)),
                        warmup = warmup,
                        probs  = c(0.05, 0.95))
        temp <- cbind(effect = dimnames(scaledMicrobeNodeEffects)[[3]][j], temp)
        rownames(temp) <- c('alphaDiversity', rownames(microbeAncestors))
        allRes <- rbind(allRes, temp)
    }
    ##

    ## summarize posterior distibutions of base-level effects
    for(m in sumconts) {
        temp <- monitor(array(baseLevelEffects[,,m,],
                              dim = c(NMCSamples, NChains * NTrees, NMicrobeNodes + 1)),
                        warmup = warmup,
                        probs  = c(0.05, 0.95))
        temp <- cbind(effect = paste0(m, levels(newermap[,m])[nlevels(newermap[,m])]), temp)
        rownames(temp) <- c('alphaDiversity', rownames(microbeAncestors))
        allRes <- rbind(allRes, temp)
    }
    ##

    ## write posterior summaries of effects to file
    cat('microbeNode\t', file = file.path(currtabledir, 'nodeEffects.txt'))
    write.table(allRes,
                file   = file.path(currtabledir, 'nodeEffects.txt'),
                sep    = '\t',
                quote  = F,
                append = T)
    ##
}

## fin
