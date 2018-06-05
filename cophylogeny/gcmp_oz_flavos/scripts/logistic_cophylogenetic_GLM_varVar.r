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


microbeTreePath <- 'raw_data/gg_constrained_fastttree.tre' #ML microbial phylogeny
hosttreePath <- 'raw_data/combined_trees.newick' #set of Bayesian draws of host species phylogeny
mapfilePath <- 'raw_data/GCMP_EMP_map_r28.txt' #mapping file
fulltablePath <- 'raw_data/reference-hit.txt' #250 bp deblur otu table output
taxAssignmentPath <- 'raw_data/reference-hit.seqs_tax_assignments.txt' #greegenes taxonomy
modelPath <- 'scripts/logistic_cophylogenetic_GLM_varVar.stan' #stan model


outdir <- file.path('output',gsub(' ', '_', Sys.time()))

## filtration options
minCountSamp <- 5 # minimum sequencing depth for a sample to be included
minPercent <- 0 # minimum percent of a sample composed of a sequence variant for it to pass the below filter
minSamps <- 1 # minimum number of samples that a sequence variant is present in at the above threshold for it to be included
##


## Stan options
NChains <- 1 ## this is per tree; since I'm doing a large number of trees in parallel i'll just do one chain for each
NIterations <- 500 ## will probably need >10,000? maybe start with 2, check convergence, double it, check, double, check, double, etc.?
max_treedepth <- 15 ## a warning will tell you if this needs to be increased
adapt_delta <- 0.9 ## increase this if you get 'divergences' - even one means your model fit sucks!
taxAveStDPriorExpect <- 1.0
NTrees <- 1 ## number of random trees to sample and to fit the model to
NSplits <- 25 ## desired number of nodes per host timeBin
##

possibleFungidGenera <- c('Fungia_','Danafungia_','Cycloseris_','Pleuractis_')
## would be good here and for samples ID'd to genus to have some kind of weights for each possible species. i.e. it's really not possible that the acropora samples are A. cervicornis or palmata because their range is wrong, so even though I don't know which Acropora species exactly, I do have some info about the relative probabilities for many of them. similarly, some species might just be super rare in general and should likely be downweighted

taxaToExclude <- c('Aiptasia_sp','Entacmaea_quadricolor','Heteractis_aurora','Macrorhynchia_philippina','Mnemiopsis_sp','Sarcophyton_sp','Sinularia_polydactyla','Stichodactyla__gigantea','Xenia_umbellata') # don't currently have a way to incorporate these into the phylogeny...

sampleTipKey <- 'host_scientific_name'

filterfunction <- function(dfin) {
    levels(dfin[,sampleTipKey]) <- gsub(' ','_',levels(dfin[,sampleTipKey]))
    df2 <- droplevels(dfin[(dfin$tissue_compartment=='T') & !grepl(paste(c('Unknown|Missing',taxaToExclude),collapse='|'),dfin[,sampleTipKey],ignore.case=T) ,]) #| dfin$tissue_compartment=='S' | dfin$tissue_compartment=='M'
    contrasts(df2$ocean) <- 'contr.sum'
    contrasts(df2$complex_robust) <- 'contr.sum'
    contrasts(df2$host_scientific_name) <- 'contr.sum'
    #contrasts(df2$tissue_compartment) <- 'contr.sum'
    contrasts(df2$reef_name) <- 'contr.sum'
    contrasts(df2$colony_name) <- 'contr.sum'
    contrasts(df2$host_genus) <- 'contr.sum'
	contrasts(df2$host_clade_sensu_fukami) <- 'contr.sum'
	df2$longitude <- as.numeric(as.character(df2$longitude))
    df2$latitude <- as.numeric(as.character(df2$latitude))
    levels(df2[,sampleTipKey])[levels(df2[,sampleTipKey]) == 'Homophyllia_hillae'] <- "Homophyllia_bowerbanki"
    levels(df2[,sampleTipKey])[levels(df2[,sampleTipKey]) == 'Pocillopora_eydouxi'] <- "Pocillopora_grandis"
    levels(df2[,sampleTipKey])[levels(df2[,sampleTipKey]) == 'Pseudosiderastrea_tayami'] <- "Pseudosiderastrea_tayamai"
    return(df2)
}

hosttree <- read.tree(hosttreePath)
hosttree <- .compressTipLabel(hosttree)



map <- read.table(mapfilePath,header=T,sep='\t',comment.char='',check.names=F)
rownames(map) <- map[,'#SampleID']


newmap <- filterfunction(map)

fulltable <- t(read.table(fulltablePath, header=T, sep='\t', skip=1, comment.char='', row.names=1, check.names=F))

idx <- rownames(fulltable)[rownames(fulltable) %in% rownames(newmap) & rowSums(fulltable) >= minCountSamp]
newermap <- newmap[idx,]



## identify unique Scleractinian species in the data, and replace spaces with underscores

study.species <- gsub(' ', '_', levels(newermap[,sampleTipKey]))

## identify the Scleractinian species in the data that do not exist in the template tree

study.species.missing <- study.species[!study.species %in% grep(paste(c(attr(hosttree, "TipLabel"),'Fungid','not_applicable'),collapse='|'), study.species, ignore.case=T, value=T)]
generaOfUnknowns <- sapply(study.species.missing, function(x) strsplit(x, '_')[[1]][[1]])


y.old <- fulltable[idx, ]

y.old.filt <- apply(y.old,2,function(x) x/sum(x) > minPercent)
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


flavos <- rownames(tax[tax[,'Family']=='f__Flavobacteriaceae' & !is.na(tax[,'Family']),])

y.binary.filtered.flavos <- y.binary.filtered[,colnames(y.binary.filtered) %in% flavos]




bacttree <- read.tree(microbeTreePath)

bacttree.Y <- drop.tip(bacttree,bacttree$tip.label[!bacttree$tip.label %in% colnames(y.binary.filtered.flavos)])




## root the tree if it's unrooted
if(is.rooted(bacttree.Y)) {
    bacttree.Y.root <- reorder(bacttree.Y, order='pruningwise')
} else {
    bacttree.Y.root <- reorder(midpoint.root(bacttree.Y), order='pruningwise')
}
##

if(is.null(bacttree.Y.root$edge.length)) {
    bacttree.Y.root$edge.length <- rep(1,length(bacttree.Y.root$edge))
}


## summarize the putative taxa to be estimated
taxa <- factor(colnames(y.binary.filtered.flavos),levels=bacttree.Y.root$tip.label)
taxnames <- levels(taxa)
NTaxa <- length(taxnames)
##

## sort the OTU table so its entries match the tree's tip labels
y <- y.binary.filtered.flavos[,taxnames]
##

esttips <- colnames(y)
NEstTips <- length(esttips)


senddat <- melt(y, varnames=c('sample','tip'), value.name='present')
samplenames <- as.numeric(factor(senddat[,1]))
tipnames <- as.numeric(factor(senddat[,2], levels=esttips))
present <- senddat[,3]
NObs <- nrow(senddat)
NEstSamples <- length(unique(samplenames))

newermap <- newermap[levels(factor(senddat[,1])),]


NEstIntTaxNodes <- bacttree.Y.root$Nnode
NEstTaxNodes <- NEstTips + NEstIntTaxNodes - 1

taxEdges <- bacttree.Y.root$edge.length

allnodes <- unique(bacttree.Y.root$edge[,1])
NSamplesRaw <- nrow(y)
ancestors <- matrix(0, NEstTaxNodes + 1, NEstTaxNodes + 1)
for(node in 1:(NEstTaxNodes + 1)) {
    ancestors[, node] <- as.numeric(1:(NEstTaxNodes + 1) %in% c(Ancestors(bacttree.Y.root, node), node))
}
##

colnames(ancestors) <- rownames(ancestors) <- paste0('i',1:(NEstTaxNodes + 1))
colnames(ancestors)[1:NEstTips] <- rownames(ancestors)[1:NEstTips] <- paste0('t',colnames(y))

ancestors <- ancestors[-(NEstTips + 1), -(NEstTips + 1)]

newermap$sequencing_depth <- rowSums(y[,esttips])
newermap$log_sequencing_depth <- log(newermap$sequencing_depth)
newermap$log_sequencing_depth_scaled <- scale(newermap$log_sequencing_depth)





modelform <- ~ ocean + reef_name + log_sequencing_depth_scaled # + tissue_compartment + colony_name
allfactors <- attr(terms.formula(modelform), "term.labels")
NFactors <- length(allfactors)
allfactorder <- sapply(allfactors, function(x) sum(gregexpr(':', x, fixed=TRUE)[[1]] > 0))

modelMat <- model.matrix(modelform, model.frame(newermap,na.action=NULL))
modelMat[is.na(modelMat)] <- 0
sumconts <- names(attr(modelMat, "contrasts")[attr(modelMat, "contrasts")=='contr.sum'])
for(j in sumconts) {
    colnames(modelMat)[grep(j,colnames(modelMat))] <- paste0(j,levels(newermap[,j])[-nlevels(newermap[,j])]) ##this will not work when there are interaction effects!!!
}
modelMat <- modelMat[,-1]
NEffects <- ncol(modelMat)


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

NTimeBins <- ceiling(length(levels(newermap[,sampleTipKey])) / NSplits)

### starting here, generate multiple random samples of the map and trees (later to be summarized to define the time Bins and also to be separately used for replicate runs of the model)
sampleMap <- list()
hosttreesSampled <- list()
splittimes <- list()
boundaries <- matrix(NA,nrow=NTrees, ncol=NTimeBins-1)
for(i in 1:NTrees) {
    sampleMap[[i]] <- newermap
    fungidSps <- grep('Fungid',levels(sampleMap[[i]][,sampleTipKey]))
    levels(sampleMap[[i]][,sampleTipKey])[fungidSps] <- sample(grep(paste0(paste(possibleFungidGenera,collapse='|'),'_'), attr(hosttree, "TipLabel"), value=T), length(fungidSps)) ##assign unidentified Fungids to a random member of the group independently for each tree
    for (j in unique(generaOfUnknowns)) {
        possibleGenera <- attr(hosttree, "TipLabel")[!attr(hosttree, "TipLabel") %in% levels(sampleMap[[i]][,sampleTipKey])]
        if(!any(grepl(j,possibleGenera))) {
            possibleGenera <- attr(hosttree, "TipLabel")
        }
        levels(sampleMap[[i]][,sampleTipKey])[levels(sampleMap[[i]][,sampleTipKey]) %in% grep(paste0(j,'_'),study.species.missing,value=T)] <- sample(grep(paste0(j,'_'), possibleGenera, value=T), sum(generaOfUnknowns == j)) ## assign unidentified species to a random member of their genus independently for each tree
    }

	#filter the tree only contain the sampled (or assigned) species
	hosttreesSampled[[i]] <- ladderize(drop.tip(hosttree[[i]],hosttree[[i]]$tip.label[!hosttree[[i]]$tip.label %in% levels(sampleMap[[i]][,sampleTipKey])]))

	#divide total evolutionary time into chunks that contain approximately equal numbers of splits
	lttHosttreesSampled <- ltt(hosttreesSampled[[i]],log.lineages=F,plot=F)
	temp <- max(nodeHeights(hosttreesSampled[[i]])) - lttHosttreesSampled$times[-length(lttHosttreesSampled$times)]
	splittimes[[i]] <- split(temp, ceiling(seq_along(temp)/NSplits))
    if(length(splittimes[[i]][[NTimeBins]]) < NSplits/2) {
        splittimes[[i]][[NTimeBins - 1]] <- c(splittimes[[i]][[NTimeBins - 1]], splittimes[[i]][[NTimeBins]])
        splittimes[[i]][[NTimeBins]] <- NULL
        NTimeBins <- NTimeBins - 1
        boundaries <- boundaries <- matrix(NA,nrow=NTrees, ncol=NTimeBins-1)
    }
    
    #cut points in each phylogeny that would result in approximately equal numbers of splits per bin of time
    boundaries[i,] <- sapply(2:NTimeBins, function(x) max(splittimes[[i]][[x]]))

}


#average cut points among all the trees so the bins of time are consistent across replicates
meanBoundaries <- apply(boundaries,2,mean)


NHostTips <- length(hosttreesSampled[[1]]$tip.label)
NIntHostNodes <- hosttreesSampled[[1]]$Nnode
NHostNodes <- NIntHostNodes + NHostTips - 1


timeBinSizes <- list()
relativeTimeBinSizes <- list()
edgetobin <- list()
for(i in 1:NTrees) {
    nh <- nodeHeights(hosttreesSampled[[i]])
    maxNH <- max(nh)
    timeBinSizes[[i]] <- sapply(2:(length(meanBoundaries)+2),function(x) c(maxNH,meanBoundaries,0)[x-1] - c(maxNH,meanBoundaries,0)[x]) #size in millions of years of each bin (consistent across replicates /except/ for the oldest bin
    relativeTimeBinSizes[[i]] <- timeBinSizes[[i]] / sum(timeBinSizes[[i]]) #proportion of total tree height belonging to each bin (differs for each replicate due to the varying sizes of the oldest bin)
    nd <- maxNH - nh
    edgetobin[[i]] <- matrix(NA,ncol=NTimeBins,nrow=nrow(hosttreesSampled[[i]]$edge)) #create a matrix for each replicate that describes the proportion of each time bin that each branch exists for
    for(j in 1:NTimeBins) {
        if(j == 1) {
            allin <- which(nd[,2] >= meanBoundaries[j])
			allout <- which(nd[,1] <= meanBoundaries[j])
            cedge <- which((nd[,1] > meanBoundaries[j]) & (nd[,2] < meanBoundaries[j]))
            edgetobin[[i]][cedge,j] <- nd[cedge,1] - meanBoundaries[j]
        } else if(j == NTimeBins) {
            allin <- which(nd[,1] <= meanBoundaries[j-1])
            allout <- which(nd[,2] >= meanBoundaries[j-1])
            cedge <- which((nd[,1] > meanBoundaries[j-1]) & (nd[,2] < meanBoundaries[j-1]))
            edgetobin[[i]][cedge,j] <- meanBoundaries[j-1] - nd[cedge,2]
        } else {
            allin <- which((nd[,1] <= meanBoundaries[j-1]) & (nd[,2] >= meanBoundaries[j]))
            allout <- which((nd[,1] <= meanBoundaries[j]) | (nd[,2] >= meanBoundaries[j-1]))
            cedge1 <- which((nd[,1] < meanBoundaries[j-1]) & (nd[,1] > meanBoundaries[j]) & (nd[,2] < meanBoundaries[j]))
            edgetobin[[i]][cedge1,j] <- nd[cedge1,1] - meanBoundaries[j]
            cedge2 <- which((nd[,1] > meanBoundaries[j-1]) & (nd[,2] < meanBoundaries[j-1]) & (nd[,2] > meanBoundaries[j]))
            edgetobin[[i]][cedge2,j] <- meanBoundaries[j-1] - nd[cedge2,2]
            cedge3 <- which((nd[,1] > meanBoundaries[j-1]) & (nd[,2] < meanBoundaries[j]))
            edgetobin[[i]][cedge3,j] <- meanBoundaries[j-1] - meanBoundaries[j]
        }
        edgetobin[[i]][allin,j] <- hosttreesSampled[[i]]$edge.length[allin]
        edgetobin[[i]][allout,j] <- 0

		edgetobin[[i]][,j] <- edgetobin[[i]][,j] / timeBinSizes[[i]][j]
    }
    rownames(edgetobin[[i]]) <- hosttreesSampled[[i]]$edge[,2]
    edgetobin[[i]] <- edgetobin[[i]][order(as.numeric(rownames(edgetobin[[i]]))),]
}





### expand this so it's unique for each tree
hostAncestors <- list()
hostAncestorsExpanded <- list()
for (i in 1:NTrees) {
	hostAncestors[[i]] <- matrix(NA, NHostNodes+1, NHostNodes+1)
	for(node in 1:(NHostNodes+1)) {
	    hostAncestors[[i]][node, ] <- as.numeric(1:(NHostNodes+1) %in% c(Ancestors(hosttreesSampled[[i]], node), node))
	}
    colnames(hostAncestors[[i]]) <- rownames(hostAncestors[[i]]) <- paste0('i',1:(NHostNodes+1))
    colnames(hostAncestors[[i]])[1:NHostTips] <- rownames(hostAncestors[[i]])[1:NHostTips] <- hosttreesSampled[[i]]$tip.label
    hostAncestors[[i]] <- hostAncestors[[i]][-(NHostTips + 1), -(NHostTips + 1)]
    hostAncestorsExpanded[[i]] <- hostAncestors[[i]][as.character(sampleMap[[i]][,sampleTipKey]),]
    rownames(hostAncestorsExpanded[[i]]) <- rownames(sampleMap[[i]])
}
##


standat <- list()
for (i in 1:NTrees) {
    standat[[i]] <- list(NEstSamples=NEstSamples, NObs=NObs, NEstTaxNodes=NEstTaxNodes, NEstTips=NEstTips, NFactors=NFactors, NEffects=NEffects, present=present, samplenames=samplenames, tipnames=tipnames, factLevelMat=factLevelMat, ancestors=ancestors, modelMat=modelMat, hostAncestors=hostAncestors[[i]], hostAncestorsExpanded=hostAncestorsExpanded[[i]], edgetobin=edgetobin[[i]], NHostNodes=NHostNodes, NTimeBins=NTimeBins, taxAveStDPriorExpect=taxAveStDPriorExpect, timeBinSizes=relativeTimeBinSizes[[i]], taxEdges=taxEdges)
}
dir.create(outdir, recursive=T)

save.image(file.path(outdir,'setup.RData'))
save(ancestors,file=file.path(outdir,'ancestors.RData'))
save(meanBoundaries,file=file.path(outdir,'meanBoundaries.RData'))
save(hostAncestors,file=file.path(outdir,'hostAncestors.RData'))
save(hostAncestorsExpanded,file=file.path(outdir,'hostAncestorsExpanded.RData'))
save(hosttreesSampled,file=file.path(outdir,'sampledHostTrees.RData'))
save(timeBinSizes,file=file.path(outdir,'timeBinSizes.RData'))


## logistic regression
##
fit <- mclapply(1:NTrees, function(i) stan(file=modelPath, data=standat[[i]], control=list(adapt_delta=adapt_delta, max_treedepth=max_treedepth), iter=NIterations, thin=max(c(floor(NIterations/2000),1)), chains=NChains ))

save(fit, file=file.path(outdir,'fit.RData'))

for(i in 1:NTrees) {
    
    currplotdir <- file.path(outdir,paste0('tree_',i),'plots')
    currtabledir <- file.path(outdir,paste0('tree_',i),'tables','nodes')
    currdatadir <- file.path(outdir,paste0('tree_',i),'data')


    dir.create(currplotdir, recursive=T)
    dir.create(currtabledir, recursive=T)
    dir.create(currdatadir, recursive=T)


    pdf(file=file.path(currplotdir,'sampledTree.pdf'), width=25, height=15)
    plot(hosttreesSampled[[i]],cex=0.75)
    for(age in max(nodeHeights(hosttreesSampled[[i]]))-meanBoundaries) {
        lines(x=c(age,age),y=c(1,length(hosttreesSampled[[i]]$tip.label)),lwd=1)
    }
    graphics.off()


    vars1 <- extract(fit[[i]], pars='stDProps')[[1]]
    colnames(vars1) <- c(paste0('ADiv.',colnames(factLevelMat)), paste0('Specificity.',colnames(factLevelMat)), 'ADiv.host', 'host.specificity', 'taxa.prevalence')

    pdf(file=file.path(currplotdir,'scalesboxes.pdf'), width=25, height=15)
    boxplot(vars1, cex.axis=0.5, las=2)
    graphics.off()

    save(vars1,file=file.path(currdatadir,'stDProps.RData'))
    
    
    scaledAlphaDivEffects <- extract(fit[[i]], pars='scaledAlphaDivEffects')[[1]]
    
    save(scaledAlphaDivEffects,file=file.path(currdatadir,'scaledAlphaDivEffects.RData'))


    timeBinProps <- extract(fit[[i]], pars='timeBinProps')[[1]]

    pdf(file=file.path(currplotdir,'timeBinProps_boxes.pdf'), width=25, height=15)
    boxplot(timeBinProps, cex.axis=0.5, las=2)
    graphics.off()

    save(timeBinProps,file=file.path(currdatadir,'timeBinProps.RData'))
    
    
	meanBoundariesRounded <- round(meanBoundaries,1)
    relativeEvolRates <- extract(fit[[i]], pars='relativeEvolRates')[[1]]
    colnames(relativeEvolRates) <- c(paste0('before ',meanBoundariesRounded[1],' mya'), paste0(meanBoundariesRounded[1],' - ',meanBoundariesRounded[2],' mya'), paste0(meanBoundariesRounded[2],' - ',meanBoundariesRounded[3],' mya'), paste0(meanBoundariesRounded[3],' - ',meanBoundariesRounded[4],' mya'), paste0(meanBoundariesRounded[4],' - present'))

    pdf(file=file.path(currplotdir,'evolRatesRelToWeightedMean.pdf'), width=25, height=15)
    boxplot(relativeEvolRates, xlab='Time Period', ylab='Rate of Evolution Relative to Weighted Mean')
    graphics.off()

	save(relativeEvolRates,file=file.path(currdatadir,'relativeEvolRates.RData'))
    
    relativeEvolRate4p5 <- (relativeEvolRates[,4] * timeBinSizes[[1]][[4]] + relativeEvolRates[,5] * timeBinSizes[[1]][[5]]) / (timeBinSizes[[1]][[4]] + timeBinSizes[[1]][[5]])
    relativeEvolRatesMerged45 <- cbind(relativeEvolRates[,-c(4,5)],relativeEvolRate4p5)
    colnames(relativeEvolRates)[[4]] <- paste0(meanBoundariesRounded[3],' - present')
    
    pdf(file=file.path(currplotdir,'evolRatesRelToWeightedMeanMerged45.pdf'), width=25, height=15)
    boxplot(relativeEvolRatesMerged45, xlab='Time Period', ylab='Rate of Evolution Relative to Weighted Mean')
    graphics.off()
    
    pdf(file=file.path(currplotdir,'logEvolRatesRelToWeightedMeanMerged45.pdf'), width=25, height=15)
    boxplot(log(relativeEvolRatesMerged45), xlab='Time Period', ylab='Log Rate of Evolution Relative to Weighted Mean')
    graphics.off()

    sums <- summary(fit[[i]], pars='phyloLogVarMultRaw', probs=c(0.05,0.95), use_cache = F)

    sums3d <- array(NA, dim=c(NHostNodes - NHostTips, NEstTaxNodes - NEstTips, ncol(sums$summary)))
    for(effect in 1:(NHostNodes - NHostTips)) {
        sums3d[effect,,] <- sums$summary[(effect-1) * (NEstTaxNodes - NEstTips) + (1:(NEstTaxNodes - NEstTips)),]
    }
    dimnames(sums3d) <- list(colnames(hostAncestors[[i]])[(NHostTips + 1):NHostNodes], rownames(ancestors)[(NEstTips + 1):NEstTaxNodes], colnames(sums$summary))
    factorfilenames <- colnames(hostAncestors[[i]])[(NHostTips + 1):NHostNodes]

    for(effect in 1:(NHostNodes - NHostTips)) {
        cat('\t', file=file.path(currtabledir, paste0(factorfilenames[[effect]],'.txt')))
        write.table(sums3d[effect,,], file=file.path(currtabledir, paste0(factorfilenames[[effect]],'.txt')), sep='\t', quote=F,append=T)
    }

}

currplotdir <- file.path(outdir,'alltrees','plots')
currtabledir <- file.path(outdir,'alltrees','tables','nodes')
currdatadir <- file.path(outdir,'alltrees','data')

dir.create(currplotdir, recursive=T)
dir.create(currtabledir, recursive=T)
dir.create(currdatadir, recursive=T)




allfit <- sflist2stanfit(fit)



vars1 <- extract(allfit, pars='stDProps')[[1]]
colnames(vars1) <- c(paste0('ADiv.',colnames(factLevelMat)), paste0('Specificity.',colnames(factLevelMat)), 'ADiv.host', 'host.specificity')

pdf(file=file.path(currplotdir,'scalesboxes.pdf'), width=25, height=15)
boxplot(vars1, cex.axis=0.5, las=2)
graphics.off()

save(vars1,file=file.path(currdatadir,'stDProps.RData'))



timeBinProps <- extract(allfit, pars='timeBinProps')[[1]]

pdf(file=file.path(currplotdir,'timeBinPropsboxes.pdf'), width=25, height=15)
boxplot(timeBinProps, cex.axis=0.5, las=2)
graphics.off()

save(timeBinProps,file=file.path(currdatadir,'timeBinProps.RData'))



relativeEvolRates <- extract(allfit, pars='relativeEvolRates')[[1]]
colnames(relativeEvolRates) <- c(paste0('before ',meanBoundariesRounded[1],' mya'), paste0(meanBoundariesRounded[1],' - ',meanBoundariesRounded[2],' mya'), paste0(meanBoundariesRounded[2],' - ',meanBoundariesRounded[3],' mya'), paste0(meanBoundariesRounded[3],' - ',meanBoundariesRounded[4],' mya'), paste0(meanBoundariesRounded[4],' - present'))

pdf(file=file.path(currplotdir,'evolRatesRelToWeightedMean.pdf'), width=25, height=15)
boxplot(relativeEvolRates, xlab='Time Period', ylab='Rate of Evolution Relative to Weighted Mean')
graphics.off()

save(relativeEvolRates,file=file.path(currdatadir,'relativeEvolRates.RData'))

relativeEvolRate4p5 <- (relativeEvolRates[,4] * timeBinSizes[[1]][[4]] + relativeEvolRates[,5] * timeBinSizes[[1]][[5]]) / (timeBinSizes[[1]][[4]] + timeBinSizes[[1]][[5]])
relativeEvolRatesMerged45 <- cbind(relativeEvolRates[,-c(4,5)],relativeEvolRate4p5)
colnames(relativeEvolRates)[[4]] <- paste0(meanBoundariesRounded[3],' - present')

pdf(file=file.path(currplotdir,'evolRatesADivRelToWeightedMeanMerged45.pdf'), width=25, height=15)
boxplot(relativeEvolRatesMerged45, xlab='Time Period', ylab='Rate of Evolution Relative to Weighted Mean')
graphics.off()

pdf(file=file.path(currplotdir,'logEvolRatesADivRelToWeightedMeanMerged45.pdf'), width=25, height=15)
boxplot(log(relativeEvolRatesMerged45), xlab='Time Period', ylab='Log Rate of Evolution Relative to Weighted Mean')
graphics.off()

