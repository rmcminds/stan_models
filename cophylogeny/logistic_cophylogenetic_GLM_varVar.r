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


spectree <- 'species_rooted_mapped_labels.tree'
mapfile <- 'GCMP_symbiodinium_map2.txt'
outdir <- 'logistic_GLM_out_test/'
fulltable <- t(read.table('species_table.txt', header=T, sep='\t', comment.char='', row.names=1, check.names=F))
modelfile <- 'logistic_cophylogenetic_GLM_varVar.stan'
hosttreef <- 'combined_trees.newick'
minsampsbalance <- 1
mincountsamp <- 5


## Stan options
nchains <- 1 ## this is per tree; since I'm doing a large number of trees in parallel i'll just do one chain for each
iterations <- 5000 ## will probably need >10,000? maybe start with 2, check convergence, double it, check, double, check, double, etc.?
max_treedepth <- 15 ## a warning will tell you if this needs to be increased
adapt_delta <- 0.9 ## increase this if you get 'divergences' - even one means your model fit sucks!
taxAveStDPriorExpect <- 1.0
nparalleltrees <- 1 ## number of random trees to sample and to fit the model to
nsplits <- 15 ## desired number of nodes per host timeBin
##

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

hosttree <- read.tree(hosttreef)

hosttree <- .compressTipLabel(hosttree)


bacttree <- read.tree(spectree)


## root the tree if it's unrooted
if(is.rooted(bacttree)) {
    bacttree.root <- reorder(bacttree, order='pruningwise')
} else {
    bacttree.root <- reorder(midpoint.root(bacttree), order='pruningwise')
}
##

if(is.null(bacttree.root$edge.length)) {
    bacttree.root.chronos <- bacttree.root
    bacttree.root.chronos$edge.length <- rep(1,length(bacttree.root.chronos$edge))
    bacttree.root.chronos <- chronos(bacttree.root.chronos, model = "discrete", control = chronos.control(nb.rate.cat = 1))
} else {
    bacttree.root.chronos <- chronos(bacttree.root)
}

class(bacttree.root.chronos) <- 'phylo'


newtable <- fulltable[-nrow(fulltable),]
mode(newtable) <- 'numeric'


map <- read.table(mapfile,header=T,sep='\t',comment.char='',check.names=F)
rownames(map) <- map[,'#SampleID']


newmap <- filterfunction(map)
idx <- rownames(newtable)[rownames(newtable) %in% rownames(newmap) & rowSums(newtable) >= mincountsamp]
y.old <- newtable[idx,colnames(newtable) %in% bacttree.root.chronos$tip.label]
newermap <- newmap[idx,]
bacttree.root.chronos.Y <- drop.tip(bacttree.root.chronos,bacttree.root.chronos$tip.label[!bacttree.root.chronos$tip.label %in% colnames(y.old)])



## summarize the putative taxa to be estimated
taxa <- factor(colnames(y.old),levels=bacttree.root.chronos.Y$tip.label)
taxnames <- levels(taxa)
NTaxa <- length(taxnames)
##

## sort the OTU table so its entries match the tree's tip labels
y <- y.old[,taxnames]
##
ybinary <- apply(y,2,function(x) x > 0)
mode(ybinary) <- 'numeric'

yb <- ybinary[,colSums(ybinary) >= minsampsbalance]
esttips <- colnames(yb)
NEstTips <- length(esttips)


senddat <- melt(yb, varnames=c('sample','tip'), value.name='present')
samplenames <- as.numeric(factor(senddat[,1]))
tipnames <- as.numeric(factor(senddat[,2], levels=esttips))
present <- senddat[,3]
NObs <- nrow(senddat)
NEstSamples <- length(unique(samplenames))

newermap <- newermap[levels(factor(senddat[,1])),]

bacttree.root.chronos.Y.filtered <- drop.tip(bacttree.root.chronos.Y, bacttree.root.chronos.Y$tip.label[!bacttree.root.chronos.Y$tip.label %in% esttips])

NEstIntTaxNodes <- bacttree.root.chronos.Y.filtered$Nnode
NEstTaxNodes <- NEstTips + NEstIntTaxNodes - 1

taxEdges <- bacttree.root.chronos.Y.filtered$edge.length

allnodes <- unique(bacttree.root.chronos.Y.filtered$edge[,1])
NSamplesRaw <- nrow(yb)
ancestors <- matrix(0, NEstTaxNodes + 1, NEstTaxNodes + 1)
for(node in 1:NEstTaxNodes) {
    ancestors[, node] <- as.numeric(1:(NEstTaxNodes + 1) %in% c(Ancestors(bacttree.root.chronos.Y.filtered, node), node))
}
##

colnames(ancestors) <- rownames(ancestors) <- paste0('i',1:(NEstTaxNodes + 1))
colnames(ancestors)[1:NEstTips] <- rownames(ancestors)[1:NEstTips] <- paste0('t',colnames(yb))

ancestors <- ancestors[-(NEstTips + 1), -(NEstTips + 1)]

newermap$sequencing_depth <- rowSums(y[,esttips])
newermap$log_sequencing_depth <- log(newermap$sequencing_depth)
newermap$log_sequencing_depth_scaled <- scale(newermap$log_sequencing_depth)





modelform <- ~ ocean + reef_name + colony_name + tissue_compartment + log_sequencing_depth_scaled
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



## identify unique Scleractinian species in the data, and replace spaces with underscores

study.species <- gsub(' ', '_', levels(newermap[,sampleTipKey]))

## identify the Scleractinian species in the data that do not exist in the template tree

study.species.missing <- study.species[!study.species %in% c(attr(hosttree, "TipLabel"),'Fungid_sp','not_applicable')]
generaOfUnknowns <- sapply(study.species.missing, function(x) strsplit(x, '_')[[1]][[1]])


NTimeBins <- ceiling(length(levels(newermap[,sampleTipKey])) / nsplits)

### starting here, generate multiple random samples of the map and trees (later to be summarized to define the time Bins and also to be separately used for replicate runs of the model)
sampleMap <- list()
hosttreesSampled <- list()
splittimes <- list()
boundaries <- matrix(NA,nrow=nparalleltrees, ncol=NTimeBins-1)
for(i in 1:nparalleltrees) {
    sampleMap[[i]] <- newermap
    levels(sampleMap[[i]][,sampleTipKey])[levels(sampleMap[[i]][,sampleTipKey]) == 'Fungid_sp'] <- sample(grep(paste0(paste(possibleFungidGenera,collapse='|'),'_'), attr(hosttree, "TipLabel"), value=T), 1) ##assign unidentified Fungids to a random member of the group independently for each tree
    for (j in unique(generaOfUnknowns)) {
        possibleGenera <- attr(hosttree, "TipLabel")[!attr(hosttree, "TipLabel") %in% levels(sampleMap[[i]][,sampleTipKey])]
        levels(sampleMap[[i]][,sampleTipKey])[levels(sampleMap[[i]][,sampleTipKey]) %in% grep(paste0(j,'_'),study.species.missing,value=T)] <- sample(grep(paste0(j,'_'), possibleGenera, value=T), sum(generaOfUnknowns == j)) ## assign unidentified species to a random member of their genus independently for each tree
    }

	#filter the tree only contain the sampled (or assigned) species
	hosttreesSampled[[i]] <- ladderize(drop.tip(hosttree[[i]],hosttree[[i]]$tip.label[!hosttree[[i]]$tip.label %in% levels(sampleMap[[i]][,sampleTipKey])]))

	#divide total evolutionary time into chunks that contain approximately equal numbers of splits
	lttHosttreesSampled <- ltt(hosttreesSampled[[i]],log.lineages=F,plot=F)
	temp <- max(nodeHeights(hosttreesSampled[[i]])) - lttHosttreesSampled$times[-length(lttHosttreesSampled$times)]
	splittimes[[i]] <- split(temp, ceiling(seq_along(temp)/nsplits))
    if(length(splittimes[[i]][[NTimeBins]]) < nsplits/2) {
        splittimes[[i]][[NTimeBins - 1]] <- c(splittimes[[i]][[NTimeBins - 1]], splittimes[[i]][[NTimeBins]])
        splittimes[[i]][[NTimeBins]] <- NULL
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
for(i in 1:nparalleltrees) {
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
for (i in 1:nparalleltrees) {
	hostAncestors[[i]] <- matrix(0, NHostNodes+1, NHostNodes+1)
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
for (i in 1:nparalleltrees) {
    standat[[i]] <- list(NEstSamples=NEstSamples, NObs=NObs, NEstTaxNodes=NEstTaxNodes, NEstTips=NEstTips, NFactors=NFactors, NEffects=NEffects, present=present, samplenames=samplenames, tipnames=tipnames, factLevelMat=factLevelMat, ancestors=ancestors, modelMat=modelMat, hostAncestors=hostAncestors[[i]], hostAncestorsExpanded=hostAncestorsExpanded[[i]], edgetobin=edgetobin[[i]], NHostNodes=NHostNodes, NTimeBins=NTimeBins, taxAveStDPriorExpect=taxAveStDPriorExpect, timeBinSizes=relativeTimeBinSizes[[i]], taxEdges=taxEdges)
}
dir.create(outdir, recursive=T)

save.image(paste0(outdir,'image.RData'))

## logistic regression
##
fit <- mclapply(1:nparalleltrees, function(i) stan(file=modelfile, data=standat[[i]], control=list(adapt_delta=adapt_delta, max_treedepth=max_treedepth), iter=iterations, thin=floor(iterations/2000), chains=nchains ))

dir.create(outdir, recursive=T)

save.image(paste0(outdir,'image.RData'))

save(ancestors,file=paste0(outdir,'ancestors.RData'))
save(meanBoundaries,file=paste0(outdir,'meanBoundaries.RData'))
save(hostAncestors,file=paste0(outdir,'hostAncestors.RData'))
save(hostAncestorsExpanded,file=paste0(outdir,'hostAncestorsExpanded.RData'))
save(hosttreesSampled,file=paste0(outdir,'sampledHostTrees.RData'))
save(timeBinSizes,file=paste0(outdir,'timeBinSizes.RData'))

effectList <- list()

for(i in 1:nparalleltrees) {
    
    currplotdir <- paste0(outdir,'tree_',i,'/plots/')
    currtabledir <- paste0(outdir,'tree_',i,'/tables/nodes/')
    currdatadir <- paste0(outdir,'tree_',i,'/data/')


    dir.create(currplotdir, recursive=T)
    dir.create(currtabledir, recursive=T)
    dir.create(currdatadir, recursive=T)


    pdf(file=paste0(currplotdir,'sampledTree.pdf'), width=25, height=15)
    plot(hosttreesSampled[[i]],cex=0.75)
    for(age in max(nodeHeights(hosttreesSampled[[i]]))-meanBoundaries) {
        lines(x=c(age,age),y=c(1,length(hosttreesSampled[[i]]$tip.label)),lwd=1)
    }
    graphics.off()


    vars1 <- extract(fit[[i]], pars='stDProps')[[1]]
    colnames(vars1) <- c(paste0('ADiv.',colnames(factLevelMat)), paste0('Specificity.',colnames(factLevelMat)), 'ADiv.host', 'cophylogeny')

    pdf(file=paste0(currplotdir,'scalesboxes.pdf'), width=25, height=15)
    boxplot(vars1, cex.axis=0.5, las=2)
    graphics.off()

    save(vars1,file=paste0(currdatadir,'stDProps.RData'))
    
    
    scaledAlphaDivEffects <- extract(fit[[i]], pars='scaledAlphaDivEffects')[[1]]
    
    save(scaledAlphaDivEffects,file=paste0(currdatadir,'scaledAlphaDivEffects.RData'))


    timeBinProps <- extract(fit[[i]], pars='timeBinProps')[[1]]

    pdf(file=paste0(currplotdir,'timeBinProps_boxes.pdf'), width=25, height=15)
    boxplot(timeBinProps, cex.axis=0.5, las=2)
    graphics.off()

    save(timeBinProps,file=paste0(currdatadir,'timeBinProps.RData'))
    
    
	meanBoundariesRounded <- round(meanBoundaries,1)
    relativeEvolRates <- extract(fit[[i]], pars='relativeEvolRates')[[1]]
    colnames(relativeEvolRates) <- c(paste0('before ',meanBoundariesRounded[1],' mya'), paste0(meanBoundariesRounded[1],' - ',meanBoundariesRounded[2],' mya'), paste0(meanBoundariesRounded[2],' - ',meanBoundariesRounded[3],' mya'), paste0(meanBoundariesRounded[3],' - ',meanBoundariesRounded[4],' mya'), paste0(meanBoundariesRounded[4],' - present'))

    pdf(file=paste0(currplotdir,'evolRatesRelToWeightedMean.pdf'), width=25, height=15)
    boxplot(relativeEvolRates, xlab='Time Period', ylab='Rate of Evolution Relative to Weighted Mean')
    graphics.off()

	save(relativeEvolRates,file=paste0(currdatadir,'relativeEvolRates.RData'))
    
    relativeEvolRate4p5 <- (relativeEvolRates[,4] * timeBinSizes[[1]][[4]] + relativeEvolRates[,5] * timeBinSizes[[1]][[5]]) / (timeBinSizes[[1]][[4]] + timeBinSizes[[1]][[5]])
    relativeEvolRatesMerged45 <- cbind(relativeEvolRates[,-c(4,5)],relativeEvolRate4p5)
    colnames(relativeEvolRates)[[4]] <- paste0(meanBoundariesRounded[3],' - present')
    
    pdf(file=paste0(currplotdir,'evolRatesRelToWeightedMeanMerged45.pdf'), width=25, height=15)
    boxplot(relativeEvolRatesMerged45, xlab='Time Period', ylab='Rate of Evolution Relative to Weighted Mean')
    graphics.off()
    
    pdf(file=paste0(currplotdir,'logEvolRatesRelToWeightedMeanMerged45.pdf'), width=25, height=15)
    boxplot(log(relativeEvolRatesMerged45), xlab='Time Period', ylab='Log Rate of Evolution Relative to Weighted Mean')
    graphics.off()

    sums <- summary(fit[[i]], pars='phyloLogVarMultRaw', probs=c(0.05,0.95), use_cache = F)

    sums3d <- array(NA, dim=c(NHostNodes - NHostTips, NEstTaxNodes - NEstTips, ncol(sums$summary)))
    for(effect in 1:(NHostNodes - NHostTips)) {
        sums3d[effect,,] <- sums$summary[(effect-1) * (NEstTaxNodes - NEstTips) + (1:(NEstTaxNodes - NEstTips)),]
    }
    dimnames(sums3d) <- list(colnames(hostAncestors[[i]]), rownames(ancestors), colnames(sums$summary))
    factorfilenames <- colnames(hostAncestors[[i]])

    for(effect in 1:(NHostNodes)) {
        cat('\t', file=paste0(currtabledir, factorfilenames[[effect]],'.txt'))
        write.table(sums3d[effect,,], file=paste0(currtabledir, factorfilenames[[effect]],'.txt'), sep='\t', quote=F,append=T)
    }

    dir.create(currtabledir, recursive=T)


    sums <- summary(fit[[i]], pars='scaledTipEffects', probs=c(0.05,0.95), use_cache = F)
    
    sums3d <- array(NA, dim=c(NEffects + NHostNodes, NEstTips, ncol(sums$summary)))
    for(effect in 1:(NEffects + NHostNodes)) {
        sums3d[effect,,] <- sums$summary[(effect-1) * NEstTips + (1:NEstTips),]
    }
    dimnames(sums3d) <- list(c(colnames(modelMat)[1:NEffects],colnames(hostAncestors[[i]])), colnames(ancestors), colnames(sums$summary))
    factorfilenames <- c(gsub(':','.',colnames(modelMat)), colnames(hostAncestors[[i]]))
    
    for(effect in 1:(NEffects + NHostNodes)) {
        cat('\t', file=paste0(currtabledir, factorfilenames[[effect]],'_summed_tips.txt'))
        write.table(sums3d[effect,,], file=paste0(currtabledir, factorfilenames[[effect]],'_summed_tips.txt'), sep='\t', quote=F,append=T)
    }
    
    ext <- extract(fit[[i]], pars='scaledTipEffects')[[1]]
    
    for(j in sumconts) {
        hi <- ext[,factLevelMat[,j] == 1,]
        sup <- -apply(hi, 1, colSums)
        rownames(sup) <- colnames(ancestors)
        sup2 <- sup
        dim(sup2) <- c(dim(sup),1)
        yeah <- monitor(aperm(sup2,c(2,3,1)), warmup=0, probs=c(0.05,0.95))
        rownames(yeah) <- colnames(ancestors)
        cat('\t', file=paste0(currtabledir,j,levels(newermap[,j])[nlevels(newermap[,j])],'_summed_tips.txt'))
        write.table(yeah, file=paste0(currtabledir,j,levels(newermap[,j])[nlevels(newermap[,j])],'_summed_tips.txt'), sep='\t', quote=F,append=T)
    }
    
    save(ext,file=paste0(currdatadir,'scaledTipEffects.RData'))

}

currplotdir <- paste0(outdir,'alltrees/plots/')
currtabledir <- paste0(outdir,'alltrees/tables/nodes/')
currdatadir <- paste0(outdir,'alltrees/data/')

dir.create(currplotdir, recursive=T)
dir.create(currtabledir, recursive=T)
dir.create(currdatadir, recursive=T)


save(effectList,file=paste0(currdatadir,'scaledADivAndTaxNodeEffects.RData'))




allfit <- sflist2stanfit(fit)



vars1 <- extract(allfit, pars='stDProps')[[1]]
colnames(vars1) <- c(paste0('ADiv.',colnames(factLevelMat)), paste0('Specificity.',colnames(factLevelMat)), 'ADiv.host', 'cophylogeny')

pdf(file=paste0(currplotdir,'scalesboxes.pdf'), width=25, height=15)
boxplot(vars1, cex.axis=0.5, las=2)
graphics.off()

save(vars1,file=paste0(currdatadir,'stDProps.RData'))



timeBinPropsADiv <- extract(allfit, pars='timeBinProps')[[1]]

pdf(file=paste0(currplotdir,'timeBinPropsboxes.pdf'), width=25, height=15)
boxplot(timeBinPropsADiv, cex.axis=0.5, las=2)
graphics.off()

save(timeBinPropsADiv,file=paste0(currdatadir,'timeBinPropsADiv.RData'))


timeBinPropsSpec <- extract(allfit, pars='timeBinPropsSpec')[[1]]

pdf(file=paste0(currplotdir,'timeBinPropsSpecboxes.pdf'), width=25, height=15)
boxplot(timeBinPropsSpec, cex.axis=0.5, las=2)
graphics.off()

save(timeBinPropsSpec,file=paste0(currdatadir,'timeBinPropsSpec.RData'))



relativeEvolRates <- extract(allfit, pars='relativeEvolRates')[[1]]
colnames(relativeEvolRates) <- c(paste0('before ',meanBoundariesRounded[1],' mya'), paste0(meanBoundariesRounded[1],' - ',meanBoundariesRounded[2],' mya'), paste0(meanBoundariesRounded[2],' - ',meanBoundariesRounded[3],' mya'), paste0(meanBoundariesRounded[3],' - ',meanBoundariesRounded[4],' mya'), paste0(meanBoundariesRounded[4],' - present'))

pdf(file=paste0(currplotdir,'evolRatesADivRelToWeightedMean.pdf'), width=25, height=15)
boxplot(relativeEvolRates, xlab='Time Period', ylab='Rate of Evolution Relative to Weighted Mean')
graphics.off()

save(relativeEvolRates,file=paste0(currdatadir,'relativeEvolRates.RData'))

relativeEvolRateADiv4p5 <- (relativeEvolRates[,4] * timeBinSizes[[1]][[4]] + relativeEvolRates[,5] * timeBinSizes[[1]][[5]]) / (timeBinSizes[[1]][[4]] + timeBinSizes[[1]][[5]])
relativeEvolRatesMerged45 <- cbind(relativeEvolRates[,-c(4,5)],relativeEvolRateADiv4p5)
colnames(relativeEvolRates)[[4]] <- paste0(meanBoundariesRounded[3],' - present')

pdf(file=paste0(currplotdir,'evolRatesADivRelToWeightedMeanMerged45.pdf'), width=25, height=15)
boxplot(relativeEvolRatesMerged45, xlab='Time Period', ylab='Rate of Evolution Relative to Weighted Mean')
graphics.off()

pdf(file=paste0(currplotdir,'logEvolRatesADivRelToWeightedMeanMerged45.pdf'), width=25, height=15)
boxplot(log(relativeEvolRatesMerged45), xlab='Time Period', ylab='Log Rate of Evolution Relative to Weighted Mean')
graphics.off()



relativeEvolRatesSpec <- extract(allfit, pars='relativeEvolRatesSpec')[[1]]
colnames(relativeEvolRatesSpec) <- c(paste0('before ',meanBoundariesRounded[1],' mya'), paste0(meanBoundariesRounded[1],' - ',meanBoundariesRounded[2],' mya'), paste0(meanBoundariesRounded[2],' - ',meanBoundariesRounded[3],' mya'), paste0(meanBoundariesRounded[3],' - ',meanBoundariesRounded[4],' mya'), paste0(meanBoundariesRounded[4],' - present'))

pdf(file=paste0(currplotdir,'evolRatesSpecRelToWeightedMean.pdf'), width=25, height=15)
boxplot(relativeEvolRatesSpec, xlab='Time Period', ylab='Rate of Evolution Relative to Weighted Mean')
graphics.off()

save(relativeEvolRatesSpec,file=paste0(currdatadir,'relativeEvolRatesSpec.RData'))

relativeEvolRateSpec4p5 <- (relativeEvolRatesSpec[,4] * timeBinSizes[[1]][[4]] + relativeEvolRatesSpec[,5] * timeBinSizes[[1]][[5]]) / (timeBinSizes[[1]][[4]] + timeBinSizes[[1]][[5]])
relativeEvolRatesSpecMerged45 <- cbind(relativeEvolRatesSpec[,-c(4,5)],relativeEvolRateSpec4p5)
colnames(relativeEvolRatesSpec)[[4]] <- paste0(meanBoundariesRounded[3],' - present')

pdf(file=paste0(currplotdir,'evolRatesSpecRelToWeightedMeanMerged45.pdf'), width=25, height=15)
boxplot(relativeEvolRatesSpecMerged45, xlab='Time Period', ylab='Rate of Evolution Relative to Weighted Mean')
graphics.off()

pdf(file=paste0(currplotdir,'logEvolRatesSpecRelToWeightedMeanMerged45.pdf'), width=25, height=15)
boxplot(log(relativeEvolRatesSpecMerged45), xlab='Time Period', ylab='Log Rate of Evolution Relative to Weighted Mean')
graphics.off()
