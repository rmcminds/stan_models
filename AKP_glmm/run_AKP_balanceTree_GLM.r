library(coda)
library(rstan)
library(ape)
library(geiger)
library(phytools)
library(phangorn)
library(shinystan)
library(nlme)
library(reshape2)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


################################## customize options in this section ###########################

## locate all your files
outdir <- 'balance_out/'
otu_table <- 'feature-table.txt' ##tsv as produced by biom convert --to-tsv
intree <- 'tree.nwk'
mapfile <- 'map.txt'
modelfile <- 'AKP_balanceTree_GLM.stan'
##

## set some thresholds for filtering your data
mincountnode <- 5 #minimum count of a node in a sample
minsampsnode <- 1 #minimum number of samples in which a node must have a count of at least 'mincountnode' to be kept
mincountparent <- 25 #minimum count of a parent node for a balance to be estimated in a given sample
minsampsbalance <- 3 #minimum number of samples in which the balance is present after filtering (i.e. if it's been filtered to a number of samples lower than this, then don't estimate the balance in the model at all)
##
## mostly for testing purposes, trim heavily to only look at the very largest balances
mincountnode <- 10
minsampsnode <- 1
mincountparent <- 100
minsampsbalance <- 5
##

## place custom code within this function to filter your mapping file to contain a subset of its samples and to recast variables into their proper formats (e.g. factors aren't continuous, base level is correct)
filterfunction <- function(df2) {
    df2$treatment <- relevel(df2$treatment, 'C')
    df2$time <- relevel(df2$time, 'january2016')
    df2$temperature <- scale(as.numeric(gsub(',','.',as.character(df2$temperature))), center=T, scale=T)
    contrasts(df2$species) <- 'contr.sum'
    df2$tag <- as.factor(df2$tag)
    contrasts(df2$tag) <- 'contr.sum'
    contrasts(df2$plot) <- 'contr.sum'
    return(df2)
}
##

## specify your model formula for location effects
modelform <- ~ treatment * species * (time + temperature) + plot + tag
allfactors <- attr(terms.formula(modelform), "term.labels")
##

## specify your model formula for dispersion effects
dispform <- ~ treatment * species * (time + temperature)
##

## factors to split plots by
plotfacts <- c('treatment','time','species')
##

## Stan options
nchains <- 4 ## should run at least 4 chains to make sure they converge on the same answer (they will be run in parallel)
iterations <- 2000 ## will probably need >10,000? maybe start with 2, check convergence, double it, check, double, check, double, etc.?
max_treedepth <- 20 ## a warning will tell you if this needs to be increased
adapt_delta <- 0.99 ## increase this if you get 'divergences' - even one means your model fit sucks!
##


######################### and then just run the following code (hopefully) ##########################

## load the data
fulltable <- t(read.table(otu_table, header=T, sep='\t', comment.char='', skip=1, row.names=1))
bacttree <- read.tree(intree)
map <- read.table(mapfile,header=T,sep='\t',comment.char='',check.names=F)
rownames(map) <- map[,'#SampleID']
##

## filter the data
newmap <- filterfunction(map)
idx <- rownames(fulltable)[rownames(fulltable) %in% rownames(newmap)]
y.old <- fulltable[idx,colnames(fulltable) %in% bacttree$tip.label]
newermap <- newmap[idx,]
y.filt <- y.old[,apply(y.old,2,function(x) sum(x >= mincountnode) >= minsampsnode)]
bacttreeY <- drop.tip(bacttree,bacttree$tip.label[!bacttree$tip.label %in% colnames(y.filt)])
##


## root the tree if it's unrooted
if(is.rooted(bacttreeY)) {
    bacttreeY.root <- reorder(bacttreeY, order='pruningwise')
} else {
    bacttreeY.root <- reorder(midpoint.root(bacttreeY), order='pruningwise')
}
##

## summarize the putative taxa to be estimated
taxa <- factor(colnames(y.filt),levels=bacttreeY.root$tip.label)
taxnames <- levels(taxa)
NTaxa <- length(taxnames)
NIntTaxNodes <- bacttreeY.root$Nnode
NTaxNodes <- NTaxa + NIntTaxNodes
##

## sort the OTU table so its entries match the tree's tip labels
y <- y.filt[,taxnames]
##

## remove some potentially large objects to free up some memory
rm(fulltable,y.old,y.filt,bacttree,bacttreeY)
gc()
##

## expand the 'OTU table' so it contains entries for every node in the tree by summing the counts of all descendant taxa
allnodes <- unique(bacttreeY.root$edge[,1])
NSamplesRaw <- nrow(y)
z <- cbind(y, matrix(NA, NSamplesRaw, NIntTaxNodes) )
for(node in allnodes) {
    z[, node] <- apply(z[, bacttreeY.root$edge[bacttreeY.root$edge[,1]==node, 2]], 1, sum)
}
##

## keep only observations where the balance's sum is above the specified threshold remove nodes from the analysis entirely if they're not observed above thresholds at least a certain number of times
colnames(z)[(ncol(y)+1):ncol(z)] <- as.character((ncol(y)+1):ncol(z))
temp <- melt(z[, (ncol(y)+1):ncol(z)], varnames=c('sample','node'), value.name='parentcount')
temp <- temp[temp$parentcount >= mincountparent,]

halfedges <- bacttreeY.root$edge[!duplicated(bacttreeY.root$edge[,1]),]
rownames(halfedges) <- halfedges[,1]
halfedges <- halfedges[,-1]

temp$childcount <- sapply(1:nrow(temp), function(x) z[temp[x,1], halfedges[[as.character(temp[x,2])]]])

estnodes <- allnodes[sapply(allnodes, function(x) sum(x==temp$node) >= minsampsbalance)]

senddat <- temp[temp$node %in% estnodes,]

## remove some more potentially large objects to free up some memory
rm(z)
gc()
##

## filter samples from the map if they got entirely filtered out of the data, and sort the map so the samples are in the order of the levels being sent to Stan
newestmap <- newermap[levels(factor(senddat[,1])),]
##

## create the model matrix for location effects
modelMat <- model.matrix(modelform, data=newestmap)
##

modelFacts <- attr(terms.formula(modelform), 'factors')
modelOrds <- attr(terms.formula(modelform), 'order')
factLevOris <- attr(modelMat, 'assign') #which column of modelFacts gave rise to each column of subModelMat
expModelOrds <- c(0,modelOrds[factLevOris])
senddatFacts <- merge(senddat,modelMat,by.x=1,by.y=0)

## do not estimate factor levels if the number of samples that correspond to that level is the same as the number of samples in a lower order of the relevant factors. for instance, if there are two nutrientsyes:algae_treatmentHalimeda samples, but there are only two algae_treatmentHalimeda samples total, then only estimate the algae_treatmentHalimeda effect, because it is impossible to estimate both (you can't see how the interaction is different from the main effect because there are no samples that are algae_treatmentHalimeda without also being nutrientsyes).
fiidx <- NULL
NFactsPerTax <- NULL
factIndices <- NULL
for(o in estnodes) {
    nsamp <- colSums(senddatFacts[senddatFacts[,2]==o, 5:ncol(senddatFacts)])
    estimate <- sapply(1:length(nsamp), function(x) {
        higher <- c(strsplit(names(nsamp[x]),':')[[1]], '(Intercept)')
        higher <- higher[higher != names(nsamp[x])]
        all( nsamp[[x]] < nsamp[names(nsamp) %in% higher] ) & nsamp[[x]] > 0
    })
    factIndices <- c(factIndices, which(estimate))
    fiidx <- c(fiidx, sum(NFactsPerTax) + 1)
    NFactsPerTax <- c(NFactsPerTax, sum(estimate))
}
##

## create the model matrix for dispersion effects
dispMat <- model.matrix(dispform, data=newestmap)[,-1]
##

## organize and summarize model parameters to pass to Stan
samplenames <- as.numeric(factor(senddat[,1]))
nodenames <- as.numeric(factor(senddat[,2], levels=estnodes))
sampledcounts <- as.numeric(senddat[,3])
datacounts <- as.numeric(senddat[,4])


NEstTaxNodes <- length(estnodes)
NObs <- length(sampledcounts)
NEstSamples <- length(unique(samplenames))
NFactors <- ncol(modelMat)
NFactorsDisp <- ncol(dispMat)







standat <- list(NObs=NObs, NEstSamples=NEstSamples, sampledcounts=sampledcounts, datacounts=datacounts, samplenames=samplenames, nodenames=nodenames, NEstTaxNodes=NEstTaxNodes, NFactors=NFactors, modelMat=modelMat, dispMat=dispMat, NFactorsDisp=NFactorsDisp, fiidx=fiidx, NFactsPerTax=NFactsPerTax, factIndices=factIndices)
##


### run the model!
fit <- stan(file=modelfile, data=standat, control=list(adapt_delta=adapt_delta, max_treedepth=max_treedepth), iter=iterations, thin=floor(iterations/2000), chains=nchains )

dir.create(outdir)

save.image(file=paste0(outdir,'balanceTree_GLM_results.RData'))

shinystanres <- as.shinystan(fit, pars=c('tax_normFacts','factorDispersionMultipliers'))

save(shinystanres,file=paste0(outdir,'balanceTree_GLM_results_shinystan.RData'))



######################### plotting functions; maybe call these beta stage #######################

obsrats <- datacounts/sampledcounts

plotdat <- data.frame(obsrats=obsrats,samplenames=samplenames,nodenames=nodenames)
for(fact in plotfacts) {
	plotdat[,fact] <- sapply(plotdat$samplenames,function(x) return(newestmap[x,fact]))
}

level <- colnames(tax)[1:7]

dir.create(paste0(outdir,'balanceTree_GLM_plots/'))

for(i in unique(nodenames)) {
    pdf(file=paste0(outdir,'balanceTree_GLM_plots/',estnodes[[i]],'.pdf'),width=10)
    boxplot(update(modelform,obsrats~.),data=plotdat[plotdat$nodenames==i,],las=2,cex.axis=0.5)
    graphics.off()
    
    for(node in Children(bacttreeY.root, estnodes[[i]])) {
        tips <- bacttreeY.root$tip.label[Descendants(bacttreeY.root, node,type='tips')[[1]]]
        cat(paste0('node ',node,':\n'),file=paste0(outdir,'balanceTree_GLM_plots/',estnodes[[i]],'.txt'),append=T)
        if(length(tips) != 1) {
            cat(paste0(sort(unique(apply(tax[tips,level],1,paste,collapse='.'))),'\n'),file=paste0(outdir,'balanceTree_GLM_plots/',estnodes[[i]],'.txt'),append=T)
            cat(paste0(sort(unique(tips)),'\n'),file=paste0(outdir,'balanceTree_GLM_plots/',estnodes[[i]],'.txt'),append=T)
        } else {
            cat(paste0(paste(tax[tips,level],collapse='.'),'\n'),file=paste0(outdir,'balanceTree_GLM_plots/',estnodes[[i]],'.txt'),append=T)
            cat(paste0(tips,'\n'),file=paste0(outdir,'balanceTree_GLM_plots/',estnodes[[i]],'.txt'),append=T)
        }
    }
    
}


######################### summary functions; maybe call these beta stage #######################

sums <- summary(fit, pars='tax_normFacts', probs = c(0.05, 0.95))

sums3d <- array(NA, dim=c(NFactors,NEstTaxNodes,ncol(sums$summary)))
for(factor in 1:NFactors) {
    sums3d[factor,,] <- sums$summary[(factor-1) * NEstTaxNodes + (1:NEstTaxNodes),]
    write.table(sums3d[factor,,], file=paste0(outdir,'balanceTree_GLM_plots/',colnames(modelMat)[[factor]],'.txt'), sep='\t', quote=F)
}
dimnames(sums3d) <- list(colnames(modelMat)[1:NFactors],estnodes,colnames(sums$summary))
factorfilenames <- gsub(':','.',colnames(modelMat))

for(factor in 1:NFactors) {
    cat('\t', file=paste0(outdir,'balanceTree_GLM_plots/',factorfilenames[[factor]],'.txt'))
    write.table(sums3d[factor,,], file=paste0(outdir,'balanceTree_GLM_plots/',factorfilenames[[factor]],'.txt'), sep='\t', quote=F,append=T)
}

sumsDisp <- summary(fit, pars='factorDispersionMultipliers', probs = c(0.05, 0.95))

rownames(sumsDisp$summary) <- colnames(dispMat)
write.table(sumsDisp$summary, file=paste0(outdir,'balanceTree_GLM_plots/factorDispersionMultipliers.txt'), sep='\t', quote=F)
