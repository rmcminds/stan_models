library(ape)
library(phangorn)
library(geosphere)
library(rstan)
library(geiger)
library(phytools)
library(shinystan)
library(nlme)
library(reshape2)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


spectree <- 'symbionttree.tree'
mapfile <- 'map.txt'
outdir <- 'logistic_GLM_out/'
fulltable <- t(read.table('species_table.txt', header=T, sep='\t', comment.char='', row.names=1, check.names=F))
modelfile <- 'logistic_GLM.stan'
minsampsbalance <- 5
mincountsamp <- 5

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
    return(df2)
}



bacttree <- read.tree(spectree)

newtable <- fulltable[-nrow(fulltable),]
mode(newtable) <- 'numeric'


map <- read.table(mapfile,header=T,sep='\t',comment.char='',check.names=F)
rownames(map) <- map[,'#SampleID']


newmap <- filterfunction(map)
idx <- rownames(newtable)[rownames(newtable) %in% rownames(newmap) & rowSums(newtable) >= mincountsamp]
y.old <- newtable[idx,colnames(newtable) %in% bacttree$tip.label]
newermap <- newmap[idx,]
bacttreeY <- drop.tip(bacttree,bacttree$tip.label[!bacttree$tip.label %in% colnames(y.old)])

## root the tree if it's unrooted
if(is.rooted(bacttreeY)) {
    bacttreeY.root <- reorder(bacttreeY, order='pruningwise')
} else {
    bacttreeY.root <- reorder(midpoint.root(bacttreeY), order='pruningwise')
}
##


## summarize the putative taxa to be estimated
taxa <- factor(colnames(y.old),levels=bacttreeY.root$tip.label)
taxnames <- levels(taxa)
NTaxa <- length(taxnames)
NIntTaxNodes <- bacttreeY.root$Nnode
NTaxNodes <- NTaxa + NIntTaxNodes
##

## sort the OTU table so its entries match the tree's tip labels
y <- y.old[,taxnames]
##

## expand the 'OTU table' so it contains entries for every node in the tree by summing the counts of all descendant taxa
allnodes <- unique(bacttreeY.root$edge[,1])
NSamplesRaw <- nrow(y)
z <- cbind(y, matrix(NA, NSamplesRaw, NIntTaxNodes) )
for(node in allnodes) {
    z[, node] <- apply(z[, bacttreeY.root$edge[bacttreeY.root$edge[,1]==node, 2]], 1, sum)
}
##

zb <- apply(z,2,function(x) x > 0)
mode(zb) <- 'numeric'

colnames(zb)[1:ncol(y)] <- paste0('t',colnames(zb)[1:ncol(y)])
colnames(zb)[(1 + ncol(y)):ncol(z)] <- paste0('i',(1 + ncol(y)):ncol(z))








## logistic regression


modelform <- ~ ocean + reef_name + complex_robust + host_clade_sensu_fukami + host_genus + host_scientific_name + colony_name + tissue_compartment + log_sequencing_depth
allfactors <- attr(terms.formula(modelform), "term.labels")




## Stan options
nchains <- 4 ## should run at least 4 chains to make sure they converge on the same answer (they will be run in parallel)
iterations <- 5000 ## will probably need >10,000? maybe start with 2, check convergence, double it, check, double, check, double, etc.?
max_treedepth <- 15 ## a warning will tell you if this needs to be increased
adapt_delta <- 0.9 ## increase this if you get 'divergences' - even one means your model fit sucks!
##

newermap$sequencing_depth <- rowSums(y)

newermap$log_sequencing_depth <- log(newermap$sequencing_depth)



tempdat <- melt(zb, varnames=c('sample','node'), value.name='present')
tempdat$parentpres <- apply(tempdat, 1, function(x) {
    zb[x[[1]], Ancestors(bacttreeY.root,match(x[[2]],colnames(zb)),'parent')]
})


estnodes <- colnames(zb)[sapply(colnames(zb), function(x) sum(x==tempdat$node) >= minsampsbalance)]

senddat <- tempdat[tempdat$node %in% estnodes,]





newestmap <- newermap[levels(factor(senddat[,1])),]




modelMat <- model.matrix(modelform, model.frame(newestmap,na.action=NULL))
modelMat[is.na(modelMat)] <- 0
##

sumconts <- names(attr(modelMat, "contrasts")[attr(modelMat, "contrasts")=='contr.sum'])

##this will not work when there are interaction effects!!!
for(i in sumconts) {
	colnames(modelMat)[grep(i,colnames(modelMat))] <- paste0(i,levels(newermap[,i])[-nlevels(newermap[,i])])
}

senddatFacts <- merge(senddat,modelMat,by.x=1,by.y=0)

## do not estimate factor levels if the number of samples that correspond to that level is the same as the number of samples in a lower order of the relevant factors. for instance, if there are two nutrientsyes:algae_treatmentHalimeda samples, but there are only two algae_treatmentHalimeda samples total, then only estimate the algae_treatmentHalimeda effect, because it is impossible to estimate both (you can't see how the interaction is different from the main effect because there are no samples that are algae_treatmentHalimeda without also being nutrientsyes).
fiidx <- NULL
NFactsPerTax <- NULL
factIndices <- NULL
for(o in estnodes) {
    nsamp <- colSums(senddatFacts[senddatFacts[,2]==o, 5:ncol(senddatFacts)] > 0, na.rm=T)
    estimate <- sapply(1:length(nsamp), function(x) {
        higher <- c(strsplit(names(nsamp[x]),':')[[1]], '(Intercept)')
        higher <- higher[higher != names(nsamp[x])]
        any(c(!modelMat[,names(nsamp)[[x]]] %in% -1:1, all(nsamp[[x]] < nsamp[names(nsamp) %in% higher]) & nsamp[[x]] > 0))
    })
    factIndices <- c(factIndices, which(estimate))
    fiidx <- c(fiidx, sum(NFactsPerTax) + 1)
    NFactsPerTax <- c(NFactsPerTax, sum(estimate))
}
##


samplenames <- as.numeric(factor(senddat[,1]))
nodenames <- as.numeric(factor(senddat[,2], levels=estnodes))
present <- senddat[,3]


NEstTaxNodes <- length(estnodes)
NObs <- nrow(senddat)
NEstSamples <- length(unique(samplenames))
NEffects <- ncol(modelMat)
NFactors <- length(allfactors)+1
sampleCount <- newestmap$sequencing_depth

allfactorder <- sapply(allfactors, function(x) sum(gregexpr(':', x, fixed=TRUE)[[1]] > 0))

factLevelMat <- matrix(NA, NEffects, NFactors)
colnames(factLevelMat) <- c('Intercept',allfactors)
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

factLevelMat[,'Intercept'] <- 0
factLevelMat[1,1] <- 1


standat <- list(NObs=NObs, NEstSamples=NEstSamples, samplenames=samplenames, nodenames=nodenames, NEstTaxNodes=NEstTaxNodes, NEffects=NEffects, modelMat=modelMat, NFactors=NFactors, factLevelMat=factLevelMat, present=present, fiidx=fiidx, NFactsPerTax=NFactsPerTax, factIndices=factIndices)
##
fit <- stan(file=modelfile, data=standat, control=list(adapt_delta=adapt_delta, max_treedepth=max_treedepth), iter=iterations, thin=floor(iterations/2000), chains=nchains )




dir.create(paste0(outdir,'balanceTree_GLM_plots'), recursive=T)

sums <- summary(fit, pars='tax_normFacts', probs=c(0.05,0.95))

sums3d <- array(NA, dim=c(NEffects,NEstTaxNodes,ncol(sums$summary)))
for(factor in 1:NEffects) {
    sums3d[factor,,] <- sums$summary[(factor-1) * NEstTaxNodes + (1:NEstTaxNodes),]
}
dimnames(sums3d) <- list(colnames(modelMat)[1:NEffects],estnodes,colnames(sums$summary))
factorfilenames <- gsub(':','.',colnames(modelMat))

for(effect in 1:NEffects) {
    cat('\t', file=paste0(outdir,'balanceTree_GLM_plots/',factorfilenames[[effect]],'.txt'))
    write.table(sums3d[effect,,], file=paste0(outdir,'balanceTree_GLM_plots/',factorfilenames[[effect]],'.txt'), sep='\t', quote=F,append=T)
}

#    factname <- names(which(factLevelMat[effect,]==1))

ext <- extract(fit, pars='tax_normFacts')[[1]]

for(i in sumconts) {
    hi <- ext[,factLevelMat[,i] == 1,]
    sup <- -apply(hi, 1, colSums)
    rownames(sup) <- estnodes
    sup2 <- sup
    dim(sup2) <- c(dim(sup),1)
    yeah <- monitor(aperm(sup2,c(2,3,1)), warmup=0, probs=c(0.05,0.95))
    rownames(yeah) <- estnodes
    cat('\t', file=paste0(outdir,'balanceTree_GLM_plots/',i,levels(newermap[,i])[nlevels(newermap[,i])],'.txt'))
    write.table(yeah, file=paste0(outdir,'balanceTree_GLM_plots/',i,levels(newermap[,i])[nlevels(newermap[,i])],'.txt'), sep='\t', quote=F,append=T)
}

save.image(paste0(outdir,'image.RData'))

sumsGlob <- summary(fit, pars='global', probs=c(0.05,0.95))
dimnames(sumsGlob$summary) <- list(colnames(modelMat)[1:NEffects], colnames(sumsGlob$summary))
cat('\t', file=paste0(outdir,'balanceTree_GLM_plots/globaleffects.txt'))
write.table(sumsGlob$summary, file=paste0(outdir,'balanceTree_GLM_plots/globaleffects.txt'), sep='\t', quote=F,append=T)

extScales <- extract(fit, pars='scales')[[1]]
vars1 <- t(apply(extScales, 1, function(x) x^2/sum(x^2)))
colnames(vars1) <- c(colnames(factLevelMat), paste0('global',colnames(factLevelMat)))

pdf(file=paste0(outdir,'balanceTree_GLM_plots/scalesboxes.pdf'))
boxplot(vars1, cex.axis=0.5, las=2)
graphics.off()

save(vars1,file=paste0(outdir,'balanceTree_GLM_plots/scales.RData'))

#bacttreeY.root$tip.label[Descendants(bacttreeY.root,157)[[1]]]
