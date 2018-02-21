library(coda)
library(rstan)
library(ape)
library(geiger)
library(phytools)
library(phangorn)
library(shinystan)
library(reshape2)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

indir <-  commandArgs(trailingOnly=TRUE)[[1]]
outdir <-  commandArgs(trailingOnly=TRUE)[[2]]

modelfile <- paste0(getwd(),'/balanceTree_sourceTracking.stan')
fulltable <- t(read.table(paste0(indir,'/otu_table_n500.txt'), header=T, sep='\t', comment.char='', skip=1, row.names=1)) # OTU-table
bacttree <- read.tree(paste0(indir,'/gg_constrained_fastttree.tre'))
tax_assignments <- paste0(indir,'/rep_set_tax_assignments.txt')
map <- read.table(paste0(indir,'/map_full_20171115.txt'),header=T,sep='\t',comment.char='',check.names=F)
rownames(map) <- map[,'#SampleID']


mincountnode <- 5
minsampsnode <- 1
mincountparent <- 25
minsampsbalance <- 3

##for testing purposes, trim heavily
mincountnode <- 5
minsampsnode <- 1
mincountparent <- 1500
minsampsbalance <- 7


## filter the data
newmap <- droplevels(map[map$date==20120816,])

newmap$nutrients <- relevel(newmap$nutrients, 'no')
newmap$dictyota_contact <- relevel(newmap$dictyota_contact, 'no')
newmap$sargassum_contact <- relevel(newmap$sargassum_contact, 'no')
newmap$halimeda_contact <- relevel(newmap$halimeda_contact, 'no')
newmap$plastic_contact <- relevel(newmap$plastic_contact, 'no')

newmap[newmap$Env != 'OFav_mucus' & newmap$Env != 'OFav_tissue', c('dictyota_contact','sargassum_contact','halimeda_contact','plastic_contact')] <- 'no'

newmap$nutrientsXdictyota[newmap$nutrients == 'yes' & newmap$dictyota_contact == 'yes'] <- 'yes'
newmap$nutrientsXdictyota[newmap$nutrients == 'no' | newmap$dictyota_contact == 'no'] <- 'no'
newmap$nutrientsXdictyota <- as.factor(newmap$nutrientsXdictyota)
newmap$nutrientsXsargassum[newmap$nutrients == 'yes' & newmap$sargassum_contact == 'yes'] <- 'yes'
newmap$nutrientsXsargassum[newmap$nutrients == 'no' | newmap$sargassum_contact == 'no'] <- 'no'
newmap$nutrientsXsargassum <- as.factor(newmap$nutrientsXsargassum)
newmap$nutrientsXhalimeda[newmap$nutrients == 'yes' & newmap$halimeda_contact == 'yes'] <- 'yes'
newmap$nutrientsXhalimeda[newmap$nutrients == 'no' | newmap$halimeda_contact == 'no'] <- 'no'
newmap$nutrientsXhalimeda <- as.factor(newmap$nutrientsXhalimeda)
newmap$nutrientsXplastic[newmap$nutrients == 'yes' & newmap$plastic_contact == 'yes'] <- 'yes'
newmap$nutrientsXplastic[newmap$nutrients == 'no' | newmap$plastic_contact == 'no'] <- 'no'
newmap$nutrientsXplastic <- as.factor(newmap$nutrientsXplastic)

modelform <- ~ nutrients + dictyota_contact + sargassum_contact + halimeda_contact + plastic_contact + nutrientsXdictyota + nutrientsXsargassum + nutrientsXhalimeda + nutrientsXplastic

idx <- rownames(fulltable)[rownames(fulltable) %in% rownames(newmap)]

y.old <- fulltable[idx,colnames(fulltable) %in% bacttree$tip.label]
newermap <- newmap[idx,]

y.filt <- y.old[,apply(y.old,2,function(x) sum(x >= mincountnode) >= minsampsnode)]
##will also want a minimum sample count. would be less important for basic analyses of the balances (i.e. a sample with only a count of 50 could still create a couple of observations to inform the top-level balances), but for mixtures of envs we will need more information?

bacttreeY <- drop.tip(bacttree,bacttree$tip.label[!bacttree$tip.label %in% colnames(y.filt)])

taxdat <- read.table(tax_assignments,sep='\t',stringsAsFactors=F, row.names=1)
x <- strsplit(taxdat[,1],'; ')
most <- max(sapply(x,length))
parsedtax <- lapply(x,function(x) {length(x) <- most; return(x)})
tax <- do.call('rbind',parsedtax)
rownames(tax) <- rownames(taxdat)
colnames(tax) <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')

archs <- rownames(tax)[tax[,'Kingdom']=='k__Archaea']


## root the tree if it's unrooted
if(is.rooted(bacttreeY)) {
    bacttreeY.root <- reorder(bacttreeY, order='pruningwise')
} else {
    bacttreeY.root <- reorder(root(bacttreeY, node=getMRCA(bacttreeY,archs[archs %in% bacttreeY$tip.label]), resolve.root = T), order='pruningwise')
    #bacttreeY.root <- midpoint.root(bacttreeY)
}
##

## summarize the putative taxa to be estimated
taxa <- factor(colnames(y.filt),levels=bacttreeY.root$tip.label)
taxnames <- levels(taxa)
NTaxa <- length(taxnames)
NIntTaxNodes <- bacttreeY.root$Nnode
NTaxNodes <- NTaxa + NIntTaxNodes

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
newestmap <- droplevels(newermap[levels(factor(senddat[,1])),])
##

envs <- newestmap$Env
envnames <- levels(envs)
NEnvs <- length(envnames)

## create the model matrix for location effects
modelMat <- model.matrix(modelform, data=newestmap)



## organize and summarize model parameters to pass to Stan
samplenames <- as.numeric(factor(senddat[,1]))
nodenames <- as.numeric(factor(senddat[,2], levels=estnodes))
sampledcounts <- as.numeric(senddat[,3])
datacounts <- as.numeric(senddat[,4])


NEstTaxNodes <- length(estnodes)
NObs <- length(sampledcounts)
NEstSamples <- length(unique(samplenames))
NFactors <- ncol(modelMat)

#####

nchains=4

init <- lapply(1:nchains, function(...) list(env_prop_normal_raw=matrix(-5.0,nrow=(NEnvs),ncol=NEnvs)))

fit <- stan(file=modelfile, data=list(NObs=NObs, NEstSamples=NEstSamples, sampledcounts=sampledcounts, datacounts=datacounts, samplenames=samplenames, nodenames=nodenames, NEstTaxNodes=NEstTaxNodes, NEnvs=NEnvs, envs=as.numeric(envs), NFactors=NFactors, modelMat=modelMat ), control=list(adapt_delta=0.8, max_treedepth=15), iter=1000, chains=nchains, thin=1, init_r=0.2, init=init )

save.image(file=paste0(outdir,'/st_out.RData'))


