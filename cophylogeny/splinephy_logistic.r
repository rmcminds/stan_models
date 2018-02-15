library(coda)
library(rstan)
library(ape)
library(geiger)
library(phytools)
library(phangorn)
library(shinystan)
library(nlme)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


################################## customize options in this section ###########################

## locate all your files
outdir <- 'splinephy_logistic_out/'
otu_table <- 'otutable.txt' ##tsv as produced by biom convert --to-tsv
tax_assignments <- 'reference-hit.seqs_tax_assignments.txt'
intree <- 'symbionttree.tre'
mapfile <- 'map.txt'
hosttreefile <- 'hosttree.newick'
modelfile <- 'splinephy_logistic.stan'
##



## set some thresholds for filtering your data
mincountnode <- 100 #minimum count of a node in a sample
minsampsnode <- 1 #minimum number of samples in which a node must have a count of at least 'mincountnode' to be kept
mincountsamp <- 50
#

## place custom code within this function to filter your mapping file to contain a subset of its samples and to recast variables into their proper formats (e.g. factors aren't continuous, base level is correct)
filterfunction <- function(dfin) {
    levels(dfin[,sampleTipKey]) <- gsub(' ','_',levels(dfin[,sampleTipKey]))
    df2 <- droplevels(dfin[dfin$tissue_compartment=='T',])
    df2[,sampleTipKey][df2[,sampleTipKey] %in% c('Missing:_Not_collected', 'unknown_unknown', 'unknown_species_A')] <- NA
    df2 <- droplevels(df2[df2[,sampleTipKey] %in% hosttree$tip.label,])
    return(df2)
}
##

sampleTipKey <- 'host_scientific_name'


## Stan options
nchains <- 4 ## should run at least 4 chains to make sure they converge on the same answer (they will be run in parallel)
iterations <- 2000 ## will probably need >10,000? maybe start with 2, check convergence, double it, check, double, check, double, etc.?
max_treedepth <- 15 ## a warning will tell you if this needs to be increased
adapt_delta <- 0.9 ## increase this if you get 'divergences' - even one means your model fit sucks!
##


######################### and then just run the following code (hopefully) ##########################

## load the data
fulltable <- t(read.table(otu_table, header=T, sep='\t', comment.char='', skip=1, row.names=1, check.names=F))
bacttree <- read.tree(intree)
map <- read.table(mapfile,header=T,sep='\t',comment.char='',check.names=F)
rownames(map) <- map[,'#SampleID']
hosttree <- read.tree(hosttreefile)

taxdat <- read.table(tax_assignments,sep='\t',stringsAsFactors=F, row.names=1)
x <- strsplit(taxdat[,1],'; ')
most <- max(sapply(x,length))
parsedtax <- lapply(x,function(x) {length(x) <- most; return(x)})
tax <- do.call('rbind',parsedtax)
rownames(tax) <- rownames(taxdat)
colnames(tax) <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')
##


## filter the data
newmap <- filterfunction(map)
idx <- rownames(fulltable)[rownames(fulltable) %in% rownames(newmap)]
subsetTaxa <- rownames(tax[tax[,'Family'] == 'f__Endozoicimonaceae',])
y.old <- fulltable[idx, colnames(fulltable) %in% bacttree$tip.label & colnames(fulltable) %in% subsetTaxa]
newermap <- newmap[idx,]
y.filt <- y.old[,apply(y.old,2,function(x) sum(x >= mincountnode) >= minsampsnode)]
y.filt2 <- y.filt[rowSums(y.filt) > mincountsamp,]
newestmap <- newermap[rownames(y.filt2),]
bacttreeY <- drop.tip(bacttree,bacttree$tip.label[!bacttree$tip.label %in% colnames(y.filt2)])
##

## root the tree if it's unrooted
if(is.rooted(bacttreeY)) {
    bacttreeY.root <- bacttreeY
} else {
	bacttreeY.root <- midpoint.root(bacttreeY)
}
##

symtreeheight <- max(nodeHeights(bacttreeY.root))
symNodeHeights <- apply(mrca(bacttreeY.root), c(1,2), function(x) nodeheight(bacttreeY.root,x)) / symtreeheight
symCophenetic <- cophenetic(bacttreeY.root) / symtreeheight

## summarize the putative taxa to be estimated
taxa <- factor(colnames(y.filt2),levels=bacttreeY.root$tip.label)
taxnames <- levels(taxa)
NTaxa <- length(taxnames)
NIntTaxNodes <- bacttreeY.root$Nnode
NTaxNodes <- NTaxa + NIntTaxNodes
##

## sort the OTU table so its entries match the tree's tip labels
y <- y.filt2[,taxnames]
##

## remove some potentially large objects to free up some memory
rm(fulltable,y.old,y.filt,bacttree,bacttreeY)
gc()
##


## parse the host phylogeny and create the phylogenetic distance matrices
hosts <- droplevels(newestmap[,sampleTipKey])
hosttreeY <- drop.tip(hosttree,hosttree$tip.label[!hosttree$tip.label %in% hosts])
hosts <- factor(hosts,levels=hosttreeY$tip.label)
sampleTips <- levels(hosts)
sampleTipMat <- sapply(sampleTips,function(x) as.numeric(hosts == x))

NTips <- ncol(sampleTipMat)

hosttreeheight <- max(nodeHeights(hosttreeY))
hostNodeHeights <- apply(mrca(hosttreeY), c(1,2), function(x) nodeheight(hosttreeY,x)) / hosttreeheight
hostCophenetic <- cophenetic(hosttreeY) / hosttreeheight
##hosttreeY <- drop.tip(hosttree,which(!hosttree$tip.label %in% sample(hosttree$tip.label,10)))
##ou <- function(coph,nh,alpha) {temp <- exp(-coph * alpha) * (1 - exp(-2 * nh * alpha)); return(diag(1/sqrt(diag(temp))) %*% temp %*% diag(1/sqrt(diag(temp))))}

## create the model matrix
modelMat <- sampleTipMat
##

NSym <- ncol(y)
NSamples <- nrow(y)
NHosts <- length(levels(hosts))
sampleCount <- rowSums(y)
presence <- y > 0
mode(presence) <- 'numeric'

NSymKnots_temp <- ceiling(1 + 3.322*(log10(NSym^2/2+NSym))) ##Sturge's rule
symCophLT <- symCophenetic[!lower.tri(symCophenetic)]
NSymCophLT <- length(symCophLT)
symKnots <- unique(round(unname(quantile(symCophLT,probs=seq(from=0, to=1, length.out = NSymKnots_temp))),7))
NSymKnots <- length(symKnots)
NHostKnots_temp <- ceiling(1 + 3.322*(log10(NHosts^2/2+NHosts))) ##Sturge's rule
hostCophLT <- hostCophenetic[!lower.tri(hostCophenetic)]
NHostCophLT <- length(hostCophLT)
hostKnots <- unique(round(unname(quantile(hostCophLT,probs=seq(from=0, to=1, length.out = NHostKnots_temp))),7))
NHostKnots <- length(hostKnots)
##http://mc-stan.org/users/documentation/case-studies/splines_in_stan.html
splineDegree <- 3

standat <- list(NSamples=NSamples, y=t(presence), NSym=NSym, NHosts=NHosts, modelMat=modelMat, hostCophenetic=hostCophenetic, symCophenetic=symCophenetic, NHostKnots=NHostKnots, NSymKnots=NSymKnots, symKnots=symKnots, hostKnots=hostKnots, splineDegree=splineDegree, hostCophLT=hostCophLT, symCophLT=symCophLT, NHostCophLT=NHostCophLT, NSymCophLT=NSymCophLT, sampleCount=sampleCount)
##


### run the model!
fit <- stan(file=modelfile, data=standat, control=list(adapt_delta=adapt_delta, max_treedepth=max_treedepth), iter=iterations, thin=floor(iterations/2000), chains=nchains )

dir.create(outdir)

save.image(file=paste0(outdir,'balanceTree_GLM_results.RData'))

shinystanres <- as.shinystan(fit, pars=c('tax_normFactsScaled','scaleFacts'))

save(shinystanres,file=paste0(outdir,'balanceTree_GLM_results_shinystan.RData'))



######################### plotting functions; maybe call these beta stage #######################

taxdat <- read.table(tax_assignments,sep='\t',stringsAsFactors=F, row.names=1)
x <- strsplit(taxdat[,1],'; ')
most <- max(sapply(x,length))
parsedtax <- lapply(x,function(x) {length(x) <- most; return(x)})
tax <- do.call('rbind',parsedtax)
rownames(tax) <- rownames(taxdat)
colnames(tax) <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')

obsrats <- datacounts/sampledcounts

plotdat <- data.frame(obsrats=obsrats,samplenames=samplenames,nodenames=nodenames)
for(fact in plotfacts) {
	plotdat[,fact] <- sapply(plotdat$samplenames,function(x) return(newestmap[x,fact]))
}

newtax <- tax[bacttreeY.root$tip.label[bacttreeY.root$tip.label %in% rownames(tax)],]

level <- colnames(newtax)[1:7]

dir.create(paste0(outdir,'balanceTree_GLM_plots/'),recursive = T)

for(i in unique(nodenames)) {
    pdf(file=paste0(outdir,'balanceTree_GLM_plots/',estnodes[[i]],'.pdf'),width=10)
    boxplot(update(modelform,obsrats~.),data=plotdat[plotdat$nodenames==i,])
    graphics.off()
    
    for(node in Children(bacttreeY.root, estnodes[[i]])) {
        tips <- bacttreeY.root$tip.label[Descendants(bacttreeY.root, node,type='tips')[[1]]]
        tips <- tips[tips %in% rownames(newtax)]
        cat(paste0('node ',node,':\n'),file=paste0(outdir,'balanceTree_GLM_plots/',estnodes[[i]],'.txt'),append=T)
        if(length(tips) != 1) {
            cat(paste0(sort(unique(apply(newtax[tips,level],1,paste,collapse='.'))),'\n'),file=paste0(outdir,'balanceTree_GLM_plots/',estnodes[[i]],'.txt'),append=T)
            cat(paste0(sort(unique(tips)),'\n'),file=paste0(outdir,'balanceTree_GLM_plots/',estnodes[[i]],'.txt'),append=T)
        } else {
            cat(paste0(paste(newtax[tips,level],collapse='.'),'\n'),file=paste0(outdir,'balanceTree_GLM_plots/',estnodes[[i]],'.txt'),append=T)
            cat(paste0(tips,'\n'),file=paste0(outdir,'balanceTree_GLM_plots/',estnodes[[i]],'.txt'),append=T)
        }
    }
    
}

