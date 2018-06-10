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

outdir <- file.path('output',gsub(':','-',gsub(' ', '_', Sys.time())))

modelPath <- 'scripts/multinomial_sourceTracking.stan' #stan model
pathd8path <- 'scripts/PATHd8'

microbeTreePath <- 'raw_data/gg_constrained_fastttree.tre' # greengenes-constrained phylogeny of 97% OTUs
fulltablePath <- 'raw_data/otu_table_n500.txt' # open reference 97% OTUs, filtered to samples that had a minimum 500 reads
taxAssignmentsPath <- 'raw_data/rep_set_tax_assignments.txt' # greengenes taxonomy assignments
mapfilePath <- 'raw_data/map_full_20171115.txt' #mapping file


mincountnode <- 5
minsampsnode <- 1
mincountsamp <- 100


map <- read.table(mapfilePath, header=T, sep='\t', comment.char='', check.names=F)
rownames(map) <- map[,'#SampleID']


## filter the data
newmap <- droplevels(map[map$date==20120816,])

newmap$nutrients <- relevel(newmap$nutrients, 'no')
newmap$dictyota_contact <- relevel(newmap$dictyota_contact, 'no')
newmap$sargassum_contact <- relevel(newmap$sargassum_contact, 'no')
newmap$halimeda_contact <- relevel(newmap$halimeda_contact, 'no')
newmap$plastic_contact <- relevel(newmap$plastic_contact, 'no')

newmap[newmap$Env != 'OFav_mucus' & newmap$Env != 'OFav_tissue', c('dictyota_contact','sargassum_contact','halimeda_contact','plastic_contact')] <- 'no'

contrasts(newmap$Env) <- 'contr.sum'

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

fulltable <- t(read.table(fulltablePath, header=T, sep='\t', comment.char='', skip=1, row.names=1))

idx <- rownames(fulltable)[rownames(fulltable) %in% rownames(newmap)]

bacttree <- read.tree(microbeTreePath)

y.old <- fulltable[idx,colnames(fulltable) %in% bacttree$tip.label]
newermap <- newmap[idx,]

cat('Starting taxa number: ', ncol(y.old), '\n')
cat('Starting sample number: ', nrow(y.old), '\n')

y.taxfilt <- y.old[,apply(y.old,2,function(x) sum(x >= mincountnode) >= minsampsnode)]

cat('Taxa remaining after first filter: ', ncol(y.taxfilt), '\n')

y.sampsfilt <- y.taxfilt[apply(y.taxfilt, 1, function(x) sum(x) >= mincountsamp),]
##will also want a minimum sample count. would be less important for basic analyses of the balances (i.e. a sample with only a count of 50 could still create a couple of observations to inform the top-level balances), but for mixtures of envs we will need more information?

cat('Samples remaining after first filter: ', nrow(y.sampsfilt), '\n')

newestmap <- newermap[rownames(y.sampsfilt),]

y.filt <- y.sampsfilt[,apply(y.sampsfilt,2,function(x) sum(x >= mincountnode) >= minsampsnode)]

cat('Taxa remaining after second filter: ', ncol(y.filt), '\n')


bacttreeY <- drop.tip(bacttree,bacttree$tip.label[!bacttree$tip.label %in% colnames(y.filt)])

taxdat <- read.table(taxAssignmentsPath, sep='\t', stringsAsFactors=F, row.names=1)
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
    #bacttreeY.root <- reorder(root(bacttreeY, node=getMRCA(bacttreeY,archs[archs %in% bacttreeY$tip.label]), resolve.root = T), order='pruningwise') #archaea appear to be so polyphyletic that this is essentially random?
    bacttreeY.root <- midpoint.root(bacttreeY)
}
##

dir.create(outdir, recursive=T)

if(is.null(bacttreeY.root$edge.length)) {
    bacttreeY.root.chronos <- bacttreeY.root
    bacttreeY.root.chronos$edge.length <- rep(1,length(bacttreeY.root.chronos$edge))
    bacttreeY.root.chronos <- chronos(bacttreeY.root.chronos, model = "discrete", control = chronos.control(nb.rate.cat = 1))
    class(bacttreeY.root.chronos) <- 'phylo'
} else {
    bacttreeY.root.chronos <- bacttreeY.root
    bacttreeY.root.chronos$node.label <- NULL
    write.tree(bacttreeY.root.chronos, file=file.path(outdir,'bacttreeY.root.tree'), digits = 20)
    system(paste(pathd8path, file.path(outdir,'bacttreeY.root.tree'), file.path(outdir,'bacttreeY.root.pathd8.tree'), sep = " ")) ##taken from geiger's utility function PATHd8.phylo
    system(paste("grep \"d8 tree\" ", file.path(outdir,'bacttreeY.root.pathd8.tree'), ">", file.path(outdir,'bacttreeY.root.pathd8.parsed.tree'),
    sep = " "))
    bacttreeY.root.chronos = read.tree(file.path(outdir,'bacttreeY.root.pathd8.parsed.tree'))
}



## summarize the putative taxa to be estimated
taxnames <- levels(factor(colnames(y.filt),levels=bacttreeY.root.chronos$tip.label))
NMicrobeTips <- length(taxnames)
NIntTaxNodes <- bacttreeY.root.chronos$Nnode
NMicrobeNodes <- NMicrobeTips + NIntTaxNodes


microbeNodeScalesRaw <- sqrt(bacttreeY.root.chronos$edge.length)
names(microbeNodeScalesRaw) <- bacttreeY.root.chronos$edge[,2]
microbeNodeScalesRaw[as.character(length(bacttreeY.root.chronos$tip.label) + 1)] <- 1

microbeNodeScales <- microbeNodeScalesRaw[as.character(1:NMicrobeNodes)]


## sort the OTU table so its entries match the tree's tip labels
y <- y.filt[,taxnames]
##

microbeAncestors <- matrix(0, NMicrobeNodes, NMicrobeTips)
for(tip in 1:NMicrobeTips) {
    microbeAncestors[, tip] <- as.numeric(1:NMicrobeNodes %in% c(Ancestors(bacttreeY.root.chronos, tip), tip))
}

colnames(microbeAncestors) <- paste0('t',colnames(y))
rownames(microbeAncestors) <- paste0('i',1:NMicrobeNodes)
rownames(microbeAncestors)[1:NMicrobeTips] <- colnames(microbeAncestors)


## remove some potentially large objects to free up some memory
rm(fulltable,y.old,y.filt,bacttree,bacttreeY)
gc()
##


envs <- newestmap$Env
envnames <- levels(envs)
NEnvs <- length(envnames)
envFact <- cbind(as.numeric(envnames == 'OFav_mucus'), as.numeric(envnames == 'OFav_tissue'))
NSepFacts <- 2

## create the model matrix for location effects
modelMat <- model.matrix(modelform, data=newestmap)[,-1]


NSamples <- nrow(y)
NFactors <- ncol(modelMat)

envMat <- matrix(0,NSamples,NEnvs)
for (env in 1:NEnvs) {
    envMat[,env] <- as.numeric(newestmap$Env == envnames[[env]])
}

envEffects <- factor(c(envnames), levels=envnames)
contrasts(envEffects) <- 'contr.sum'
envTaxMat <- model.matrix(~env, data=list(env=envEffects))

NSampsPerEnv <- sapply(envnames, function(x) sum(envs==x))
isBaseSample <- as.numeric(sapply(1:length(envs), function(x) sum(envs[x]==envs[1:x]) == NSampsPerEnv[envs[x]])) ## MUST be last sample of a particular env!

#####

nchains=4


save.image(file=file.path(outdir, 'stMultinomial_vb_setup.RData'))


fit <- vb(stan_model(file=modelPath), data=list(y=t(y), NSamples=NSamples, NMicrobeNodes=NMicrobeNodes, NMicrobeTips=NMicrobeTips, NEnvs=NEnvs, envs=as.numeric(envs), envTaxMat=envTaxMat, NFactors=NFactors, modelMat=modelMat, envFact=envFact, NSepFacts=NSepFacts, microbeAncestors=microbeAncestors, taxAveStDPriorExpect=1.0, envPropAveStDPriorExpect=1.0, microbeNodeScales=microbeNodeScales ), iter=10000, pars=c('env_props','samp_props') ) #adapt_engaged=F, eta=0.1

save(fit, file=file.path(outdir, 'stMultinomial_vb_fit.RData'))

ext <- extract(fit, pars=c('env_props[1,1]','env_props[2,2]','env_props[3,3]','env_props[4,4]','env_props[5,5]','env_props[6,6]','env_props[7,7]','env_props[1,8]','env_props[4,8]','samp_props[7,4]','samp_props[4,7]'))

pdf(file=file.path(outdir, 'stMultinomial_vb_boxes.pdf'))
boxplot(ext)
graphics.off()

