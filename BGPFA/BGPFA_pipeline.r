##want to modify such that script can be separated from actual data
# make utility functions to reshape input data? simplification of code is necessary anyway
library(rstan)
library(cmdstanr)
library(FactoMineR)
library(ape)
library(phangorn)
library(phytools)
library(MASS)
library(geosphere)
library(DESeq2)
utilities_dir <- file.path(Sys.getenv('HOME'), 'scripts/stan_models/utility/')
source(file.path(utilities_dir, 'read_stan_csv_subset.r'))
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
logit <- function(p) log(p/(1-p))
inv_logit <- function(x) { 1 / (1 + exp(-x)) }

nMicrobeKeep <- 500
K_linear <- 10
K_gp <- 20
KG <- 2
K <- K_linear + KG * K_gp
global_scale_prior = 2.5
rate_gamma_fact = 10
shape_gamma_fact = 2
site_smoothness <- 2
nu_residuals <- 25

input_prefix <- file.path(Sys.getenv('HOME'), 'data/tara_unsupervised_analyses')
if(exists('myargs')) {if(length(myargs)==1) {input_prefix <- myargs[[1]]}} else if(length(commandArgs(trailingOnly=TRUE)) > 0) {input_prefix <- commandArgs(trailingOnly=TRUE)[[1]]}
preprocess_prefix <- paste0(Sys.getenv('HOME'), '/outputs/tara/intermediate/')
model_dir <- file.path(Sys.getenv('HOME'), 'scripts/stan_models/BGPFA/')
model_name <- 'BGPFA'
engine <- 'advi'
opencl <- FALSE
output_prefix <- paste0(Sys.getenv('HOME'), '/outputs/tara/BGPFA_', nMicrobeKeep)

dir.create(output_prefix, recursive = TRUE)

sampling_commands <- list(sampling = paste(paste0('./', model_name),
                                          paste0('data file=', file.path(output_prefix, 'data.json')),
                                          paste0('init=', file.path(output_prefix, 'inits.json')),
                                          'output',
                                          paste0('file=', file.path(output_prefix, 'samples_hmc_quick.txt')),
                                          paste0('refresh=', 1),
                                          'method=sample algorithm=hmc',
                                          'stepsize=0.01',
                                          'engine=nuts',
                                          'max_depth=7',
                                          'num_warmup=300',
                                          'num_samples=100',
                                          ('opencl platform=0 device=1')[opencl],
                                          sep=' '),
                          advi = paste(paste0('./', model_name),
                                       paste0('data file=', file.path(output_prefix, 'data.json')),
                                       paste0('init=', file.path(output_prefix, 'inits.json')),
                                       'output',
                                       paste0('file=', file.path(output_prefix, 'samples_advi.txt')),
                                       paste0('refresh=', 100),
                                       'method=variational algorithm=meanfield',
                                       'grad_samples=1',
                                       'elbo_samples=100',
                                       'iter=30000',
                                       'eta=1',
                                       'adapt engaged=1',
                                       'tol_rel_obj=0.001',
                                       'output_samples=200',
                                       ('opencl platform=0 device=1')[opencl],
                                       sep=' '))

sample_data <- read.table(file.path(input_prefix, 'Stephane/20201008/TARA-PACIFIC_samples-provenance_20200731d.txt'), sep='\t', skip=1, header=T, comment.char='', stringsAsFactors=TRUE)
sample_data$species <- NA
sample_data$species[grepl('Pocillopora',sample_data$sample.material_taxonomy_host)] <-'Pocillopora'
sample_data$species[grepl('Porites',sample_data$sample.material_taxonomy_host)] <-'Porites'
sample_data$species[grepl('Millepora',sample_data$sample.material_taxonomy_host)] <-'Millepora'
sample_data$species[grepl('Heliopora',sample_data$sample.material_taxonomy_host)] <-'Heliopora'
sample_data <- sample_data[!is.na(sample_data$species),]
sample_data$myname <- sapply(as.character(sample_data$sampling.design_label),function(x) paste0(strsplit(x,'-')[[1]][2:4],collapse=''))

## import han/shini's ASV data (pre-processed by me)
load(file.path(preprocess_prefix, '20201231/TARA_PACIFIC_16S_asv_CO-0-filtered_combined_filtered.RData')) ##loads 'bacteria'
load(file.path(preprocess_prefix, '20201231/TARA_PACIFIC_16S_asv_CO-0-filtered_combined_filtered_disps.RData'))
load(file.path(preprocess_prefix, '20201231/TARA_PACIFIC_16S_asv_CO-0-filtered_combined_filtered_binary_disps.RData'))
load(file.path(preprocess_prefix, '20201231/TARA_PACIFIC_16S_asv_CO-0-filtered_combined_filtered_tax.RData'))

keepers <- c()
binar <- FALSE
abo <- names(shrink_list$ql_disp_shrunken[order(shrink_list$ql_disp_shrunken,decreasing=TRUE)])
bio <- names(binary_disps[order(binary_disps,decreasing=TRUE)])
for(i in 1:nMicrobeKeep) {
    if(!binar) {keepers <- c(keepers,abo[!abo %in% keepers][[1]])}
    else {keepers <- c(keepers,bio[!bio %in% keepers][[1]])}
    binar <- !binar
}
keepers <- sort(keepers)

bacteriaFilt <- t(bacteria[keepers,])
rownames(bacteriaFilt) <- sample_data[match(rownames(bacteriaFilt), sample_data$sample.storage_container.label), 'myname']
bacteriaFilt <-bacteriaFilt[order(rownames(bacteriaFilt)),]

bacttree <- read.tree(file.path(preprocess_prefix, '20201231/TARA_PACIFIC_16S_asv_CO-0-filtered_combined_filtered_aligned.tree'))
bacttreeY <- drop.tip(bacttree, bacttree$tip.label[!bacttree$tip.label %in% colnames(bacteriaFilt)])

#colnames(tax) <- c('Root','Domain','Phylum','Class','Order','Family','Genus')

archs <- rownames(tax)[tax[,2]=='Archaea' & !is.na(tax[,2])]
myarchs <- archs[archs %in% bacttreeY$tip.label]
## root the tree
if(length(myarchs) > 1) {
    bacttreeY.root <- reorder(root(bacttreeY, node = getMRCA(bacttreeY, myarchs), resolve.root = TRUE), order='pruningwise')
} else if(length(myarchs) == 1) {
    bacttreeY.root <- reorder(root(bacttreeY, outgroup = myarchs, resolve.root = TRUE), order='pruningwise')
} else {
    bacttreeY.root <- MidpointRooter:::midpoint.root2(bacttreeY)
}

bacteriaFilt <- bacteriaFilt[,bacttreeY.root$tip.label]

NTips <- length(bacttreeY.root$tip.label)
NNodes <- bacttreeY.root$Nnode + NTips - 1


## import Nicolas' ASV data (pre-processed by me)
load(file.path(preprocess_prefix, '20210102/TARA_PACIFIC_18SV9_4191_samples_v202004.OTU.filtered_CO-0-filtered_combined_filtered.RData')) ##loads 'bacteria'
load(file.path(preprocess_prefix, '20210102/TARA_PACIFIC_18SV9_4191_samples_v202004.OTU.filtered_CO-0-filtered_combined_filtered_disps.RData'))
load(file.path(preprocess_prefix, '20210102/TARA_PACIFIC_18SV9_4191_samples_v202004.OTU.filtered_CO-0-filtered_combined_filtered_binary_disps.RData'))
load(file.path(preprocess_prefix, '20210102/TARA_PACIFIC_18SV9_4191_samples_v202004.OTU.filtered_CO-0-filtered_combined_filtered_tax.RData'))

keepers_euk <- c()
binar <- FALSE
abo_euks <- names(shrink_list_euks$ql_disp_shrunken[order(shrink_list_euks$ql_disp_shrunken,decreasing=TRUE)])
bio_euks <- names(binary_disps_euks[order(binary_disps_euks,decreasing=TRUE)])
for(i in 1:nMicrobeKeep) {
    if(!binar) {keepers_euk <- c(keepers_euk, abo_euks[!abo_euks %in% keepers_euk][[1]])}
    else {keepers_euk <- c(keepers_euk, bio_euks[!bio_euks %in% keepers_euk][[1]])}
    binar <- !binar
}
keepers_euk <- sort(keepers_euk)

euksFilt <- t(euks[keepers_euk,])
rownames(euksFilt) <- sample_data[match(rownames(euksFilt), sample_data$sample.storage_container.label), 'myname']
euksFilt <-euksFilt[order(rownames(euksFilt)),]

euktree <- read.tree(file.path(preprocess_prefix, '20210102/TARA_PACIFIC_18SV9_4191_samples_v202004.OTU.filtered_CO-0-filtered_combined_filtered_aligned.tree'))
euktreeY <- drop.tip(euktree, euktree$tip.label[!euktree$tip.label %in% colnames(euksFilt)])
euktreeY.root <- MidpointRooter:::midpoint.root2(euktreeY)

euksFilt <- euksFilt[,euktreeY.root$tip.label]

NTipsEuks <- length(euktreeY.root$tip.label)
NNodesEuks <- euktreeY.root$Nnode + NTipsEuks - 1


###

transcr <- as.matrix(read.csv(file.path(input_prefix, 'Alice/20190924/transcriptomics/tara_I10_trans.csv'), sep=';', row.names=1, stringsAsFactors=TRUE)[,1:48])
rownames(transcr) <- paste(substr(rownames(transcr), 1, 7), 0, substr(rownames(transcr), 8, nchar(rownames(transcr))), sep = "")
transcr <- transcr[,colSums(transcr) > 0]
transcr <- t(apply(transcr,1,function(x) round(10 * x)))
transcr <- cbind(transcr, apply(transcr,1,function(x) 10000000 - sum(x))) #data are normalized out of a million reads and have accuracy of 1/100th. to treat them as counts i am going to assume all samples had exactly 10 million reads originally - because i suspect the siginificant digits are optimistic and this number was not actually sequenced. i'm hoping more were sequenced, but assuming a lower count makes them less confident and this is therefore conservative
colnames(transcr)[ncol(transcr)] <- 'X.Other'
colnames(transcr) <- sub('^X.','transcr_',colnames(transcr))

## import telomere data

t2Raw <- read.csv(file.path(input_prefix, 'Alice/20200527/T2_NOTNorm_Table_130520_final.csv'), sep=';', row.names=1, stringsAsFactors=TRUE,)
colnames(t2Raw)[1:5] <- sub('T2','T2_raw', colnames(t2Raw)[1:5])
rownames(t2Raw) <- sapply(t2Raw$ircan, function(x) paste0(c('I',substr(x,2,3),'S0',substr(x,4,4),'C0',substr(x,6,7)),collapse=''))

t2Norm <- read.csv(file.path(input_prefix, 'Alice/20200527/T2_Norm_Table_130520_final.csv'), sep=';', row.names=1, stringsAsFactors=TRUE,)
colnames(t2Norm)[1:5] <- sub('T2','T2_norm', colnames(t2Norm)[1:5])
rownames(t2Norm) <- sapply(t2Norm$ircan, function(x) paste0(c('I',substr(x,2,3),'S0',substr(x,4,4),'C0',substr(x,6,7)),collapse=''))

t2in <- merge(t2Raw[,2:4], t2Norm[,2:4], by=0)
rownames(t2in) <- t2in[,1]
t2in <- t2in[,-1]

t3Raw <- read.csv(file.path(input_prefix, 'Alice/20200527/T3_NOTNorm_Table_260520_final.csv'), sep=';', row.names=1, stringsAsFactors=TRUE,)
colnames(t3Raw)[1:5] <- sub('T3','T3_raw', colnames(t3Raw)[1:5])
t3Raw$myname <- sapply(t3Raw$ircan, function(x) paste0(c('I',substr(x,2,3),'S0',substr(x,4,4),'C0',substr(x,6,7)),collapse=''))
t3Raw <- t3Raw[!t3Raw$myname %in% t3Raw$myname[duplicated(t3Raw$myname)],]
rownames(t3Raw) <- sapply(t3Raw$ircan, function(x) paste0(c('I',substr(x,2,3),'S0',substr(x,4,4),'C0',substr(x,6,7)),collapse=''))

t3Norm <- read.csv(file.path(input_prefix, 'Alice/20200527/T3_Norm_Table_260520_final.csv'), sep=';', row.names=1, stringsAsFactors=TRUE,)
colnames(t3Norm)[1:5] <- sub('T3','T3_norm', colnames(t3Norm)[1:5])
t3Norm$myname <- sapply(t3Norm$ircan, function(x) paste0(c('I',substr(x,2,3),'S0',substr(x,4,4),'C0',substr(x,6,7)),collapse=''))
t3Norm <- t3Norm[!t3Norm$myname %in% t3Norm$myname[duplicated(t3Norm$myname)],]
rownames(t3Norm) <- sapply(t3Norm$ircan, function(x) paste0(c('I',substr(x,2,3),'S0',substr(x,4,4),'C0',substr(x,6,7)),collapse=''))

t3in <- merge(t3Raw[,2:4], t3Norm[,2:4], by=0)
rownames(t3in) <- t3in[,1]
t3in <- t3in[,-1]

## import ITS2 profiles table
its <- as.matrix(read.table(file=file.path(input_prefix, 'Ben/20200212/its2_type_profiles/90_20200125_DBV_2020-01-28_05-39-25.985991.profiles.absolute.abund_only_named.txt'), row.names=1, header=T, sep='\t', comment.char='', stringsAsFactors=TRUE))
rownames(its) <- sub('TARA_','',rownames(its),fixed=T)
colnames(its) <- sub('^X','ITS2_',colnames(its))
its <- its[grepl('CO-', rownames(its)), ]
its <- its[,apply(its,2,function(x) sum(x>0)) > 0]
itsFilt <- its
rownames(itsFilt) <- sample_data[match(rownames(itsFilt), sample_data$sample.storage_container.label), 'myname']

itsMeta <- read.table(file=file.path(input_prefix, 'Ben/20200212/its2_type_profiles/90_20200125_DBV_2020-01-28_05-39-25.985991.profiles.absolute.meta_only.txt'), row.names=1, header=T, sep='\t', comment.char='', stringsAsFactors=TRUE)
rownames(itsMeta) <- paste0('ITS2_', rownames(itsMeta))

## import paola's biomarker data
biomarkers <- read.table(file.path(input_prefix, 'Paola/20200218/paola_readable.txt'), sep='\t', header=T, row.names=1, stringsAsFactors=TRUE)
rownames(biomarkers) <- gsub('-','',sub('OA000-','',rownames(biomarkers)))

snps <- read.csv(file.path(input_prefix, 'Didier/geneAlex_LD02.csv'), skip=2, row.names=2, stringsAsFactors=TRUE)[,-1]

rownames(snps) <- gsub('o',0,sapply(rownames(snps), function(x) paste0(c('I', substr(x,2,3), 'S', substr(x,5,6), 'C', substr(x,8,10)), collapse='')))

nSNPs <- ncol(snps) / 2

snpmat <- matrix(0, nrow = nrow(snps), ncol = 4 * nSNPs)
colnames(snpmat) <- 1:(4 * nSNPs)
rownames(snpmat) <- rownames(snps)
Msnp <- vector('numeric')
newsnpnames <- vector()
snpFilter <- vector()
indSNP <- 1
indSNPMat <- 1
for(i in 1:nSNPs) {
    snpmat[, indSNPMat:(indSNPMat + 3)] <- t(apply(snps[, indSNP:(indSNP+1)], 1, function(x) {
        sapply(1:4, function(y) sum(y == x))
    }))
    colnames(snpmat)[indSNPMat:(indSNPMat + 3)] <- paste(colnames(snps)[indSNP], 1:4, sep='.')
    newFilt <- apply(snpmat[, indSNPMat:(indSNPMat + 3)], 2, function(x) any(x > 0) & sd(x) > 0)
    if(sum(newFilt) == 2) {
        snpPoss <- sort(unique(snpmat[,(indSNPMat:(indSNPMat + 3))[which(newFilt)[[2]]]]))
        if(length(snpPoss) == 2) {

            poss1 <- snpmat[,(indSNPMat:(indSNPMat + 3))[which(newFilt)[[1]]]] == snpPoss[[1]]
            poss2 <- snpmat[,(indSNPMat:(indSNPMat + 3))[which(newFilt)[[1]]]] == snpPoss[[2]]

            if(sum(poss1) > 1 & sum(poss2) > 1) {
                newsnpnames <- c(newsnpnames, paste(names(newFilt)[which(newFilt)[[2]]], names(newFilt)[which(newFilt)[[1]]], sep='_vs_'))
                newFilt[which(newFilt)[[1]]] <- FALSE

                snpmat[poss1,(indSNPMat:(indSNPMat + 3))[newFilt]] <- 0
                snpmat[poss2,(indSNPMat:(indSNPMat + 3))[newFilt]] <- 1

                Msnp <- c(Msnp, 1)
            } else {
                newFilt[newFilt] <- FALSE
            }
        }
        else if(length(snpPoss) > 2) {

            het <- snpmat[,(indSNPMat:(indSNPMat + 3))[which(newFilt)[[1]]]] == 1
            hom1 <- snpmat[,(indSNPMat:(indSNPMat + 3))[which(newFilt)[[1]]]] == 2
            hom2 <- snpmat[,(indSNPMat:(indSNPMat + 3))[which(newFilt)[[1]]]] == 0

            homname <- paste(names(newFilt)[which(newFilt)[[2]]], names(newFilt)[which(newFilt)[[1]]], sep='_vs_')
            snpmat[het, (indSNPMat:(indSNPMat + 3))[which(newFilt)[[1]]]] <- NA
            snpmat[hom1, (indSNPMat:(indSNPMat + 3))[which(newFilt)[[1]]]] <- 0
            snpmat[hom2, (indSNPMat:(indSNPMat + 3))[which(newFilt)[[1]]]] <- 1

            hetname <- paste(names(newFilt)[which(newFilt)[[2]]], names(newFilt)[which(newFilt)[[1]]], sep='_AND_')
            snpmat[het, (indSNPMat:(indSNPMat + 3))[which(newFilt)[[2]]]] <- 1
            snpmat[hom1, (indSNPMat:(indSNPMat + 3))[which(newFilt)[[2]]]] <- 0
            snpmat[hom2, (indSNPMat:(indSNPMat + 3))[which(newFilt)[[2]]]] <- 0

            if(sum(hom1) > 1 & sum(hom2) > 1) {
                newsnpnames <- c(newsnpnames, homname)
                Msnp <- c(Msnp, 1)
            } else {
                newFilt[which(newFilt)[[1]]] <- FALSE
            }

            if(sum(het) > 1 & sum(!het) > 1) {
                newsnpnames <- c(newsnpnames, hetname)
                Msnp <- c(Msnp, 1)
            } else {
                newFilt[which(newFilt)[[2]]] <- FALSE
            }
        }
    } else if(sum(newFilt) > 2){
        newsnpnames <- c(newsnpnames, paste(colnames(snps)[indSNP], 1:4, sep='.'))
        Msnp <- c(Msnp, sum(newFilt))
    }
    snpFilter <- c(snpFilter, newFilt)
    indSNP <- indSNP + 2
    indSNPMat <- indSNPMat + 4
}
snpmat <- snpmat[,snpFilter]
colnames(snpmat) <- newsnpnames

### photos
photodat <- read.table(file.path(input_prefix,'Ryan/20200406/results_202004061414.txt'), header=T, sep='\t', quote = '', stringsAsFactors=FALSE)
photodat <- photodat[!is.na(photodat[,1]) & !is.na(photodat[,2]) & !is.na(photodat[,3]),]
photodat <- photodat[as.logical(photodat[,2]) & as.logical(photodat[,3]),]
rownames(photodat) <- paste0(rownames(photodat),'_',sapply(photodat[,1],function(x) strsplit(x,split='_',fixed=TRUE)[[1]][[4]]))

envvars <- c('bleached','parrotfish_scars','bivalve','Spirobranchus','Tridacna','boringurchin','other_polychaete','sponge','Halimeda','Turbinaria','Dictyota','Lobophora','cca','Galaxaura','Sargassum','unhealthy','turf_over','sediment','pigmentation','trematodiasis','gallcrabs','ascidians')

photonames <- gsub('-','',sub('^.*_OA000-','',rownames(photodat)))
photosamples <- unique(photonames)

photodatagg <- t(sapply(photosamples, function(x) {
    photosub <- photodat[photonames == x, c('porites_lumps','porites_ridges',envvars[envvars != 'bleached'])]
    if(sum(photonames == x) > 1) {
        apply(photosub,
              2,
              function(y) as.numeric(if(all(is.na(as.logical(y)))) NA else any(as.logical(y),na.rm=TRUE)))
    } else {
        return(as.numeric(as.logical(photosub)))
    }
}))

## host
allgenera <- unique(photodat$target)
ngen <- length(allgenera)

allmorphs <- sort(unique(paste(photodat$target,photodat$morphotype,sep='_'))[unique(photodat$morphotype) != ''])
nlevs <- length(allmorphs)
morphgroups <- sapply(allmorphs,function(x) strsplit(x,'_')[[1]][[1]])

hostphotomat <- sapply(allmorphs, function(y) y == paste(photodat$target,photodat$morphotype,sep='_'))
colnames(hostphotomat) <- paste('morphotype', allmorphs, sep='_')
Mhostphoto <- nlevs

hostphotomat2 <- t(sapply(photosamples, function(x) if(sum(photonames == x) > 1) {apply(hostphotomat[photonames == x,],2,function(y) sum(y,na.rm=TRUE)/length(y))} else {hostphotomat[photonames == x,]}))

hostphotomat2 <- cbind(hostphotomat2, photodatagg[rownames(hostphotomat2),c('porites_lumps','porites_ridges')])
Mhostphoto <- c(Mhostphoto, 1, 1)

## env
bleachLevs <- unique(photodat[,'bleached'])[!is.na(unique(photodat[,'bleached'])) & unique(photodat[,'bleached']) != '']
Menvphoto <- length(bleachLevs)

envphotomat <- sapply(bleachLevs, function(y) y == photodat$bleached)
colnames(envphotomat) <- paste('bleached', bleachLevs, sep='_')
envphotomat2 <- t(sapply(photosamples, function(x) if(sum(photonames == x) > 1) {apply(envphotomat[photonames == x,],2,function(y) sum(y,na.rm=TRUE)/length(y))} else {envphotomat[photonames == x,]}))

envphotomat2 <- cbind(envphotomat2, photodatagg[,envvars[envvars != 'bleached']])
envphotomat2 <- envphotomat2[,apply(envphotomat2,2,sd,na.rm=TRUE) > 0]
Menvphoto <- c(Menvphoto, rep(1,length(envvars[envvars %in% colnames(envphotomat2,2)])))

########

repeatedsamples <- c(rownames(bacteriaFilt), rownames(euksFilt), rownames(transcr), rownames(itsFilt), rownames(biomarkers), rownames(t2in), rownames(t3in), rownames(snpmat), rownames(hostphotomat2), rownames(envphotomat2))
allsamples <- unique(repeatedsamples)[sapply(unique(repeatedsamples), function(x) sum(x==repeatedsamples)>1)]

bacteriaFilt <- bacteriaFilt[allsamples[allsamples %in% rownames(bacteriaFilt)],]
euksFilt <- euksFilt[allsamples[allsamples %in% rownames(euksFilt)],]
transcr <- transcr[allsamples[allsamples %in% rownames(transcr)],]
itsFilt <- itsFilt[allsamples[allsamples %in% rownames(itsFilt)],]
itsFilt <- itsFilt[,apply(itsFilt,2,function(x) sum(x>0)) > 0]

t2log <- as.matrix(log(t2in[allsamples[allsamples %in% rownames(t2in)],]))
t3log <- as.matrix(log(t3in[allsamples[allsamples %in% rownames(t3in)],]))

biomarkersLog <- log(biomarkers[allsamples[allsamples %in% rownames(biomarkers)],])
biomarkersLog_inits <- biomarkersLog
for(i in 1:ncol(biomarkersLog_inits)){
    biomarkersLog_inits[is.infinite(biomarkersLog_inits[,i]),i] <- min(biomarkersLog_inits[!is.infinite(biomarkersLog_inits[,i]),i], na.rm=TRUE)
}
biomarkersLog_inits <- as.matrix(biomarkersLog_inits)

snpmat <- snpmat[allsamples[allsamples %in% rownames(snpmat)],]

hostphotomat2 <- hostphotomat2[allsamples[allsamples %in% rownames(hostphotomat2)],]
keepcol <- rep(TRUE,ncol(hostphotomat2))
MFilt <- Mhostphoto
for(u in 1:length(Mhostphoto)) {
    for(p in 1:Mhostphoto[u]) {
        if(sum(hostphotomat2[,sum(Mhostphoto[0:(u-1)])+p], na.rm=TRUE) == 0) {
            keepcol[sum(Mhostphoto[0:(u-1)])+p] <- FALSE
            MFilt[u] <- MFilt[u] - 1
        }
    }
}
hostphotomat2 <- hostphotomat2[,keepcol]
Mhostphoto <- MFilt

envphotomat2 <- envphotomat2[allsamples[allsamples %in% rownames(envphotomat2)],]
keepcol <- rep(TRUE,ncol(envphotomat2))
MFilt <- Menvphoto
for(u in 1:length(Menvphoto)) {
    for(p in 1:Menvphoto[u]) {
        if(sum(envphotomat2[,sum(Menvphoto[0:(u-1)])+p], na.rm=TRUE) == 0) {
            keepcol[sum(Menvphoto[0:(u-1)])+p] <- FALSE
            MFilt[u] <- MFilt[u] - 1
        }
    }
}
envphotomat2 <- envphotomat2[,keepcol]
Menvphoto <- MFilt


filtData <- sample_data[match(allsamples[allsamples %in% sample_data$myname], sample_data$myname),]
filtData <- droplevels(filtData)
rownames(filtData) <- filtData$myname
filtData$site <- as.factor(substr(filtData$myname,1,6))
filtData$island <- as.factor(substr(filtData$myname,1,3))

mmSpec <- model.matrix(~0 + species, data = filtData)
rownames(mmSpec) <- filtData$myname[!is.na(filtData$species)]
mmSpec <- mmSpec[allsamples[allsamples %in% rownames(mmSpec)],]
mSpec <- ncol(mmSpec)


### didier's pop gen analysis

snpSVDs <- read.table(file.path(input_prefix,'Didier/20200419/snp_clusters.txt'), header=T, sep='\t', quote = '', row.names = 1, stringsAsFactors=FALSE)
snpSVDmat <- model.matrix(~0 + clade, data = snpSVDs)
snpSVDmat <- snpSVDmat[allsamples[allsamples %in% rownames(snpSVDmat)],]
MSVD <- ncol(snpSVDmat)

###

fabT1 <- read.table(file.path(input_prefix, 'Fabien/20200127/TotalTable_meteo_nav_tsg_acs_par_cdom.txt'), sep='\t', header=T, comment.char='', stringsAsFactors=TRUE)
sites <- substr(fabT1$Var1, 21,27)
keep <- sapply(sites, function(x) strsplit(x,'')[[1]][[1]] == 'I' & (strsplit(x,'')[[1]][[6]] != '0' | strsplit(x,'')[[1]][[7]] != '0')) & (grepl('3x10', fabT1$Var1) | grepl('3X10', fabT1$Var1))
fabT1 <- fabT1[,10:ncol(fabT1)]
fabT1 <- fabT1[,!grepl('lat|lon|dire|unixtime|current_speed|Course_Over_Ground|Speed_over_Ground|apparent',colnames(fabT1))]
fabT1 <- as.data.frame(apply(fabT1,2,function(x) if(any(x<=0,na.rm=T)) x else log(x)))

fabT1sites <- unique(sites[keep])

fabT1agg <- array(0, dim = c(length(fabT1sites), ncol(fabT1)), dimnames = list(fabT1sites, colnames(fabT1)))
for(f in fabT1sites) {
    fabT1agg[f,] <- apply(fabT1[sites==f & keep,], 2, mean, na.rm = TRUE)
}
repeats <- substr(grep('Q50',colnames(fabT1agg),ignore.case = TRUE,value=TRUE),4,100)
fabT1agg <- fabT1agg[,!colnames(fabT1agg) %in% c(repeats, paste0('Q25',repeats), paste0('Q75',repeats))]
fabT1agg <- fabT1agg[,apply(fabT1agg, 2, function(x) {temp <- sd(x,na.rm=T); temp > 0 & !is.na(temp)})]

fabT2 <- read.table(file.path(input_prefix, 'Fabien/20190329/TotalTable_Satellites_2S.txt'), sep='\t', header=T, comment.char='', stringsAsFactors=TRUE)
islands <- substr(fabT2$st_ID, 1,3)
keep <- sapply(islands, function(x) strsplit(x,'')[[1]][[2]] != '0' | strsplit(x,'')[[1]][[2]] != '0')
fabT2 <- fabT2[,25:(ncol(fabT2)-1)]
fabT2 <- fabT2[,apply(fabT2, 2, function(x) {temp <- sd(x,na.rm=T); temp > 0 & !is.na(temp)})]
fabT2 <- log(fabT2[,!grepl('_SD',colnames(fabT2))])

fabT2isls <- unique(islands[keep])
fabT2agg <- array(0, dim = c(length(fabT2isls), ncol(fabT2)), dimnames = list(fabT2isls, colnames(fabT2)))
for(f in fabT2isls) {
    fabT2agg[f,] <- apply(fabT2[islands==f & keep,], 2, mean, na.rm = TRUE)
}


allVarGroups <- unique(c(fabT1sites,fabT2isls))
nVarGroups <- length(allVarGroups)
nVarsWGroups <- ncol(fabT1agg) + ncol(fabT2agg)

fabT1agg <- fabT1agg[allVarGroups[allVarGroups %in% rownames(fabT1agg)],]
fabT2agg <- fabT2agg[allVarGroups[allVarGroups %in% rownames(fabT2agg)],]


### site matrix (special treatment in model)
longlat <- matrix(nrow=0,ncol=2)
for(i in levels(filtData$site)) {
    lo <- as.numeric(as.character(filtData$sampling.event_longitude_start_ddd.dddddd[filtData$site == i][[1]]))
    la <- as.numeric(as.character(filtData$sampling.event_latitude_start_dd.dddddd[filtData$site == i][[1]]))
    longlat <- rbind(longlat,c(lo,la))
    rownames(longlat)[nrow(longlat)] <- i
}
dupsites <- rownames(longlat)[duplicated(longlat)]
for(i in dupsites) {
    matchsites <- apply(longlat,1,function(x) all(x == longlat[i,]))
    standardsite <-rownames(longlat)[matchsites][[1]]
    levels(filtData$site)[levels(filtData$site) %in% rownames(longlat)[matchsites]] <- standardsite
}
longlat <- longlat[levels(filtData$site),]

missinglonglat <- apply(longlat,1,function(x) any(is.na(x)))
dist_sites <- geosphere:::distm(longlat[!missinglonglat,], fun = function(x,y) geosphere:::distVincentySphere(x,y,1))
for(i in which(missinglonglat)) {
    dist_sites <- rbind(dist_sites[1:(i-1),],1,dist_sites[i:nrow(dist_sites),])
    dist_sites <- cbind(dist_sites[,1:(i-1)],1,dist_sites[,i:ncol(dist_sites)])
    dist_sites[i,i] <- 0
}

siteMat <- model.matrix(~0+site, data=filtData)
siteMat <- siteMat[allsamples[allsamples %in% rownames(siteMat)],]
nSites <- ncol(siteMat)
##

N <- length(allsamples)
D <- 4
R <- 5
C <- 6
M <- c(ncol(bacteriaFilt), ncol(euksFilt), ncol(transcr), ncol(itsFilt), ncol(biomarkersLog), ncol(t2log), ncol(t3log), ncol(fabT1agg), ncol(fabT2agg), ncol(snpmat), ncol(hostphotomat2), ncol(envphotomat2), ncol(mmSpec), ncol(snpSVDmat), ncol(siteMat))
VOB <- sum(M)
samp2group <- cbind(sapply(fabT1sites, function(x) as.integer(x == filtData[allsamples,'site'])), sapply(fabT2isls, function(x) as.integer(x == filtData[allsamples,'island'])))
samp2group <- samp2group %*% diag(1.0 / colSums(samp2group, na.rm=TRUE))
samp2group[is.na(samp2group)] <- 0
Mc <- c(Msnp, Mhostphoto, Menvphoto, mSpec, MSVD, nSites)
C_vars <- length(Mc)
I <- matrix(0, nrow=D, ncol=(N+nVarGroups))
I[1,1:N] <- as.integer(allsamples %in% rownames(bacteriaFilt))
I[2,1:N] <- as.integer(allsamples %in% rownames(euksFilt))
I[3,1:N] <- as.integer(allsamples %in% rownames(transcr))
I[4,1:N] <- as.integer(allsamples %in% rownames(itsFilt))
colnames(I) <- c(allsamples, paste0('varGroup',1:nVarGroups))

O <- ncol(biomarkersLog) + ncol(t2log) + ncol(t3log) + ncol(fabT1agg) + ncol(fabT2agg)
IR <- matrix(0, nrow=O, ncol=(N+nVarGroups)) #analogous to I above but per-variable rather than per-dataset
IR[1:ncol(biomarkersLog),1:N] <- t(apply(biomarkersLog, 2, function(x) as.integer(allsamples %in% rownames(biomarkersLog)[!is.na(x)])))
IR[ncol(biomarkersLog) + (1:ncol(t2log)),1:N] <- t(apply(t2log, 2, function(x) as.integer(allsamples %in% rownames(t2log)[!is.na(x)])))
IR[ncol(biomarkersLog) + ncol(t2log) + (1:ncol(t3log)),1:N] <- t(apply(t3log, 2, function(x) as.integer(allsamples %in% rownames(t3log)[!is.na(x)])))
IR[ncol(biomarkersLog) + ncol(t2log) + ncol(t3log) + (1:ncol(fabT1agg)),(N+1):(N+nVarGroups)] <- t(apply(fabT1agg, 2, function(x) as.integer(allVarGroups %in% rownames(fabT1agg)[!is.na(x)])))
IR[ncol(biomarkersLog) + ncol(t2log) + ncol(t3log) + ncol(fabT1agg) + (1:ncol(fabT2agg)),(N+1):(N+nVarGroups)] <- t(apply(fabT2agg, 2, function(x) as.integer(allVarGroups %in% rownames(fabT2agg)[!is.na(x)])))

IC <- matrix(0, nrow=C_vars, ncol=(N+nVarGroups))
for(var in 1:length(Msnp)) {
    if(Msnp[var] > 1) {
        IC[var,1:N] <- as.integer(allsamples %in% rownames(snpmat)[apply(!is.na(snpmat[,(sum(Msnp[0:(var-1)])+1):sum(Msnp[0:var])]), 1, all)])
    } else {
        IC[var,1:N] <- as.integer(allsamples %in% rownames(snpmat)[!is.na(snpmat[,sum(Msnp[0:var])])])
    }
}
for(var in 1:length(Mhostphoto)) {
    if(Mhostphoto[var] > 1) {
        IC[length(Msnp) + var,1:N] <- as.integer(allsamples %in% rownames(hostphotomat2)[apply(!is.na(hostphotomat2[,(sum(Mhostphoto[0:(var-1)])+1):sum(Mhostphoto[0:var])]), 1, all)])
    } else {
        IC[length(Msnp) + var,1:N] <- as.integer(allsamples %in% rownames(hostphotomat2)[!is.na(hostphotomat2[,sum(Mhostphoto[0:var])])])
    }
}
for(var in 1:length(Menvphoto)) {
    if(Menvphoto[var] > 1) {
        IC[length(Msnp) + length(Mhostphoto) + var,1:N] <- as.integer(allsamples %in% rownames(envphotomat2)[apply(!is.na(envphotomat2[,(sum(Menvphoto[0:(var-1)])+1):sum(Menvphoto[0:var])]), 1, all)])
    } else {
        IC[length(Msnp) + length(Mhostphoto) + var,1:N] <- as.integer(allsamples %in% rownames(envphotomat2)[!is.na(envphotomat2[,sum(Menvphoto[0:var])])])
    }
}
for(var in 1:length(mSpec)) {
    if(mSpec[var] > 1) {
        IC[length(Msnp) + length(Mhostphoto) + length(Menvphoto) + var,1:N] <- as.integer(allsamples %in% rownames(mmSpec)[apply(!is.na(mmSpec[,(sum(mSpec[0:(var-1)])+1):sum(mSpec[0:var])]), 1, all)])
    } else {
        IC[length(Msnp) + length(Mhostphoto) + length(Menvphoto) + var,1:N] <- as.integer(allsamples %in% rownames(mmSpec)[!is.na(mmSpec[,sum(mSpec[0:var])])])
    }
}
for(var in 1:length(MSVD)) {
    if(MSVD[var] > 1) {
        IC[length(Msnp) + length(Mhostphoto) + length(Menvphoto) + length(mSpec) + var,1:N] <- as.integer(allsamples %in% rownames(snpSVDmat)[apply(!is.na(snpSVDmat[,(sum(MSVD[0:(var-1)])+1):sum(MSVD[0:var])]), 1, all)])
    } else {
        IC[length(Msnp) + length(Mhostphoto) + length(Menvphoto) + length(mSpec) + var,1:N] <- as.integer(allsamples %in% rownames(snpSVDmat)[!is.na(snpSVDmat[,sum(MSVD[0:var])])])
    }
}
for(var in 1:length(nSites)) {
    if(nSites[var] > 1) {
        IC[length(Msnp) + length(Mhostphoto) + length(Menvphoto) + length(mSpec) + length(MSVD) + var,1:N] <- as.integer(allsamples %in% rownames(siteMat)[apply(!is.na(siteMat[,(sum(nSites[0:(var-1)])+1):sum(nSites[0:var])]), 1, all)])
    } else {
        IC[length(Msnp) + length(Mhostphoto) + length(Menvphoto) + length(mSpec) + length(MSVD) + var,1:N] <- as.integer(allsamples %in% rownames(siteMat)[!is.na(siteMat[,sum(nSites[0:var])])])
    }
}
colnames(IC) <- c(allsamples, paste0('varGroup',1:nVarGroups))

I_cs <- t(apply(I, 1, cumsum))
IR_cs <- t(apply(IR, 1, cumsum))
IC_cs <- t(apply(IC, 1, cumsum))

X <- unlist(c(sapply(1:(N+nVarGroups), function(x) if(I[1,x]) bacteriaFilt[I_cs[1,x],]),
              sapply(1:(N+nVarGroups), function(x) if(I[2,x]) euksFilt[I_cs[2,x],]),
              sapply(1:(N+nVarGroups), function(x) if(I[3,x]) transcr[I_cs[3,x],]),
              sapply(1:(N+nVarGroups), function(x) if(I[4,x]) itsFilt[I_cs[4,x],])))

P <- unlist(c(sapply(1:(N+nVarGroups), function(x) unlist(sapply(1:ncol(biomarkersLog), function(y) {
                                                          temp <- biomarkersLog[!is.na(biomarkersLog[,y]),y]
                                                          names(temp) <- rep(colnames(biomarkersLog)[y], length(temp))
                                                          if(IR[y,x]) temp[IR_cs[y,x]]
                                                   }))),
              sapply(1:(N+nVarGroups), function(x) unlist(sapply(1:ncol(t2log), function(y) {
                                                          temp <- t2log[!is.na(t2log[,y]),y]
                                                          names(temp) <- rep(colnames(t2log)[y], length(temp))
                                                          if(IR[ncol(biomarkersLog) + y,x]) temp[IR_cs[ncol(biomarkersLog) + y,x]]
                                                   }))),
              sapply(1:(N+nVarGroups), function(x) unlist(sapply(1:ncol(t3log), function(y) {
                                                          temp <- t3log[!is.na(t3log[,y]),y]
                                                          names(temp) <- rep(colnames(t3log)[y], length(temp))
                                                          if(IR[ncol(biomarkersLog) + ncol(t2log) + y,x]) temp[IR_cs[ncol(biomarkersLog) + ncol(t2log) + y,x]]
                                                   }))),
              sapply(1:(N+nVarGroups), function(x) unlist(sapply(1:ncol(fabT1agg), function(y) {
                                                          temp <- fabT1agg[!is.na(fabT1agg[,y]),y]
                                                          names(temp) <- rep(colnames(fabT1agg)[y], length(temp))
                                                          if(IR[ncol(biomarkersLog) + ncol(t2log) + ncol(t3log) + y,x]) temp[IR_cs[ncol(biomarkersLog) + ncol(t2log) + ncol(t3log) + y,x]]
                                                   }))),
              sapply(1:(N+nVarGroups), function(x) unlist(sapply(1:ncol(fabT2agg), function(y) {
                                                          temp <- fabT2agg[!is.na(fabT2agg[,y]),y]
                                                          names(temp) <- rep(colnames(fabT2agg)[y], length(temp))
                                                          if(IR[ncol(biomarkersLog) + ncol(t2log) + ncol(t3log) + ncol(fabT1agg) + y,x]) temp[IR_cs[ncol(biomarkersLog) + ncol(t2log) + ncol(t3log) + ncol(fabT1agg) + y,x]]
                                                   })))))

Y <- unlist(c(sapply(1:length(Msnp), function(var) unlist(sapply(1:(N+nVarGroups), function(x) {
                                                          if(IC[var, x]) return(snpmat[allsamples[[x]], (sum(Msnp[0:(var-1)])+1):sum(Msnp[0:var])])
                                                   }))),
              sapply(1:length(Mhostphoto), function(var) unlist(sapply(1:(N+nVarGroups), function(x) {
                                                          if(IC[length(Msnp) + var, x]) return(hostphotomat2[allsamples[[x]], (sum(Mhostphoto[0:(var-1)])+1):sum(Mhostphoto[0:var])])
                                                   }))),
              sapply(1:length(Menvphoto), function(var) unlist(sapply(1:(N+nVarGroups), function(x) {
                                                          if(IC[length(Msnp) + length(Mhostphoto) + var, x]) return(envphotomat2[allsamples[[x]], (sum(Menvphoto[0:(var-1)])+1):sum(Menvphoto[0:var])])
                                                   }))),
              sapply(1:length(mSpec), function(var) unlist(sapply(1:(N+nVarGroups), function(x) {
                                                          if(IC[length(Msnp) + length(Mhostphoto) + length(Menvphoto) + var, x]) return(mmSpec[allsamples[[x]], (sum(mSpec[0:(var-1)])+1):sum(mSpec[0:var])])
                                                   }))),
              sapply(1:length(MSVD), function(var) unlist(sapply(1:(N+nVarGroups), function(x) {
                                                          if(IC[length(Msnp) + length(Mhostphoto) + length(Menvphoto) + length(mSpec) + var, x]) return(snpSVDmat[allsamples[[x]], (sum(MSVD[0:(var-1)])+1):sum(MSVD[0:var])])
                                                   }))),
              sapply(1:length(nSites), function(var) unlist(sapply(1:(N+nVarGroups), function(x) {
                                                          if(IC[length(Msnp) + length(Mhostphoto) + length(Menvphoto) + length(mSpec) + length(MSVD) + var, x]) return(siteMat[allsamples[[x]], (sum(nSites[0:(var-1)])+1):sum(nSites[0:var])])
                                                   })))))

indMP <- which(is.infinite(P))
P_max <- sapply(names(indMP), function(x) {
    temp <- biomarkersLog[,x]
    min(temp[!is.infinite(temp)], na.rm=TRUE)
})
nMP <- length(P_max)

P[is.infinite(P)] <- 0


F <- length(bacteriaFilt) + length(euksFilt) + length(transcr) + length(itsFilt)
H <- sum(!is.na(biomarkersLog)) + sum(!is.na(t2log)) + sum(!is.na(t3log)) + sum(!is.na(fabT1agg)) + sum(!is.na(fabT2agg))
G <- length(Y)


cophenbact <- cophenetic(bacttreeY.root)[colnames(bacteriaFilt),colnames(bacteriaFilt)]
bacteriaFilt_inits <- log(bacteriaFilt)
bacteriaFilt_inits[bacteriaFilt == 0] <- NA
bacteriaFilt_nuisance_inits <- numeric()
for(x in 1:nrow(bacteriaFilt)) {
    nuisance <- mean(bacteriaFilt_inits[x,],na.rm=TRUE)
    bacteriaFilt_nuisance_inits <- c(bacteriaFilt_nuisance_inits,nuisance)
    bacteriaFilt_inits[x,] <- bacteriaFilt_inits[x,] - nuisance
}
for(x in 1:nrow(bacteriaFilt)) {
    wObs <- which(!is.na(bacteriaFilt_inits[x,]))
    for(y in which(is.na(bacteriaFilt_inits[x,]))) {
        closestRel <- wObs[which.min(cophenbact[y,wObs])]
        ypres <- bacteriaFilt[,y] > 0
        closepres <- bacteriaFilt[,closestRel] > 0
        bothpres <- ypres & closepres
        if(any(bothpres))
            bacteriaFilt_inits[x,y] <- bacteriaFilt_inits[x,closestRel] + mean(bacteriaFilt_inits[bothpres,y] - bacteriaFilt_inits[bothpres,closestRel], na.rm=TRUE)
        else
            bacteriaFilt_inits[x,y] <- bacteriaFilt_inits[x,closestRel] + mean(bacteriaFilt_inits[ypres,y],na.rm=TRUE) - mean(bacteriaFilt_inits[closepres,closestRel],na.rm=TRUE)
    }
}
bacteria_intercept_inits <- apply(bacteriaFilt_inits,2,mean)
bacteriaFilt_inits <- bacteriaFilt_inits - mean(bacteria_intercept_inits)
bacteriaFilt_nuisance_inits <- bacteriaFilt_nuisance_inits + mean(bacteria_intercept_inits)
bacteria_intercept_inits <- bacteria_intercept_inits - mean(bacteria_intercept_inits)
##might want to modify the replacement by doing something like the ITS version, where the mean of multiple close relatives is used (or better some kind of ancestral state estimate)

copheneuk <- cophenetic(euktreeY.root)[colnames(euksFilt),colnames(euksFilt)]
euksFilt_inits <- log(euksFilt)
euksFilt_inits[euksFilt == 0] <- NA
euksFilt_nuisance_inits <- numeric()
for(x in 1:nrow(euksFilt)) {
    nuisance <- mean(euksFilt_inits[x,],na.rm=TRUE)
    euksFilt_nuisance_inits <- c(euksFilt_nuisance_inits,nuisance)
    euksFilt_inits[x,] <- euksFilt_inits[x,] - nuisance
}
for(x in 1:nrow(euksFilt)) {
    wObs <- which(!is.na(euksFilt_inits[x,]))
    for(y in which(is.na(euksFilt_inits[x,]))) {
        closestRel <- wObs[which.min(copheneuk[y,wObs])]
        ypres <- euksFilt[,y] > 0
        closepres <- euksFilt[,closestRel] > 0
        bothpres <- ypres & closepres
        if(any(bothpres))
            euksFilt_inits[x,y] <- euksFilt_inits[x,closestRel] + mean(euksFilt_inits[bothpres,y] - euksFilt_inits[bothpres,closestRel], na.rm=TRUE)
        else
            euksFilt_inits[x,y] <- euksFilt_inits[x,closestRel] + mean(euksFilt_inits[ypres,y],na.rm=TRUE) - mean(euksFilt_inits[closepres,closestRel],na.rm=TRUE)
    }
}
euks_intercept_inits <- apply(euksFilt_inits,2,mean)
euksFilt_inits <- euksFilt_inits - mean(euks_intercept_inits)
euksFilt_nuisance_inits <- euksFilt_nuisance_inits + mean(euks_intercept_inits)
euks_intercept_inits <- euks_intercept_inits - mean(euks_intercept_inits)


transcr_inits <- log(transcr)
transcr_inits[transcr == 0] <- NA
transcr_nuisance_inits <- numeric()
for(x in 1:nrow(transcr)) {
    nuisance <- mean(transcr_inits[x,],na.rm=TRUE)
    transcr_nuisance_inits <- c(transcr_nuisance_inits,nuisance)
    transcr_inits[x,] <- transcr_inits[x,] - nuisance
}
for(x in 1:nrow(transcr)) {
    for(y in 1:ncol(transcr)) {
        if(is.na(transcr_inits[x,y])) {
            transcr_inits[x,y] <- mean(transcr_inits[transcr[,y] > 0, y],na.rm=TRUE)
        }
    }
}
transcr_intercept_inits <- apply(transcr_inits,2,mean)
transcr_inits <- transcr_inits - mean(transcr_intercept_inits)
transcr_nuisance_inits <- transcr_nuisance_inits + mean(transcr_intercept_inits)
transcr_intercept_inits <- transcr_intercept_inits - mean(transcr_intercept_inits)


allclades <- as.character(unique(itsMeta[colnames(itsFilt), 'Clade']))
allmaj <- as.character(unique(itsMeta[colnames(itsFilt), 'Majority.ITS2.sequence']))
itsMat <- cbind(sapply(unique(itsMeta[colnames(itsFilt), 'Clade']), function(x) x == itsMeta[colnames(itsFilt), 'Clade']),
                sapply(allmaj, function(x) x == itsMeta[colnames(itsFilt), 'Majority.ITS2.sequence']))
colnames(itsMat) <- paste0('ITS2_', c(allclades,allmaj))
itsMat <- itsMat[,apply(itsMat,2,function(x) sum(x) > 1)]

covits <- tcrossprod(itsMat)
colnames(covits) <- rownames(covits) <- colnames(itsFilt)
itsFilt_inits <- log(itsFilt)
itsFilt_inits[itsFilt == 0] <- NA
itsFilt_nuisance_inits <- numeric()
for(x in 1:nrow(itsFilt)) {
    nuisance <- mean(itsFilt_inits[x,],na.rm=TRUE)
    itsFilt_nuisance_inits <- c(itsFilt_nuisance_inits,nuisance)
    itsFilt_inits[x,] <- itsFilt_inits[x,] - nuisance
}
for(x in 1:nrow(itsFilt)) {
    wObs <- which(!is.na(itsFilt_inits[x,]))
    for(y in which(is.na(itsFilt_inits[x,]))) {
        closestRel <- wObs[covits[y,wObs] == max(covits[y,wObs])]
        ypres <- itsFilt[,y] > 0
        if(length(closestRel) == 1) {
            closepres <- itsFilt[,closestRel] > 0
            bothpres <- ypres & closepres
            if(any(bothpres))
                itsFilt_inits[x,y] <- itsFilt_inits[x,closestRel] + mean(itsFilt_inits[bothpres,y] - itsFilt_inits[bothpres,closestRel], na.rm=TRUE)
            else
                itsFilt_inits[x,y] <- itsFilt_inits[x,closestRel] + mean(itsFilt_inits[ypres,y],na.rm=TRUE) - mean(itsFilt_inits[closepres,closestRel],na.rm=TRUE)
        } else {
            itsFilt_inits[x,y] <- mean(itsFilt_inits[x,closestRel]) + mean(itsFilt_inits[ypres,y],na.rm=TRUE) - mean(sapply(closestRel, function(z) mean(itsFilt_inits[itsFilt[,z] > 0, z])), na.rm=TRUE)
        }
    }
}
its_intercept_inits <- apply(itsFilt_inits,2,mean)
itsFilt_inits <- itsFilt_inits - mean(its_intercept_inits)
itsFilt_nuisance_inits <- itsFilt_nuisance_inits + mean(its_intercept_inits)
its_intercept_inits <- its_intercept_inits - mean(its_intercept_inits)


bactPhyMat <- matrix(0, NNodes + 1, NNodes + 1)
for(node in 1:(NNodes + 1)) {
    bactPhyMat[node, ] <- as.numeric(1:(NNodes + 1) %in% c(Ancestors(bacttreeY.root, node), node))
}
colnames(bactPhyMat) <- rownames(bactPhyMat) <- paste0('16S_i', 1:(NNodes + 1))
colnames(bactPhyMat)[1:NTips] <- rownames(bactPhyMat)[1:NTips] <- bacttreeY.root$tip.label
bactPhyMat <- bactPhyMat[colnames(bacteriaFilt), (NTips+2):ncol(bactPhyMat)]

eukPhyMat <- matrix(0, NNodesEuks + 1, NNodesEuks + 1)
for(node in 1:(NNodesEuks + 1)) {
    eukPhyMat[node, ] <- as.numeric(1:(NNodesEuks + 1) %in% c(Ancestors(euktreeY.root, node), node))
}
colnames(eukPhyMat) <- rownames(eukPhyMat) <- paste0('euk_i', 1:(NNodesEuks + 1))
colnames(eukPhyMat)[1:NTipsEuks] <- rownames(eukPhyMat)[1:NTipsEuks] <- euktreeY.root$tip.label
eukPhyMat <- eukPhyMat[colnames(euksFilt), (NTipsEuks+2):ncol(eukPhyMat)]

biomarkermatNames <- c('symbiont_biomass','ubiquitin')
biomarkermat <- sapply(biomarkermatNames, function(x) grepl(x,colnames(biomarkersLog),ignore.case = TRUE))

t2Mat <- cbind(1, as.numeric(grepl('raw', colnames(t2log))), as.numeric(grepl('norm', colnames(t2log))), as.numeric(grepl('Q1', colnames(t2log))), as.numeric(grepl('Q3', colnames(t2log))))
colnames(t2Mat) <- c('T2_all', 'T2_raw_all', 'T2_norm_all', 'T2_Q1', 'T2_Q3')
t3Mat <- cbind(1, as.numeric(grepl('raw', colnames(t3log))), as.numeric(grepl('norm', colnames(t3log))), as.numeric(grepl('Q1', colnames(t3log))), as.numeric(grepl('Q3', colnames(t3log))))
colnames(t3Mat) <- c('T3_all', 'T3_raw_all', 'T3_norm_all', 'T3_Q1', 'T3_Q3')

fabT1aggMatNames <- list(windspeed_group = 'speed',
                         windspeed05_group = 'Q05.*speed',
                         windspeed95_group = 'Q95.*speed',
                         true_windspeed_group = 'true_wind_speed',
                         windspeed_knots_group = 'wind_speed_knots',
                         wind_direction_group = 'true_wind_dir|eastward_wind|orthward_wind',
                         true_wind_dir_group = 'true_wind_dir',
                         barometer_hp_group = 'barometer_hp',
                         temperature_group = 'temp|sst|T_acs|T_tsg',
                         air_temperature_group = 'airTemp',
                         water_temperature_group = 'waterTemp|T_acs|sst_batos|T_tsg',
                         waterTemp_group = 'waterTemp',
                         T_acs_group = 'T_acs',
                         sst_batos_group = 'sst_batos',
                         T_tsg_group = 'T_tsg',
                         temperature05_group = 'Q05.*(temp|sst|T_acs|T_tsg)',
                         temperature95_group = 'Q95.*(temp|sst|T_acs|T_tsg)',
                         water_temperature05_group = 'Q05.*(waterTemp|T_acs|sst_batos|T_tsg)',
                         water_temperature95_group = 'Q95.*(waterTemp|T_acs|sst_batos|T_tsg)',
                         dewPoint_group = 'dewPoint',
                         salinity_group = 'S_tsg|S_acs|sal',
                         salinity05_group = 'Q05.*(S_tsg|S_acs|sal)',
                         salinity95_group = 'Q95.*(S_tsg|S_acs|sal)',
                         S_tsg_group = 'S_tsg',
                         S_acs_group = 'S_acs',
                         S_tsg_uncorr_group = 'S_tsg_uncorr',
                         S_tsg_corr_group = 'S_tsg_corr',
                         chlorophyll_group = 'Chl',
                         chlorophyll05_group = 'Q05.*Chl',
                         chlorophyll95_group = 'Q95.*Chl',
                         Chl_acs_group = 'Chl_acs$',
                         Chl_acs_final_group = 'Chl_acs_final',
                         POC_group = 'POC',
                         POC05_group = 'Q05.*POC',
                         POC95_group = 'Q95.*POC',
                         POC_acs_group = 'POC_acs$',
                         POC_acs_final_group = 'POC_acs_final',
                         gamma_group = 'gamma',
                         gamma05_group = 'Q05.*gamma',
                         gamma95_group = 'Q95.*gamma',
                         gamma_acs_group = 'gamma_acs$',
                         gamma_acs_final_group = 'gamma_acs_final',
                         fCDOM_group = 'fCDOM',
                         PAR_group = 'PAR',
                         PAR05_group = 'Q05.*PAR',
                         PAR95_group = 'Q95.*PAR',
                         PAR_bincount_group = 'PAR_bincount',
                         PAR_sd_group = 'PAR_sd',
                         PAR_higher_group = 'PAR$',
                         Twilight_group = 'Twilight',
                         AstronauticalTwilight_group = 'AstronauticalTwilight',
                         CivilTwilight_group = 'CivilTwilight',
                         NauticalTwilight_group = 'NauticalTwilight',
                         astronomical_group = 'rise|set',
                         rise_group = 'rise',
                         set_group = 'set',
                         zenith_group = 'zenith',
                         azimuth_group = 'azimuth',
                         moon_group = 'moon',
                         moon_rise_set_group = 'moon_(rise|set)')
fabT1aggMatNames <- fabT1aggMatNames[sapply(fabT1aggMatNames, function(x) sum(grepl(x,colnames(fabT1agg),ignore.case = TRUE))) > 1]
fabT1aggMat <- sapply(fabT1aggMatNames, function(x) grepl(x,colnames(fabT1agg),ignore.case = TRUE))


hostphotoHigherMatNames <- c('morphotype_Millepora','morphotype_Pocillopora','morphotype_Porites')
hostphotoHigherMat <- sapply(hostphotoHigherMatNames, function(x) grepl(x,colnames(hostphotomat2),ignore.case = TRUE))


envphotoHigherMatNames <- list(bleached_group='bleached_light|bleached_bleached', borers_group='bivalve|Spirobranchus|Tridacna|boringurchin|other_polychaete|sponge|gallcrabs|ascidians', polychaetes_group='Spirobranchus|other_polychaete', bivalve_group='bivalve|Tridacna', algae_contact_group='Halimeda|Turbinaria|Dictyota|Lobophora|Galaxaura|Sargassum', sargassaceae_group='Turbinaria|Sargassum')
envphotoHigherMatNames <- envphotoHigherMatNames[sapply(envphotoHigherMatNames, function(x) sum(grepl(x,colnames(envphotomat2),ignore.case = TRUE))) > 1]
envphotoHigherMat <- sapply(envphotoHigherMatNames, function(x) grepl(x,colnames(envphotomat2),ignore.case = TRUE))

snpSVDHigherMatNames <- c('cladePOC','cladePOR')
snpSVDHigherMat <- sapply(snpSVDHigherMatNames, function(x) grepl(x,colnames(snpSVDmat),ignore.case = TRUE))

siteHigherMatNames <- levels(filtData$island)
siteHigherMat <- sapply(siteHigherMatNames, function(x) grepl(x,colnames(siteMat),ignore.case = TRUE))

intercepts_inits <- c(bacteria_intercept_inits,
                      euks_intercept_inits,
                      transcr_intercept_inits,
                      its_intercept_inits,
                      apply(biomarkersLog_inits,2,mean,na.rm=TRUE),
                      apply(t2log,2,mean,na.rm=TRUE),
                      apply(t3log,2,mean,na.rm=TRUE),
                      apply(fabT1agg,2,mean,na.rm=TRUE),
                      apply(fabT2agg,2,mean,na.rm=TRUE),
                      unlist(sapply(1:length(Msnp), function(x) if(Msnp[[x]] == 1) {logit(mean(snpmat[,sum(Msnp[0:x])],na.rm=TRUE))} else {temp <- log(apply(snpmat[,(sum(Msnp[0:(x-1)])+1):sum(Msnp[0:x])],2,mean,na.rm=TRUE));temp - mean(temp)})),
                      unlist(sapply(1:length(Mhostphoto), function(x) if(Mhostphoto[[x]] == 1) {logit(mean(hostphotomat2[,sum(Mhostphoto[0:x])],na.rm=TRUE))} else {temp <- log(apply(hostphotomat2[,(sum(Mhostphoto[0:(x-1)])+1):sum(Mhostphoto[0:x])],2,mean,na.rm=TRUE));temp - mean(temp)})),
                      unlist(sapply(1:length(Menvphoto), function(x) if(Menvphoto[[x]] == 1) {logit(mean(envphotomat2[,sum(Menvphoto[0:x])],na.rm=TRUE))} else {temp <- log(apply(envphotomat2[,(sum(Menvphoto[0:(x-1)])+1):sum(Menvphoto[0:x])],2,mean,na.rm=TRUE));temp - mean(temp)})),
                      unlist(sapply(1:length(mSpec), function(x) if(mSpec[[x]] == 1) {logit(mean(mmSpec[,sum(mSpec[0:x])],na.rm=TRUE))} else {temp <- log(apply(mmSpec[,(sum(mSpec[0:(x-1)])+1):sum(mSpec[0:x])],2,mean,na.rm=TRUE));temp - mean(temp)})),
                      unlist(sapply(1:length(MSVD), function(x) if(MSVD[[x]] == 1) {logit(mean(snpSVDmat[,sum(MSVD[0:x])],na.rm=TRUE))} else {temp <- log(apply(snpSVDmat[,(sum(MSVD[0:(x-1)])+1):sum(MSVD[0:x])],2,mean,na.rm=TRUE));temp - mean(temp)})),
                      {temp <- log(apply(siteMat,2,mean,na.rm=TRUE));temp - mean(temp)})

intercept_fun <- function(x) {
    if(sum(x,na.rm=TRUE) != sum(!is.na(x)))
        mean(x,na.rm=TRUE)
    else
        sum(x,na.rm=TRUE)/(sum(!is.na(x))+1)
}
binary_count_intercepts_inits <- logit(c(apply(bacteriaFilt > 0,2,intercept_fun),
                                         apply(euksFilt > 0,2,intercept_fun),
                                         apply(transcr > 0,2,intercept_fun),
                                         apply(itsFilt > 0,2,intercept_fun)))

multinomial_nuisance_inits <- c(bacteriaFilt_nuisance_inits,
                                euksFilt_nuisance_inits,
                                transcr_nuisance_inits,
                                itsFilt_nuisance_inits)


biomarkersInv <- t(ginv(cbind(diag(ncol(biomarkersLog_inits)),biomarkermat)))
fabT1aggMatInv <- t(ginv(cbind(diag(ncol(fabT1agg)),fabT1aggMat)))

prior_scales <- c(apply(bacteriaFilt_inits %*% t(ginv(cbind(1,diag(ncol(bacteriaFilt_inits)),bactPhyMat))[-1,]),2,sd),
                  apply(euksFilt_inits %*% t(ginv(cbind(1,diag(ncol(euksFilt_inits)),eukPhyMat))[-1,]),2,sd),
                  apply(transcr_inits, 2, sd),
                  apply(itsFilt_inits %*% t(ginv(cbind(1,diag(ncol(itsFilt_inits)),itsMat))[-1,]),2,sd),
                  apply(sapply(1:ncol(biomarkersInv), function(x) matrix(biomarkersLog_inits[,biomarkersInv[,x]!=0],nrow=nrow(biomarkersLog_inits)) %*% biomarkersInv[biomarkersInv[,x]!=0,x]),2,sd,na.rm=TRUE),
                  apply(t2log %*% t(ginv(cbind(diag(ncol(t2log)),t2Mat))),2,sd),
                  apply(t3log %*% t(ginv(cbind(diag(ncol(t3log)),t3Mat))),2,sd),
                  apply(sapply(1:ncol(fabT1aggMatInv), function(x) matrix(fabT1agg[,fabT1aggMatInv[,x]!=0],nrow=nrow(fabT1agg)) %*% fabT1aggMatInv[fabT1aggMatInv[,x]!=0,x]),2,sd,na.rm=TRUE),
                  apply(fabT2agg,2,sd,na.rm=TRUE),
                  rep(1,ncol(snpmat)),
                  rep(1,ncol(hostphotomat2)),
                  rep(1,ncol(hostphotoHigherMat)),
                  rep(1,ncol(envphotomat2)),
                  rep(1,ncol(envphotoHigherMat)),
                  rep(1,ncol(mmSpec)),
                  rep(1,ncol(snpSVDmat)),
                  rep(1,ncol(snpSVDHigherMat)),
                  rep(1,ncol(siteMat)),
                  rep(1,ncol(siteHigherMat)))

prior_intercept_scales <- c(apply(bacteriaFilt_inits,2,sd),
                            apply(euksFilt_inits,2,sd),
                            apply(transcr_inits,2, sd),
                            apply(itsFilt_inits,2,sd),
                            apply(biomarkersLog_inits,2,sd,na.rm=TRUE),
                            apply(t2log,2,sd),
                            apply(t3log,2,sd),
                            apply(fabT1agg,2,sd,na.rm=TRUE),
                            apply(fabT2agg,2,sd,na.rm=TRUE),
                            rep(1,ncol(snpmat)),
                            rep(1,ncol(hostphotomat2)),
                            rep(1,ncol(envphotomat2)),
                            rep(1,ncol(mmSpec)),
                            rep(1,ncol(snpSVDmat)),
                            rep(1,ncol(siteMat)))

prior_intercept_centers <- intercepts_inits
binary_count_intercept_centers <- binary_count_intercepts_inits

Mplus <- M + c(ncol(bactPhyMat), ncol(eukPhyMat), 0, ncol(itsMat), ncol(biomarkermat), ncol(t2Mat), ncol(t3Mat), ncol(fabT1aggMat), 0, 0, ncol(hostphotoHigherMat), ncol(envphotoHigherMat), 0, ncol(snpSVDHigherMat), ncol(siteHigherMat))
VOBplus <- sum(Mplus)
mm <- c(bactPhyMat, eukPhyMat, itsMat, biomarkermat, t2Mat, t3Mat, fabT1aggMat, hostphotoHigherMat, envphotoHigherMat, snpSVDHigherMat, siteHigherMat)
sizeMM <- length(mm)

varlabs <- c(colnames(bacteriaFilt), colnames(bactPhyMat), colnames(euksFilt), colnames(eukPhyMat), colnames(transcr), colnames(itsFilt), colnames(itsMat), colnames(biomarkersLog), biomarkermatNames, colnames(t2log), colnames(t2Mat), colnames(t3log), colnames(t3Mat), colnames(fabT1agg), names(fabT1aggMatNames), colnames(fabT2agg), colnames(snpmat), colnames(hostphotomat2), hostphotoHigherMatNames, colnames(envphotomat2), names(envphotoHigherMatNames), colnames(mmSpec), colnames(snpSVDmat), colnames(snpSVDHigherMat), colnames(siteMat), colnames(siteHigherMat))
varlabsM <- c(colnames(bacteriaFilt), colnames(euksFilt), colnames(transcr), colnames(itsFilt), colnames(biomarkersLog), colnames(t2log), colnames(t3log), colnames(fabT1agg), colnames(fabT2agg), colnames(snpmat), colnames(hostphotomat2), colnames(envphotomat2), colnames(mmSpec), colnames(snpSVDmat), colnames(siteMat))

varsWGroupsInds <- c(sapply(colnames(fabT1agg), function(x) which(x==varlabsM)), sapply(colnames(fabT2agg), function(x) which(x==varlabsM)))

inv_log_max_contam <- 1 / (log(2) + log(c(max(apply(diag(1/rowSums(bacteriaFilt)) %*% bacteriaFilt, 2, function(x) max(x[x>0]) / min(x[x>0]))),
                                          max(apply(diag(1/rowSums(euksFilt)) %*% euksFilt, 2, function(x) max(x[x>0]) / min(x[x>0]))),
                                          max(apply(diag(1/rowSums(transcr)) %*% transcr, 2, function(x) max(x[x>0]) / min(x[x>0]))),
                                          max(apply(diag(1/rowSums(itsFilt)) %*% itsFilt, 2, function(x) max(x[x>0]) / min(x[x>0]))))))

data <- list(N = N,
             D = D,
             R = R,
             C = C,
             M = M,
             Mplus = Mplus,
             I = I,
             O = O,
             IR = IR,
             IC = IC,
             F = F,
             X = X,
             G = G,
             Y = Y,
             H = H,
             P = P,
             C_vars = C_vars,
             Mc = Mc,
             prior_scales                   = prior_scales,
             prior_intercept_scales         = prior_intercept_scales,
             prior_intercept_centers        = prior_intercept_centers,
             binary_count_intercept_centers = binary_count_intercept_centers,
             sizeMM = sizeMM,
             mm     = mm,
             global_scale_prior = global_scale_prior,
             K_linear = K_linear,
             K_gp     = K_gp,
             KG       = KG,
             nMP      = nMP,
             indMP    = indMP,
             P_max    = P_max,
             shape_gamma_fact    = shape_gamma_fact,
             rate_gamma_fact     = rate_gamma_fact,
             dist_sites          = dist_sites[lower.tri(dist_sites)],
             rho_sites_prior     = mean(dist_sites[lower.tri(dist_sites)]),
             nVarGroups          = nVarGroups,
             nVarsWGroups        = nVarsWGroups,
             samp2group          = samp2group,
             varsWGroupsInds     = varsWGroupsInds,
             ms                  = site_smoothness,
             nu_residuals        = nu_residuals,
             inv_log_max_contam  = inv_log_max_contam)

#### create initiliazations
abundance_true_vector_inits <- unlist(c(sapply(1:N, function(x) if(I[1,x]) bacteriaFilt_inits[I_cs[1,x],]),
                                        sapply(1:N, function(x) if(I[2,x]) euksFilt_inits[I_cs[2,x],]),
                                        sapply(1:N, function(x) if(I[3,x]) transcr_inits[I_cs[3,x],]),
                                        sapply(1:N, function(x) if(I[4,x]) itsFilt_inits[I_cs[4,x],])))

Z_raw <- matrix(rnorm((K_linear+KG*K_gp)*N) * 0.001, nrow=K_linear+KG*K_gp)
Z_raw <- diag(sqrt(colSums(t(Z_raw)^2))) %*% svd(t(Z_raw))$v %*% t(svd(t(Z_raw))$u)

W_raw <- matrix(rnorm((VOBplus+sum(Mplus[1:D])+D) * K) * 0.001, ncol=K)
W_raw <- svd(W_raw)$u %*% t(svd(W_raw)$v) %*% diag(sqrt(colSums(W_raw^2)))

init <- list(abundance_true_vector = abundance_true_vector_inits,
             intercepts = intercepts_inits,
             binary_count_intercepts = binary_count_intercepts_inits,
             binary_count_dataset_intercepts = rep(0,D),
             multinomial_nuisance = multinomial_nuisance_inits,
             global_effect_scale = global_scale_prior,
             latent_scales = rep(global_scale_prior * 0.05,K),
             sds = rep(1, VOBplus+sum(Mplus[1:D])+D),
             dataset_scales = rep(1, 2*D+R+C),
             nu_factors_raw = matrix(10, nrow=2*D+R+C, ncol=K),
             weight_scales = matrix(global_scale_prior * 0.05, nrow=2*D+R+C, ncol=K),
             rho_sites = as.array(rep(mean(dist_sites[lower.tri(dist_sites)]) * 0.001, K)),
             site_prop = as.array(rep(0.001, K)),
             Z_raw = Z_raw,
             W_raw = W_raw,
             P_missing = rep(-1,nMP),
             rho_Z = matrix(0.0001, nrow = K_linear, ncol = KG),
             inv_log_less_contamination = -inv_log_max_contam,
             contaminant_overdisp = rep(1,D))

save.image(file.path(output_prefix, 'setup.RData'))

write_stan_json(init, file.path(output_prefix, 'inits.json'))
write_stan_json(data, file.path(output_prefix, 'data.json'))

setwd(cmdstan_path())
system(paste0(c('make ', 'make STAN_OPENCL=true ')[opencl+1], file.path(model_dir, model_name)))

setwd(model_dir)
print(sampling_commands[[engine]])
print(date())
system(sampling_commands[[engine]])

importparams <- c('W_norm','Z','sds','latent_scales','global_effect_scale','log_post','dataset_scales','var_scales','weight_scales','nu_factors','rho_sites', 'site_prop', 'cov_sites', 'binary_count_dataset_intercepts','log_less_contamination','contaminant_overdisp','rho_Z')

stan.fit.vb <- read_stan_csv_subset(file.path(output_prefix, paste0('samples_',engine,'.txt')),
                                    params = importparams)

save.image(file.path(output_prefix, paste0('res_',engine,'.RData')))

source(file.path(model_dir, 'GPBFA_extract.r'))

save.image(file.path(output_prefix, paste0('res_',engine,'_processed.RData')))
