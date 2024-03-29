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
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
logit <- function(p) log(p/(1-p))
inv_logit <- function(x) { 1 / (1 + exp(-x)) }

nMicrobeKeep <- 100
K_linear <- 35#10
K_gp <- 15
KG <- 0#3
K <- K_linear + KG * K_gp
global_scale_prior = 2.5
rate_gamma_fact = 10
shape_gamma_fact = 2
site_smoothness <- 2
nu_residuals <- 5
shape_gnorm <- 7

input_prefix <- file.path(Sys.getenv('HOME'), 'data/tara_unsupervised_analyses')
if(exists('myargs')) {if(length(myargs)==1) {input_prefix <- myargs[[1]]}} else if(length(commandArgs(trailingOnly=TRUE)) > 0) {input_prefix <- commandArgs(trailingOnly=TRUE)[[1]]}
preprocess_prefix <- paste0(Sys.getenv('HOME'), '/outputs/tara/intermediate/')
include_path <- file.path(Sys.getenv('HOME'), 'scripts/stan_models/utility/')
model_dir <- file.path(Sys.getenv('HOME'), 'scripts/stan_models/BGPFA/')
source(file.path(model_dir, 'functions.r'))
model_name <- 'BGPFA'
engine <- c('sampling','advi')[1]
mass_slow <- 1 / 100 # multiplied by select raw parameters so that they initially are explored slower during sampling
opencl <- FALSE
output_prefix <- paste0(Sys.getenv('HOME'), '/outputs/tara/BGPFA_', nMicrobeKeep)

dir.create(output_prefix, recursive = TRUE)

sampling_commands <- list(sampling = paste(paste0('./', model_name),
                                          paste0('data file=', file.path(output_prefix, 'data.json')),
                                          paste0('init=', file.path(output_prefix, 'inits.json')),
                                          'output',
                                          paste0('file=', file.path(output_prefix, 'samples_sampling.csv')),
                                          paste0('refresh=', 1),
                                          'method=sample algorithm=hmc',
                                          'stepsize=1',
                                          'engine=nuts',
                                          'max_depth=10',
                                          'adapt t0=10',
                                          'delta=0.8',
                                          'kappa=0.75',
                                          #'init_buffer=10',
                                          #'window=2',
                                          'term_buffer=50',
                                          'num_warmup=1000',
                                          'num_samples=1000',
                                          ('opencl platform=0 device=0')[opencl],
                                          sep=' '),
                          advi = paste(paste0('./', model_name),
                                       paste0('data file=', file.path(output_prefix, 'data.json')),
                                       paste0('init=', file.path(output_prefix, 'inits.json')),
                                       'output',
                                       paste0('file=', file.path(output_prefix, 'samples_advi.csv')),
                                       paste0('refresh=', 100),
                                       'method=variational algorithm=meanfield',
                                       'grad_samples=1',
                                       'elbo_samples=1',
                                       'iter=20000',
                                       'eta=1',
                                       'adapt engaged=0',
                                       'tol_rel_obj=0.0001',
                                       'eval_elbo=1',
                                       'output_samples=200',
                                       ('opencl platform=0 device=0')[opencl],
                                       sep=' '))

dataset_names <- c('mb16S','mb18S','rna','its2','biomarkers','t2','t3','fab_T1_agg','fab_T2_agg','snps','hphoto','ephoto','species','svd','sites')
in_data <- list()
inits_abundance <- list()
inits_nuisance <- list()
inits_intercepts <- list()
mm_h <- list()
names_mm_h <- list()

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

in_data$mb16S <- t(bacteria[keepers,])
rownames(in_data$mb16S) <- sample_data[match(rownames(in_data$mb16S), sample_data$sample.storage_container.label), 'myname']
in_data$mb16S <-in_data$mb16S[order(rownames(in_data$mb16S)),]

bacttree <- read.tree(file.path(preprocess_prefix, '20201231/TARA_PACIFIC_16S_asv_CO-0-filtered_combined_filtered_aligned.tree'))
bacttree$tip.label <- sub('^_R_','',bacttree$tip.label)
bacttreeY <- drop.tip(bacttree, bacttree$tip.label[!bacttree$tip.label %in% colnames(in_data$mb16S)])

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
bactsnotintree <- colnames(in_data$mb16S)[!colnames(in_data$mb16S) %in% bacttreeY.root$tip.label]
bacttreeY.root <- add.tips(bacttreeY.root, bactsnotintree, rep(length(bacttreeY.root$tip.label)+1,length(bactsnotintree)), rep(1,length(bactsnotintree)))

in_data$mb16S <- in_data$mb16S[,bacttreeY.root$tip.label]

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
} ## filter based on shapiro-wilk divergence from normality rather than on variance? non-normal things likely to have small number of explanatory factors? need deseq2 shrunk data?
keepers_euk <- sort(keepers_euk)

in_data$mb18S <- t(euks[keepers_euk,])
rownames(in_data$mb18S) <- sample_data[match(rownames(in_data$mb18S), sample_data$sample.storage_container.label), 'myname']
in_data$mb18S <-in_data$mb18S[order(rownames(in_data$mb18S)),]

euktree <- read.tree(file.path(preprocess_prefix, '20210102/TARA_PACIFIC_18SV9_4191_samples_v202004.OTU.filtered_CO-0-filtered_combined_filtered_aligned.tree'))
euktree$tip.label <- sub('^_R_','',euktree$tip.label)
euktreeY <- drop.tip(euktree, euktree$tip.label[!euktree$tip.label %in% colnames(in_data$mb18S)])
euktreeY.root <- MidpointRooter:::midpoint.root2(euktreeY)
euksnotintree <- colnames(in_data$mb18S)[!colnames(in_data$mb18S) %in% euktree$tip.label]
euktreeY.root <- add.tips(euktreeY.root, euksnotintree, rep(length(euktreeY.root$tip.label)+1,length(euksnotintree)), rep(1,length(euksnotintree)))

in_data$mb18S <- in_data$mb18S[,euktreeY.root$tip.label]

NTipsEuks <- length(euktreeY.root$tip.label)
NNodesEuks <- euktreeY.root$Nnode + NTipsEuks - 1


###

in_data$rna <- as.matrix(read.csv(file.path(input_prefix, 'Alice/20190924/transcriptomics/tara_I10_trans.csv'), sep=';', row.names=1, stringsAsFactors=TRUE)[,1:48])
rownames(in_data$rna) <- paste(substr(rownames(in_data$rna), 1, 7), 0, substr(rownames(in_data$rna), 8, nchar(rownames(in_data$rna))), sep = "")
in_data$rna <- in_data$rna[,colSums(in_data$rna) > 0]
in_data$rna <- t(apply(in_data$rna,1,function(x) round(10 * x)))
in_data$rna <- cbind(in_data$rna, apply(in_data$rna,1,function(x) 10000000 - sum(x))) #data are normalized out of a million reads and have accuracy of 1/100th. to treat them as counts i am going to assume all samples had exactly 10 million reads originally - because i suspect the siginificant digits are optimistic and this number was not actually sequenced. i'm hoping more were sequenced, but assuming a lower count makes them less confident and this is therefore conservative
colnames(in_data$rna)[ncol(in_data$rna)] <- 'X.Other'
colnames(in_data$rna) <- sub('^X.','transcr_',colnames(in_data$rna))

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
in_data$its2 <- its
rownames(in_data$its2) <- sample_data[match(rownames(in_data$its2), sample_data$sample.storage_container.label), 'myname']

itsMeta <- read.table(file=file.path(input_prefix, 'Ben/20200212/its2_type_profiles/90_20200125_DBV_2020-01-28_05-39-25.985991.profiles.absolute.meta_only.txt'), row.names=1, header=T, sep='\t', comment.char='', stringsAsFactors=TRUE)
rownames(itsMeta) <- paste0('ITS2_', rownames(itsMeta))

## import paola's biomarker data
biomarkers <- read.table(file.path(input_prefix, 'Paola/20200218/paola_readable.txt'), sep='\t', header=T, row.names=1, stringsAsFactors=TRUE)
rownames(biomarkers) <- gsub('-','',sub('OA000-','',rownames(biomarkers)))

## import didier's snp data
snps <- read.csv(file.path(input_prefix, 'Didier/geneAlex_LD02.csv'), skip=2, row.names=2, stringsAsFactors=TRUE)[,-1]

rownames(snps) <- gsub('o',0,sapply(rownames(snps), function(x) paste0(c('I', substr(x,2,3), 'S', substr(x,5,6), 'C', substr(x,8,10)), collapse='')))

nSNPs <- ncol(snps) / 2

in_data$snps <- matrix(0, nrow = nrow(snps), ncol = 4 * nSNPs)
colnames(in_data$snps) <- 1:(4 * nSNPs)
rownames(in_data$snps) <- rownames(snps)
M_snp <- vector('numeric')
newsnpnames <- vector()
snpFilter <- vector()
indSNP <- 1
indSNPMat <- 1
for(i in 1:nSNPs) {
    in_data$snps[, indSNPMat:(indSNPMat + 3)] <- t(apply(snps[, indSNP:(indSNP+1)], 1, function(x) {
        sapply(1:4, function(y) sum(y == x))
    }))
    colnames(in_data$snps)[indSNPMat:(indSNPMat + 3)] <- paste(colnames(snps)[indSNP], 1:4, sep='.')
    newFilt <- apply(in_data$snps[, indSNPMat:(indSNPMat + 3)], 2, function(x) any(x > 0) & sd(x) > 0)
    if(sum(newFilt) == 2) {
        snpPoss <- sort(unique(in_data$snps[,(indSNPMat:(indSNPMat + 3))[which(newFilt)[[2]]]]))
        if(length(snpPoss) == 2) {

            poss1 <- in_data$snps[,(indSNPMat:(indSNPMat + 3))[which(newFilt)[[1]]]] == snpPoss[[1]]
            poss2 <- in_data$snps[,(indSNPMat:(indSNPMat + 3))[which(newFilt)[[1]]]] == snpPoss[[2]]

            if(sum(poss1) > 1 & sum(poss2) > 1) {
                newsnpnames <- c(newsnpnames, paste(names(newFilt)[which(newFilt)[[2]]], names(newFilt)[which(newFilt)[[1]]], sep='_vs_'))
                newFilt[which(newFilt)[[1]]] <- FALSE

                in_data$snps[poss1,(indSNPMat:(indSNPMat + 3))[newFilt]] <- 0
                in_data$snps[poss2,(indSNPMat:(indSNPMat + 3))[newFilt]] <- 1

                M_snp <- c(M_snp, 1)
            } else {
                newFilt[newFilt] <- FALSE
            }
        }
        else if(length(snpPoss) > 2) {

            het <- in_data$snps[,(indSNPMat:(indSNPMat + 3))[which(newFilt)[[1]]]] == 1
            hom1 <- in_data$snps[,(indSNPMat:(indSNPMat + 3))[which(newFilt)[[1]]]] == 2
            hom2 <- in_data$snps[,(indSNPMat:(indSNPMat + 3))[which(newFilt)[[1]]]] == 0

            homname <- paste(names(newFilt)[which(newFilt)[[2]]], names(newFilt)[which(newFilt)[[1]]], sep='_vs_')
            in_data$snps[het, (indSNPMat:(indSNPMat + 3))[which(newFilt)[[1]]]] <- NA
            in_data$snps[hom1, (indSNPMat:(indSNPMat + 3))[which(newFilt)[[1]]]] <- 0
            in_data$snps[hom2, (indSNPMat:(indSNPMat + 3))[which(newFilt)[[1]]]] <- 1

            hetname <- paste(names(newFilt)[which(newFilt)[[2]]], names(newFilt)[which(newFilt)[[1]]], sep='_AND_')
            in_data$snps[het, (indSNPMat:(indSNPMat + 3))[which(newFilt)[[2]]]] <- 1
            in_data$snps[hom1, (indSNPMat:(indSNPMat + 3))[which(newFilt)[[2]]]] <- 0
            in_data$snps[hom2, (indSNPMat:(indSNPMat + 3))[which(newFilt)[[2]]]] <- 0

            if(sum(hom1) > 1 & sum(hom2) > 1) {
                newsnpnames <- c(newsnpnames, homname)
                M_snp <- c(M_snp, 1)
            } else {
                newFilt[which(newFilt)[[1]]] <- FALSE
            }

            if(sum(het) > 1 & sum(!het) > 1) {
                newsnpnames <- c(newsnpnames, hetname)
                M_snp <- c(M_snp, 1)
            } else {
                newFilt[which(newFilt)[[2]]] <- FALSE
            }
        }
    } else if(sum(newFilt) > 2){
        newsnpnames <- c(newsnpnames, paste(colnames(snps)[indSNP], 1:4, sep='.'))
        M_snp <- c(M_snp, sum(newFilt))
    }
    snpFilter <- c(snpFilter, newFilt)
    indSNP <- indSNP + 2
    indSNPMat <- indSNPMat + 4
}
in_data$snps <- in_data$snps[,snpFilter]
colnames(in_data$snps) <- newsnpnames

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
M_hphoto <- nlevs

in_data$hphoto <- t(sapply(photosamples, function(x) if(sum(photonames == x) > 1) {apply(hostphotomat[photonames == x,],2,function(y) sum(y,na.rm=TRUE)/length(y))} else {hostphotomat[photonames == x,]}))

in_data$hphoto <- cbind(in_data$hphoto, photodatagg[rownames(in_data$hphoto),c('porites_lumps','porites_ridges')])
M_hphoto <- c(M_hphoto, 1, 1)

## env
bleachLevs <- unique(photodat[,'bleached'])[!is.na(unique(photodat[,'bleached'])) & unique(photodat[,'bleached']) != '']
M_ephoto <- length(bleachLevs)

envphotomat <- sapply(bleachLevs, function(y) y == photodat$bleached)
colnames(envphotomat) <- paste('bleached', bleachLevs, sep='_')
in_data$ephoto <- t(sapply(photosamples, function(x) if(sum(photonames == x) > 1) {apply(envphotomat[photonames == x,],2,function(y) sum(y,na.rm=TRUE)/length(y))} else {envphotomat[photonames == x,]}))

in_data$ephoto <- cbind(in_data$ephoto, photodatagg[,envvars[envvars != 'bleached']])
in_data$ephoto <- in_data$ephoto[,apply(in_data$ephoto,2,sd,na.rm=TRUE) > 0]
M_ephoto <- c(M_ephoto, rep(1,length(envvars[envvars %in% colnames(in_data$ephoto,2)])))

########

repeatedsamples <- c(rownames(in_data$mb16S), rownames(in_data$mb18S), rownames(in_data$rna), rownames(in_data$its2), rownames(biomarkers), rownames(t2in), rownames(t3in), rownames(in_data$snps), rownames(in_data$hphoto), rownames(in_data$ephoto))
allsamples <- unique(repeatedsamples)[sapply(unique(repeatedsamples), function(x) sum(x==repeatedsamples)>1)]

in_data$mb16S <- in_data$mb16S[allsamples[allsamples %in% rownames(in_data$mb16S)],]
in_data$mb18S <- in_data$mb18S[allsamples[allsamples %in% rownames(in_data$mb18S)],]
in_data$rna <- in_data$rna[allsamples[allsamples %in% rownames(in_data$rna)],]
in_data$its2 <- in_data$its2[allsamples[allsamples %in% rownames(in_data$its2)],]
in_data$its2 <- in_data$its2[,apply(in_data$its2,2,function(x) sum(x>0)) > 0]

in_data$t2 <- as.matrix(log(t2in[allsamples[allsamples %in% rownames(t2in)],]))
in_data$t3 <- as.matrix(log(t3in[allsamples[allsamples %in% rownames(t3in)],]))

in_data$biomarkers <- log(biomarkers[allsamples[allsamples %in% rownames(biomarkers)],])
inits_biomarkers <- in_data$biomarkers
for(i in 1:ncol(inits_biomarkers)){
    inits_biomarkers[is.infinite(inits_biomarkers[,i]),i] <- min(inits_biomarkers[!is.infinite(inits_biomarkers[,i]),i], na.rm=TRUE)
}
inits_biomarkers <- as.matrix(inits_biomarkers)

in_data$snps <- in_data$snps[allsamples[allsamples %in% rownames(in_data$snps)],]

in_data$hphoto <- in_data$hphoto[allsamples[allsamples %in% rownames(in_data$hphoto)],]
keepcol <- rep(TRUE,ncol(in_data$hphoto))
MFilt <- M_hphoto
for(u in 1:length(M_hphoto)) {
    for(p in 1:M_hphoto[u]) {
        if(sum(in_data$hphoto[,sum(M_hphoto[0:(u-1)])+p], na.rm=TRUE) == 0) {
            keepcol[sum(M_hphoto[0:(u-1)])+p] <- FALSE
            MFilt[u] <- MFilt[u] - 1
        }
    }
}
in_data$hphoto <- in_data$hphoto[,keepcol]
M_hphoto <- MFilt

in_data$ephoto <- in_data$ephoto[allsamples[allsamples %in% rownames(in_data$ephoto)],]
keepcol <- rep(TRUE,ncol(in_data$ephoto))
MFilt <- M_ephoto
for(u in 1:length(M_ephoto)) {
    for(p in 1:M_ephoto[u]) {
        if(sum(in_data$ephoto[,sum(M_ephoto[0:(u-1)])+p], na.rm=TRUE) == 0) {
            keepcol[sum(M_ephoto[0:(u-1)])+p] <- FALSE
            MFilt[u] <- MFilt[u] - 1
        }
    }
}
in_data$ephoto <- in_data$ephoto[,keepcol]
M_ephoto <- MFilt


filtData <- sample_data[match(allsamples[allsamples %in% sample_data$myname], sample_data$myname),]
filtData <- droplevels(filtData)
rownames(filtData) <- filtData$myname
filtData$site <- as.factor(substr(filtData$myname,1,6))
filtData$island <- as.factor(substr(filtData$myname,1,3))

in_data$species <- model.matrix(~0 + species, data = filtData)
rownames(in_data$species) <- filtData$myname[!is.na(filtData$species)]
in_data$species <- in_data$species[allsamples[allsamples %in% rownames(in_data$species)],]
M_spec <- ncol(in_data$species)


### didier's pop gen analysis

snpSVDs <- read.table(file.path(input_prefix,'Didier/20200419/snp_clusters.txt'), header=T, sep='\t', quote = '', row.names = 1, stringsAsFactors=FALSE)
in_data$svd <- model.matrix(~0 + clade, data = snpSVDs)
in_data$svd <- in_data$svd[allsamples[allsamples %in% rownames(in_data$svd)],]
M_svd <- ncol(in_data$svd)

###

fabT1 <- read.table(file.path(input_prefix, 'Fabien/20200127/TotalTable_meteo_nav_tsg_acs_par_cdom.txt'), sep='\t', header=T, comment.char='', stringsAsFactors=TRUE)
sites <- substr(fabT1$Var1, 21,27)
keep <- sapply(sites, function(x) strsplit(x,'')[[1]][[1]] == 'I' & (strsplit(x,'')[[1]][[6]] != '0' | strsplit(x,'')[[1]][[7]] != '0')) & (grepl('3x10', fabT1$Var1) | grepl('3X10', fabT1$Var1))
fabT1 <- fabT1[,10:ncol(fabT1)]
fabT1 <- fabT1[,!grepl('lat|lon|dire|unixtime|current_speed|Course_Over_Ground|Speed_over_Ground|apparent',colnames(fabT1))]
fabT1 <- as.data.frame(apply(fabT1,2,function(x) if(any(x<=0,na.rm=T)) x else log(x)))

fabT1sites <- unique(sites[keep])

in_data$fab_T1_agg <- array(0, dim = c(length(fabT1sites), ncol(fabT1)), dimnames = list(fabT1sites, colnames(fabT1)))
for(f in fabT1sites) {
    in_data$fab_T1_agg[f,] <- apply(fabT1[sites==f & keep,], 2, mean, na.rm = TRUE)
}
repeats <- substr(grep('Q50',colnames(in_data$fab_T1_agg),ignore.case = TRUE,value=TRUE),4,100)
in_data$fab_T1_agg <- in_data$fab_T1_agg[,!colnames(in_data$fab_T1_agg) %in% c(repeats, paste0('Q25',repeats), paste0('Q75',repeats))]
in_data$fab_T1_agg <- in_data$fab_T1_agg[,apply(in_data$fab_T1_agg, 2, function(x) {temp <- sd(x,na.rm=T); temp > 0 & !is.na(temp)})]

fabT2 <- read.table(file.path(input_prefix, 'Fabien/20190329/TotalTable_Satellites_2S.txt'), sep='\t', header=T, comment.char='', stringsAsFactors=TRUE)
islands <- substr(fabT2$st_ID, 1,3)
keep <- sapply(islands, function(x) strsplit(x,'')[[1]][[2]] != '0' | strsplit(x,'')[[1]][[2]] != '0')
fabT2 <- fabT2[,25:(ncol(fabT2)-1)]
fabT2 <- fabT2[,apply(fabT2, 2, function(x) {temp <- sd(x,na.rm=T); temp > 0 & !is.na(temp)})]
fabT2 <- log(fabT2[,!grepl('_SD',colnames(fabT2))])

fabT2isls <- unique(islands[keep])
in_data$fab_T2_agg <- array(0, dim = c(length(fabT2isls), ncol(fabT2)), dimnames = list(fabT2isls, colnames(fabT2)))
for(f in fabT2isls) {
    in_data$fab_T2_agg[f,] <- apply(fabT2[islands==f & keep,], 2, mean, na.rm = TRUE)
}


allVarGroups <- unique(c(fabT1sites,fabT2isls))
N_var_groups <- length(allVarGroups)

in_data$fab_T1_agg <- in_data$fab_T1_agg[allVarGroups[allVarGroups %in% rownames(in_data$fab_T1_agg)],]
in_data$fab_T2_agg <- in_data$fab_T2_agg[allVarGroups[allVarGroups %in% rownames(in_data$fab_T2_agg)],]


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

in_data$sites <- model.matrix(~0+site, data=filtData)
in_data$sites <- in_data$sites[allsamples[allsamples %in% rownames(in_data$sites)],]
N_sites <- ncol(in_data$sites)
##

in_data <- in_data[dataset_names]
N <- length(allsamples)
N_all <- N + N_var_groups
D <- 4
R <- 5
C <- 6
M <- sapply(in_data, ncol)
VOB <- sum(M)
samp2group <- cbind(sapply(fabT1sites, function(x) as.integer(x == filtData[allsamples,'site'])), sapply(fabT2isls, function(x) as.integer(x == filtData[allsamples,'island'])))
samp2group <- samp2group %*% diag(1.0 / colSums(samp2group, na.rm=TRUE))
samp2group[is.na(samp2group)] <- 0
Mc <- c(M_snp, M_hphoto, M_ephoto, M_spec, M_svd, N_sites)
C_vars <- length(Mc)

I <- matrix(0, nrow=D, ncol=N_all)
I[1:D,1:N] <- t(sapply(1:D, function(d) as.integer(allsamples %in% rownames(in_data[[d]]))))
dimnames(I) <- list(dataset_names[1:D], c(allsamples, paste0('varGroup', 1:N_var_groups)))

O <- sum(M[(D+1):(D+R)])
IR <- matrix(0, nrow=0, ncol=N_all) #analogous to I above but per-variable rather than per-dataset
for(y in (D+1):(D+R)) {
    IR <- rbind(IR,
                t(apply(in_data[[y]], 2, function(x) as.integer(c(allsamples,allVarGroups) %in% rownames(in_data[[y]])[!is.na(x)]))))
}
colnames(IR) <- c(allsamples, paste0('varGroup', 1:N_var_groups))

B <- sum(Mc)
IC <- matrix(0, nrow=0, ncol=N_all)
for(y in (D+R+1):(D+R+C)) {
    IC <- rbind(IC,
                t(apply(in_data[[y]], 2, function(x) as.integer(c(allsamples,allVarGroups) %in% rownames(in_data[[y]])[!is.na(x)]))))
}
colnames(IC) <- c(allsamples, paste0('varGroup', 1:N_var_groups))


ICv <- matrix(0, nrow=C_vars, ncol=N_all)
for(var in 1:length(M_snp)) {
    if(M_snp[var] > 1) {
        ICv[var,1:N] <- as.integer(allsamples %in% rownames(in_data$snps)[apply(!is.na(in_data$snps[,(sum(M_snp[0:(var-1)])+1):sum(M_snp[0:var])]), 1, all)])
    } else {
        ICv[var,1:N] <- as.integer(allsamples %in% rownames(in_data$snps)[!is.na(in_data$snps[,sum(M_snp[0:var])])])
    }
}
for(var in 1:length(M_hphoto)) {
    if(M_hphoto[var] > 1) {
        ICv[length(M_snp) + var,1:N] <- as.integer(allsamples %in% rownames(in_data$hphoto)[apply(!is.na(in_data$hphoto[,(sum(M_hphoto[0:(var-1)])+1):sum(M_hphoto[0:var])]), 1, all)])
    } else {
        ICv[length(M_snp) + var,1:N] <- as.integer(allsamples %in% rownames(in_data$hphoto)[!is.na(in_data$hphoto[,sum(M_hphoto[0:var])])])
    }
}
for(var in 1:length(M_ephoto)) {
    if(M_ephoto[var] > 1) {
        ICv[length(M_snp) + length(M_hphoto) + var,1:N] <- as.integer(allsamples %in% rownames(in_data$ephoto)[apply(!is.na(in_data$ephoto[,(sum(M_ephoto[0:(var-1)])+1):sum(M_ephoto[0:var])]), 1, all)])
    } else {
        ICv[length(M_snp) + length(M_hphoto) + var,1:N] <- as.integer(allsamples %in% rownames(in_data$ephoto)[!is.na(in_data$ephoto[,sum(M_ephoto[0:var])])])
    }
}
for(var in 1:length(M_spec)) {
    if(M_spec[var] > 1) {
        ICv[length(M_snp) + length(M_hphoto) + length(M_ephoto) + var,1:N] <- as.integer(allsamples %in% rownames(in_data$species)[apply(!is.na(in_data$species[,(sum(M_spec[0:(var-1)])+1):sum(M_spec[0:var])]), 1, all)])
    } else {
        ICv[length(M_snp) + length(M_hphoto) + length(M_ephoto) + var,1:N] <- as.integer(allsamples %in% rownames(in_data$species)[!is.na(in_data$species[,sum(M_spec[0:var])])])
    }
}
for(var in 1:length(M_svd)) {
    if(M_svd[var] > 1) {
        ICv[length(M_snp) + length(M_hphoto) + length(M_ephoto) + length(M_spec) + var,1:N] <- as.integer(allsamples %in% rownames(in_data$svd)[apply(!is.na(in_data$svd[,(sum(M_svd[0:(var-1)])+1):sum(M_svd[0:var])]), 1, all)])
    } else {
        ICv[length(M_snp) + length(M_hphoto) + length(M_ephoto) + length(M_spec) + var,1:N] <- as.integer(allsamples %in% rownames(in_data$svd)[!is.na(in_data$svd[,sum(M_svd[0:var])])])
    }
}
for(var in 1:length(N_sites)) {
    if(N_sites[var] > 1) {
        ICv[length(M_snp) + length(M_hphoto) + length(M_ephoto) + length(M_spec) + length(M_svd) + var,1:N] <- as.integer(allsamples %in% rownames(in_data$sites)[apply(!is.na(in_data$sites[,(sum(N_sites[0:(var-1)])+1):sum(N_sites[0:var])]), 1, all)])
    } else {
        ICv[length(M_snp) + length(M_hphoto) + length(M_ephoto) + length(M_spec) + length(M_svd) + var,1:N] <- as.integer(allsamples %in% rownames(in_data$sites)[!is.na(in_data$sites[,sum(N_sites[0:var])])])
    }
}
colnames(ICv) <- c(allsamples, paste0('varGroup',1:N_var_groups))

I_cs <- t(apply(I, 1, cumsum))
IR_cs <- t(apply(IR, 1, cumsum))
IC_cs <- t(apply(ICv, 1, cumsum))

X <- unlist(sapply(1:D, function(d) sapply(1:N_all, function(x) if(I[d,x]) in_data[[d]][I_cs[d,x],])))

P <- unlist(c(sapply(1:N_all, function(x) unlist(sapply(1:ncol(in_data$biomarkers), function(y) {
                                                          temp <- in_data$biomarkers[!is.na(in_data$biomarkers[,y]),y]
                                                          names(temp) <- rep(colnames(in_data$biomarkers)[y], length(temp))
                                                          if(IR[y,x]) temp[IR_cs[y,x]]
                                                   }))),
              sapply(1:N_all, function(x) unlist(sapply(1:ncol(in_data$t2), function(y) {
                                                          temp <- in_data$t2[!is.na(in_data$t2[,y]),y]
                                                          names(temp) <- rep(colnames(in_data$t2)[y], length(temp))
                                                          if(IR[ncol(in_data$biomarkers) + y,x]) temp[IR_cs[ncol(in_data$biomarkers) + y,x]]
                                                   }))),
              sapply(1:N_all, function(x) unlist(sapply(1:ncol(in_data$t3), function(y) {
                                                          temp <- in_data$t3[!is.na(in_data$t3[,y]),y]
                                                          names(temp) <- rep(colnames(in_data$t3)[y], length(temp))
                                                          if(IR[ncol(in_data$biomarkers) + ncol(in_data$t2) + y,x]) temp[IR_cs[ncol(in_data$biomarkers) + ncol(in_data$t2) + y,x]]
                                                   }))),
              sapply(1:N_all, function(x) unlist(sapply(1:ncol(in_data$fab_T1_agg), function(y) {
                                                          temp <- in_data$fab_T1_agg[!is.na(in_data$fab_T1_agg[,y]),y]
                                                          names(temp) <- rep(colnames(in_data$fab_T1_agg)[y], length(temp))
                                                          if(IR[ncol(in_data$biomarkers) + ncol(in_data$t2) + ncol(in_data$t3) + y,x]) temp[IR_cs[ncol(in_data$biomarkers) + ncol(in_data$t2) + ncol(in_data$t3) + y,x]]
                                                   }))),
              sapply(1:N_all, function(x) unlist(sapply(1:ncol(in_data$fab_T2_agg), function(y) {
                                                          temp <- in_data$fab_T2_agg[!is.na(in_data$fab_T2_agg[,y]),y]
                                                          names(temp) <- rep(colnames(in_data$fab_T2_agg)[y], length(temp))
                                                          if(IR[ncol(in_data$biomarkers) + ncol(in_data$t2) + ncol(in_data$t3) + ncol(in_data$fab_T1_agg) + y,x]) temp[IR_cs[ncol(in_data$biomarkers) + ncol(in_data$t2) + ncol(in_data$t3) + ncol(in_data$fab_T1_agg) + y,x]]
                                                   })))))

Y <- unlist(c(sapply(1:length(M_snp), function(var) unlist(sapply(1:N_all, function(x) {
                                                          if(ICv[var, x]) return(in_data$snps[allsamples[[x]], (sum(M_snp[0:(var-1)])+1):sum(M_snp[0:var])])
                                                   }))),
              sapply(1:length(M_hphoto), function(var) unlist(sapply(1:N_all, function(x) {
                                                          if(ICv[length(M_snp) + var, x]) return(in_data$hphoto[allsamples[[x]], (sum(M_hphoto[0:(var-1)])+1):sum(M_hphoto[0:var])])
                                                   }))),
              sapply(1:length(M_ephoto), function(var) unlist(sapply(1:N_all, function(x) {
                                                          if(ICv[length(M_snp) + length(M_hphoto) + var, x]) return(in_data$ephoto[allsamples[[x]], (sum(M_ephoto[0:(var-1)])+1):sum(M_ephoto[0:var])])
                                                   }))),
              sapply(1:length(M_spec), function(var) unlist(sapply(1:N_all, function(x) {
                                                          if(ICv[length(M_snp) + length(M_hphoto) + length(M_ephoto) + var, x]) return(in_data$species[allsamples[[x]], (sum(M_spec[0:(var-1)])+1):sum(M_spec[0:var])])
                                                   }))),
              sapply(1:length(M_svd), function(var) unlist(sapply(1:N_all, function(x) {
                                                          if(ICv[length(M_snp) + length(M_hphoto) + length(M_ephoto) + length(M_spec) + var, x]) return(in_data$svd[allsamples[[x]], (sum(M_svd[0:(var-1)])+1):sum(M_svd[0:var])])
                                                   }))),
              sapply(1:length(N_sites), function(var) unlist(sapply(1:N_all, function(x) {
                                                          if(ICv[length(M_snp) + length(M_hphoto) + length(M_ephoto) + length(M_spec) + length(M_svd) + var, x]) return(in_data$sites[allsamples[[x]], (sum(N_sites[0:(var-1)])+1):sum(N_sites[0:var])])
                                                   })))))

idx_Pm <- which(is.infinite(P))
idx_Pnm <- which(!is.infinite(P))
P_max <- sapply(names(idx_Pm), function(x) {
    temp <- in_data$biomarkers[,x]
    min(temp[!is.infinite(temp)], na.rm=TRUE)
})
N_Pm <- length(P_max)

P[is.infinite(P)] <- P_max

F <- sum(sapply(in_data[1:D],length))
H <- sum(sapply((D+1):(D+R), function(x) sum(!is.na(in_data[[x]]))))
G <- length(Y)


mm_h$mb16S <- matrix(0, NNodes + 1, NNodes + 1)
for(node in 1:(NNodes + 1)) {
    mm_h$mb16S[node, ] <- as.numeric(1:(NNodes + 1) %in% c(Ancestors(bacttreeY.root, node), node))
}
colnames(mm_h$mb16S) <- rownames(mm_h$mb16S) <- paste0('16S_i', 1:(NNodes + 1))
colnames(mm_h$mb16S)[1:NTips] <- rownames(mm_h$mb16S)[1:NTips] <- bacttreeY.root$tip.label
mm_h$mb16S <- mm_h$mb16S[colnames(in_data$mb16S), (NTips+2):ncol(mm_h$mb16S)]

cophenbact <- cophenetic(bacttreeY.root)[colnames(in_data$mb16S),colnames(in_data$mb16S)]
inits_abundance$mb16S <- log(in_data$mb16S)
inits_abundance$mb16S[in_data$mb16S == 0] <- NA
inits_nuisance$mb16S <- numeric()
for(x in 1:nrow(in_data$mb16S)) {
    nuisance <- mean(inits_abundance$mb16S[x,],na.rm=TRUE)
    inits_nuisance$mb16S <- c(inits_nuisance$mb16S,nuisance)
    inits_abundance$mb16S[x,] <- inits_abundance$mb16S[x,] - nuisance
}
inits_intercepts$mb16S <- apply(inits_abundance$mb16S,2,mean,na.rm=TRUE)
inits_abundance$mb16S <- inits_abundance$mb16S - mean(inits_intercepts$mb16S)
inits_nuisance$mb16S <- inits_nuisance$mb16S + mean(inits_intercepts$mb16S)
inits_intercepts$mb16S <- inits_intercepts$mb16S - mean(inits_intercepts$mb16S)
prior_intercept_scales_mb16S <- apply(inits_abundance$mb16S,2,sd,na.rm=TRUE)
prior_intercept_scales_mb16S[prior_intercept_scales_mb16S == 0] <- min(prior_intercept_scales_mb16S[prior_intercept_scales_mb16S > 0],na.rm=TRUE)
prior_intercept_scales_mb16S[is.na(prior_intercept_scales_mb16S)] <- mean(prior_intercept_scales_mb16S, na.rm=TRUE)
for(x in 1:nrow(in_data$mb16S)) {
    wObs <- which(!is.na(inits_abundance$mb16S[x,]))
    for(y in which(is.na(inits_abundance$mb16S[x,]))) {
        closestRel <- wObs[which.min(cophenbact[y,wObs])]
        ypres <- in_data$mb16S[,y] > 0
        closepres <- in_data$mb16S[,closestRel] > 0
        bothpres <- ypres & closepres
        if(any(bothpres))
            inits_abundance$mb16S[x,y] <- inits_abundance$mb16S[x,closestRel] + mean(inits_abundance$mb16S[bothpres,y] - inits_abundance$mb16S[bothpres,closestRel], na.rm=TRUE)
        else
            inits_abundance$mb16S[x,y] <- inits_abundance$mb16S[x,closestRel] + mean(inits_abundance$mb16S[ypres,y],na.rm=TRUE) - mean(inits_abundance$mb16S[closepres,closestRel],na.rm=TRUE)
    }
}

mm_h$mb18S <- matrix(0, NNodesEuks + 1, NNodesEuks + 1)
for(node in 1:(NNodesEuks + 1)) {
    mm_h$mb18S[node, ] <- as.numeric(1:(NNodesEuks + 1) %in% c(Ancestors(euktreeY.root, node), node))
}
colnames(mm_h$mb18S) <- rownames(mm_h$mb18S) <- paste0('euk_i', 1:(NNodesEuks + 1))
colnames(mm_h$mb18S)[1:NTipsEuks] <- rownames(mm_h$mb18S)[1:NTipsEuks] <- euktreeY.root$tip.label
mm_h$mb18S <- mm_h$mb18S[colnames(in_data$mb18S), (NTipsEuks+2):ncol(mm_h$mb18S)]

copheneuk <- cophenetic(euktreeY.root)[colnames(in_data$mb18S),colnames(in_data$mb18S)]
inits_abundance$mb18S <- log(in_data$mb18S)
inits_abundance$mb18S[in_data$mb18S == 0] <- NA
inits_nuisance$mb18S <- numeric()
for(x in 1:nrow(in_data$mb18S)) {
    nuisance <- mean(inits_abundance$mb18S[x,],na.rm=TRUE)
    inits_nuisance$mb18S <- c(inits_nuisance$mb18S,nuisance)
    inits_abundance$mb18S[x,] <- inits_abundance$mb18S[x,] - nuisance
}
inits_intercepts$mb18S <- apply(inits_abundance$mb18S,2,mean,na.rm=TRUE)
inits_abundance$mb18S <- inits_abundance$mb18S - mean(inits_intercepts$mb18S)
inits_nuisance$mb18S <- inits_nuisance$mb18S + mean(inits_intercepts$mb18S)
inits_intercepts$mb18S <- inits_intercepts$mb18S - mean(inits_intercepts$mb18S)
prior_intercept_scales_mb18S <- apply(inits_abundance$mb18S,2,sd,na.rm=TRUE)
prior_intercept_scales_mb18S[prior_intercept_scales_mb18S == 0] <- min(prior_intercept_scales_mb18S[prior_intercept_scales_mb18S > 0],na.rm=TRUE)
prior_intercept_scales_mb18S[is.na(prior_intercept_scales_mb18S)] <- mean(prior_intercept_scales_mb18S,na.rm=TRUE)
for(x in 1:nrow(in_data$mb18S)) {
    wObs <- which(!is.na(inits_abundance$mb18S[x,]))
    for(y in which(is.na(inits_abundance$mb18S[x,]))) {
        closestRel <- wObs[which.min(copheneuk[y,wObs])]
        ypres <- in_data$mb18S[,y] > 0
        closepres <- in_data$mb18S[,closestRel] > 0
        bothpres <- ypres & closepres
        if(any(bothpres))
            inits_abundance$mb18S[x,y] <- inits_abundance$mb18S[x,closestRel] + mean(inits_abundance$mb18S[bothpres,y] - inits_abundance$mb18S[bothpres,closestRel], na.rm=TRUE)
        else
            inits_abundance$mb18S[x,y] <- inits_abundance$mb18S[x,closestRel] + mean(inits_abundance$mb18S[ypres,y],na.rm=TRUE) - mean(inits_abundance$mb18S[closepres,closestRel],na.rm=TRUE)
    }
}

inits_abundance$rna <- log(in_data$rna)
inits_abundance$rna[in_data$rna == 0] <- NA
inits_nuisance$rna <- numeric()
for(x in 1:nrow(in_data$rna)) {
    nuisance <- mean(inits_abundance$rna[x,],na.rm=TRUE)
    inits_nuisance$rna <- c(inits_nuisance$rna,nuisance)
    inits_abundance$rna[x,] <- inits_abundance$rna[x,] - nuisance
}
inits_intercepts$rna <- apply(inits_abundance$rna,2,mean,na.rm=TRUE)
inits_abundance$rna <- inits_abundance$rna - mean(inits_intercepts$rna)
inits_nuisance$rna <- inits_nuisance$rna + mean(inits_intercepts$rna)
inits_intercepts$rna <- inits_intercepts$rna - mean(inits_intercepts$rna)
prior_intercept_scales_rna <- apply(inits_abundance$rna,2,sd,na.rm=TRUE)
prior_intercept_scales_rna[prior_intercept_scales_rna == 0] <- min(prior_intercept_scales_rna[prior_intercept_scales_rna > 0],na.rm=TRUE)
prior_intercept_scales_rna[is.na(prior_intercept_scales_rna)] <- mean(prior_intercept_scales_rna, na.rm=TRUE)
for(x in 1:nrow(in_data$rna)) {
    for(y in 1:ncol(in_data$rna)) {
        if(is.na(inits_abundance$rna[x,y])) {
            inits_abundance$rna[x,y] <- mean(inits_abundance$rna[in_data$rna[,y] > 0, y],na.rm=TRUE)
        }
    }
}

allclades <- as.character(unique(itsMeta[colnames(in_data$its2), 'Clade']))
allmaj <- as.character(unique(itsMeta[colnames(in_data$its2), 'Majority.ITS2.sequence']))
mm_h$its2 <- cbind(sapply(unique(itsMeta[colnames(in_data$its2), 'Clade']), function(x) x == itsMeta[colnames(in_data$its2), 'Clade']),
                sapply(allmaj, function(x) x == itsMeta[colnames(in_data$its2), 'Majority.ITS2.sequence']))
colnames(mm_h$its2) <- paste0('ITS2_', c(allclades,allmaj))
mm_h$its2 <- mm_h$its2[,apply(mm_h$its2,2,function(x) sum(x) > 1)]

covits <- tcrossprod(mm_h$its2)
colnames(covits) <- rownames(covits) <- colnames(in_data$its2)
inits_abundance$its2 <- log(in_data$its2)
inits_abundance$its2[in_data$its2 == 0] <- NA
inits_nuisance$its2 <- numeric()
for(x in 1:nrow(in_data$its2)) {
    nuisance <- mean(inits_abundance$its2[x,],na.rm=TRUE)
    inits_nuisance$its2 <- c(inits_nuisance$its2,nuisance)
    inits_abundance$its2[x,] <- inits_abundance$its2[x,] - nuisance
}
inits_intercepts$its2 <- apply(inits_abundance$its2,2,mean,na.rm=TRUE)
inits_abundance$its2 <- inits_abundance$its2 - mean(inits_intercepts$its2)
inits_nuisance$its2 <- inits_nuisance$its2 + mean(inits_intercepts$its2)
inits_intercepts$its2 <- inits_intercepts$its2 - mean(inits_intercepts$its2)
prior_intercept_scales_its2 <- apply(inits_abundance$its2,2,sd,na.rm=TRUE)
prior_intercept_scales_its2[prior_intercept_scales_its2 == 0] <- min(prior_intercept_scales_its2[prior_intercept_scales_its2 > 0],na.rm=TRUE)
prior_intercept_scales_its2[is.na(prior_intercept_scales_its2)] <- mean(prior_intercept_scales_its2,na.rm=TRUE)
for(x in 1:nrow(in_data$its2)) {
    wObs <- which(!is.na(inits_abundance$its2[x,]))
    for(y in which(is.na(inits_abundance$its2[x,]))) {
        closestRel <- wObs[covits[y,wObs] == max(covits[y,wObs])]
        ypres <- in_data$its2[,y] > 0
        if(length(closestRel) == 1) {
            closepres <- in_data$its2[,closestRel] > 0
            bothpres <- ypres & closepres
            if(any(bothpres))
                inits_abundance$its2[x,y] <- inits_abundance$its2[x,closestRel] + mean(inits_abundance$its2[bothpres,y] - inits_abundance$its2[bothpres,closestRel], na.rm=TRUE)
            else
                inits_abundance$its2[x,y] <- inits_abundance$its2[x,closestRel] + mean(inits_abundance$its2[ypres,y],na.rm=TRUE) - mean(inits_abundance$its2[closepres,closestRel],na.rm=TRUE)
        } else {
            inits_abundance$its2[x,y] <- mean(inits_abundance$its2[x,closestRel]) + mean(inits_abundance$its2[ypres,y],na.rm=TRUE) - mean(sapply(closestRel, function(z) mean(inits_abundance$its2[in_data$its2[,z] > 0, z])), na.rm=TRUE)
        }
    }
}

names_mm_h$biomarkers <- c('symbiont_biomass','ubiquitin')
mm_h$biomarkers <- sapply(names_mm_h$biomarkers, function(x) grepl(x,colnames(in_data$biomarkers),ignore.case = TRUE))

mm_h$t2 <- cbind(1, as.numeric(grepl('raw', colnames(in_data$t2))), as.numeric(grepl('norm', colnames(in_data$t2))), as.numeric(grepl('Q1', colnames(in_data$t2))), as.numeric(grepl('Q3', colnames(in_data$t2))))
colnames(mm_h$t2) <- c('T2_all', 'T2_raw_all', 'T2_norm_all', 'T2_Q1', 'T2_Q3')
mm_h$t3 <- cbind(1, as.numeric(grepl('raw', colnames(in_data$t3))), as.numeric(grepl('norm', colnames(in_data$t3))), as.numeric(grepl('Q1', colnames(in_data$t3))), as.numeric(grepl('Q3', colnames(in_data$t3))))
colnames(mm_h$t3) <- c('T3_all', 'T3_raw_all', 'T3_norm_all', 'T3_Q1', 'T3_Q3')

names_mm_h$fab_t1_agg <- list(windspeed_group = 'speed',
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
names_mm_h$fab_t1_agg <- names_mm_h$fab_t1_agg[sapply(names_mm_h$fab_t1_agg, function(x) sum(grepl(x,colnames(in_data$fab_T1_agg),ignore.case = TRUE))) > 1]
mm_h$fab_T1_agg <- sapply(names_mm_h$fab_t1_agg, function(x) grepl(x,colnames(in_data$fab_T1_agg),ignore.case = TRUE))

names_mm_h$hphoto <- c('morphotype_Millepora','morphotype_Pocillopora','morphotype_Porites')
mm_h$hphoto <- sapply(names_mm_h$hphoto, function(x) grepl(x,colnames(in_data$hphoto),ignore.case = TRUE))
prior_scales_hphoto <- apply(matrix(as.numeric(in_data$hphoto > 0),ncol=ncol(in_data$hphoto)) %*% t(MASS:::ginv(cbind(1,diag(ncol(in_data$hphoto)),mm_h$hphoto))[-1,]),2,sd,na.rm=TRUE)
prior_scales_hphoto[prior_scales_hphoto==0] <- min(c(prior_scales_hphoto[prior_scales_hphoto>0],1)) / 2

names_mm_h$ephoto <- list(bleached_group='bleached_light|bleached_bleached', borers_group='bivalve|Spirobranchus|Tridacna|boringurchin|other_polychaete|sponge|gallcrabs|ascidians', polychaetes_group='Spirobranchus|other_polychaete', bivalve_group='bivalve|Tridacna', algae_contact_group='Halimeda|Turbinaria|Dictyota|Lobophora|Galaxaura|Sargassum', sargassaceae_group='Turbinaria|Sargassum')
names_mm_h$ephoto <- names_mm_h$ephoto[sapply(names_mm_h$ephoto, function(x) sum(grepl(x,colnames(in_data$ephoto),ignore.case = TRUE))) > 1]
mm_h$ephoto <- sapply(names_mm_h$ephoto, function(x) grepl(x,colnames(in_data$ephoto),ignore.case = TRUE))
prior_scales_ephoto <- apply(matrix(as.numeric(in_data$ephoto > 0),ncol=ncol(in_data$ephoto)) %*% t(MASS:::ginv(cbind(1,diag(ncol(in_data$ephoto)),mm_h$ephoto))[-1,]),2,sd,na.rm=TRUE)
prior_scales_ephoto[prior_scales_ephoto==0] <- min(c(prior_scales_ephoto[prior_scales_ephoto>0],1)) / 2

names_mm_h$svd <- c('cladePOC','cladePOR')
mm_h$svd <- sapply(names_mm_h$svd, function(x) grepl(x,colnames(in_data$svd),ignore.case = TRUE))
prior_scales_svd <- apply(matrix(as.numeric(in_data$svd > 0),ncol=ncol(in_data$svd)) %*% t(MASS:::ginv(cbind(1,diag(ncol(in_data$svd)),mm_h$svd))[-1,]),2,sd,na.rm=TRUE)
prior_scales_svd[prior_scales_svd==0] <- min(c(prior_scales_svd[prior_scales_svd>0],1)) / 2

names_mm_h$sites <- levels(filtData$island)
mm_h$sites <- sapply(names_mm_h$sites, function(x) grepl(x,colnames(in_data$sites),ignore.case = TRUE))
prior_scales_sites <- apply(matrix(as.numeric(in_data$sites > 0),ncol=ncol(in_data$sites)) %*% t(MASS:::ginv(cbind(1,diag(ncol(in_data$sites)),mm_h$sites))[-1,]),2,sd,na.rm=TRUE)
prior_scales_sites[prior_scales_sites==0] <- min(c(prior_scales_sites[prior_scales_sites>0],1)) / 2

intercept_fun <- function(x) {
    if(sum(x,na.rm=TRUE) != sum(!is.na(x)))
        mean(x,na.rm=TRUE)
    else
        sum(x,na.rm=TRUE)/(sum(!is.na(x))+1)
}
binary_count_intercepts_inits_log_list <- sapply(names(in_data)[1:D], function(d) log(apply(in_data[[d]] > 0, 2, intercept_fun)))
binary_count_intercepts_inits <- logit(exp(unlist(binary_count_intercepts_inits_log_list)))

ii_fun <- function(x,M,mm) {
    if(M[[x]] == 1) {
        logit(mean(mm[,sum(M[0:x])], na.rm=TRUE))
    } else {
        temp <- log(apply(mm[,(sum(M[0:(x-1)])+1):sum(M[0:x])],
                          2,
                          mean,
                          na.rm=TRUE))
        return(temp - mean(temp))
    }
}
intercepts_inits <- c(inits_intercepts$mb16S,
                      inits_intercepts$mb18S,
                      inits_intercepts$rna,
                      inits_intercepts$its2,
                      apply(inits_biomarkers,2,mean,na.rm=TRUE),
                      apply(in_data$t2,2,mean,na.rm=TRUE),
                      apply(in_data$t3,2,mean,na.rm=TRUE),
                      apply(in_data$fab_T1_agg,2,mean,na.rm=TRUE),
                      apply(in_data$fab_T2_agg,2,mean,na.rm=TRUE),
                      unlist(sapply(1:length(M_snp), ii_fun, M_snp, in_data$snps)),
                      unlist(sapply(1:length(M_hphoto), ii_fun, M_hphoto, in_data$hphoto)),
                      unlist(sapply(1:length(M_ephoto), ii_fun, M_ephoto, in_data$ephoto)),
                      unlist(sapply(1:length(M_spec), ii_fun, M_spec, in_data$species)),
                      unlist(sapply(1:length(M_svd), ii_fun, M_svd, in_data$svd)),
                      {temp <- log(apply(in_data$sites,2,mean,na.rm=TRUE)); temp - mean(temp)},
                      binary_count_intercepts_inits,
                      rep(0,4))

prior_binary_count_intercept_scales <- sapply(names(in_data)[1:D], function(d) {out <- apply(in_data[[d]] > 0, 2, sd); out[out==0] <- min(c(out[out>0],1)) / 2; return(out)})
prior_binary_count_dataset_intercept_scales <- sapply(names(prior_binary_count_intercept_scales), function(d) exp(mean(log(prior_binary_count_intercept_scales[[d]]))))

multinomial_nuisance_inits <- c(inits_nuisance$mb16S,
                                inits_nuisance$mb18S,
                                inits_nuisance$rna,
                                inits_nuisance$its2)


biomarkersInv <- t(MASS:::ginv(cbind(diag(ncol(inits_biomarkers)),mm_h$biomarkers)))
fabT1aggMatInv <- t(MASS:::ginv(cbind(diag(ncol(in_data$fab_T1_agg)),mm_h$fab_T1_agg)))

inits_abundance_solved <- lapply(names(inits_abundance), function(d) {
    inits_abundance[[d]] %*% t(MASS:::ginv(cbind(1,diag(ncol(inits_abundance[[d]])),mm_h[[d]]))[-1,])
})
names(inits_abundance_solved) <- names(inits_abundance)

inits_prevalence_solved <- lapply(names(inits_abundance), function(d) {
    matrix(as.numeric(in_data[[d]] > 0),ncol=ncol(in_data[[d]])) %*% t(MASS:::ginv(cbind(1,diag(ncol(inits_abundance[[d]])),mm_h[[d]]))[-1,])
})
names(inits_prevalence_solved) <- names(inits_abundance)

prior_scales_abundance <- sapply(names(inits_abundance_solved), function(d) apply(inits_abundance_solved[[d]], 2, sd))
prior_scales_prevalence <- sapply(names(inits_abundance_solved), function(d) {out <- apply(inits_prevalence_solved[[d]], 2, sd); out[out==0] <- min(c(out[out>0],1)) / 2; return(out)})

prior_scales <- unlist(c(prior_scales_abundance,
                         apply(sapply(1:ncol(biomarkersInv), function(x) matrix(inits_biomarkers[,biomarkersInv[,x]!=0],nrow=nrow(inits_biomarkers)) %*% biomarkersInv[biomarkersInv[,x]!=0,x]),2,sd,na.rm=TRUE),
                         apply(in_data$t2 %*% t(MASS:::ginv(cbind(diag(ncol(in_data$t2)),mm_h$t2))),2,sd),
                         apply(in_data$t3 %*% t(MASS:::ginv(cbind(diag(ncol(in_data$t3)),mm_h$t3))),2,sd),
                         apply(sapply(1:ncol(fabT1aggMatInv), function(x) matrix(in_data$fab_T1_agg[,fabT1aggMatInv[,x]!=0],nrow=nrow(in_data$fab_T1_agg)) %*% fabT1aggMatInv[fabT1aggMatInv[,x]!=0,x]),2,sd,na.rm=TRUE),
                         apply(in_data$fab_T2_agg,2,sd,na.rm=TRUE),
                         apply(in_data$snps,2,sd,na.rm=TRUE),
                         prior_scales_hphoto,
                         prior_scales_ephoto,
                         apply(in_data$species,2,sd,na.rm=TRUE),
                         prior_scales_svd,
                         prior_scales_sites,
                         prior_scales_prevalence,
                         prior_binary_count_dataset_intercept_scales))

prior_intercept_scales <- c(prior_intercept_scales_mb16S,
                            prior_intercept_scales_mb18S,
                            prior_intercept_scales_rna,
                            prior_intercept_scales_its2,
                            apply(inits_biomarkers,2,sd,na.rm=TRUE),
                            apply(in_data$t2,2,sd),
                            apply(in_data$t3,2,sd),
                            apply(in_data$fab_T1_agg,2,sd,na.rm=TRUE),
                            apply(in_data$fab_T2_agg,2,sd,na.rm=TRUE),
                            apply(in_data$snps,2,sd,na.rm=TRUE),
                            apply(in_data$hphoto,2,sd,na.rm=TRUE),
                            apply(in_data$ephoto,2,sd,na.rm=TRUE),
                            apply(in_data$species,2,sd,na.rm=TRUE),
                            apply(in_data$svd,2,sd,na.rm=TRUE),
                            apply(in_data$sites,2,sd,na.rm=TRUE),
                            prior_binary_count_intercept_scales,
                            prior_binary_count_dataset_intercept_scales)

inv_log_max_contam <- 1 / log(sapply(names(in_data)[1:D], function(d) max(apply(diag(1/rowSums(in_data[[d]])) %*% in_data[[d]], 2, function(x) max(x[x>0]) / min(x[x>0])))))

inits_contamination <- lapply(names(inits_abundance), function(d) {
    y <- sapply(1:nrow(inits_abundance[[d]]), function(x) {
        z <- inits_intercepts[[d]] + binary_count_intercepts_inits_log_list[[d]] - 1 / inv_log_max_contam[d]
        return(z)
    })
    return(t(y))
})
names(inits_contamination) <- names(inits_abundance)

inits_counts <- lapply(names(inits_abundance), function(d) {
    sapply(1:ncol(inits_abundance[[d]]), function(y) {
        sapply(1:nrow(inits_abundance[[d]]), function(x) {
            if(in_data[[d]][x,y] > 0) {
                return(inits_abundance[[d]][x,y])
            } else {
                return(inits_contamination[[d]][x,y] / 4 + inits_abundance[[d]][x,y] / 4 + (log(0.5 / sum(in_data[[d]][x,])) - inits_nuisance[[d]][x]) / 2)
            }
        }) # initiate with log of observed counts if positive, or split the difference between mean of log of observed counts and initial contamination estimate and a count of 1 if zero
    })
})
names(inits_counts) <- names(inits_abundance)

scales_observed <- sapply(names(inits_counts), function(d) apply(inits_counts[[d]], 2, sd))

prior_intercept_centers <- intercepts_inits

M_higher <- c(ncol(mm_h$mb16S), ncol(mm_h$mb18S), 0, ncol(mm_h$its2), ncol(mm_h$biomarkers), ncol(mm_h$t2), ncol(mm_h$t3), ncol(mm_h$fab_T1_agg), 0, 0, ncol(mm_h$hphoto), ncol(mm_h$ephoto), 0, ncol(mm_h$svd), ncol(mm_h$sites))
M_all <- M + M_higher
VOBplus <- sum(M_all)
mm <- unlist(mm_h)
size_mm <- length(mm)

F_higher <- sapply(1:D, function(x) sum(I[x,]*(M_all[x]-M[x])))

IR_higher <- rbind(t(mm_h$biomarkers) %*% IR[1:M[D+1],],
                   t(mm_h$t2) %*% IR[(M[D+1] + 1):sum(M[(D+1):(D+2)]),],
                   t(mm_h$t3) %*% IR[(sum(M[(D+1):(D+2)]) + 1):sum(M[(D+1):(D+3)]),],
                   t(mm_h$fab_T1_agg) %*% IR[(sum(M[(D+1):(D+3)]) + 1):sum(M[(D+1):(D+4)]),])
IR_higher[IR_higher > 0] <- 1
O_higher = nrow(IR_higher)
H_higher <- sum(IR_higher > 0)

IC_higher <- rbind(t(mm_h$hphoto) %*% IC[1:M[D+R+2],],
                   t(mm_h$ephoto) %*% IC[(M[D+R+2] + 1):sum(M[(D+R+2):(D+R+3)]),],
                   t(mm_h$svd) %*% IC[(sum(M[(D+R+2):(D+R+4)]) + 1):sum(M[(D+R+2):(D+R+5)]),],
                   t(mm_h$sites) %*% IC[(sum(M[(D+R+2):(D+R+5)]) + 1):sum(M[(D+R+2):(D+R+6)]),])
IC_higher[IC_higher > 0] <- 1
B_higher <-nrow(IC_higher)
G_higher <- sapply(1:length(Mc), function(x) sum(ICv[x,]*Mc[x]))

varlabs <- c(colnames(in_data$mb16S), colnames(mm_h$mb16S), colnames(in_data$mb18S), colnames(mm_h$mb18S), colnames(in_data$rna), colnames(in_data$its2), colnames(mm_h$its2), colnames(in_data$biomarkers), names_mm_h$biomarkers, colnames(in_data$t2), colnames(mm_h$t2), colnames(in_data$t3), colnames(mm_h$t3), colnames(in_data$fab_T1_agg), names(names_mm_h$fab_t1_agg), colnames(in_data$fab_T2_agg), colnames(in_data$snps), colnames(in_data$hphoto), names_mm_h$hphoto, colnames(in_data$ephoto), names(names_mm_h$ephoto), colnames(in_data$species), colnames(in_data$svd), colnames(mm_h$svd), colnames(in_data$sites), colnames(mm_h$sites))
varlabsM <- unlist(sapply(in_data, colnames))

data <- list(N            = N,
             N_var_groups = N_var_groups,
             N_all        = N_all,
             D = D,
             R = R,
             C = C,
             M = M,
             M_higher = M_higher,
             M_all = M_all,
             I = I,
             O = O,
             O_higher = O_higher,
             IR = IR,
             IR_higher = IR_higher,
             IC = IC,
             B_higher = B_higher,
             IC_higher = IC_higher,
             ICv = ICv,
             F = F,
             F_higher = sum(F_higher),
             X = X,
             G = G,
             Y = Y,
             H = H,
             H_higher = H_higher,
             P = P,
             B = B,
             C_vars = C_vars,
             Mc = Mc,
             G_higher = sum(G_higher),
             scales_observed                = unlist(scales_observed) * mass_slow,
             prior_scales                   = prior_scales,
             prior_intercept_centers        = prior_intercept_centers,
             size_mm = size_mm,
             mm     = mm,
             global_scale_prior = global_scale_prior,
             K_linear = K_linear,
             K_gp     = K_gp,
             KG       = KG,
             N_Pm     = N_Pm,
             idx_Pm   = idx_Pm,
             idx_Pnm  = idx_Pnm,
             dist_sites          = dist_sites[lower.tri(dist_sites)],
             rho_sites_prior     = mean(dist_sites[lower.tri(dist_sites)]),
             samp2group          = samp2group,
             site_smoothness     = site_smoothness,
             nu_residuals        = nu_residuals,
             inv_log_max_contam  = inv_log_max_contam,
             shape_gnorm         = shape_gnorm,
             mass_slow           = mass_slow)

#### create initializations
inits_counts_norm <- lapply(names(inits_counts), function(d) {
    inits_counts[[d]] %*% diag(1 / scales_observed[[d]])
})
names(inits_counts_norm) <- names(inits_counts)

abundance_observed_vector_inits <- unlist(c(sapply(1:N, function(x) if(I[1,x]) inits_counts_norm$mb16S[I_cs[1,x],]),
                                            sapply(1:N, function(x) if(I[2,x]) inits_counts_norm$mb18S[I_cs[2,x],]),
                                            sapply(1:N, function(x) if(I[3,x]) inits_counts_norm$rna[I_cs[3,x],]),
                                            sapply(1:N, function(x) if(I[4,x]) inits_counts_norm$its2[I_cs[4,x],])))

inits_counts_solved <- lapply(names(inits_counts), function(d) {
    inits_counts[[d]] %*% t(MASS:::ginv(cbind(1,diag(ncol(inits_counts[[d]])),mm_h[[d]]))[-1,]) %*% diag(1 / prior_scales_abundance[[d]])
})
names(inits_counts_solved) <- names(inits_counts)

abundance_higher_vector_inits <- unlist(c(sapply(1:N, function(x) if(I[1,x]) inits_counts_solved$mb16S[I_cs[1,x],-(1:M[[1]])]),
                                          sapply(1:N, function(x) if(I[2,x]) inits_counts_solved$mb18S[I_cs[2,x],-(1:M[[2]])]),
                                          sapply(1:N, function(x) if(I[3,x]) inits_counts_solved$rna[I_cs[3,x],-(1:M[[3]])]),
                                          sapply(1:N, function(x) if(I[4,x]) inits_counts_solved$its2[I_cs[4,x],-(1:M[[4]])])))

inits_prevalence_solved_norm <- lapply(names(inits_counts), function(d) {
    inits_prevalence_solved[[d]] %*% diag(1 / prior_scales_prevalence[[d]])
})
names(inits_prevalence_solved_norm) <- names(inits_prevalence_solved)

prevalence_higher_vector_inits <- unlist(c(sapply(1:N, function(x) if(I[1,x]) inits_prevalence_solved_norm$mb16S[I_cs[1,x],-(1:M[[1]])]),
                                           sapply(1:N, function(x) if(I[2,x]) inits_prevalence_solved_norm$mb18S[I_cs[2,x],-(1:M[[2]])]),
                                           sapply(1:N, function(x) if(I[3,x]) inits_prevalence_solved_norm$rna[I_cs[3,x],-(1:M[[3]])]),
                                           sapply(1:N, function(x) if(I[4,x]) inits_prevalence_solved_norm$its2[I_cs[4,x],-(1:M[[4]])])))

inits_P_higher <- list(biomarkers = {temp <- inits_biomarkers %*% t(MASS:::ginv(cbind(1,diag(ncol(inits_biomarkers)),mm_h$biomarkers)))[,-(1:(M['biomarkers']+1))]; for(j in 1:ncol(temp)) temp[is.na(temp[,j]),j] <- mean(temp[,j],na.rm=TRUE); temp},
                       t2 = in_data$t2 %*% t(MASS:::ginv(cbind(1,diag(ncol(in_data$t2)),mm_h$t2)))[,-(1:(M['t2']+1))],
                       t3 = in_data$t3 %*% t(MASS:::ginv(cbind(1,diag(ncol(in_data$t3)),mm_h$t3)))[,-(1:(M['t3']+1))],
                       fab_T1_agg = {temp <- in_data$fab_T1_agg %*% t(MASS:::ginv(cbind(1,diag(ncol(in_data$fab_T1_agg)),mm_h$fab_T1_agg)))[,-(1:(M['fab_T1_agg']+1))]; for(j in 1:ncol(temp)) temp[is.na(temp[,j]),j] <- mean(temp[,j],na.rm=TRUE); temp})

inits_P_higher <- unlist(sapply(1:length(inits_P_higher), function(r) {
    sapply(1:N_all, function(n) {
        sapply(1:M_higher[D+r], function(m) {
            if(IR_higher[sum(M_higher[1:(D+r-1)]) - sum(M_higher[1:D]) + m, n]) mean(inits_P_higher[[r]][,m])
        })
    })
}))

inits_Y_higher <- list(matrix(as.numeric(in_data$hphoto > 0),ncol=ncol(in_data$hphoto)) %*% t(MASS:::ginv(cbind(1,diag(ncol(in_data$hphoto)),mm_h$hphoto))[-(1:(M['hphoto']+1)),]),
                       matrix(as.numeric(in_data$ephoto > 0),ncol=ncol(in_data$ephoto)) %*% t(MASS:::ginv(cbind(1,diag(ncol(in_data$ephoto)),mm_h$ephoto))[-(1:(M['ephoto']+1)),]))

Z_init <- matrix(rnorm((K_linear+KG*K_gp)*N), nrow=K_linear+KG*K_gp) * 1e-5

latent_var_raw_init <- c(2*log(global_scale_prior), rep(0, K_linear))
if(KG > 0) {
    latent_var_raw_init <- c(latent_var_raw_init, rep(0, K_gp))
}
if(KG > 1) {
    for(g in 2:KG) {
        latent_var_raw_init <- c(latent_var_raw_init, rep(0, K_gp))
    }
}

init <- list(abundance_observed_vector       = abundance_observed_vector_inits / mass_slow,
             intercepts                      = prior_intercept_centers,
             multinomial_nuisance            = multinomial_nuisance_inits,
             latent_var_raw      = latent_var_raw_init,
             weight_scales_raw  = matrix(1 / mass_slow, nrow=2*D+R+C, ncol=K),
             sds_raw            = rep(0, VOBplus+sum(M_all[1:D])+D),
             dataset_scales_raw = rep(0, 2*D+R+C),
             rho_sites          = as.array(rep(mean(dist_sites[lower.tri(dist_sites)]), K)),
             site_prop          = as.array(rep(0.5, K)),
             abundance_higher_vector  = abundance_higher_vector_inits,
             prevalence_higher_vector = prevalence_higher_vector_inits,
             P_higher        = inits_P_higher,
             Y_higher_vector = rep(0,sum(G_higher)),
             Z_linear_raw    = Z_init[1:K_linear,],
             Z_gp_raw        = {if(K > K_linear) Z_init[(K_linear+1):K,] else 1},
             W_raw     = matrix(0, nrow=VOBplus+sum(M_all[1:D])+D, ncol=K),
             rho_Z     = matrix(2 / K_gp, nrow = K_linear, ncol = KG),
             inv_log_less_contamination  = -inv_log_max_contam,
             contaminant_overdisp_raw    = rep(0,D),
             alpha_Z_raw                 = rep(0,K),
             beta_Z_prop_raw             = rep(-1 / mass_slow,K))

save.image(file.path(output_prefix, 'setup.RData'))

write_stan_json(data, file.path(output_prefix, 'data.json'))
write_stan_json(init, file.path(output_prefix, 'inits.json'))

setwd(cmdstan_path())
system(paste0(c('make ', 'make STAN_OPENCL=true ')[opencl+1], 'STANCFLAGS="--include-paths=', include_path, '" ', file.path(model_dir, model_name)))

setwd(model_dir)
print(sampling_commands[[engine]])
print(date())
system(sampling_commands[[engine]])

importparams <- c('W_norm','Z','sds','latent_scales','dataset_scales','var_scales','weight_scales','rho_sites', 'site_prop', 'corr_sites', 'intercepts', 'log_less_contamination', 'contaminant_overdisp', 'rho_Z'[if(KG > 0) TRUE else NULL],'alpha_Z','beta_Z_prop')

filter_large_stan_csv(file.path(output_prefix, paste0('samples_',engine,'.csv')), file.path(output_prefix, paste0('samples_',engine,'_filtered.csv')), importparams)

stan.fit <- read_cmdstan_csv(file.path(output_prefix, paste0('samples_',engine,'_filtered.csv')),
                             variables = importparams,
                             format = 'draws_array')

save.image(file.path(output_prefix, paste0('res_',engine,'.RData')))

source(file.path(model_dir, 'BGPFA_extract.r'))

save.image(file.path(output_prefix, paste0('res_',engine,'_processed.RData')))
