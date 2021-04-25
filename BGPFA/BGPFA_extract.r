library(rstan)
library(phangorn)
library(RColorBrewer)
source(file.path(model_dir, 'my_triplot.r'))

monteCarloP <- function(x, pn='p') {
    if(pn == 'n') {
        res <- (1 + sum(x >= 0, na.rm=TRUE)) / (1 + length(x))
    } else if(pn == 'p') {
        res <- (1 + sum(x <= 0, na.rm=TRUE)) / (1 + length(x))
    }
} # https://arxiv.org/pdf/1603.05766.pdf

sumfunc <- median
#sumfunc <- mean
#sumfunc <- function(x) x[length(x)]
#sumfunc <- function(x) x[maxind]

labs <- c(varlabs, paste0(varlabs[1:sum(M_all[1:D])],'.binary'), paste0('binary_count_dataset_int_',1:D))

dataset_scales <- extract(stan.fit, pars='dataset_scales', permuted=FALSE)
dataset_scales <- apply(dataset_scales,3,sumfunc)

var_scales <- extract(stan.fit, pars='var_scales', permuted=FALSE)
var_scales <- apply(var_scales,3,sumfunc)
names(var_scales) <- labs

global_effect_scale <- extract(stan.fit, pars='global_effect_scale', permuted=FALSE)
global_effect_scale <- apply(global_effect_scale, 3, sumfunc)

latent_scales <- extract(stan.fit, pars='latent_scales', permuted=FALSE)
latent_scales <- apply(latent_scales, 3, sumfunc)
axisOrder <- order(latent_scales, decreasing=TRUE)
latent_scales <- latent_scales[axisOrder]

weight_scales <- extract(stan.fit, pars='weight_scales', permuted=FALSE)
weight_scales <- apply(weight_scales,3,sumfunc)
weight_scales <- array(weight_scales,dim=c(length(weight_scales)/K,K))
weight_scales <- array(weight_scales[,axisOrder], dim=c(length(weight_scales)/K,K))

if('sds' %in% importparams) {
    sds <- extract(stan.fit, pars='sds', permuted=FALSE)
    sds <- apply(sds,3,sumfunc)
}

nu_factors <- extract(stan.fit, pars='nu_factors', permuted=FALSE)
nu_factors <- matrix(apply(nu_factors, 3, sumfunc), ncol=K)
nu_factors <- nu_factors[,axisOrder]

if('log_less_contamination' %in% importparams) {
    log_less_contamination <- extract(stan.fit, pars='log_less_contamination', permuted=FALSE)
    log_less_contamination <- apply(log_less_contamination, 3, sumfunc)
}

if('contaminant_overdisp' %in% importparams) {
    contaminant_overdisp <- extract(stan.fit, pars='contaminant_overdisp', permuted=FALSE)
    contaminant_overdisp <- apply(contaminant_overdisp, 3, sumfunc)
}

Z_all <- extract(stan.fit, pars='Z', permuted=FALSE)

Z <- apply(Z_all,3,sumfunc)
Z <- t(array(Z,dim=c(K,length(Z)/K)))
Z <- array(Z[,axisOrder],dim=c(length(Z)/K,K))
rownames(Z) <- allsamples

Z_all <- array(Z_all, dim=c(dim(Z_all)[[1]],N,K))[,,axisOrder,drop=FALSE]

W_norm_all <- extract(stan.fit, pars='W_norm', permuted=FALSE)
W_norm <- apply(W_norm_all,3,sumfunc)
W_norm <- array(W_norm,dim=c(length(W_norm)/K,K))[,axisOrder]

binary_count_dataset_intercepts <- extract(stan.fit, pars='binary_count_dataset_intercepts', permuted=FALSE)
binary_count_dataset_intercepts <- apply(binary_count_dataset_intercepts, 3, sumfunc)

DRC <- D+R+C

rho_sites <- extract(stan.fit, pars='rho_sites', permuted=FALSE)
rho_sites <- apply(rho_sites, 3, sumfunc)

cov_sites <- extract(stan.fit, pars='cov_sites', permuted=FALSE)
cov_sites <- apply(cov_sites, 3, sumfunc)
dim(cov_sites) <- c(K,N_sites,N_sites)
cov_sites <- cov_sites[axisOrder,,]

if('rho_Z' %in% importparams) {
    rho_Z <- extract(stan.fit, pars='rho_Z', permuted=FALSE)
    if(exists('KG')) {
        rho_Z <- array(apply(rho_Z, 3, sumfunc),dim=c(K_linear,KG))[axisOrder[axisOrder <= K_linear],]
    } else {
        rho_Z <- apply(rho_Z, 3, sumfunc)[axisOrder[axisOrder <= K_linear]]
    }
}

pos1 <- array(apply(W_norm_all,3,monteCarloP,pn='p'),dim=c(dim(W_norm_all)[[3]]/K,K))[,axisOrder]
possig <- pos1 < 0.05
neg1 <- array(apply(W_norm_all,3,monteCarloP,pn='n'),dim=c(dim(W_norm_all)[[3]]/K,K))[,axisOrder]
negsig <- neg1 < 0.05
anysig <- possig | negsig
rownames(pos1) <- rownames(neg1) <- labs
gc()

nullfunc <- function() {

    gc()

    drivers <- mytriplot(Z, W_norm, Z, 1,2, as.factor(filtData[allsamples,]$species), labs, 50, TRUE, NULL, NULL, FALSE, TRUE)

    drivers <- mytriplot(Z, W_norm, Z, 1,2, as.factor(filtData[allsamples,]$species), labs, 50, TRUE, NULL, NULL, FALSE, TRUE, anysig)

    i <- 1
    hist((var_scales / prior_scales / dataset_scales[i])[(sum(M_all[1:(i-1)])+1):sum(M_all[1:i])])

    labs[negsig[,1]]
    labs[possig[,1]]
    list(positive=labs[possig[,1]],negative=labs[negsig[,1]])

}
