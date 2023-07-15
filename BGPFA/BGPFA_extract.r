library(rstan)
library(phangorn)
library(RColorBrewer)
model_dir <- file.path(Sys.getenv('HOME'), 'scripts/stan_models/BGPFA/')
source(file.path(model_dir, 'functions.r'))

if(length(stan.fit) == 2) {
    stan.fit.draws <- stan.fit[[2]]
} else {
    stan.fit.draws <- stan.fit$post_warmup_draws
}

sumfunc <- median
#sumfunc <- mean
#sumfunc <- function(x) x[length(x)]
#sumfunc <- function(x) x[maxind]

labs <- c(varlabs, paste0(varlabs[1:sum(M_all[1:D])],'.binary'), paste0('binary_count_dataset_int_',dataset_names[1:D]))
dataset_names_expanded <- c(dataset_names, paste0(dataset_names[1:D],'.binary'))

dataset_scales <- stan.fit.draws[,,grep('dataset_scales\\[.*',dimnames(stan.fit.draws)[[3]])]
dataset_scales <- apply(dataset_scales,3,sumfunc)
names(dataset_scales) <- dataset_names_expanded

var_scales <- stan.fit.draws[,,grep('var_scales\\[.*',dimnames(stan.fit.draws)[[3]])]
var_scales <- apply(var_scales,3,sumfunc)
names(var_scales) <- labs

global_effect_scale <- stan.fit.draws[,,grep('global_effect_scale',dimnames(stan.fit.draws)[[3]])]
global_effect_scale <- apply(global_effect_scale, 3, sumfunc)

W_norm_all <- stan.fit.draws[,,grep('W_norm\\[.*',dimnames(stan.fit.draws)[[3]]),drop = FALSE]
W_norm_all <- array(W_norm_all, dim = c(dim(W_norm_all)[[1]], dim(W_norm_all)[[3]]/K, K))
Z_all <- stan.fit.draws[,,grep('Z\\[.*',dimnames(stan.fit.draws)[[3]])]
Z_all <- aperm(array(Z_all, dim = c(dim(Z_all)[[1]], K, N)), c(1,3,2))
WZ_all <- abind:::abind(W_norm_all,Z_all,along=2)
WZ <- apply(WZ_all,c(2,3),sumfunc,na.rm=TRUE)
W_norm <- WZ[1:dim(W_norm_all)[2],]
Z <- WZ[(dim(W_norm_all)[2]+1):nrow(WZ),]
rownames(Z) <- allsamples

latent_scales <- stan.fit.draws[,,grep('latent_scales\\[.*',dimnames(stan.fit.draws)[[3]])]
latent_scales <- apply(latent_scales, 3, sumfunc)

weight_scales <- stan.fit.draws[,,grep('weight_scales\\[.*',dimnames(stan.fit.draws)[[3]])]
weight_scales <- apply(weight_scales,3,sumfunc)
weight_scales <- array(weight_scales,dim=c(length(weight_scales)/K,K))
weight_scales <- array(weight_scales, dim=c(length(weight_scales)/K,K))
rownames(weight_scales) <- dataset_names_expanded

if('sds' %in% importparams) {
    sds <- stan.fit.draws[,,grep('sds\\[.*',dimnames(stan.fit.draws)[[3]])]
    sds <- apply(sds,3,sumfunc)
}

if('log_less_contamination' %in% importparams) {
    log_less_contamination <- stan.fit.draws[,,grep('log_less_contamination\\[.*',dimnames(stan.fit.draws)[[3]])]
    log_less_contamination <- apply(log_less_contamination, 3, sumfunc)
}

if('contaminant_overdisp' %in% importparams) {
    contaminant_overdisp <- stan.fit.draws[,,grep('contaminant_overdisp\\[.*',dimnames(stan.fit.draws)[[3]])]
    contaminant_overdisp <- apply(contaminant_overdisp, 3, sumfunc)
}

if('binary_count_dataset_intercepts' %in% importparams) {
    binary_count_dataset_intercepts <- stan.fit.draws[,,grep('binary_count_dataset_intercepts\\[.*',dimnames(stan.fit.draws)[[3]])]
    binary_count_dataset_intercepts <- apply(binary_count_dataset_intercepts, 3, sumfunc)
    names(binary_count_dataset_intercepts) <- dataset_names[1:D]
}

DRC <- D+R+C

rho_sites <- stan.fit.draws[,,grep('rho_sites',dimnames(stan.fit.draws)[[3]])]
rho_sites <- apply(rho_sites, 3, sumfunc)

if('corr_sites' %in% importparams) {
    corr_sites <- stan.fit.draws[,,grep('corr_sites\\[.*',dimnames(stan.fit.draws)[[3]])]
    corr_sites <- apply(corr_sites, 3, sumfunc)
    dim(corr_sites) <- c(K,N_sites,N_sites)
}

if('rho_Z' %in% importparams & KG > 0) {
    rho_Z <- stan.fit.draws[,,grep('rho_Z\\[.*',dimnames(stan.fit.draws)[[3]])]
    rho_Z <- array(apply(rho_Z, 3, sumfunc),dim=c(K_linear,KG))
}

if('skew_Z' %in% importparams) {
    skew_Z <- stan.fit.draws[,,grep('skew_Z\\[.*',dimnames(stan.fit.draws)[[3]])]
    skew_Z <- apply(skew_Z, 3, sumfunc)
}

if('alpha_Z' %in% importparams) {
    alpha_Z <- stan.fit.draws[,,grep('alpha_Z\\[.*',dimnames(stan.fit.draws)[[3]]), drop=FALSE]
    alpha_Z <- apply(alpha_Z, 3, sumfunc)
}

if('beta_Z_prop' %in% importparams) {
    beta_Z_prop <- stan.fit.draws[,,grep('beta_Z_prop\\[.*',dimnames(stan.fit.draws)[[3]]), drop=FALSE]
    beta_Z_prop <- apply(beta_Z_prop, 3, sumfunc)
}

pos1 <- apply(WZ_all,c(2,3),monteCarloP,pn='p')
possig <- pos1 < 0.05
neg1 <- apply(WZ_all,c(2,3),monteCarloP,pn='n')
negsig <- neg1 < 0.05
anysig <- possig | negsig
rownames(pos1) <- rownames(neg1) <- c(labs,allsamples)
gc()

nullfunc <- function() {

    gc()

    drivers <- mytriplot(Z, W_norm, Z, 1,2, as.factor(filtData[allsamples,]$species), labs, 50, TRUE, NULL, NULL, FALSE, TRUE, anysig[1:nrow(W_norm),])

    mshape_pcs <- princomp(WZ)
    W_pc <- mshape_pcs$scores[1:nrow(W_norm),]
    Z_pc <- mshape_pcs$scores[(dim(W_norm_all)[2]+1):nrow(WZ),]
    drivers <- mytriplot(Z_pc, W_pc, Z_pc, 1,2, as.factor(filtData[allsamples,]$species), labs, 50, TRUE, NULL, NULL, FALSE, TRUE)

    c(labs,allsamples)[negsig[,1]]
    c(labs,allsamples)[possig[,1]]
    list(positive=c(labs,allsamples)[possig[,1]],negative=c(labs,allsamples)[negsig[,1]])

    euktax[Descendants(euktreeY.root,137,'tips')[[1]],]

    i <- 1; hist((var_scales / prior_scales / dataset_scales[i])[(sum(M_all[1:(i-1)])+1):sum(M_all[1:i])])

}
