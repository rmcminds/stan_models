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

orient_axes <- function(extracted, reference = NULL) {
    K <- dim(extracted)[[2]]
    if(is.null(reference)) reference <- extracted[,,sample(dim(extracted)[[3]],1)]
    oriented <- array(NA, dim=dim(extracted))
    which_match <- array(NA,dim(extracted)[c(2,3)])
    reflected <- array(NA,dim(extracted)[c(2,3)])
    for(draw in 1:dim(extracted)[[3]]) {
        available <- rep(TRUE,K)
        for(axis in 1:K) {
            idx <- (1:K)[available]
            cors <- sapply(idx, function(k) cor(extracted[,k,draw], reference[,axis]))
            which_match[axis,draw] <- idx[which.max(cors)]
            available[which_match[axis,draw]] <- FALSE
            reflected[axis,draw] <- sign(cors[which.max(cors)])
            oriented[,axis,draw] <- reflected[axis,draw] * extracted[,which_match[axis,draw],draw]
        }
    }
    return(list(oriented=oriented, which_match=which_match, reflected=reflected))
} # swap axes such that they best match a reference (by default a random draw). If axes are degenerate during model fitting, this should orient them. May have issues if there are multiple very highly correlated axes?

sumfunc <- median
#sumfunc <- mean
#sumfunc <- function(x) x[length(x)]
#sumfunc <- function(x) x[maxind]

labs <- c(varlabs, paste0(varlabs[1:sum(M_all[1:D])],'.binary'), paste0('binary_count_dataset_int_',dataset_names[1:D]))
dataset_names_expanded <- c(dataset_names, paste0(dataset_names[1:D],'.binary'))

dataset_scales <- extract(stan.fit, pars='dataset_scales', permuted=FALSE)
dataset_scales <- apply(dataset_scales,3,sumfunc)
names(dataset_scales) <- dataset_names_expanded

var_scales <- extract(stan.fit, pars='var_scales', permuted=FALSE)
var_scales <- apply(var_scales,3,sumfunc)
names(var_scales) <- labs

global_effect_scale <- extract(stan.fit, pars='global_effect_scale', permuted=FALSE)
global_effect_scale <- apply(global_effect_scale, 3, sumfunc)

W_norm_all <- extract(stan.fit, pars='W_norm', permuted=FALSE)
W_norm_all <- array(W_norm_all, dim = c(dim(W_norm_all)[[1]], dim(W_norm_all)[[3]]/K, K))
Z_all <- extract(stan.fit, pars='Z', permuted=FALSE)
Z_all <- aperm(array(Z_all, dim = c(dim(Z_all)[[1]], K, dim(Z_all)[[3]]/K)), c(1,3,2))
WZ_all_raw <- aperm(abind:::abind(W_norm_all,Z_all,along=2), c(2,3,1))
mutated <- Morpho:::procSym(WZ_all_raw, CSinit = FALSE, scale = FALSE, reflect = TRUE, orp = FALSE, bending = FALSE, pcAlign = FALSE)
mutated_projected <- shapes:::procOPA(WZ_all_raw[,,1], mutated$mshape, scale=FALSE)
oriented <- orient_axes(WZ_all_raw, mutated_projected$Bhat)
WZ_all <- aperm(oriented$oriented, c(3,1,2))
axisOrder_all <- oriented$which_match

latent_scales <- extract(stan.fit, pars='latent_scales', permuted=FALSE)
latent_scales <- sapply(1:dim(latent_scales)[[1]], function(x) latent_scales[x,,axisOrder_all[,x]])
latent_scales <- apply(latent_scales, 1, sumfunc)
axisOrder <- order(latent_scales, decreasing=TRUE)
latent_scales <- latent_scales[axisOrder]

WZ <- apply(WZ_all,c(2,3),sumfunc)
WZ <- array(WZ, dim = c(length(WZ)/K, K))[,axisOrder]
W_norm <- WZ[1:dim(W_norm_all)[2],]
Z <- WZ[(dim(W_norm_all)[2]+1):nrow(WZ),]
rownames(Z) <- allsamples

weight_scales <- extract(stan.fit, pars='weight_scales', permuted=FALSE)
##need to reorder each draw according to axisOrder_all- this is currently inaccurate
weight_scales <- apply(weight_scales,3,sumfunc)
weight_scales <- array(weight_scales,dim=c(length(weight_scales)/K,K))
weight_scales <- array(weight_scales[,axisOrder], dim=c(length(weight_scales)/K,K))
rownames(weight_scales) <- dataset_names_expanded

if('sds' %in% importparams) {
    sds <- extract(stan.fit, pars='sds', permuted=FALSE)
    sds <- apply(sds,3,sumfunc)
}

nu_factors <- extract(stan.fit, pars='nu_factors', permuted=FALSE)
##need to reorder each draw according to axisOrder_all- this is currently inaccurate
nu_factors <- matrix(apply(nu_factors, 3, sumfunc), ncol=K)
nu_factors <- nu_factors[,axisOrder]
rownames(nu_factors) <- dataset_names_expanded

if('log_less_contamination' %in% importparams) {
    log_less_contamination <- extract(stan.fit, pars='log_less_contamination', permuted=FALSE)
    log_less_contamination <- apply(log_less_contamination, 3, sumfunc)
}

if('contaminant_overdisp' %in% importparams) {
    contaminant_overdisp <- extract(stan.fit, pars='contaminant_overdisp', permuted=FALSE)
    contaminant_overdisp <- apply(contaminant_overdisp, 3, sumfunc)
}

binary_count_dataset_intercepts <- extract(stan.fit, pars='binary_count_dataset_intercepts', permuted=FALSE)
binary_count_dataset_intercepts <- apply(binary_count_dataset_intercepts, 3, sumfunc)
names(binary_count_dataset_intercepts) <- dataset_names[1:D]

DRC <- D+R+C

rho_sites <- extract(stan.fit, pars='rho_sites', permuted=FALSE)
rho_sites <- apply(rho_sites, 3, sumfunc)

if('corr_sites' %in% importparams) {
    corr_sites <- extract(stan.fit, pars='corr_sites', permuted=FALSE)
    corr_sites <- apply(corr_sites, 3, sumfunc)
    dim(corr_sites) <- c(K,N_sites,N_sites)
    corr_sites <- corr_sites[axisOrder,,]
}

if('rho_Z' %in% importparams & KG > 0) {
    rho_Z <- extract(stan.fit, pars='rho_Z', permuted=FALSE)
    ##need to reorder each draw according to axisOrder_all- this is currently inaccurate
    rho_Z <- array(apply(rho_Z, 3, sumfunc),dim=c(K_linear,KG))[axisOrder[axisOrder <= K_linear],]
}

if('skew_Z' %in% importparams) {
    skew_Z <- extract(stan.fit, pars='skew_Z', permuted=FALSE)
    ##need to reorder each draw according to axisOrder_all- this is currently inaccurate
    skew_Z <- apply(skew_Z, 3, sumfunc)[axisOrder]
}

pos1 <- apply(WZ_all,c(2,3),monteCarloP,pn='p')[,axisOrder]
possig <- pos1 < 0.05
neg1 <- apply(WZ_all,c(2,3),monteCarloP,pn='n')[,axisOrder]
negsig <- neg1 < 0.05
anysig <- possig | negsig
rownames(pos1) <- rownames(neg1) <- c(labs,allsamples)
gc()

nullfunc <- function() {

    gc()

    drivers <- mytriplot(Z, W_norm, Z, 1,2, as.factor(filtData[allsamples,]$species), labs, 50, TRUE, NULL, NULL, FALSE, TRUE)

    drivers <- mytriplot(Z, W_norm, Z, 1,2, as.factor(filtData[allsamples,]$species), labs, 50, TRUE, NULL, NULL, FALSE, TRUE, anysig)

    i <- 1
    hist((var_scales / prior_scales / dataset_scales[i])[(sum(M_all[1:(i-1)])+1):sum(M_all[1:i])])

    c(labs,allsamples)[negsig[,1]]
    c(labs,allsamples)[possig[,1]]
    list(positive=c(labs,allsamples)[possig[,1]],negative=c(labs,allsamples)[negsig[,1]])

}
