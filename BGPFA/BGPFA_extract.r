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

mybiplot <- function(object1, object2, a1 = 1, a2 = 2, fact, labs = rep('',nrow(object2)), nv = nrow(object2), normalize = FALSE, pca = FALSE, center=FALSE, effsBoth = FALSE) {

    if(pca) {
        pcob <- prcomp(object1, center=center)
        ob1 <- pcob$x[,c(a1,a2)]
        ob2 <- (object2 %*% pcob$rotation)[,c(a1,a2)]
    } else {
        ob1 <- object1[,ncol(object1) - c(a1,a2) + 1]
        ob2 <- object2[,ncol(object2) - c(a1,a2) + 1]
    }

    if(normalize) {
        ob1[,1] <- ob1[,1] / max(abs(ob1[,1]))
        ob1[,2] <- ob1[,2] / max(abs(ob1[,2]))
        ob2[,1] <- ob2[,1] / max(abs(ob2[,1]))
        ob2[,2] <- ob2[,2] / max(abs(ob2[,2]))
    }

    plot(ob1,
         col  = fact,
         pch  = 16,
         xlim = range(c(ob1[,1],ob2[,1])),
         ylim = range(c(ob1[,2],ob2[,2])),
         cex  = 0.5)
    legend("bottomright", levels(fact), col=1:nlevels(fact), pch=16)

    effs <- setNames(apply(ob2, 1, function(x) sqrt(sum(x^2))), labs)

    if(effsBoth) {
        keep <- tail(order(effs),nv)
        ob2 <- ob2[keep, ]
        res <- effs[rev(labs[keep])]
    } else {
        keep1 <- tail(order(abs(ob2[,1])),nv)
        keep2 <- tail(order(abs(ob2[,2])),nv)
        keep <- union(keep1, keep2)
        res <- list(rev(setNames(ob2[keep1,1], labs[keep1])), rev(setNames(ob2[keep2,2], labs[keep2])))
        ob2 <- ob2[keep, ]
    }

    text(ob2, labels=labs[keep], cex=0.4)
    return(res)
}

myplotres <- function(object1, a1, a2, labs, nv) {

    ob1 <- object1[,c(a1,a2)]
    effs <- apply(ob1, 1, function(x) sqrt(sum(x^2)))
    keep <- tail(order(effs), nv)
    ob1 <- ob1[keep, ]
    newlabs = labs[keep]

    plot(ob1,
         pch = 16,
         cex = 0.1)
    text(ob1,
         labels = newlabs,
         cex = 0.4)
}

sumfunc <- median
#sumfunc <- mean
#sumfunc <- function(x) x[length(x)]
#sumfunc <- function(x) x[maxind]

labs <- c(varlabs, paste0(varlabs[1:sum(Mplus[1:D])],'.binary'), paste0('countDatasetBinaryInt_',1:D))

lps <- extract(stan.fit.vb, pars='log_post', permuted=FALSE)[,1,1]
maxind <- min(which(lps == max(lps)))

dataset_scales.vb <- extract(stan.fit.vb, pars='dataset_scales', permuted=FALSE)
dataset_scales.vb <- apply(dataset_scales.vb,3,sumfunc)

var_scales <- extract(stan.fit.vb, pars='var_scales', permuted=FALSE)
var_scales.vb <- apply(var_scales,3,sumfunc)
names(var_scales.vb) <- labs

global_effect_scale.vb <- extract(stan.fit.vb, pars='global_effect_scale', permuted=FALSE)
global_effect_scale.vb <- apply(global_effect_scale.vb, 3, sumfunc)

latent_scales.vb <- extract(stan.fit.vb, pars='latent_scales', permuted=FALSE)
latent_scales.vb <- apply(latent_scales.vb, 3, sumfunc)
axisOrder <- order(latent_scales.vb, decreasing=TRUE)
latent_scales.vb <- latent_scales.vb[axisOrder]

weight_scales.vb <- extract(stan.fit.vb, pars='weight_scales', permuted=FALSE)
weight_scales.vb <- apply(weight_scales.vb,3,sumfunc)
weight_scales.vb <- array(weight_scales.vb,dim=c(length(weight_scales.vb)/K,K))
weight_scales.vb <- array(weight_scales.vb[,axisOrder], dim=c(length(weight_scales.vb)/K,K))

if('sds' %in% importparams) {
    sds.vb <- extract(stan.fit.vb, pars='sds', permuted=FALSE)
    sds.vb <- apply(sds.vb,3,sumfunc)
}

nu_factors.vb <- extract(stan.fit.vb, pars='nu_factors', permuted=FALSE)
nu_factors.vb <- matrix(apply(nu_factors.vb, 3, sumfunc), ncol=K)
nu_factors.vb <- nu_factors.vb[,axisOrder]

if('log_less_contamination' %in% importparams) {
    log_less_contamination.vb <- extract(stan.fit.vb, pars='log_less_contamination', permuted=FALSE)
    log_less_contamination.vb <- apply(log_less_contamination.vb, 3, sumfunc)
}

if('contaminant_overDisp' %in% importparams) {
    contaminant_overDisp.vb <- extract(stan.fit.vb, pars='contaminant_overDisp', permuted=FALSE)
    contaminant_overDisp.vb <- apply(contaminant_overDisp.vb, 3, sumfunc)
}

Z_all <- extract(stan.fit.vb, pars='Z', permuted=FALSE)

Z.vb <- apply(Z_all,3,sumfunc)
Z.vb <- t(array(Z.vb,dim=c(K,length(Z.vb)/K)))
Z.vb <- array(Z.vb[,axisOrder],dim=c(length(Z.vb)/K,K))
rownames(Z.vb) <- allsamples

Z_all <- array(Z_all, dim=c(dim(Z_all)[[1]],N,K))[,,axisOrder,drop=FALSE]

W_norm <- extract(stan.fit.vb, pars='W_norm', permuted=FALSE)
W_norm.vb <- apply(W_norm,3,sumfunc)
W_norm.vb <- array(W_norm.vb,dim=c(length(W_norm.vb)/K,K))[,axisOrder]

binary_count_dataset_intercepts.vb <- extract(stan.fit.vb, pars='binary_count_dataset_intercepts', permuted=FALSE)
binary_count_dataset_intercepts.vb <- apply(binary_count_dataset_intercepts.vb, 3, sumfunc)

DRC <- D+R+C

rhoSites.vb <- extract(stan.fit.vb, pars='rhoSites', permuted=FALSE)
rhoSites.vb <- apply(rhoSites.vb, 3, sumfunc)

covSites.vb <- extract(stan.fit.vb, pars='covSites', permuted=FALSE)
covSites.vb <- apply(covSites.vb, 3, sumfunc)
dim(covSites.vb) <- c(K,nSites,nSites)
covSites.vb <- covSites.vb[axisOrder,,]

if('rhoZ' %in% importparams) {
    rhoZ.vb <- extract(stan.fit.vb, pars='rhoZ', permuted=FALSE)
    if(exists('KG')) {
        rhoZ.vb <- array(apply(rhoZ.vb, 3, sumfunc),dim=c(K_linear,KG))[axisOrder[axisOrder <= K_linear],]
    } else {
        rhoZ.vb <- apply(rhoZ.vb, 3, sumfunc)[axisOrder[axisOrder <= K_linear]]
    }
}

pos1 <- array(apply(W_norm,3,monteCarloP,pn='p'),dim=c(dim(W_norm)[[3]]/K,K))[,axisOrder]
possig <- pos1 < 0.05
neg1 <- array(apply(W_norm,3,monteCarloP,pn='n'),dim=c(dim(W_norm)[[3]]/K,K))[,axisOrder]
negsig <- neg1 < 0.05
anysig <- possig | negsig
rownames(pos1) <- rownames(neg1) <- labs
gc()

nullfunc <- function() {

    gc()

    drivers <- mytriplot(Z.vb, W_norm.vb, Z.vb, 1,2, as.factor(filtData[allsamples,]$species), labs, 50, TRUE, NULL, NULL, FALSE, TRUE)

    drivers <- mytriplot(Z.vb, W_norm.vb, Z.vb, 1,2, as.factor(filtData[allsamples,]$species), labs, 50, TRUE, NULL, NULL, FALSE, TRUE, anysig)

    i <- 1
    hist((var_scales.vb / priorScales/dataset_scales.vb[i])[(sum(Mplus[1:(i-1)])+1):sum(Mplus[1:i])])

    labs[negsig[,1]]
    labs[possig[,1]]
    list(positive=labs[possig[,1]],negative=labs[negsig[,1]])

}
