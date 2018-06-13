library(ape)
library(phangorn)
library(geosphere)
library(rstan)
library(geiger)
library(phytools)
library(shinystan)
library(nlme)
library(reshape2)
library(paleotree)
library(parallel)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

outdir <- '/raid1/home/micro/mcmindsr/ryan/20180131_stansymbio/logistic_GLM_RMcMinds_nodeEffects_host16/'

load(paste0(outdir,'setup.RData'))
load(paste0(outdir,'fit.RData'))

outdir <- '/raid1/home/micro/mcmindsr/ryan/20180131_stansymbio/logistic_GLM_RMcMinds_nodeEffects_host16/reprocessed_3/'

for(i in 1:nparalleltrees) {
    
    check_hmc_diagnostics(fit[[i]])

    currdatadir <- paste0(outdir,'tree_',i,'/data/')
    currtabledir <- paste0(outdir,'tree_',i,'/tables/')

    dir.create(currdatadir, recursive=T)
    dir.create(paste0(currtabledir, 'nodes/'), recursive=T)
    dir.create(paste0(currtabledir, 'nodes_summed/'), recursive=T)

    scaledTaxNodeEffects <- array(c(scaledTaxNodeEffects, extract(fit[[i]], pars='scaledTaxNodeEffects', permuted=F, inc_warmup=T)), dim=c(NMCSamples, nchains, NEffects + NHostNodes, NEstTaxNodes, i), dimnames=list(sample=NULL, chain=NULL, effect=c(colnames(modelMat)[1:NEffects],colnames(hostAncestors[[i]])), taxnode=colnames(microbeAncestors), hostTree=i))
    
    
    summedScaledTaxNodeEffects <- array(NA, dim=c(NMCSamples, nchains, NEffects + NHostNodes, NEstTaxNodes, i), dimnames=list(sample=NULL, chain=NULL, effect=c(colnames(modelMat)[1:NEffects],colnames(hostAncestors[[i]])), taxnode=colnames(microbeAncestors), hostTree=i))
    baseLevelEffects <- array(NA, dim=c(NMCSamples, nchains, NFactors, NEstTaxNodes, i))
    summedBaseLevelEffects <- array(NA, dim=c(NMCSamples, nchains, NFactors, NEstTaxNodes, i))
    for(j in 1:NMCSamples) {
        for(k in 1:nchains) {
            summedScaledTaxNodeEffects[j,k,,,i] <- scaledTaxNodeEffects[j,k,,,i] %*% microbeAncestors
            for(m in sumconts) {
                baseLevelEffects[j,k,m,,i] <- -colSums(scaledTaxNodeEffects[j,k,factLevelMat[,m] == 1,,i])
                summedBaseLevelEffects[j,k,m,,i] <- -colSums(summedScaledTaxNodeEffects[j,k,factLevelMat[,m] == 1,,i])
            }
        }
    }
    
    for(l in 1:(NEffects + NHostNodes)) {
        yeah <- monitor(scaledTaxNodeEffects[,,l,i], warmup=warmup, probs=c(0.05,0.95))
        cat('\t', file = paste0(currtabledir, 'nodes/', dimnames(scaledTaxNodeEffects)[[3]][l], '.txt'))
        write.table(yeah, file = paste0(currtabledir, 'nodes/', dimnames(scaledTaxNodeEffects)[[3]][l], '.txt'), sep='\t', quote=F,append=T)
        
        yeah <- monitor(summedScaledTaxNodeEffects[,,l,i], warmup=warmup, probs=c(0.05,0.95))
        cat('\t', file = paste0(currtabledir, 'nodes_summed/', dimnames(summedScaledTaxNodeEffects)[[3]][l], '.txt'))
        write.table(yeah, file = paste0(currtabledir, 'nodes_summed/', dimnames(summedScaledTaxNodeEffects)[[3]][l], '.txt'), sep='\t', quote=F,append=T)
    }
    
    for(m in sumconts) {
        yeah <- monitor(baseLevelEffects[,,m,i], warmup=warmup, probs=c(0.05,0.95))
        cat('\t', file = paste0(currtabledir, 'nodes/', m, levels(newermap[,m])[nlevels(newermap[,m])], '.txt'))
        write.table(yeah, file = paste0(currtabledir, 'nodes/', m, levels(newermap[,m])[nlevels(newermap[,m])], '.txt'), sep='\t', quote=F,append=T)
        
        yeah <- monitor(summedBaseLevelEffects[,,m,i], warmup=warmup, probs=c(0.05,0.95))
        cat('\t', file = paste0(currtabledir, 'nodes_summed/', m, levels(newermap[,m])[nlevels(newermap[,m])], '.txt'))
        write.table(yeah, file = paste0(currtabledir, 'nodes_summed/', m, levels(newermap[,m])[nlevels(newermap[,m])], '.txt'), sep='\t', quote=F,append=T)
    }
    

}


currdatadir <- paste0(outdir,'alltrees/data/')
dir.create(currdatadir, recursive=T)

save(scaledTaxNodeEffects,file=paste0(currdatadir,'scaledTaxNodeEffects.RData'))
save(summedScaledTaxNodeEffects,file=paste0(currdatadir,'summedScaledTaxNodeEffects.RData'))
save(baseLevelEffects,file=paste0(currdatadir,'baseLevelEffects.RData'))
save(summedBaseLevelEffects,file=paste0(currdatadir,'summedBaseLevelEffects.RData'))

