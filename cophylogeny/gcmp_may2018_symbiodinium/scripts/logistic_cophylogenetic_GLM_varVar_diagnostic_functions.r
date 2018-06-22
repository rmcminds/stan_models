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

outdir <- file.path('output','2018-06-19_10-31-37')

##
load(file=file.path(outdir,'setup.RData'))
load(file=file.path(outdir,'fit.RData'))
##

allfit <- c(fit, sflist2stanfit(fit))

for(i in 1:(NTrees + 1)) {
    currdiagnosticdir <- file.path(outdir, 'diagnostics', if(i <= NTrees) {paste0('tree_', i)} else {'allfit'})
    dir.create(currdiagnosticdir, recursive=T)
    
    pars <- c('stDProps', 'hostMetaVarProps', 'timeBinMetaVar', 'metaVarProps', 'aveStD')
    
    pdf(file=file.path(currdiagnosticdir,'pairs_grabBag.pdf'), width=50, height=50)
    pairs(allfit[[i]], pars=pars)
    graphics.off()
    
    pdf(file=file.path(currdiagnosticdir,'pairs_grabBagLog.pdf'), width=50, height=50)
    pairs(allfit[[i]], pars=pars, log=T)
    graphics.off()
    
    pars <- paste0('phyloLogVarMultPrev[',1:20,']')
    
    pdf(file=file.path(currdiagnosticdir,'pairs_phyloLogVarMultPrev.pdf'), width=50, height=50)
    pairs(allfit[[i]], pars=pars, log=T)
    graphics.off()
    
    pars <- c(paste0('phyloLogVarMultADiv[',1:10,']'), paste0('rawAlphaDivEffects[',1:10,']'))
    
    pdf(file=file.path(currdiagnosticdir,'pairs_phyloLogVarMultADiv.pdf'), width=100, height=100)
    pairs(allfit[[i]], pars=pars, log=T)
    graphics.off()
}

## fin
