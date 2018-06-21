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

outdir <- file.path('output','2018-06-05_17-12-31')

##
load(file=file.path(outdir,'setup.RData'))
load(file=file.path(outdir,'fit.RData'))
##

## summarize results for parameters that can be interpretted across all sampled host trees
currdiagnosticdir <- file.path(outdir,'diagnostics')
dir.create(currdiagnosticdir, recursive=T)

allfit <- sflist2stanfit(fit)

pdf(file=file.path(currdiagnosticdir,'pairs_StDProps.pdf'), width=50, height=50)
pairs(allfit, pars=c('stDProps', 'timeBinProps', 'globalIntercept', 'aveStD'))
graphics.off()

pdf(file=file.path(currdiagnosticdir,'pairs_StDPropsLog.pdf'), width=50, height=50)
pairs(allfit, pars=c('stDProps', 'timeBinProps', 'globalIntercept', 'aveStD'), log=T)
graphics.off()

pars <- paste0('phyloLogVarMultPrev[',1:20,']')

pdf(file=file.path(currdiagnosticdir,'pairs_phyloLogVarMultPrev.pdf'), width=50, height=50)
pairs(allfit, pars=pars, log=T)
graphics.off()

pars <- c(paste0('phyloLogVarMultADiv[',1:10,']'), paste0('rawAlphaDivEffects[',1:10,']'))

pdf(file=file.path(currdiagnosticdir,'pairs_phyloLogVarMultADiv.pdf'), width=100, height=100)
pairs(allfit, pars=pars, log=T)
graphics.off()


## fin
