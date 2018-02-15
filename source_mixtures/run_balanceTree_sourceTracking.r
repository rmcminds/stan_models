library(coda)
library(rstan)
library(shinystan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

load('preparation.RData')

nchains=4



init <- lapply(1:nchains, function(...) list(env_prop_normal_raw=matrix(-5.0,nrow=(NEnvs),ncol=NEnvs)))

fit <- stan(file='balanceTree_sourceTracking.stan', data=list(NObs=NObs, NSamples=NSamples, sampledcounts=sampledcounts, datacounts=datacounts, samplenames=samplenames, nodenames=nodenames, NEstTaxNodes=NEstTaxNodes, NEnvs=NEnvs, envs=as.numeric(envs) ), control=list(adapt_delta=0.9, max_treedepth=20), iter=5000, chains=nchains, thin=2, init=init )

save.image(file='results.RData')

pdf(file='pairs.pdf')
pairs(fit,pars=c('env_props[1,1]','env_props[2,2]','env_props[3,3]','env_props[4,4]','lp__'))
graphics.off()



outdir <- 'sourcetracking_balances_out/'
dir.create(outdir)

taxdat <- read.table('rep_set_tax_assignments.txt',sep='\t',stringsAsFactors=F, row.names=1)
x <- strsplit(taxdat[,1],'; ')
most <- max(sapply(x,length))
parsedtax <- lapply(x,function(x) {length(x) <- most; return(x)})
tax <- do.call('rbind',parsedtax)
rownames(tax) <- rownames(taxdat)
colnames(tax) <- c('Kingdom','Phylum','Class','Order','Family','Genus','Species')

obsrats <- datacounts/sampledcounts

plotdat <- data.frame(obsrats=obsrats,samplenames=samplenames,nodenames=nodenames)
plotdat$envs <- sapply(plotdat$samplenames,function(x) return(envs[[x]]))

level <- colnames(tax)[1:7]

for(i in unique(nodenames)) {
    pdf(file=paste0(outdir,estnodes[[i]],'.pdf'),width=10)
    boxplot(obsrats~envs,data=plotdat[plotdat$nodenames==i,])
    graphics.off()
    
    
    for(node in Children(bacttreeY.root, estnodes[[i]])) {
        tips <- bacttreeY.root$tip.label[Descendants(bacttreeY.root, node,type='tips')[[1]]]
        cat(paste0('node ',node,':\n'),file=paste0(outdir,estnodes[[i]],'.txt'),append=T)
        if(length(tips) != 1) {
            cat(paste0(sort(unique(apply(tax[tips,level],1,paste,collapse='.'))),'\n'),file=paste0(outdir,estnodes[[i]],'.txt'),append=T)
            cat(paste0(sort(unique(tips)),'\n'),file=paste0(outdir,estnodes[[i]],'.txt'),append=T)
        } else {
            cat(paste0(paste(tax[tips,level],collapse='.'),'\n'),file=paste0(outdir,estnodes[[i]],'.txt'),append=T)
            cat(paste0(tips,'\n'),file=paste0(outdir,estnodes[[i]],'.txt'),append=T)
        }
    }
    
}
