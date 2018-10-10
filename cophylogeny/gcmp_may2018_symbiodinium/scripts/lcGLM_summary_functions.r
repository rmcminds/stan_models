makeParentMat <- function(NNodes, tree, NTips, tipNames) {
    parentMat <- matrix(0, NNodes + 1, NNodes + 1)
    for(node in 1:(NNodes + 1)) {
        parentMat[node, ] <- as.numeric(1:(NNodes + 1) %in% c(Ancestors(tree, node, 'parent'), node))
    }
    colnames(parentMat) <- rownames(parentMat) <- paste0('i', 1:(NNodes + 1))
    colnames(parentMat)[1:NTips] <- rownames(parentMat)[1:NTips] <- paste0('t', tipNames)
    parentMat <- cbind(1, parentMat[-(NTips + 1), -(NTips + 1)])
    return(parentMat)
}

summarizeLcGLM <- function(...) {
    
    fitModes <- sapply(fit, function(x) x@mode)
    
    ## create parental ancestry matrix for microbes (only sums effect of a given node and its direct parent, not all ancestors
    microbeParentsT <- t(makeParentMat(NMicrobeNodes,
                                       finalMicrobeTree,
                                       NMicrobeTips,
                                       microbeTips))

    hostParents <- list()
    microbeAncestorsT <- t(cbind(1, microbeAncestors))

    colorpal <- colorRampPalette(brewer.pal(9, 'Blues'))
    plotcolors <- c('white', colorpal(100), 'black')
    
    baseNames <- sapply(sumconts, function(m) paste0(m, levels(newermap[,m])[nlevels(newermap[,m])]))

    ## summarize the results separately for each sampled host tree
    for(i in 1:NTrees) {
        
        if (fitModes[[i]] != 0) next
        
        check_hmc_diagnostics(fit[[i]])

        currplotdir <- file.path(outdir, paste0('tree_', i), 'plots')
        currtabledir <- file.path(outdir, paste0('tree_', i), 'tables')
        currdatadir <- file.path(outdir, paste0('tree_', i), 'data')

        dir.create(currplotdir, recursive = T)
        dir.create(currtabledir, recursive = T)
        dir.create(currdatadir, recursive = T)

        ## plot the sampled tree with the time bins marked
        pdf(file = file.path(currplotdir, 'sampledHostTree.pdf'), width = 25, height = 15)
        plot(hostTreesSampled[[i]], cex = 0.75)
        for(age in hostTreeDetails[[i]]$maxNHs - meanHostBoundaries) {
            lines(x = c(age, age), y = c(1, NHostTips), lwd = 1)
        }
        graphics.off()
        ##
        
        ## plot the sampled tree with the time bins marked
        pdf(file = file.path(currplotdir, 'sampledMicrobeTree.pdf'), width = 25, height = 15)
        plot(finalMicrobeTree, cex = 0.75)
        for(age in microbeTreeDetails$maxNHs - microbeCutPoints$boundaries) {
            lines(x = c(age, age), y = c(1, NMicrobeTips), lwd = 1)
        }
        graphics.off()
        ##
        
        ## summarize the mean branch lengths of the microbes
        sums <- summary(fit[[i]], pars = 'microbeScales', probs = c(0.05,0.95), use_cache = F)
        newEdges <- sums$summary[,'mean']^2
        finalMicrobeTree.newEdges <- finalMicrobeTree
        finalMicrobeTree.newEdges$edge.length <- newEdges[order(microbeTreeDetails$edgeOrder)]
        pdf(file = file.path(currplotdir, 'microbeTreeWEstimatedEdgeLengths.pdf'), width = 25, height = 15)
        plot(finalMicrobeTree.newEdges, cex = 0.5)
        graphics.off()
        ##
        
        ## summarize the mean branch lengths of the hosts
        sums <- summary(fit[[i]], pars = 'hostScales', probs = c(0.05,0.95), use_cache = F)
        newEdges <- sums$summary[,'mean']^2
        hostTreesSampled.newEdges <- hostTreesSampled[[i]]
        hostTreesSampled.newEdges$edge.length <- newEdges[order(hostTreeDetails[[i]]$edgeOrder)]
        pdf(file = file.path(currplotdir, 'hostTreeWEstimatedEdgeLengths.pdf'), width = 25, height = 15)
        plot(hostTreesSampled.newEdges, cex = 0.5)
        graphics.off()
        ##
        
        ## plot heatmap of cophylogenetic patterns
        plotmicrobetree <- ladderize(multi2di(finalMicrobeTree))
        hclmicrobetree <- as.hclust(plotmicrobetree)
        # for each sample in the data, assign its mitotype to a vector
        hostvect <- sampleMap[[i]]$host_scientific_name
        # name the vector with the sample IDs
        names(hostvect) <- rownames(sampleMap[[i]])
        # sort the samples in the vector
        temp <- sampleMap[[i]][order(sampleMap[[i]]$concatenated_date),]
        temp <- temp[order(temp$reef_name),]
        temp <- temp[order(temp$host_scientific_name),]
        for(comp in c('T', 'S', 'M')) {
            temp2 <- temp[temp$tissue_compartment == comp,]
            hostvectTemp <- hostvect[rownames(temp2)]
            # expand the tips into polytomies containing a tip for each sample
            hosttree <- expandTaxonTree(hostTreesSampled[[i]], hostvectTemp, keepBrLen = T)
            hosttree <- drop.tip(hosttree,
                                 hosttree$tip.label[!hosttree$tip.label %in% names(hostvectTemp)])
            # convert polytomies into randomly-split, binary subtrees with 0-length branches, and ladderize the whole tree
            hosttree.dichotomous <- ladderize(multi2di(hosttree, random = F), right = F)
            # convert the phylogenetic tree into an hclust object
            hclhosttree <- as.hclust(force.ultrametric(hosttree.dichotomous))
            plotFilt <- as.matrix(t(y)[plotmicrobetree$tip.label, hosttree.dichotomous$tip.label])
            pdf(file   = file.path(currplotdir,
                                   paste0('cophylogeny_heatmap_',
                                          comp,
                                          '.pdf')),
                width  = 10,
                height = 10)
            heatmap(plotFilt,
                    Rowv   = as.dendrogram(hclmicrobetree),
                    Colv   = as.dendrogram(hclhosttree),
                    col    = plotcolors,
                    cexCol = 0.2,
                    cexRow = 0.1,
                    scale  = 'none')
            graphics.off()
        }
        ##

        ## variance partitioning
        stDProps <- extract(fit[[i]], pars = 'stDProps')[[1]]
        colnames(stDProps) <- c(paste0('ADiv.', colnames(factLevelMat)),
                                paste0('Specificity.', colnames(factLevelMat)),
                                'ADiv.host',
                                'host.specificity',
                                'microbe.prevalence')

        pdf(file = file.path(currplotdir,'scalesboxes.pdf'), width = 25, height = 15)
        boxplot(stDProps, cex.axis = 0.5, las = 2)
        graphics.off()

        save(stDProps, file = file.path(currdatadir, 'stDProps.RData'))
        ##
        
        ## proportion of variance explained by each time bin
        metaVarProps <- extract(fit[[i]], pars = 'metaVarProps')[[1]]
        colnames(metaVarProps) <- c('prevalence', 'adiv', 'specificty')
        pdf(file = file.path(currplotdir, 'metaVarProps_boxes.pdf'), width = 25, height = 15)
        boxplot(metaVarProps, cex.axis = 0.5, las = 2)
        graphics.off()
        save(metaVarProps, file = file.path(currdatadir, 'metaVarProps.RData'))
        
        metaScales <- array(extract(fit[[i]],
                            pars = 'metaScales',
                            permuted = F,
                            inc_warmup = T),
                    dim = c(NMCSamples,
                            NChains,
                            3))
                            
        pdf(file = file.path(currplotdir, 'metaScales_boxes.pdf'), width = 25, height = 15)
        boxplot(t(apply(metaScales, c(1,3), rbind)), cex.axis = 0.5, las = 2)
        graphics.off()
                            
        save(metaScales, file = file.path(currdatadir, 'metaScales.RData'))
        
        
        ## compare rates of host evolution in each time bin
        relativeHostEvolRates <- exp(extract(fit[[i]], pars = 'logRelativeHostEvolRates')[[1]])
        colnames(relativeHostEvolRates) <- sapply(1:(length(meanHostBoundariesRounded) + 1),
                                                  function(x) paste0(c('before', meanHostBoundariesRounded)[[x]],
                                                                     ' - ',
                                                                     c(meanHostBoundariesRounded, 'present')[[x]],
                                                                     ' mya'))

        pdf(file = file.path(currplotdir, 'evolRatesHost.pdf'), width = 25, height = 15)
        boxplot(relativeHostEvolRates, xlab = 'Time Period', ylab = 'Rate of Evolution Relative to Mean')
        graphics.off()
        
        pdf(file = file.path(currplotdir, 'logEvolRatesHost.pdf'), width = 25, height = 15)
        boxplot(log(relativeHostEvolRates), xlab = 'Time Period', ylab = 'Log Rate of Evolution Relative to Mean')
        graphics.off()

        save(relativeHostEvolRates, file = file.path(currdatadir, 'relativeHostEvolRates.RData'))
        ##
        
        ## compare rates of microbe evolution in each time bin
        relativeMicrobeEvolRates <- exp(extract(fit[[i]], pars = 'logRelativeMicrobeEvolRates')[[1]])
        
        pdf(file = file.path(currplotdir, 'evolRatesMicrobe.pdf'), width = 25, height = 15)
        boxplot(relativeMicrobeEvolRates, xlab = 'Time Period', ylab = 'Rate of Evolution Relative to Mean')
        graphics.off()
        
        pdf(file = file.path(currplotdir, 'logEvolRatesMicrobe.pdf'), width = 25, height = 15)
        boxplot(log(relativeMicrobeEvolRates), xlab = 'Time Period', ylab = 'Log Rate of Evolution Relative to Mean')
        graphics.off()
        
        save(relativeMicrobeEvolRates, file = file.path(currdatadir, 'relativeMicrobeEvolRates.RData'))
        ##
        
        ## ornstein-uhlenbeck parameters
        hostOUAlpha <- extract(fit[[i]], pars = 'hostOUAlpha')[[1]]
        microbeOUAlpha <- extract(fit[[i]], pars = 'microbeOUAlpha')[[1]]
        OUAlphas <- cbind(hostOUAlpha, microbeOUAlpha)
        colnames(OUAlphas) <- c('host', 'microbe')
        
        pdf(file = file.path(currplotdir, 'OUAlphas.pdf'), width = 25, height = 15)
        boxplot(OUAlphas, xlab = 'Host or microbe', ylab = 'alpha')
        graphics.off()
        
        save(OUAlphas, file = file.path(currdatadir, 'OUAlphas.RData'))
        ##
        
        ## extract effects from model
        scaledMicrobeNodeEffects <- array(extract(fit[[i]],
                                                  pars       = 'scaledMicrobeNodeEffects',
                                                  permuted   = F,
                                                  inc_warmup = T),
                                          dim = c(NMCSamples,
                                                  NChains,
                                                  NEffects + NHostNodes + 1,
                                                  NMicrobeNodes + 1),
                                          dimnames = list(sample  = NULL,
                                                          chain   = NULL,
                                                          effect  = c('microbePrevalence',
                                                                      colnames(modelMat)[1:NEffects],
                                                                      paste0('host_', colnames(hostAncestors[[i]]))),
                                                          taxnode = c('alphaDiversity', colnames(microbeAncestors))))
                                                        
        save(scaledMicrobeNodeEffects, file = file.path(currdatadir, 'scaledMicrobeNodeEffects.RData'))
        ##
        
        ## calculate base-level effects (negative sum of all others in category)
        baseLevelEffects <- array(NA,
                                  dim = c(NMCSamples,
                                          NChains,
                                          length(sumconts),
                                          NMicrobeNodes + 1),
                                  dimnames = list(sample  = NULL,
                                                  chain   = NULL,
                                                  effect  = sumconts,
                                                  taxnode = c('alphaDiversity', colnames(microbeAncestors))))
        for(j in 1:NMCSamples) {
            for(k in 1:NChains) {
                for(m in sumconts) {
                    baseLevelEffects[j,k,m,] <- -colSums(scaledMicrobeNodeEffects[j, k, rownames(factLevelMat)[factLevelMat[,m] == 1],])
                }
            }
        }
        
        save(baseLevelEffects, file = file.path(currdatadir, 'baseLevelEffects.RData'))
        ##

        ## summarize posterior distributions of effects
        allRes <- NULL
        for(j in 1:(NEffects + NHostNodes + 1)) {
            temp <- monitor(array(scaledMicrobeNodeEffects[,,j,],
                                  dim = c(NMCSamples, NChains, NMicrobeNodes + 1)),
                            warmup = warmup,
                            probs  = c(0.05, 0.95))
            temp <- cbind(effect = dimnames(scaledMicrobeNodeEffects)[[3]][j], temp)
            rownames(temp) <- c('alphaDiversity', rownames(microbeAncestors))
            allRes <- rbind(allRes, temp)
        }
        ##
        
        ## summarize posterior distibutions of base-level effects
        for(m in sumconts) {
            temp <- monitor(array(baseLevelEffects[,,m,],
                                  dim = c(NMCSamples, NChains, NMicrobeNodes + 1)),
                            warmup = warmup,
                            probs = c(0.05, 0.95))
            temp <- cbind(effect = baseNames[[m]], temp)
            rownames(temp) <- c('alphaDiversity', rownames(microbeAncestors))
            allRes <- rbind(allRes, temp)
        }
        ##
        
        ## write posterior summaries of effects to file
        cat('microbeNode\t', file = file.path(currtabledir, 'nodeEffects.txt'))
        write.table(allRes,
                    file   = file.path(currtabledir, 'nodeEffects.txt'),
                    sep    = '\t',
                    quote  = F,
                    append = T)
        ##
        
        ## extract variance modifier terms from the fit model
        phyloLogVarMultPrev <- array(extract(fit[[i]],
                                             pars = 'phyloLogVarMultPrev',
                                             permuted = F,
                                             inc_warmup = T),
                                     dim = c(NMCSamples,
                                             NChains,
                                             NMicrobeNodes),
                                     dimnames = list(sample = NULL,
                                                     chain = NULL,
                                                     taxnode = colnames(microbeAncestors)))
        save(phyloLogVarMultPrev, file = file.path(currdatadir, 'phyloLogVarMultPrev.RData'))
        
        phyloLogVarMultADiv <- array(extract(fit[[i]],
                                             pars = 'phyloLogVarMultADiv',
                                             permuted = F,
                                             inc_warmup = T),
                                     dim = c(NMCSamples,
                                             NChains,
                                             NHostNodes),
                                     dimnames = list(sample = NULL,
                                                     chain = NULL,
                                                     taxnode = colnames(hostAncestors[[i]])))
        save(phyloLogVarMultADiv, file = file.path(currdatadir, 'phyloLogVarMultADiv.RData'))
        
        phyloLogVarMultRaw <- array(extract(fit[[i]],
                                            pars = 'phyloLogVarMultRaw',
                                            permuted = F,
                                            inc_warmup = T),
                                    dim = c(NMCSamples,
                                            NChains,
                                            NHostNodes,
                                            NMicrobeNodes),
                                    dimnames = list(sample = NULL,
                                                    chain = NULL,
                                                    hostnode = colnames(hostAncestors[[i]]),
                                                    microbenode = colnames(microbeAncestors)))
        save(phyloLogVarMultRaw, file = file.path(currdatadir, 'phyloLogVarMultRaw.RData'))
        
        matMult <- array(NA,
                         dim = c(NMCSamples,
                                 NChains,
                                 NHostNodes,
                                 NMicrobeNodes),
                         dimnames = list(sample = NULL,
                                         chain = NULL,
                                         hostnode = colnames(hostAncestors[[i]]),
                                         microbenode = colnames(microbeAncestors)))
        ##
        
        ## see if any clades have higher variance among their descendants (maybe suggesting codiversification)
        allRes <- NULL
        for(j in 1:NHostNodes) {
            temp <- monitor(array(phyloLogVarMultRaw[,,j,],
                                  dim = c(NMCSamples, NChains, NMicrobeNodes)),
                            warmup = warmup,
                            probs = c(0.05, 0.95))
            temp <- cbind(hostNode = colnames(hostAncestors[[i]])[[j]], temp)
            rownames(temp) <- rownames(microbeAncestors)
            allRes <- rbind(allRes, temp)
        }
        
        cat('microbeNode\t', file = file.path(currtabledir, 'individualPhyloVarianceEffects.txt'))
        write.table(allRes,
                    file   = file.path(currtabledir, 'individualPhyloVarianceEffects.txt'),
                    sep    = '\t',
                    quote  = F,
                    append = T)
        ##
        
        ## since the model has many degrees of freedom, it's possible that certain regions of the tree have 'significant' effects but that those effects are 'implemented' via different, nearby nodes in each iteration. We can sum the effects of all ancestral nodes to see if a given node is consistently influenced by all of itself plus its ancestral terms, even if all of those terms are individually inconsistent
        for(j in 1:NMCSamples) {
            for(k in 1:NChains) {
                matMult[j,k,,] <- cbind(1,
                                        hostAncestors[[i]]) %*%
                                  rbind(c(0,
                                          phyloLogVarMultPrev[j,k,] * metaScales[j,k,1]),
                                        cbind(phyloLogVarMultADiv[j,k,] * metaScales[j,k,2],
                                              phyloLogVarMultRaw[j,k,,] * metaScales[j,k,3])) %*%
                                  microbeAncestorsT
            }
        }
        
        allRes <- NULL
        for(j in 1:NHostNodes) {
            temp <- monitor(array(matMult[,,j,],
                                  dim = c(NMCSamples, NChains, NMicrobeNodes)),
                            warmup = warmup,
                            probs = c(0.05, 0.95))
            temp <- cbind(hostNode = colnames(hostAncestors[[i]])[[j]], temp)
            rownames(temp) <- rownames(microbeAncestors)
            allRes <- rbind(allRes, temp)
        }
        
        cat('microbeNode\t', file = file.path(currtabledir, 'allSummedPhyloVarianceEffects.txt'))
        write.table(allRes,
                    file   = file.path(currtabledir, 'allSummedPhyloVarianceEffects.txt'),
                    sep    = '\t',
                    quote  = F,
                    append = T)
        ##
        
        ## we can also sum the effects of just a node and its parent, as a bit of a middle-ground approach to see whether a given node is significantly different from its grandparent instead of its parent or instead of the overall mean
        hostParents[[i]] <- makeParentMat(NHostNodes,
                                          hostTreesSampled[[i]],
                                          NHostTips,
                                          hostTreesSampled[[i]]$tip.label)

        for(j in 1:NMCSamples) {
            for(k in 1:NChains) {
                matMult[j,k,,] <- hostParents[[i]] %*%
                                  rbind(c(0,
                                          phyloLogVarMultPrev[j,k,] * metaScales[j,k,1]),
                                        cbind(phyloLogVarMultADiv[j,k,] * metaScales[j,k,2],
                                              phyloLogVarMultRaw[j,k,,] * metaScales[j,k,3])) %*%
                                  microbeParentsT
            }
        }
        
        allRes <- NULL
        for(j in 1:NHostNodes) {
            temp <- monitor(array(matMult[,,j,],
                                  dim = c(NMCSamples, NChains, NMicrobeNodes)),
                            warmup = warmup,
                            probs = c(0.05, 0.95))
            temp <- cbind(hostNode = colnames(hostAncestors[[i]])[[j]], temp)
            rownames(temp) <- rownames(microbeAncestors)
            allRes <- rbind(allRes, temp)
        }
        
        cat('microbeNode\t', file = file.path(currtabledir, 'parentalSummedPhyloVarianceEffects.txt'))
        write.table(allRes,
                    file   = file.path(currtabledir, 'parentalSummedPhyloVarianceEffects.txt'),
                    sep    = '\t',
                    quote  = F,
                    append = T)
        ##
        
        ## see if any host clades have higher variance among their descendants
        sums <- summary(fit[[i]], pars = 'phyloLogVarMultADiv', probs = c(0.05,0.95), use_cache = F)
        rownames(sums$summary) <- colnames(hostAncestors[[i]])
        
        cat('hostNode\t', file = file.path(currtabledir, 'phylogeneticADivEffects.txt'))
        write.table(sums$summary,
                    file   = file.path(currtabledir, 'phylogeneticADivEffects.txt'),
                    sep    = '\t',
                    quote  = F,
                    append = T)
        ##
        
        ## see if any microbe clades have higher variance among their descendants
        sums <- summary(fit[[i]], pars = 'phyloLogVarMultPrev', probs = c(0.05,0.95), use_cache = F)
        rownames(sums$summary) <- rownames(microbeAncestors)
        
        cat('microbeNode\t', file = file.path(currtabledir, 'phylogeneticPrevalenceEffects.txt'))
        write.table(sums$summary,
                    file   = file.path(currtabledir, 'phylogeneticPrevalenceEffects.txt'),
                    sep    = '\t',
                    quote  = F,
                    append = T)
        ##
        
        ##
        effectNames <- c('microbePrevalence',
                         colnames(modelMat)[1:NEffects],
                         paste0('host_', colnames(hostAncestors[[i]])),
                         baseNames)
        
        matMult <- array(NA,
                         dim = c(NMCSamples,
                                 NChains,
                                 NHostNodes + NEffects + length(sumconts) + 1,
                                 NMicrobeNodes + 1),
                         dimnames = list(sample = NULL,
                                         chain = NULL,
                                         hostnode = effectNames,
                                         microbenode = c('alphaDiversity',
                                                         colnames(microbeAncestors))))
        ##
        
        ## build a temporary model matrix
        hostMat <- matrix(0, nrow = NHostNodes + NEffects + length(sumconts) + 1, ncol = NHostNodes + NEffects + length(sumconts) + 1)
        hostMat[1:(NEffects + 1), 1:(NEffects + 1)] <- diag(1, nrow = NEffects + 1)
        hostMat[(NEffects + 2):(NEffects + NHostNodes + 1), (NEffects + 2):(NEffects + NHostNodes + 1)] <- hostAncestors[[i]]
        hostMat[(NEffects + NHostNodes + 2):(NEffects + NHostNodes + length(sumconts) + 1), (NEffects + NHostNodes + 2):(NEffects + NHostNodes + length(sumconts) + 1)] <- diag(1, nrow = length(sumconts))
        ##
        
        ## since the model has many degrees of freedom, it's possible that certain regions of the tree have 'significant' effects but that those effects are 'implemented' via different, nearby nodes in each iteration. We can sum the effects of all ancestral nodes to see if a given node is consistently influenced by all of itself plus its ancestral terms, even if all of those terms are individually inconsistent
        for(j in 1:NMCSamples) {
            for(k in 1:NChains) {
                matMult[j,k,,] <- hostMat %*% rbind(scaledMicrobeNodeEffects[j,k,,], baseLevelEffects[j,k,,]) %*% cbind(1, microbeAncestorsT)
            }
        }
        
        allRes <- NULL
        for(j in 1:(NHostNodes + NEffects + length(sumconts) + 1)) {
            temp <- monitor(array(matMult[,,j,],
                                  dim = c(NMCSamples, NChains, NMicrobeNodes + 1)),
                            warmup = warmup,
                            probs = c(0.05, 0.95))
            temp <- cbind(hostEffect = effectNames[[j]], temp)
            rownames(temp) <- c('alphaDiversity', rownames(microbeAncestors))
            allRes <- rbind(allRes, temp)
        }
        
        cat('microbeNode\t', file = file.path(currtabledir, 'allSummedNodeEffects.txt'))
        write.table(allRes,
                    file   = file.path(currtabledir, 'allSummedNodeEffects.txt'),
                    sep    = '\t',
                    quote  = F,
                    append = T)
        ##
        
        ## we can also sum the effects of just a node and its parent, as a bit of a middle-ground approach to see whether a given node is significantly different from its grandparent instead of its parent or instead of the overall mean
        hostMat[(NEffects + 2):(NEffects + NHostNodes + 1), (NEffects + 2):(NEffects + NHostNodes + 1)] <- hostParents[[i]][,-1]
        
        for(j in 1:NMCSamples) {
            for(k in 1:NChains) {
                matMult[j,k,,] <- hostMat %*% rbind(scaledMicrobeNodeEffects[j,k,,], baseLevelEffects[j,k,,]) %*% cbind(1, microbeParentsT)
            }
        }
        
        allRes <- NULL
        for(j in 1:(NHostNodes + NEffects + length(sumconts) + 1)) {
            temp <- monitor(array(matMult[,,j,],
                                  dim = c(NMCSamples, NChains, NMicrobeNodes + 1)),
                            warmup = warmup,
                            probs = c(0.05, 0.95))
            temp <- cbind(hostNode = effectNames[[j]], temp)
            rownames(temp) <- c('alphaDiversity', rownames(microbeAncestors))
            allRes <- rbind(allRes, temp)
        }
        
        cat('microbeNode\t', file = file.path(currtabledir, 'parentalSummedNodeEffects.txt'))
        write.table(allRes,
                    file   = file.path(currtabledir, 'parentalSummedNodeEffects.txt'),
                    sep    = '\t',
                    quote  = F,
                    append = T)
        ##
    }
    ##

    ## summarize results for parameters that can be interpretted across all sampled host trees
    currplotdir <- file.path(outdir,'alltrees','plots')
    currtabledir <- file.path(outdir,'alltrees','tables')
    currdatadir <- file.path(outdir,'alltrees','data')

    dir.create(currplotdir, recursive=T)
    dir.create(currtabledir, recursive=T)
    dir.create(currdatadir, recursive=T)
    
    NSuccessTrees <- sum(fitModes == 0)
    
    if (NSuccessTrees > 1) {

        allfit <- sflist2stanfit(fit[fitModes == 0])

        stDProps <- extract(allfit, pars = 'stDProps')[[1]]
        colnames(stDProps) <- c(paste0('ADiv.', colnames(factLevelMat)),
                                paste0('Specificity.', colnames(factLevelMat)),
                                'ADiv.host',
                                'host.specificity',
                                'microbe.prevalence')

        pdf(file = file.path(currplotdir, 'scalesboxes.pdf'), width = 25, height = 15)
        boxplot(stDProps, cex.axis = 0.5, las = 2)
        graphics.off()

        save(stDProps, file = file.path(currdatadir, 'stDProps.RData'))

        relativeHostEvolRates <- exp(extract(allfit, pars='logRelativeHostEvolRates')[[1]])
        colnames(relativeHostEvolRates) <- sapply(1:(length(meanHostBoundariesRounded) + 1),
                                                  function(x) paste0(c('before', meanHostBoundariesRounded)[[x]],
                                                                     ' - ',
                                                                     c(meanHostBoundariesRounded, 'present')[[x]],
                                                                     ' mya'))

        pdf(file = file.path(currplotdir, 'evolRatesHost.pdf'), width = 25, height = 15)
        boxplot(relativeHostEvolRates, xlab = 'Time Period', ylab = 'Rate of Evolution Relative to Mean')
        graphics.off()

        pdf(file = file.path(currplotdir, 'logEvolRatesHost.pdf'), width = 25, height = 15)
        boxplot(log(relativeHostEvolRates), xlab = 'Time Period', ylab = 'Log Rate of Evolution Relative to Mean')
        graphics.off()

        save(relativeHostEvolRates, file = file.path(currdatadir, 'relativeHostEvolRates.RData'))

        ## compare rates of microbe evolution in each time bin
        relativeMicrobeEvolRates <- exp(extract(allfit, pars = 'logRelativeMicrobeEvolRates')[[1]])

        pdf(file = file.path(currplotdir, 'evolRatesMicrobe.pdf'), width = 25, height = 15)
        boxplot(relativeMicrobeEvolRates, xlab = 'Time Period', ylab = 'Rate of Evolution Relative to Mean')
        graphics.off()

        pdf(file = file.path(currplotdir, 'logEvolRatesMicrobe.pdf'), width = 25, height = 15)
        boxplot(log(relativeMicrobeEvolRates), xlab = 'Time Period', ylab = 'Log Rate of Evolution Relative to Mean')
        graphics.off()

        save(relativeMicrobeEvolRates, file = file.path(currdatadir, 'relativeMicrobeEvolRates.RData'))
        ##

        ## ornstein-uhlenbeck parameters
        hostOUAlpha <- extract(allfit, pars = 'hostOUAlpha')[[1]]
        microbeOUAlpha <- extract(allfit, pars = 'microbeOUAlpha')[[1]]
        OUAlphas <- cbind(hostOUAlpha, microbeOUAlpha)
        colnames(OUAlphas) <- c('host', 'microbe')

        pdf(file = file.path(currplotdir,'OUAlphas.pdf'), width = 25, height = 15)
        boxplot(OUAlphas, xlab = 'Host or microbe', ylab = 'alpha')
        graphics.off()

        save(OUAlphas, file = file.path(currdatadir, 'OUAlphas.RData'))
        ##

        ## summarize the mean branch lengths of the microbes
        sums <- summary(allfit, pars = 'microbeScales', probs = c(0.05,0.95), use_cache = F)
        newEdges <- sums$summary[,'mean']^2
        finalMicrobeTree.newEdges <- finalMicrobeTree
        finalMicrobeTree.newEdges$edge.length <- newEdges[order(microbeTreeDetails$edgeOrder)]
        pdf(file = file.path(currplotdir, 'microbeTreeWEstimatedEdgeLengths.pdf'), width = 25, height = 15)
        plot(finalMicrobeTree.newEdges, cex = 0.5)
        graphics.off()
        ##

        ## extract effects from model
        scaledMicrobeNodeEffects <- array(extract(allfit, pars='scaledMicrobeNodeEffects', permuted=F, inc_warmup=T),
                                          dim = c(NMCSamples,
                                                  NChains * NTrees,
                                                  NEffects + NHostNodes + 1,
                                                  NMicrobeNodes + 1),
                                          dimnames = list(sample  = NULL,
                                                          chain   = NULL,
                                                          effect  = c('microbePrevalence',
                                                                      colnames(modelMat)[1:NEffects],
                                                                      paste0('host_', colnames(hostAncestors[[i]]))),
                                                          taxnode = c('alphaDiversity', colnames(microbeAncestors))))
                                                        
        save(scaledMicrobeNodeEffects, file = file.path(currdatadir, 'scaledMicrobeNodeEffects.RData'))
        ##

        ## calculate base-level effects (negative sum of all others in category)
        baseLevelEffects <- array(NA,
                                  dim = c(NMCSamples,
                                          NChains * NTrees,
                                          length(sumconts),
                                          NMicrobeNodes + 1),
                                  dimnames = list(sample  = NULL,
                                                  chain   = NULL,
                                                  effect  = sumconts,
                                                  taxnode = c('alphaDiversity', colnames(microbeAncestors))))
        for(j in 1:NMCSamples) {
            for(k in 1:(NChains * NTrees)) {
                for(m in sumconts) {
                    baseLevelEffects[j,k,m,] <- -colSums(scaledMicrobeNodeEffects[j, k, rownames(factLevelMat)[factLevelMat[,m] == 1],])
                }
            }
        }

        save(baseLevelEffects, file = file.path(currdatadir, 'baseLevelEffects.RData'))
        ##

        ## summarize posterior distributions of effects
        allRes <- NULL
        for(j in 1:(NEffects + NHostNodes + 1)) {
            temp <- monitor(array(scaledMicrobeNodeEffects[,,j,],
                                  dim = c(NMCSamples, NChains * NTrees, NMicrobeNodes + 1)),
                            warmup = warmup,
                            probs  = c(0.05, 0.95))
            temp <- cbind(effect = dimnames(scaledMicrobeNodeEffects)[[3]][j], temp)
            rownames(temp) <- c('alphaDiversity', rownames(microbeAncestors))
            allRes <- rbind(allRes, temp)
        }
        ##

        ## summarize posterior distibutions of base-level effects
        for(m in sumconts) {
            temp <- monitor(array(baseLevelEffects[,,m,],
                                  dim = c(NMCSamples, NChains * NTrees, NMicrobeNodes + 1)),
                            warmup = warmup,
                            probs  = c(0.05, 0.95))
            temp <- cbind(effect = baseNames[[m]], temp)
            rownames(temp) <- c('alphaDiversity', rownames(microbeAncestors))
            allRes <- rbind(allRes, temp)
        }
        ##

        ## write posterior summaries of effects to file
        cat('microbeNode\t', file = file.path(currtabledir, 'nodeEffects.txt'))
        write.table(allRes,
                    file   = file.path(currtabledir, 'nodeEffects.txt'),
                    sep    = '\t',
                    quote  = F,
                    append = T)
        ##

        ##
        effectNames <- c('microbePrevalence',
                         colnames(modelMat)[1:NEffects],
                         paste0('host_', colnames(hostAncestors[[1]])),
                         baseNames)

        matMult <- array(NA,
                         dim = c(NMCSamples,
                                 NChains,
                                 NHostNodes + NEffects + length(sumconts) + 1,
                                 NMicrobeNodes + 1),
                         dimnames = list(sample = NULL,
                                         chain = NULL,
                                         hostnode = effectNames,
                                         microbenode = c('alphaDiversity',
                                                         colnames(microbeAncestors))))
        ##

        ## build a temporary model matrix
        hostMat <- matrix(0, nrow = NHostNodes + NEffects + length(sumconts) + 1, ncol = NHostNodes + NEffects + length(sumconts) + 1)
        hostMat[1:(NEffects + 1), 1:(NEffects + 1)] <- diag(1, nrow = NEffects + 1)
        hostMat[(NEffects + 2):(NEffects + NHostNodes + 1), (NEffects + 2):(NEffects + NHostNodes + 1)] <- hostAncestors[[1]]
        hostMat[(NEffects + NHostNodes + 2):(NEffects + NHostNodes + length(sumconts) + 1), (NEffects + NHostNodes + 2):(NEffects + NHostNodes + length(sumconts) + 1)] <- diag(1, nrow = length(sumconts))
        ##

        ## since the model has many degrees of freedom, it's possible that certain regions of the tree have 'significant' effects but that those effects are 'implemented' via different, nearby nodes in each iteration. We can sum the effects of all ancestral nodes to see if a given node is consistently influenced by all of itself plus its ancestral terms, even if all of those terms are individually inconsistent
        for(j in 1:NMCSamples) {
            for(k in 1:NChains) {
                matMult[j,k,,] <- hostMat %*%
                                  rbind(scaledMicrobeNodeEffects[j,k,,],
                                        baseLevelEffects[j,k,,]) %*%
                                  cbind(1, microbeAncestorsT)
            }
        }

        allRes <- NULL
        for(j in 1:(NHostNodes + NEffects + length(sumconts) + 1)) {
            temp <- monitor(array(matMult[,,j,],
                                  dim = c(NMCSamples, NChains, NMicrobeNodes + 1)),
                            warmup = warmup,
                            probs = c(0.05, 0.95))
            temp <- cbind(hostEffect = effectNames[[j]], temp)
            rownames(temp) <- c('alphaDiversity', rownames(microbeAncestors))
            allRes <- rbind(allRes, temp)
        }

        cat('microbeNode\t', file = file.path(currtabledir, 'allSummedNodeEffects.txt'))
        write.table(allRes,
                    file   = file.path(currtabledir, 'allSummedNodeEffects.txt'),
                    sep    = '\t',
                    quote  = F,
                    append = T)
        ##

        ## we can also sum the effects of just a node and its parent, as a bit of a middle-ground approach to see whether a given node is significantly different from its grandparent instead of its parent or instead of the overall mean
        hostMat[(NEffects + 2):(NEffects + NHostNodes + 1), (NEffects + 2):(NEffects + NHostNodes + 1)] <- hostParents[[1]][,-1]

        for(j in 1:NMCSamples) {
            for(k in 1:NChains) {
                matMult[j,k,,] <- hostMat %*%
                                  rbind(scaledMicrobeNodeEffects[j,k,,],
                                        baseLevelEffects[j,k,,]) %*%
                                  cbind(1, microbeParentsT)
            }
        }

        allRes <- NULL
        for(j in 1:(NHostNodes + NEffects + length(sumconts) + 1)) {
            temp <- monitor(array(matMult[,,j,],
                                  dim = c(NMCSamples, NChains, NMicrobeNodes + 1)),
                            warmup = warmup,
                            probs = c(0.05, 0.95))
            temp <- cbind(hostNode = effectNames[[j]], temp)
            rownames(temp) <- c('alphaDiversity', rownames(microbeAncestors))
            allRes <- rbind(allRes, temp)
        }

        cat('microbeNode\t', file = file.path(currtabledir, 'parentalSummedNodeEffects.txt'))
        write.table(allRes,
                    file   = file.path(currtabledir, 'parentalSummedNodeEffects.txt'),
                    sep    = '\t',
                    quote  = F,
                    append = T)
        ##

        allfit <- c(fit, allfit)
        
    }
}

makeDiagnosticPlots <- function(...) {
    
    ## make some diagnostic plots
    for(i in 1:(NTrees + 1)) {
        currdiagnosticdir <- file.path(outdir, 'diagnostics', if(i <= NTrees) {paste0('tree_', i)} else {'allfit'})
        dir.create(currdiagnosticdir, recursive = T)
        
        pars <- c('stDProps', 'hostMetaVarProps', 'metaVarProps', 'aveStD')
        
        pdf(file=file.path(currdiagnosticdir,'pairs_grabBag.pdf'), width = 50, height = 50)
        pairs(allfit[[i]], pars = pars)
        graphics.off()
        
        pdf(file=file.path(currdiagnosticdir,'pairs_grabBagLog.pdf'), width = 50, height = 50)
        pairs(allfit[[i]], pars = pars, log = T)
        graphics.off()
        
        pars <- paste0('phyloLogVarMultPrev[',1:20,']')
        
        pdf(file=file.path(currdiagnosticdir,'pairs_phyloLogVarMultPrev.pdf'), width = 50, height = 50)
        pairs(allfit[[i]], pars = pars, log = T)
        graphics.off()
        
    }
    ##
}
## fin
