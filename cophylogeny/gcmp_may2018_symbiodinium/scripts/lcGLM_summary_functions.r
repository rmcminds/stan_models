makeParentMat <- function(NNodes, tree, NTips, tipNames) {
    parentMat <- matrix(0, NNodes + 1, NNodes + 1)
    grandparentMat <- matrix(0, NNodes + 1, NNodes + 1)
    for(node in 1:(NNodes + 1)) {
        p <- Ancestors(tree, node, 'parent')
        parentMat[node, ] <- as.numeric(1:(NNodes + 1) %in% c(p, node))
        grandparentMat[node, ] <- as.numeric(1:(NNodes + 1) %in% c(Ancestors(tree, p, 'parent'), p, node))
    }
    colnames(parentMat) <- rownames(parentMat) <- paste0('i', 1:(NNodes + 1))
    colnames(parentMat)[1:NTips] <- rownames(parentMat)[1:NTips] <- paste0('t', tipNames)
    parentMat <- cbind(0, parentMat[-(NTips + 1), -(NTips + 1)])
    parentMat <- rbind(c(1, rep(0, ncol(parentMat) - 1)), parentMat)
    
    colnames(grandparentMat) <- rownames(grandparentMat) <- paste0('i', 1:(NNodes + 1))
    colnames(grandparentMat)[1:NTips] <- rownames(grandparentMat)[1:NTips] <- paste0('t', tipNames)
    grandparentMat <- cbind(0, grandparentMat[-(NTips + 1), -(NTips + 1)])
    grandparentMat <- rbind(c(1, rep(0, ncol(grandparentMat) - 1)), grandparentMat)
    
    return(list(parentMat      = parentMat,
                grandparentMat = grandparentMat))
}

statusUpdate <- function(iter, N) {
    if(iter == round(N / 4)) {
        cat('\t25%...')
    } else if(iter == round(N / 2)) {
        cat('50%...')
    } else if (iter == round(3 * N / 4)) {
        cat('75%...')
    } else if (iter == N) {
        cat('100%\n')
    }
}

summarizeLcGLM <- function(combineTrees = T, separateTrees = T, ...) {
    
    fitModes <- sapply(fit, function(x) x@mode)
    
    ## create parental ancestry matrix for microbes (only sums effect of a given node and its direct parent, not all ancestors
    microbeP.GP <- makeParentMat(NMicrobeNodes,
                                 finalMicrobeTree,
                                 NMicrobeTips,
                                 microbeTips)
                                 
    microbeParentsT <- t(microbeP.GP$parentMat)
    microbeGrandparentsT <- t(microbeP.GP$grandparentMat)

    hostP.GP <- list()
    hostP.GP[[1]] <- makeParentMat(NHostNodes,
                                   hostTreesSampled[[1]],
                                   NHostTips,
                                   hostTreesSampled[[1]]$tip.label)
    
    hostParents <- list()
    hostParents[[1]] <- hostP.GP[[1]]$parentMat
    hostGrandparents <- list()
    hostGrandparents[[1]] <- hostP.GP[[1]]$grandparentMat
    
    microbeAncestorsT <- t(microbeAncestors)
    microbeMat <- matrix(0, nrow = NMicrobeNodes + 1, ncol = NMicrobeNodes + 1)
    microbeMat[1, 1] <- 1
    microbeMat[2:(NMicrobeNodes + 1), 2:(NMicrobeNodes + 1)] <- microbeAncestorsT

    colorpal <- colorRampPalette(brewer.pal(9, 'Blues'))
    plotcolors <- c('white', colorpal(100), 'black')
    
    baseNames <- sapply(sumconts, function(m) paste0(m, levels(newermap[,m])[nlevels(newermap[,m])]))

    ## summarize results for parameters that can be interpretted across all sampled host trees
    NSuccessTrees <- sum(fitModes == 0)

    ## summarize the results separately for each sampled host tree and/or for all trees combined
    for(i in 1:(separateTrees * NTrees + combineTrees)) {
        
        if (NSuccessTrees > 1 & i == (separateTrees * NTrees + combineTrees)) {
    
            currplotdir <- file.path(subdir, 'alltrees', 'plots')
            currtabledir <- file.path(subdir, 'alltrees', 'tables')
            currdatadir <- file.path(subdir, 'alltrees', 'data')
            
            dir.create(currplotdir, recursive = T)
            dir.create(currtabledir, recursive = T)
            dir.create(currdatadir, recursive = T)
            
            cat('\nSummarizing effects of all trees combined:\n')
            cat(paste0(as.character(Sys.time()), '\n'))

            i <- 1
            fit[[i]] <- sflist2stanfit(fit[fitModes == 0])
            NChains <- NSuccessTrees * NChains
            
            fit[-i] <- NA
            gc()
        
        } else if (fitModes[[i]] != 0) {
            
            cat(paste0('\nTree ', i, ' did not complete successfully\n'))
            cat(paste0(as.character(Sys.time()), '\n'))
            next
            
        } else {
            
            currplotdir <- file.path(subdir, paste0('tree_', i), 'plots')
            currtabledir <- file.path(subdir, paste0('tree_', i), 'tables')
            currdatadir <- file.path(subdir, paste0('tree_', i), 'data')
            
            dir.create(currplotdir, recursive = T)
            dir.create(currtabledir, recursive = T)
            dir.create(currdatadir, recursive = T)
            
            cat(paste0('\nSummarizing Tree ', i, ':\n'))
            cat(paste0(as.character(Sys.time()), '\n'))
            
        }
        
        sink(stdout(), type = "message")
        check_hmc_diagnostics(fit[[i]])
        sink(NULL, type = "message")

        ## plot the sampled tree with the time bins marked
        pdf(file   = file.path(currplotdir, 'sampledHostTree.pdf'),
            width  = 25,
            height = 15)
        plot(hostTreesSampled[[i]],
             cex = 0.75)
        graphics.off()
        ##
        
        ## plot the sampled tree with the time bins marked
        pdf(file   = file.path(currplotdir, 'sampledMicrobeTree.pdf'),
            width  = 25,
            height = 15)
        plot(finalMicrobeTree,
             cex = 0.75)
        graphics.off()
        ##
        
        ## summarize the mean branch lengths of the microbes
        newEdges <- apply(array(extract(fit[[i]],
                                        pars       = 'microbeScales',
                                        permuted   = F),
                                dim = c(NMCSamples - warmup,
                                        NChains,
                                        NMicrobeNodes))^2,
                          3,
                          mean)
        finalMicrobeTree.newEdges <- finalMicrobeTree
        finalMicrobeTree.newEdges$edge.length <- newEdges[order(microbeTreeDetails$edgeOrder)]
        pdf(file   = file.path(currplotdir, 'microbeTreeWEstimatedEdgeLengths.pdf'),
            width  = 25,
            height = 15)
        plot(finalMicrobeTree.newEdges,
             cex = 0.5)
        graphics.off()
        ##
        
        ## summarize the mean branch lengths of the hosts
        newEdges <- apply(array(extract(fit[[i]],
                                        pars       = 'hostScales',
                                        permuted   = F),
                                dim = c(NMCSamples - warmup,
                                        NChains,
                                        NHostNodes))^2,
                          3,
                          mean)
        hostTreesSampled.newEdges <- hostTreesSampled[[i]]
        hostTreesSampled.newEdges$edge.length <- newEdges[order(hostTreeDetails[[i]]$edgeOrder)]
        pdf(file   = file.path(currplotdir, 'hostTreeWEstimatedEdgeLengths.pdf'),
            width  = 25,
            height = 15)
        plot(hostTreesSampled.newEdges,
             cex = 0.5)
        graphics.off()
        ##
        
        ## plot heatmap of cophylogenetic patterns
        plotmicrobetree <- ladderize(multi2di(finalMicrobeTree))
        hclmicrobetree <- as.hclust(plotmicrobetree)
        # for each sample in the data, assign its mitotype to a vector
        hostvect <- sampleMap[[i]][,sampleTipKey]
        # name the vector with the sample IDs
        names(hostvect) <- rownames(sampleMap[[i]])
        # sort the samples in the vector
        temp <- sampleMap[[i]][order(sampleMap[[i]]$concatenated_date),]
        temp <- temp[order(temp$reef_name),]
        temp <- temp[order(temp[,sampleTipKey]),]
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

        cat('\n\tRaw variance estimates\n\t')
        cat(paste0(as.character(Sys.time()), '\n'))
        
        ## variance partitioning
        stDProps <- array(extract(fit[[i]],
                                  pars       = 'stDProps',
                                  permuted   = F,
                                  inc_warmup = T),
                          dim = c(NMCSamples,
                                  NChains,
                                  2 * NFactors + 3),
                          dimnames = list(sample  = NULL,
                                          chain   = NULL,
                                          factor  = c(paste0('ADiv.', colnames(factLevelMat)),
                                                      paste0('Specificity.', colnames(factLevelMat)),
                                                      'ADiv.host',
                                                      'host.specificity',
                                                      'microbe.prevalence')))
        save(stDProps, file = file.path(currdatadir, 'stDProps.RData'))
        
        stDPropsPlot <- NULL
        for(j in 1:NChains) {
            stDPropsPlot <- rbind(stDPropsPlot, stDProps[(warmup + 1):NMCSamples, j,])
        }
        colnames(stDPropsPlot) <- gsub('_',
                                       ' ',
                                       c(paste0(colnames(factLevelMat), ' (ADiv)'),
                                         paste0(colnames(factLevelMat), ' (Specificity)'),
                                         'Host (ADiv)',
                                         'Host (Specificity)',
                                         'Microbe prevalence'))
                                
        ADivInd <- c(1:NFactors, 2 * NFactors + 1)
        specInd <- c((NFactors + 1):(2 * NFactors), 2 * NFactors + 2)
        meds <- apply(stDPropsPlot, 2, median)
        stDPropsPlot <- stDPropsPlot[, c(2 * NFactors + 3,
                                         ADivInd[order(meds[ADivInd], decreasing = T)],
                                         specInd[order(meds[specInd], decreasing = T)])]
                                   
        pdf(file   = file.path(currplotdir, 'stDProps_boxes.pdf'),
            width  = 7,
            height = 7)
        boxplot(stDPropsPlot,
                at       = c(1, 3:(NFactors + 3), (NFactors + 5):(2 * NFactors + 5)),
                cex.axis = 0.5,
                las      = 2,
                xlab     = 'Factor',
                ylab     = 'Percent of total variance')
        graphics.off()
        ##
        
        cat('\n\tRaw meta-variance estimates\n\t')
        cat(paste0(as.character(Sys.time()), '\n'))
        
        ## meta-variance partitioning
        metaVarProps <- array(extract(fit[[i]],
                                      pars       = 'metaVarProps',
                                      permuted   = F,
                                      inc_warmup = T),
                              dim = c(NMCSamples,
                                      NChains,
                                      3),
                              dimnames = list(sample  = NULL,
                                              chain   = NULL,
                                              effect  = c('Prevalence',
                                                          'ADiv',
                                                          'Specificty')))
        metaVarPropsPlot <- NULL
        for(j in 1:NChains) {
            metaVarPropsPlot <- rbind(metaVarPropsPlot, metaVarProps[(warmup + 1):NMCSamples, j,])
        }
        pdf(file   = file.path(currplotdir, 'metaVarProps_boxes.pdf'),
            width  = 7,
            height = 7)
        boxplot(metaVarPropsPlot,
                cex.axis = 0.5,
                las      = 2)
        graphics.off()
        
        save(metaVarProps, file = file.path(currdatadir, 'metaVarProps.RData'))
        ##
        
        ## actual scales of metavariance
        metaScales <- array(extract(fit[[i]],
                                    pars       = 'metaScales',
                                    permuted   = F,
                                    inc_warmup = T),
                            dim = c(NMCSamples,
                                    NChains,
                                    3),
                            dimnames = list(sample  = NULL,
                                            chain   = NULL,
                                            effect  = c('Prevalence',
                                                        'ADiv',
                                                        'Specificty')))
        metaScalesPlot <- NULL
        for(j in 1:NChains) {
            metaScalesPlot <- rbind(metaScalesPlot, metaScales[(warmup + 1):NMCSamples, j,])
        }
        pdf(file   = file.path(currplotdir, 'metaScales_boxes.pdf'),
            width  = 7,
            height = 7)
        boxplot(metaScalesPlot,
                cex.axis = 0.5,
                las      = 2)
        graphics.off()
        
        save(metaScales, file = file.path(currdatadir, 'metaScales.RData'))
        ##
        
        ## ornstein-uhlenbeck parameters
        OUAlphas <- array(extract(fit[[i]],
                                  pars       = c('hostOUAlpha', 'microbeOUAlpha'),
                                  permuted   = F,
                                  inc_warmup = T),
                          dim = c(NMCSamples,
                                  NChains,
                                  2),
                          dimnames = list(sample = NULL,
                                          chain  = NULL,
                                          alpha  = c('host', 'microbe')))
        OUAlphaPlot <- NULL
        for(j in 1:NChains) {
            OUAlphaPlot <- rbind(OUAlphaPlot, OUAlphas[(warmup + 1):NMCSamples, j,])
        }
        pdf(file   = file.path(currplotdir, 'OUAlphas.pdf'),
            width  = 7,
            height = 7)
        boxplot(OUAlphaPlot,
                xlab = 'Host or microbe',
                ylab = 'alpha')
        graphics.off()
        
        save(OUAlphas, file = file.path(currdatadir, 'OUAlphas.RData'))
        ##
        
        cat('\n\tRaw variance modifiers\n\t')
        cat(paste0(as.character(Sys.time()), '\n'))
        
        ## extract variance modifier terms from the fit model
        phyloLogVarMultPrev <- array(extract(fit[[i]],
                                             pars       = 'phyloLogVarMultPrev',
                                             permuted   = F,
                                             inc_warmup = T),
                                     dim = c(NMCSamples,
                                             NChains,
                                             NMicrobeNodes),
                                     dimnames = list(sample  = NULL,
                                                     chain   = NULL,
                                                     taxnode = colnames(microbeAncestors)))
        
        phyloLogVarMultADiv <- array(extract(fit[[i]],
                                             pars       = 'phyloLogVarMultADiv',
                                             permuted   = F,
                                             inc_warmup = T),
                                     dim = c(NMCSamples,
                                             NChains,
                                             NHostNodes),
                                     dimnames = list(sample  = NULL,
                                                     chain   = NULL,
                                                     taxnode = colnames(hostAncestors[[i]])))
        
        phyloLogVarMultRaw <- array(extract(fit[[i]],
                                            pars       = 'phyloLogVarMultRaw',
                                            permuted   = F,
                                            inc_warmup = T),
                                    dim = c(NMCSamples,
                                            NChains,
                                            NHostNodes,
                                            NMicrobeNodes),
                                    dimnames = list(sample      = NULL,
                                                    chain       = NULL,
                                                    hostnode    = colnames(hostAncestors[[i]]),
                                                    microbenode = colnames(microbeAncestors)))
        ##
        
        phyloLogVarMultScaled <- array(NA,
                                       dim = c(NMCSamples,
                                               NChains,
                                               NHostNodes + 1,
                                               NMicrobeNodes + 1))
        for(j in 1:NMCSamples) {
            for(k in 1:NChains) {
                phyloLogVarMultScaled[j,k,,] <- rbind(c(0,
                                                       phyloLogVarMultPrev[j,k,] * metaScales[j,k,1]),
                                                     cbind(phyloLogVarMultADiv[j,k,] * metaScales[j,k,2],
                                                           phyloLogVarMultRaw[j,k,,] * metaScales[j,k,3]))
            }
        }
        
        ## see if any clades have higher variance among their descendants (maybe suggesting codiversification)
        allRes <- NULL
        for(j in 1:(NHostNodes + 1)) {
            temp <- monitor(array(phyloLogVarMultScaled[,,j,],
                                  dim = c(NMCSamples,
                                          NChains,
                                          NMicrobeNodes + 1)),
                            warmup = warmup,
                            probs  = c(0.05, 0.95),
                            print  = F)
            temp <- cbind(hostNode = c('microbePrevalence',
                                       paste0('host_', colnames(hostAncestors[[i]])))[[j]], temp)
            rownames(temp) <- c('alphaDiversity', rownames(microbeAncestors))
            allRes <- rbind(allRes, temp)
        }
        
        cat('microbeNode\t', file = file.path(currtabledir, 'individualPhyloVarianceEffects.txt'))
        write.table(allRes,
                    file   = file.path(currtabledir, 'individualPhyloVarianceEffects.txt'),
                    sep    = '\t',
                    quote  = F,
                    append = T)
        ##
        
        cat('\n\tSummed variance modifiers\n\t')
        cat(paste0(as.character(Sys.time()), '\n'))
        
        ## since the model has many degrees of freedom, it's possible that certain regions of the tree have 'significant' effects but that those effects are 'implemented' via different, nearby nodes in each iteration. We can sum the effects of all ancestral nodes to see if a given node is consistently influenced by all of itself plus its ancestral terms, even if all of those terms are individually inconsistent
        ## first build a temporary model matrix
        hostMat <- matrix(0,
                          nrow = NHostNodes + NEffects + NSumTo0 + 1,
                          ncol = NHostNodes + NEffects + NSumTo0 + 1)
        hostMat[1:(NEffects + 1),
                1:(NEffects + 1)] <- diag(1, nrow = NEffects + 1)
        hostMat[(NEffects + 2):(NEffects + NHostNodes + 1),
                (NEffects + 2):(NEffects + NHostNodes + 1)] <- hostAncestors[[i]]
        hostMat[(NEffects + NHostNodes + 2):(NEffects + NHostNodes + NSumTo0 + 1),
                (NEffects + NHostNodes + 2):(NEffects + NHostNodes + NSumTo0 + 1)] <- diag(1, nrow = NSumTo0)
        ##
        
        ## sum the effects
        matMult <- array(NA,
             dim = c(NMCSamples,
                     NChains,
                     NHostNodes + 1,
                     NMicrobeNodes + 1),
             dimnames = list(sample      = NULL,
                             chain       = NULL,
                             hostnode    = c('microbePrevalence', colnames(hostAncestors[[i]])),
                             microbenode = c('alphaDiversity', colnames(microbeAncestors))))
                             
        for(j in 1:NMCSamples) {
            for(k in 1:NChains) {
                matMult[j,k,,] <- hostMat[c(1, (NEffects + 2):(NEffects + NHostNodes + 1)),
                                          c(1, (NEffects + 2):(NEffects + NHostNodes + 1))] %*%
                                  phyloLogVarMultScaled[j,k,,] %*%
                                  microbeMat
            }
        }
        
        allRes <- NULL
        for(j in 1:(NHostNodes + 1)) {
            temp <- monitor(array(matMult[,,j,],
                                  dim = c(NMCSamples,
                                          NChains,
                                          NMicrobeNodes + 1)),
                            warmup = warmup,
                            probs  = c(0.05, 0.95),
                            print  = F)
            temp <- cbind(hostNode = c('microbePrevalence',
                                       paste0('host_', colnames(hostAncestors[[i]])))[[j]], temp)
            rownames(temp) <- c('alphaDiversity', rownames(microbeAncestors))
            allRes <- rbind(allRes, temp)
            statusUpdate(j, NHostNodes)
        }
        
        cat('microbeNode\t', file = file.path(currtabledir, 'allSummedPhyloVarianceEffects.txt'))
        write.table(allRes,
                    file   = file.path(currtabledir, 'allSummedPhyloVarianceEffects.txt'),
                    sep    = '\t',
                    quote  = F,
                    append = T)
        ##
        
        cat('\n\tSums of node and parent node variance modifiers\n\t')
        cat(paste0(as.character(Sys.time()), '\n'))
        
        ## we can also sum the effects of just a node and its parent, as a bit of a middle-ground approach to see whether a given node is significantly different from its grandparent instead of its parent or instead of the overall mean
        hostP.GP[[i]] <- makeParentMat(NHostNodes,
                                       hostTreesSampled[[i]],
                                       NHostTips,
                                       hostTreesSampled[[i]]$tip.label)
                                       
        hostParents[[i]] <- hostP.GP[[i]]$parentMat

        for(j in 1:NMCSamples) {
            for(k in 1:NChains) {
                matMult[j,k,,] <- hostParents[[i]] %*%
                                  phyloLogVarMultScaled[j,k,,] %*%
                                  microbeParentsT
            }
        }
        
        allRes <- NULL
        for(j in 1:(NHostNodes + 1)) {
            temp <- monitor(array(matMult[,,j,],
                                  dim = c(NMCSamples,
                                          NChains,
                                          NMicrobeNodes + 1)),
                            warmup = warmup,
                            probs  = c(0.05, 0.95),
                            print  = F)
            temp <- cbind(hostNode = c('microbePrevalence',
                                       paste0('host_', colnames(hostAncestors[[i]])))[[j]], temp)
            rownames(temp) <- c('alphaDiversity', rownames(microbeAncestors))
            allRes <- rbind(allRes, temp)
            statusUpdate(j, NHostNodes)
        }
        
        cat('microbeNode\t', file = file.path(currtabledir, 'parentalSummedPhyloVarianceEffects.txt'))
        write.table(allRes,
                    file   = file.path(currtabledir, 'parentalSummedPhyloVarianceEffects.txt'),
                    sep    = '\t',
                    quote  = F,
                    append = T)
        ##
        
        cat('\n\tSums of node, parent, and grandparent node variance modifiers\n\t')
        cat(paste0(as.character(Sys.time()), '\n'))
        
        ## or a node, its parent, and its grandparent
        hostGrandparents[[i]] <- hostP.GP[[i]]$grandparentMat
        
        for(j in 1:NMCSamples) {
            for(k in 1:NChains) {
                matMult[j,k,,] <- hostGrandparents[[i]] %*%
                                  phyloLogVarMultScaled[j,k,,] %*%
                                  microbeGrandparentsT
            }
        }
        
        allRes <- NULL
        for(j in 1:(NHostNodes + 1)) {
            temp <- monitor(array(matMult[,,j,],
                                  dim = c(NMCSamples,
                                          NChains,
                                          NMicrobeNodes + 1)),
                            warmup = warmup,
                            probs  = c(0.05, 0.95),
                            print  = F)
            temp <- cbind(hostNode = c('microbePrevalence',
                                       paste0('host_', colnames(hostAncestors[[i]])))[[j]], temp)
            rownames(temp) <- c('alphaDiversity', rownames(microbeAncestors))
            allRes <- rbind(allRes, temp)
            statusUpdate(j, NHostNodes)
        }
        
        cat('microbeNode\t', file = file.path(currtabledir, 'grandparentalSummedPhyloVarianceEffects.txt'))
        write.table(allRes,
                    file   = file.path(currtabledir, 'grandparentalSummedPhyloVarianceEffects.txt'),
                    sep    = '\t',
                    quote  = F,
                    append = T)
        ##
        
        cat('\n\tRaw node effects\n\t')
        cat(paste0(as.character(Sys.time()), '\n'))
        
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
                                                                      colnames(modelMat)[2:(NEffects + 1)],
                                                                      paste0('host_', colnames(hostAncestors[[i]]))),
                                                          taxnode = c('alphaDiversity', colnames(microbeAncestors))))
                                                        
        save(scaledMicrobeNodeEffects, file = file.path(currdatadir, 'scaledMicrobeNodeEffects.RData'))
        ##
        
        ## calculate base-level effects (negative sum of all others in category)
        baseLevelEffects <- array(extract(fit[[i]],
                                          pars       = 'baseLevelEffects',
                                          permuted   = F,
                                          inc_warmup = T),
                                  dim = c(NMCSamples,
                                          NChains,
                                          NSumTo0,
                                          NMicrobeNodes + 1),
                                  dimnames = list(sample  = NULL,
                                                  chain   = NULL,
                                                  effect  = baseNames,
                                                  taxnode = c('alphaDiversity', colnames(microbeAncestors))))
        
        save(baseLevelEffects, file = file.path(currdatadir, 'baseLevelEffects.RData'))
        ##

        ## summarize posterior distributions of effects
        allRes <- NULL
        for(j in 1:(NEffects + NHostNodes + 1)) {
            temp <- monitor(array(scaledMicrobeNodeEffects[,,j,],
                                  dim = c(NMCSamples,
                                          NChains,
                                          NMicrobeNodes + 1)),
                            warmup = warmup,
                            probs  = c(0.05, 0.95),
                            print  = F)
            temp <- cbind(effect = dimnames(scaledMicrobeNodeEffects)[[3]][j], temp)
            rownames(temp) <- c('alphaDiversity', rownames(microbeAncestors))
            allRes <- rbind(allRes, temp)
            statusUpdate(j, NEffects + NHostNodes + 1)
        }
        ##
        
        ## summarize posterior distibutions of base-level effects
        for(m in baseNames) {
            temp <- monitor(array(baseLevelEffects[,,m,],
                                  dim = c(NMCSamples,
                                          NChains,
                                          NMicrobeNodes + 1)),
                            warmup = warmup,
                            probs  = c(0.05, 0.95),
                            print  = F)
            temp <- cbind(effect = m, temp)
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
                         colnames(modelMat)[2:(NEffects + 1)],
                         paste0('host_', colnames(hostAncestors[[i]])),
                         baseNames)
        
        matMult <- array(NA,
                         dim = c(NMCSamples,
                                 NChains,
                                 NHostNodes + NEffects + NSumTo0 + 1,
                                 NMicrobeNodes + 1),
                         dimnames = list(sample = NULL,
                                         chain = NULL,
                                         hostnode = effectNames,
                                         microbenode = c('alphaDiversity',
                                                         colnames(microbeAncestors))))
        ##
        
        cat('\n\tSummed node effects\n\t')
        cat(paste0(as.character(Sys.time()), '\n'))
        
        ## since the model has many degrees of freedom, it's possible that certain regions of the tree have 'significant' effects but that those effects are 'implemented' via different, nearby nodes in each iteration. We can sum the effects of all ancestral nodes to see if a given node is consistently influenced by all of itself plus its ancestral terms, even if all of those terms are individually inconsistent
        for(j in 1:NMCSamples) {
            for(k in 1:NChains) {
                matMult[j,k,,] <- hostMat %*%
                                  rbind(scaledMicrobeNodeEffects[j,k,,],
                                        baseLevelEffects[j,k,,]) %*%
                                  microbeMat
            }
        }
        
        allRes <- NULL
        for(j in 1:(NHostNodes + NEffects + NSumTo0 + 1)) {
            temp <- monitor(array(matMult[,,j,],
                                  dim = c(NMCSamples,
                                          NChains,
                                          NMicrobeNodes + 1)),
                            warmup = warmup,
                            probs  = c(0.05, 0.95),
                            print  = F)
            temp <- cbind(hostEffect = effectNames[[j]], temp)
            rownames(temp) <- c('alphaDiversity', rownames(microbeAncestors))
            allRes <- rbind(allRes, temp)
            statusUpdate(j, NHostNodes + NEffects + NSumTo0 + 1)
        }
        
        cat('microbeNode\t', file = file.path(currtabledir, 'allSummedNodeEffects.txt'))
        write.table(allRes,
                    file   = file.path(currtabledir, 'allSummedNodeEffects.txt'),
                    sep    = '\t',
                    quote  = F,
                    append = T)
        ##
        
        cat('\n\tSums of node and parent node effects\n\t')
        cat(paste0(as.character(Sys.time()), '\n'))
        
        ## we can also sum the effects of just a node and its parent, as a bit of a middle-ground approach to see whether a given node is significantly different from its grandparent instead of its parent or instead of the overall mean
        hostParentMat <- hostMat
        hostParentMat[(NEffects + 2):(NEffects + NHostNodes + 1),
                      (NEffects + 2):(NEffects + NHostNodes + 1)] <- hostParents[[i]][-1, -1]
        
        for(j in 1:NMCSamples) {
            for(k in 1:NChains) {
                matMult[j,k,,] <- hostParentMat %*%
                                  rbind(scaledMicrobeNodeEffects[j,k,,],
                                        baseLevelEffects[j,k,,]) %*%
                                  microbeParentsT
            }
        }
        
        allRes <- NULL
        for(j in 1:(NHostNodes + NEffects + NSumTo0 + 1)) {
            temp <- monitor(array(matMult[,,j,],
                                  dim = c(NMCSamples,
                                          NChains,
                                          NMicrobeNodes + 1)),
                            warmup = warmup,
                            probs  = c(0.05, 0.95),
                            print  = F)
            temp <- cbind(hostNode = effectNames[[j]], temp)
            rownames(temp) <- c('alphaDiversity', rownames(microbeAncestors))
            allRes <- rbind(allRes, temp)
            statusUpdate(j, NHostNodes + NEffects + NSumTo0 + 1)
        }
        
        cat('microbeNode\t', file = file.path(currtabledir, 'parentalSummedNodeEffects.txt'))
        write.table(allRes,
                    file   = file.path(currtabledir, 'parentalSummedNodeEffects.txt'),
                    sep    = '\t',
                    quote  = F,
                    append = T)
        ##
        
        cat('\n\tSums of node, parent, and grandparent node effects\n\t')
        cat(paste0(as.character(Sys.time()), '\n'))
        
        ## or a node, its parent, and its grandparent
        hostGrandparentMat <- hostMat
        hostGrandparentMat[(NEffects + 2):(NEffects + NHostNodes + 1),
                           (NEffects + 2):(NEffects + NHostNodes + 1)] <- hostGrandparents[[i]][-1, -1]
        
        for(j in 1:NMCSamples) {
            for(k in 1:NChains) {
                matMult[j,k,,] <- hostGrandparentMat %*%
                                  rbind(scaledMicrobeNodeEffects[j,k,,],
                                        baseLevelEffects[j,k,,]) %*%
                                  microbeGrandparentsT
            }
        }
        
        allRes <- NULL
        for(j in 1:(NHostNodes + NEffects + NSumTo0 + 1)) {
            temp <- monitor(array(matMult[,,j,],
                                  dim = c(NMCSamples,
                                          NChains,
                                          NMicrobeNodes + 1)),
                            warmup = warmup,
                            probs  = c(0.05, 0.95),
                            print  = F)
            temp <- cbind(hostNode = effectNames[[j]], temp)
            rownames(temp) <- c('alphaDiversity', rownames(microbeAncestors))
            allRes <- rbind(allRes, temp)
            statusUpdate(j, NHostNodes + NEffects + NSumTo0 + 1)
        }
        
        cat('microbeNode\t', file = file.path(currtabledir, 'grandparentalSummedNodeEffects.txt'))
        write.table(allRes,
                    file   = file.path(currtabledir, 'grandparentalSummedNodeEffects.txt'),
                    sep    = '\t',
                    quote  = F,
                    append = T)
        ##
    }
}

makeDiagnosticPlots <- function(...) {
    
    allfit <- c(fit, allfit)
    
    ## make some diagnostic plots
    for(i in 1:(NTrees + 1)) {
        currdiagnosticdir <- file.path(subdir, 'diagnostics', if(i <= NTrees) {paste0('tree_', i)} else {'allfit'})
        dir.create(currdiagnosticdir, recursive = T)
        
        pars <- c('stDProps', 'metaVarProps', 'aveStD')
        
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

runStanModel <- function(noData = F, shuffleData = F, ...) {
    
    if(sum(noData, shuffleData) > 1) {
        cat('\at most one of noData and shuffleData can be TRUE.\n')
        q()
    } else if(noData) {

        ## erase data from stan list
        for (i in 1:NTrees) {
            standat[[i]]$NObs <- 0
            standat[[i]]$present <- standat[[i]]$sampleNames <- standat[[i]]$microbeTipNames <- vector()
        }
        
        subdir <- file.path(outdir, 'resampled_from_prior')
        dir.create(subdir, recursive = T)
        save.image(file.path(subdir, 'setup.RData'))

        cat('\nSampling from prior\n')
    } else if(shuffleData) {

        ## shuffle data from stan list
        for (i in 1:NTrees) {
            standat[[i]]$present <- sample(present)
        }
        
        subdir <- file.path(outdir, 'resampled_with_randomObs')
        dir.create(subdir, recursive = T)
        save.image(file.path(subdir, 'setup.RData'))

        cat('\nSampling from shuffled data\n')
    } else {
        
        subdir <- file.path(outdir, 'primary_sampling')
        dir.create(subdir, recursive = T)
        save.image(file.path(subdir, 'setup.RData'))
        
        cat('\nSampling from model\n')
    }
    
    cat(paste0(as.character(Sys.time()), '\n'))

    ## run the model!
    fit <- mclapply(1:NTrees,
        function(i) {
            setTimeLimit(timeLimit)
            tryCatch({
                stan(file     = modelPath,
                     data     = standat[[i]],
                     control  = list(adapt_delta   = adapt_delta,
                                     max_treedepth = max_treedepth),
                     iter     = NIterations,
                     thin     = thin,
                     chains   = NChains,
                     seed     = seed,
                     chain_id = (NChains * (i - 1) + (1:NChains)),
                     pars     = c('rawMicrobeNodeEffects'),
                     include  = FALSE,
                     init_r   = init_r)
            }, error = function(e) NA)
        }, mc.preschedule = F,
           mc.cores       = NCores)

    cat('\nSaving results\n')
    cat(paste0(as.character(Sys.time()), '\n'))

    save(fit, file = file.path(subdir, 'fit.RData'))
    ##

    ## summarize the model fit
    summarizeLcGLM()
}
## fin
