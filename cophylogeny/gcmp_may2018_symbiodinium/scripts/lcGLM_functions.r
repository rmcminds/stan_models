logit <- function(x) { log(x / (1 - x)) }
inv_logit <- function(x) { 1 / (1 + exp(-x)) }

createAncestryMat <- function(NNodes, tree, NTips, tipNames) {
    ancestryMat <- matrix(0, NNodes + 1, NNodes + 1)
    for(node in 1:(NNodes + 1)) {
        ancestryMat[node, ] <- as.numeric(1:(NNodes + 1) %in% c(Ancestors(tree, node), node))
    }
    colnames(ancestryMat) <- rownames(ancestryMat) <- paste0('i', 1:(NNodes + 1))
    colnames(ancestryMat)[1:NTips] <- rownames(ancestryMat)[1:NTips] <- tipNames
    ancestryMat <- ancestryMat[-(NTips + 1), -(NTips + 1)]
    return(ancestryMat)
}

getTreeDetails <- function(tree) {
    edgeOrder <- order(tree$edge[,2])
    edgeLengths <- tree$edge.length[edgeOrder]
    NHs <- nodeHeights(tree)
    maxNHs <- max(NHs)
    NHRel <- NHs / maxNHs
    rownames(NHRel) <- tree$edge[,2]
    NHRel <- NHRel[edgeOrder,]
    logitNH <- logit(apply(NHRel[(length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode - 1),], 1, function(x) (x[[2]] - x[[1]]) / (1 - x[[1]])))
    pm <- makeParentMat(tree$Nnode - 1 + length(tree$tip.label), tree, length(tree$tip.label), tree$tip.label)[-1, -1]

    return(list(edgeOrder   = edgeOrder,
                NHs         = NHs,
                maxNHs      = maxNHs,
                NHRel       = NHRel,
                edgeLengths = edgeLengths,
                logitNH     = logitNH,
                pm          = pm))
}

makeIdentityMat <- function(NNodes, ...) {
    idMat <- diag(NNodes + 1)
    return(idMat)
}

makeParentMat <- function(NNodes, tree, NTips, tipNames, ...) {
    parentMat <- matrix(0, NNodes + 1, NNodes + 1)
    for(node in 1:(NNodes + 1)) {
        p <- Ancestors(tree, node, 'parent')
        parentMat[node, ] <- as.numeric(1:(NNodes + 1) %in% c(p, node))
    }
    colnames(parentMat) <- rownames(parentMat) <- paste0('i', 1:(NNodes + 1))
    colnames(parentMat)[1:NTips] <- rownames(parentMat)[1:NTips] <- paste0('t', tipNames)
    parentMat <- cbind(0, parentMat[-(NTips + 1), -(NTips + 1)])
    parentMat <- rbind(c(1, rep(0, ncol(parentMat) - 1)), parentMat)
    
    return(parentMat)
}

makeGrandparentMat <- function(NNodes, tree, NTips, tipNames, ...) {
    grandparentMat <- matrix(0, NNodes + 1, NNodes + 1)
    for(node in 1:(NNodes + 1)) {
        p <- Ancestors(tree, node, 'parent')
        grandparentMat[node, ] <- as.numeric(1:(NNodes + 1) %in% c(Ancestors(tree, p, 'parent'), p, node))
    }

    colnames(grandparentMat) <- rownames(grandparentMat) <- paste0('i', 1:(NNodes + 1))
    colnames(grandparentMat)[1:NTips] <- rownames(grandparentMat)[1:NTips] <- paste0('t', tipNames)
    grandparentMat <- cbind(0, grandparentMat[-(NTips + 1), -(NTips + 1)])
    grandparentMat <- rbind(c(1, rep(0, ncol(grandparentMat) - 1)), grandparentMat)
    
    return(grandparentMat)
}

makeRootMat <- function(NNodes, tree, NTips, tipNames, ...) {
    idMat <- matrix(0, nrow = NNodes + 1, ncol = NNodes + 1)
    idMat[1, 1] <- 1
    idMat[2:(NNodes + 1), 2:(NNodes + 1)] <- createAncestryMat(NNodes, tree, NTips, tipNames)
    return(idMat)
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

summarizeLcGLM <- function(combineTrees  = T,
                           separateTrees = T,
                           whichTrees    = NULL,
                           plotTrees     = T,
                           sumScales     = T,
                           sumVarMods    = T,
                           sumEffects    = T,
                           contrastFuns  = list(vsParent                     = list(host    = makeIdentityMat,
                                                                                    microbe = makeIdentityMat),
                                                vsGrandparent                = list(host    = makeParentMat,
                                                                                    microbe = makeParentMat),
                                                vsGreatgrandparent           = list(host    = makeGrandparentMat,
                                                                                    microbe = makeGrandparentMat),
                                                vsRoot                       = list(host    = makeRootMat,
                                                                                    microbe = makeRootMat),
                                                microbeVsRootHostVsParent    = list(host    = makeIdentityMat,
                                                                                    microbe = makeRootMat)),
                           ...) {
    
    fitModes <- sapply(fit, function(x) x@mode)
    
    colorpal <- colorRampPalette(brewer.pal(9, 'Blues'))
    plotcolors <- c('white', colorpal(100), 'black')
    
    baseNames <- sapply(sumconts, function(m) paste0(m, levels(newermap[,m])[nlevels(newermap[,m])]))

    ## summarize results for parameters that can be interpretted across all sampled host trees
    NSuccessTrees <- sum(fitModes == 0)

    ## summarize the results separately for each sampled host tree and/or for all trees combined
    for(i in 1:(separateTrees * NTrees + combineTrees * (NTrees > 1))) {
        
        if (NSuccessTrees > 1 & i == (separateTrees * NTrees + 1)) {
    
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
        
        } else if (!is.null(whichTrees) & !(i %in% whichTrees)) {
            
            next
            
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
        
        if(!variational) {
            sink(stdout(), type = "message")
            check_hmc_diagnostics(fit[[i]])
            sink(NULL, type = "message")
        }

        if(plotTrees) {
        
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
                                            pars     = 'microbeScales',
                                            permuted = F),
                                    dim = c(NMCSamples - warmup,
                                            NChains,
                                            NMicrobeNodes))^2,
                              3,
                              mean)
            finalMicrobeTree.newEdges <- finalMicrobeTree
            finalMicrobeTree.newEdges$tip.label <- substr(finalMicrobeTree.newEdges$tip.label, 1, 30)
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
                                            pars     = 'hostScales',
                                            permuted = F),
                                    dim = c(NMCSamples - warmup,
                                            NChains,
                                            NHostNodes))^2,
                              3,
                              mean)
            hostTreesSampled.newEdges <- hostTreesSampled[[i]]
            hostTreesSampled.newEdges$tip.label <- substr(hostTreesSampled.newEdges$tip.label, 1, 30)
            hostTreesSampled.newEdges$edge.length <- newEdges[order(hostTreeDetails[[i]]$edgeOrder)]
            pdf(file   = file.path(currplotdir, 'hostTreeWEstimatedEdgeLengths.pdf'),
                width  = 25,
                height = 15)
            plot(hostTreesSampled.newEdges,
                 cex = 0.5)
            graphics.off()
            ##
            
            ## plot heatmap of cophylogenetic patterns
            plotmicrobetree <- force.ultrametric(ladderize(multi2di(finalMicrobeTree)))
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
                
                temp2 <- temp[temp$tissue_compartment == comp & !is.na(temp$tissue_compartment),]
                hostvectTemp <- hostvect[rownames(temp2)]
                
                # expand the tips into polytomies containing a tip for each sample
                hosttree <- expandTaxonTree(hostTreesSampled[[i]],
                                            hostvectTemp,
                                            keepBrLen = T)
                hosttree <- drop.tip(hosttree,
                                     hosttree$tip.label[!hosttree$tip.label %in% names(hostvectTemp)])
                                     
                # convert polytomies into randomly-split, binary subtrees with 0-length branches,
                # and ladderize the whole tree
                hosttree.dichotomous <- force.ultrametric(ladderize(multi2di(hosttree, random = F), right = F))
                
                # convert the phylogenetic tree into an hclust object
                hclhosttree <- as.hclust(hosttree.dichotomous)
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
        
        }

        if(sumScales) {
            
            cat('\n\tRaw variance estimates\n\t')
            cat(paste0(as.character(Sys.time()), '\n'))
            
            ## variance partitioning (top level factors)
            stDProps <- array(extract(fit[[i]],
                                      pars       = 'stDProps',
                                      permuted   = F,
                                      inc_warmup = T),
                              dim = c(NMCSamples,
                                      NChains,
                                      2 * NFactors + 3),
                              dimnames = list(sample  = NULL,
                                              chain   = NULL,
                                              factor  = c(paste0('ADiv.', names(groupedFactors)),
                                                          paste0('Specificity.', names(groupedFactors)),
                                                          'ADiv.host',
                                                          'host.specificity',
                                                          'microbe.prevalence')))
            save(stDProps, file = file.path(currdatadir, 'stDProps.RData'))
            
            allRes <- monitor(stDProps,
                              warmup = warmup,
                              probs  = c(0.025, 0.5, 0.975),
                              print  = F)
            rownames(allRes) <- c(paste0('ADiv.', names(groupedFactors)),
                              paste0('Specificity.', names(groupedFactors)),
                              'ADiv.host',
                              'host.specificity',
                              'microbe.prevalence')
            cat('factor\t', file = file.path(currtabledir, 'stDProps.txt'))
            write.table(allRes,
                        file   = file.path(currtabledir, 'stDProps.txt'),
                        sep    = '\t',
                        quote  = F,
                        append = T)
            
            stDPropsPlot <- NULL
            for(j in 1:NChains) {
                stDPropsPlot <- rbind(stDPropsPlot, stDProps[(warmup + 1):NMCSamples, j,])
            }
            colnames(stDPropsPlot) <- gsub('_',
                                           ' ',
                                           c(paste0(names(groupedFactors), ' (ADiv)'),
                                             paste0(names(groupedFactors), ' (Specificity)'),
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
            
            ## variance partitioning (broken down by subfactors)
            NSubfactorGammas <- 0
            gammaNames <- NULL
            for(j in 1:NFactors) {
                if(NSubPerFactor[[j]] > 1) {
                    NSubfactorGammas <- NSubfactorGammas + NSubPerFactor[[j]]
                    gammaNames <- c(gammaNames, groupedFactors[[j]])
                }
            }
            subfactProps <- array(extract(fit[[i]],
                                          pars       = 'subfactProps',
                                          permuted   = F,
                                          inc_warmup = T),
                                  dim = c(NMCSamples,
                                          NChains,
                                          2 * sum(NSubPerFactor) + 3),
                                  dimnames = list(sample  = NULL,
                                                  chain   = NULL,
                                                  factor  = c(paste0('ADiv.', unlist(groupedFactors)),
                                                              paste0('Specificity.', unlist(groupedFactors)),
                                                              'ADiv.host',
                                                              'host.specificity',
                                                              'microbe.prevalence')))
            start <- 1
            for (j in 1:NFactors) {
                if(NSubPerFactor[[j]] > 1) {
                    
                    subfactPropsPlot <- NULL
                    for(k in 1:NChains) {
                        subfactPropsADivPlot <- rbind(subfactPropsPlot, subfactProps[(warmup + 1):NMCSamples,
                                                                                     k,
                                                                                     start:(start - 1 + NSubPerFactor[[j]])])
                        subfactPropsSpecPlot <- rbind(subfactPropsPlot, subfactProps[(warmup + 1):NMCSamples,
                                                                                     k,
                                                                                     (sum(NSubPerFactor) + start):(sum(NSubPerFactor) + start - 1 + NSubPerFactor[[j]])])
                    }
                    subfactPropsADivPlot <- t(apply(subfactPropsADivPlot, 1, function(x) x / sum(x)))
                    subfactPropsSpecPlot <- t(apply(subfactPropsSpecPlot, 1, function(x) x / sum(x)))
                    save(subfactPropsADivPlot, file = file.path(currdatadir, paste0('subfactProps_ADiv_', names(groupedFactors)[[j]], '.RData')))
                    save(subfactPropsSpecPlot, file = file.path(currdatadir, paste0('subfactProps_Spec_', names(groupedFactors)[[j]], '.RData')))
                    
                    allres <- monitor(array(subfactPropsADivPlot,
                                            dim = c(nrow(subfactPropsADivPlot),
                                                    1,
                                                    ncol(subfactPropsADivPlot))),
                                      warmup = 0,
                                      probs  = c(0.025, 0.5, 0.975),
                                      print  = F)
                    rownames(allres) <- groupedFactors[[j]]
                    cat('subfactor\t', file = file.path(currtabledir, paste0('subfactProps_ADiv_', names(groupedFactors)[[j]], '.txt')))
                    write.table(allres,
                                file   = file.path(currtabledir, paste0('subfactProps_ADiv_', names(groupedFactors)[[j]], '.txt')),
                                sep    = '\t',
                                quote  = F,
                                append = T)
                                
                    allres <- monitor(array(subfactPropsSpecPlot,
                                            dim = c(nrow(subfactPropsSpecPlot),
                                                    1,
                                                    ncol(subfactPropsSpecPlot))),
                                      warmup = 0,
                                      probs  = c(0.025, 0.5, 0.975),
                                      print  = F)
                    rownames(allres) <- groupedFactors[[j]]
                    cat('subfactor\t', file = file.path(currtabledir, paste0('subfactProps_Spec_', names(groupedFactors)[[j]], '.txt')))
                    write.table(allres,
                                file   = file.path(currtabledir, paste0('subfactProps_Spec_', names(groupedFactors)[[j]], '.txt')),
                                sep    = '\t',
                                quote  = F,
                                append = T)

                    pdf(file   = file.path(currplotdir, paste0('subfactProps_ADiv_', names(groupedFactors)[[j]], '_boxes.pdf')),
                        width  = 7,
                        height = 7)
                    boxplot(subfactPropsADivPlot,
                            cex.axis = 0.5,
                            las      = 2,
                            xlab     = 'Subfactor',
                            ylab     = paste0('Percent of ', names(groupedFactors)[[j]], ' (Adiv) variance'))
                    graphics.off()
                    
                    pdf(file   = file.path(currplotdir, paste0('subfactProps_Spec_', names(groupedFactors)[[j]], '_boxes.pdf')),
                        width  = 7,
                        height = 7)
                    boxplot(subfactPropsSpecPlot,
                            cex.axis = 0.5,
                            las      = 2,
                            xlab     = 'Subfactor',
                            ylab     = paste0('Percent of ', names(groupedFactors)[[j]], ' (Specificity) variance'))
                    graphics.off()
                    
                }
                start <- start + NSubPerFactor[[j]]
            }
            
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
        
        }
        
        if(sumVarMods) {
            
            if(!exists('metaScales')) {
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
            }
            
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
            
            dir.create(file.path(currtabledir, 'phyloVarianceEffects'), recursive = T)
                            
            for(contrast in names(contrastFuns)) {
                
                cat(paste0('\n\tVariance modifiers (', contrast, ')\n\t'))
                cat(paste0(as.character(Sys.time()), '\n'))
                
                ##
                hostMat <- contrastFuns[[contrast]][['host']](NHostNodes,
                                                              hostTreesSampled[[i]],
                                                              NHostTips,
                                                              hostTreesSampled[[i]]$tip.label)
                microbeMat <- t(contrastFuns[[contrast]][['microbe']](NMicrobeNodes,
                                                                      finalMicrobeTree,
                                                                      NMicrobeTips,
                                                                      microbeTips))
                ##
                
                ## sum the effects
                matMult <- array(NA,
                                 dim = c(NMCSamples,
                                         NChains,
                                         NHostNodes + 1,
                                         NMicrobeNodes + 1),
                                 dimnames = list(sample      = NULL,
                                                 chain       = NULL,
                                                 hostnode    = c('microbePrevalence',
                                                                 colnames(hostAncestors[[i]])),
                                                 microbenode = c('alphaDiversity',
                                                                 colnames(microbeAncestors))))
                                     
                for(j in 1:NMCSamples) {
                    for(k in 1:NChains) {
                        matMult[j,k,,] <- hostMat %*%
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
                
                cat('microbeNode\t', file = file.path(currtabledir, 'phyloVarianceEffects', paste0(contrast, '.txt')))
                write.table(allRes,
                            file   = file.path(currtabledir, 'phyloVarianceEffects', paste0(contrast, '.txt')),
                            sep    = '\t',
                            quote  = F,
                            append = T)
                ##
                
            }
        }
        
        if(sumEffects) {
            
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
            
            ##
            effectNames <- c('microbePrevalence',
                             colnames(modelMat)[2:(NEffects + 1)],
                             paste0('host_', colnames(hostAncestors[[i]])),
                             baseNames)
                             
            dir.create(file.path(currtabledir, 'nodeEffects'), recursive = T)
            
            for(contrast in names(contrastFuns)) {
                
                cat(paste0('\n\tNode effects (', contrast, ')\n\t'))
                cat(paste0(as.character(Sys.time()), '\n'))
                
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
                
                ##build a temporary model matrix
                hostMat <- diag(NHostNodes + NEffects + NSumTo0 + 1)
                hostMat[(NEffects + 2):(NEffects + NHostNodes + 1),
                        (NEffects + 2):(NEffects + NHostNodes + 1)] <- contrastFuns[[contrast]][['host']](NHostNodes,
                                                                                                          hostTreesSampled[[i]],
                                                                                                          NHostTips,
                                                                                                          hostTreesSampled[[i]]$tip.label)[-1, -1]
                microbeMat <- t(contrastFuns[[contrast]][['microbe']](NMicrobeNodes,
                                                                      finalMicrobeTree,
                                                                      NMicrobeTips,
                                                                      microbeTips))
                ##
                
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
                
                cat('microbeNode\t', file = file.path(currtabledir, 'nodeEffects', paste0(contrast, '.txt')))
                write.table(allRes,
                            file   = file.path(currtabledir, 'nodeEffects', paste0(contrast, '.txt')),
                            sep    = '\t',
                            quote  = F,
                            append = T)
                ##
            }
        }
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

runStanModel <- function(noData = F, shuffleData = F, shuffleSamples = F, variational = F, ...) {
    
    if(variational) {
        NMCSamples <- 1000
        warmup <- 0
    }
    
    if(sum(noData, shuffleData, shuffleSamples) > 1) {
        cat('\at most one of noData, shuffleData, and shuffleSamples can be TRUE.\n')
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
        presentShuffled <- sample(present)
        for (i in 1:NTrees) {
            standat[[i]]$present <- presentShuffled
        }
        
        subdir <- file.path(outdir, 'resampled_with_randomObs')
        dir.create(subdir, recursive = T)
        save.image(file.path(subdir, 'setup.RData'))

        cat('\nSampling from shuffled data\n')
    } else if(shuffleSamples) {
        
        ## shuffle data from stan list
        senddatShuffled <- melt(y[sample(1:nrow(y)),], varnames = c('sample', 'tip'), value.name = 'present')
        presentSamplesShuffled <- senddatShuffled[,3]
        for (i in 1:NTrees) {
            standat[[i]]$present <- presentSamplesShuffled
        }
        
        subdir <- file.path(outdir, 'resampled_with_randomSamples')
        dir.create(subdir, recursive = T)
        save.image(file.path(subdir, 'setup.RData'))
        
        cat('\nSampling from shuffled samples\n')
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
                if(!variational) {
                stan(file     = modelPath,
                     data     = standat[[i]],
                     control  = list(adapt_delta   = adapt_delta,
                                     max_treedepth = max_treedepth),
                     iter     = NIterations,
                     thin     = thin,
                     chains   = NChains,
                     seed     = seed,
                     chain_id = (NChains * (i - 1) + (1:NChains)),
                     pars     = c('rawMicrobeNodeEffects', 'sampleTipEffects'),
                     include  = FALSE,
                     init_r   = init_r)
                } else {
                    vb(stan_model(file = modelPath),
                       data     = standat[[i]],
                       iter     = 25000,
                       seed     = seed,
                       pars     = c('rawMicrobeNodeEffects', 'sampleTipEffects'),
                       include  = FALSE,
                       init_r   = init_r,
                       sample_file = file.path(subdir, 'samples.csv'))
                }
            }, error = function(e) {
                   print(e)
                   return(NA)
               })
        }, mc.preschedule = F,
           mc.cores       = NCores)

    cat('\nSaving results\n')
    cat(paste0(as.character(Sys.time()), '\n'))

    save(fit, file = file.path(subdir, 'fit.RData'))
    ##

    ## summarize the model fit
    environment(summarizeLcGLM) <- environment()
    summarizeLcGLM()
}
## fin
