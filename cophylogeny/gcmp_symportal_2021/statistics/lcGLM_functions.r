logit <- function(x) { log(x / (1 - x)) }
inv_logit <- function(x) { 1 / (1 + exp(-x)) }

getTreeDetails <- function(tree) {
    edgeOrder <- order(tree$edge[,2])
    edgeLengths <- tree$edge.length[edgeOrder]
    NHs <- nodeHeights(tree)
    maxNHs <- max(NHs)
    NHRel <- NHs / maxNHs
    rownames(NHRel) <- tree$edge[,2]
    NHRel <- NHRel[edgeOrder,]
    logitNH <- logit(apply(NHRel[(length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode - 1),],
                           1,
                           function(x) (x[[2]] - x[[1]]) / (1 - x[[1]])))
    pm <- makeContrastMat(tree$Nnode - 1 + length(tree$tip.label),
                          tree,
                          length(tree$tip.label),
                          tree$tip.label,
                          1)[-1, -1]

    return(list(edgeOrder   = edgeOrder,
                NHs         = NHs,
                maxNHs      = maxNHs,
                NHRel       = NHRel,
                edgeLengths = edgeLengths,
                logitNH     = logitNH,
                pm          = pm))
}

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

makeContrastMat <- function(NNodes, tree = NULL, NTips = NULL, tipNames = NULL, nLevels = 0, ...) {
    ## creates a matrix that, when multiplied by a matrix of raw per-branch effects, sums the effects of parental branches to some specified degree in order to draw contrasts between different parts of the tree
    if (nLevels == 0) {
        contrastMat <- diag(NNodes + 1)
    } else {
        contrastMat <- matrix(0, nrow = NNodes + 1, ncol = NNodes + 1)
        if(nLevels == 'vsRoot') {
            contrastMat[1, 1] <- 1
            contrastMat[2:(NNodes + 1), 2:(NNodes + 1)] <- createAncestryMat(NNodes, tree, NTips, tipNames)
        } else {
            for(node in 1:(NNodes + 1)) {
                x <- node
                for(i in 1:nLevels) {
                    x <- c(Ancestors(tree, x[[1]], 'parent'), x)
                }
                contrastMat[node, ] <- as.numeric(1:(NNodes + 1) %in% x)
            }
            colnames(contrastMat) <- rownames(contrastMat) <- paste0('i', 1:(NNodes + 1))
            colnames(contrastMat)[1:NTips] <- rownames(contrastMat)[1:NTips] <- paste0('t', tipNames)
            contrastMat <- cbind(0, contrastMat[-(NTips + 1), -(NTips + 1)])
            contrastMat <- rbind(c(1, rep(0, ncol(contrastMat) - 1)), contrastMat)
        }
    }
    return(contrastMat)
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

bridge_sampler_diag <- function(sf, df, warmup=0, repetitions = 1, method = "normal", cores = 1, use_neff = TRUE, maxiter = 1000, silent = FALSE, verbose = FALSE, ...) {
    #concepts largely copied from bridge_sampler.stanreg
    samples <- coda::as.mcmc.list(lapply(df, FUN = function(f) {
        cat('Loading data\n')
        nupars <- rstan::get_num_upars(sf)
        #skip first cols ("lp__", "accept_stat__", "stepsize__", "treedepth__", "n_leapfrog__", "divergent__", "energy__") and all cols after the number of unconstrained params
        coda::as.mcmc(as.matrix(data.table::fread(sep=',', cmd=paste0('grep -v "^#" "', f, '" | tail -n +', warmup+1), select = 8:(nupars+7))))
    }))
    lb <- rep(-Inf, ncol(samples[[1]]))
    ub <- rep(Inf, ncol(samples[[1]]))
    names(lb) <- names(ub) <- colnames(samples[[1]])
    cat('Bridge sampling\n')
    bridge_output <- bridge_sampler(samples = samples, log_posterior = bridgesampling:::.stan_log_posterior,
                                    data = list(stanfit = sf), lb = lb, ub = ub, repetitions = repetitions,
                                    method = method, cores = cores, use_neff = use_neff,
                                    packages = "rstan", maxiter = maxiter, silent = silent,
                                    verbose = verbose)
    return(bridge_output)
}

summarizeLcGLM <- function(combineTrees    = T,
                           separateTrees   = T,
                           whichTrees      = NULL,
                           plotTrees       = T,
                           sumScales       = T,
                           sumVarMods      = T,
                           sumEffects      = T,
                           tipNamesAreSeqs = F,
                           contrastLevels  = list(vsParent                     = list(host    = 0,
                                                                                      microbe = 0),
                                                  vsGrandparent                = list(host    = 1,
                                                                                      microbe = 1),
                                                  vsGreatgrandparent           = list(host    = 2,
                                                                                      microbe = 2),
                                                  vsGreatgreatgrandparent      = list(host    = 3,
                                                                                      microbe = 3),
                                                  vsRoot                       = list(host    = 'vsRoot',
                                                                                      microbe = 'vsRoot'),
                                                  microbeVsRootHostVsParent    = list(host    = 0,
                                                                                      microbe = 'vsRoot'),
                                                  hostVsRootMicrobeVsParent    = list(host    = 'vsRoot',
                                                                                      microbe = 0)),
                           ...) {

    fitModes <- sapply(fit, function(x) x@mode)

    colorpal <- colorRampPalette(brewer.pal(9, 'Blues'))
    plotcolors <- c('white', colorpal(100), 'black')

    mycols <- rev(brewer.pal(11, 'RdYlBu'))
    mycols[[6]] < 'white'
    colorpal <- colorRampPalette(mycols)
    plotcolorsVar <- c(colorpal(101))

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

            csvfiles <- file.path(subdir, paste0('samples_chain', (1:NTrees)[fitModes == 0], '.csv'))

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

            csvfiles <- file.path(subdir, paste0('samples_chain', i, '.csv'))

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
                                            pars     = 'microbeNewEdges',
                                            permuted = F),
                                    dim = c(NMCSamples - warmup,
                                            NChains,
                                            NMicrobeNodes)),
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
                                            pars     = 'hostNewEdges',
                                            permuted = F),
                                    dim = c(NMCSamples - warmup,
                                            NChains,
                                            NHostNodes)),
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
            plotmicrobetree <- ladderize(multi2di(force.ultrametric(finalMicrobeTree), random = F))
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
                tryCatch({

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
                    hosttree.dichotomous <- ladderize(multi2di(force.ultrametric(hosttree), random = F), right = F)

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

                }, error = function(e) print(e))
            }
            ##

        }

        if(sumScales) {
            tryCatch({

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
                                              NFactors + 3),
                                      dimnames = list(sample  = NULL,
                                                      chain   = NULL,
                                                      effect  = c(paste0('Specificity.',
                                                                         names(groupedFactors)),
                                                                  'Prevalence',
                                                                  'ADiv',
                                                                  'host.specificty')))
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
                                            NSubfactors + 3),
                                    dimnames = list(sample  = NULL,
                                                    chain   = NULL,
                                                    effect  = c(paste0('Specificity.',
                                                                       unlist(groupedFactors)),
                                                                'Prevalence',
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

            }, error = function(e) print(e))
        }

        if(sumVarMods) {

            cat('\n\tLoading raw meta-variance\n\t')
            cat(paste0(as.character(Sys.time()), '\n'))

            if(!exists('metaScales')) {
                metaScales <- array(extract(fit[[i]],
                            pars       = 'metaScales',
                            permuted   = F,
                            inc_warmup = T),
                    dim = c(NMCSamples,
                            NChains,
                            NSubfactors + 3),
                    dimnames = list(sample  = NULL,
                                    chain   = NULL,
                                    effect  = c(paste0('Specificity.',
                                                       unlist(groupedFactors)),
                                                'Prevalence',
                                                'ADiv',
                                                'Specificty')))
            }

            ## extract variance modifier terms from the fit model
            csvBuffer <- read_stan_csv_subset(csvfiles,
                                              params = c('phyloLogVarMultPrev','phyloLogVarMultADiv','phyloLogVarMultRaw','phyloLogVarMultFacts'))

            phyloLogVarMultPrev <- array(extract(csvBuffer,
                                                 pars       = 'phyloLogVarMultPrev',
                                                 permuted   = F,
                                                 inc_warmup = T),
                                         dim = c(NMCSamples,
                                                 NChains,
                                                 NMicrobeNodes),
                                         dimnames = list(sample  = NULL,
                                                         chain   = NULL,
                                                         taxnode = colnames(microbeAncestors)))

            phyloLogVarMultADiv <- array(extract(csvBuffer,
                                                 pars       = 'phyloLogVarMultADiv',
                                                 permuted   = F,
                                                 inc_warmup = T),
                                         dim = c(NMCSamples,
                                                 NChains,
                                                 NHostNodes),
                                         dimnames = list(sample  = NULL,
                                                         chain   = NULL,
                                                         taxnode = colnames(hostAncestors[[i]])))

            phyloLogVarMultRaw <- array(extract(csvBuffer,
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

            phyloLogVarMultFacts <- array(extract(csvBuffer,
                                                  pars       = 'phyloLogVarMultFacts',
                                                  permuted   = F,
                                                  inc_warmup = T),
                                          dim = c(NMCSamples,
                                                  NChains,
                                                  NSubfactors,
                                                  NMicrobeNodes),
                                          dimnames = list(sample      = NULL,
                                                          chain       = NULL,
                                                          factor      = unlist(groupedFactors),
                                                          microbenode = colnames(microbeAncestors)))

            rm('csvBuffer')
            gc()
            ##

            phyloLogVarMultScaled <- array(NA,
                                           dim = c(NMCSamples,
                                                   NChains,
                                                   NSubfactors + NHostNodes + 1,
                                                   NMicrobeNodes + 1))

            for(j in 1:NMCSamples) {
                for(k in 1:NChains) {
                    outerFactors <- outer(metaScales[j,k, 1:NSubfactors], sqrt(microbeTreeDetails$edgeLengths))
                    sqrtMicrobeEdges <- sqrt(microbeTreeDetails$edgeLengths) * metaScales[j,k, NSubfactors + 1]
                    sqrtHostEdges <- sqrt(hostTreeDetails[[i]]$edgeLengths) * metaScales[j,k, NSubfactors + 2]
                    outerEdges <- sqrt(outer(hostTreeDetails[[i]]$edgeLengths, microbeTreeDetails$edgeLengths)) * metaScales[j,k, NSubfactors + 3]

                    phyloLogVarMultScaled[j,k, 1:(NSubfactors + 1), 1] <- 0
                    phyloLogVarMultScaled[j,k, (NSubfactors + 2):(NSubfactors + NHostNodes + 1), 1] <- sqrtHostEdges * phyloLogVarMultADiv[j,k,]
                    phyloLogVarMultScaled[j,k, 1, 2:(NMicrobeNodes + 1)] <- sqrtMicrobeEdges * phyloLogVarMultPrev[j,k,]
                    phyloLogVarMultScaled[j,k, 2:(NSubfactors + 1), 2:(NMicrobeNodes + 1)] <- outerFactors * phyloLogVarMultFacts[j,k,,]
                    phyloLogVarMultScaled[j,k, (NSubfactors + 2):(NSubfactors + NHostNodes + 1), 2:(NMicrobeNodes + 1)] <- outerEdges * phyloLogVarMultRaw[j,k,,]
                }
            }

            rm('phyloLogVarMultPrev', 'phyloLogVarMultADiv', 'phyloLogVarMultRaw', 'phyloLogVarMultFacts')
            gc()

            ##

            plotmicrobetree <- ladderize(multi2di(force.ultrametric(finalMicrobeTree), random = F))
            if(tipNamesAreSeqs) {
                plotmicrobetree2 <- plotmicrobetree
                newMicrobeNames <- plotmicrobetree$tip.label
                names(newMicrobeNames) <- paste0('t', 1:length(plotmicrobetree$tip.label))
                plotmicrobetree2$tip.label <- names(newMicrobeNames)
                write.dna(newMicrobeNames,
                          file   = file.path(currdatadir, 'seqnames.fasta'),
                          format = 'fasta',
                          nbcol  = -1,
                          colw   = 10000)
                hclmicrobetree <- as.hclust(plotmicrobetree2)
            } else {
                hclmicrobetree <- as.hclust(plotmicrobetree)
            }
            plothosttree <- ladderize(multi2di(force.ultrametric(hostTreesSampled[[i]]), random = F), right = F)
            hclhosttree <- as.hclust(plothosttree)

            dir.create(file.path(currtabledir, 'phyloVarianceEffects'), recursive = T)

            for(contrast in names(contrastLevels)) {
                tryCatch({

                    cat(paste0('\n\tVariance modifiers (', contrast, ')\n\t'))
                    cat(paste0(as.character(Sys.time()), '\n'))

                    ##
                    hostMat <- diag(NHostNodes + NSubfactors + 1)
                    hostMat[(NSubfactors + 2):(NSubfactors + NHostNodes + 1),
                            (NSubfactors + 2):(NSubfactors + NHostNodes + 1)] <- makeContrastMat(NHostNodes,
                                                                                                 hostTreesSampled[[i]],
                                                                                                 NHostTips,
                                                                                                 hostTreesSampled[[i]]$tip.label,
                                                                                                 contrastLevels[[contrast]][['host']])[-1, -1]
                    hostMat <- Matrix(hostMat)

                    microbeMat <- Matrix(t(makeContrastMat(NMicrobeNodes,
                                                           finalMicrobeTree,
                                                           NMicrobeTips,
                                                           microbeTips,
                                                           contrastLevels[[contrast]][['microbe']])))
                    ##

                    ## sum the effects
                    matMult <- array(NA,
                                     dim = c(NMCSamples,
                                             NChains,
                                             NSubfactors + NHostNodes + 1,
                                             NMicrobeNodes + 1),
                                     dimnames = list(sample      = NULL,
                                                     chain       = NULL,
                                                     hostnode    = c('microbePrevalence',
                                                                     unlist(groupedFactors),
                                                                     colnames(hostAncestors[[i]])),
                                                     microbenode = c('alphaDiversity',
                                                                     colnames(microbeAncestors))))

                    for(j in 1:NMCSamples) {
                        for(k in 1:NChains) {
                            matMult[j,k,,] <- as.matrix(hostMat %*%
                                                        phyloLogVarMultScaled[j,k,,] %*%
                                                        microbeMat)
                        }
                    }

                    allRes <- NULL
                    for(j in 1:(NSubfactors + NHostNodes + 1)) {
                        temp <- monitor(array(matMult[,,j,],
                                              dim = c(NMCSamples,
                                                      NChains,
                                                      NMicrobeNodes + 1)),
                                        warmup = warmup,
                                        probs  = c(0.05, 0.95),
                                        print  = F)
                        temp <- cbind(hostNode = c('microbePrevalence',
                                                   unlist(groupedFactors),
                                                   paste0('host_', colnames(hostAncestors[[i]])))[[j]], temp)
                        temp <- cbind(microbeNode = c('alphaDiversity',
                                                      rownames(microbeAncestors)),
                                      temp)
                        allRes <- rbind(allRes, temp)
                        statusUpdate(j, NHostNodes)
                    }

                    write.table(allRes,
                                file   = file.path(currtabledir, 'phyloVarianceEffects', paste0(contrast, '.txt')),
                                sep    = '\t',
                                quote  = F,
                                row.names = F,
                                append = T)
                    ##

                    allRes <- matrix(NA, nrow = NMicrobeNodes + 1, ncol = NSubfactors + NHostNodes + 1)
                    for(j in 1:(NSubfactors + NHostNodes + 1)) {
                        for(k in 1:(NMicrobeNodes + 1)) {
                            allRes[k,j] <- median(matMult[,,j,k])
                        }
                    }
                    rownames(allRes) <- c('alphaDiversity',
                                          substr(colnames(microbeAncestors),
                                                 2,
                                                 nchar(colnames(microbeAncestors))))
                    colnames(allRes) <- c('microbePrevalence', unlist(groupedFactors), colnames(hostAncestors[[i]]))
                    ##

                    save(allRes, file = file.path(currdatadir, paste0('variance_heatmap_', contrast, '.RData')))

                    plotFilt <- as.matrix(allRes[plotmicrobetree$tip.label, plothosttree$tip.label])
                    if(tipNamesAreSeqs) {
                        rownames(plotFilt) <- names(newMicrobeNames)
                    }

                    pdf(file   = file.path(currplotdir,
                                           paste0('variance_heatmap_',
                                                  contrast,
                                                  '.pdf')),
                        width  = 10,
                        height = 10)
                    heatmap(plotFilt,
                            Rowv   = as.dendrogram(hclmicrobetree),
                            Colv   = as.dendrogram(hclhosttree),
                            col    = plotcolorsVar,
                            cexCol = 0.2,
                            cexRow = 0.1,
                            scale  = 'none')
                    graphics.off()

                }, error = function(e) print(e))
            }

            rm('phyloLogVarMultScaled', 'matMult', 'allRes')
            gc()
        }

        if(sumEffects) {

            ## extract effects from model
            scaledMicrobeNodeEffects <- array(extract(read_stan_csv_subset(csvfiles,
                                                                           params = 'scaledMicrobeNodeEffects'),
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
            baseLevelEffects <- array(extract(read_stan_csv_subset(csvfiles,
                                                                   params = 'baseLevelEffects'),
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

            for(contrast in names(contrastLevels)) {
                tryCatch({

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
                            (NEffects + 2):(NEffects + NHostNodes + 1)] <- makeContrastMat(NHostNodes,
                                                                                           hostTreesSampled[[i]],
                                                                                           NHostTips,
                                                                                           hostTreesSampled[[i]]$tip.label,
                                                                                           contrastLevels[[contrast]][['host']])[-1, -1]
                    hostMat <- Matrix(hostMat)
                    microbeMat <- Matrix(t(makeContrastMat(NMicrobeNodes,
                                                           finalMicrobeTree,
                                                           NMicrobeTips,
                                                           microbeTips,
                                                           contrastLevels[[contrast]][['microbe']])))
                    ##

                    for(j in 1:NMCSamples) {
                        for(k in 1:NChains) {
                            matMult[j,k,,] <- as.matrix(hostMat %*%
                                                        rbind(scaledMicrobeNodeEffects[j,k,,],
                                                              baseLevelEffects[j,k,,]) %*%
                                                        microbeMat)
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
                        temp <- cbind(microbeNode = c('alphaDiversity', rownames(microbeAncestors)), temp)
                        allRes <- rbind(allRes, temp)
                        statusUpdate(j, NHostNodes + NEffects + NSumTo0 + 1)
                    }

                    write.table(allRes,
                                file      = file.path(currtabledir, 'nodeEffects', paste0(contrast, '.txt')),
                                sep       = '\t',
                                quote     = F,
                                row.names = F,
                                append    = T)
                    ##
                }, error = function(e) print(e))
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

    if(exists('fit')) {
        rm(fit)
        gc()
    }

    if(variational) {
        NIterations <- NMCSamples <- round(1000 / NTrees)
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

        cat('\nSampling from prior\n')
    } else if(shuffleData) {

        ## shuffle data from stan list
        presentShuffled <- sample(present)
        for (i in 1:NTrees) {
            standat[[i]]$present <- presentShuffled
        }

        subdir <- file.path(outdir, 'resampled_with_randomObs')

        cat('\nSampling from shuffled data\n')
    } else if(shuffleSamples) {

        ## shuffle data from stan list
        senddatShuffled <- melt(y[sample(1:nrow(y)),], varnames = c('sample', 'tip'), value.name = 'present')
        presentSamplesShuffled <- senddatShuffled[,3]
        for (i in 1:NTrees) {
            standat[[i]]$present <- presentSamplesShuffled
        }

        subdir <- file.path(outdir, 'resampled_with_randomSamples')

        cat('\nSampling from shuffled samples\n')
    } else {

        subdir <- file.path(outdir, 'primary_sampling')

        cat('\nSampling from model\n')
    }

    dir.create(subdir, recursive = T)
    save(list = c(ls(envir = sys.frame(0)),
                  ls(envir = sys.frame(1))),
         file = file.path(subdir, 'setup.RData'))

    cat(paste0(as.character(Sys.time()), '\n'))

    ## run the model!
    fit <- mclapply(1:NTrees,
        function(i) {
            setTimeLimit(timeLimit)
            tryCatch({
                if(!variational) {
                sampling(object          = sm,
                         data            = standat[[i]],
                         control         = list(adapt_delta   = adapt_delta,
                                                adapt_t0      = adapt_t0,
                                                adapt_kappa   = adapt_kappa,
                                                adapt_gamma   = adapt_gamma,
                                                max_treedepth = max_treedepth),
                         iter            = NIterations,
                         thin            = thin,
                         chains          = NChains,
                         seed            = seed,
                         chain_id        = (NChains * (i - 1) + (1:NChains)),
                         pars            = c('aveStD', 'stDProps', 'aveStDMeta', 'metaScales', 'metaVarProps', 'subfactProps', 'subfactMetaProps', 'microbeNewEdges', 'hostNewEdges'),
                         include         = TRUE,
                         init_r          = init_r,
                         sample_file     = file.path(subdir, paste0('samples_chain', i, '.csv')),
                         diagnostic_file = file.path(subdir, paste0('diagnostics_chain', i, '.csv')))
                } else {
                    Sys.sleep((i - 1) * 3 * 60)
                    vb(stan_model(file = modelPath),
                       data           = standat[[i]],
                       iter           = 25000,
                       seed           = seed,
                       pars           = c('rawMicrobeNodeEffects'),
                       include        = FALSE,
                       init_r         = init_r,
                       output_samples = NIterations,
                       sample_file    = file.path(subdir, paste0('samples_chain', i, '.csv')))
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
    tryCatch(summarizeLcGLM(), error = function(e) print(e))
}

read_stan_csv_subset <- function (csvfiles, col_major = TRUE, params = NULL, keep = TRUE) {
    ##mostly copied from rstan's read_stan_csv, but with additional code to enable the loading of just a subset of parameters
    if (length(csvfiles) < 1)
        stop("csvfiles does not contain any CSV file name")
    g_skip <- 10
    ss_lst <- vector("list", length(csvfiles))
    cs_lst2 <- vector("list", length(csvfiles))
    for (i in seq_along(csvfiles)) {
        header <- rstan:::read_csv_header(csvfiles[i])
        lineno <- attr(header, "lineno")
        vnames <- strsplit(header, ",")[[1]]
        ##
        if(!is.null(params)) {
            oldnames <- vnames
            alwayskeep <- 'lp__|accept_stat__|stepsize__|treedepth__|n_leapfrog__|divergent__|energy__'
            if(keep)
                vnames <- grep(paste(c(params,alwayskeep),collapse='|'), vnames, value=T)
            else
                vnames <- grep(paste(c(params,alwayskeep),collapse='|'), vnames, value=T, invert=T)
            inds <- which(oldnames %in% vnames)
        }
        ##
        iter.count <- attr(header, "iter.count")
        variable.count <- length(vnames)
        df <- structure(replicate(variable.count, list(numeric(iter.count))),
            names = vnames, row.names = c(NA, -iter.count), class = "data.frame")
        comments = character()
        con <- file(csvfiles[[i]], "rb")
        buffer.size <- min(ceiling(1e+06/variable.count), iter.count)
        row.buffer <- matrix(ncol = variable.count, nrow = buffer.size)
        row <- 1
        buffer.pointer <- 1
        while (length(char <- readBin(con, "int", size = 1L)) > 0) {
            seek(con, origin = "current", -1)
            if (char == 35) {
                line <- readLines(con, n = 1)
                comments <- c(comments, line)
                next
            }
            if (char == 108) {
                readLines(con, n = 1)
                next
            }
            if (char == 10) {
                readLines(con, n = 1)
                next
            }
            row.buffer[buffer.pointer, ] <- scan(con, nlines = 1,
                sep = ",", quiet = TRUE)[inds]
            if (buffer.pointer == buffer.size) {
                df[row:(row + buffer.size - 1), ] <- row.buffer
                row <- row + buffer.size
                buffer.pointer <- 0
            }
            buffer.pointer <- buffer.pointer + 1
        }
        if (buffer.pointer > 1) {
            df[row:(row + buffer.pointer - 2), ] <- row.buffer[1:(buffer.pointer -
                1), ]
            df <- df[1:(row + buffer.pointer - 2),]
        } else {
            df <- df[1:(row - 1),]
        }
        close(con)
        cs_lst2[[i]] <- rstan:::parse_stancsv_comments(comments)
        if ("output_samples" %in% names(cs_lst2[[i]]))
            df <- df[-1, ]
        ss_lst[[i]] <- df
    }
    m_name <- sub("(_\\d+)*$", "", rstan:::filename_rm_ext(basename(csvfiles[1])))
    sdate <- do.call(max, lapply(csvfiles, function(csv) file.info(csv)$mtime))
    sdate <- format(sdate, "%a %b %d %X %Y")
    chains <- length(ss_lst)
    fnames <- names(ss_lst[[1]])
    n_save <- nrow(ss_lst[[1]])
    paridx <- rstan:::paridx_fun(fnames)
    lp__idx <- attr(paridx, "meta")["lp__"]
    par_fnames <- c(fnames[paridx], "lp__")
    pars_oi <- rstan:::unique_par(par_fnames)
    dims_oi <- lapply(pars_oi, function(i) {
        pat <- paste("^", i, "(\\.\\d+)*$", sep = "")
        i_fnames <- par_fnames[grepl(pat, par_fnames)]
        rstan:::get_dims_from_fnames(i_fnames, i)
    })
    names(dims_oi) <- pars_oi
    midx <- if (!col_major)
        rstan:::multi_idx_row2colm(dims_oi)
    else 1:length(par_fnames)
    if (chains > 1) {
        if (!all(sapply(ss_lst[-1], function(i) identical(names(i),
            fnames))))
            stop("the CSV files do not have same parameters")
        if (!all(sapply(ss_lst[-1], function(i) identical(length(i[[1]]),
            n_save))))
            stop("the number of iterations are not the same in all CSV files")
    }
    mode <- 0L
    samples <- lapply(ss_lst, function(df) {
        ss <- df[c(paridx, lp__idx)[midx]]
        attr(ss, "sampler_params") <- df[setdiff(attr(paridx,
            "meta"), lp__idx)]
        ss
    })
    par_fnames <- par_fnames[midx]
    for (i in seq_along(samples)) {
        attr(samples[[i]], "adaptation_info") <- cs_lst2[[i]]$adaptation_info
        attr(samples[[i]], "args") <- list(sampler_t = cs_lst2[[i]]$sampler_t,
            chain_id = cs_lst2[[i]]$chain_id)
        if (cs_lst2[[i]]$has_time)
            attr(samples[[i]], "elapsed_time") <- rstan:::get_time_from_csv(cs_lst2[[i]]$time_info)
    }
    save_warmup <- sapply(cs_lst2, function(i) i$save_warmup)
    warmup <- sapply(cs_lst2, function(i) i$warmup)
    thin <- sapply(cs_lst2, function(i) i$thin)
    iter <- sapply(cs_lst2, function(i) i$iter)
    if (!rstan:::all_int_eq(warmup) || !rstan:::all_int_eq(thin) || !rstan:::all_int_eq(iter))
        stop("not all iter/warmups/thin are the same in all CSV files")
    n_kept0 <- 1 + (iter - warmup - 1)%/%thin
    warmup2 <- 0
    if (max(save_warmup) == 0L) {
        n_kept <- n_save
    }
    else if (min(save_warmup) == 1L) {
        warmup2 <- 1 + (warmup[1] - 1)%/%thin[1]
        n_kept <- n_save - warmup2
    }
    if (n_kept0[1] != n_kept) {
        warning("the number of iterations after warmup found (",
            n_kept, ") does not match iter/warmup/thin from CSV comments (",
            paste(n_kept0, collapse = ","), ")")
        if (n_kept < 0) {
            warmup <- warmup + n_kept
            n_kept <- 0
            mode <- 2L
        }
        n_kept0 <- n_save
        iter <- n_save
        for (i in 1:length(cs_lst2)) {
            cs_lst2[[i]]$warmup <- warmup
            cs_lst2[[i]]$iter <- iter
        }
    }
    idx_kept <- if (warmup2 == 0)
        1:n_kept
    else -(1:warmup2)
    for (i in seq_along(samples)) {
        m <- vapply(samples[[i]], function(x) mean(x[idx_kept]),
            numeric(1))
        attr(samples[[i]], "mean_pars") <- m[-length(m)]
        attr(samples[[i]], "mean_lp__") <- m["lp__"]
    }
    perm_lst <- lapply(1:chains, function(id) sample.int(n_kept))
    sim = list(samples = samples, iter = iter[1], thin = thin[1],
        warmup = warmup[1], chains = chains, n_save = rep(n_save,
            chains), warmup2 = rep(warmup2, chains), permutation = perm_lst,
        pars_oi = pars_oi, dims_oi = dims_oi, fnames_oi = rstan:::dotfnames_to_sqrfnames(par_fnames),
        n_flatnames = length(par_fnames))
    null_dso <- new("cxxdso", sig = list(character(0)), dso_saved = FALSE,
        dso_filename = character(0), modulename = character(0),
        system = R.version$system, cxxflags = character(0), .CXXDSOMISC = new.env(parent = emptyenv()))
    null_sm <- new("stanmodel", model_name = m_name, model_code = character(0),
        model_cpp = list(), dso = null_dso)
    nfit <- new("stanfit", model_name = m_name, model_pars = pars_oi,
        par_dims = dims_oi, mode = mode, sim = sim, inits = list(),
        stan_args = cs_lst2, stanmodel = null_sm, date = sdate,
        .MISC = new.env(parent = emptyenv()))
    return(nfit)
}

## fin
