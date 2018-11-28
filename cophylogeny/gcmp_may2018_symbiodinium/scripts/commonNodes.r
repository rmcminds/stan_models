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
library(ggplot2)
library(RColorBrewer)
library(Matrix)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

outdir <- 'output/2018-11-25_13-52-47'
subdir <- file.path(outdir, 'primary_sampling')

load(file.path(subdir, 'setup.RData'))
load(file.path(subdir, 'fit.RData'))

source(file.path('scripts', 'lcGLM_functions.r'))

## consider the effects of nodes from different trees 'the same' if all of their descendants are the same.

currtabledir <- file.path(subdir, 'alltrees', 'tables')
currdatadir <- file.path(subdir, 'alltrees', 'data')
dir.create(file.path(currtabledir, 'nodeEffects'), recursive = T)
dir.create(file.path(currtabledir, 'phyloVarianceEffects'), recursive = T)
dir.create(currdatadir, recursive = T)

contrastLevels  = list(vsParent                = 0,
                       vsGrandparent           = 1,
                       vsGreatgrandparent      = 2,
                       vsGreatgreatgrandparent = 3)

contrastLevels  = list(vsRoot = 'vsRoot')

for(contrast in names(contrastLevels)) {

    ##build a temporary model matrix
    microbeMat <- Matrix(t(makeContrastMat(NMicrobeNodes,
                                           finalMicrobeTree,
                                           NMicrobeTips,
                                           microbeTips,
                                           contrastLevels[[contrast]])))
    ##

    hostNodes <- sort(unique(hostTreesSampled[[1]]$edge[,2]))

    if(contrastLevels[[contrast]] == 'vsRoot') {
        allMatches <- NULL
        for(tree1 in 1:(NTrees - 1)) {
            for(node1 in 1:length(hostNodes)) {
                for(tree2 in (tree1 + 1):NTrees) {
                    treeDiff <- setdiff(union(hostTreesSampled[[tree1]]$tip.label,
                                              hostTreesSampled[[tree2]]$tip.label),
                                        intersect(hostTreesSampled[[tree1]]$tip.label,
                                                  hostTreesSampled[[tree2]]$tip.label))
                    
                    isMatch <- sapply(hostNodes, function(node2) {
                        nodeUnion <- union(hostTreesSampled[[tree1]]$tip.label[Descendants(hostTreesSampled[[tree1]], hostNodes[[node1]])[[1]]],
                                           hostTreesSampled[[tree2]]$tip.label[Descendants(hostTreesSampled[[tree2]], node2)[[1]]])
                        nodeIntersect <- intersect(hostTreesSampled[[tree1]]$tip.label[Descendants(hostTreesSampled[[tree1]], hostNodes[[node1]])[[1]]],
                                                   hostTreesSampled[[tree2]]$tip.label[Descendants(hostTreesSampled[[tree2]], node2)[[1]]])
                        nodeDiffs <- setdiff(nodeUnion,
                                             nodeIntersect)
                                             
                        daughters1 <- Children(hostTreesSampled[[tree1]], hostNodes[[node1]])
                        if(length(daughters1) > 0) {
                            daughter1Miss <- sum(!sapply(daughters1, function(x) {
                                all(hostTreesSampled[[tree1]]$tip.label[Descendants(hostTreesSampled[[tree1]], x)[[1]]] %in% treeDiff)
                            })) > 1
                        } else {
                            daughter1Miss <- TRUE
                        }
                        
                        daughters2 <- Children(hostTreesSampled[[tree2]], node2)
                        if(length(daughters2) > 0) {
                            daughter2Miss <- sum(!sapply(daughters2, function(x) {
                                all(hostTreesSampled[[tree2]]$tip.label[Descendants(hostTreesSampled[[tree2]], x)[[1]]] %in% treeDiff)
                            })) > 1
                        } else {
                            daughter2Miss <- TRUE
                        }
                        
                        return(length(nodeIntersect) > 0 & all(nodeDiffs %in% treeDiff) & daughter1Miss & daughter2Miss)
                        # there is at least some overlap in the descendants of each node and some overlap in each node's sister clade; all the differences that do exist are due to species that are simply not present in one of the trees; and the allowable differences due to absent species does not include the case where the extra species form a basal clade under the node of interest
                    })
                    if(sum(isMatch) > 0) {
                        allMatches <- rbind(allMatches,
                                            cbind(paste0('tree', tree1),
                                                  colnames(hostAncestors[[tree1]])[[node1]],
                                                  paste0('tree', tree2),
                                                  colnames(hostAncestors[[tree2]])[isMatch]))
                    }
                }
            }
        }
    } else {
        allMatches <- NULL
        for(tree1 in 1:(NTrees - 1)) {
            for(node1 in 1:length(hostNodes)) {
                x <- hostNodes[[node1]]
                for(i in 1:(contrastLevels[[contrast]] + 1)) {
                    x <- c(Ancestors(hostTreesSampled[[tree1]], x[[1]], 'parent'), x)
                }
                containingNodes1 <- min(x)
                for(tree2 in (tree1 + 1):NTrees) {
                    treeDiff <- setdiff(union(hostTreesSampled[[tree1]]$tip.label,
                                              hostTreesSampled[[tree2]]$tip.label),
                                        intersect(hostTreesSampled[[tree1]]$tip.label,
                                                  hostTreesSampled[[tree2]]$tip.label))
                                                  
                    daughters1 <- Children(hostTreesSampled[[tree1]], hostNodes[[node1]])
                    if(length(daughters1) > 0) {
                        daughter1Miss <- sum(!sapply(daughters1, function(x) {
                            all(hostTreesSampled[[tree1]]$tip.label[Descendants(hostTreesSampled[[tree1]], x)[[1]]] %in% treeDiff)
                        })) > 1
                    } else {
                        daughter1Miss <- TRUE
                    }
                    
                    isMatch <- sapply(hostNodes, function(node2) {
                        y <- node2
                        for(i in 1:(contrastLevels[[contrast]] + 1)) {
                            y <- c(Ancestors(hostTreesSampled[[tree2]], y[[1]], 'parent'), y)
                        }
                        containingNodes2 <- min(y)
                        
                        nodeUnion <- union(hostTreesSampled[[tree1]]$tip.label[Descendants(hostTreesSampled[[tree1]], hostNodes[[node1]])[[1]]],
                                           hostTreesSampled[[tree2]]$tip.label[Descendants(hostTreesSampled[[tree2]], node2)[[1]]])
                        nodeIntersect <- intersect(hostTreesSampled[[tree1]]$tip.label[Descendants(hostTreesSampled[[tree1]], hostNodes[[node1]])[[1]]],
                                                   hostTreesSampled[[tree2]]$tip.label[Descendants(hostTreesSampled[[tree2]], node2)[[1]]])
                        nodeDiffs <- setdiff(nodeUnion,
                                             nodeIntersect)
                                             
                        containingUnion <- union(hostTreesSampled[[tree1]]$tip.label[Descendants(hostTreesSampled[[tree1]],
                                                                                                 containingNodes1)[[1]]],
                                                 hostTreesSampled[[tree2]]$tip.label[Descendants(hostTreesSampled[[tree2]],
                                                                                                 containingNodes2)[[1]]])
                        containingIntersect <- intersect(hostTreesSampled[[tree1]]$tip.label[Descendants(hostTreesSampled[[tree1]],
                                                                                                         containingNodes1)[[1]]],
                                                         hostTreesSampled[[tree2]]$tip.label[Descendants(hostTreesSampled[[tree2]],
                                                                                                         containingNodes2)[[1]]])
                        containingDiffs <- setdiff(containingUnion,
                                                   containingIntersect)
                        
                        daughters2 <- Children(hostTreesSampled[[tree2]], node2)
                        if(length(daughters2) > 0) {
                            daughter2Miss <- sum(!sapply(daughters2, function(x) {
                                all(hostTreesSampled[[tree2]]$tip.label[Descendants(hostTreesSampled[[tree2]], x)[[1]]] %in% treeDiff)
                            })) > 1
                        } else {
                            daughter2Miss <- TRUE
                        }
                        
                        return(length(nodeIntersect) > 0 & length(containingIntersect) > 0 & all(c(nodeDiffs, containingDiffs) %in% treeDiff) & daughter1Miss & daughter2Miss)
                        # there is at least some overlap in the descendants of each node and some overlap in each node's containing clade; all the differences that do exist are due to species that are simply not present in one of the trees; and the allowable differences due to absent species does not include the case where the extra species form a basal clade under the node of interest
                    })
                    if(sum(isMatch) > 0) {
                        allMatches <- rbind(allMatches,
                                            cbind(paste0('tree', tree1),
                                                  colnames(hostAncestors[[tree1]])[[node1]],
                                                  paste0('tree', tree2),
                                                  colnames(hostAncestors[[tree2]])[isMatch]))
                    }
                }
            }
        }
    }


    tree1Comps <- allMatches[allMatches[,1] == 'tree1', -1]
    tree1Nodes <- unique(tree1Comps[,1])

    inAll <- sapply(tree1Nodes, function(node) {
        sum(tree1Comps[,1] == node) == NTrees - 1
    })
    NNodesInAll <- sum(inAll)

    combinedEffects <- array(NA,
                             dim = c(NMCSamples,
                                     NChains * NTrees,
                                     NNodesInAll,
                                     NMicrobeNodes + 1),
                             dimnames = list(sample = NULL,
                                             chain = NULL,
                                             hostnode = tree1Nodes[inAll],
                                             microbenode = c('alphaDiversity',
                                                             colnames(microbeAncestors))))

    #location effects
    for(tree in 1:NTrees) {
        load(file.path(subdir, paste0('tree_', tree), 'data', 'scaledMicrobeNodeEffects.RData'))
        hostMat <- diag(NHostNodes + NEffects + 1)
        hostMat[(NEffects + 2):(NEffects + NHostNodes + 1),
                (NEffects + 2):(NEffects + NHostNodes + 1)] <- makeContrastMat(NHostNodes,
                                                                               hostTreesSampled[[tree]],
                                                                               NHostTips,
                                                                               hostTreesSampled[[tree]]$tip.label,
                                                                               contrastLevels[[contrast]])[-1, -1]
        hostMat <- Matrix(hostMat)
        effectNames <- c('microbePrevalence',
                         colnames(modelMat)[2:(NEffects + 1)],
                         paste0('host_', colnames(hostAncestors[[tree]])))
        matMult <- array(NA,
                         dim = c(NMCSamples,
                                 NChains,
                                 NEffects + NHostNodes + 1,
                                 NMicrobeNodes + 1),
                         dimnames = list(sample = NULL,
                                         chain = NULL,
                                         hostnode = effectNames,
                                         microbenode = c('alphaDiversity',
                                                         colnames(microbeAncestors))))
        for(j in 1:NMCSamples) {
            for(k in 1:NChains) {
                matMult[j,k,,] <- as.matrix(hostMat %*%
                                            scaledMicrobeNodeEffects[j,k,,] %*%
                                            microbeMat)
            }
        }
        for(node in 1:length(tree1Nodes)) {
            if(inAll[[node]]) {
                if(tree == 1) {
                    nodeName <- tree1Nodes[[node]]
                } else {
                    nodeName <- tree1Comps[tree1Comps[,1] == tree1Nodes[[node]] &
                                           tree1Comps[,2] == paste0('tree', tree),
                                           3]
                }
                combinedEffects[,
                                NChains * (tree - 1) + 1:NChains,
                                tree1Nodes[[node]],
                                ] <- matMult[,, paste0('host_', nodeName),]
            }
        }
    }

    save(combinedEffects, file = file.path(currdatadir, paste0(contrast, '_combinedEffects.RData')))

    allRes <- NULL
    for(j in 1:NNodesInAll) {
        temp <- monitor(array(combinedEffects[,,j,],
                              dim = c(NMCSamples,
                                      NChains * NTrees,
                                      NMicrobeNodes + 1)),
                        warmup = warmup,
                        probs  = c(0.05, 0.95),
                        print  = F)
        temp <- cbind(hostEffect = tree1Nodes[inAll][[j]], temp)
        rownames(temp) <- c('alphaDiversity', rownames(microbeAncestors))
        allRes <- rbind(allRes, temp)
        statusUpdate(j, NNodesInAll)
    }

    cat('microbeNode\t', file = file.path(currtabledir, 'nodeEffects', paste0(contrast, '_matchedNodes.txt')))
    write.table(allRes,
                file   = file.path(currtabledir, 'nodeEffects', paste0(contrast, '_matchedNodes.txt')),
                sep    = '\t',
                quote  = F,
                append = T)
    ##

    ##variance effects
    for(tree in 1:NTrees) {
        metaScales <- array(extract(fit[[tree]],
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
        phyloLogVarMultPrev <- array(extract(fit[[tree]],
                                                     pars       = 'phyloLogVarMultPrev',
                                                     permuted   = F,
                                                     inc_warmup = T),
                                             dim = c(NMCSamples,
                                                     NChains,
                                                     NMicrobeNodes),
                                             dimnames = list(sample  = NULL,
                                                             chain   = NULL,
                                                             taxnode = colnames(microbeAncestors)))
        phyloLogVarMultADiv <- array(extract(fit[[tree]],
                                             pars       = 'phyloLogVarMultADiv',
                                             permuted   = F,
                                             inc_warmup = T),
                                     dim = c(NMCSamples,
                                             NChains,
                                             NHostNodes),
                                     dimnames = list(sample  = NULL,
                                                     chain   = NULL,
                                                     taxnode = colnames(hostAncestors[[tree]])))
        
        phyloLogVarMultRaw <- array(extract(fit[[tree]],
                                            pars       = 'phyloLogVarMultRaw',
                                            permuted   = F,
                                            inc_warmup = T),
                                    dim = c(NMCSamples,
                                            NChains,
                                            NHostNodes,
                                            NMicrobeNodes),
                                    dimnames = list(sample      = NULL,
                                                    chain       = NULL,
                                                    hostnode    = colnames(hostAncestors[[tree]]),
                                                    microbenode = colnames(microbeAncestors)))
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
        hostMat <- Matrix(makeContrastMat(NHostNodes,
                                          hostTreesSampled[[tree]],
                                          NHostTips,
                                          hostTreesSampled[[tree]]$tip.label,
                                          contrastLevels[[contrast]]))
        effectNames <- c('microbePrevalence',
                         paste0('host_', colnames(hostAncestors[[tree]])))
        matMult <- array(NA,
                         dim = c(NMCSamples,
                                 NChains,
                                 NHostNodes + 1,
                                 NMicrobeNodes + 1),
                         dimnames = list(sample = NULL,
                                         chain = NULL,
                                         hostnode = effectNames,
                                         microbenode = c('alphaDiversity',
                                                         colnames(microbeAncestors))))
        for(j in 1:NMCSamples) {
            for(k in 1:NChains) {
                matMult[j,k,,] <- as.matrix(hostMat %*%
                                            phyloLogVarMultScaled[j,k,,] %*%
                                            microbeMat)
            }
        }
        for(node in 1:length(tree1Nodes)) {
            if(inAll[[node]]) {
                if(tree == 1) {
                    nodeName <- tree1Nodes[[node]]
                } else {
                    nodeName <- tree1Comps[tree1Comps[,1] == tree1Nodes[[node]] &
                                           tree1Comps[,2] == paste0('tree', tree),
                                           3]
                }
                combinedEffects[,
                                NChains * (tree - 1) + 1:NChains,
                                tree1Nodes[[node]],
                                ] <- matMult[,, paste0('host_', nodeName),]
            }
        }
    }

    save(combinedEffects, file = file.path(currdatadir, paste0(contrast, '_combined_phyloVar.RData')))

    allRes <- NULL
    for(j in 1:NNodesInAll) {
        temp <- monitor(array(combinedEffects[,,j,],
                              dim = c(NMCSamples,
                                      NChains * NTrees,
                                      NMicrobeNodes + 1)),
                        warmup = warmup,
                        probs  = c(0.05, 0.95),
                        print  = F)
        temp <- cbind(hostEffect = tree1Nodes[inAll][[j]], temp)
        rownames(temp) <- c('alphaDiversity', rownames(microbeAncestors))
        allRes <- rbind(allRes, temp)
        statusUpdate(j, NNodesInAll)
    }

    cat('microbeNode\t', file = file.path(currtabledir, 'phyloVarianceEffects', paste0(contrast, '_matchedNodes.txt')))
    write.table(allRes,
                file   = file.path(currtabledir, 'phyloVarianceEffects', paste0(contrast, '_matchedNodes.txt')),
                sep    = '\t',
                quote  = F,
                append = T)
}
