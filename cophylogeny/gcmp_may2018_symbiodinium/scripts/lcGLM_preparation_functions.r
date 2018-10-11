findCutPoints <- function(tree, maxNH, NSplits, NTimeBins) {
    
    #divide total evolutionary time into chunks that contain approximately equal numbers of splits
    lttTree <- ltt(tree, log.lineages = F, plot = F)
    temp <- maxNH - lttTree$times[-length(lttTree$times)]
    splitTimes <- split(temp, ceiling(seq_along(temp) / NSplits))
    if(length(splitTimes[[NTimeBins]]) < NSplits / 2) {
        splitTimes[[NTimeBins - 1]] <- c(splitTimes[[NTimeBins - 1]],
                                         splitTimes[[NTimeBins]])
        splitTimes[[NTimeBins]] <- NULL
        NTimeBins <- NTimeBins - 1
    }
    
    #cut points in each phylogeny that would result in approximately equal numbers of splits per bin of time
    boundaries <- sapply(2:NTimeBins, function(x) max(splitTimes[[x]]))
    boundariesRounded <- round(boundaries, 1)
    
    return(list(NTimeBins         = NTimeBins,
                boundaries        = boundaries,
                boundariesRounded = boundariesRounded))
    
}

createTimeBins <- function(cutPoints, maxNH, NHs, tree, edgeOrder) {
    
    #size of each bin
    timeBinSizes <- sapply(2:(length(cutPoints$boundaries) + 2),
                           function(x) c(maxNH,
                                         cutPoints$boundaries,
                                         0)[x-1] - c(maxNH,
                                                     cutPoints$boundaries,
                                                     0)[x])
    #proportion of total tree height belonging to each bin (differs for each replicate due to the varying sizes of the oldest bin)
    relativeTimeBinSizes <- timeBinSizes / sum(timeBinSizes)
    
    nd <- maxNH - NHs
    
    #create a matrix for each replicate that describes the proportion of each time bin that each branch exists for
    edgeToBin <- matrix(NA, ncol = cutPoints$NTimeBins, nrow = nrow(tree$edge))
    for(j in 1:cutPoints$NTimeBins) {
        if(j == 1) {
            allin  <- which(nd[,2] >= cutPoints$boundaries[j])
            allout <- which(nd[,1] <= cutPoints$boundaries[j])
            cedge  <- which((nd[,1] > cutPoints$boundaries[j]) &
                            (nd[,2] < cutPoints$boundaries[j]))
            edgeToBin[cedge, j] <- nd[cedge, 1] - cutPoints$boundaries[j]
        } else if(j == cutPoints$NTimeBins) {
            allin  <- which(nd[,1] <= cutPoints$boundaries[j - 1])
            allout <- which(nd[,2] >= cutPoints$boundaries[j - 1])
            cedge  <- which((nd[,1] > cutPoints$boundaries[j - 1]) &
                            (nd[,2] < cutPoints$boundaries[j - 1]))
            edgeToBin[cedge, j] <- cutPoints$boundaries[j - 1] - nd[cedge, 2]
        } else {
            allin  <- which((nd[,1] <= cutPoints$boundaries[j - 1]) &
                            (nd[,2] >= cutPoints$boundaries[j]))
            allout <- which((nd[,1] <= cutPoints$boundaries[j]) |
                            (nd[,2] >= cutPoints$boundaries[j - 1]))
            cedge1 <- which((nd[,1] <= cutPoints$boundaries[j - 1]) &
                            (nd[,1] > cutPoints$boundaries[j]) &
                            (nd[,2] < cutPoints$boundaries[j]))
            cedge2 <- which((nd[,1] > cutPoints$boundaries[j - 1]) &
                            (nd[,2] < cutPoints$boundaries[j - 1]) &
                            (nd[,2] >= cutPoints$boundaries[j]))
            cedge3 <- which((nd[,1] > cutPoints$boundaries[j - 1]) &
                            (nd[,2] < cutPoints$boundaries[j]))
                            
            edgeToBin[cedge1, j] <- nd[cedge1, 1] - cutPoints$boundaries[j]
            edgeToBin[cedge2, j] <- cutPoints$boundaries[j - 1] - nd[cedge2, 2]
            edgeToBin[cedge3, j] <- cutPoints$boundaries[j - 1] - cutPoints$boundaries[j]
        }
        edgeToBin[allin, j] <- tree$edge.length[allin]
        edgeToBin[allout, j] <- 0
    }
    edgeToBin <- edgeToBin / maxNH
    rownames(edgeToBin) <- tree$edge[,2]
    edgeToBin <- edgeToBin[edgeOrder,]

    return(list(edgeToBin            = edgeToBin,
                timeBinSizes         = timeBinSizes,
                relativeTimeBinSizes = relativeTimeBinSizes))
    
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

getTreeDetails <- function(tree) {
    edgeOrder <- order(tree$edge[,2])
    NHs <- nodeHeights(tree)
    maxNHs <- max(NHs)
    NHRel <- NHs / maxNHs
    rownames(NHRel) <- tree$edge[,2]
    NHRel <- NHRel[edgeOrder,]
    edgeLengths <- tree$edge.length[edgeOrder]
    return(list(edgeOrder   = edgeOrder,
                NHs         = NHs,
                maxNHs      = maxNHs,
                NHRel       = NHRel,
                edgeLengths = edgeLengths))
}
