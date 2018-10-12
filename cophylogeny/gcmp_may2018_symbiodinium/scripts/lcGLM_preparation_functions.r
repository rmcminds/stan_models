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
