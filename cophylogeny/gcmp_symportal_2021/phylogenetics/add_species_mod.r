## slightly modified version of add.species.to.genus() from phytools
add.species.to.genus.mod <- function (tree, species, genus = NULL, tips = NULL, where = c("root", "random"))
{
    if (!inherits(tree, "phylo"))
    stop("tree should be an object of class \"phylo\".")
    if (!is.ultrametric(tree))
    warning("this code has only been tested with ultrametric tree\n  your tree may be returned without edge lengths")
    where <- where[1]
    if (is.null(genus)) {
        x <- strsplit(species, "")[[1]]
        i <- 1
        while (x[i] != "_" && x[i] != " ") i <- i + 1
        genus <- paste(x[2:i - 1], collapse = "")
    }
    if (is.null(tips)) {
        ii <- grep(paste(genus, "_", sep = ""), tree$tip.label)
    } else {ii <- grep(paste(tips,collapse='|'),tree$tip.label)}
    if (length(ii) > 1) {
        if (!is.monophyletic(tree, tree$tip.label[ii]))
        warning(paste(genus, "may not be monophyletic\n  attaching to the most inclusive group containing members of this genus"))
        nn <- findMRCA(tree, tree$tip.label[ii])
        if (where == "root")
        tree <- bind.tip(tree, gsub(" ", "_", species), where = nn)
        else if (where == "random") {
            tt <- splitTree(tree, list(node = nn, bp = tree$edge.length[which(tree$edge[,
            2] == nn)]))
            tt[[2]] <- add.random(tt[[2]], tips = gsub(" ", "_",
            species))
            tree <- paste.tree(tt[[1]], tt[[2]])
        }
        else stop("option 'where' not recognized")
    }
    else if (length(ii) == 1) {
        nn <- ii
        if (where == "root")
        tree <- bind.tip(tree, gsub(" ", "_", species), where = nn,
        position = 0.5 * tree$edge.length[which(tree$edge[,
        2] == nn)])
        else if (where == "random")
        tree <- bind.tip(tree, gsub(" ", "_", species), where = nn,
        position = runif(n = 1) * tree$edge.length[which(tree$edge[,
        2] == nn)])
        else stop("option 'where' not recognized")
    }
    else warning("could not match your species to a genus\n  check spelling, including case")
    tree
}
