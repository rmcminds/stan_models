library(ape)

args <- commandArgs(trailingOnly = TRUE)
collapseMultis <- args[[1]]
model <- args[[2]]
input <- args[[3]]
output <- args[[4]]

tree <- read.tree(input)
if(collapseMultis) {tree <- di2multi(tree)}

originalmin <- min(tree$edge.length[tree$edge.length > 0])
tree$edge.length[tree$edge.length <= 0] <- originalmin / 2

newtree <- tree
newtree$edge.length <- tree$edge.length * 0.001 / originalmin

caltree <- ladderize(chronos(newtree, model=model))

write.tree(caltree, file=output)
