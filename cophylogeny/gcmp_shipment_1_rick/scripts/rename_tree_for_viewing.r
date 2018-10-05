microbeTree.Y.root.temp <- microbeTree.Y.root
microbeTree.Y.root.temp$tip.label <- tax[microbeTree.Y.root.temp$tip.label,'Family']
microbeTree.Y.root.temp$tip.label[is.na(microbeTree.Y.root.temp$tip.label)] <- 'NA'

for(i in (length(microbeTree.Y.root.temp$tip.label) + 1):(length(microbeTree.Y.root.temp$tip.label) + microbeTree.Y.root.temp$Nnode)) {
    desc <- Descendants(microbeTree.Y.root.temp, i)[[1]]
    
    if(all(microbeTree.Y.root.temp$tip.label[desc] == microbeTree.Y.root.temp$tip.label[desc[[1]]])) {
        microbeTree.Y.root.temp$tip.label[desc[[1]]] <- paste0(microbeTree.Y.root.temp$tip.label[desc[[1]]], '.', length(desc))
        microbeTree.Y.root.temp <- drop.tip(microbeTree.Y.root.temp, desc[-1])
    }
}

microbeTree.Y.root.temp$tip.label <- paste(1:length(microbeTree.Y.root.temp$tip.label), microbeTree.Y.root.temp$tip.label)

write.tree(microbeTree.Y.root.temp, file='~/Dropbox/0.tree.tree')
