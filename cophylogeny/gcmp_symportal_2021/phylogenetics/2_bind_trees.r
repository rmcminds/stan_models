library(ape)
library(phytools)
library(phangorn)
library(worms)

HRmolTree <- read.nexus('raw_data/host_phylo/Huang&Roy_Molecular.tre')
coralCorallimorphSplit <- max(nodeHeights(HRmolTree))

nsamples <- 25

outgroupTrees <- read.nexus('output/host_phylo/linked_tree.trees')
outgroupTreesBurnedIn <- outgroupTrees[50002:100001]
outgroupTreesSampled <- sample(outgroupTreesBurnedIn,nsamples)

HRtrees <- read.nexus('raw_data/host_phylo/Huang&Roy_Supertree.tre')

oldnames <- sub('_',
                ' ',
                sub('..._', '', HRtrees$tip.label[[1]]))
newnames <- sub(' .* ',
                ' ',
                wormsbynames(oldnames, marine_only = F)$valid_name)

fixedByLike <- cbind(oldnames[is.na(newnames)],
                     sub(' .* ',
                         ' ',
                         wormsbynames(sub(' ', '%', oldnames[is.na(newnames)]))$valid_name))

newnames[is.na(newnames)] <- fixedByLike[,2]
fixedByFuzzy <- cbind(oldnames[is.na(newnames)],
                      sapply(oldnames[is.na(newnames)],
                             function(x) tryCatch(wormsbymatchnames(x)$valid_name,
                                                  error = function(e) NA)))
newnames[is.na(newnames)] <- fixedByFuzzy[,2]
noMatches <- oldnames[is.na(newnames)]
newnames[is.na(newnames)] <- noMatches

HRtreesRelabeled <- HRtrees

newnames <- sub(' .* ', '_', newnames)
attr(HRtreesRelabeled, "TipLabel") <- sub(' ', '_', newnames)

HRtreesConsolidated <- lapply(HRtreesRelabeled,
                              function(x) drop.tip(x, which(duplicated(attr(HRtreesRelabeled, "TipLabel")))))
class(HRtreesConsolidated) <- "multiPhylo"



HRtreesSampled <- sample(HRtreesConsolidated,nsamples)

source('gcmp_stan_symportal/phylogenetics/add_species_mod.r')


combined.trees <- lapply(1:nsamples, function(i) {
    
    tempHR <- HRtreesSampled[[i]]
    HR_height <- max(branching.times(tempHR))
    trimmed <- drop.tip(outgroupTreesSampled[[i]], 'Pocillopora_damicornis')
    trimmed$tip.label[trimmed$tip.label=='Acropora_palmata'] <- 'Corals'
    
    bind_to <- Ancestors(trimmed, which(trimmed$tip.label == 'Corals'), type = "parent")
    bind_to_height <- dist.nodes(trimmed)[which(trimmed$tip.label == 'Corals'), bind_to]
    tempHR$root.edge <- bind_to_height - HR_height
    
    full.tree <- bind.tree(trimmed, tempHR, where = bind_to)
    full.tree.trimmed <- drop.tip(full.tree, 'Corals')
    
    bind_toTrimmed <- getMRCA(full.tree.trimmed, tempHR$tip.label)
    bind_toTrimmedHeight <- max(nodeHeights(full.tree.trimmed)) - nodeheight(full.tree.trimmed, bind_toTrimmed)
    bind_toTrimmedParentHeight <- max(nodeHeights(full.tree.trimmed)) - nodeheight(full.tree.trimmed, Ancestors(full.tree.trimmed, bind_toTrimmed, type = "parent"))
    
    if (bind_toTrimmedHeight > coralCorallimorphSplit | bind_toTrimmedParentHeight < coralCorallimorphSplit) {
        position <- runif(1,0, full.tree.trimmed$edge.length[bind_toTrimmed])
    } else {
    	position <- coralCorallimorphSplit - bind_toTrimmedHeight
    }
    
    corallimorphs <- subtrees(HRmolTree)[[getMRCA(HRmolTree,
                                                  c('COR_Ricordea_florida', 'COR_Discosoma')) - length(HRmolTree$tip.label)]]
    corallimorphs$root.edge <- bind_toTrimmedHeight + position - max(branching.times(corallimorphs))
    
    newtree <- bind.tree(full.tree.trimmed, corallimorphs, where = bind_toTrimmed, position = position)
    newtree$tip.label[newtree$tip.label %in% c('COR_Ricordea_florida', 'COR_Discosoma')] <- c('Discosoma_sp', 'Ricordea_florida')
    fullertree <- force.ultrametric(add.species.to.genus.mod(newtree, 'Rhodactis_sp', tips = c('Discosoma_sp', 'Ricordea_florida'), where = 'random'))
    fullertree <- force.ultrametric(add.species.to.genus.mod(fullertree, 'Xenia_umbellata', tips = 'Lobophytum_sp', where = 'random'))
    fullertree <- force.ultrametric(add.species.to.genus.mod(fullertree, 'Sarcophyton_sp', tips = 'Lobophytum_sp', where = 'random'))
    fullertree <- force.ultrametric(add.species.to.genus.mod(fullertree, 'Sinularia_polydactyla', tips = c('Lobophytum_sp', 'Sarcophyton_sp'), where = 'random'))
    fullertree <- force.ultrametric(add.species.to.genus.mod(fullertree, 'Heteractis_aurora', tips = 'Anemonia_viridis', where = 'random'))
    fullertree <- force.ultrametric(add.species.to.genus.mod(fullertree, 'Entacmaea_quadricolor', tips = 'Anemonia_viridis', where = 'random'))
    fullertree <- force.ultrametric(add.species.to.genus.mod(fullertree, 'Aiptasia_sp', tips = 'Nematostella_vectensis', where = 'random'))
    fullertree <- force.ultrametric(add.species.to.genus.mod(fullertree, 'Stylaster_sp', tips = 'Millepora_sp', where = 'random'))
    fullertree <- force.ultrametric(add.species.to.genus.mod(fullertree, 'Macrorhynchia_philippina', tips = c('Millepora_sp', 'Hydra_vulgaris'), where = 'random'))
    fullertree$tip.label[fullertree$tip.label == 'Sycon_ciliatum'] <- 'Mnemiopsis_sp'
	return(fullertree)
})

class(combined.trees) <- "multiPhylo"
combined.trees <- .compressTipLabel(combined.trees)

write.tree(combined.trees,'output/host_phylo/combined_trees.newick')

