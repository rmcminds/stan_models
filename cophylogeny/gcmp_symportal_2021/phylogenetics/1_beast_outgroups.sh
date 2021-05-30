#calibrated_augmented_host_phylogeny.xml created in beauti using results of align_12s.sh and pratlong alignment

beast -beagle -threads 10 output/host_phylo/calibrated_augmented_host_phylogeny.xml
treeannotator -b 50 -heights ca -lowMem linked_tree.trees final_tree.tree
