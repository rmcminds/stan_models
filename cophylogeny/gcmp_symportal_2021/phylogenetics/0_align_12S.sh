# fasta to nexus converter from https://sites.google.com/site/shannonhedtke/Scripts

align_seqs.py -i raw_data/host_phylo/combined_mitochondria_minimal.fasta -o output/host_phylo/combined_mitochondria_minimal_aligned -m mafft

awk '1 {print $1}' output/host_phylo/combined_mitochondria_minimal_aligned/*.fasta | sed -E '/>/!s/[^ATCG]/?/g' > output/host_phylo/combined_mitochondria_minimal_aligned/combined_mitochondria_minimal_aligned_clean.fasta

convertfasta2nex.pl output/host_phylo/combined_mitochondria_minimal_aligned/combined_mitochondria_minimal_aligned_clean.fasta > output/host_phylo/combined_mitochondria_minimal_aligned/combined_mitochondria_minimal_aligned_clean.nex

