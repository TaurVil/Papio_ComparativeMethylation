# Papio_ComparativeMethylation
Re-analysis of DNA methylation patterns in papio baboons, published in Molecular Biology and Evolution 2019 (https://doi.org/10.1093/molbev/msy227) and adjusted to the Panubis1 baboon genome. Please contact me (taur.vil at gmail.com) with any questions. 

## Goal
Redo mapping and analysis of published results for the new genome version. 
Useful for followup work in other studies to identify overlap in differentially methylated sites

## Processing Scripts

set.01.sh : get and map data, quality metrics

set.02.sh : call genotypes from bam file. Produces the file "cleaned.genotypes_matrix.txt"

set.03.sh : get count data for the autosomes. Remove sites that are variable in the bisSNP calls or based on the methratio files. 

