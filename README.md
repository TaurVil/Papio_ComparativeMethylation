# Papio_ComparativeMethylation
Re-analysis of DNA methylation patterns in papio baboons, published in Molecular Biology and Evolution 2019 (https://doi.org/10.1093/molbev/msy227) and adjusted to the Panubis1 baboon genome. Please contact me (taur.vil at gmail.com) with any questions. 

## Goal
Redo mapping and analysis of published results for the new genome version. 
Useful for followup work in other studies to identify overlap in differentially methylated sites

## Processing Scripts

set.01.sh : get and map data, quality metrics

set.02.sh : call genotypes from bam file. Produces the file "cleaned.genotypes_matrix.txt" for 012 matrices and the merged.vcf.gz file. Contains all 20 autosomes and the 2 sex chromosomes. 

set.03.sh : get count data for the autosomes. Remove sites that are variable in the bisSNP calls or based on the methratio files. n44 info and count files, saved within n44.raw_count_data.RData

After step 3, I did a lot of cleaning. Statistics and intermediate count files and bams were stored in "archived" and moved to a local hard drive. Fastq files (trimmed and raw) were removed and can be re-downloaded from SRA following set.01.sh. Genotypes by chormosome were removed and can be extracted from the total genotype file. Methratio folder was removed after zipping into a tar.gz file called "methratio_files.tar.gz" and stored in the archive. 

set.04.sh : hardac based models, namely regressing counts out, ANOVA, and macau

set.05.sh : define genomic contexts

## Analysis Scripts

GenotypeStructure.Rmd : plots but no additional files


