# Papio_ComparativeMethylation
Re-analysis of DNA methylation patterns in papio baboons, published in Molecular Biology and Evolution 2019 (https://doi.org/10.1093/molbev/msy227) and adjusted to the Panubis1 baboon genome. Please contact me (taur.vil at gmail.com) with any questions. 

## Goal
Redo mapping and analysis of published results for the new genome version. 
Useful for followup work in other studies to identify overlap in differentially methylated sites

## Processing Scripts

**set.01.sh**: get and map data, quality metrics

**set.02.sh**: call genotypes from bam file. Produces the file "cleaned.genotypes_matrix.txt" for 012 matrices and the merged.vcf.gz file. Contains all 20 autosomes and the 2 sex chromosomes. 

**set.03.sh**: get count data for the autosomes. Remove sites that are variable in the bisSNP calls or based on the methratio files. n44 info and count files, saved within **n44.raw_count_data.RData**.

After step 3, I did a lot of cleaning. Statistics and intermediate count files and bams were stored in "archived" and moved to a local hard drive. Fastq files (trimmed and raw) were removed and can be re-downloaded from SRA following set.01.sh. Genotypes by chromosome were removed and can be extracted from the total genotype file. Methratio folder was removed after zipping into a tar.gz file called "methratio_files.tar.gz" and stored in the archive. 

**set.04.sh**: hardac based models, namely regressing counts out, ANOVA, and macau

**set.05.sh**: define genomic contexts: CpG islands, CpG shores, genes, promoters, enhancers, regions where the local phylogeny matches the global phylogeny. 


## Local Analysis Scripts

**01_GenotypeStructure.Rmd** : Fig S8. Includes other QC plots, but no output files. 

**02_Convert_Raw_Counts_to_Continuous_Data.R** : Calculates ratios, counts regressed out, and mvalues from n44.raw_count_data.Rmd. Here we discover a subset of sites with methylation level > 1, which we remove and refilter the info and count file for. We are left with 979,099 sites with >= 5x coverage in each species; we had already filtered for sites covered in 1/2 the dataset. The outcome file is saved as **n44.continuous_data.RData**; we're done with the raw data file (n44.raw_count_data.RData) and that's moved to the archive. 

**03_MethylationStructure.Rmd** : Fig 1A & B. Fig S2. Procrustes analysis. 

**04_Anova.R** : ANOVA analysis for sites with intermediate methylation (with 10 permutations). Note that our requirement of 5x coverage per species retains ~250k sites with intermediate methylation vs ~700k in the original paper (which required 5x coverage across all samples only). Saves output (ANOVA results only) as **n44.anova.RData**.
