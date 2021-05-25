# Call genotypes from bam files 
	
  ## Could use a dbsnp file of known baboon variants, but that's not available yet (will be with Vilgalys, Fogel et al.)
  
  ## Prepare bam files for genotype calling
  sbatch --array=1-44 --mem=8G run.a1.sort_convert_RG.sh
  ## Run BisSNP for each chromosome. 
  sbatch --array=1-22 --mem=16G run.a2.bisSNP.sh
  
  ## Bis SNP is called with minimum read quality score of 30, minimum base quality score of 0, and standard call confidence of 20
  ## Includes filtering for 5% minor allele frequency, called in at least 10% of samples, minQ > 25, thined of clusters (3 SNPs within 35bp), and filtered following GATK hard filters (FS > 30, QD < 2)
  ## The 4th line from BisSNP neeeds to be removed to run GATK for variant filtration
  
  ## Merge genotype calls across chromosomes
  ls genotype_by_chrom/*gz > a1_vcf_by_chrom.list
  java -jar ~/picard.jar GatherVcfs I=a1_vcf_by_chrom.list O=merged.vcf.gz
  ## Get 012 genotype matrix to accompany the actual genotype calls
  vcftools --gzvcf merged.vcf.gz --012 --out merged

# Return vcf file and 012 genotype matrix

# Fix 012 genotype matrix. 
  ## Original output is sample-by-site rather than site-by-sample, and includes a placeholder row which we'll get rid of. 
cd /data/tunglab/tpv/CM/genotype_all
module load R; R
	library(data.table)
	fread("merged.012") -> d; t(d) -> d
	fread("merged.012.pos") -> pos; pos$site <- paste(pos$V1 , pos$V2, sep="_")
	as.data.frame(d) -> d; d[-1,] -> d
	pos$site -> row.names(d)
	write.table(d, "cleaned.genotype_matrix.txt", row.names=T, col.names=T, sep="\t", quote=F)

sed -i 's/-1/NA/g' cleaned.genotype_matrix.txt

## at this point there are 963,160 variants in the "genotype_all" folder
