# Call genotypes from bam files 
	
  ## Could use a dbsnp file of known baboon variants, but that's not available yet (will be with Vilgalys, Fogel et al.)
	
  sbatch --array=1-44 --mem=8G run.a1.sort_convert_RG.sh
	sbatch --array=1-22 --mem=16G run.a2.bisSNP.sh
	ls genotype_by_chrom/*gz > a1_vcf_by_chrom.list
	java -jar ~/picard.jar GatherVcfs I=a1_vcf_by_chrom.list O=merged.vcf.gz
	vcftools --gzvcf merged.vcf.gz --012 --out merged

# Return vcf file and 012 genotype matrix

# Fix 012 genotype matrix
cd /data/tunglab/tpv/CM/genotype_all
module load R
library(data.table)
fread("merged.012") -> d; t(d) -> d
fread("merged.012.pos") -> pos; pos$site <- paste(pos$V1 , pos$V2, sep="_")
as.data.frame(d) -> d; d[-1,] -> d
pos$site -> row.names(d)
write.table(d, "cleaned.genotype_matrix.txt", row.names=T, col.names=T, sep="\t", quote=F)
