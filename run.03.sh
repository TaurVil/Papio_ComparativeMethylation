# Count data

# Convert individual mratio files into a single matrix for methylated and unmethylated counts
  cd /data/tunglab/tpv/CM
  ls methratio/*methratio.txt > methratio/methratio_files.txt
  module load bedtools2; for f in `cat 00_chroms.txt`; do cat run.c1.get_counts_part1.sh | sed -e s/CHROMNAME/$f/g > get_data.$f.sh; sbatch -N 1 --mem=10000 -p all --nice ./get_data.$f.sh; done
  conda activate r-env 
  cd methratio; for f in `cat ../00_chroms.txt`; do cat ../run.c2.get_counts_part2.R | sed -e s/CHROMNAME/$f/g > $f.get_data.R; sbatch -N 1 --mem=10000 -p all --nice ./$f.get_data.R; done; cd ..

# Merge into one file for counts, mcounts, and info
  head -1 info_chr1.txt > n44.info.txt
  for f in `seq 1 20`; do cat counts_table_chr$f.txt >> n44.counts.txt;  
    sed '1d' info_chr$f.txt | sed 's/"//g' >> n44.info.txt; 
    sed '1d' mcounts_table_chr$f.txt >> n44.mcounts.txt; done

# only keeping the autosomes

# Remove CpG sites where the C or G is variable
  # Because C and G counts are condensed, we'll need to remove sites where the variant overlaps the CpG site value or the CpG site value plus 1 (disrupting the C or G)

  # Variable sites determined using BisSNP genotype calls for these samples, rather than an externally validated file set
  # found in /data/tunglab/tpv/CM/genotype_all/merged.012.pos
  
  # Remove sites where the methylation data implies polymorphism (reverse G counts are not the same as reverse GA counts)
  module load R; R
  options(scipen = 20)
  library(data.table); to_remove <- NULL; files <- list.files(pattern="methratio.txt")
  for (i in files) {
    fread(i) -> data; subset(data, ! data$rev_G_count == data$rev_GA_count) -> d2
    rbind(d2, to_remove) -> to_remove
  }
  to_remove[!duplicated(paste(to_remove$chr, to_remove$pos, sep="_")),] -> to_rem
  to_rem[,3] <- to_rem[,2] + 1
  write.table(to_rem[,1:3], "./variable_CT.bed", row.names=F, col.names=F, sep="\t", quote=F)
  
  read.delim("/data/tunglab/tpv/CM/genotype_all/merged.012.pos", header=F) -> snps
  snps[,3] <- snps[,2] + 1; snps[,2] <- snps[,2] - 1
  write.table(snps[,1:3], "./snps.bed", row.names=F, col.names=F, sep="\t", quote=F)
  
  module load bedtools2; cat snps.bed variable_CT.bed > remove.bed; bedtools sort -i remove.bed > tmp.bed; bedtools merge -i tmp.bed > to_exclude.bed; rm tmp.bed
  sed '1d' n44.info.txt > tmp; bedtools sort -i tmp > tmp2; bedtools subtract -a tmp2 -b to_exclude.bed | cut -f 1-4 > n44.nonvariable.txt; rm tmp; rm tmp2

  # Other option is to identify them from known variants. We'll not do this. 
	# vcftools --gzvcf known_variants.vcf.gz --012 --out known_variants
	# vcftools --gzvcf known_variants.vcf.gz --remove-filtered-all --mac 2 --max-alleles 2 --minQ 100 --012 --out filtered_variants
			
## Export R data object to use for future analyses
## Remove non-variable sites

  library(data.table); fread("./n44.nonvariable.txt") -> keep

  fread("./n44.info.txt") -> info
  fread("./n44.counts.txt") -> c
  fread("./n44.mcounts.txt") -> m

  m[,-1] -> m; colnames(m) <- colnames(c)

  c[info$site %in% keep$V4,] -> c
  m[info$site %in% keep$V4,] -> m
  info[info$site %in% keep$V4,] -> info

  rm(keep); save.image("../n44.raw_count_data.RData")

  ## Specially pull out sites for yellow-anubis analyses 
  info$n_anu <- 9-rowSums(c[,1:9] == 0)
  info$n_yel <- 6-rowSums(c[,34:39] == 0)

  ## retains 1,014,362 sites out of 1,015,557 
  
### Zip up methratio folder, we should be done with it
tar -cvzf methratio_files.tar.gz methratio/
