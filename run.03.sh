# Count data

# Convert individual mratio files into a single matrix for methylated and unmethylated counts
  cd /data/tunglab/tpv/CM
  ls methratio/*methratio.txt > methratio/methratio_files.txt
  for f in `cat 00_chroms.txt`; do cat run.c1.get_counts_part1.sh | sed -e s/CHROMNAME/$f/g > get_data.$f.sh; sbatch -N 1 --mem=10000 -p all --nice ./get_data.$f.sh; done
  cd methratio; module load R; for f in `cat ../00_chroms.txt`; do cat ../run.c2.get_counts_part2.R | sed -e s/CHROMNAME/$f/g > $f.get_data.R; sbatch -N 1 --mem=10000 -p all --nice ./$f.get_data.R; done; cd ..

# Merge into one file for counts, mcounts, and info

# Remove CpG sites where the C or G is variable
  # Variable sites determined using BisSNP genotype calls for these samples, rather than an externally validated file set
  # found in /data/tunglab/tpv/CM/genotype_all/merged.012.pos
  
  # Other option is to identify them from known variants. We'll not do this. 
	# vcftools --gzvcf known_variants.vcf.gz --012 --out known_variants
	# vcftools --gzvcf known_variants.vcf.gz --remove-filtered-all --mac 2 --max-alleles 2 --minQ 100 --012 --out filtered_variants
			
  # Because C and G counts are condensed, we'll need to remove sites where the variant overlaps the CpG site value or the CpG site value plus 1 (disrupting the C or G)
  
  
