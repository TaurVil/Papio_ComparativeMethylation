#!/bin/bash
path_genome=/data/tunglab/shared/genomes/panubis1/Panubis_1.0.fa

index=${SLURM_ARRAY_TASK_ID}


## This will work as the chromosomes are indexed 1-20,X,Y,scafoldZZZ rather than "chr1", "chr2", etc.
chrom=$index

module load java/1.8.0_45-fasrc01
module load tabix
module load samtools
module load vcftools
export PATH=$PATH:/data/tunglab/tpv/Programs/cmake/bin/
module load gcc

# Get genotype calls that in either reference panel, then identify invariant sites and prepare to remove singletons
vcftools --gzvcf /data/tunglab/asf40/wgs_data/tmp2/vcf_files/baboon1k_v1_snpEff_chr$chrom'.vcf.gz' --keep 00_allspecies.list --max-alleles 2 --remove-indels --recode --out ./chrom_vcfs/01.raw_biallelic_snps.$chrom.unadmixed --recode-INFO-all

vcftools --gzvcf ./chrom_vcfs/01.raw_biallelic_snps.$chrom.unadmixed.recode.vcf --max-missing 0.5 --minQ 30 --max-alleles 2 --remove-indels --maf 0.01 --recode --out ./chrom_vcfs/01.$chrom.missing_minQ_biallelic.unadmixed --recode-INFO-all

vcftools --gzvcf ./chrom_vcfs/01.$chrom.missing_minQ_biallelic.unadmixed.recode.vcf --singletons --out ./chrom_vcfs/01.$chrom.unadmixed


# GATK variant filtering
java -jar /data/tunglab/tpv/Programs/GenomeAnalysisTK.jar -T VariantFiltration -R ./Panubis1.nochromname.fa -V ./chrom_vcfs/01.$chrom.missing_minQ_biallelic.unadmixed.recode.vcf -filterName "FS" --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" -o ./chrom_vcfs/01b.gatk_filtered.unadmixed.$chrom.vcf.gz

# Get separate vcf for each species, retaining only singleton sites with adequate coverage
vcftools --gzvcf ./chrom_vcfs/01b.gatk_filtered.unadmixed.$chrom.vcf.gz --exclude-positions ./chrom_vcfs/01.$chrom.unadmixed.singletons --recode --out ./chrom_vcfs/02.$chrom --recode-INFO-all --remove-filtered-all

bgzip ./chrom_vcfs/01.$chrom.missing_minQ_biallelic.unadmixed.recode.vcf
bgzip ./chrom_vcfs/02.$chrom.recode.vcf
bgzip ./chrom_vcfs/01.raw_biallelic_snps.$chrom.unadmixed.recode.vcf

vcftools --gzvcf ./chrom_vcfs/01.$chrom.missing_minQ_biallelic.unadmixed.recode.vcf.gz --max-missing 1 --thin 50 --recode --out ./chrom_vcfs/02.no_missing.thin50.$chrom
bgzip ./chrom_vcfs/02.no_missing.thin50.$chrom.recode.vcf
tabix ./chrom_vcfs/02.no_missing.thin50.$chrom.recode.vcf.gz 

