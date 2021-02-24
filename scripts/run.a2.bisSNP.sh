#!/bin/bash

module load java
module load samtools
module load tabix

#ls bams/*.bam > a1_bam.list

path_genome=/data/tunglab/shared/genomes/panubis1/Panubis1.0.fa
index=${SLURM_ARRAY_TASK_ID}

chrom=`head -$index 00_chroms.txt | tail -1`


#java -Xmx4g -jar /data/tunglab/tpv/Programs/BisSNP-0.82.2.jar -R $path_genome -T BisulfiteGenotyper -I a1_bam.list -vfn1 cpg_1.raw.$chrom.vcf -vfn2 snp_1.$chrom.raw.vcf -stand_call_conf 20 -mmq 30 -mbq 0 -L $chrom
                        #Lower stand_call_conf for low coverage samples
                        #mmq is the minimum read quality score necessary
                        #mbq is the minimum base quality score necessary

#java -Xmx4g -jar /data/tunglab/tpv/Programs/BisSNP-0.82.2.jar -R $path_genome -T VCFpostprocess -oldVcf snp_1.$chrom.raw.vcf -newVcf snp_1.filtered.$chrom.vcf -snpVcf snp_1.$chrom.raw.vcf -o snp_1.filter.$chrom.summary.txt

##$path_java -Xmx4g -jar /data/tunglab/tpv/SGE/RRBS/bams/BisSNP-0.82.2.jar -R $path_genome -T VCFpostprocess -oldVcf cpg_1.raw.CHROM.vcf -newVcf cpg_1.filtered.CHROM.vcf -snpVcf snp_1.raw.CHROM.vcf -L CHROM -o cpg_1.filter.CHROM.summary.txt
#

#sed '4d' snp_1.filtered.$chrom.vcf > tmp.$chrom; mv tmp.$chrom snp_1.filtered.$chrom.vcf
#bgzip cpg_1.raw.CHROM.vcf; tabix cpg_1.raw.CHROM.vcf.gz
#bgzip snp_1.filtered.$chrom.vcf; tabix snp_1.filtered.$chrom.vcf.gz

module load java/1.8.0_45-fasrc01
module load vcftools
#java -Xmx4g -jar /data/tunglab/tpv/Programs/GenomeAnalysisTK.jar -R $path_genome -T VariantFiltration -V snp_1.filtered.$chrom.vcf.gz -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" -o snp_1.2xfilt.$chrom.vcf

vcftools --vcf snp_1.2xfilt.$chrom.vcf --maf 0.05 --minQ 25 --max-missing 0.10 --recode --out snp_1.3xfil.$chrom.vcf

bgzip snp_1.3xfil.$chrom.vcf.recode.vcf
mv snp_1.3xfil.$chrom.vcf.recode.vcf.gz snp.3xfil.maf05.maxmissing10.minQ25.$chrom.vcf.gz

tabix snp.3xfil.maf05.maxmissing10.minQ25.$chrom.vcf.gz

#bgzip snp_1.3xfil.CHROM.vcf
#bgzip snp_1.2xfilt.CHROM.vcf

