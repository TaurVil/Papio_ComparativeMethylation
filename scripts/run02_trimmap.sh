#!/bin/bash
#SBATCH --get-user-env

# run trim galore
module load TrimGalore

index=${SLURM_ARRAY_TASK_ID}
id=`head -$index 00_libnames.txt | tail -1`
fastq=$id.fq.gz

#trim_galore -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --gzip --rrbs --length 15 --stringency 4 $fastq

ls $id.fq.gz
#mv "$id"_trimmed.fq.gz trimmed/$id.trimmed.fq.gz

module load BSMAP
module load samtools

path_fastq=trimmed/$id.trimmed.fq.gz
path_genome=/data/tunglab/shared/genomes/panubis1/Panubis1.0.fa
path_sam=$id.sam

bsmap -a $path_fastq -d $path_genome -o $path_sam -3 -v 0.1 -r 0

# get mratios
module load python/2.7.6-fasrc01
module load samtools/1.3.1-gcb01
path_out=methratio/$id.methratio.txt
#python /data/tunglab/tpv/Programs/methratio.py -o $path_out -d $path_genome --combine-CpG --context=CG -u $path_sam
samtools view -bS $path_sam > $id.bam
#rm $path_sam

echo "$id done!"

