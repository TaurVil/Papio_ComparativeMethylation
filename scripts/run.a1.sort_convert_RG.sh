#!/bin/bash

module load java
module load samtools

path_genome=/data/tunglab/shared/genomes/panubis1/Panubis1.0.fa
index=${SLURM_ARRAY_TASK_ID}
id=`head -$index 00_libnames.txt | tail -1`

in_bam=$id.bam


samtools sort -o sort.$id.bam $id.bam
# rm $id.bam

java -Xmx4g -jar ~/picard.jar AddOrReplaceReadGroups VALIDATION_STRINGENCY=LENIENT INPUT=./sort.$id.bam OUTPUT=./rg.$id.bam RGPL=illumina RGLB=reseq RGPU=reseq RGSM=$id

samtools index rg.$id.bam

#rm sort.$id.bam

