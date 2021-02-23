#!/bin/bash
#SBATCH --get-user-env
module load samtools
module load sratoolkit

index=${SLURM_ARRAY_TASK_ID}
SRR=`head -$index 00_SRAfiles.txt | tail -1`
nom=`head -$index 00_libnames.txt | tail -1`

fastq-dump  $SRR --gzip

mv $SRR.fastq.gz $nom.fq.gz

#--split-files used for paired-end data
