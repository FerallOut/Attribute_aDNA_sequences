#!/bin/bash -l
#SBATCH -A snic2019-8-150
#SBATCH -t 22:00:00
#SBATCH -p core

# load modules
module load bioinfo-tools samtools

# run mpileup to get the necessary positional information
# samtools mpileup -f reference_file bam_file > mpileup_output_file

samtools mpileup -f /home/miba8458/private/workingFolder/raw_internal/hs37d5.fa.gz /home/miba8458/private/workingFolder/raw_internal/jar001-b1e1l1p1_CATGCTC_L004_merged.170314_ST-E00201_0191_AHH7YYALXX.hs37d5.fa.cons.90perc.bam > /home/miba8458/private/workingFolder/raw_internal/jar001.mpileup

