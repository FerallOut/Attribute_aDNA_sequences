#!/bin/bash -l

module load bioinfo-tools samtools

# run mpileup to get the necessary positional information  
samtools view -u -s 83.10 /home/miba8458/LnktoWorkingFolder/UppmaxDataLnk/data/raw_external/I1300.hs37d5.fa.cons.90perc.bam| samtools mpileup -f initial_files/hs37d5.fa.gz - > output/I1300.mpileup
samtools view -u -s 83.10 /home/miba8458/LnktoWorkingFolder/UppmaxDataLnk/data/raw_external/I1302.hs37d5.fa.cons.90perc.bam| samtools mpileup -f initial_files/hs37d5.fa.gz - > output/I1302.mpileup
samtools view -u -s 83.10 /home/miba8458/LnktoWorkingFolder/UppmaxDataLnk/data/raw_external/I1303.hs37d5.fa.cons.90perc.bam| samtools mpileup -f initial_files/hs37d5.fa.gz - > output/I1303.mpileup
samtools view -u -s 83.10 /home/miba8458/LnktoWorkingFolder/UppmaxDataLnk/data/raw_external/I1314.hs37d5.fa.cons.90perc.bam| samtools mpileup -f initial_files/hs37d5.fa.gz -  > output/I1324.mpileup
samtools view -u -s 83.10 /home/miba8458/LnktoWorkingFolder/UppmaxDataLnk/data/raw_external/Vi33.16-alllib.final.bam| samtools mpileup -f initial_files/hs37d5.fa.gz -  > output/Vi33.16.mpileup
samtools view -u -s 83.10 /home/miba8458/LnktoWorkingFolder/UppmaxDataLnk/data/raw_external/Vi33.25-alllib.final.bam| samtools mpileup -f initial_files/hs37d5.fa.gz -  > output/Vi33.25.mpileup
samtools view -u -s 83.10 /home/miba8458/LnktoWorkingFolder/UppmaxDataLnk/data/raw_external/Vi33.26-alllib.final.bam| samtools mpileup -f initial_files/hs37d5.fa.gz -  > output/Vi33.26.mpileup


# run the python script to count the Neanderthal vs hu genomic positions

module load python3
python3 -u count_ND_ext.py > ext_NDresults    # unbuffered log
