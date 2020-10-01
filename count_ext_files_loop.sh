#!/bin/bash -l

module load bioinfo-tools samtools

# run mpileup to get the necessary positional information  
# for loop goes through all the '.bam' files of the folder

extension=bam
source_in=/proj/snic2019-8-150/nobackup/private/Mirela/neanderthal_project/data/raw_external
source_out=/home/miba8458/private/workingFolder/output
file=${source_in}/*.${extension}


for i in $file
do
filename=$(basename -a $i)
if [[ $filename = *"-"* ]] 
then
   samtools view -u -s 83.10 $i | samtools mpileup -f initial_files/hs37d5.fa.gz - > $source_out/"${filename%-*}".mpileup 
else
   samtools view -u -s 83.10 $i | samtools mpileup -f initial_files/hs37d5.fa.gz - > $source_out/${filename%%.*}.mpileup
fi
done


# run the python script to count the Neanderthal vs hu genomic positions

module load python3
python3 -u count_ND_ext.py > ext_NDresults    # unbuffered log

