|Date  	|Activity 	|What was learned 	|Trigger 	|
|---	|---	|---	|---	|
|Tue, 2020.01.21	|* getting the genomic locations known to be different in modern human and in Neanderthals > needed for counting these locations in newly sequenced aDNA   	| * If the data is too large and can't be put in the supplementary data section of the paper, they are linked in the paper somewhere towards the end. They are not necessarily in an establishe database and accessed through an accession nr. They can just be stored on the website of the institution they work at.  	| * I didn't find them, Torsten showed me that they are usually in a link at the end of the paper.  |
|   	|* basic coding to look at the structure of a file: <br> - unique values in a column<br> -how many times the unique value appears|* **Always verify what you get from your code by getting to the same outcome using at least 2 different approaches**  <br> * The bash commands I used to find out some general information about the files I was looking at, were giving different results. `sort \| uniq` works a bit different than I realised. I chose an `awk` code.| * I wasted a few hours because I was too tired to think straight.|
||* understanding how the numerotation of the positions in the chromosomes correspond to the position numbers given in the sequencing files|* yes, they correspond|* need to know what to extract|
|Thu, 2020.01.23|* learning to use `mpileup` from `samtools` to extract the DNA information for each location in the chromosome. Later will query only the necessary locations, described in the file from the paper.|   	|   	|
|Fri, 2020.01.24|* coding|* I understood that I don't need either `bcftools` or a .vcf file. A .bam file passed through `mpileup` is enough|* I checked what I get after using all these programs to work on the .bam files that I received|
|   |   	|* finaly made it to the actual coding task ':D   	|   	|
|Mon, 2020.01.27|* sending the jobs to slurm - had a problem|* keep your eyes on the syntax, damn it!!|* help from Torsten|


