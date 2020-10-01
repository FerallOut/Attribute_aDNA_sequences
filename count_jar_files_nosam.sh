#!/bin/bash -l


# run the python script to count the Neanderthal vs hu genomic positions

module load python3
python3 -u count_ND_jar.py > jar_NDresults
