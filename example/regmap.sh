#!/bin/bash
#$ -l h_vmem=10G


module load python/x86_64/3.6.5


python3 ../bin/MINIEX_regmap.py -c INPUTS/miniexExample_cells2clusters.txt \
                                -i INPUTS/miniexExample_identities.txt \
                                -r OUTPUTS/regulons_output/miniexExample_rankedRegulons.xlsx \
                                -m INPUTS/miniexExample_matrix.txt \
                                -t 10,25,50,100,150 \
                                -o OUTPUTS/figures \
                                -d miniexExample
