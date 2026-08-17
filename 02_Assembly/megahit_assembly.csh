#!/bin/bash
#$ -o $5/Output_MegaHit
#$ -e $5/Error_MegaHit

source activate --stack /gpfs0/tals/users/androsiu/.conda/envs/megahit

megahit -t 16 -1 $1 -2 $2 --continue -o $3  # 1 paired-end library

megahit_core contig2fastg 59 $3/intermediate_contigs/k59.contigs.fa > $4/assembly_graph.fastg # get FASTG from the intermediate contigs of k=59

