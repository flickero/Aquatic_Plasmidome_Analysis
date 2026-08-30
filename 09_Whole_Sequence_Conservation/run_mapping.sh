#!/bin/bash

source activate plasmid_mapping_env


bowtie2 -a \
  -x /gpfs0/tals/projects/Analysis/Oded_Project/test/SRR23967291/filtered_plasmids_unique_headers \
  -1 /gpfs0/tals/projects/Analysis/Oded_Project/test/SRR23967291/SRR23967291_val_1.fq.gz \
  -2 /gpfs0/tals/projects/Analysis/Oded_Project/test/SRR23967291/SRR23967291_val_2.fq.gz \
  | samtools view -bS - > /gpfs0/tals/projects/Analysis/Oded_Project/test/SRR23967291/SRR23967291_double.bam

