#!/bin/bash
#$ -cwd
#$ -q tals.q
#$ -V
#$ -o /gpfs0/tals/projects/Analysis/Oded_Project/test/SRR23967295/Output_mapping
#$ -e /gpfs0/tals/projects/Analysis/Oded_Project/test/SRR23967295/Error_mapping

source activate plasmid_mapping_env

bowtie2 -a \
  -x /gpfs0/tals/projects/Analysis/Oded_Project/test/SRR23967291/drep_plasmids_index \
  -1 /gpfs0/tals/projects/Analysis/Oded_Project/Alexander_est/reads/SAMN33870684/SRR23967295_val_1.fq.gz \
  -2 /gpfs0/tals/projects/Analysis/Oded_Project/Alexander_est/reads/SAMN33870684/SRR23967295_val_2.fq.gz \
  | samtools view -bS - > /gpfs0/tals/projects/Analysis/Oded_Project/test/SRR23967295/SRR23967295_vs_drep_double.bam
samtools sort /gpfs0/tals/projects/Analysis/Oded_Project/test/SRR23967295/SRR23967295_vs_drep_double.bam -o /gpfs0/tals/projects/Analysis/Oded_Project/test/SRR23967295/SRR23967295_vs_drep_double.sorted.bam
samtools index SRR23967295_vs_drep_double.sorted.bam