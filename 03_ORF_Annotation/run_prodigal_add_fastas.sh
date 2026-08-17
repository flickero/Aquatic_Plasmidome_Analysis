#!/bin/bash
#$ -N prodigal_run
#$ -q tals.q
#$ -cwd
#$ -V
#$ -o prodigal_run.out
#$ -e prodigal_run.err
#$ -l h_rt=01:00:00
#$ -l h_vmem=4G

# Load conda and activate the environment with Prodigal
source ~/.bashrc
conda activate prodigal_env

prodigal \
  -i /gpfs0/tals/projects/Analysis/Oded_Project/drep_output/all_dereplicated_plasmids_doubled.fasta \
  -o all_ORFs.gff \
  -a all_ORFs.faa \
  -d all_ORFs.fna \
  -p meta

