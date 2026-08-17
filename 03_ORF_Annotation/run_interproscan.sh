#!/bin/bash
#$ -cwd
#$ -q tals.q
#$ -V
#$ -o ./Output_interproscan_tsv
#$ -e ./Error_interproscan_tsv

/gpfs0/tals/projects/software/interproscan/interproscan-5.64-96.0/interproscan.sh \
  -i /gpfs0/tals/projects/Analysis/Oded_Project/prodigal/all_ORFs_clean.faa \
  -f tsv \
  -iprlookup \
  -goterms \
  --disable-precalc \
  --tempdir /gpfs0/tals/projects/Analysis/Oded_Project/interproscan_temp_tsv \
  -cpu 8 \
  -o /gpfs0/tals/projects/Analysis/Oded_Project/prodigal/interproscan_results_tsv
