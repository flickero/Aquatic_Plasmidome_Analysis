#!/bin/bash
#$ -N blast_plasmidome_PLSDB
#$ -cwd
#$ -V
#$ -l h_vmem=8G
#$ -pe ompi 4
#$ -q tals.q
#$ -o /gpfs0/tals/projects/Analysis/Oded_Project/BLAST_results/blast_2.15_plasmidome_PLSDB.o$JOB_ID
#$ -e /gpfs0/tals/projects/Analysis/Oded_Project/BLAST_results/blast_2.15_plasmidome_PLSDB.e$JOB_ID

# Activate environment
source /gpfs0/tals/projects/software/Anaconda3-2025.06/bin/activate \
    /gpfs0/tals/projects/software/Anaconda3-2025.06/envs/blast_2.15

blastn \
  -query /gpfs0/tals/projects/Analysis/Oded_Project/drep_output/all_dereplicated_plasmids.fasta \
  -db /gpfs0/tals/projects/Analysis/Oded_Project/DBs/PLSDB/plsdb \
  -task megablast \
  -evalue 1e-5 \
  -perc_identity 50 \
  -num_threads 1 \
  -outfmt "6 qseqid sseqid pident length evalue bitscore qcovs qstart qend sstart send slen stitle" \
  -out /gpfs0/tals/projects/Analysis/Oded_Project/BLAST_results/blast_2.15_plasmidome_vs_PLSDB_exploratory.tab

