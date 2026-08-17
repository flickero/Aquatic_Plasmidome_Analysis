#!/bin/bash
#$ -cwd
#$ -q tals.q
#$ -V
#$ -N eggnog_ORFs
#$ -o eggnog/eggnog.o$JOB_ID
#$ -e eggnog/eggnog.e$JOB_ID

# Load conda inside SGE job
source /gpfs0/tals/projects/software/Anaconda3-2025.06/etc/profile.d/conda.sh

# Activate EggNOG-mapper environment
conda activate /gpfs0/tals/users/flickero/.conda/envs/emapper_env

# Input ORFs
INPUT="/gpfs0/tals/projects/Analysis/Oded_Project/prodigal/all_ORFs_clean.faa"

# Output prefix
OUT="/gpfs0/tals/projects/Analysis/Oded_Project/eggnog/all_ORFs"

# Run EggNOG-mapper
emapper.py \
    -i "$INPUT" \
    --itype proteins \
    --cpu 8 \
    --data_dir /gpfs0/tals/projects/data/Datasets/Oded_data/eggnog_db \
    --override \
    --go_evidence non-electronic \
    --target_orthologs all \
    --report_orthologs \
    --excel \
    --pfam_realign realign \
    --decorate_gff no \
    --output $(basename "$OUT") \
    --output_dir $(dirname "$OUT")

echo "EggNOG-mapper annotation completed."
