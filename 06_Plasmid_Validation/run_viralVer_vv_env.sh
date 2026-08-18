#!/bin/bash
#$ -cwd
#$ -q tals.q
#$ -V
#$ -o ./Output_viralVer_vv_env
#$ -e ./Error_viralVer_vv_env

# Activate environment
source /gpfs0/tals/projects/software/Anaconda3-2025.06/etc/profile.d/conda.sh
conda activate vv_env

# Run ViralVerify
viralverify \
  -f /gpfs0/tals/projects/Analysis/Oded_Project/drep_output/all_dereplicated_plasmids.fasta \
  -o viralverify_results \
  --hmm /gpfs0/tals/users/flickero/.conda/envs/vv_env/db/viralverify_db/viral.hmm \
  -p

