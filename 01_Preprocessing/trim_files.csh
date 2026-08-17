#!/bin/bash
#$ -o ../logs/${4}_Output_Trim
#$ -e ../logs/${4}_Error_Trim

source activate --stack /gpfs0/tals/projects/software/anaconda3/envs/cutadaptenv

sample=$1

# Define the file patterns
# trim galore output
file_pattern1="${4}_val_1.fq.gz"
file_pattern2="${4}_val_2.fq.gz"

#input files
in1="$sample/$file_pattern1"
in2="$sample/$file_pattern2"

# Check if both files exist and are non-empty
if [ -s "$in1" ] && [ -s "$in2" ]; then
    echo "Reads are trimmed"
else
    echo "One or both input files do not exist or are empty. Running trim_galore."

    # Run trim_galore
    trim_galore -j 2 --illumina --paired --clip_R1 5 --clip_R2 5 --three_prime_clip_R1 5 --three_prime_clip_R2 5 --basename "$4" -o $sample $2 $3
    
    # Update input file paths
    in1="$sample/$file_pattern1"
    in2="$sample/$file_pattern2"
fi
