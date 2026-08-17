#!/bin/bash
#$ -o $5/Output_SCAPP
#$ -e $5/Error_SCAPP

source activate --stack /gpfs0/tals/users/androsiu/.conda/envs/scapp

scapp -g $1 -o $2 -k 59 -r1 $3 -r2 $4

