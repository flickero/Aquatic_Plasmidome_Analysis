# Plasmid Assembly and Dereplication

This directory contains the scripts used to reconstruct plasmid candidates from metagenomic sequencing data.

Pipeline

Raw paired-end reads

↓

MEGAHIT assembly

↓

Assembly graph generation (FASTG)

↓

SCAPP plasmid detection

↓

Aggregation of plasmid candidates

↓

dRep dereplication

Scripts

- PlasmidomeAssembly.py
- megahit_assembly.csh
- scapp.csh
- run_drep.sh

Notes

The master pipeline includes an optional MetaPlasmidSPAdes branch that was evaluated during pipeline development but was not used in the analyses presented in the thesis.