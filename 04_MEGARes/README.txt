# Workflow 4 - Plasmid-associated resistance profiling

Purpose
-------
Identify plasmid-encoded resistance determinants and quantify their ecological
distribution across environmental samples.

Pipeline

Filtered ORFs
    │
    ▼
tblastn against MEGARes v3.0
    │
    ▼
ORFs_vs_MEGARes.tblastn.tsv
    │
    ▼
build_plasmid_AMR_resistance_targets_presence.py
    │
    ▼
plasmid_AMR_target_presence.tsv
    │
    ▼
build_sample_AMR_target_richness.py
    │
    ▼
sample_AMR_target_richness.tsv
    │
    ├──► Mres_Final.R
    │         Figure 4.7
    │
    └──► AMR_statistical_analysis.R
              PERMANOVA
              Pairwise PERMANOVA
              Supplementary statistics