# Workflow 3 – ORF Prediction and Functional Annotation

## Overview

This workflow predicts open reading frames (ORFs) from the dereplicated plasmidome and performs functional annotation using multiple complementary databases.

To avoid duplicate ORF predictions arising from circular plasmids, all plasmid sequences were first represented as doubled sequences (P+P). ORFs whose start coordinates fell within the central copy of the doubled sequence were retained for downstream analyses.

The resulting filtered protein sequences were subsequently annotated using:

- InterProScan (protein domains)
- eggNOG-mapper (orthology and COG functional categories)
- MEGARes (antimicrobial resistance determinants)

The annotation outputs were finally integrated into a unified plasmid-level functional profile used throughout the thesis.

---

## Input

Output from Workflow 2:

- `all_dereplicated_plasmids_doubled.fasta`

---

## Workflow

```
Doubled plasmid FASTA
        │
        ▼
run_prodigal_add_fastas.sh
        │
        ▼
Predicted ORFs
(.faa / .fna / .gff)
        │
        ▼
filter_centered_orfs.py
        │
        ▼
Filtered ORFs
(all_ORFs_clean.*)
        │
        ├────────► run_interproscan.sh
        │
        ├────────► run_eggnog.sh
        │
        └────────► tblastn against MEGARes
        │
        ▼
build_plasmidome_functional_profiles.py
        │
        ▼
Integrated plasmid functional profile
```

---

## Scripts

### run_prodigal_add_fastas.sh

Predicts protein-coding genes from doubled plasmid sequences using Prodigal in metagenomic mode.

Outputs:

- `all_ORFs.faa`
- `all_ORFs.fna`
- `all_ORFs.gff`

---

### filter_centered_orfs.py

Retains only ORFs whose start coordinate falls within the central copy of each doubled plasmid sequence, preventing duplicate ORF predictions introduced by sequence doubling.

---

### run_interproscan.sh

Annotates filtered protein sequences with conserved protein domains using InterProScan.

Output:

- `interproscan_results.filtered.tsv`

---

### run_eggnog.sh

Assigns orthology, COG functional categories and functional descriptions using eggNOG-mapper.

Output:

- `all_ORFs.emapper.annotations.filtered`

---

### extract_plasmid_interpro.sh

Identifies plasmid-associated InterPro domains and generates plasmid-level support tables.

Outputs include:

- `interproscan_plasmid_hits_filtered.tsv`
- `plasmid_support_counts_filtered.txt`
- `unique_plasmid_candidates_filtered.txt`

---

### build_plasmidome_functional_profiles.py

Integrates Prodigal coordinates, EggNOG annotations, InterPro domains and MEGARes annotations into a single plasmid-level functional annotation table.

Output:

- `plasmid_functional_profile.tsv`

---

### analyze_orfs_from_gff_fixed.py

Generates summary statistics describing ORF count, ORF length and coding density distributions used for Figures 4.2–4.4.

---

## Software

- Prodigal
- InterProScan 5.64-96.0
- eggNOG-mapper 2.1.13
- Python 3
- Biopython
- pandas

---

## Thesis outputs generated from this workflow

Results chapter:

- Figure 4.2 – ORF count distribution
- Figure 4.3 – ORF length distribution
- Figure 4.4 – Coding density distribution
- Figure 4.5 – Global COG functional composition
- Figure 4.6 – Sample-level functional composition
- Section 4.2.4 – Conserved domain and plasmid hallmark signatures
- Functional annotation tables used throughout the resistance analyses
