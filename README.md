# Aquatic Plasmidome Analysis Pipeline

## Repository accompanying MSc Thesis

**Author:** Oded Flicker

**Institution:** Ben-Gurion University of the Negev

**Degree:** M.Sc. in Marine Biotechnology (Bioinformatics)

---

## Overview

This repository accompanies the MSc thesis:

> *Identification and Characterization of Aquatic Plasmids from Variable Environmental Conditions*

The repository contains the complete computational workflow used to reconstruct, validate, annotate and analyze plasmids recovered from aquatic metagenomes originating from two independent environmental systems:

- Alexander River Estuary (Israel)
- Bogotá River Basin (Colombia)

The analysis integrates plasmid assembly, dereplication, functional annotation, antimicrobial resistance profiling, ecological analyses and figure generation.

---

# Repository structure

```
01_Assembly/
02_Dereplication/
03_ORF_Annotation/
04_AMR_Analysis/
05_Coverage_Analysis/
06_Plasmid_Validation/
07_Metadata/
08_Figure_Scripts/
supplementary_data/
```

Each workflow contains:

- README describing the workflow
- scripts used in the thesis
- required input files
- expected outputs

---

# Analysis workflow

The complete analysis proceeds as follows:

```
Raw metagenomic reads
        │
        ▼
Assembly (MEGAHIT + SCAPP)
        │
        ▼
Dereplication (dRep)
        │
        ▼
ORF prediction (Prodigal)
        │
        ▼
Functional annotation
(EggNOG + InterProScan)
        │
        ▼
AMR annotation (MEGARes)
        │
        ▼
Coverage mapping
(Bowtie2 + SAMtools)
        │
        ▼
Ecological analyses
        │
        ▼
Figures and statistics
```

---

# Software

Major software used in this study

| Software | Version |
|-----------|----------|
| MEGAHIT | (see manuscript) |
| SCAPP | (see manuscript) |
| dRep | (see manuscript) |
| Prodigal | (see manuscript) |
| EggNOG-mapper | 2.1.13 |
| InterProScan | 5.64-96.0 |
| BLAST+ | 2.15.0 |
| Bowtie2 | (see manuscript) |
| SAMtools | (see manuscript) |
| ViralVerify | (see manuscript) |
| R | 4.x |

Complete software versions are reported in the thesis Methods section.

---

# Repository contents

The repository contains the scripts required to reproduce the analyses presented in the thesis.

Large intermediate datasets, raw sequencing data and reference databases are **not** included.

Examples include:

- MEGARes database
- PLSDB database
- EggNOG database
- InterProScan databases
- Raw FASTQ files

These should be downloaded separately.

---

# Supplementary datasets

Machine-readable supplementary datasets are provided under

```
supplementary_data/
```

including for example

- plasmid validation matrix
- functional ORF annotations
- sample AMR richness matrix
- cleaned metadata tables

---

# Reproducibility

All scripts are provided exactly as used during the thesis unless otherwise noted.

Several HPC scripts contain hard-coded project paths corresponding to the Ben-Gurion University HPC environment. Users will need to adapt file paths and scheduler directives for their own systems.

---

# Citation

If you use this repository, please cite the accompanying MSc thesis.

---

# Contact

Oded Flicker

Department of Marine Biotechnology

Ben-Gurion University of the Negev

Israel
