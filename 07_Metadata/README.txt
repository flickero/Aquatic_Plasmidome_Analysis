# Workflow 7 - Metadata Preparation

## Summary

This workflow standardizes metadata from the Alexander River and Bogotá River studies into harmonized tables used throughout the downstream analyses. Environmental variables, sampling locations, seasons, and sequencing statistics are reformatted into consistent metadata tables suitable for visualization and statistical analyses.

---

## Description

The two public datasets used in this study were accompanied by heterogeneous metadata formats. This workflow harmonizes both datasets into standardized metadata tables that can be consumed directly by downstream R and Python scripts.

### Alexander River

The preprocessing script:

- extracts run accession identifiers
- assigns station numbers
- converts station names into standardized identifiers
- groups stations into ecological regions
- converts sampling dates into seasonal categories
- exports the final cleaned metadata table

The resulting metadata includes:

- sample ID
- station number
- station name
- ecological region
- sampling date
- season
- sequencing depth
- water temperature

### Bogotá River

The preprocessing script:

- standardizes sample identifiers
- assigns ecological regions
- extracts sampling year
- merges antibiotic concentration measurements
- merges physicochemical measurements
- averages technical replicates where required
- exports the final cleaned metadata table

The resulting metadata includes:

- sample ID
- station
- ecological region
- sampling year
- antibiotic concentrations
- BOD
- COD
- TSS
- TN
- sequencing information

Additional metadata files containing sequencing depth and ORF counts are included for figure annotation.

---

## Scripts

| Script | Purpose |
|---------|---------|
| `block2_metadata_preprocess.R` | Generates harmonized Alexander and Bogotá metadata tables from the original supplementary spreadsheets. |

---

## Input

Original supplementary metadata tables associated with the Alexander River and Bogotá River studies.

---

## Main outputs

- `alex_clean_metadata.csv`
- `alex_clean_metadata.xlsx`
- `bogo_clean_metadata.csv`
- `bogo_clean_metadata.xlsx`
- `read_per_sample.csv`
- `ORF_count_per_sample_filtered.tsv`

---

## Downstream usage

These metadata tables are used throughout the repository, including:

- AMR heatmaps
- Functional heatmaps
- Coverage heatmaps
- PERMANOVA analyses
- PCoA analyses
- Figure annotations
- Sample ordering