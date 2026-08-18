# Workflow 6 – Plasmid Validation

## Summary

This workflow integrates multiple independent validation approaches to assess the confidence of dereplicated plasmid candidates. Evidence from graph-based reconstruction, plasmid classification, conserved plasmid-associated domains, and similarity to curated reference plasmids is combined into a unified validation matrix that forms the basis for the validation statistics reported in the thesis.

---

## Description

Plasmid candidates reconstructed after dereplication were evaluated using four complementary validation strategies.

### 1. SCAPP confidence

Plasmids identified as confident cyclic assemblies by SCAPP were retained as one source of structural support.

### 2. ViralVerify classification

Dereplicated plasmids were classified using ViralVerify to distinguish plasmid sequences from viral or chromosomal elements based on sequence composition and Hidden Markov Models (HMMs).

### 3. InterProScan plasmid-associated domains

Predicted ORFs were scanned using InterProScan. Plasmids containing conserved plasmid-associated domains (e.g., replication initiation, partitioning, mobilization) were identified as receiving functional support.

### 4. PLSDB similarity search

Dereplicated plasmids were compared against the PLSDB reference plasmid database using BLASTn. High-confidence matches (≥70% reciprocal coverage, as defined in the thesis) provided homology-based support from previously reported plasmids.

Finally, results from all validation approaches were merged into a unified plasmid validation matrix indicating which independent methods supported each plasmid candidate. This matrix was subsequently used to summarize validation statistics reported in the Results and Discussion.

---

## Scripts

| Script | Purpose |
|---------|---------|
| `run_viralVer_vv_env.sh` | Runs ViralVerify on the dereplicated plasmidome. |
| `blast_2.15_plsdb.sh` | Performs BLASTn comparison of dereplicated plasmids against the PLSDB database. |
| `build_plasmid_verification_table.sh` | Integrates ViralVerify, SCAPP, InterProScan, and PLSDB evidence into a unified plasmid validation matrix. |

---

## Input

- Dereplicated plasmid FASTA (`all_dereplicated_plasmids.fasta`)
- ViralVerify HMM database
- PLSDB BLAST database
- SCAPP confident plasmid list
- InterProScan plasmid-support table

---

## Main outputs

- `viralverify_plasmid_predictions.csv`
- `known_plasmid_pairs_70.tsv`
- `unique_plasmid_candidates_filtered.txt`
- `plasmid_validation_matrix.tsv`

---

## Software

- ViralVerify
- BLAST+ (blastn)
- InterProScan
- SCAPP

---

## Notes

The validation framework intentionally combines independent sources of evidence rather than relying on a single prediction method. Structural reconstruction (SCAPP), sequence composition (ViralVerify), conserved functional domains (InterProScan), and homology to curated reference plasmids (PLSDB) provide complementary support for plasmid identity. The resulting validation matrix serves as the basis for confidence assessment throughout the thesis.