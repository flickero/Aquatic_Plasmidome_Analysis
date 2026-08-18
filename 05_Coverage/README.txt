# Workflow 5 - Ecological Distribution of Dereplicated Plasmids

## Summary

This workflow characterizes the ecological distribution of dereplicated plasmids across all metagenomic samples. Reads from each sample are mapped to the dereplicated plasmid catalog, coverage is calculated for every plasmid, and individual coverage reports are merged into a global plasmid × sample coverage matrix. The resulting matrix is used for ecological visualization and statistical analyses, including heatmaps, hierarchical clustering, principal coordinate analysis (PCoA), and PERMANOVA.

---

## Description

The workflow consists of five major steps:

1. **Read mapping**
   - Paired-end metagenomic reads are aligned against the dereplicated plasmid reference using Bowtie2.
   - Alignments are converted into sorted and indexed BAM files with SAMtools.
   - The repository contains an example mapping script; the same procedure was executed independently for every sample.

2. **Coverage calculation**
   - Coverage statistics are calculated for every plasmid using the central 50% of each doubled plasmid sequence to avoid circular edge artifacts.
   - Per-sample reports include coverage percentage, depth statistics, and mapping summaries.

3. **Coverage matrix generation**
   - Individual coverage reports are merged into a single plasmid × sample coverage matrix (`merged.tsv`).
   - Coverage values below 5% are treated as background noise and set to zero.

4. **Ecological visualization**
   - Coverage matrices are used to generate publication-quality heatmaps for the Alexander River and Global datasets.
   - Plasmids are clustered using Pearson correlation with UPGMA hierarchical clustering.
   - Sample annotations are incorporated from study metadata.

5. **Community statistical analyses**
   - Plasmid presence is defined using a ≥70% coverage threshold.
   - Presence–absence matrices are analyzed using Jaccard distances.
   - Community structure is evaluated using PCoA, PERMANOVA, and homogeneity-of-dispersion tests (betadisper).
   - The Global analysis additionally includes pairwise PERMANOVA comparisons.

---

## Scripts

| Script | Purpose |
|---------|---------|
| `run_mapping.sh` | Example Bowtie2/SAMtools mapping pipeline for one sample. |
| `coverage_analysis_lucy.sh` | Calculates plasmid coverage statistics from BAM files. |
| `plasmid_coverage_matrix_generator.sh` | Merges individual coverage reports into a global coverage matrix. |
| `Alexander_stage2_environmental_column_order_w_numbering_season.R` | Generates the final Alexander River coverage heatmap. |
| `Alexander_PERMANOVA.R` | Community ecology analyses for the Alexander dataset. |
| `Bogota_heatmap_statistical_analysis.R` | Community ecology and environmental analyses for the Bogotá dataset. |
| `Global_new_heatmap_generator_annotation_added.R` | Generates the final Global plasmid coverage heatmap. |
| `Global_community_analysis.R` | Global community ecology analyses, including PCoA and PERMANOVA. |

---

## Main outputs

- `merged.tsv`
- Coverage heatmaps
- PCoA coordinates and figures
- PERMANOVA reports
- Community statistics