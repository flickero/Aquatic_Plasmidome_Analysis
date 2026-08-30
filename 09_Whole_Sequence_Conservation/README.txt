# Workflow 09 - Whole-Sequence Conservation of Cosmopolitan Plasmids

## Summary

This workflow investigates whether cosmopolitan plasmids detected across multiple metagenomic samples are conserved throughout their entire sequence or only within specific genomic regions.

Rather than evaluating plasmid presence solely at the whole-sequence level, this workflow performs nucleotide-resolution coverage analysis by mapping metagenomic reads against doubled plasmid references, extracting contiguous covered genomic intervals, and comparing these intervals with annotated plasmid features.

The resulting visualizations reveal which functional modules are consistently conserved across environments and which regions exhibit reduced conservation.

---

## Biological objective

Whole-plasmid coverage thresholds indicate whether a plasmid is present within a sample but provide no information regarding which genomic regions are actually conserved.

Coverage intervals are generated solely from mapped reads across the complete plasmid sequence. Functional annotations are added only after interval extraction to facilitate biological interpretation without influencing interval detection.

---

## Pipeline overview

```
Representative cosmopolitan plasmids
                │
                ▼
Bowtie2 mapping against doubled plasmid references
                │
                ▼
Sorted BAM files
                │
                ▼
samtools depth (-aa)
                │
                ▼
Per-base coverage
                │
                ▼
Coverage interval extraction
                │
                ▼
coverage_intervals.tsv
                │
                ▼
Feature overlap analysis
                │
                ▼
feature_summary.tsv
                │
                ▼
Publication-quality conservation plots
```

---

## Workflow steps

### 1. Read mapping

Metagenomic reads are aligned against doubled plasmid reference sequences using Bowtie2.

Doubling the plasmid sequence prevents artificial fragmentation of genes and coverage intervals caused by the arbitrary linearization point of circular plasmids.

---

### 2. Coverage calculation

Per-base sequencing depth is extracted from each BAM file using

```
samtools depth -aa
```

The `-aa` option reports coverage for every nucleotide position, including positions with zero coverage.

---

### 3. Coverage interval extraction

Continuous regions with sequencing depth greater than zero are merged into genomic coverage intervals.

Coordinates are converted from the doubled reference back into the original plasmid coordinate system, yielding complete whole-plasmid coverage maps.

Output:

```
coverage_intervals.tsv
```

---

### 4. Feature overlap

Coverage intervals are intersected with annotated plasmid features, including:

- replication genes
- maintenance genes
- mobility genes
- recombination genes
- toxin-antitoxin systems
- resistance genes
- OriV regions
- hypothetical proteins

Importantly, annotations are **not** used to define conserved regions.

They are overlaid only after coverage intervals have been identified.

Output:

```
feature_summary.tsv
```

---

### 5. Visualization

Coverage intervals and annotated plasmid features are combined to generate plasmid conservation maps.

Each figure displays:

- complete plasmid architecture
- functional annotation
- sample-specific conserved regions
- conservation patterns across multiple environments

---

## Directory structure

```
01_Mapping/
02_Coverage/
03_Interval_Extraction/
04_Feature_Overlap/
05_Visualization/
06_Input_Examples/
07_Output_Examples/
```

---

## Input

- Representative plasmid FASTA sequences
- Sorted BAM alignment files
- Curated plasmid feature annotations
- Feature label dictionary

---

## Output

- `coverage_intervals.tsv`
- `feature_summary.tsv`
- Publication-quality whole-sequence conservation figures

---

## Software

- Bowtie2
- SAMtools
- Python 3
- pandas
- R
- ggplot2

---

## Main scripts

| Script | Purpose |
|---------|----------|
| `run_mapping.sh` | Maps reads against doubled plasmid references |
| `coverage_analysis_lucy.sh` | Calculates nucleotide-level coverage |
| `plasmid_coverage_matrix_generator.sh` | Generates coverage matrices |
| `analyze_plasmid_intervals.py` | Converts per-base coverage into genomic intervals |
| `build_feature_summary.py` | Computes overlap between intervals and annotated features |
| `plot_cosmopolitan_plasmids.R` | Produces publication-quality conservation figures |

---

## Notes

Coverage intervals are derived exclusively from nucleotide-level sequencing coverage across the complete plasmid sequence.

Annotated genes and functional elements are incorporated only during the interpretation stage and are **not** used to define conserved regions. This ensures that the conservation analysis remains independent of prior biological annotation while enabling direct visualization of conserved functional modules.