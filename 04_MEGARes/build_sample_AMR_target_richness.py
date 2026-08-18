#!/usr/bin/env python3

import pandas as pd

# ==============================
# INPUT FILES
# ==============================

PLASMID_TARGET_TSV = (
    "/gpfs0/tals/projects/Analysis/Oded_Project/MEGAres/tblastn/"
    "plasmid_AMR_target_presence.tsv"
)

MERGED_TSV = (
    "/gpfs0/tals/projects/Analysis/Oded_Project/test/"
    "merged_coverage_output/merged.tsv"
)

OUTPUT_TSV = (
    "/gpfs0/tals/projects/Analysis/Oded_Project/MEGAres/tblastn/"
    "sample_AMR_target_richness.tsv"
)

COVERAGE_THRESHOLD = 70.0

# ==============================
# 1. LOAD PLASMID × TARGET TABLE
# ==============================

plasmid_target = pd.read_csv(PLASMID_TARGET_TSV, sep="\t", index_col=0)

plasmid_target_long = (
    plasmid_target
    .reset_index()
    .melt(
        id_vars="plasmid_id",
        var_name="resistance_target",
        value_name="presence"
    )
)

plasmid_target_long = plasmid_target_long[
    plasmid_target_long["presence"] == 1
]

# ==============================
# 2. LOAD MERGED COVERAGE MATRIX
# ==============================

merged = pd.read_csv(MERGED_TSV, sep="\t")

if "Contig" not in merged.columns:
    raise ValueError("merged.tsv must contain a 'Contig' column")

coverage_long = merged.melt(
    id_vars="Contig",
    var_name="sample_id",
    value_name="coverage"
).rename(columns={"Contig": "plasmid_id"})

coverage_present = coverage_long[
    coverage_long["coverage"] >= COVERAGE_THRESHOLD
]

# ==============================
# 3. JOIN SAMPLE ← PLASMID ← TARGET
# ==============================

sample_target = (
    coverage_present
    .merge(plasmid_target_long, on="plasmid_id", how="inner")
    [["sample_id", "plasmid_id", "resistance_target"]]
)

# ==============================
# 4. COMPUTE TARGET RICHNESS
# ==============================

sample_target_richness = (
    sample_target
    .groupby(["sample_id", "resistance_target"])["plasmid_id"]
    .nunique()
    .reset_index(name="richness")
)

sample_target_matrix = (
    sample_target_richness
    .pivot_table(
        index="sample_id",
        columns="resistance_target",
        values="richness",
        fill_value=0
    )
    .sort_index()
)

sample_target_matrix.to_csv(OUTPUT_TSV, sep="\t")

print("✔ Sample × resistance_target richness table written.")
print(f"Samples: {sample_target_matrix.shape[0]}")
print(f"Resistance targets: {sample_target_matrix.shape[1]}")
