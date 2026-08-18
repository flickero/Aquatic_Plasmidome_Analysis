#!/usr/bin/env python3

import pandas as pd
import re

# ==============================
# INPUT / OUTPUT
# ==============================

INPUT_TSV = "/gpfs0/tals/projects/Analysis/Oded_Project/MEGAres/tblastn/ORFs_vs_MEGARes.tblastn.tsv"
FILTERED_ORFS = "/gpfs0/tals/projects/Analysis/Oded_Project/prodigal/filtered_orf_ids.txt"
OUTPUT_TSV = "/gpfs0/tals/projects/Analysis/Oded_Project/MEGAres/tblastn/plasmid_AMR_target_presence.tsv"

# ==============================
# COLUMN NAMES
# ==============================

columns = [
    "qseqid", "sseqid", "stitle", "evalue", "length", "nident", "pident",
    "mismatch", "score", "qcovs", "qstart", "qend",
    "sstart", "send", "sframe", "slen", "qseq", "sseq"
]

# ==============================
# LOAD DATA
# ==============================

df = pd.read_csv(INPUT_TSV, sep="\t", header=None, names=columns)

# Filter ORFs
with open(FILTERED_ORFS) as f:
    valid_orfs = set(line.strip() for line in f if line.strip())

df = df[df["qseqid"].isin(valid_orfs)]

# Extract plasmid ID
def extract_plasmid_id(orf_id):
    return re.sub(r"_\d+$", "", orf_id)

df["plasmid_id"] = df["qseqid"].apply(extract_plasmid_id)

def extract_resistance_target(stitle):
    parts = stitle.split("|")
    return parts[2] if len(parts) >= 3 else "UNKNOWN"

df["resistance_target"] = df["stitle"].apply(extract_resistance_target)

# Build plasmid × resistance target presence matrix
df["presence"] = 1

presence_matrix = (
    df
    .groupby(["plasmid_id", "resistance_target"])["presence"]
    .max()
    .unstack(fill_value=0)
    .sort_index()
)

presence_matrix.to_csv(OUTPUT_TSV, sep="\t")

print("✔ Plasmid × resistance_target presence table written.")
print(f"Plasmids: {presence_matrix.shape[0]}")
print(f"Resistance targets: {presence_matrix.shape[1]}")
