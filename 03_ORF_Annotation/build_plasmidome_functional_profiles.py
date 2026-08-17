#!/usr/bin/env python3

import pandas as pd

# ==========================================================
# INPUT FILES
# ==========================================================

PRODIGAL_FAA = "/gpfs0/tals/projects/Analysis/Oded_Project/prodigal/all_ORFs_clean.faa"
FILTERED_ORFS = "/gpfs0/tals/projects/Analysis/Oded_Project/prodigal/filtered_orf_ids.txt"

EGGNOG_FILE = "/gpfs0/tals/projects/Analysis/Oded_Project/eggnog/all_ORFs.emapper.annotations.filtered"
INTERPRO_FILE = "/gpfs0/tals/projects/Analysis/Oded_Project/prodigal/interproscan_results.filtered.tsv"
MEGARES_FILE = "/gpfs0/tals/projects/Analysis/Oded_Project/MEGAres/tblastn/ORFs_vs_MEGARes.tblastn.tsv"

OUTPUT_FILE = "plasmid_functional_profile.tsv"

# ==========================================================
# 1️⃣ LOAD FILTERED ORF LIST
# ==========================================================

print("Loading filtered ORF IDs...")

relevant_orfs = set(
    pd.read_csv(FILTERED_ORFS, header=None)[0]
    .astype(str)
    .str.strip()
)

# ==========================================================
# 2️⃣ PARSE PRODIGAL COORDINATES
# ==========================================================

print("Parsing Prodigal coordinates...")

records = []

with open(PRODIGAL_FAA) as f:
    for line in f:
        if line.startswith(">"):
            header = line.strip()[1:]
            parts = header.split("#")
            id_part = parts[0].strip()

            if id_part not in relevant_orfs:
                continue

            start = int(parts[1].strip())
            end = int(parts[2].strip())
            strand = parts[3].strip()

            plasmid_id = "_".join(id_part.split("_")[:-1])

            records.append({
                "plasmid_id": plasmid_id,
                "orf_id": id_part,
                "start": start,
                "end": end,
                "strand": strand
            })

orf_coords = pd.DataFrame(records)

# ==========================================================
# 3️⃣ PARSE EGGNOG
# ==========================================================

print("Parsing EggNOG annotations...")

with open(EGGNOG_FILE) as f:
    lines = f.readlines()

header_line_idx = None
for i, line in enumerate(lines):
    if line.startswith("#query"):
        header_line_idx = i
        break

header = lines[header_line_idx].strip().lstrip("#").split("\t")

eggnog_df = pd.read_csv(
    EGGNOG_FILE,
    sep="\t",
    skiprows=header_line_idx + 1,
    header=None
)

eggnog_df.columns = header

eggnog_df = eggnog_df[["query", "COG_category", "Description"]]
eggnog_df = eggnog_df.rename(columns={
    "query": "orf_id",
    "Description": "COG_description"
})

eggnog_df["orf_id"] = eggnog_df["orf_id"].astype(str).str.strip()
eggnog_df = eggnog_df[eggnog_df["orf_id"].isin(relevant_orfs)]

# ==========================================================
# 4️⃣ PARSE INTERPRO
# ==========================================================

print("Parsing InterPro annotations...")

interpro_cols = [
    "orf_id","md5","length","analysis","signature_accession",
    "signature_description","start","end","evalue","status",
    "date","ipr_accession","ipr_description","extra"
]

interpro_df = pd.read_csv(
    INTERPRO_FILE,
    sep="\t",
    header=None,
    names=interpro_cols
)

interpro_df["orf_id"] = interpro_df["orf_id"].astype(str).str.strip()
interpro_df = interpro_df[interpro_df["orf_id"].isin(relevant_orfs)]

interpro_df["domain"] = (
    interpro_df["signature_accession"].astype(str) + "|" +
    interpro_df["signature_description"].astype(str)
)

interpro_grouped = (
    interpro_df.groupby("orf_id")["domain"]
    .apply(lambda x: ";".join(sorted(set(x))))
    .reset_index()
    .rename(columns={"domain": "Domains"})
)

# ==========================================================
# 5️⃣ PARSE MEGARES (WHITESPACE-ROBUST)
# ==========================================================

print("Parsing MEGARes hits...")

meg_cols = [
    "orf_id","subject","duplicate_subject",
    "evalue","score","length","pident",
    "qstart","qend","sstart","send",
    "gap","bitscore","qlen","slen",
    "qseq","sseq"
]

# Use whitespace separator because file is mixed tab/space separated
meg_df = pd.read_csv(
    MEGARES_FILE,
    sep=r"\s+",
    header=None,
    engine="python"
)

# Keep only first 17 columns
meg_df = meg_df.iloc[:, :17]
meg_df.columns = meg_cols

meg_df["orf_id"] = meg_df["orf_id"].astype(str).str.strip()

# === DEBUG BLOCK ===
print("Example MEGARes ORF ID:", repr(meg_df["orf_id"].iloc[0]))
print("Example filtered ORF ID:", repr(next(iter(relevant_orfs))))
print("Is first MEGARes ORF in filtered set?",
      meg_df["orf_id"].iloc[0] in relevant_orfs)
# ===================

meg_df = meg_df[meg_df["orf_id"].isin(relevant_orfs)]

if meg_df.empty:
    print("No MEGARes hits found among filtered ORFs.")
    meg_grouped = pd.DataFrame(columns=[
        "orf_id",
        "AMR_class",
        "AMR_mechanism",
        "AMR_gene",
        "AMR_identity"
    ])
else:
    def parse_meg_subject(subject):
        parts = str(subject).split("|")
        if len(parts) >= 5:
            return pd.Series({
                "AMR_class": parts[2],
                "AMR_mechanism": parts[3],
                "AMR_gene": parts[4]
            })
        return pd.Series({
            "AMR_class": None,
            "AMR_mechanism": None,
            "AMR_gene": None
        })

    parsed = meg_df["subject"].apply(parse_meg_subject)

    meg_parsed = pd.concat(
        [meg_df.reset_index(drop=True),
         parsed.reset_index(drop=True)],
        axis=1
    )

    meg_grouped = (
        meg_parsed.groupby("orf_id")
        .agg({
            "AMR_class": lambda x: ";".join(sorted(set(x.dropna()))),
            "AMR_mechanism": lambda x: ";".join(sorted(set(x.dropna()))),
            "AMR_gene": lambda x: ";".join(sorted(set(x.dropna()))),
            "pident": "max"
        })
        .reset_index()
        .rename(columns={"pident": "AMR_identity"})
    )

# ==========================================================
# 6️⃣ MERGE DATA
# ==========================================================

print("Merging datasets...")

full_df = orf_coords.merge(eggnog_df, on="orf_id", how="left")
full_df = full_df.merge(interpro_grouped, on="orf_id", how="left")
full_df = full_df.merge(meg_grouped, on="orf_id", how="left")

full_df = full_df.sort_values(["plasmid_id", "start"])

# ==========================================================
# 7️⃣ EXPORT
# ==========================================================

print("Writing output...")
full_df.to_csv(OUTPUT_FILE, sep="\t", index=False)

print("Done.")
print("Total filtered ORFs processed:", len(full_df))
print("Unique plasmids:", full_df["plasmid_id"].nunique())