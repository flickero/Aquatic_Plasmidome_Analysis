#!/bin/bash
set -euo pipefail

###############################################################################
# INPUT FILES
###############################################################################

VIRALVERIFY="/gpfs0/tals/projects/Analysis/Oded_Project/viralVerify/viralverify_results/viralverify_plasmid_predictions.csv"
SCAPP="/gpfs0/tals/projects/Analysis/Oded_Project/SCAPP_validation/scapp_confirmed_plasmids.txt"
BLAST="/gpfs0/tals/projects/Analysis/Oded_Project/BLAST_results/known_plasmid_pairs_70.tsv"
INTERPRO="/gpfs0/tals/projects/Analysis/Oded_Project/interproscan/unique_plasmid_candidates_filtered.txt"

OUT="plasmid_validation_matrix.tsv"

###############################################################################
# TEMP FILES
###############################################################################

TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

VV_IDS="$TMPDIR/viralverify.ids"
SCAPP_IDS="$TMPDIR/scapp.ids"
BLAST_IDS="$TMPDIR/blast.ids"
IPS_IDS="$TMPDIR/interproscan.ids"
ALL_IDS="$TMPDIR/all.ids"

###############################################################################
# 1. EXTRACT CANONICAL PLASMID IDS (SAMPLE_RNODE_X)
###############################################################################

# ViralVerify (CSV, ID in column A)
awk -F',' 'NR>1 {
  match($1, /(.*_RNODE_[0-9]+)/, a);
  if (a[1]) print a[1];
}' "$VIRALVERIFY" | sort -u > "$VV_IDS"

# SCAPP
awk '{
  match($0, /(.*_RNODE_[0-9]+)/, a);
  if (a[1]) print a[1];
}' "$SCAPP" | sort -u > "$SCAPP_IDS"

# BLAST vs PLSDB (first field until tab)
awk -F'\t' '{
  match($1, /(.*_RNODE_[0-9]+)/, a);
  if (a[1]) print a[1];
}' "$BLAST" | sort -u > "$BLAST_IDS"

# InterProScan (already RNODE-based)
awk '{
  match($0, /(.*_RNODE_[0-9]+)/, a);
  if (a[1]) print a[1];
}' "$INTERPRO" | sort -u > "$IPS_IDS"

###############################################################################
# 2. UNION OF ALL PLASMID IDS
###############################################################################

cat "$VV_IDS" "$SCAPP_IDS" "$BLAST_IDS" "$IPS_IDS" \
  | sort -u > "$ALL_IDS"

###############################################################################
# 3. BUILD BOOLEAN VERIFICATION MATRIX
###############################################################################

{
  echo -e "plasmid_id\tViralVerify\tSCAPP\tBLAST_PLSDB\tInterProScan"

  while read -r pid; do
    vv=$(grep -Fxq "$pid" "$VV_IDS" && echo TRUE || echo FALSE)
    sc=$(grep -Fxq "$pid" "$SCAPP_IDS" && echo TRUE || echo FALSE)
    bl=$(grep -Fxq "$pid" "$BLAST_IDS" && echo TRUE || echo FALSE)
    ip=$(grep -Fxq "$pid" "$IPS_IDS" && echo TRUE || echo FALSE)

    echo -e "$pid\t$vv\t$sc\t$bl\t$ip"
  done < "$ALL_IDS"

} > "$OUT"

echo "✔ Verification table written to: $OUT"
