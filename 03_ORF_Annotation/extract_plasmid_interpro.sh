#!/bin/bash
#$ -cwd
#$ -q tals.q
#$ -V

###############################################################################
# PURPOSE
# Identify plasmid candidates supported by InterProScan functional annotation
# by detecting "plasmid" terms in InterPro description (column M),
# collapsing annotation hits → ORFs → plasmids (RNODE-level).
###############################################################################

set -euo pipefail

###############################################################################
# PATHS
###############################################################################

INPUT="/gpfs0/tals/projects/Analysis/Oded_Project/prodigal/interproscan_results.filtered.tsv"
OUTDIR="/gpfs0/tals/projects/Analysis/Oded_Project/interproscan/"

mkdir -p "$OUTDIR"

PLASMID_HITS="$OUTDIR/interproscan_plasmid_hits_filtered.tsv"
PLASMID_COUNTS="$OUTDIR/plasmid_support_counts_filtered.txt"
UNIQUE_PLASMIDS="$OUTDIR/unique_plasmid_candidates_filtered.txt"

###############################################################################
# SANITY CHECK
###############################################################################

if [[ ! -f "$INPUT" ]]; then
    echo "ERROR: Input InterProScan TSV not found:"
    echo "  $INPUT"
    exit 1
fi

###############################################################################
# STEP 1
# Extract InterProScan entries whose InterPro description (column M, field 13)
# contains the word "plasmid" (case-insensitive).
###############################################################################

awk -F'\t' '
BEGIN { OFS="\t" }
NR==1 { print; next }
tolower($13) ~ /plasmid/ { print }
' "$INPUT" > "$PLASMID_HITS"

###############################################################################
# STEP 2
# Deduplicate at ORF level, collapse to plasmid level (RNODE),
# and count UNIQUE ORFs per plasmid.
###############################################################################

cut -f1 "$PLASMID_HITS" \
| sort -u \
| awk '
{
    if (match($0, /^[^_]+_RNODE_[0-9]+/, a)) {
        print a[0]
    }
}
' \
| sort \
| uniq -c \
| sort -nr > "$PLASMID_COUNTS"

###############################################################################
# STEP 3
# Extract clean plasmid-level candidate list
###############################################################################

awk '{print $2}' "$PLASMID_COUNTS" > "$UNIQUE_PLASMIDS"

###############################################################################
# STEP 4
# FINAL VALIDATION CHECKS
###############################################################################

ORF_COUNT=$(cut -f1 "$PLASMID_HITS" | sort -u | wc -l)
SUM_COUNTS=$(awk '{sum += $1} END {print sum}' "$PLASMID_COUNTS")
PLASMID_COUNT=$(wc -l < "$UNIQUE_PLASMIDS")

echo "============================================================"
echo "InterProScan plasmid validation summary"
echo "============================================================"
echo "Unique ORFs with plasmid-associated annotation : $ORF_COUNT"
echo "Summed ORF counts across plasmids             : $SUM_COUNTS"
echo "Unique plasmid candidates (RNODE-level)       : $PLASMID_COUNT"
echo
echo "Outputs:"
echo "  $PLASMID_HITS"
echo "  $PLASMID_COUNTS"
echo "  $UNIQUE_PLASMIDS"
echo "============================================================"

###############################################################################
# CONSISTENCY GUARANTEE
###############################################################################

if [[ "$ORF_COUNT" -ne "$SUM_COUNTS" ]]; then
    echo "WARNING: ORF count mismatch — investigate immediately"
    exit 2
fi
