#!/usr/bin/env bash

# -------------------------
# CONFIGURATION
# -------------------------
MASTER_DIR="/gpfs0/tals/projects/Analysis/Oded_Project/test/"
OUTPUT_DIR="${MASTER_DIR}/merged_coverage_output"
mkdir -p "${OUTPUT_DIR}"

echo "Scanning ${MASTER_DIR} for plasmid coverage reports..."

# -------------------------
# STEP 1: Extract Contig + Coverage from each sample report
# -------------------------
tmp_dir="${OUTPUT_DIR}/tmp_extract"
mkdir -p "${tmp_dir}"

for file in ${MASTER_DIR}/*/plasmidome_coverage_by_*.tsv; do
    sample=$(basename "${file}" | sed 's/plasmidome_coverage_by_//' | sed 's/.tsv//')
    echo "Processing sample: ${sample}"

    # Extract Contig and Coverage columns and apply 5% noise filter
    awk -v s="${sample}" '
        BEGIN { OFS="\t" }
        NR==1 { print "Contig", s; next }
        {
            cov = $11 + 0;             # Coverage column
            if (cov < 5.0) cov = 0.0;  # Noise filter
            print $1, cov
        }
    ' "${file}" > "${tmp_dir}/${sample}.tsv"
done

echo "Finished extracting individual sample tables."

# -------------------------
# STEP 2: Build list of sample tables
# -------------------------
sample_tables=(${tmp_dir}/*.tsv)

# -------------------------
# STEP 3: Merge all tables by 'Contig'
#         (since all files share identical Contig rows in identical order)
# -------------------------
echo "Merging all tables into matrix..."

# Start with the first sample table
cp "${sample_tables[0]}" "${OUTPUT_DIR}/merged.tsv"

# Join remaining tables column-wise
for ((i=1; i<${#sample_tables[@]}; i++)); do
    paste "${OUTPUT_DIR}/merged.tsv" <(cut -f2 "${sample_tables[$i]}") \
        > "${OUTPUT_DIR}/merged.tmp"
    mv "${OUTPUT_DIR}/merged.tmp" "${OUTPUT_DIR}/merged.tsv"
done

echo "Merged matrix created."

# -------------------------
# STEP 4: Create CSV version
# -------------------------
tr '\t' ',' < "${OUTPUT_DIR}/merged.tsv" > "${OUTPUT_DIR}/merged.csv"
echo "CSV file created."

# -------------------------
# STEP 5: Create RDS loader R script
# -------------------------
cat <<EOF > "${OUTPUT_DIR}/load_matrix_in_R.R"
library(readr)

# Load TSV
mat <- read_tsv("${OUTPUT_DIR}/merged.tsv")

# Save as RDS
saveRDS(mat, "${OUTPUT_DIR}/merged.rds")

EOF

echo "R loader script created (load_matrix_in_R.R)."

echo "DONE."
