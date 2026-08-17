#!/bin/bash
#$ -N drep_run
#$ -cwd
#$ -j y
#$ -o drep_run.log
#$ -q tals.q
#$ -l mem_free=64G
#$ -pe ompi 8

# ===============================
#  Environment setup
# ===============================

# Limit multi-threaded BLAS libraries to prevent oversubscription
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# Point to the shared dRep environment
export CONDA_ENV_PATH="/gpfs0/tals/projects/software/Anaconda3-2025.06/envs/drep_copy"
export PATH="$CONDA_ENV_PATH/bin:$PATH"
export LD_LIBRARY_PATH="$CONDA_ENV_PATH/lib:$LD_LIBRARY_PATH"

# Temporary directory for sorting and intermediate files
export TMPDIR=/gpfs0/tals/projects/Analysis/Oded_Project/tmp

# ===============================
#  Directories and inputs
# ===============================
WORKDIR="/gpfs0/tals/projects/Analysis/Oded_Project"
INPUT_LIST="$WORKDIR/Output/All_plasmids/plasmid_list.txt"
INPUT_DIR="$WORKDIR/Output/All_plasmids"
OUTDIR="$WORKDIR/drep_output"

# ===============================
#  Step 0: Clean previous output
# ===============================
if [ -d "$OUTDIR" ]; then
    echo "🧹 Cleaning old dRep output..."
    mv "$OUTDIR" "${OUTDIR}_old_$(date +%Y%m%d_%H%M)" 2>/dev/null
fi
mkdir -p "$OUTDIR"

# ===============================
#  Step 1: Sanity check
# ===============================
list_count=$(wc -l < "$INPUT_LIST")
file_count=$(ls -1 "$INPUT_DIR"/*.fasta 2>/dev/null | wc -l)
echo "🧬 FASTA files listed: $list_count | Found: $file_count"

missing=$(grep -L "^>" "$INPUT_DIR"/*.fasta | wc -l)
if [ "$missing" -gt 0 ]; then
    echo "❌ Warning: $missing FASTA files missing headers."
    grep -L "^>" "$INPUT_DIR"/*.fasta > "$OUTDIR/missing_headers.txt"
fi

# ===============================
#  Step 2: Run dRep
# ===============================
echo "🚀 Launching dRep dereplication (shared env)..."
"$CONDA_ENV_PATH/bin/dRep" dereplicate \
    "$OUTDIR" \
    -g "$INPUT_LIST" \
    --ignoreGenomeQuality \
    --S_algorithm ANImf \
    -nc 0.5 \
    -l 1000 \
    -N50W 0 \
    -sizeW 1 \
    --clusterAlg single \
    -p 2

run_status=$?

# ===============================
#  Step 3: Check for failed genomes
# ===============================
echo "------------------------------------------"
if [ "$run_status" -ne 0 ]; then
    echo "⚠️  dRep run exited with errors (code $run_status). Checking for failed genomes..."
    grep -i "CRITICAL" drep_run.log | awk '{print $NF}' | sort | uniq > "$OUTDIR/failed_genomes.txt"
    fail_count=$(wc -l < "$OUTDIR/failed_genomes.txt")
    if [ "$fail_count" -gt 0 ]; then
        echo "💥 $fail_count failed genomes detected — attempting isolated re-clustering..."
        "$CONDA_ENV_PATH/bin/dRep" cluster \
            -g "$OUTDIR/failed_genomes.txt" \
            -p 4 \
            --S_algorithm fastANI \
            --clusterAlg single \
            --multiround_primary_clustering \
            -l 1000 \
            --ignoreGenomeQuality \
            -o "$OUTDIR/recluster_failed"
    else
        echo "✅ No specific failed genomes found — check log manually."
    fi
else
    echo "✅ dRep completed successfully."
fi

# ===============================
#  Step 4: Summary
# ===============================
if [ -f "$OUTDIR/data_tables/Cdb.csv" ]; then
    total_clusters=$(($(wc -l < "$OUTDIR/data_tables/Cdb.csv") - 1))
    echo "🏷️  Total clusters detected: $total_clusters"
else
    echo "⚠️  No Cdb.csv found — run may have failed early."
fi

echo "Run finished at $(date)"
