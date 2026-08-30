#!/bin/bash
#$ -cwd
#$ -q tals.q
#$ -V
#$ -o /gpfs0/tals/projects/Analysis/Oded_Project/test/SRR23967296/Output_coverage
#$ -e /gpfs0/tals/projects/Analysis/Oded_Project/test/SRR23967296/Error_coverage

BAM="$1"
OUTFILE="$2"

LENGTHS="/gpfs0/tals/projects/Analysis/Oded_Project/test/SRR23967291/lengths_double_all.txt"

# Output header
echo -e "Contig\tSingle_Length\tStart_Pos\tEnd_Pos\tMean_Depth\tMedian_Depth\tMin_Depth\tMax_Depth\tNum_Pos_≥Min\t%_Pos_≥Min\tCoverage" > "$OUTFILE"

if [ ! -f "${BAM}.bai" ]; then
    echo "Indexing BAM..."
    samtools index -@ 8 "$BAM"
fi

# Process each contig
tail -n +2 "$LENGTHS" | while IFS=$'\t' read -r contig length; do
    echo "Processing $contig (length $length)"

    half=$((length / 2))
    quarter=$((half / 2))
    start=$((half - quarter))
    end=$((half + quarter))
    region="${contig}:${start}-${end}"

    echo "  ROI: $start to $end"

    # Get depth values for the region
    depths=$(samtools depth -aa -g SECONDARY -r "$region" "$BAM" | awk '{ print $3 }')

    if [ -z "$depths" ]; then
        echo "Warning: No depth found for $contig in region $region"
        echo -e "${contig}\t${half}\t${start}\t${end}\t0\t0\t0\t0\t0\t0\t0" >> "$OUTFILE"
        continue
    fi

    mapfile -t depth_array <<< "$depths"
    total_positions=${#depth_array[@]}

    stats=$(printf "%s\n" "${depth_array[@]}" | awk -v n="$total_positions" '
    {
        depth[NR] = $1
        sum += $1
        if ($1 > max) max = $1
        if (NR == 1 || $1 < min) min = $1
        if ($1 >= 1) covered++
    }
    END {
        # Sort for median
        npos = asort(depth)
        if (n % 2 == 1) {
            median = depth[(n + 1) / 2]
        } else {
            median = (depth[n / 2] + depth[n / 2 + 1]) / 2
        }

        # Count positions ≥ min
        count_ge_min = 0
        for (i = 1; i <= n; i++) {
            if (depth[i] >= min) count_ge_min++
        }

        mean = sum / n
        perc_ge_min = (count_ge_min * 100.0) / n
        coverage = (covered * 100.0) / n

        printf "%.2f\t%.2f\t%d\t%d\t%d\t%.2f\t%.2f", mean, median, min, max, count_ge_min, perc_ge_min, coverage
    }')

    echo -e "${contig}\t${half}\t${start}\t${end}\t${stats}" >> "$OUTFILE"
done
