#!/usr/bin/env python3
"""
analyze_plasmid_intervals.py

Scan a merged plasmid coverage matrix, identify interesting plasmids present
above a user-defined coverage threshold, and extract contiguous covered
intervals from indexed BAM files using samtools depth.

Requirements
------------
- Python 3.8+
- samtools available in PATH
- Indexed BAMs (.bam + .bai)

Usage
-----
python analyze_plasmid_intervals.py 70

Output
------
interesting_plasmid_intervals.tsv
"""

import csv
import os
import sys
import subprocess
from collections import OrderedDict

PROJECT = "/gpfs0/tals/projects/Analysis/Oded_Project"

PLASMID_FASTA = os.path.join(PROJECT, "interesting_plasmids.fasta")
MERGED = os.path.join(PROJECT, "test", "merged_coverage_output", "merged.tsv")
BAM_ROOT = os.path.join(PROJECT, "test")
OUTPUT = os.path.join(PROJECT, "interesting_plasmid_intervals.tsv")

SAMTOOLS = "/gpfs0/tals/projects/software/samtools-1.21/bin/samtools"


def read_interesting_plasmids(fasta):
    plasmids = []
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                plasmids.append(line[1:].strip())
    return plasmids


def read_coverage_table(tsv):
    coverage = OrderedDict()
    with open(tsv, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        samples = header[1:]
        for row in reader:
            plasmid = row[0]
            vals = {}
            for sample, value in zip(samples, row[1:]):
                try:
                    vals[sample] = float(value)
                except ValueError:
                    vals[sample] = 0.0
            coverage[plasmid] = vals
    return coverage


def qualifying_samples(coverage_dict, plasmid, threshold):
    if plasmid not in coverage_dict:
        return []
    out = []
    for sample, cov in coverage_dict[plasmid].items():
        if cov >= threshold:
            out.append((sample, cov))
    return out


def bam_path(sample):
    return os.path.join(
        BAM_ROOT,
        sample,
        f"{sample}_vs_drep_double.sorted.bam"
    )


def samtools_depth(contig, bam):
    cmd = [
        SAMTOOLS,
        "depth",
        "-aa",
        "-r",
        contig,
        bam
    ]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        raise RuntimeError(result.stderr.strip())

    depths = []
    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) != 3:
            continue
        pos = int(fields[1])
        depth = int(fields[2])
        depths.append((pos, depth))

    return depths


def collapse_intervals(depths):
    intervals = []

    inside = False
    start = end = None
    sum_depth = 0
    max_depth = 0
    length = 0

    for pos, depth in depths:

        if depth > 0:

            if not inside:
                inside = True
                start = pos
                end = pos
                sum_depth = 0
                max_depth = depth
                length = 0

            end = pos
            sum_depth += depth
            max_depth = max(max_depth, depth)
            length += 1

        else:

            if inside:
                intervals.append({
                    "start": start,
                    "end": end,
                    "length": length,
                    "mean_depth": sum_depth / length,
                    "max_depth": max_depth
                })
                inside = False

    if inside:
        intervals.append({
            "start": start,
            "end": end,
            "length": length,
            "mean_depth": sum_depth / length,
            "max_depth": max_depth
        })

    return intervals


def main():

    if len(sys.argv) != 2:
        print("Usage: python analyze_plasmid_intervals.py <coverage_threshold>")
        sys.exit(1)

    threshold = float(sys.argv[1])

    print(f"Coverage threshold: {threshold:.1f}%")

    plasmids = read_interesting_plasmids(PLASMID_FASTA)
    coverage = read_coverage_table(MERGED)

    with open(OUTPUT, "w", newline="") as out:

        writer = csv.writer(out, delimiter="\t")

        writer.writerow([
            "Plasmid",
            "Sample",
            "CoveragePercent",
            "Interval",
            "Start",
            "End",
            "Length",
            "MeanDepth",
            "MaxDepth"
        ])

        total = len(plasmids)

        for idx, plasmid in enumerate(plasmids, start=1):

            print("=" * 60)
            print(f"[{idx}/{total}] {plasmid}")

            qs = qualifying_samples(coverage, plasmid, threshold)

            if not qs:
                print("  No qualifying samples.")
                continue

            for sample, cov in qs:

                bam = bam_path(sample)

                if not os.path.isfile(bam):
                    print(f"  WARNING: BAM missing: {bam}")
                    continue

                print(f"  {sample} ({cov:.1f}%)")

                try:
                    depths = samtools_depth(plasmid, bam)
                except Exception as e:
                    print(f"    samtools failed: {e}")
                    continue

                intervals = collapse_intervals(depths)

                if not intervals:
                    writer.writerow([
                        plasmid,
                        sample,
                        cov,
                        0,
                        "",
                        "",
                        0,
                        0,
                        0
                    ])
                    continue

                for n, iv in enumerate(intervals, start=1):

                    writer.writerow([
                        plasmid,
                        sample,
                        cov,
                        n,
                        iv["start"],
                        iv["end"],
                        iv["length"],
                        round(iv["mean_depth"], 2),
                        iv["max_depth"]
                    ])

    print()
    print(f"Finished.")
    print(f"Output written to:\n{OUTPUT}")


if __name__ == "__main__":
    main()
