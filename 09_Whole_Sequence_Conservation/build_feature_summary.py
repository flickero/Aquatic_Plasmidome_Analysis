#!/usr/bin/env python3

"""
build_feature_summary.py

Purpose
-------
Summarize how well each annotated feature is represented in each sample.

Inputs
------
coverage_intervals.tsv
master_annotations.tsv

Output
------
feature_summary.tsv

Coordinates are assumed to already be in the ORIGINAL plasmid coordinate
system.
"""

import csv
import os

###############################################################################
# Project paths
###############################################################################

PROJECT = "/gpfs0/tals/projects/Analysis/Oded_Project"

INTERVALS = os.path.join(
    PROJECT,
    "plasmid_interval_analysis",
    "coverage",
    "coverage_intervals.tsv"
)

ANNOTATIONS = os.path.join(
    PROJECT,
    "master_annotations.tsv"
)

OUTPUT = os.path.join(
    PROJECT,
    "feature_summary.tsv"
)

###############################################################################
# Read coverage intervals
###############################################################################

def read_intervals():

    intervals = []

    with open(INTERVALS) as fh:

        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )

        for row in reader:

            intervals.append({

                "plasmid": row["Plasmid"],

                "sample": row["Sample"],

                "coverage": float(row["CoveragePercent"]),

                "start": int(row["Original_Start"]),

                "end": int(row["Original_End"])

            })

    return intervals


###############################################################################
# Read annotations
###############################################################################

def read_annotations():

    annotations = []

    with open(ANNOTATIONS) as fh:

        reader = csv.DictReader(
            fh,
            delimiter="\t"
        )

        for row in reader:

            annotations.append({

                "plasmid": row["Plasmid"],

                "feature": row["Feature"],

                "type": row["Type"],

                "source": row["Source"],

                "start": int(row["Start"]),

                "end": int(row["End"]),

                "strand": row["Strand"],

                "COG_category": row["COG_category"],

                "COG_description": row["COG_description"],

                "Domains": row["Domains"],

                "AMR_class": row["AMR_class"],

                "AMR_mechanism": row["AMR_mechanism"],

                "AMR_gene": row["AMR_gene"],

                "AMR_identity": row["AMR_identity"]

            })

    return annotations


###############################################################################
# Calculate overlap
###############################################################################

def overlap(start1, end1,
            start2, end2):

    start = max(start1, start2)

    end = min(end1, end2)

    if start > end:

        return None

    return (

        start,

        end,

        end - start + 1

    )


###############################################################################
# Merge overlapping segments
###############################################################################

def merge_segments(segments):
    """
    Merge overlapping or adjacent overlap segments.

    Example

        100-200
        180-250
        251-300

    becomes

        100-300
    """

    if len(segments) == 0:

        return []

    segments.sort()

    merged = [segments[0]]

    for start, end in segments[1:]:

        last_start, last_end = merged[-1]

        if start <= last_end + 1:

            merged[-1] = (

                last_start,

                max(last_end, end)

            )

        else:

            merged.append(

                (start, end)

            )

    return merged
###############################################################################
# Create display name for plotting
###############################################################################

def get_display_name(feature):

    # OriV features
    if feature["type"] == "OriV":
        return feature["feature"]

    # AMR genes have highest priority
    if feature["AMR_gene"].strip():
        return feature["AMR_gene"]

    # Then COG description
    if feature["COG_description"].strip():

        desc = feature["COG_description"].strip()

        if desc != "-":
            return desc

    # Then first InterPro/Pfam domain
    if feature["Domains"].strip():

        first = feature["Domains"].split(";")[0]

        if "|" in first:

            return first.split("|", 1)[1]

        return first

    # Otherwise fall back to ORF ID
    return feature["feature"]
###############################################################################
# Build feature summary
###############################################################################

def build_summary(intervals, annotations):

    summary = []

    for feature in annotations:

        feature_length = (

            feature["end"]

            - feature["start"]

            + 1

        )

        #######################################################################
        # Find samples containing this plasmid
        #######################################################################

        samples = sorted({

            interval["sample"]

            for interval in intervals

            if interval["plasmid"] == feature["plasmid"]

        })

        for sample in samples:

            matching = [

                interval

                for interval in intervals

                if interval["plasmid"] == feature["plasmid"]

                and interval["sample"] == sample

            ]

            segments = []

            contributing = 0

            coverage_percent = matching[0]["coverage"]

            ###############################################################
            # Collect overlap segments
            ###############################################################

            for interval in matching:

                ov = overlap(

                    interval["start"],
                    interval["end"],

                    feature["start"],
                    feature["end"]

                )

                if ov is None:

                    continue

                contributing += 1

                segments.append(

                    (

                        ov[0],

                        ov[1]

                    )

                )

            ###############################################################
            # Merge overlap segments
            ###############################################################

            merged = merge_segments(

                segments

            )

            covered_bp = 0

            segment_strings = []

            for start, end in merged:

                covered_bp += (

                    end

                    - start

                    + 1

                )

                segment_strings.append(

                    f"{start}-{end}"

                )

            covered_percent = round(

                covered_bp

                * 100.0

                / feature_length,

                2

            )

            ###############################################################
            # Coverage pattern
            ###############################################################

            if covered_bp == 0:

                pattern = "Absent"

            elif covered_bp == feature_length:

                pattern = "Complete"

            else:

                pattern = "Fragmented"

            ###############################################################
            # Save row
            ###############################################################

            summary.append({

                "Plasmid": feature["plasmid"],

                "Sample": sample,

                "CoveragePercent": coverage_percent,

                "Feature": feature["feature"],
                
                "Display_Name": get_display_name(feature),

                "FeatureType": feature["type"],

                "Source": feature["source"],

                "Start": feature["start"],

                "End": feature["end"],

                "Length": feature_length,

                "CoveredBP": covered_bp,

                "FeatureCoveredPercent": covered_percent,

                "CoveragePattern": pattern,

                "IntervalsContributing": len(merged),

                "CoveredSegments": ";".join(segment_strings),

                "COG_category": feature["COG_category"],

                "COG_description": feature["COG_description"],

                "Domains": feature["Domains"],

                "AMR_class": feature["AMR_class"],

                "AMR_mechanism": feature["AMR_mechanism"],

                "AMR_gene": feature["AMR_gene"],

                "AMR_identity": feature["AMR_identity"]

            })

    return summary


###############################################################################
# Write output
###############################################################################

def write_summary(summary):

    with open(

        OUTPUT,

        "w",

        newline=""

    ) as out:

        writer = csv.DictWriter(

            out,

            fieldnames=summary[0].keys(),

            delimiter="\t"

        )

        writer.writeheader()

        writer.writerows(summary)


###############################################################################
# Main
###############################################################################

def main():

    print("=" * 70)

    print("Building feature summary")

    print("=" * 70)

    print()

    print("Reading coverage intervals...")

    intervals = read_intervals()

    print(f"  Loaded {len(intervals)} intervals.")

    print()

    print("Reading annotations...")

    annotations = read_annotations()

    print(f"  Loaded {len(annotations)} annotations.")

    print()

    print("Summarizing feature coverage...")

    summary = build_summary(

        intervals,

        annotations

    )

    print(f"  Generated {len(summary)} feature summaries.")

    print()

    print("Writing output...")

    write_summary(

        summary

    )

    print()

    print("=" * 70)

    print("Finished successfully")

    print("=" * 70)

    print()

    print(f"Output written to:\n{OUTPUT}")

    print()


###############################################################################
# Entry point
###############################################################################

if __name__ == "__main__":

    main()
