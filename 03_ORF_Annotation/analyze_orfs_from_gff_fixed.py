
import argparse
import pandas as pd

def parse_gff(input_file):
    plasmid_lengths = {}
    cds_entries = []

    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith("# Sequence Data:"):
                parts = line.strip().split(";")
                seqlen = int(parts[1].split("=")[1])
                seqhdr = parts[2].split("=")[1].strip('"')
                plasmid_lengths[seqhdr] = seqlen // 2  # original length
            elif not line.startswith("#") and "\tCDS\t" in line:
                parts = line.strip().split("\t")
                seqhdr = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                cds_entries.append((seqhdr, start, end, strand))
    return plasmid_lengths, cds_entries

def filter_and_analyze(plasmid_lengths, cds_entries):
    filtered = []
    orf_counts = {}
    orf_lengths = []
    total_orf_lengths = {}

    for seqhdr, start, end, strand in cds_entries:
        if seqhdr not in plasmid_lengths:
            continue

        seq_len = plasmid_lengths[seqhdr] * 2
        lower = plasmid_lengths[seqhdr] // 2
        upper = lower + plasmid_lengths[seqhdr]

        if lower <= start <= upper:
            orf_len = end - start + 1
            filtered.append((seqhdr, start, end, strand, orf_len))
            orf_counts[seqhdr] = orf_counts.get(seqhdr, 0) + 1
            orf_lengths.append((seqhdr, orf_len))
            total_orf_lengths[seqhdr] = total_orf_lengths.get(seqhdr, 0) + orf_len

    orf_counts_df = pd.DataFrame(list(orf_counts.items()), columns=["Plasmid", "ORF_Count"])
    orf_lengths_df = pd.DataFrame(orf_lengths, columns=["Plasmid", "ORF_Length"])
    orf_ratio_df = pd.DataFrame([
        (pid, round(total_orf_lengths[pid] / plasmid_lengths[pid], 2))
        for pid in total_orf_lengths
    ], columns=["Plasmid", "ORF_to_Plasmid_Ratio"])

    return orf_counts_df, orf_lengths_df, orf_ratio_df

def write_to_excel(orf_counts_df, orf_lengths_df, orf_ratio_df, output_file):
    with pd.ExcelWriter(output_file) as writer:
        orf_counts_df.to_excel(writer, sheet_name="ORF_Counts", index=False)
        orf_lengths_df.to_excel(writer, sheet_name="ORF_Lengths", index=False)
        orf_ratio_df.to_excel(writer, sheet_name="ORF_to_Plasmid_Ratio", index=False)

def main():
    parser = argparse.ArgumentParser(description="Analyze Prodigal GFF output for ORFs in plasmid sequences")
    parser.add_argument("--input", required=True, help="Input GFF file from Prodigal")
    args = parser.parse_args()

    plasmid_lengths, cds_entries = parse_gff(args.input)
    orf_counts_df, orf_lengths_df, orf_ratio_df = filter_and_analyze(plasmid_lengths, cds_entries)
    write_to_excel(orf_counts_df, orf_lengths_df, orf_ratio_df, "predictions.xlsx")

if __name__ == "__main__":
    main()
