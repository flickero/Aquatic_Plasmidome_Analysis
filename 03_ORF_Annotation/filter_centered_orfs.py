from Bio import SeqIO
import re

input_file = "all_ORFs.fna"       # Replace with your actual .fna file
output_file = "filtered_orfs_centered.fna"

def extract_length_and_start(header):
    # Extract original plasmid length (before doubling)
    length_match = re.search(r"length_(\d+)", header)
    if not length_match:
        return None, None
    orig_len = int(length_match.group(1))

    # Extract start coordinate from header
    parts = header.split("#")
    if len(parts) < 2:
        return None, None
    start = int(parts[1].strip())

    return orig_len, start

kept = 0
total = 0

with open(input_file, "r") as infile, open(output_file, "w") as outfile:
    for record in SeqIO.parse(infile, "fasta"):
        total += 1
        header = record.description
        orig_len, start = extract_length_and_start(header)
        if orig_len is None or start is None:
            continue

        doubled_len = orig_len * 2
        lower = doubled_len // 4
        upper = (3 * doubled_len) // 4

        if lower <= start <= upper:
            SeqIO.write(record, outfile, "fasta")
            kept += 1

print(f"✅ Filtered {kept} of {total} ORFs (with start inside 25–75% of doubled plasmid length)")
