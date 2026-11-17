#!/usr/bin/env python3

import os
from Bio import SeqIO
import csv

# ---------------------------
# Define input directories
# ---------------------------
with_gap_dir = "/home/naren/Desktop/ANAND/PAN_GENOME/SELECTION_PRESSURE/4_pal2nal/codon_aligned_output_without_no_gap_option"
no_gap_dir = "/home/naren/Desktop/ANAND/PAN_GENOME/SELECTION_PRESSURE/4_pal2nal/codon_aligned_output_with_gap_removal_option"
output_csv = "alignment_length_comparison.csv"

# ---------------------------
# Helper function to get alignment length
# ---------------------------
def get_alignment_length(fasta_file):
    """
    Returns alignment length (all sequences assumed equal length)
    """
    for record in SeqIO.parse(fasta_file, "fasta"):
        return len(record.seq)
    return 0

# ---------------------------
# Collect lengths
# ---------------------------
data = []

for file_name in sorted(os.listdir(with_gap_dir)):
    if file_name.endswith(".fasta"):
        orthogroup = os.path.splitext(file_name)[0]

        with_gap_path = os.path.join(with_gap_dir, file_name)
        no_gap_path = os.path.join(no_gap_dir, file_name)

        # Calculate lengths
        len_with_gap = get_alignment_length(with_gap_path)
        len_no_gap = get_alignment_length(no_gap_path) if os.path.exists(no_gap_path) else "NA"

        # Calculate % gap content
        if isinstance(len_no_gap, int) and len_with_gap > 0:
            gap_percentage = round(100 * (len_with_gap - len_no_gap) / len_with_gap, 2)
        else:
            gap_percentage = "NA"

        data.append([orthogroup, len_with_gap, len_no_gap, gap_percentage])

# ---------------------------
# Write CSV report
# ---------------------------
with open(output_csv, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Orthogroup", "Alignment_Length_With_Gap", "Alignment_Length_No_Gap", "Gap_Percentage"])
    writer.writerows(data)

print(f"✅ Report generated: {output_csv}")

