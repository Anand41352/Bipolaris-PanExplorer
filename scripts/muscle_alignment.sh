#!/bin/bash
# -----------------------------
# Batch Protein Alignment using MUSCLE v5
# High accuracy mode
# -----------------------------

# Input and output folders
PROTEIN_DIR="/home/naren/Desktop/ANAND/PAN_GENOME/SELECTION_PRESSURE/3_MACASE/10093_PROTEIN"
OUT_DIR="/home/naren/Desktop/ANAND/PAN_GENOME/SELECTION_PRESSURE/3_MACASE/aligned_PROTEIN_10093"

# Create output folder if it doesn't exist
mkdir -p "$OUT_DIR"

# Loop through all protein FASTA files
for file in "$PROTEIN_DIR"/*.fa; do
    base=$(basename "$file" .fa)
    echo "Aligning protein $base ..."

    # MUSCLE v5 alignment (default is high accuracy)
    muscle -align "$file" -output "$OUT_DIR/${base}_aa_aln.fa"
done

echo "All protein alignments completed."
