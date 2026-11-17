#!/bin/bash
# ==========================================
# DeepFRI batch prediction script
# Author: Anand
# ==========================================

# Input and output folders
INPUT_DIR="input_sequences_for_annotation"
OUTPUT_DIR="results_go_prediction"

# Make sure output directory exists
mkdir -p "$OUTPUT_DIR"

# List of ontologies
ONTOLOGIES=("mf" "bp" "cc")

# Loop through all FASTA files in input folder
for fasta in $INPUT_DIR/*.fasta; do
    # Get base name without path and extension
    base=$(basename "$fasta" .fasta)
    echo "=========================================="
    echo "🔍 Processing $base ..."
    echo "=========================================="

    # Make output subfolder for this dataset
    mkdir -p "$OUTPUT_DIR/$base"

    # Run DeepFRI for each ontology
    for ont in "${ONTOLOGIES[@]}"; do
        echo "▶️ Running ontology: $ont on $base ..."
        python predict.py \
            --fasta_fn "$fasta" \
            --ont "$ont" \
            --output_fn_prefix "$OUTPUT_DIR/$base/DeepFRI_${ont}_${base}" \
            --verbose
    done

    echo "✅ Finished $base"
    echo
done

echo "🎯 All predictions completed successfully!"

