#!/bin/bash

# ===== SETTINGS =====
BASE_DIR="/home/naren/Desktop/ANAND/PAN_GENOME/SELECTION_PRESSURE/7_busted/BATCH_1/3775/6_FAILED"
THREADS_PER_JOB=10
MAX_PARALLEL_JOBS=8
# =====================

export OMP_NUM_THREADS=$THREADS_PER_JOB

run_analysis () {

    OG_DIR=$1
    OG_NAME=$(basename "$OG_DIR")

    ALIGN="$OG_DIR/$OG_NAME.fasta"
    TREE="$OG_DIR/$OG_NAME.tree"

    echo "======================================"
    echo "Running analyses for $OG_NAME"
    echo "======================================"

    # ---------------- BUSTED ----------------
    hyphy busted \
        --alignment "$ALIGN" \
        --tree "$TREE" \
        --output "$OG_DIR/${OG_NAME}_BUSTED.json" \
        CPU=$THREADS_PER_JOB \
        > "$OG_DIR/${OG_NAME}_BUSTED.log" 2>&1

    # ---------------- MEME ----------------
    hyphy meme \
        --alignment "$ALIGN" \
        --tree "$TREE" \
        --output "$OG_DIR/${OG_NAME}_MEME.json" \
        CPU=$THREADS_PER_JOB \
        > "$OG_DIR/${OG_NAME}_MEME.log" 2>&1

    # ---------------- FEL ----------------
    hyphy fel \
        --alignment "$ALIGN" \
        --tree "$TREE" \
        --output "$OG_DIR/${OG_NAME}_FEL.json" \
        CPU=$THREADS_PER_JOB \
        > "$OG_DIR/${OG_NAME}_FEL.log" 2>&1


    # ---------------- aBSREL ----------------
    hyphy absrel \
        --alignment "$ALIGN" \
        --tree "$TREE" \
        --output "$OG_DIR/${OG_NAME}_aBSREL.json" \
        CPU=$THREADS_PER_JOB \
        > "$OG_DIR/${OG_NAME}_aBSREL.log" 2>&1

    echo "Finished $OG_NAME"
}

export -f run_analysis

find "$BASE_DIR" -maxdepth 1 -type d -name "OG*" | \
    xargs -I {} -P $MAX_PARALLEL_JOBS bash -c 'run_analysis "$@"' _ {}
