#!/usr/bin/env bash
set -euo pipefail

#######################
# User-editable config
#######################
IQTREE="/home/naren/Desktop/ANAND/PAN_GENOME/SELECTION_PRESSURE/5_iqtree/iqtree-3.0.1-Linux/bin/iqtree3"
ALN_DIR="/home/naren/Desktop/ANAND/PAN_GENOME/SELECTION_PRESSURE/5_iqtree/6723_file_input_for_iqtree"              # directory with codon alignments
PATTERN="*_codon_aln.fasta"
THREADS="AUTO"
MERIT="AICc"
OUT_ROOT="iqtree_codon_runs_6723_FAILED_SOLVE"
KEEP_EXISTING=false

#######################
# Prepare absolute paths
#######################
OUT_ROOT="$(realpath "$OUT_ROOT")"
FAILED_LOG="$OUT_ROOT/failed_alignments.txt"
FAILED_DIR="$OUT_ROOT/failed_alignments_folder"

#######################
# Safety checks
#######################
if [ ! -x "$IQTREE" ]; then
  echo "ERROR: iqtree binary not found or not executable at: $IQTREE"
  exit 1
fi

if [ ! -d "$ALN_DIR" ]; then
  echo "ERROR: alignment directory not found: $ALN_DIR"
  exit 1
fi

# Ensure output dirs and log exist
mkdir -p "$OUT_ROOT" "$FAILED_DIR"
touch "$FAILED_LOG"

echo "Outputs will be placed in: $OUT_ROOT"
echo "Failed files will be logged to: $FAILED_LOG and stored in: $FAILED_DIR"
echo

#######################
# Loop over alignments
#######################
shopt -s nullglob
for aln in "$ALN_DIR"/$PATTERN; do
  base=$(basename "$aln")
  prefix="${base%.*}"
  run_outdir="$OUT_ROOT/$prefix"
  mkdir -p "$run_outdir"

  if [ -f "$run_outdir/$prefix.iqtree" ] && [ "$KEEP_EXISTING" = true ]; then
    echo "Skipping $base (already has results in $run_outdir)."
    continue
  fi

  echo "-------------------------------------------"
  echo "Running IQ-TREE on: $base"
  echo " -> output prefix: $prefix"
  echo " -> output dir: $run_outdir"
  echo

  cp "$aln" "$run_outdir/$base"
  cd "$run_outdir" || continue

  #########################
  # Detect variable sites
  #########################
  total_sites=$(grep -v ">" "$base" | tr -d "\n" | awk '{print length}')
  var_sites=$(awk '
    BEGIN{FS=""; OFS=""; nseq=0}
    /^>/{next}
    {seq[nseq++]=$0}
    END{
      L=length(seq[0]);
      vars=0;
      for(i=1;i<=L;i++){
        split("", a);
        for(j=0;j<nseq;j++){a[substr(seq[j],i,1)]++}
        if(length(a)>1) vars++;
      }
      print vars;
    }' "$base")

  perc_const=$(awk -v v=$var_sites -v t=$total_sites 'BEGIN{printf "%.2f", 100*(1 - v/t)}')

  echo "Alignment stats for $base: $var_sites variable of $total_sites sites (${perc_const}% constant)"

  #########################
  # Choose run mode
  #########################
  if [ "$var_sites" -lt 10 ] || (( $(echo "$perc_const > 90" | bc -l) )); then
    echo "⚠️  Low variation detected: switching to SAFE mode (simple model)"
    CMD=( "$IQTREE" -s "$base" -st CODON -m GY -safe -T "$THREADS" -pre "$prefix" )
  else
    CMD=( "$IQTREE" -s "$base" -st CODON -m MFP -merit "$MERIT" -bnni -T "$THREADS" -pre "$prefix" )
  fi

  echo "Running command: ${CMD[*]}"
  if ! "${CMD[@]}"; then
    echo "❌ ERROR: IQ-TREE failed for $base"
    echo "$base" >> "$FAILED_LOG"
    mv "$aln" "$FAILED_DIR/" 2>/dev/null || true
    cd - > /dev/null || true
    continue
  fi

  echo "✅ Finished successfully: $base"
  cd - > /dev/null || true
done

echo
echo "==========================================="
echo "All IQ-TREE analyses completed."
echo "Results directory: $OUT_ROOT"
echo "Failed alignments logged in: $FAILED_LOG"
echo "Failed files stored in: $FAILED_DIR"
echo "==========================================="
