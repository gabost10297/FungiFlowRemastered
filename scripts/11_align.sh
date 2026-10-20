#!/usr/bin/env bash
# Align all consensus sequences of each top organism and compare them with its representative.
# Usage: bash scripts/11_align.sh
#
# Output in results/11_alignments/:
#   NN_<organism>.aln.fasta   MAFFT alignment (open in AliView or Jalview)
#   vs_representative.tsv     every sequence compared with its organism's representative
#   summary.tsv               one row per organism
set -euo pipefail
source config.sh

IN_DIR="results/10_organisms/top"
OUT_DIR="results/11_alignments"
LOG="logs/11_align.log"

shopt -s nullglob
files=("$IN_DIR"/*.fasta)
if (( ${#files[@]} == 0 )); then
    echo "Error: no FASTA files in $IN_DIR. Run script 10 first." >&2
    exit 1
fi
rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR" logs

{
    for fa in "${files[@]}"; do
        name=$(basename "$fa" .fasta)
        mafft --auto --thread "$THREADS" --quiet "$fa" > "$OUT_DIR/$name.aln.fasta"
        echo "Aligned $name: $(grep -c '>' "$fa") sequences"
    done
    echo
    python3 scripts/py/compare_alignment.py \
        results/10_organisms/organisms.tsv "$OUT_DIR" \
        --same "$SAME_ITS_IDENTITY" --length-diff "$LENGTH_DIFF"
} 2>&1 | tee "$LOG"