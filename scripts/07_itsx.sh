#!/usr/bin/env bash
# Trim every consensus to the full ITS region (ITS1-5.8S-ITS2) with ITSx.
# Usage: bash scripts/07_itsx.sh

# Output in results/07_itsx/:
#   <sample>.its.fasta   full ITS of every consensus in which ITSx found it
#   itsx/<sample>.*      all ITSx output: single regions, positions, no detections
#   summary.tsv          consensuses and reads kept per sample

set -euo pipefail
source config.sh

IN_DIR="results/06_clusters"
OUT_DIR="results/07_itsx"
LOG="logs/07_itsx.log"

shopt -s nullglob

FASTAS=("$IN_DIR"/*.consensus.fasta)
if (( ${#FASTAS[@]} == 0 )); then
    echo "Error: no consensus files in $IN_DIR. Run script 06 first." >&2
    exit 1
fi

mkdir -p "$OUT_DIR/itsx" logs

# Print the number of sequences and the total reads behind them (from ";size=") in a FASTA file, or "0 0" if the file does not exist

seq_stats() {
    if [[ -e "$1" ]]; then
        awk -F'size=' '/^>/ { n++; r += $2 } END { print n + 0, r + 0 }' "$1"
    else
        echo "0 0"
    fi
}

{
    for fa in "${FASTAS[@]}"; do
        sample=$(basename "$fa" .consensus.fasta)
        prefix="$OUT_DIR/itsx/$sample"
        its="$OUT_DIR/$sample.its.fasta"
        echo ">>> $sample"
        if [[ -s "$its" ]]; then
            echo "already done, skipping"
            continue
        fi
        ITSx -i "$fa" \
             -o "$prefix" \
             --cpu "$THREADS" \
             -E "$ITSX_EVALUE" \
             --preserve T \
             --save_regions all \
             --graphical F
        if [[ -e "$prefix.full.fasta" ]]; then
            cp "$prefix.full.fasta" "$its"
        else
            : > "$its"
        fi
    done

    echo "Summary:"
    printf "sample\tconsensuses\tfull_ITS\tITS1_found\tITS2_found\tno_ITS\treads\treads_full_ITS\tpct_reads_full_ITS\n" \
        > "$OUT_DIR/summary.tsv"
    for fa in "${FASTAS[@]}"; do
        sample=$(basename "$fa" .consensus.fasta)
        prefix="$OUT_DIR/itsx/$sample"
        read -r n_in r_in     < <(seq_stats "$fa")
        read -r n_full r_full < <(seq_stats "$prefix.full.fasta")
        read -r n_its1 _      < <(seq_stats "$prefix.ITS1.fasta")
        read -r n_its2 _      < <(seq_stats "$prefix.ITS2.fasta")
        read -r n_none _      < <(seq_stats "${prefix}_no_detections.fasta")
        pct=$(awk -v a="$r_full" -v b="$r_in" 'BEGIN { printf "%.1f", 100 * a / b }')
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$sample" "$n_in" "$n_full" "$n_its1" "$n_its2" "$n_none" "$r_in" "$r_full" "$pct" \
            >> "$OUT_DIR/summary.tsv"
    done
    cat "$OUT_DIR/summary.tsv"
    echo "Length of the full ITS sequences:"
    seqkit stats --all "$OUT_DIR"/*.its.fasta
} 2>&1 | tee "$LOG"

