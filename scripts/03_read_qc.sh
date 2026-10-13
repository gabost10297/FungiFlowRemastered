#!/usr/bin/env bash
# Summarise read counts, lengths and quality for every raw FASTQ file (barcodes)

# Usage: bash scripts/03_read_qc.sh [input_dir]   (default: raw_data)

set -euo pipefail

IN_DIR="${1:-raw_data}"
OUT_DIR="results/03_read_qc"
THREADS="${THREADS:-$(nproc)}"
LOG="logs/03_read_qc.log"


shopt -s nullglob   # set nullglob so that the array is empty if no files match

FASTQS=("$IN_DIR"/*.fastq.gz)
if (( ${#FASTQS[@]} == 0 )); then
    echo "Error: no .fastq.gz files in $IN_DIR" >&2
    exit 1
fi

mkdir -p "$OUT_DIR" logs

{
    echo "Read statistics for ${#FASTQS[@]} samples"

    seqkit stats --all --tabular --threads "$THREADS" "${FASTQS[@]}" \
        > "$OUT_DIR/read_stats.tsv"
    
    # Print a few columns, looked up by name because their position may change in future versions of seqkit

    awk -F'\t' '
        NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i
                    printf "%-12s %9s %6s %6s %6s %6s %6s %8s\n",
                            "sample", "reads", "min", "Q1", "median", "Q3", "max", "avg_qual"
                    next }
        { sample = $col["file"]; sub(".*/", "", sample); sub(".fastq.gz", "", sample)
          printf "%-12s %9s %6s %6s %6s %6s %6s %8s\n",
                 sample, $col["num_seqs"], $col["min_len"], $col["Q1"], $col["Q2"],
                 $col["Q3"], $col["max_len"], $col["AvgQual"] }
    ' "$OUT_DIR/read_stats.tsv"


    echo "Percentage of reads per length and quality class"
    printf "sample\treads\tlen_lt300\tlen_300_499\tlen_500_899\tlen_ge900\tq_lt10\tq_10_14\tq_15_19\tq_ge20\n" \
        > "$OUT_DIR/length_quality_classes.tsv"
    

    for fq in "${FASTQS[@]}"; do
        sample=$(basename "$fq" .fastq.gz)
        # One line per read: id, length, average quality
        seqkit fx2tab --name --only-id --length --avg-qual "$fq" |
        awk -F'\t' -v s="$sample" '
            { n++
              len = $2; q = $3
              if      (len < 300) l1++; else if (len < 500) l2++
              else if (len < 900) l3++; else                l4++
              if      (q < 10)    q1++; else if (q < 15)    q2++
              else if (q < 20)    q3++; else                q4++ }
            END { printf "%s\t%d\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\n",
                         s, n, 100*l1/n, 100*l2/n, 100*l3/n, 100*l4/n,
                         100*q1/n, 100*q2/n, 100*q3/n, 100*q4/n }'
    done >> "$OUT_DIR/length_quality_classes.tsv"

    cat "$OUT_DIR/length_quality_classes.tsv"
} 2>&1 | tee "$LOG"
