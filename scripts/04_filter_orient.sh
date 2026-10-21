#!/usr/bin/env bash
# Filter reads by quality and length, then turn them all to the same strand.

# Usage: bash scripts/04_filter_orient.sh [input_dir]   (default: raw_data)
#
# Output per sample, in results/04_filter_orient/:
#   <sample>.filtered.fastq.gz    passed quality/length filter (input for Kraken2)
#   <sample>.oriented.fastq.gz    filtered and oriented (input for clustering)
#   <sample>.unoriented.fastq.gz  filtered, but no similar UNITE sequence found

set -euo pipefail
source config.sh

IN_DIR="${1:-raw_data}"
OUT_DIR="results/04_filter_orient"
LOG="logs/04_filter_orient.log"

shopt -s nullglob  # set nullglob so that the array is empty if no files match

FASTQS=("$IN_DIR"/*.fastq.gz)
if (( ${#FASTQS[@]} == 0 )); then
    echo "Error: no .fastq.gz files in $IN_DIR" >&2
    exit 1
fi
if [[ ! -s "$UNITE_FASTA" ]]; then
    echo "Error: reference $UNITE_FASTA not found or empty." >&2
    exit 1
fi

mkdir -p "$OUT_DIR" logs

count_reads() {
    seqkit stats --tabular "$1" | awk 'NR == 2 { print $4 }' # return number of reads in a FASTQ file
}

{
    echo "Settings: quality >= $MIN_QUAL, length $MIN_LEN-$MAX_LEN bp, $THREADS threads"
    echo "Orientation reference: $UNITE_FASTA"

    for fq in "${FASTQS[@]}"; do
        sample=$(basename "$fq" .fastq.gz)
        filtered="$OUT_DIR/$sample.filtered.fastq.gz"
        oriented="$OUT_DIR/$sample.oriented.fastq.gz"
        unoriented="$OUT_DIR/$sample.unoriented.fastq.gz"

        echo ">>> $sample"
        if [[ -s "$oriented" ]]; then
            echo "already done, skipping"
            continue
        fi

        echo "filtering..."
        zcat "$fq" |
            chopper --quality "$MIN_QUAL" \
                    --minlength "$MIN_LEN" \
                    --maxlength "$MAX_LEN" \
                    --threads "$THREADS" |
            gzip > "$filtered.tmp"
        mv "$filtered.tmp" "$filtered"

        echo "orienting..."
        # Compare each read with UNITE in both directions; reverse-complement
        # reads that match the opposite strand.
        # (an old Illumina limit), but Nanopore reads can go higher.
        vsearch --orient "$filtered" \
                --db "$UNITE_FASTA" \
                --threads "$THREADS" \
                --fastqout "$OUT_DIR/$sample.oriented.fastq" \
                --notmatched "$OUT_DIR/$sample.unoriented.fastq" 
        gzip -f "$OUT_DIR/$sample.unoriented.fastq"
        gzip -f "$OUT_DIR/$sample.oriented.fastq"
    done

    echo "Summary:"
    printf "sample\traw\tfiltered\toriented\tunoriented\tpct_kept\n" > "$OUT_DIR/summary.tsv"
    for fq in "${FASTQS[@]}"; do
        sample=$(basename "$fq" .fastq.gz)
        raw=$(count_reads "$fq")
        filt=$(count_reads "$OUT_DIR/$sample.filtered.fastq.gz")
        ori=$(count_reads "$OUT_DIR/$sample.oriented.fastq.gz")
        unori=$(count_reads "$OUT_DIR/$sample.unoriented.fastq.gz")
        pct=$(awk -v a="$ori" -v b="$raw" 'BEGIN { printf "%.1f", 100 * a / b }')
        printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$sample" "$raw" "$filt" "$ori" "$unori" "$pct" \
            >> "$OUT_DIR/summary.tsv"
    done
    cat "$OUT_DIR/summary.tsv"
} 2>&1 | tee "$LOG"