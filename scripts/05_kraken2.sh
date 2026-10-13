#!/usr/bin/env bash
# Classify every filtered read with Kraken2 and summarise the results per sample.
# Usage: bash scripts/05_kraken2.sh
#
# Output in results/05_kraken2/:
#   <sample>.report.txt     Kraken2 report: reads per taxon
#   <sample>.reads.tsv.gz   one line per read: read id and assigned taxid
#   classification.tsv      % classified / Fungi / plants / animals per sample
#   genus_reads.tsv         reads per genus and sample

set -euo pipefail
source config.sh

IN_DIR="results/04_filter_orient"
OUT_DIR="results/05_kraken2"
LOG="logs/05_kraken2.log"

shopt -s nullglob
FASTQS=("$IN_DIR"/*.filtered.fastq.gz)
if (( ${#FASTQS[@]} == 0 )); then
    echo "Error: no filtered reads in $IN_DIR. Run script 04 first." >&2
    exit 1
fi
if [[ ! -s "$KRAKEN_DB/hash.k2d" ]]; then
    echo "Error: no Kraken2 database in $KRAKEN_DB. Run script 02 first." >&2
    exit 1
fi

mkdir -p "$OUT_DIR" logs

{
    echo "Settings: confidence $KRAKEN_CONFIDENCE, minimum hit groups $KRAKEN_MIN_HITS, $THREADS threads"

    for fq in "${FASTQS[@]}"; do
        sample=$(basename "$fq" .filtered.fastq.gz)
        report="$OUT_DIR/$sample.report.txt"

        echo ">>> $sample"
        if [[ -s "$report" ]]; then
            echo "already done, skipping"
            continue
        fi

        kraken2 --db "$KRAKEN_DB" \
                --threads "$THREADS" \
                --confidence "$KRAKEN_CONFIDENCE" \
                --minimum-hit-groups "$KRAKEN_MIN_HITS" \
                --gzip-compressed \
                --report "$report.tmp" \
                "$fq" |
            gzip > "$OUT_DIR/$sample.reads.tsv.gz"
        mv "$report.tmp" "$report"
    done

    echo "Summary:"
    python3 scripts/py/kraken_summary.py "$OUT_DIR"
} 2>&1 | tee "$LOG"