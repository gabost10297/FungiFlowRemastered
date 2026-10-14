#!/usr/bin/env bash
# Experiment: how many clusters does vsearch make at different identity thresholds?
# Usage: bash explore/06_test_clustering.sh [sample] [n_reads]
#        defaults: barcode25, 10000 random reads

set -euo pipefail
source config.sh

SAMPLE="${1:-barcode25}"
N_READS="${2:-10000}"
IN="results/04_filter_orient/$SAMPLE.oriented.fastq.gz"
OUT_DIR="results/explore/06_clustering/$SAMPLE"
MIN_SIZE=20

if [[ ! -s "$IN" ]]; then
    echo "Error: $IN not found. Run script 04 first." >&2
    exit 1
fi
mkdir -p "$OUT_DIR"

# A fixed random subset keeps the test fast and repeatable.
seqkit shuffle --rand-seed 1 "$IN" -o "$OUT_DIR/shuffled.fastq.gz" 2>/dev/null
seqkit head -n "$N_READS" "$OUT_DIR/shuffled.fastq.gz" |
    seqkit fq2fa > "$OUT_DIR/subset.fasta"
rm "$OUT_DIR/shuffled.fastq.gz"

total=$(grep -c '^>' "$OUT_DIR/subset.fasta")
echo "Clustering $total reads of $SAMPLE"

SUMMARY="$OUT_DIR/summary.tsv"
printf "identity\tclusters\tclusters_ge%d\tpct_reads_in_ge%d\tlargest\n" \
    "$MIN_SIZE" "$MIN_SIZE" > "$SUMMARY"

for id in 0.80 0.85 0.88 0.90 0.95 0.98; do
    echo "identity $id..."
    uc="$OUT_DIR/id_$id.uc"
    vsearch --cluster_fast "$OUT_DIR/subset.fasta" \
            --id "$id" \
            --strand plus \
            --threads "$THREADS" \
            --uc "$uc" \
            --quiet

    # In a .uc file, every cluster has one line starting with "C";
    # column 3 is the number of reads in that cluster.
    awk -F'\t' -v id="$id" -v min="$MIN_SIZE" -v total="$total" '
        $1 == "C" { n++
                    if ($3 >= min) { big++; in_big += $3 }
                    if ($3 > largest) largest = $3 }
        END { printf "%s\t%d\t%d\t%.1f\t%d\n", id, n, big, 100 * in_big / total, largest }
    ' "$uc" >> "$SUMMARY"
done

cat "$SUMMARY"