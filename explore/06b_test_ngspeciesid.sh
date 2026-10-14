#!/usr/bin/env bash
# Experiment: cluster one sample with NGSpeciesID and polish the consensuses with medaka.
# Runs in the fungiflow-ngspeciesid container.
# Usage: bash explore/06b_test_ngspeciesid.sh [sample] [n_reads]
#        defaults: barcode25, 10000 random reads

set -euo pipefail
source config.sh

SAMPLE="${1:-barcode25}"
N_READS="${2:-10000}"
IN="results/04_filter_orient/$SAMPLE.oriented.fastq.gz"
OUT_DIR="results/explore/06b_ngspeciesid/$SAMPLE"
LOG="$OUT_DIR/run.log"

if [[ ! -s "$IN" ]]; then
    echo "Error: $IN not found. Run script 04 first." >&2
    exit 1
fi
if [[ -e "$OUT_DIR" ]]; then
    echo "Error: $OUT_DIR already exists. Delete it first to rerun." >&2
    exit 1
fi
mkdir -p "$OUT_DIR"

# NGSpeciesID keeps clusters larger than (abundance_ratio x number of reads).
# Convert our minimum cluster size in reads into that fraction.
ratio=$(awk -v m="$MIN_CLUSTER_SIZE" -v n="$N_READS" 'BEGIN { printf "%.6f", m / n }')

# NGSpeciesID reads uncompressed FASTQ only.
zcat "$IN" > "$OUT_DIR/reads.fastq"

{
    echo "Settings: $N_READS reads, clusters >= $MIN_CLUSTER_SIZE reads (ratio $ratio)"
    echo "Medaka model: $MEDAKA_MODEL, $THREADS threads"

    NGSpeciesID --ont --consensus --medaka \
                --medaka_model "$MEDAKA_MODEL" \
                --fastq "$OUT_DIR/reads.fastq" \
                --outfolder "$OUT_DIR/ngspeciesid" \
                --sample_size "$N_READS" \
                --seed 1 \
                --q "$MIN_QUAL" \
                --abundance_ratio "$ratio" \
                --max_seqs_for_consensus "$MAX_SEQS_CONSENSUS" \
                --sample_name "$SAMPLE" \
                --t "$THREADS"

    echo "Finished in $SECONDS seconds"
} 2>&1 | tee "$LOG"

rm "$OUT_DIR/reads.fastq"