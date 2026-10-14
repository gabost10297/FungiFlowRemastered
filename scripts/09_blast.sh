#!/usr/bin/env bash
# Identify every OTU by BLAST against UNITE and give it a name and confidence level.
# Usage: bash scripts/09_blast.sh
#
# Output in results/09_blast/:
#   blast_hits.tsv   up to 10 hits per OTU
set -euo pipefail
source config.sh

OTUS="results/08_otus/otus.fasta"
TABLE="results/08_otus/otu_table.tsv"
OUT_DIR="results/09_blast"
LOG="logs/09_blast.log"

if [[ ! -s "$OTUS" || ! -s "$TABLE" ]]; then
    echo "Error: OTU files not found. Run script 08 first." >&2
    exit 1
fi
mkdir -p "$OUT_DIR" logs

{
    echo "Thresholds: high >= ${BLAST_HIGH_PIDENT}% identity and ${BLAST_HIGH_QCOV}% coverage,"
    echo "            low >= ${BLAST_LOW_PIDENT}% and ${BLAST_LOW_QCOV}%,"
    echo "            genus >= ${BLAST_GENUS_PIDENT}% and ${BLAST_GENUS_QCOV}%"

    blastn -query "$OTUS" \
           -db "$BLAST_DB" \
           -max_target_seqs 10 \
           -max_hsps 1 \
           -evalue 1e-10 \
           -num_threads "$THREADS" \
           -outfmt "6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue bitscore qcovhsp" \
           -out "$OUT_DIR/blast_hits.tsv"
    echo "BLAST hits: $(wc -l < "$OUT_DIR/blast_hits.tsv")"

    python3 scripts/py/assign_taxonomy.py \
        "$OUT_DIR/blast_hits.tsv" "$TABLE" "$OUT_DIR/taxonomy.tsv" \
        --high-pident "$BLAST_HIGH_PIDENT" --high-qcov "$BLAST_HIGH_QCOV" \
        --low-pident "$BLAST_LOW_PIDENT" --low-qcov "$BLAST_LOW_QCOV" \
        --genus-pident "$BLAST_GENUS_PIDENT" --genus-qcov "$BLAST_GENUS_QCOV" \
        --ambiguous-gap "$BLAST_AMBIGUOUS_GAP" --min-ref-cov "$BLAST_MIN_REF_COV"
} 2>&1 | tee "$LOG"