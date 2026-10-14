#!/usr/bin/env bash
# Group identified OTUs into organisms, count them per sample and pick the top ones.
# Usage: bash scripts/10_organisms.sh
#
# Output in results/10_organisms/:
#   organisms.tsv           one row per organism: reads per sample, representative OTU
#   otu_groups.tsv          the group (and organism) of every OTU
#   per_sample.tsv          number of organisms per sample
#   representatives.fasta   best OTU of every organism (for the tree)
#   top/NN_<name>.fasta     all consensus sequences of the top organisms (for alignment)

set -euo pipefail
source config.sh

TAXONOMY="results/09_blast/taxonomy.tsv"
OTU_DIR="results/08_otus"
OUT_DIR="results/10_organisms"
LOG="logs/10_organisms.log"

if [[ ! -s "$TAXONOMY" ]]; then
    echo "Error: $TAXONOMY not found. Run script 09 first." >&2
    exit 1
fi
rm -rf "$OUT_DIR/top"
mkdir -p "$OUT_DIR" logs

{
    python3 scripts/py/organisms.py \
        "$TAXONOMY" \
        "$OTU_DIR/otu_table.tsv" \
        "$OTU_DIR/otus.fasta" \
        "$OTU_DIR/all_its.fasta" \
        "$OTU_DIR/members.tsv" \
        "$OUT_DIR" \
        --top "$TOP_ORGANISMS"
} 2>&1 | tee "$LOG"