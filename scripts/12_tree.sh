#!/usr/bin/env bash
# Build a phylogenetic tree from the representatives of the most abundant organisms.
# Usage: bash scripts/12_tree.sh
#
# Output in results/12_tree/:
#   tree_input.fasta       representatives of the top organisms
#   tree_input.aln.fasta   MAFFT alignment
#   tree.treefile          maximum likelihood tree with bootstrap support (Newick)
#   tree.contree           bootstrap consensus tree
#   tree.iqtree            IQ-TREE report: chosen model, likelihood, text drawing of the tree
#   tree.log               full IQ-TREE log
set -euo pipefail
source config.sh

REPS="results/10_organisms/representatives.fasta"
ORGANISMS="results/10_organisms/organisms.tsv"
OUT_DIR="results/12_tree"
LOG="logs/12_tree.log"

if [[ ! -s "$REPS" || ! -s "$ORGANISMS" ]]; then
    echo "Error: results of script 10 not found. Run script 10 first." >&2
    exit 1
fi
rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR" logs

{
    # representatives.fasta is sorted by reads, so the first N are the top organisms.
    seqkit head -n "$TREE_ORGANISMS" "$REPS" | seqkit seq -i > "$OUT_DIR/tree_input.fasta"
    echo "Organisms in the tree: $(grep -c '>' "$OUT_DIR/tree_input.fasta")"

    mafft --localpair --maxiterate 1000 --thread "$THREADS" --quiet \
        "$OUT_DIR/tree_input.fasta" > "$OUT_DIR/tree_input.aln.fasta"

    # Root the tree between the phyla: all Basidiomycota in the tree are the outgroup.
    outgroup=$(awk -F'\t' -v n="$TREE_ORGANISMS" '
        NR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
        NR <= n + 1 && $col["phylum"] == "Basidiomycota" {
            id = $col["organism"]; gsub(" ", "_", id); gsub("\\.", "", id)
            printf "%s%s", sep, id; sep = ","
        }' "$ORGANISMS")
    root=()
    if [[ -n "$outgroup" ]]; then
        root=(-o "$outgroup")
        echo "Rooted on Basidiomycota: $outgroup"
    fi

    iqtree3 -s "$OUT_DIR/tree_input.aln.fasta" \
            -m MFP \
            -B "$TREE_BOOTSTRAP" \
            -T AUTO \
            --seed 1 \
            "${root[@]}" \
            --prefix "$OUT_DIR/tree" \
            --quiet

    grep -E "Best-fit model|Log-likelihood of the tree" "$OUT_DIR/tree.iqtree"
    grep "WARNING" "$OUT_DIR/tree.log" || true
    echo "Tree: $OUT_DIR/tree.treefile"
} 2>&1 | tee "$LOG"