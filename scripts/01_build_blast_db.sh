#!/usr/bin/env bash
# Builds a BLAST database from a FASTA file.
# How to use: bash scripts/01_build_blast_db.sh <unite.fasta>

# Example: 
#   docker run --rm --user "$(id -u):$(id -g)" -v "$PWD":/data fungiflow-remastered \
#   bash scripts/01_build_blast_db.sh database/unite/sh_general_release_dynamic_all_19.02.2025.fasta

set -euo pipefail #stop at first error

UNITE_FASTA="${1:?Usage: $0 <unite.fasta>}" #take the fasta file from command line
OUT_DIR="database/blast" 
DB_NAME="unite"
LOG="logs/01_build_blast_db.log"

mkdir -p "$OUT_DIR" logs

if [[ ! -f "$UNITE_FASTA" ]]; then
    echo "Error: FASTA file '$UNITE_FASTA' not found." >&2
    exit 1
fi

{

    echo "INPUT FASTA: $UNITE_FASTA"
    echo "seqkit stats (log how many sequences and how long they are): "
    seqkit stats "$UNITE_FASTA"

    echo "Building BLAST database..."
    makeblastdb \
        -in "$UNITE_FASTA" \
        -dbtype nucl \
        -out "$OUT_DIR/$DB_NAME" \
        -title "$(basename "$UNITE_FASTA")"

    echo "BLAST database built successfully. Results: "
    blastdbcmd -db "$OUT_DIR/$DB_NAME" -info  #read the database and print info about it

} 2>&1 | tee "$LOG"