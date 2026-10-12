#!/usr/bin/env bash
# Build a Kraken2 database from a UNITE FASTA file.
# Usage: bash scripts/02_build_kraken_db.sh <unite.fasta>
#Example:
# docker run --rm --user "$(id -u):$(id -g)" -v "$PWD":/data fungiflow-remastered bash scripts/02_build_kraken_db.sh database/unite/sh_general_release_dynamic_all_19.02.2025.fasta

set -euo pipefail

UNITE_FASTA="${1:?Usage: $0 <unite.fasta>}"
DB_DIR="database/kraken2"
THREADS="${THREADS:-$(nproc)}"      # Use all available CPU threads by default
LOG="logs/02_build_kraken_db.log" 

if [[ ! -s "$UNITE_FASTA" ]]; then
    echo "Error: FASTA file '$UNITE_FASTA' not found or empty." >&2
    exit 1
fi
if [[ -e "$DB_DIR" ]]; then
    echo "Error: $DB_DIR already exists. Delete it first to rebuild." >&2
    exit 1
fi

mkdir -p "$DB_DIR" logs

{
    echo "Converting UNITE to Kraken2 format..."
    python3 scripts/py/unite_to_kraken.py "$UNITE_FASTA" "$DB_DIR"
    echo "Adding sequences to the Kraken2 library..."
    kraken2-build --add-to-library "$DB_DIR/library.fasta" --db "$DB_DIR"
    echo "Building the database with $THREADS threads..."
    kraken2-build --build --db "$DB_DIR" --threads "$THREADS"
    echo "Results:"
    ls -lh "$DB_DIR"/*.k2d
    kraken2-inspect --db "$DB_DIR" > "$DB_DIR/inspect.txt"
    head -n 15 "$DB_DIR/inspect.txt"

} 2>&1 | tee "$LOG"