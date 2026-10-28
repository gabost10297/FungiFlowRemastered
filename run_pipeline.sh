#!/usr/bin/env bash
# Run steps 03-12 for one run folder, every step in its own Docker container.
# Usage: bash run_pipeline.sh <run name>
#
# The barcodes must be in runs/<run name>/raw_data/*.fastq.gz. Each run folder
# works like a small copy of the project: its own raw_data/, results/, logs/ and
# config.sh, with scripts/ and database/ linked to the shared ones.
# Steps 01-02 (the reference databases) run first if they are missing.
#
# Writes runs/<run name>/status.tsv: one line per step with step, state, time.
# Started from the app it runs inside fungiflow-viz; from a terminal it works the same.

set -euo pipefail
cd "$(dirname "$0")"

RUN="${1:?Usage: $0 <run name>}"
RUN_DIR="runs/$RUN"
STATUS="$RUN_DIR/status.tsv"

if [[ ! "$RUN" =~ ^[A-Za-z0-9_-]+$ ]]; then
    echo "Error: run name may only contain letters, digits, _ and -" >&2
    exit 1
fi
shopt -s nullglob
FASTQS=("$RUN_DIR"/raw_data/*.fastq.gz)
if (( ${#FASTQS[@]} == 0 )); then
    echo "Error: no .fastq.gz files in $RUN_DIR/raw_data" >&2
    exit 1
fi

# How the step containers see the project at /data. The desktop app sets
# DATA_CONTAINER, a stopped container that holds the mount, so no host path is
# needed. Otherwise "docker run -v" needs the host path, which the Docker engine reads.
if [[ -n "${DATA_CONTAINER:-}" ]]; then
    MOUNT=(--volumes-from "$DATA_CONTAINER")
else
    MOUNT=(-v "${HOST_PROJECT:-$PWD}":/data)
fi
RUN_AS="${RUN_AS:-$(id -u):$(id -g)}"

mkdir -p "$RUN_DIR/results" "$RUN_DIR/logs"
ln -sfn ../../scripts "$RUN_DIR/scripts"
ln -sfn ../../database "$RUN_DIR/database"
[[ -f "$RUN_DIR/config.sh" ]] || cp config.sh "$RUN_DIR/config.sh"

source config.sh

mark() {
    printf "%s\t%s\t%s\n" "$1" "$2" "$(date '+%Y-%m-%d %H:%M:%S')" >> "$STATUS"
}

# step: <name> <image> <working dir inside /data> <script and arguments...>
step() {
    local name="$1" image="$2" workdir="$3"
    shift 3
    echo "=== $name ($(date '+%H:%M:%S'))"
    mark "$name" running
    if docker run --rm --name "fungiflow_${RUN}_${name}" \
            --user "$RUN_AS" -e HOME=/tmp -e THREADS="$THREADS" \
            "${MOUNT[@]}" -w "$workdir" \
            "$image" bash "$@"; then
        mark "$name" done
    else
        mark "$name" failed
        echo "Step $name failed. Its log is in $RUN_DIR/logs/." >&2
        exit 1
    fi
}

: > "$STATUS"
echo "Run $RUN: ${#FASTQS[@]} barcodes, $THREADS threads"

if [[ ! -s "$UNITE_FASTA" ]]; then
    echo "Unpacking UNITE"
    tar -xzf database/unite/*.tgz -C database/unite
fi
if [[ ! -s "$BLAST_DB.nsq" ]]; then
    step 01 fungiflow-remastered /data scripts/01_build_blast_db.sh "$UNITE_FASTA"
fi
if [[ ! -s "$KRAKEN_DB/hash.k2d" ]]; then
    rm -rf "$KRAKEN_DB"
    step 02 fungiflow-remastered /data scripts/02_build_kraken_db.sh "$UNITE_FASTA"
fi

WORK="/data/$RUN_DIR"
step 03 fungiflow-remastered  "$WORK" scripts/03_read_qc.sh
step 04 fungiflow-remastered  "$WORK" scripts/04_filter_orient.sh
step 05 fungiflow-remastered  "$WORK" scripts/05_kraken2.sh
step 06 fungiflow-ngspeciesid "$WORK" scripts/06_cluster_consensus.sh
step 07 fungiflow-remastered  "$WORK" scripts/07_itsx.sh
step 08 fungiflow-remastered  "$WORK" scripts/08_otus.sh
step 09 fungiflow-remastered  "$WORK" scripts/09_blast.sh
step 10 fungiflow-remastered  "$WORK" scripts/10_organisms.sh
step 11 fungiflow-remastered  "$WORK" scripts/11_align.sh
step 12 fungiflow-remastered  "$WORK" scripts/12_tree.sh

mark all done
echo "Finished. Results in $RUN_DIR/results/"
