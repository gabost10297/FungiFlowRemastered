#!/usr/bin/env bash
# Cluster each sample's oriented reads with NGSpeciesID and polish one consensus per cluster with medaka. Runs in the fungiflow-ngspeciesid container.
# Usage: bash scripts/06_cluster_consensus.sh
#
# Output in results/06_clusters/:
#   <sample>.consensus.fasta   one consensus per cluster, header: >sample_clID;size=READS
#   <sample>.ngspeciesid.log   NGSpeciesID's own log (read error rates etc.)
#   <sample>.reads.tsv         reads available and reads used
#   summary.tsv                all samples together

set -euo pipefail
source config.sh

IN_DIR="results/04_filter_orient"
OUT_DIR="results/06_clusters"
LOG="logs/06_cluster_consensus.log"

shopt -s nullglob

FASTQS=("$IN_DIR"/*.oriented.fastq.gz)
if (( ${#FASTQS[@]} == 0 )); then
    echo "Error: no oriented reads in $IN_DIR. Run script 04 first." >&2
    exit 1
fi
mkdir -p "$OUT_DIR" logs

{
    echo "Settings: $CLUSTER_READS reads per sample (0 = all), clusters >= $MIN_CLUSTER_SIZE reads"
    echo "Medaka model: $MEDAKA_MODEL, $THREADS threads"

    for fq in "${FASTQS[@]}"; do
        sample=$(basename "$fq" .oriented.fastq.gz)
        final="$OUT_DIR/$sample.consensus.fasta"
        work="$OUT_DIR/work_$sample"

        echo ">>> $sample"
        if [[ -s "$final" ]]; then
            echo "already done, skipping"
            continue
        fi

        # Start clean: remove leftovers of an interrupted run.
        rm -rf "$work"
        mkdir -p "$work"
        zcat "$fq" > "$work/reads.fastq"

        available=$(( $(wc -l < "$work/reads.fastq") / 4 ))
        if (( CLUSTER_READS > 0 && CLUSTER_READS < available )); then
            used=$CLUSTER_READS
        else
            used=$available
        fi
        printf "%s\t%s\t%s\n" "$sample" "$available" "$used" > "$OUT_DIR/$sample.reads.tsv"

        # NGSpeciesID keeps clusters larger than (abundance_ratio x reads used).
        ratio=$(awk -v m="$MIN_CLUSTER_SIZE" -v n="$used" 'BEGIN { printf "%.6f", m / n }')
        echo "clustering $used of $available reads (abundance ratio $ratio)"

        NGSpeciesID --ont --consensus --medaka \
                    --medaka_model "$MEDAKA_MODEL" \
                    --fastq "$work/reads.fastq" \
                    --outfolder "$work/ngspeciesid" \
                    --sample_size "$used" \
                    --seed 1 \
                    --q "$MIN_QUAL" \
                    --abundance_ratio "$ratio" \
                    --max_seqs_for_consensus "$MAX_SEQS_CONSENSUS" \
                    --sample_name "$sample" \
                    --t "$THREADS"

        # Collect all consensuses into one file with short headers, e.g.
        #   >barcode25_consensus_cl_id_7_total_supporting_reads_1150
        # becomes
        #   >barcode25_cl7;size=1150
        # (";size=" is the format vsearch reads cluster sizes from).
        # /dev/null keeps cat from waiting for keyboard input if no consensus exists.
        
        consensus=("$work"/ngspeciesid/medaka_cl_id_*/consensus.fasta)
        echo "consensus sequences: ${#consensus[@]}"
        cat "${consensus[@]}" /dev/null |
            awk -v s="$sample" '
                /^>/ { match($0, /cl_id_[0-9]+/)
                       id = substr($0, RSTART + 6, RLENGTH - 6)
                       size = $0
                       sub(/.*_total_supporting_reads_/, "", size)
                       print ">" s "_cl" id ";size=" size
                       next }
                { print }' > "$final.tmp"
        mv "$final.tmp" "$final"

        cp "$work/ngspeciesid/logfile.txt" "$OUT_DIR/$sample.ngspeciesid.log"
        rm -rf "$work"
        echo "done in $SECONDS seconds since start"
    done

    echo "Summary:"
    printf "sample\treads_available\treads_used\tconsensuses\treads_in_consensuses\tpct_in_consensuses\n" \
        > "$OUT_DIR/summary.tsv"
    for fq in "${FASTQS[@]}"; do
        sample=$(basename "$fq" .oriented.fastq.gz)
        [[ -e "$OUT_DIR/$sample.reads.tsv" ]] || continue
        read -r _ available used < "$OUT_DIR/$sample.reads.tsv"
        awk -F'size=' -v s="$sample" -v a="$available" -v u="$used" '
            /^>/ { n++; r += $2 }
            END { printf "%s\t%d\t%d\t%d\t%d\t%.1f\n", s, a, u, n, r, 100 * r / u }
        ' "$OUT_DIR/$sample.consensus.fasta" >> "$OUT_DIR/summary.tsv"
    done
    cat "$OUT_DIR/summary.tsv"
} 2>&1 | tee "$LOG"