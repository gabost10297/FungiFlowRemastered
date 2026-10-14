#!/usr/bin/env bash
# Remove chimeras, then group the full ITS sequences of all samples into OTUs
# (one per organism) and count the reads of every OTU in every sample.
# Usage: bash scripts/08_otus.sh
#
# Output in results/08_otus/:
#   chimeras/        per sample: chimeras found, sequences kept, summary.tsv
#   all_its.fasta    chimera-free ITS of all samples, header: >sample_clID;size=READS;sample=SAMPLE
#   otus.fasta       one representative ITS per OTU (its most abundant consensus)
#   otu_table.tsv    reads per OTU (rows) and sample (columns)
#   members.tsv      which consensus belongs to which OTU

set -euo pipefail
source config.sh

IN_DIR="results/07_itsx"
OUT_DIR="results/08_otus"
LOG="logs/08_otus.log"

shopt -s nullglob
FASTAS=("$IN_DIR"/*.its.fasta)
if (( ${#FASTAS[@]} == 0 )); then
    echo "Error: no ITS files in $IN_DIR. Run script 07 first." >&2
    exit 1
fi
mkdir -p "$OUT_DIR/chimeras" logs

{
    echo "Settings: OTU identity $OTU_IDENTITY"

    # Remove chimeras sample by sample (chimeras form within one PCR), then
    # combine all samples, adding ";sample=..." so vsearch knows the origin.
    : > "$OUT_DIR/all_its.fasta"
    printf "sample\tsequences\tborderline\tchimeras\treads\treads_borderline\treads_in_chimeras\n" > "$OUT_DIR/chimeras/summary.tsv"
    for fa in "${FASTAS[@]}"; do
        sample=$(basename "$fa" .its.fasta)
        if [[ ! -s "$fa" ]]; then
            echo "$sample: no ITS sequences, skipped"
            continue
        fi

        vsearch --uchime_denovo "$fa" \
        --sizein \
        --nonchimeras "$OUT_DIR/chimeras/$sample.ok.fasta" \
        --borderline "$OUT_DIR/chimeras/$sample.borderline.fasta" \
        --chimeras "$OUT_DIR/chimeras/$sample.chimeras.fasta" \
        --uchimeout "$OUT_DIR/chimeras/$sample.uchime.tsv" \
        --quiet

        # Keep non-chimeras and borderline sequences.
        cat "$OUT_DIR/chimeras/$sample.ok.fasta" "$OUT_DIR/chimeras/$sample.borderline.fasta" |
            awk -v s="$sample" '/^>/ { print $0 ";sample=" s; next } { print }' \
            >> "$OUT_DIR/all_its.fasta"

        # Count sequences and reads per file: 1 = kept, 2 = borderline, 3 = chimeras.
        count() { awk -F'size=' '/^>/ { n++; r += $2 } END { print n + 0, r + 0 }' "$1"; }
        read -r n_ok r_ok     < <(count "$OUT_DIR/chimeras/$sample.ok.fasta")
        read -r n_bord r_bord < <(count "$OUT_DIR/chimeras/$sample.borderline.fasta")
        read -r n_chim r_chim < <(count "$OUT_DIR/chimeras/$sample.chimeras.fasta")
        printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\n" "$sample" \
            $(( n_ok + n_bord + n_chim )) "$n_bord" "$n_chim" \
            $(( r_ok + r_bord + r_chim )) "$r_bord" "$r_chim" \
            >> "$OUT_DIR/chimeras/summary.tsv"

    done

    echo "Chimera check:"
    cat "$OUT_DIR/chimeras/summary.tsv"
    echo "ITS sequences after chimera removal: $(grep -c '^>' "$OUT_DIR/all_its.fasta")"

    vsearch --cluster_size "$OUT_DIR/all_its.fasta" \
            --id "$OTU_IDENTITY" \
            --strand plus \
            --sizein \
            --sizeout \
            --no_progress \
            --threads "$THREADS" \
            --centroids "$OUT_DIR/otus.fasta" \
            --uc "$OUT_DIR/otus.uc" \
            --otutabout "$OUT_DIR/otu_table.tsv"

    # From the .uc file: "S" lines are the OTU representatives, "H" lines the
    # sequences joined to one. Column 9 is the sequence, column 10 its OTU (or "*").
    printf "otu\tconsensus\tsample\treads\tidentity_to_otu\n" > "$OUT_DIR/members.tsv"
    awk -F'\t' '
        $1 == "S" || $1 == "H" {
            otu = ($1 == "S") ? $9 : $10
            ident = ($1 == "S") ? "100.0" : $4
            split($9, f, ";"); split(otu, o, ";")
            reads = f[2]; sub("size=", "", reads)
            sample = f[3]; sub("sample=", "", sample)
            print o[1] "\t" f[1] "\t" sample "\t" reads "\t" ident }
    ' "$OUT_DIR/otus.uc" >> "$OUT_DIR/members.tsv"

    echo "OTUs: $(grep -c '^>' "$OUT_DIR/otus.fasta")"
    echo "Most abundant OTUs (reads, number of samples):"
    awk -F'\t' '
        NR == 1 { next }
        { total = 0; n = 0
          for (i = 2; i <= NF; i++) { total += $i; if ($i > 0) n++ }
          print total "\t" n "\t" $1 }
    ' "$OUT_DIR/otu_table.tsv" | sort -k1,1nr | awk 'NR <= 20'
} 2>&1 | tee "$LOG"