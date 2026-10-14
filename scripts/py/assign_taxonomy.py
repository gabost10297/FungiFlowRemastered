#!/usr/bin/env python3
"""Assign a name and a confidence level to every OTU from its BLAST hits.

Usage: python3 scripts/assign_taxonomy.py <blast_hits.tsv> <otu_table.tsv> <out.tsv> [thresholds]

Confidence levels:
  high        identity and query coverage pass the high thresholds  -> species
  low         identity and query coverage pass the low thresholds   -> species
  genus       identity and query coverage pass the genus thresholds -> genus
  unassigned  hits too weak for any of the above
  no_hit      no BLAST hit at all
If another species is almost as good a match (ambiguous), the name is moved up
to the deepest rank both hits share.
"""
import argparse
import csv
from collections import defaultdict

RANK_CODES = {"k": "kingdom", "p": "phylum", "c": "class", "o": "order",
              "f": "family", "g": "genus", "s": "species"}
RANKS = list(RANK_CODES.values())
BLAST_COLUMNS = ["qseqid", "sseqid", "pident", "length", "qstart", "qend", "qlen",
                 "sstart", "send", "slen", "evalue", "bitscore", "qcovhsp"]


def parse_lineage(sseqid):
    """'Name|Acc|SH|reftype|k__Fungi;...;s__X_y' -> {'kingdom': 'Fungi', ..., 'species': 'X_y'}"""
    lineage = {}
    for part in sseqid.split("|")[-1].split(";"):
        code, _, name = part.partition("__")
        if code in RANK_CODES and name:
            lineage[RANK_CODES[code]] = name
    return lineage


def read_hits(path):
    """Return {otu: [hit, ...]} with each OTU's hits sorted from best to worst."""
    hits = defaultdict(list)
    with open(path) as f:
        for line in f:
            hit = dict(zip(BLAST_COLUMNS, line.rstrip("\n").split("\t")))
            for key in ("pident", "bitscore", "qcovhsp"):
                hit[key] = float(hit[key])
            for key in ("sstart", "send", "slen"):
                hit[key] = int(hit[key])
            # How much of the reference sequence the alignment covers.
            hit["ref_cov"] = 100 * (abs(hit["send"] - hit["sstart"]) + 1) / hit["slen"]
            hit["lineage"] = parse_lineage(hit["sseqid"])
            otu = hit["qseqid"].split(";")[0]
            hits[otu].append(hit)
    for otu_hits in hits.values():
        otu_hits.sort(key=lambda h: -h["bitscore"])
    return hits


def read_otu_table(path):
    """Return {otu: (total reads, number of samples with reads)}."""
    otus = {}
    with open(path) as f:
        next(f)   # skip the header line
        for line in f:
            fields = line.rstrip("\n").split("\t")
            counts = [int(x) for x in fields[1:]]
            otus[fields[0]] = (sum(counts), sum(1 for c in counts if c > 0))
    return otus


def confidence_level(hit, t):
    """Apply the thresholds to the best hit: return (confidence, rank)."""
    pid, qcov = hit["pident"], hit["qcovhsp"]
    if pid >= t.high_pident and qcov >= t.high_qcov:
        return "high", "species"
    if pid >= t.low_pident and qcov >= t.low_qcov:
        return "low", "species"
    if pid >= t.genus_pident and qcov >= t.genus_qcov:
        return "genus", "genus"
    return "unassigned", None


def shared_rank(a, b):
    """Deepest rank at which two lineages still have the same name, or None."""
    shared = None
    for rank in RANKS:
        if a.get(rank) and a.get(rank) == b.get(rank):
            shared = rank
        else:
            break
    return shared


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("blast_hits")
    p.add_argument("otu_table")
    p.add_argument("out")
    p.add_argument("--high-pident", type=float, default=99)
    p.add_argument("--high-qcov", type=float, default=92)
    p.add_argument("--low-pident", type=float, default=97)
    p.add_argument("--low-qcov", type=float, default=85)
    p.add_argument("--genus-pident", type=float, default=90)
    p.add_argument("--genus-qcov", type=float, default=75)
    p.add_argument("--ambiguous-gap", type=float, default=0.5)
    p.add_argument("--min-ref-cov", type=float, default=85)
    t = p.parse_args()

    hits = read_hits(t.blast_hits)
    otus = read_otu_table(t.otu_table)

    rows = []
    for otu, (reads, n_samples) in otus.items():
        row = {"otu": otu, "reads": reads, "samples": n_samples}
        otu_hits = hits.get(otu, [])
        if not otu_hits:
            row.update(confidence="no_hit", rank="", name="")
            rows.append(row)
            continue

        best = otu_hits[0]
        best_species = best["lineage"].get("species", "")
        # The best hit to a *different* species, if any.
        second = next((h for h in otu_hits[1:]
                       if h["lineage"].get("species", "") != best_species), None)

        confidence, rank = confidence_level(best, t)
        ambiguous = (second is not None
                     and best["pident"] - second["pident"] < t.ambiguous_gap)
        if ambiguous and rank:
            common = shared_rank(best["lineage"], second["lineage"])
            if common is None or RANKS.index(common) < RANKS.index(rank):
                rank = common

        row.update(
            confidence=confidence,
            rank=rank or "",
            name=best["lineage"].get(rank, "") if rank else "",
            best_species=best_species,
            pident=f"{best['pident']:.2f}",
            qcov=f"{best['qcovhsp']:.0f}",
            ref_cov=f"{best['ref_cov']:.0f}",
            second_species=second["lineage"].get("species", "") if second else "",
            second_pident=f"{second['pident']:.2f}" if second else "",
            ambiguous="yes" if ambiguous else "no",
            partial_ref="yes" if best["ref_cov"] < t.min_ref_cov else "no",
            **{r: best["lineage"].get(r, "") for r in RANKS[:-1]},
        )
        rows.append(row)

    rows.sort(key=lambda r: -r["reads"])
    columns = ["otu", "reads", "samples", "confidence", "rank", "name", "best_species",
               "pident", "qcov", "ref_cov", "second_species", "second_pident",
               "ambiguous", "partial_ref"] + RANKS[:-1]
    with open(t.out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t", restval="")
        writer.writeheader()
        writer.writerows(rows)

    total_reads = sum(r["reads"] for r in rows)
    print(f"{'confidence':<12}{'OTUs':>6}{'reads':>9}{'% reads':>9}")
    for level in ["high", "low", "genus", "unassigned", "no_hit"]:
        group = [r for r in rows if r["confidence"] == level]
        n_reads = sum(r["reads"] for r in group)
        print(f"{level:<12}{len(group):>6}{n_reads:>9}{100 * n_reads / total_reads:>8.1f}%")
    print(f"ambiguous: {sum(r.get('ambiguous') == 'yes' for r in rows)} OTUs, "
          f"partial reference: {sum(r.get('partial_ref') == 'yes' for r in rows)} OTUs")

    print("\nMost abundant OTUs:")
    for r in rows[:25]:
        print(f"  {r['reads']:>7}  {r['samples']} samples  {r['confidence']:<10} "
              f"{r['rank']:<8} {r['name']}")


if __name__ == "__main__":
    main()