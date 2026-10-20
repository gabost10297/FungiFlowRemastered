#!/usr/bin/env python3
"""Compare every sequence of an organism with the organism's representative.

Only the stretch both sequences cover is compared, so a sequence that is
longer at one end is not counted as different; its extra length is reported
separately. Every column in that stretch where the two differ (another base,
or a gap in one of them) counts as one difference.
"""
import argparse
import csv
from pathlib import Path


def read_fasta(path):
    """Return {id: (header, sequence)}; the id is the header up to the first space."""
    seqs = {}
    name = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = (line[1:], [])
            elif name:
                seqs[name][1].append(line.upper())
    return {name: (header, "".join(parts)) for name, (header, parts) in seqs.items()}


def covered(seq):
    """First and last alignment column that is not a gap."""
    bases = [i for i, c in enumerate(seq) if c != "-"]
    return bases[0], bases[-1]


def compare(a, b):
    """Return (columns compared, differences) for two aligned sequences."""
    start = max(covered(a)[0], covered(b)[0])
    end = min(covered(a)[1], covered(b)[1])
    compared = differences = 0
    for x, y in zip(a[start:end + 1], b[start:end + 1]):
        if x == "-" and y == "-":
            continue
        compared += 1
        if x != y:
            differences += 1
    return compared, differences


def header_value(header, key):
    """'id otu=X reads=12' -> value of key"""
    for field in header.split()[1:]:
        k, _, v = field.partition("=")
        if k == key:
            return v
    return ""


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("organisms", help="organisms.tsv from script 10")
    parser.add_argument("aln_dir", help="folder with NN_<organism>.aln.fasta; tables are written here")
    parser.add_argument("--same", type=float, default=98.5,
                        help="identity treated as the same ITS (default: the OTU threshold)")
    parser.add_argument("--length-diff", type=int, default=20,
                        help="report sequences this much longer or shorter than the representative")
    args = parser.parse_args()
    out = Path(args.aln_dir)

    with open(args.organisms) as f:
        representative = {int(r["position"]): (r["organism"], r["representative"])
                          for r in csv.DictReader(f, delimiter="\t")}

    seq_rows, summary_rows = [], []
    for path in sorted(out.glob("*.aln.fasta")):
        organism, rep = representative[int(path.name[:2])]
        seqs = read_fasta(path)
        rep_seq = seqs[rep][1]
        rep_len = len(rep_seq.replace("-", ""))

        rows = []
        for name, (header, seq) in seqs.items():
            compared, differences = compare(rep_seq, seq)
            length = len(seq.replace("-", ""))
            rows.append({"name": name, "sample": name.split("_")[0],
                         "otu": header_value(header, "otu"),
                         "reads": int(header_value(header, "reads")),
                         "length": length, "length_diff": length - rep_len,
                         "differences": differences,
                         "identity": 100 * (compared - differences) / compared})
        rows.sort(key=lambda r: -r["reads"])
        for r in rows:
            seq_rows.append([organism, r["name"], r["sample"], r["otu"], r["reads"],
                             "yes" if r["name"] == rep else "no", r["length"], r["length_diff"],
                             r["differences"], f"{r['identity']:.2f}"])

        others = [r for r in rows if r["name"] != rep]
        reads = sum(r["reads"] for r in rows)
        same_reads = sum(r["reads"] for r in rows if r["identity"] >= args.same)
        summary_rows.append([
            organism, rep, len(rows), len({r["sample"] for r in rows}),
            f"{min(r['identity'] for r in others):.2f}" if others else "",
            sum(1 for r in others if r["identity"] < args.same),
            f"{100 * same_reads / reads:.1f}",
            sum(1 for r in others if abs(r["length_diff"]) > args.length_diff),
        ])

    with open(out / "vs_representative.tsv", "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["organism", "sequence", "sample", "otu", "reads", "representative",
                         "length", "length_diff", "differences", "identity"])
        writer.writerows(seq_rows)

    with open(out / "summary.tsv", "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["organism", "representative", "sequences", "samples", "min_identity",
                         f"below_{args.same:g}", f"pct_reads_at_least_{args.same:g}",
                         f"length_diff_over_{args.length_diff}"])
        writer.writerows(summary_rows)

    print(f"Each sequence compared with its organism's representative "
          f"(same ITS: >= {args.same:g}% identity).\n")
    print(f"{'organism':<26}{'seqs':>5}{'samples':>8}{'min %':>8}{'< same':>7}"
          f"{'% reads same':>13}{'length off':>11}")
    for organism, _, n, n_samples, low, below, pct, length_off in summary_rows:
        print(f"{organism:<26}{n:>5}{n_samples:>8}{low:>8}{below:>7}{pct:>13}{length_off:>11}")


if __name__ == "__main__":
    main()