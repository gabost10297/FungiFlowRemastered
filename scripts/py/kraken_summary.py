#!/usr/bin/env python3

# Summarise Kraken2 reports into two tables.

#       classification.tsv  per sample: reads, % classified, % Fungi, % plants, % animals
#       genus_reads.tsv     reads per genus (rows) and sample (columns), most abundant first

# Usage: python3 scripts/kraken_summary.py <report_dir>

import sys
from pathlib import Path

KINGDOMS = ["Fungi", "Viridiplantae", "Metazoa"]

def read_report(path):
    # Yield (rank, name, clade_reads) for every line of a Kraken2 report.
    with open(path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            yield fields[-3], fields[-1].strip(), int(fields[1])



def main():

    if len(sys.argv) != 2:
        sys.exit("Usage: python3 kraken_summary.py <report_dir>")

    report_dir = Path(sys.argv[1])
    reports = sorted(report_dir.glob("*.report.txt"))

    if not reports:
        sys.exit(f"No *.report.txt files in {report_dir}")

    samples = []
    genus_reads = {}   # genus > {sample: reads}

    with open(report_dir / "classification.tsv", "w") as out:

        out.write("sample\treads\tpct_classified\t"
                  + "\t".join(f"pct_{k}" for k in KINGDOMS) + "\n")

        for report in reports:

            sample = report.name.removesuffix(".report.txt")
            samples.append(sample)

            counts = {}

            for rank, name, reads in read_report(report):
                if rank in ("U", "R"):
                    counts[rank] = reads
                elif rank == "K":
                    counts[name] = reads
                elif rank == "G":
                    genus_reads.setdefault(name, {})[sample] = reads

            total = counts.get("U", 0) + counts.get("R", 0)

            row = [sample, str(total), f"{100 * counts.get('R', 0) / total:.1f}"]
            row += [f"{100 * counts.get(k, 0) / total:.1f}" for k in KINGDOMS]

            out.write("\t".join(row) + "\n")

    genera = sorted(genus_reads, key=lambda g: -sum(genus_reads[g].values()))

    with open(report_dir / "genus_reads.tsv", "w") as out:

        out.write("genus\ttotal\t" + "\t".join(samples) + "\n")

        for genus in genera:

            row = [genus_reads[genus].get(s, 0) for s in samples]
            out.write(f"{genus}\t{sum(row)}\t" + "\t".join(map(str, row)) + "\n")

    print((report_dir / "classification.tsv").read_text())
    print("Most abundant genera (reads, all samples):")

    for genus in genera[:15]:
        print(f"  {sum(genus_reads[genus].values()):>8}  {genus}")

if __name__ == "__main__":
    main()