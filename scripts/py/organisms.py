#!/usr/bin/env python3
"""Group identified OTUs into organisms and count them per sample.

Only high and low confidence fungal OTUs count as organisms:
  - a species name (Cladosporium_herbarum) joins all OTUs with that name
  - a genus name (BLAST could not decide between species) becomes "Genus sp."
  - an unnamed reference (Pleosporales_sp) or a higher rank stays one organism
    per OTU, because two such OTUs may be two different unnamed species
"""
import argparse
import csv
from collections import defaultdict
from pathlib import Path

RANKS = ["kingdom", "phylum", "class", "order", "family", "genus"]
GROUPS = ["organism", "genus_only", "non_fungal", "unassigned"]


def read_fasta(path):
    """Return {id: sequence}; the id is the header up to the first ';'."""
    seqs = {}
    name = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split(";")[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line)
    return {name: "".join(parts) for name, parts in seqs.items()}


def read_tsv(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def read_otu_table(path):
    """Return the sample names and {otu: {sample: reads}}."""
    with open(path) as f:
        rows = list(csv.reader(f, delimiter="\t"))
    samples = rows[0][1:]
    table = {row[0]: dict(zip(samples, map(int, row[1:]))) for row in rows[1:]}
    return samples, table


def group_of(row):
    """Which of the four groups an OTU belongs to."""
    if row["confidence"] in ("unassigned", "no_hit"):
        return "unassigned"
    if row["kingdom"] != "Fungi":
        return "non_fungal"
    if row["confidence"] in ("high", "low"):
        return "organism"
    return "genus_only"


def organism_key(row):
    """Return (key, label). OTUs with the same key are one organism.

    Unnamed organisms get label None and are numbered later.
    """
    name, rank = row["name"], row["rank"]
    if rank == "species" and not name.endswith("_sp"):
        return name, name.replace("_", " ")
    if rank == "genus":
        return name, f"{name} sp."
    return row["otu"], None


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("taxonomy", help="taxonomy.tsv from script 09")
    parser.add_argument("otu_table", help="otu_table.tsv from script 08")
    parser.add_argument("otus_fasta", help="otus.fasta from script 08")
    parser.add_argument("its_fasta", help="all_its.fasta from script 08")
    parser.add_argument("members", help="members.tsv from script 08")
    parser.add_argument("out_dir")
    parser.add_argument("--top", type=int, default=10, help="how many organisms to align")
    args = parser.parse_args()

    out = Path(args.out_dir)
    (out / "top").mkdir(parents=True, exist_ok=True)
    taxonomy = read_tsv(args.taxonomy)
    samples, table = read_otu_table(args.otu_table)

    # 1. Put every OTU in a group; collect the organism OTUs by key.
    groups = {g: [0, 0] for g in GROUPS}   # group -> [OTUs, reads]
    organisms = {}
    for row in taxonomy:
        row["group"] = group_of(row)
        groups[row["group"]][0] += 1
        groups[row["group"]][1] += int(row["reads"])
        if row["group"] == "organism":
            key, label = organism_key(row)
            organisms.setdefault(key, {"label": label, "otus": []})["otus"].append(row)

    # 2. Reads per sample, and the representative: best identity, then coverage, then reads.
    for org in organisms.values():
        org["per_sample"] = {s: sum(table[o["otu"]][s] for o in org["otus"]) for s in samples}
        org["reads"] = sum(org["per_sample"].values())
        org["samples"] = sum(1 for n in org["per_sample"].values() if n > 0)
        org["rep"] = max(org["otus"], key=lambda o: (float(o["pident"]), float(o["qcov"]), int(o["reads"])))

    # 3. Number unnamed organisms per name, most reads first: "Pleosporales sp. 1".
    numbers = defaultdict(int)
    for org in sorted(organisms.values(), key=lambda o: -o["reads"]):
        if org["label"] is None:
            base = org["otus"][0]["name"].removesuffix("_sp")
            numbers[base] += 1
            org["label"] = f"{base} sp. {numbers[base]}"
        org["id"] = org["label"].replace(" ", "_").replace(".", "")

    ranked = sorted(organisms.values(), key=lambda o: (-o["reads"], o["label"]))

    # 4. Tables.
    write_tsv(out / "organisms.tsv",
              ["position", "organism", "otus", "reads", "samples", *samples,
               "representative", "pident", "qcov", "ref_cov", "confidence", *RANKS],
              ([i, o["label"], len(o["otus"]), o["reads"], o["samples"],
                *(o["per_sample"][s] for s in samples),
                o["rep"]["otu"], o["rep"]["pident"], o["rep"]["qcov"], o["rep"]["ref_cov"],
                o["rep"]["confidence"], *(o["rep"][r] for r in RANKS)]
               for i, o in enumerate(ranked, 1)))

    label_of = {otu["otu"]: o["label"] for o in ranked for otu in o["otus"]}
    write_tsv(out / "otu_groups.tsv",
              ["otu", "reads", "group", "organism", "confidence", "name"],
              ([r["otu"], r["reads"], r["group"], label_of.get(r["otu"], ""),
                r["confidence"], r["name"]] for r in taxonomy))

    per_sample = []
    for s in samples:
        present = [o for o in ranked if o["per_sample"][s] > 0]
        per_sample.append([s,
                           sum(counts[s] for counts in table.values()),
                           sum(o["per_sample"][s] for o in present),
                           len(present),
                           sum(1 for o in present if o["samples"] == 1)])
    write_tsv(out / "per_sample.tsv",
              ["sample", "reads", "organism_reads", "organisms", "only_in_this_sample"],
              per_sample)

    # 5. Sequences: one representative per organism, and all consensuses of the top ones.
    otu_seqs = read_fasta(args.otus_fasta)
    with open(out / "representatives.fasta", "w") as f:
        for o in ranked:
            rep = o["rep"]
            f.write(f">{o['id']} otu={rep['otu']} pident={rep['pident']} qcov={rep['qcov']}\n"
                    f"{otu_seqs[rep['otu']]}\n")

    its_seqs = read_fasta(args.its_fasta)
    members = read_tsv(args.members)
    for i, o in enumerate(ranked[:args.top], 1):
        otus = {otu["otu"] for otu in o["otus"]}
        with open(out / "top" / f"{i:02d}_{o['id']}.fasta", "w") as f:
            for m in members:
                if m["otu"] in otus:
                    f.write(f">{m['consensus']} otu={m['otu']} reads={m['reads']}\n"
                            f"{its_seqs[m['consensus']]}\n")

    # 6. Summary.
    total = sum(reads for _, reads in groups.values())
    print(f"{'group':<12}{'OTUs':>6}{'reads':>9}{'% reads':>9}")
    for g in GROUPS:
        n, reads = groups[g]
        print(f"{g:<12}{n:>6}{reads:>9}{100 * reads / total:>8.1f}%")

    single = sum(1 for o in ranked if o["samples"] == 1)
    print(f"\nOrganisms: {len(ranked)} ({single} found in one sample only)")

    print(f"\n{'sample':<12}{'reads':>8}{'organism reads':>16}{'organisms':>11}{'only here':>11}")
    for s, reads, org_reads, n, only in per_sample:
        print(f"{s:<12}{reads:>8}{org_reads:>16}{n:>11}{only:>11}")

    print(f"\nTop {args.top} organisms:")
    for i, o in enumerate(ranked[:args.top], 1):
        print(f"{i:>3}  {o['label']:<30}{o['reads']:>7} reads{o['samples']:>3} samples"
              f"{len(o['otus']):>3} OTUs")


if __name__ == "__main__":
    main()