#!/usr/bin/env python3

# Convert UNITE database to Kraken2 format

# Writes into OUT_DIR:
#    taxonomy/nodes.dmp     taxid > parent taxid and rank
#    taxonomy/names.dmp     taxid > names
#    library.fasta          'kraken:taxid|N' in the header 

# Usage: python3 scripts/unite_to_kraken.py <unite.fasta> <out_dir>

import sys
from pathlib import Path

RANKS = {"k": "kingdom", "p": "phylum", "c": "class", "o": "order",
         "f": "family", "g": "genus", "s": "species"}


def is_placeholder(name):
    # Check if the name is a placeholder ("unidentified", "..._sp" etc.)
    return name == "unidentified" or name.endswith("_sp")

def parse_lineage(taxonomy):
    # 'k__Fungi;p__Ascomycota;...' > [('kingdom', 'Fungi'), ('phylum', 'Ascomycota'), ...]

    lineage = []

    for part in taxonomy.split(";"):
        code, _, name = part.partition("__")
        if code in RANKS and name and not is_placeholder(name):
            lineage.append((RANKS[code], name))
    return lineage


def main():

    if len(sys.argv) != 3:                                                          
        sys.exit("Usage: python3 unite_to_kraken.py <unite.fasta> <out_dir>")   #Expect two: unite fasta file and folder 

    unite_fasta = Path(sys.argv[1])
    out_dir = Path(sys.argv[2])
    (out_dir / "taxonomy").mkdir(parents=True, exist_ok=True)


    # The tree starts with only the root, which has taxid 1 and is its own parent
    nodes = [(1, 1, "no rank", "root")]  # taxid, parent_taxid, rank, name
    taxid_of = {}                        # parent_taxid, rank, name) > taxid
    n_seqs = 0                           # how many sequences written
    n_to_species = 0                     # how many of them reached species level


    # Go through UNITE line by line
    # fin  = UNITE input, read-only
    # fout = library.fasta - the kraken2 formatted
    with open(unite_fasta) as fin, open(out_dir / "library.fasta", "w") as fout:
        
        # Seqs not starting with '>' are copied in uppercase and cleaned
        for line in fin:
            if not line.startswith(">"):
                fout.write(line.strip().upper() + "\n")
                continue
            
            # Split on '|'
            fields = line[1:].strip().split("|")
            if len(fields) < 5:
                sys.exit(f"Unexpected header: {line.strip()}")
            accession, taxonomy = fields[1], fields[-1]

    
            taxid = 1
            rank = "no rank"
            for rank, name in parse_lineage(taxonomy):
                    key = (taxid, rank, name)
                    # First time the taxon is under this parent - give it a taxid:
                    if key not in taxid_of:
                        taxid_of[key] = len(nodes) + 1
                        nodes.append((taxid_of[key], taxid, rank, name))
                    # Move one level down: this taxon is the parent of the next one.
                    taxid = taxid_of[key] # returns the number saved the first time
            
            n_seqs += 1
            if rank == "species":
                n_to_species += 1
            fout.write(f">{accession}_{n_seqs}|kraken:taxid|{taxid}\n")

    # Write taxonomy tree in NCBI style
    # Columns are separated by "\t|\t".
    # nodes.dmp: taxid | parent taxid | rank
    with open(out_dir / "taxonomy" / "nodes.dmp", "w") as f:
        for taxid, parent, rank, _ in nodes:
            f.write(f"{taxid}\t|\t{parent}\t|\t{rank}\t|\t\t|\n")
    with open(out_dir / "taxonomy" / "names.dmp", "w") as f:
        for taxid, _, _, name in nodes:
            f.write(f"{taxid}\t|\t{name}\t|\t\t|\tscientific name\t|\n")

    print(f"Sequences written:           {n_seqs}")
    print(f"Sequences labelled to species: {n_to_species}")
    print(f"Taxa in the tree:            {len(nodes)}")
    for rank in RANKS.values():
        print(f"  {rank:<10} {sum(1 for n in nodes if n[2] == rank)}")

if __name__ == "__main__":
    main()