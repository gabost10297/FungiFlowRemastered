# FungiFlow settings. Every script reads this with: source config.sh

THREADS="${THREADS:-$(nproc)}"

# Reference databases
UNITE_FASTA="database/unite/sh_general_release_dynamic_all_19.02.2025.fasta"
BLAST_DB="database/blast/unite"
KRAKEN_DB="database/kraken2"

# Read filtering (script 04)
MIN_QUAL=10      # minimum average read quality
MIN_LEN=400      # shortest read kept, roughly a complete ITS1-5.8S-ITS2
MAX_LEN=1200     # longest read kept, removes artefacts

# Kraken2 classification (script 05)
KRAKEN_CONFIDENCE=0 # minimum confidence score for Kraken2 classification (0 = no extra filtering)
KRAKEN_MIN_HITS=3 # # minimum number of separate k-mer hits needed to classify a read

# Clustering and consensus (script 06)
MEDAKA_MODEL="r1041_e82_400bps_hac_v5.2.0"   # must match the basecalling model
MIN_CLUSTER_SIZE=20        # smallest cluster (in reads) that gets a consensus
MAX_SEQS_CONSENSUS=300     # reads used per consensus; more adds time, not accuracy
CLUSTER_READS=0            # reads per sample to cluster; 0 = all reads

# ITSx (script 07)
ITSX_EVALUE=1e-2   # default 1e-5 misses the short LSU stretch at the end of many consensuses