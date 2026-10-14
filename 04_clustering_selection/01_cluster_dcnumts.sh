#!/usr/bin/env bash

# ============================================================
# Script Name: cluster_dcnumts.sh
# Description: 
#   Identifies groups of diverged coding NUMTs (dcNUMTs) that likely 
#   share a common origin, so that shared-ancestry insertions can be
#   analysed together as clades rather than as independent observations
#   in downstream selection analysis.
# Author: Aure Kylmänen
# Created: 24.08.2026
# ============================================================

set -euo pipefail
IFS=$'\n\t'

########################################################
# CONFIG
########################################################
BASE_DIR="${BASE_DIR:-$PWD}"   # working directory; override with BASE_DIR=...
export BASE_DIR
SPECIES_DIR="$BASE_DIR/species"
MITO_DB="/vol/storage/mt_database/mito_db"
KAKS_BIN="/vol/storage/software/KaKs_Calculator2.0/bin/Linux/KaKs_Calculator"
OUT_DIR="$BASE_DIR/clustering_numt"
THREADS=12
MAFFT_THREADS=4
KS_BETWEEN=0.1
IDENTITY_THRESH=98

mkdir -p "$OUT_DIR"/{clusters,singles,logs}

########################################################
# STEP 1: Collect dcNUMT sequences per gene
########################################################
echo "[1] Collecting dcNUMT sequences per gene across species..."

# Clear per_gene/ - collection below appends, so stale content would duplicate
rm -rf "$OUT_DIR/per_gene"
mkdir -p "$OUT_DIR/per_gene"

python3 << 'EOF'
import os, glob, re

BASE_DIR = os.environ["BASE_DIR"]
SPECIES_DIR = os.path.join(BASE_DIR, "species")
OUT_DIR = os.path.join(BASE_DIR, "clustering_numt")
KS_THRESH = 0.1

GENE_SYNONYMS = {f"nd{i}": f"nad{i}" for i in range(1, 10)}
GENE_SYNONYMS["nd4l"] = "nad4l"
GENE_SYNONYMS["atpa"] = "atp1"   # plastid-style name used for atp1 in some mitogenome annotations


def normalize_gene_name(raw):
    g = raw.lower()
    return GENE_SYNONYMS.get(g, g)


def parse_axt(axt_path):
    """Parse an AXT file into {sequence_id: fragment_sequence}.

    AXT block format written by the pipeline:
        <id line>
        <gene sequence>       <- plastid reference gene
        <fragment sequence>   <- NUMT insertion (this is what we want)
        <blank line>
    """
    seqs = {}
    with open(axt_path) as f:
        block = []
        for line in f:
            line = line.rstrip("\n")
            if line.strip() == "":
                if len(block) >= 3:
                    seq_id = block[0].strip()
                    frag = block[2].strip()
                    if seq_id and frag:
                        seqs[seq_id] = frag
                block = []
            else:
                block.append(line)
        if len(block) >= 3:
            seq_id = block[0].strip()
            frag = block[2].strip()
            if seq_id and frag:
                seqs[seq_id] = frag
    return seqs


kaks_files = sorted(glob.glob(
    os.path.join(SPECIES_DIR, "*/results/*_numt_kaks.csv")))
print(f"Found {len(kaks_files)} NUMT KaKs files")

n_collected = 0
n_missing_seq = 0
n_missing_axt = 0
n_redundant = 0
n_orf = 0
# Unnamed reading frames (orf..., hypothetical, locus tags such as gp005 or
# ZemaCp066, "unknown") cannot be matched across species, so they are left out
# of clustering and branch models (they are counted as coding sequence in Table 1).
ORF_TAG = re.compile(r"^(orf|hypothetical|gp\d|p\d{3}$|locus|unknown$|[A-Za-z]{2,6}(cp|mp)\d{2,}$)", re.I)
species_contributing = set()

seen_by_species_gene = {}

for kaks_file in kaks_files:
    species = kaks_file.split(os.sep)[-3]
    axt_file = os.path.join(SPECIES_DIR, species, "results", "numt_input.axt")

    if not os.path.exists(axt_file):
        n_missing_axt += 1
        continue

    seqs = parse_axt(axt_file)

    with open(kaks_file) as f:
        next(f)
        for line in f:
            parts = line.rstrip("\n").split('\t')
            if len(parts) < 4:
                continue
            seq_id_raw = parts[0].strip()
            ks = parts[3].strip()
            if ks == 'NA':
                continue
            try:
                ks_val = float(ks)
            except ValueError:
                continue
            if ks_val <= KS_THRESH:
                continue

            if seq_id_raw not in seqs:
                n_missing_seq += 1
                continue

            seq_id_clean = seq_id_raw.replace('.mafft', '').replace('_REV', '')
            gene = seq_id_clean.split('_')[-1]
            if not gene:
                continue
            gene = normalize_gene_name(gene)
            if ORF_TAG.match(gene):
                n_orf += 1
                continue

            key = (species, gene)
            seq_str = seqs[seq_id_raw]
            if seq_str in seen_by_species_gene.setdefault(key, set()):
                n_redundant += 1
                continue
            seen_by_species_gene[key].add(seq_str)

            gene_dir = os.path.join(OUT_DIR, "per_gene", gene)
            os.makedirs(gene_dir, exist_ok=True)
            out_fa = os.path.join(gene_dir, f"{gene}_all_dcnumts.fa")

            header = f"{species}__{seq_id_raw.replace('.mafft', '')}"
            with open(out_fa, 'a') as out:
                out.write(f">{header}\n{seq_str}\n")

            n_collected += 1
            species_contributing.add(species)

print(f"Collection complete.")
print(f"  Sequences collected: {n_collected}")
print(f"  Species contributing: {len(species_contributing)}")
print(f"  Redundant copies collapsed: {n_redundant}")
print(f"  Unnamed reading frames skipped: {n_orf}")
print(f"  Skipped (id not found in AXT): {n_missing_seq}")
print(f"  Species with no numt_input.axt: {n_missing_axt}")
EOF

echo "Genes collected: $(ls "$OUT_DIR/per_gene/" 2>/dev/null | wc -l)"
echo "Sequences per gene (top 10):"
for f in "$OUT_DIR"/per_gene/*/; do
    gene=$(basename "$f")
    n=$(grep -c '^>' "$f/${gene}_all_dcnumts.fa" 2>/dev/null || echo 0)
    echo "  $gene: $n"
done | sort -t: -k2 -rn | head -10

########################################################
# STEP 2: Pairwise Ks per gene
########################################################
echo "[2] Running pairwise Ks between dcNUMTs per gene..."

pairwise_ks_gene() {
    local gene_fa="$1"
    local gene
    gene=$(basename "$(dirname "$gene_fa")")
    local work_dir
    work_dir=$(dirname "$gene_fa")
    local n_seqs
    n_seqs=$(grep -c '^>' "$gene_fa" 2>/dev/null || echo 0)

    if (( n_seqs < 2 )); then
        echo "  [SKIP] $gene: only $n_seqs sequence(s)"
        return
    fi

    echo "  $gene: $n_seqs sequences"

    mafft --quiet --thread "$MAFFT_THREADS" --auto "$gene_fa" \
        > "$work_dir/${gene}_aligned.fa" 2>/dev/null

    seqkit seq -w 0 "$work_dir/${gene}_aligned.fa" \
        > "$work_dir/${gene}_aligned_1l.fa"

    python3 - "$work_dir/${gene}_aligned_1l.fa" \
              "$work_dir/${gene}_pairwise.axt" << 'PYEOF'
import sys
from itertools import combinations

in_fa, out_axt = sys.argv[1], sys.argv[2]

seqs = {}
header = None
with open(in_fa) as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            header = line[1:]
            seqs[header] = ''
        else:
            seqs[header] += line.upper()

def trim_pair(s1, s2):
    kept1, kept2 = [], []
    for a, b in zip(s1, s2):
        if a == '-' or b == '-':
            continue
        kept1.append(a)
        kept2.append(b)
    s1, s2 = ''.join(kept1), ''.join(kept2)
    trim = len(s1) % 3
    if trim:
        s1, s2 = s1[:-trim], s2[:-trim]
    return s1, s2

written = 0
with open(out_axt, 'w') as out:
    for n1, n2 in combinations(seqs.keys(), 2):
        s1, s2 = trim_pair(seqs[n1], seqs[n2])
        if len(s1) < 30:
            continue
        out.write(f"{n1}_vs_{n2}\n{s1}\n{s2}\n\n")
        written += 1

print(f"  Wrote {written} pairs for {out_axt}", flush=True)
PYEOF

    if [[ ! -s "$work_dir/${gene}_pairwise.axt" ]]; then
        echo "  [SKIP] $gene: no valid pairs"
        return
    fi

    "$KAKS_BIN" \
        -i "$work_dir/${gene}_pairwise.axt" \
        -o "$work_dir/${gene}_pairwise_kaks.csv" \
        -m LWL -c 1 2>/dev/null \
        || echo "  WARNING: KaKs failed for $gene"
}

export -f pairwise_ks_gene
export MAFFT_THREADS KAKS_BIN

find "$OUT_DIR/per_gene" -name "*_all_dcnumts.fa" \
    | parallel -j "$THREADS" pairwise_ks_gene {}

########################################################
# STEP 3: Cluster by pairwise Ks + plastid identity
#         then connected components
########################################################
echo "[3] Clustering..."

python3 << 'EOF'
import os, glob, subprocess
from collections import defaultdict

BASE_DIR = os.environ["BASE_DIR"]
OUT_DIR = os.path.join(BASE_DIR, "clustering_numt")
MITO_DB = "/vol/storage/mt_database/mito_db"
KS_THRESH = 0.1
IDENTITY_THRESH = 98.0

def connected_components(nodes, pairs):
    parent = {n: n for n in nodes}
    def find(x):
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]
    def union(x, y):
        parent[find(x)] = find(y)
    for a, b in pairs:
        union(a, b)
    groups = defaultdict(set)
    for node in nodes:
        groups[find(node)].add(node)
    return list(groups.values())

def get_mito_identities(fa_file):
    result = subprocess.run([
        'blastn', '-query', fa_file,
        '-db', MITO_DB,
        '-outfmt', '6 qseqid pident',
        '-max_target_seqs', '1',
        '-max_hsps', '1',
        '-evalue', '1e-10',
        '-num_threads', '4'
    ], capture_output=True, text=True)
    identities = {}
    for line in result.stdout.strip().split('\n'):
        if not line:
            continue
        parts = line.split('\t')
        if len(parts) >= 2:
            qid, pident = parts[0], float(parts[1])
            if qid not in identities or pident > identities[qid]:
                identities[qid] = pident
    return identities

os.makedirs(os.path.join(OUT_DIR, 'clusters'), exist_ok=True)
os.makedirs(os.path.join(OUT_DIR, 'singles'), exist_ok=True)

total_clusters = 0
total_singles = 0

for gene_dir in sorted(glob.glob(os.path.join(OUT_DIR, 'per_gene', '*'))):
    gene = os.path.basename(gene_dir)
    kaks_file = os.path.join(gene_dir, f'{gene}_pairwise_kaks.csv')
    fa_file   = os.path.join(gene_dir, f'{gene}_all_dcnumts.fa')

    if not os.path.exists(kaks_file) or not os.path.exists(fa_file):
        continue

    seqs = {}
    header = None
    with open(fa_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                header = line[1:]
                seqs[header] = ''
            elif header:
                seqs[header] += line
    all_seqs = set(seqs.keys())

    identities = get_mito_identities(fa_file)

    candidate_pairs = []
    with open(kaks_file) as f:
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            ks = parts[3]
            if ks == 'NA':
                continue
            try:
                ks_val = float(ks)
            except ValueError:
                continue
            if ks_val >= KS_THRESH:
                continue
            if '_vs_' not in parts[0]:
                continue

            seqA, seqB = parts[0].split('_vs_', 1)
            if seqA not in all_seqs or seqB not in all_seqs:
                continue

            idA = identities.get(seqA, 0)
            idB = identities.get(seqB, 0)
            if idA >= IDENTITY_THRESH or idB >= IDENTITY_THRESH:
                continue

            candidate_pairs.append((seqA, seqB))

    components = connected_components(all_seqs, candidate_pairs)
    clustered = set()

    for i, comp in enumerate(components):
        if len(comp) < 2:
            continue
        cluster_id = f"{gene}_cluster_{i+1}"
        out_fa = os.path.join(OUT_DIR, 'clusters', f'{cluster_id}.fa')
        with open(out_fa, 'w') as out:
            for sid in sorted(comp):
                if sid in seqs:
                    out.write(f'>{sid}\n{seqs[sid]}\n')
        clustered.update(comp)
        total_clusters += 1
        print(f"  Cluster {cluster_id}: {len(comp)} seqs")

    for sid in sorted(all_seqs - clustered):
        if sid not in seqs:
            continue
        out_fa = os.path.join(OUT_DIR, 'singles', f'{sid}.fa')
        with open(out_fa, 'w') as out:
            out.write(f'>{sid}\n{seqs[sid]}\n')
        total_singles += 1

print(f"\nTotal clusters: {total_clusters}")
print(f"Total singles:  {total_singles}")
EOF

echo "Done."
echo "Clusters: $(ls "$OUT_DIR/clusters/" | wc -l)"
echo "Singles:  $(ls "$OUT_DIR/singles/"  | wc -l)"
