#!/usr/bin/env bash
#
# 02_prepare_background_alignments_nupt.sh
#
# Stages 1-2 of the selection analysis:
#   (1) Build a background protein alignment for each plastid gene across
#       all reference accessions.
#   (2) Add each dcNUPT cluster or singleton into its gene's background
#       alignment with MAFFT --add.
#
# Tree building, subtree extraction, codon alignment, codeml and result
# parsing are handled by 03_build_trees_iqtree.py, 04_rebuild_subtrees_
# and_alignments.py and 05_run_codeml_stage.py in this folder (see README).
#
# Usage: ./02_prepare_background_alignments_nupt.sh

set -euo pipefail
IFS=$'\n\t'

BASE_DIR="${BASE_DIR:-$PWD}"   # working directory; override with BASE_DIR=...
export BASE_DIR
PLASTID_GENES_DIR="$BASE_DIR/downloads/plastid_genes"
CLUSTER_DIR="$BASE_DIR/clustering_nupt/clusters"
SINGLES_DIR="$BASE_DIR/clustering_nupt/singles"
OUT_DIR="$BASE_DIR/selection_nupt"
THREADS=5
MAFFT_THREADS=4
GENETIC_CODE=11   # NCBI table 11, plant plastid

mkdir -p "$OUT_DIR"/background
mkdir -p "$OUT_DIR"/prank_alignments
mkdir -p "$OUT_DIR"/trees
mkdir -p "$OUT_DIR"/subtrees
mkdir -p "$OUT_DIR"/codon_alignments
mkdir -p "$OUT_DIR"/codeml
mkdir -p "$OUT_DIR"/results

########################################################
# STEP 2: Translate + align plastid genes across species
#         Remove sequences >20% shorter than longest
########################################################
echo "[2] Building background plastid gene alignments..."

python3 << 'EOF'
import os, glob, subprocess
from Bio import SeqIO
from Bio.Seq import Seq

PLASTID_DIR = os.path.join(os.environ["BASE_DIR"], "downloads/plastid_genes")
OUT_DIR = os.path.join(os.environ["BASE_DIR"], "selection_nupt/background")
MAFFT_THREADS = 4
GENETIC_CODE = 11
# Records labelled plastid in GenBank that are not plastid genomes (identified in
# the accession audit; PP476489.3 is a mitochondrial record and attached to the
# rpl16 trees on a ~1 substitution/site branch). Never used as background.
EXCLUDE_ACC = {"PP476489", "NC_016005", "OP649456", "OP154037"}
excluded = lambda d: os.path.basename(d).split(".")[0] in EXCLUDE_ACC

all_genes = set()
for acc_dir in glob.glob(os.path.join(PLASTID_DIR, "*")):
    if excluded(acc_dir):
        continue
    for f in glob.glob(os.path.join(acc_dir, "*.fasta")):
        all_genes.add(os.path.basename(f).replace(".fasta", ""))

print(f"Found {len(all_genes)} unique genes")

for gene in sorted(all_genes):
    nuc_out  = os.path.join(OUT_DIR, f"{gene}_nuc.fa")
    prot_out = os.path.join(OUT_DIR, f"{gene}_prot.fa")
    aln_out  = os.path.join(OUT_DIR, f"{gene}_prot_aln.fa")
    filt_out = os.path.join(OUT_DIR, f"{gene}_prot_aln_filt.fa")
    nuc_filt = os.path.join(OUT_DIR, f"{gene}_nuc_filt.fa")

    nuc_seqs = {}
    for acc_dir in glob.glob(os.path.join(PLASTID_DIR, "*")):
        if excluded(acc_dir):
            continue
        acc = os.path.basename(acc_dir)
        gene_file = os.path.join(acc_dir, f"{gene}.fasta")
        if not os.path.exists(gene_file):
            continue
        for rec in SeqIO.parse(gene_file, "fasta"):
            seq_str = str(rec.seq).upper().replace("-", "")
            trim = len(seq_str) % 3
            if trim:
                seq_str = seq_str[:-trim]
            if len(seq_str) < 30:
                continue
            nuc_seqs[acc] = seq_str

    if len(nuc_seqs) < 2:
        continue

    with open(nuc_out, 'w') as f:
        for acc, seq in nuc_seqs.items():
            f.write(f">{acc}\n{seq}\n")

    with open(prot_out, 'w') as f:
        for acc, seq in nuc_seqs.items():
            prot = str(Seq(seq).translate(table=GENETIC_CODE))
            f.write(f">{acc}\n{prot}\n")

    subprocess.run([
        "mafft", "--quiet", "--auto",
        "--thread", str(MAFFT_THREADS),
        prot_out
    ], stdout=open(aln_out, 'w'), stderr=subprocess.DEVNULL)

    import statistics

    if gene.isupper() or gene == "all_plastid" or gene.startswith("orf") \
        or gene.startswith("unknown") or " " in gene:
         continue

    seqs = {r.id: str(r.seq) for r in SeqIO.parse(aln_out, "fasta")}
    ungapped_lens = [len(s.replace("-", "")) for s in seqs.values()]
    median_len = statistics.median(ungapped_lens)
    cutoff = median_len * 0.8

    kept = {k: v for k, v in seqs.items()
            if len(v.replace("-", "")) >= cutoff}

    if len(kept) < 2:
        continue

    with open(filt_out, 'w') as f:
        for acc, seq in kept.items():
            f.write(f">{acc}\n{seq}\n")

    with open(nuc_filt, 'w') as f:
        for acc in kept:
            if acc in nuc_seqs:
                f.write(f">{acc}\n{nuc_seqs[acc]}\n")

    print(f"  {gene}: {len(nuc_seqs)} accessions -> {len(kept)} after filtering")

print("Background alignments done.")
EOF

########################################################
# STEP 3: Add NUPT sequences into background alignment
#         using MAFFT --add
########################################################
echo "[3] Aligning clusters/singles into background with MAFFT --add..."

# Background alignments are named from the plastid gene FASTA filenames,
# which keep their original capitalisation (atpA, rbcL, psbB...). Cluster
# names, however, come from 01_cluster_dcnupts.sh, which normalises gene names
# to lowercase. Resolve case-insensitively so clusters find their background.
resolve_bg() {
    local gene="$1"
    local suffix="$2"
    local exact="$OUT_DIR/background/${gene}_${suffix}.fa"
    if [[ -f "$exact" ]]; then
        echo "$exact"
        return
    fi
    find "$OUT_DIR/background" -maxdepth 1 -iname "${gene}_${suffix}.fa" 2>/dev/null | head -1
}
export -f resolve_bg

align_group() {
    local fa="$1"
    local type="$2"
    local filename
    filename=$(basename "$fa" .fa)

    if [[ "$type" == "cluster" ]]; then
        gene=$(echo "$filename" | sed 's/_cluster_.*//')
    else
        gene=$(echo "$filename" | sed 's/_REV$//' | awk -F'_' '{print $NF}')
    fi

    local bg_aln
    bg_aln=$(resolve_bg "$gene" "prot_aln_filt")
    local out_aln="$OUT_DIR/prank_alignments/${filename}"

    if [[ -z "$bg_aln" || ! -f "$bg_aln" ]]; then
        echo "  [SKIP] $filename: no background alignment for gene $gene"
        return
    fi

    mkdir -p "$out_aln"

    python3 - "$fa" "$out_aln/nupt_prot.fa" << 'PYEOF'
import sys
from Bio import SeqIO
from Bio.Seq import Seq

in_fa, out_fa = sys.argv[1], sys.argv[2]
written = 0
with open(out_fa, 'w') as out:
    for rec in SeqIO.parse(in_fa, "fasta"):
        seq = str(rec.seq).upper().replace("-", "")
        trim = len(seq) % 3
        if trim:
            seq = seq[:-trim]
        if len(seq) < 30:
            continue
        prot = str(Seq(seq).translate(table=11))
        if prot.endswith('*'):
            prot = prot[:-1]
        if '*' in prot:
            continue
        if len(prot) < 10:
            continue
        out.write(f">{rec.id}\n{prot}\n")
        written += 1
PYEOF

    if [[ ! -s "$out_aln/nupt_prot.fa" ]]; then
        echo "  [SKIP] $filename: no valid translations (all have stop codons)"
        return
    fi

    mafft --add "$out_aln/nupt_prot.fa" \
          --reorder \
          --thread "$MAFFT_THREADS" \
          --quiet \
          "$bg_aln" \
          > "$out_aln/combined_aln.fa" 2>/dev/null \
    || { echo "  WARNING: MAFFT --add failed for $filename"; return; }

    if [[ ! -s "$out_aln/combined_aln.fa" ]]; then
        echo "  WARNING: empty alignment for $filename"
        return
    fi

    echo "  Done: $filename ($gene)"
}

export -f align_group
export OUT_DIR MAFFT_THREADS

find "$CLUSTER_DIR" -name "*.fa" \
    | parallel -j 3 align_group {} cluster

find "$SINGLES_DIR" -name "*.fa" \
    | parallel -j 3 align_group {} single

echo "Alignments done: $(find "$OUT_DIR/prank_alignments" -name "combined_aln.fa" | wc -l)"

