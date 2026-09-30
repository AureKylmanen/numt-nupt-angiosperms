#!/usr/bin/env python3
"""
04_build_gene_lengths.py

Extracts canonical reference gene lengths directly from the per-species
gene FASTA files (downloads/plastid_genes/*/*.fasta and
downloads/mito_genes/*/*.fasta), averaging across all species/accessions
that have a copy of each gene. This gives the actual gene length (a
fixed property of the gene), as opposed to the average matched-fragment
length (which varies per insertion and is NOT what we want here).

Usage:
    04_python3 build_gene_lengths.py
Output:
    gene_lengths.tsv   (gene, length_bp, n_copies_averaged, source)
"""

import os
import glob
import re

PLASTID_DIR = "downloads/plastid_genes"
MITO_DIR = "downloads/mito_genes"
OUTPUT = "gene_lengths.tsv"

# Known synonym mapping: nd1-nd9 and nad1-nad9 refer to the same genes
# under two different standard mitochondrial gene naming conventions
# (NCBI-style vs MITOS-style). This is NOT a typo/case issue - it needs
# explicit mapping. Extend this dict if you spot other known synonyms.
GENE_SYNONYMS = {f"nd{i}": f"nad{i}" for i in range(1, 10)}
GENE_SYNONYMS["nd4l"] = "nad4l"


def normalize_gene_name(raw_name):
    """Lowercase and apply known synonym mapping. Does NOT collapse
    suffix/isoform variants (e.g. atp1-1 vs atp1-2, nad4l-a vs nad4l-b) -
    those may represent genuinely distinct gene copies and should be
    reviewed manually before deciding whether to merge them."""
    name = raw_name.lower()
    name = GENE_SYNONYMS.get(name, name)
    return name


def fasta_seq_length(path):
    """Return the length of the (single) sequence in a fasta file."""
    length = 0
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                continue
            length += len(line.strip())
    return length


def is_clean_gene_name(raw_name):
    """Reject filenames that aren't clean gene symbols: entries with
    spaces, apostrophes, or explicit 'pseudo' labels (e.g. "pseudo
    rpl23", "pseudo 3'-rps12a", "rps 7") are annotation artifacts, not
    standard gene names, and would corrupt per-gene averages if pooled
    under inconsistent keys."""
    if " " in raw_name or "'" in raw_name:
        return False
    if "pseudo" in raw_name.lower():
        return False
    return True


def collect_gene_lengths(base_dir, exclude_name, source_label):
    gene_lengths = {}  # gene -> list of lengths across species
    skipped = []
    dirs = glob.glob(f"{base_dir}/*/")
    for d in dirs:
        for fasta in glob.glob(os.path.join(d, "*.fasta")):
            fname = os.path.basename(fasta)
            if fname == exclude_name:
                continue
            gene_raw = fname.replace(".fasta", "")
            if not is_clean_gene_name(gene_raw):
                skipped.append(gene_raw)
                continue
            gene = normalize_gene_name(gene_raw)
            length = fasta_seq_length(fasta)
            if length > 0:
                gene_lengths.setdefault(gene, []).append(length)
    if skipped:
        unique_skipped = sorted(set(skipped))
        print(f"  [{source_label}] Skipped {len(skipped)} files with non-standard names "
              f"({len(unique_skipped)} unique): {unique_skipped[:15]}"
              f"{'...' if len(unique_skipped) > 15 else ''}")
    return gene_lengths, source_label


def main():
    plastid_lengths, _ = collect_gene_lengths(PLASTID_DIR, "all_plastid.fasta", "plastid")
    mito_lengths, _ = collect_gene_lengths(MITO_DIR, "all_mito.fasta", "mito")

    rows = []
    for gene, lengths in plastid_lengths.items():
        avg_len = sum(lengths) / len(lengths)
        rows.append((gene, round(avg_len, 1), len(lengths), "plastid"))

    for gene, lengths in mito_lengths.items():
        avg_len = sum(lengths) / len(lengths)
        rows.append((gene, round(avg_len, 1), len(lengths), "mito"))

    with open(OUTPUT, "w") as f:
        f.write("gene\tlength_bp\tn_copies_averaged\tsource\n")
        for gene, length, n, source in sorted(rows, key=lambda x: -x[2]):
            f.write(f"{gene}\t{length}\t{n}\t{source}\n")

    print(f"Wrote {len(rows)} genes to {OUTPUT}")
    print(f"  Plastid genes: {len(plastid_lengths)}")
    print(f"  Mito genes: {len(mito_lengths)}")
    print("\nSample rows:")
    with open(OUTPUT) as f:
        for i, line in enumerate(f):
            if i > 10:
                break
            print(line.strip())


if __name__ == "__main__":
    main()
