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

import glob
import os
import sys

import pandas as pd

SPECIES_TSV = "final_species_summary_v2.tsv"
DATASET = "combined_dataset.csv"
ROOTS = ["downloads", "multi_organelles/downloads"]
OUTPUT = "gene_lengths.tsv"

GENE_SYNONYMS = {f"nd{i}": f"nad{i}" for i in range(1, 10)}
GENE_SYNONYMS["nd4l"] = "nad4l"
GENE_SYNONYMS["cytb"] = "cob"

CHECK_GENES = ["rps3", "rpl2", "rps4", "rps14", "ycf2", "psba", "rbcl",
               "psbd", "atpa", "atp1"]


def normalize_gene_name(raw):
    name = raw.lower()
    # product-named files: "ribulose-1,5-bisphosphate carboxylase_oxygenase
    # large subunit" (the "/" became "_") is rbcL
    if name.endswith("oxygenase large subunit"):
        return "rbcl"
    return GENE_SYNONYMS.get(name, name)


def fasta_len(path):
    n = 0
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                n += len(line.strip())
    return n


def accession_dirs(kind, acc):
    """Folders for an accession, versioned or not, in every root."""
    acc = str(acc).strip()
    base = acc.split(".")[0]
    out = []
    for r in ROOTS:
        for pat in (f"{r}/{kind}/{base}", f"{r}/{kind}/{base}.*"):
            out += [p for p in glob.glob(pat) if os.path.isdir(p)]
    return sorted(set(out))


def species_dirs(kind, species):
    return [p for r in ROOTS for p in glob.glob(f"{r}/{kind}/{species}")
            if os.path.isdir(p)]


def gene_records(dirs, species, source):
    rec = []
    for d in dirs:
        for f in glob.glob(os.path.join(d, "*.fasta")):
            fname = os.path.basename(f)
            if fname.startswith("all_"):
                continue
            length = fasta_len(f)
            if length > 0:
                gene = normalize_gene_name(fname[: -len(".fasta")])
                rec.append((species, source, gene, length))
    return rec


def main():
    for p in (SPECIES_TSV, DATASET):
        if not os.path.exists(p):
            sys.exit(f"missing {p}; run from the organelles base directory")

    species = set(pd.read_csv(SPECIES_TSV, sep="\t")["species"])
    d = pd.read_csv(DATASET)
    d = d[d.iloc[:, 0].isin(species)]

    records, no_pl, no_mt = [], [], []
    for sp, _, pl, mt in d.iloc[:, :4].itertuples(index=False):
        pl_dirs = [p for acc in str(pl).split(";") for p in accession_dirs("plastid_genes", acc)]
        mt = str(mt)
        if ";" in mt:
            mt_dirs = species_dirs("mito_genes_merged", sp)
        else:
            mt_dirs = accession_dirs("mito_genes", mt) or species_dirs("mito_genes_merged", sp)

        r_pl = gene_records(pl_dirs, sp, "plastid")
        r_mt = gene_records(mt_dirs, sp, "mito")
        if not r_pl:
            no_pl.append(sp)
        if not r_mt:
            no_mt.append(sp)
        records += r_pl + r_mt

    df = pd.DataFrame(records, columns=["species", "source", "gene", "length"])
    per_sp = df.groupby(["source", "gene", "species"])["length"].mean().reset_index()
    out = (per_sp.groupby(["source", "gene"])["length"]
           .agg(length_bp="mean", median_bp="median", n_species="size")
           .reset_index())
    out["length_bp"] = out["length_bp"].round(1)
    out["median_bp"] = out["median_bp"].round(1)
    out = out.sort_values(["source", "n_species"], ascending=[False, False])
    out[["gene", "source", "length_bp", "median_bp", "n_species"]].to_csv(
        OUTPUT, sep="\t", index=False)

    n_sp = len(d)
    print(f"species in dataset: {n_sp}")
    print(f"  with plastid genes: {n_sp - len(no_pl)}   with mito genes: {n_sp - len(no_mt)}")
    print(f"genes: plastid {(out.source == 'plastid').sum()}, mito {(out.source == 'mito').sum()}")
    shared = set(out.gene[out.source == "plastid"]) & set(out.gene[out.source == "mito"])
    print(f"gene names present in both organelles: {len(shared)}")
    print(f"wrote {OUTPUT}\n")

    chk = out[out.gene.isin(CHECK_GENES)].pivot(index="gene", columns="source",
                                                  values=["length_bp", "n_species"])
    print("check genes (length_bp and n_species per organelle):")
    print(chk.to_string())


if __name__ == "__main__":
    main()
