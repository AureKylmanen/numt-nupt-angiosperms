#!/usr/bin/env python3
"""
build_species_summary.py -- one row per species, from the per-species result files.

Reads, for every species in the species CSV (minus EXCLUDE):
  species/<SP>/results/<SP>.<name>.txt     every single-value file written by
                                           numt_nupt_detection.sh becomes column <name>
                                           (genome_size, numt_count, numt_merged_count,
                                           numt_unmerged_prop, numt_merged_prop, numt_union_bp,
                                           numt_nupt_overlap_bp/_count, numt_nupt_shared_bp, ...)
  numts/ nupts/ merged BEDs                -> numt_merged_bp, nupt_merged_bp
  organelle FASTAs (same lookup as the
  detection script)                        -> plastid_length, mito_length, n_mito_seqs
  combined_dataset.csv                     -> assembly, plastid_accession, mito_accessions
  all_species_{nupt,numt}_kaks.tsv         -> cnupt_n, dcnupt_n, cnumt_n, dcnumt_n (if present;
                                           02_cnuxt_kaks/03_combine_kaks.py)

Writes final_species_summary_v2.tsv, the table every downstream script reads.

EXCLUDE: species whose nuclear assembly is substantially incomplete (Methods).

Usage (from the working directory):
  python3 01_detection/build_species_summary.py combined_dataset.csv
"""
import glob
import os
import sys

import pandas as pd

CSV = sys.argv[1] if len(sys.argv) > 1 else "combined_dataset.csv"
OUT = "final_species_summary_v2.tsv"
MERGE = 10000
KS_MIN = 0.1
EXCLUDE = {"Tilia_amurensis", "Tripidium_arundinaceum"}   # partial nuclear assemblies
DL = ["downloads", "multi_organelles/downloads"]


def find_fa(kind, acc):
    for d in DL:
        p = f"{d}/{kind}_fasta/{acc}.fasta"
        if os.path.getsize(p) if os.path.isfile(p) else 0:
            return p
    base = acc.rsplit(".", 1)[0]
    for d in DL:
        hits = sorted(glob.glob(f"{d}/{kind}_fasta/{base}.[0-9]*.fasta"))
        if hits:
            return hits[0]
    return None


def fasta_len(path):
    n = 0
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                n += len(line.strip())
    return n


def read_value(path):
    try:
        v = open(path).read().split()[0]
    except (OSError, IndexError):
        return None
    try:
        f = float(v)
        return int(f) if f.is_integer() and "." not in v and "e" not in v.lower() else f
    except ValueError:
        return v


def bed_bp(path):
    if not os.path.isfile(path):
        return None
    s = 0
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) >= 3:
                s += int(f[2]) - int(f[1])
    return s


ds = pd.read_csv(CSV, dtype=str).fillna("")
ds = ds[~ds.SPECIES.isin(EXCLUDE)]
rows = []
for _, r in ds.iterrows():
    sp = r.SPECIES.strip()
    res = f"species/{sp}/results"
    row = {"species": sp, "assembly": r.ASSEMBLY, "plastid_accession": r.PLASTID,
           "mito_accessions": r.MITO}
    for f in sorted(glob.glob(f"{res}/{sp}.*.txt")):
        key = os.path.basename(f)[len(sp) + 1:-4]
        if key.endswith("organelle_annotated_headers"):
            continue
        row[key] = read_value(f)
    for c in ("numt", "nupt"):
        row[f"{c}_merged_bp"] = bed_bp(f"species/{sp}/{c}s/{c}s.clean.out_AllFrags_merged{MERGE}.bed")
    pf = find_fa("plastid", r.PLASTID)
    row["plastid_length"] = fasta_len(pf) if pf else None
    mfs = [find_fa("mito", a) for a in r.MITO.split(";") if a]
    row["n_mito_seqs"] = len(mfs)
    row["mito_length"] = sum(fasta_len(m) for m in mfs) if mfs and all(mfs) else None
    rows.append(row)

s = pd.DataFrame(rows)
# names used by the downstream scripts
s = s.rename(columns={"numt_unmerged_prop": "numt_prop", "nupt_unmerged_prop": "nupt_prop",
                      "numt_nupt_overlap_bp": "overlap_bp", "numt_nupt_overlap_count": "overlap_count",
                      "numt_nupt_shared_bp": "shared_bp"})

for c in ("nupt", "numt"):
    f = f"all_species_{c}_kaks.tsv"
    if os.path.isfile(f):
        k = pd.read_csv(f, sep="\t", usecols=["Species", "Ks"])
        k["Ks"] = pd.to_numeric(k["Ks"], errors="coerce")
        s[f"c{c}_n"] = s.species.map(k.groupby("Species").size()).fillna(0).astype(int)
        s[f"dc{c}_n"] = s.species.map(k[k.Ks > KS_MIN].groupby("Species").size()).fillna(0).astype(int)

front = ["species", "assembly", "genome_size", "plastid_accession", "plastid_length",
         "mito_accessions", "n_mito_seqs", "mito_length",
         "numt_count", "numt_merged_count", "numt_merged_bp", "numt_prop", "numt_merged_prop",
         "nupt_count", "nupt_merged_count", "nupt_merged_bp", "nupt_prop", "nupt_merged_prop"]
cols = [c for c in front if c in s.columns] + [c for c in s.columns if c not in front]
s[cols].to_csv(OUT, sep="\t", index=False)

print(f"wrote {OUT}: {len(s)} species ({len(EXCLUDE)} excluded), {len(cols)} columns")
for c in ("numt", "nupt"):
    print(f"  {c.upper()}: {int(s[f'{c}_merged_count'].fillna(0).sum()):,} merged insertions, "
          f"{s[f'{c}_merged_bp'].fillna(0).sum() / 1e6:,.0f} Mb; median per species "
          f"{s[f'{c}_merged_count'].median():.0f} insertions, {s[f'{c}_merged_bp'].median() / 1e3:,.0f} kb")
miss = s.species[s.genome_size.isna()].tolist() if "genome_size" in s else s.species.tolist()
if miss:
    print(f"  no detection results for {len(miss)} species: {', '.join(miss[:10])}")
print(f"  missing organelle lengths: plastid {int(s.plastid_length.isna().sum())}, "
      f"mito {int(s.mito_length.isna().sum())}")
