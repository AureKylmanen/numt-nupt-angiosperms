#!/usr/bin/env python3
"""
03_combine_kaks.py -- combined Ka/Ks tables of all species, one per class.

Reads species/<SP>/results/<SP>_{nupt,numt}_kaks.csv (01_*_gene_alignment.sh,
deduplicated by 02_dedup_kaks_rows.py) and writes

  all_species_nupt_kaks.tsv   all_species_numt_kaks.tsv
      Species + the KaKs_Calculator columns, calculable rows only (numeric Ks).
      Column names are normalised: Ka/Ks -> Ka_Ks, P-Value(Fisher) ->
      P_Value_Fisher, parenthetical units dropped (Fold-Sites(0:2:4) -> Fold_Sites).
  coding_species_status.tsv
      species, class, n_genes in the organelle gene set, status

Species whose organelle gene set has fewer than MIN_GENES annotated genes are
left out of the coding analyses (status 'no annotated gene set'): with only a
handful of genes the Ka/Ks search is not comparable to the other species. The
gene set is the folder the alignment scripts used (plastid_genes/<PLASTID>,
mito_genes/<MITO>, or mito_genes_merged/<SPECIES> for multi-chromosome
mitogenomes).

Usage (from the working directory):
  python3 03_combine_kaks.py combined_dataset.csv
"""
import glob
import os
import re
import sys

import pandas as pd

CSV = sys.argv[1] if len(sys.argv) > 1 else "combined_dataset.csv"
MIN_GENES = 10
DROP_PAREN = ("Fold-Sites", "Fold-S-Substitutions", "Fold-N-Substitutions",
              "Substitution-Rate-Ratio", "GC(")


def norm(c):
    if any(c.startswith(x) for x in DROP_PAREN):
        c = re.sub(r"\(.*\)", "", c)
    return re.sub(r"[^A-Za-z0-9]+", "_", c).strip("_")


def gene_dir(cls, sp, plastid, mito):
    if cls == "nupt":
        hits = sorted(glob.glob(f"downloads/plastid_genes/{plastid}*"))
        return hits[0] if hits else None
    if ";" not in mito:
        hits = sorted(glob.glob(f"downloads/mito_genes/{mito}*"))
        if hits:
            return hits[0]
    return f"downloads/mito_genes_merged/{sp}"


def n_genes(d):
    if not d or not os.path.isdir(d):
        return 0
    return sum(1 for f in os.listdir(d) if f.endswith(".fasta") and not f.startswith("all_"))


ds = pd.read_csv(CSV, dtype=str).fillna("")
status = []
for cls in ("nupt", "numt"):
    parts = []
    for _, r in ds.iterrows():
        sp = r["SPECIES"].strip()
        ng = n_genes(gene_dir(cls, sp, r["PLASTID"].strip(), r["MITO"].strip()))
        if ng < MIN_GENES:
            status.append((sp, cls, ng, "no annotated gene set"))
            continue
        p = f"species/{sp}/results/{sp}_{cls}_kaks.csv"
        if not os.path.isfile(p):
            status.append((sp, cls, ng, "no Ka/Ks table"))
            continue
        d = pd.read_csv(p, sep="\t")
        d.columns = [norm(c) for c in d.columns]
        d = d[pd.to_numeric(d.get("Ks"), errors="coerce").notna()] if not d.empty else d
        status.append((sp, cls, ng, f"{len(d)} coding rows"))
        if d.empty:
            continue
        d.insert(0, "Species", sp)
        parts.append(d)
    if not parts:
        x = pd.DataFrame(status, columns=["species", "class", "n_genes", "status"])
        x = x[x["class"] == cls]
        sys.exit(f"{cls}: no Ka/Ks rows found. Per species:\n" + x.to_string(index=False)
                 + f"\nExpected species/<SP>/results/<SP>_{cls}_kaks.csv (01_{cls}_gene_alignment.sh) "
                 f"and >= {MIN_GENES} gene FASTAs in downloads/{{plastid,mito}}_genes/<ACC>/ "
                 "or downloads/mito_genes_merged/<SP>/ (00_data/build_gene_sets.py)")
    out = pd.concat(parts, ignore_index=True)
    out.to_csv(f"all_species_{cls}_kaks.tsv", sep="\t", index=False)
    print(f"{cls}: {len(out)} coding rows from {out.Species.nunique()} species -> all_species_{cls}_kaks.tsv")

st = pd.DataFrame(status, columns=["species", "class", "n_genes", "status"])
st.to_csv("coding_species_status.tsv", sep="\t", index=False)
for cls in ("nupt", "numt"):
    x = st[st["class"] == cls]
    print(f"{cls}: {int((x.status == 'no annotated gene set').sum())} species without an annotated gene set "
          f"(< {MIN_GENES} genes), {int((x.status == 'no Ka/Ks table').sum())} without a Ka/Ks table")
