#!/usr/bin/env python3
"""
gene_level_summary.py

Builds the gene-level summary of transfer and retention for plastid- and
mitochondrial-derived insertions, from the corrected combined Ka/Ks tables
and the BH-corrected dcNUxT tables.

Reports per gene:
    cNUxTs      coding insertions with a calculable Ka/Ks
    n_spp_c     species contributing cNUxTs
    dcNUxTs     of those, diverged (Ks > 0.1)
    ret%        dcNUxTs / cNUxTs - how often a transferred copy is old
    sel_dc      dcNUxTs under significant purifying selection (BH q < 0.05)
    sel%        sel_dc / dcNUxTs - how often an old copy is constrained
    n_spp_sel   independent species contributing selected dcNUxTs
    mean_len    mean aligned length of the cNUxTs (bp)
    gene_len    reference gene length, if gene_lengths.tsv is available

Usage:
    python3 gene_level_summary.py [--min-sel 1] [--latex]
"""

import os
import re
import sys
import pandas as pd

BASE = "/vol/storage/5T_disk/organelles"

NUPT_ALL = os.path.join(BASE, "all_species_nupt_kaks_fixed_combined.tsv")
NUMT_ALL = os.path.join(BASE, "all_species_numt_kaks_fixed_combined.tsv")
DCNUPT = "dcnupt_final_bh.tsv"
DCNUMT = "dcnumt_final_bh.tsv"
GENE_LENGTHS = "gene_lengths.tsv"

GENE_SYNONYMS = {f"nd{i}": f"nad{i}" for i in range(1, 10)}
GENE_SYNONYMS["nd4l"] = "nad4l"
GENE_SYNONYMS["cytb"] = "cob"      # same gene under two standard names


def extract_gene(seq_id):
    base = str(seq_id).replace(".mafft", "").replace("_REV", "")
    gene = base.rsplit("_", 1)[-1].lower()
    return GENE_SYNONYMS.get(gene, gene)


def summarise(all_path, dc_path, label, gene_lengths=None):
    if not os.path.exists(all_path):
        print(f"  [MISSING] {all_path}", file=sys.stderr)
        return None
    if not os.path.exists(dc_path):
        print(f"  [MISSING] {dc_path}", file=sys.stderr)
        return None

    alldf = pd.read_csv(all_path, sep="\t")
    dc = pd.read_csv(dc_path, sep="\t")

    alldf["gene"] = alldf["Sequence"].apply(extract_gene)
    dc["gene"] = dc["Sequence"].apply(extract_gene)
    alldf["Length"] = pd.to_numeric(alldf["Length"], errors="coerce")

    c = alldf.groupby("gene").agg(
        cNUxTs=("Sequence", "size"),
        n_spp_c=("Species", "nunique"),
        mean_len=("Length", "mean"),
    )

    sig_col = "Significant_BH_q05" if "Significant_BH_q05" in dc.columns else None
    if sig_col is None:
        raise SystemExit("dcNUxT table lacks Significant_BH_q05; rerun "
                         "bh_correction_dcnuxt.py")

    d = dc.groupby("gene").agg(
        dcNUxTs=("Sequence", "size"),
        sel_dc=(sig_col, "sum"),
        n_spp_dc=("Species", "nunique"),
    )
    sel_only = dc[dc[sig_col]]
    d["n_spp_sel"] = sel_only.groupby("gene")["Species"].nunique()
    d["n_spp_sel"] = d["n_spp_sel"].fillna(0).astype(int)

    t = c.join(d, how="left")
    for col in ("dcNUxTs", "sel_dc", "n_spp_dc", "n_spp_sel"):
        t[col] = t[col].fillna(0).astype(int)

    t["ret_pct"] = 100 * t["dcNUxTs"] / t["cNUxTs"]
    t["sel_pct"] = 100 * t["sel_dc"] / t["dcNUxTs"].replace(0, pd.NA)

    if gene_lengths is not None:
        t = t.join(gene_lengths, how="left")

    t = t.sort_values(["sel_dc", "dcNUxTs", "cNUxTs"], ascending=False)
    t.insert(0, "class", label)
    return t.reset_index()


def latex_block(t, label, min_sel):
    rows = t[t["sel_dc"] >= min_sel]
    out = []
    for _, r in rows.iterrows():
        gl = "---"
        if "length_bp" in r and pd.notna(r.get("length_bp")):
            gl = f"{r['length_bp']:.0f}"
        selpct = "---" if pd.isna(r["sel_pct"]) else f"{r['sel_pct']:.1f}"
        out.append(
            f"\\textit{{{r['gene']}}} & {r['cNUxTs']:.0f} & {r['n_spp_c']:.0f} & "
            f"{r['dcNUxTs']:.0f} & {r['ret_pct']:.1f} & {r['sel_dc']:.0f} & "
            f"{selpct} & {r['n_spp_sel']:.0f} & {r['mean_len']:.0f} & {gl} \\\\"
        )
    return out


def main():
    argv = sys.argv[1:]
    latex = "--latex" in argv
    if latex:
        argv.remove("--latex")
    min_sel = 1
    if "--min-sel" in argv:
        i = argv.index("--min-sel")
        min_sel = int(argv[i + 1])

    gene_lengths = None
    if os.path.exists(GENE_LENGTHS):
        gl = pd.read_csv(GENE_LENGTHS, sep="\t")
        gl["gene"] = gl["gene"].str.lower().replace(GENE_SYNONYMS)
        gene_lengths = gl.groupby("gene")["length_bp"].mean()

    nupt = summarise(NUPT_ALL, DCNUPT, "NUPT", gene_lengths)
    numt = summarise(NUMT_ALL, DCNUMT, "NUMT", gene_lengths)

    for t, name in ((nupt, "NUPT"), (numt, "NUMT")):
        if t is None:
            continue
        t.to_csv(f"gene_summary_{name.lower()}.tsv", sep="\t", index=False)

        print(f"\n=== {name} ===")
        print(f"genes with at least one cNUxT: {len(t)}")
        print(f"genes with at least one dcNUxT: {(t['dcNUxTs'] > 0).sum()}")
        print(f"genes with at least one selected dcNUxT: {(t['sel_dc'] > 0).sum()}")
        print(f"total cNUxTs: {t['cNUxTs'].sum()}  "
              f"dcNUxTs: {t['dcNUxTs'].sum()}  selected: {t['sel_dc'].sum()}")

        show = t.head(20)[["gene", "cNUxTs", "n_spp_c", "dcNUxTs",
                           "ret_pct", "sel_dc", "sel_pct", "n_spp_sel",
                           "mean_len"]]
        print(show.to_string(index=False,
                             float_format=lambda x: f"{x:.1f}"))

        # transfer concentration: how much of the total sits in the top genes
        top5 = t.nlargest(5, "cNUxTs")["cNUxTs"].sum()
        print(f"top 5 genes account for {100*top5/t['cNUxTs'].sum():.1f}% "
              f"of all cNUxTs")

    if latex:
        print("\n" + "=" * 60)
        print("LATEX ROWS (genes with >= "
              f"{min_sel} selected dcNUxT)")
        print("=" * 60)
        for t, name in ((nupt, "NUPT"), (numt, "NUMT")):
            if t is None:
                continue
            print(f"\n%% --- {name} ---")
            for line in latex_block(t, name, min_sel):
                print(line)

    print("\nWrote gene_summary_nupt.tsv and gene_summary_numt.tsv")


if __name__ == "__main__":
    main()
