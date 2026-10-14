#!/usr/bin/env python3
"""
04_build_funnel_table.py -- Table 1, the insertion funnel, and the unique dcNUxT set

  All merged insertions                 sum of <cls>_merged_count (species summary)
  containing organellar coding seq.     cNUxT  = rows of all_species_<cls>_kaks.tsv
                                        (calculable Ks)                        % of row 1
  diverged (Ks > 0.1)                   dcNUxT                                 % of row 2
  unique after collapsing identical     dcNUxTs collapsed within species x gene
    copies                              on identical fragment sequence (AXT)   % of row 3
  under significant purifying selection Ka/Ks < 1 and BH q < 0.05, BH applied to
                                        the unique set, separately per class   % of row 4

Writes dcnuxt_unique.tsv (one row per unique dcNUxT, with q-values). This set
is used by 05_bh_correction_dcnuxt.py and by
05_dcnuxt_integration/dcnuxt_annotation_overlap.py --final-set.

Usage (from the working directory):
  python3 04_build_funnel_table.py final_species_summary_v2.tsv
"""
import os
import re
import sys

import numpy as np
import pandas as pd

KS_MIN, ALPHA = 0.1, 0.05
AXT = {"nupt": "nupt_input.axt", "numt": "numt_input.axt"}


def bh(p):
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    n = ok.sum()
    if n:
        order = np.argsort(p[ok])
        ranked = p[ok][order] * n / np.arange(1, n + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        tmp = np.empty(n); tmp[order] = np.minimum(ranked, 1)
        q[ok] = tmp
    return q


def read_axt(path):
    """id -> fragment sequence (second sequence of each AXT block)"""
    out = {}
    if not os.path.isfile(path):
        return out
    lines = [l.rstrip("\n") for l in open(path)]
    i = 0
    while i + 2 < len(lines):
        if lines[i].strip():
            out[lines[i].strip()] = lines[i + 2].strip().upper()
            i += 4
        else:
            i += 1
    return out


def gene_of(seq_id):
    """..._<gene>[_N][_REV].mafft -> gene (divergent-copy suffix _N dropped)"""
    s = re.sub(r"\.mafft$", "", seq_id)
    s = re.sub(r"_REV$", "", s)
    parts = s.split("_")
    if len(parts) > 2 and parts[-1].isdigit() and not re.search(r"\d+-\d+$", parts[-2]):
        parts = parts[:-1]                  # atp1_2 -> atp1
    return parts[-1]


def load_class(cls):
    allc = pd.read_csv(f"all_species_{cls}_kaks.tsv", sep="\t")
    allc["Ks"] = pd.to_numeric(allc["Ks"], errors="coerce")
    dc = allc[allc["Ks"] > KS_MIN].copy()
    seqs, miss = [], 0
    cache = {}
    for sp, q in zip(dc["Species"], dc["Sequence"]):
        if sp not in cache:
            cache[sp] = read_axt(os.path.join("species", sp, "results", AXT[cls]))
            if not cache[sp]:
                print(f"  {cls}: {sp}: {AXT[cls]} missing or empty")
        s = cache[sp].get(q)
        miss += s is None
        seqs.append(s if s is not None else q)      # no AXT -> never collapsed
    if miss:
        print(f"  {cls}: {miss} dcNUxT ids not found in an AXT file (kept, not collapsed)")
    dc["gene"] = dc["Sequence"].map(gene_of)
    dc["frag_seq"] = seqs
    return len(allc), dc


def main():
    summ = pd.read_csv(sys.argv[1] if len(sys.argv) > 1 else "final_species_summary_v2.tsv", sep="\t")
    table, keep = {}, []
    for cls in ("numt", "nupt"):
        total = pd.to_numeric(summ.get(f"{cls}_merged_count"), errors="coerce").sum()
        coding, dc = load_class(cls)
        uniq = dc.drop_duplicates(["Species", "gene", "frag_seq"]).copy()
        p = pd.to_numeric(uniq.get("P_Value_Fisher"), errors="coerce")
        uniq["q_BH"] = bh(p)
        kaks = pd.to_numeric(uniq.get("Ka_Ks"), errors="coerce")
        purif = (kaks < 1) & (uniq["q_BH"] < ALPHA)
        uniq["class"] = cls.upper()
        uniq["purifying_sig"] = purif
        keep.append(uniq)
        table[cls] = dict(total=total, coding=coding, dc=len(dc), uniq=len(uniq),
                          testable=int(p.notna().sum()), purif=int(purif.sum()))

    def pct(a, b, nd=1):
        return f"{100 * a / b:.{nd}f}" if b else "NA"

    m, p = table["numt"], table["nupt"]
    lines = [
        ("All merged insertions", f"{m['total']:,.0f}", "—", f"{p['total']:,.0f}", "—"),
        ("containing organellar coding sequence", f"{m['coding']:,.0f}", pct(m['coding'], m['total'], 2),
         f"{p['coding']:,.0f}", pct(p['coding'], p['total'], 2)),
        ("diverged (Ks > 0.1)", f"{m['dc']:,}", pct(m['dc'], m['coding']), f"{p['dc']:,}", pct(p['dc'], p['coding'])),
        ("unique after collapsing identical copies", f"{m['uniq']:,}", pct(m['uniq'], m['dc']),
         f"{p['uniq']:,}", pct(p['uniq'], p['dc'])),
        ("under significant purifying selection", f"{m['purif']:,}", pct(m['purif'], m['uniq']),
         f"{p['purif']:,}", pct(p['purif'], p['uniq'])),
    ]
    print(f"\n{'Category':42s} {'Mito n':>9s} {'%':>6s} {'Plastid n':>10s} {'%':>6s}")
    for r in lines:
        print(f"{r[0]:42s} {r[1]:>9s} {r[2]:>6s} {r[3]:>10s} {r[4]:>6s}")
    print(f"\n(testable = non-NA Fisher p among unique: NUMT {m['testable']}/{m['uniq']}, "
          f"NUPT {p['testable']}/{p['uniq']}; purifying % of testable: "
          f"{pct(m['purif'], m['testable'])} / {pct(p['purif'], p['testable'])})")
    print("\nLaTeX rows:")
    for r in lines:
        print(" & ".join(r).replace("%", "\\%") + r" \\")

    out = pd.concat(keep, ignore_index=True)
    cols = ["class", "Species", "gene", "Sequence", "Ka", "Ks", "Ka_Ks", "P_Value_Fisher", "q_BH", "purifying_sig"]
    out[[c for c in cols if c in out.columns]].to_csv("dcnuxt_unique.tsv", sep="\t", index=False)
    print("\nrows -> dcnuxt_unique.tsv")


if __name__ == "__main__":
    main()
