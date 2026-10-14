#!/usr/bin/env python3
"""
nonhost_summary.py -- final non-host dcNUxT calls.

Reads (working directory):
  all_{nupt,numt}_query_categories_final.tsv   categories restricted to the unique
                                               dcNUxT set (nonhost_scaffold_composition.sh)
  organellar_fractions.tsv                     class, species, contig, length, fraction
                                               (nonhost_scaffold_composition.sh)

Exclusion rules (Methods):
  dcNUPT  excluded if its contig's organellar fraction >= 0.95
  dcNUMT  excluded if its contig has any organellar coverage (fraction > 0)
  contigs > 300 kb were recorded with fraction 0 (treated as nuclear)
  contigs that could not be tested (fraction NA) are reported as UNTESTED

Writes nonhost_candidates_final.tsv (one row per HIGH_CONF_NONHOST candidate with
contig, length, fraction and decision) and prints per-class category counts.

Usage: python3 06_nonhost/nonhost_summary.py
"""
import re
import sys

import pandas as pd

RULE = {"nupt": lambda f: f >= 0.95, "numt": lambda f: f > 0}
CONTIG = re.compile(r"_\d+-\d+_[^_]+(_REV)?$")

fr = pd.read_csv("organellar_fractions.tsv", sep="\t", header=None,
                 names=["class", "species", "contig", "contig_length", "organellar_fraction"],
                 na_values="NA")
out = []
for c in ("nupt", "numt"):
    try:
        cat = pd.read_csv(f"all_{c}_query_categories_final.tsv", sep="\t")
    except FileNotFoundError:
        sys.exit(f"missing all_{c}_query_categories_final.tsv; run nonhost_scaffold_composition.sh")
    n = cat.category.value_counts()
    print(f"dc{c.upper()} (unique set, n = {len(cat):,}): "
          + ", ".join(f"{k} {v:,}" for k, v in n.items()))
    hc = cat[cat.category == "HIGH_CONF_NONHOST"].copy()
    hc["contig"] = hc.qseqid.map(lambda q: CONTIG.sub("", q))
    hc = hc.merge(fr[fr["class"] == c].drop(columns="class"), on=["species", "contig"], how="left")
    hc.insert(0, "class", c.upper())
    hc["decision"] = [
        "UNTESTED" if pd.isna(f) else ("EXCLUDED_ORGANELLAR_CONTIG" if RULE[c](f) else "RETAINED")
        for f in hc.organellar_fraction]
    d = hc.decision.value_counts()
    print(f"  high-confidence non-host: {len(hc)}; retained {d.get('RETAINED', 0)}, "
          f"excluded {d.get('EXCLUDED_ORGANELLAR_CONTIG', 0)}, untested {d.get('UNTESTED', 0)}"
          f"; species with retained: {hc.loc[hc.decision == 'RETAINED', 'species'].nunique()}")
    out.append(hc)

pd.concat(out, ignore_index=True).to_csv("nonhost_candidates_final.tsv", sep="\t", index=False)
print("wrote nonhost_candidates_final.tsv")
