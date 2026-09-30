#!/usr/bin/env python3
"""
02_bh_correction_dcnuxt.py

Filters the combined NUPT and NUMT Ka/Ks tables (from merge_kaks.py) to the
dcNUxT threshold (Ks > 0.1), applies Benjamini-Hochberg correction
to the Fisher's exact p-values SEPARATELY for NUPTs and NUMTs
(treated as two independent test families), and writes final
annotated tables with both raw and BH-adjusted (q-value) columns.

Usage:
    python3 02_bh_correction_dcnuxt.py \
        all_species_nupt_kaks.tsv \
        all_species_numt_kaks.tsv

Outputs:
    dcnupt_final_bh.tsv
    dcnumt_final_bh.tsv
    summary_evolutionary_pressure.tsv   (matches tab:evolutionary_pressure)
"""

import sys
import pandas as pd
import numpy as np

KS_THRESHOLD = 0.1
ALPHA = 0.05


def benjamini_hochberg(pvals):
    """Return BH-adjusted q-values for an array of p-values (NaNs preserved)."""
    pvals = np.asarray(pvals, dtype=float)
    n = np.sum(~np.isnan(pvals))
    qvals = np.full_like(pvals, np.nan)

    valid_idx = np.where(~np.isnan(pvals))[0]
    valid_pvals = pvals[valid_idx]

    order = np.argsort(valid_pvals)
    ranked_pvals = valid_pvals[order]
    ranks = np.arange(1, n + 1)

    raw_q = ranked_pvals * n / ranks
    # enforce monotonicity (standard BH step-up correction)
    q_monotone = np.minimum.accumulate(raw_q[::-1])[::-1]
    q_monotone = np.clip(q_monotone, 0, 1)

    # scatter back to original order
    qvals[valid_idx[order]] = q_monotone
    return qvals


def process_class(path, label):
    df = pd.read_csv(path, sep="\t")

    # Column names as written by merge_kaks.py (Ka_Ks, P_Value_Fisher, ...)
    df.columns = [c.strip() for c in df.columns]

    n_total = len(df)
    df["Ks"] = pd.to_numeric(df["Ks"], errors="coerce")
    df["Ka_Ks"] = pd.to_numeric(df["Ka_Ks"], errors="coerce")
    df["P_Value_Fisher"] = pd.to_numeric(df["P_Value_Fisher"], errors="coerce")

    # dcNUxT filter: Ks > 0.1 AND a calculable Ka/Ks (already true post-merge)
    dcnuxt = df[df["Ks"] > KS_THRESHOLD].copy()
    n_dcnuxt = len(dcnuxt)

    dcnuxt["P_BH"] = benjamini_hochberg(dcnuxt["P_Value_Fisher"].values)
    dcnuxt["Significant_raw_p05"] = dcnuxt["P_Value_Fisher"] < ALPHA
    dcnuxt["Significant_BH_q05"] = dcnuxt["P_BH"] < ALPHA

    n_sig_raw = dcnuxt["Significant_raw_p05"].sum()
    n_sig_bh = dcnuxt["Significant_BH_q05"].sum()
    n_testable = dcnuxt["P_Value_Fisher"].notna().sum()

    print(f"--- {label} ---")
    print(f"  Total calculable rows (pre-filter): {n_total}")
    print(f"  dcNUxT rows (Ks > {KS_THRESHOLD}):        {n_dcnuxt}")
    print(f"  Of these, testable (non-NA p-value): {n_testable}")
    print(f"  Significant at raw p<0.05:            {n_sig_raw} ({100*n_sig_raw/n_testable:.1f}%)")
    print(f"  Significant at BH q<0.05:              {n_sig_bh} ({100*n_sig_bh/n_testable:.1f}%)")
    print(f"  Mean Ka/Ks (dcNUxT):                   {dcnuxt['Ka_Ks'].mean():.3f}")
    print()

    return dcnuxt, {
        "label": label,
        "n_total": n_total,
        "n_dcnuxt": n_dcnuxt,
        "n_testable": n_testable,
        "n_sig_raw": n_sig_raw,
        "n_sig_bh": n_sig_bh,
        "mean_kaks": dcnuxt["Ka_Ks"].mean(),
    }


def main():
    if len(sys.argv) != 3:
        print("Usage: python3 bh_correction_dcnuxt.py <nupt_combined.tsv> <numt_combined.tsv>")
        sys.exit(1)

    nupt_path, numt_path = sys.argv[1], sys.argv[2]

    dcnupt, nupt_summary = process_class(nupt_path, "dcNUPT")
    dcnumt, numt_summary = process_class(numt_path, "dcNUMT")

    dcnupt.to_csv("dcnupt_final_bh.tsv", sep="\t", index=False)
    dcnumt.to_csv("dcnumt_final_bh.tsv", sep="\t", index=False)

    summary = pd.DataFrame([nupt_summary, numt_summary])
    summary.to_csv("summary_evolutionary_pressure.tsv", sep="\t", index=False)

    print("Wrote: dcnupt_final_bh.tsv, dcnumt_final_bh.tsv, summary_evolutionary_pressure.tsv")


if __name__ == "__main__":
    main()
