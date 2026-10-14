#!/usr/bin/env python3
"""
06_mannwhitney_kaks_comparison.py

Compares the Ka/Ks distributions of dcNUPTs vs dcNUMTs using a
Mann-Whitney U test.

Usage:
    python3 06_mannwhitney_kaks_comparison.py dcnupt_final_bh.tsv dcnumt_final_bh.tsv
"""

import sys
import pandas as pd
from scipy.stats import mannwhitneyu, median_test
import numpy as np


def main():
    if len(sys.argv) != 3:
        print("Usage: python3 mannwhitney_kaks_comparison.py <dcnupt.tsv> <dcnumt.tsv>")
        sys.exit(1)

    nupt = pd.read_csv(sys.argv[1], sep="\t")
    numt = pd.read_csv(sys.argv[2], sep="\t")

    nupt_kaks = nupt["Ka_Ks"].dropna().values
    numt_kaks = numt["Ka_Ks"].dropna().values

    print(f"dcNUPT: n={len(nupt_kaks)}, median={np.median(nupt_kaks):.4f}, mean={np.mean(nupt_kaks):.4f}")
    print(f"dcNUMT: n={len(numt_kaks)}, median={np.median(numt_kaks):.4f}, mean={np.mean(numt_kaks):.4f}")
    print()

    # Two-sided Mann-Whitney U test
    stat, p_two_sided = mannwhitneyu(nupt_kaks, numt_kaks, alternative="two-sided")
    print(f"Mann-Whitney U (two-sided): U={stat:.1f}, p={p_two_sided:.4g}")

    # One-sided tests to establish direction, if two-sided is significant
    stat_less, p_less = mannwhitneyu(nupt_kaks, numt_kaks, alternative="less")
    stat_greater, p_greater = mannwhitneyu(nupt_kaks, numt_kaks, alternative="greater")
    print(f"  H1: NUPT Ka/Ks < NUMT Ka/Ks : p={p_less:.4g}")
    print(f"  H1: NUPT Ka/Ks > NUMT Ka/Ks : p={p_greater:.4g}")

    # Effect size: rank-biserial correlation
    n1, n2 = len(nupt_kaks), len(numt_kaks)
    r = 1 - (2 * stat) / (n1 * n2)
    print(f"\nRank-biserial correlation (effect size): r={r:.4f}")
    print("(|r| < 0.1 negligible, 0.1-0.3 small, 0.3-0.5 medium, >0.5 large)")

    alpha = 0.05
    print()
    if p_two_sided < alpha:
        direction = "higher" if np.median(nupt_kaks) > np.median(numt_kaks) else "lower"
        print(f"RESULT: Significant difference (p={p_two_sided:.4g} < {alpha}).")
        print(f"        dcNUPT Ka/Ks is significantly {direction} than dcNUMT Ka/Ks.")
    else:
        print(f"RESULT: No significant difference detected (p={p_two_sided:.4g} >= {alpha}).")
        print("        Cannot support a directional claim about NUPT vs NUMT selective constraint.")


if __name__ == "__main__":
    main()
