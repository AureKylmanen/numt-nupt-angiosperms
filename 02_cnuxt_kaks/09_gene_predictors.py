#!/usr/bin/env python3
"""
09_gene_predictors.py

Tests what predicts how often an organellar gene is transferred to the
nucleus, and how often a transferred copy is retained under purifying
selection - and whether these relationships differ between plastid- and
mitochondrial-derived insertions.

Reads gene_summary_nupt.tsv and gene_summary_numt.tsv from
08_gene_level_summary.py.

Tests run per organelle:
  1. gene length vs cNUxT count            (does length drive transfer?)
  2. gene length vs species breadth        (widespread or concentrated?)
  3. gene length vs ret%                   (does length predict age?)
  4. gene length vs sel%                   (does length predict constraint?)
  5. length-normalised transfer vs sel%    (does transfer frequency predict
                                            retention once length is
                                            accounted for?)

Test 5 is the substantive one: gene length is expected to drive transfer
simply because a longer gene is a larger target, so the question is
whether genes that transfer more than their length predicts are also more
often retained.

Between organelles, correlations are compared with Fisher's r-to-z, which
tests whether two independent correlation coefficients differ.

Usage:
    python3 09_gene_predictors.py [--min-cnuxt 5] [--min-dcnuxt 3]
"""

import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, pearsonr, norm

NUPT_TSV = "gene_summary_nupt.tsv"
NUMT_TSV = "gene_summary_numt.tsv"

# Annotation-description artifacts and genes shared between organelles.
# Plastid gene names appearing in mitochondrial annotation may be genuine
# plastid-to-mitochondrion transfers rather than mitochondrial genes, and
# carry plastid reference lengths, so they are excluded from mitochondrial
# length-based tests.
DROP_ALWAYS = {"oxygenase large subunit", "unknown", "all_plastid", "all_mito"}
PLASTID_IN_MITO = {"ycf1", "ycf2", "ycf3", "ycf4", "psaa", "psab", "psba",
                   "psbb", "psbc", "psbd", "rbcl", "atpa", "atpb", "atpe",
                   "atph", "atpi", "peta", "petb", "petd", "petg", "petl",
                   "petn", "ndhb", "ndhe", "ndhh", "infa", "cema"}


def fisher_r_to_z(r1, n1, r2, n2):
    """Test whether two independent correlations differ."""
    if min(n1, n2) < 4:
        return np.nan, np.nan
    z1 = np.arctanh(np.clip(r1, -0.999999, 0.999999))
    z2 = np.arctanh(np.clip(r2, -0.999999, 0.999999))
    se = np.sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
    z = (z1 - z2) / se
    return z, 2 * (1 - norm.cdf(abs(z)))


def load(path, label, drop_plastid):
    t = pd.read_csv(path, sep="\t")
    gcol = "gene" if "gene" in t.columns else t.columns[0]
    orf = t[gcol].astype(str).str.lower().str.match(r"^(orf|hypothetical|gp\d|p\d{3}$|locus)")
    if orf.any():
        print(f"  [{label}] excluded {int(orf.sum())} ORF/hypothetical genes: "
              f"{', '.join(t.loc[orf, gcol].astype(str).head(8))}{' ...' if orf.sum() > 8 else ''}")
    t = t[~orf].copy()
    t = t[~t["gene"].isin(DROP_ALWAYS)]
    if drop_plastid:
        n_before = len(t)
        t = t[~t["gene"].isin(PLASTID_IN_MITO)]
        dropped = n_before - len(t)
        if dropped:
            print(f"  [{label}] excluded {dropped} plastid-named genes from "
                  f"mitochondrial set")
    if "length_bp" not in t.columns:
        raise SystemExit("gene_summary tables lack length_bp; run "
                         "07_build_gene_lengths.py before 08_gene_level_summary.py")
    t = t[t["length_bp"].notna()]
    return t


def report(name, r, p, n):
    stars = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
    print(f"  {name:<48} r={r:+.3f}  p={p:.3g} {stars}  n={n}")
    return r, n


def analyse(t, label, min_cnuxt, min_dcnuxt):
    print(f"\n=== {label} ===")
    res = {}

    a = t[t["cNUxTs"] >= min_cnuxt]
    r, p = spearmanr(a["length_bp"], a["cNUxTs"])
    res["len_vs_transfer"] = report("gene length vs cNUxT count", r, p, len(a))

    r, p = spearmanr(a["length_bp"], a["n_spp_c"])
    res["len_vs_breadth"] = report("gene length vs species breadth", r, p, len(a))

    r, p = spearmanr(a["length_bp"], a["ret_pct"])
    res["len_vs_ret"] = report("gene length vs ret% (dcNUxT/cNUxT)", r, p, len(a))

    b = t[(t["dcNUxTs"] >= min_dcnuxt) & t["sel_pct"].notna()]
    if len(b) >= 4:
        r, p = spearmanr(b["length_bp"], b["sel_pct"])
        res["len_vs_sel"] = report("gene length vs sel% (selected/dcNUxT)",
                                   r, p, len(b))

        # length-normalised transfer: residual of log(cNUxT) on log(length)
        x = np.log10(b["length_bp"].values)
        y = np.log10(b["cNUxTs"].values)
        slope, intercept = np.polyfit(x, y, 1)
        resid = y - (slope * x + intercept)
        r, p = spearmanr(resid, b["sel_pct"])
        res["normtransfer_vs_sel"] = report(
            "length-normalised transfer vs sel%", r, p, len(b))

        # permutation test, since n is small and sel% is bounded
        obs = abs(r)
        rng = np.random.default_rng(0)
        perm = sum(abs(spearmanr(resid, rng.permutation(b["sel_pct"].values))[0])
                   >= obs for _ in range(10000))
        print(f"  {'  permutation p (10,000 shuffles)':<48} "
              f"p={(perm + 1) / 10001:.4f}")
    else:
        print(f"  too few genes with >= {min_dcnuxt} dcNUxTs for "
              f"retention tests (n={len(b)})")

    return res


def main():
    argv = sys.argv[1:]

    def opt(flag, default):
        if flag in argv:
            i = argv.index(flag)
            return int(argv[i + 1])
        return default

    min_cnuxt = opt("--min-cnuxt", 5)
    min_dcnuxt = opt("--min-dcnuxt", 3)

    print(f"Thresholds: cNUxT >= {min_cnuxt} for transfer tests, "
          f"dcNUxT >= {min_dcnuxt} for retention tests")

    nupt = load(NUPT_TSV, "NUPT", drop_plastid=False)
    numt = load(NUMT_TSV, "NUMT", drop_plastid=True)

    rn = analyse(nupt, "NUPT (plastid-derived)", min_cnuxt, min_dcnuxt)
    rm = analyse(numt, "NUMT (mitochondrial-derived)", min_cnuxt, min_dcnuxt)

    print("\n=== NUPT vs NUMT: do the relationships differ? ===")
    print("  (Fisher r-to-z on independent correlations)")
    labels = {
        "len_vs_transfer": "gene length vs transfer",
        "len_vs_breadth": "gene length vs species breadth",
        "len_vs_ret": "gene length vs ret%",
        "len_vs_sel": "gene length vs sel%",
        "normtransfer_vs_sel": "length-normalised transfer vs sel%",
    }
    for k, name in labels.items():
        if k in rn and k in rm:
            r1, n1 = rn[k]
            r2, n2 = rm[k]
            z, p = fisher_r_to_z(r1, n1, r2, n2)
            if np.isnan(z):
                continue
            stars = ("***" if p < 0.001 else "**" if p < 0.01
                     else "*" if p < 0.05 else "ns")
            print(f"  {name:<40} NUPT r={r1:+.3f} (n={n1})  "
                  f"NUMT r={r2:+.3f} (n={n2})  z={z:+.2f} p={p:.3g} {stars}")

    print("\nNote: correlations are across genes, not insertions, so n is the "
          "number of genes.\nGenes below the thresholds are excluded because "
          "ret% and sel% are unstable\nwhen computed from a handful of "
          "sequences.")


if __name__ == "__main__":
    main()
