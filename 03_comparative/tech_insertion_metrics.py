#!/usr/bin/env python3
"""
Per-species insertion-level metrics for the sequencing-technology test.

Prediction being tested: if short-read assemblers collapse or drop recent,
near-identical organellar insertions, short-read assemblies should be depleted
specifically in (a) high-identity hits and (b) long merged insertions, not just
in total content.

Reads, for each species:
  species/$SP/numts/numts.clean.out   (sseqid pident length sstart send)
  species/$SP/nupts/nupts.clean.out
  species/$SP/numts/numts.clean.out_AllFrags_merged10000.bed
  species/$SP/nupts/nupts.clean.out_AllFrags_merged10000.bed

Writes insertion_metrics.tsv (one row per species).

Usage:
  python3 tech_insertion_metrics.py [BASE_DIR] [OUT]
"""
import os
import sys

BASE = sys.argv[1] if len(sys.argv) > 1 else os.environ.get("BASE_DIR", os.getcwd())
OUT = sys.argv[2] if len(sys.argv) > 2 else "insertion_metrics.tsv"
MIN_LEN = 200            # same as numt_nupt_detection.sh
MERGE = 10000
HI_ID = (99.0, 95.0)     # identity thresholds for "recent" hits
LONG = (5000, 10000)     # merged-insertion length thresholds (bp)


def hit_metrics(path):
    n = 0
    hi = {t: 0 for t in HI_ID}
    hi_bp = {t: 0 for t in HI_ID}
    bp = 0
    if not os.path.isfile(path):
        return None
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                pid = float(f[1]); ln = int(f[2])
            except ValueError:
                continue
            if ln < MIN_LEN:
                continue
            n += 1
            bp += ln
            for t in HI_ID:
                if pid >= t:
                    hi[t] += 1
                    hi_bp[t] += ln
    if n == 0:
        return None
    out = {"n_hits": n}
    for t in HI_ID:
        k = f"{int(t)}"
        out[f"frac_hits_ge{k}"] = hi[t] / n
        out[f"frac_bp_ge{k}"] = hi_bp[t] / bp
    return out


def bed_metrics(path):
    if not os.path.isfile(path):
        return None
    lens = []
    with open(path) as fh:
        for line in fh:
            f = line.split("\t")
            if len(f) < 3:
                continue
            try:
                lens.append(int(f[2]) - int(f[1]))
            except ValueError:
                continue
    if not lens:
        return None
    lens.sort()
    tot = sum(lens)
    out = {"n_merged": len(lens), "median_len": lens[len(lens) // 2],
           "max_len": lens[-1]}
    for t in LONG:
        big = [x for x in lens if x >= t]
        out[f"n_ge{t // 1000}kb"] = len(big)
        out[f"frac_bp_ge{t // 1000}kb"] = sum(big) / tot
    return out


def main():
    sp_dir = os.path.join(BASE, "species")
    species = sorted(d for d in os.listdir(sp_dir)
                     if os.path.isdir(os.path.join(sp_dir, d)))
    cols = None
    rows = []
    missing = 0
    for sp in species:
        rec = {"species": sp}
        for cls in ("numt", "nupt"):
            d = os.path.join(sp_dir, sp, f"{cls}s")
            h = hit_metrics(os.path.join(d, f"{cls}s.clean.out"))
            b = bed_metrics(os.path.join(
                d, f"{cls}s.clean.out_AllFrags_merged{MERGE}.bed"))
            if h is None and b is None:
                missing += 1
            for src in (h, b):
                if src:
                    for k, v in src.items():
                        rec[f"{cls}_{k}"] = v
        rows.append(rec)
        if cols is None or len(rec) > len(cols):
            cols = list(rec)
    # union of keys, stable order
    keys = ["species"]
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(OUT, "w") as fh:
        fh.write("\t".join(keys) + "\n")
        for r in rows:
            vals = []
            for k in keys:
                v = r.get(k, "NA")
                vals.append(f"{v:.6g}" if isinstance(v, float) else str(v))
            fh.write("\t".join(vals) + "\n")
    print(f"Wrote {OUT}: {len(rows)} species "
          f"({missing} species x class combinations with no files)")


if __name__ == "__main__":
    main()
