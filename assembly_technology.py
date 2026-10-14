#!/usr/bin/env python3
"""
assembly_technology.py

Computes assembly contiguity statistics for every species from the .fai
indices, and tests whether organellar contamination tracks assembly
quality.

Contiguity is computed locally and needs no network access. Declared
sequencing technology can optionally be fetched from NCBI with
--fetch-tech, which requires internet and the `datasets` CLI or curl.

Contiguity is the more directly relevant measure: what affects the
analysis is whether organellar contigs survive in the assembly, and that
depends on assembly practice, of which read length is only one component.

Usage:
    python3 assembly_technology.py --fetch-tech     (after raw_scaffold_composition.py)
"""

import os
import sys
import glob
import subprocess

BASE = os.environ.get("BASE_DIR", os.getcwd())
GENOME_DIRS = [
    f"{BASE}/downloads/genome_fasta",
    f"{BASE}/multi_organelles/downloads/genome_fasta",
]
DATASETS = [f"{BASE}/combined_dataset.csv"]
COMPOSITION = "raw_scaffold_composition.tsv"


def load_assembly_map():
    sp2asm = {}
    for path in DATASETS:
        if not os.path.exists(path):
            continue
        with open(path) as f:
            for line in f:
                p = line.strip().split(",")
                if len(p) < 2 or p[0] in ("SPECIES", "species"):
                    continue
                sp2asm[p[0]] = p[1]
    return sp2asm


def nx(lengths, total, frac):
    """Nx and Lx: shortest scaffold at which cumulative length passes frac."""
    cum = 0
    for i, L in enumerate(lengths, 1):
        cum += L
        if cum >= total * frac:
            return L, i
    return 0, 0


def stats_from_fai(path):
    lengths = []
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                lengths.append(int(p[1]))
    if not lengths:
        return None
    lengths.sort(reverse=True)
    total = sum(lengths)
    n50, l50 = nx(lengths, total, 0.5)
    n90, _ = nx(lengths, total, 0.9)
    return {
        "n_scaffolds": len(lengths),
        "total_len": total,
        "max_scaffold": lengths[0],
        "n50": n50,
        "l50": l50,
        "n90": n90,
        "pct_in_scaffolds_over_1mb": 100 * sum(x for x in lengths if x >= 1_000_000) / total,
    }



# Submitters enter sequencing technology as free text, so classification
# is by keyword. Anything naming both a long- and a short-read platform is
# called hybrid.
LONG_READ = ["pacbio", "pac bio", "hifi", "sequel", "revio", "nanopore",
             "oxford nanopore", "ont", "promethion", "minion", "gridion",
             "clr", "rsii", "rs ii"]
SHORT_READ = ["illumina", "hiseq", "novaseq", "miseq", "nextseq", "genome analyzer",
              "bgiseq", "dnbseq", "mgiseq", "mgi", "solid", "454", "ion torrent"]


def classify_tech(s):
    if not s or s.upper() in ("NA", "NONE", ""):
        return "unknown"
    t = s.lower()
    has_long = any(k in t for k in LONG_READ)
    has_short = any(k in t for k in SHORT_READ)
    if has_long and has_short:
        return "hybrid"
    if has_long:
        return "long_read"
    if has_short:
        return "short_read"
    return "other"


def fetch_tech(accessions):
    """Declared sequencing technology from the NCBI Datasets API."""
    out = {}
    try:
        import json
        import urllib.request
    except ImportError:
        return out
    for i in range(0, len(accessions), 20):
        batch = accessions[i:i + 20]
        url = ("https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/"
               + ",".join(batch) + "/dataset_report")
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                data = json.load(r)
            for rep in data.get("reports", []):
                acc = rep.get("accession")
                info = rep.get("assembly_info", {})
                out[acc] = info.get("sequencing_tech", "NA")
        except Exception as e:
            print(f"  [warn] batch {i//20}: {e}", file=sys.stderr)
        print(f"  fetched {min(i+20, len(accessions))}/{len(accessions)}",
              file=sys.stderr)
    return out


def main():
    do_tech = "--fetch-tech" in sys.argv

    sp2asm = load_assembly_map()
    asm2sp = {v: k for k, v in sp2asm.items()}

    rows = {}
    for gdir in GENOME_DIRS:
        for f in glob.glob(f"{gdir}/*.fai"):
            parts = os.path.basename(f).split("_")
            if len(parts) < 2:
                continue
            asm = parts[0] + "_" + parts[1]
            sp = asm2sp.get(asm)
            if sp is None:
                continue
            s = stats_from_fai(f)
            if s:
                s["species"] = sp
                s["assembly"] = asm
                rows[sp] = s

    print(f"assemblies with stats: {len(rows)}")

    # join contamination if available
    comp = {}
    if os.path.exists(COMPOSITION):
        with open(COMPOSITION) as f:
            hdr = f.readline().rstrip("\n").split("\t")
            for line in f:
                p = line.rstrip("\n").split("\t")
                d = dict(zip(hdr, p))
                comp[d["species"]] = d
        print(f"contamination data for {len(comp)} species")
    else:
        print(f"[note] {COMPOSITION} not found; skipping contamination join")

    tech = {}
    if do_tech:
        print("\nfetching sequencing technology from NCBI...")
        tech = fetch_tech(sorted({r["assembly"] for r in rows.values()}))
        print(f"  got technology for {len(tech)} assemblies")

    out = "assembly_technology.tsv"
    cols = ["species", "assembly", "n_scaffolds", "total_len", "max_scaffold",
            "n50", "l50", "n90", "pct_in_scaffolds_over_1mb"]
    extra = []
    if comp:
        extra += ["numt_pct_suspect_bp", "nupt_pct_suspect_bp"]
    if tech:
        extra += ["sequencing_tech", "tech_class"]

    with open(out, "w") as fh:
        fh.write("\t".join(cols + extra) + "\n")
        for sp in sorted(rows):
            r = rows[sp]
            vals = [str(r[c]) if not isinstance(r[c], float) else f"{r[c]:.2f}"
                    for c in cols]
            if comp:
                d = comp.get(sp, {})
                vals += [d.get("numt_pct_suspect_bp", "NA"),
                         d.get("nupt_pct_suspect_bp", "NA")]
            if tech:
                raw = tech.get(r["assembly"], "NA")
                vals += [raw.replace("\t", " "), classify_tech(raw)]
            fh.write("\t".join(vals) + "\n")

    # contiguity distribution
    n50s = sorted(r["n50"] for r in rows.values())
    n = len(n50s)
    print(f"\nScaffold N50 across {n} assemblies:")
    for label, v in [("min", n50s[0]), ("25th", n50s[n//4]),
                     ("median", n50s[n//2]), ("75th", n50s[3*n//4]),
                     ("max", n50s[-1])]:
        print(f"  {label:<8} {v:>15,}")
    for thr, name in [(10_000, "<10 kb"), (100_000, "<100 kb"),
                      (1_000_000, "<1 Mb"), (10_000_000, "<10 Mb")]:
        k = sum(1 for x in n50s if x < thr)
        print(f"  N50 {name:<8} {k:>4} assemblies ({100*k/n:.1f}%)")

    # does contamination track contiguity?
    if comp:
        try:
            from scipy.stats import spearmanr
        except ImportError:
            print("\n[note] scipy unavailable; skipping correlation")
            return
        for key, label in [("nupt_pct_suspect_bp", "NUPT"),
                           ("numt_pct_suspect_bp", "NUMT")]:
            xs, ys = [], []
            for sp, r in rows.items():
                d = comp.get(sp)
                if not d or d.get(key) in (None, "NA"):
                    continue
                xs.append(r["n50"])
                ys.append(float(d[key]))
            if len(xs) > 10:
                rho, p = spearmanr(xs, ys)
                print(f"\nScaffold N50 vs {label} contamination: "
                      f"Spearman r={rho:+.3f}, p={p:.3g}, n={len(xs)}")

    if tech:
        from collections import Counter
        cls = Counter(classify_tech(tech.get(r["assembly"], "NA"))
                      for r in rows.values())
        n_tot = sum(cls.values())
        print(f"\nSequencing technology across {n_tot} assemblies:")
        for k in ["long_read", "hybrid", "short_read", "other", "unknown"]:
            if cls.get(k):
                print(f"  {k:<12} {cls[k]:>4} ({100*cls[k]/n_tot:5.1f}%)")

        if comp:
            print("\nMedian NUPT contamination by technology:")
            groups = {}
            for sp, r in rows.items():
                d = comp.get(sp)
                if not d or d.get("nupt_pct_suspect_bp") in (None, "NA"):
                    continue
                k = classify_tech(tech.get(r["assembly"], "NA"))
                groups.setdefault(k, []).append(float(d["nupt_pct_suspect_bp"]))
            for k in ["long_read", "hybrid", "short_read", "other", "unknown"]:
                v = sorted(groups.get(k, []))
                if v:
                    print(f"  {k:<12} n={len(v):>4}  median={v[len(v)//2]:6.2f}%")

    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
