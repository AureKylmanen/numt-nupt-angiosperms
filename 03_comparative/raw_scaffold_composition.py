#!/usr/bin/env python3
"""
raw_scaffold_composition.py

Applies the scaffold-composition check to ALL merged NUMT and NUPT
intervals, not just the diverged coding subset.

This matters for the per-species proportions used in the phylogenetic and
genome-size analyses. Organellar contigs left in an assembly are detected
as wall-to-wall "insertions", so any species whose assembly retains them
has its proportion inflated. Because retention of such contigs depends on
assembly practice rather than biology, the inflation is uneven across
species and makes assembly contiguity a potential confound.

Each scaffold is classified once, then intervals are attributed to their
scaffold's class:

  ORGANELLAR_SELF   >80% covered by intervals of the same organelle
  CROSS_ORGANELLE   >80% covered by the other organelle
  ORGANELLE_RICH    combined organellar coverage above --max-organellar
  SHORT_SCAFFOLD    shorter than --min-scaffold
  NUCLEAR           none of the above

Short scaffolds are reported but NOT counted as suspect. Many assemblies
are scaffold-level with N50 well below any useful threshold, so flagging
every short scaffold measures assembly fragmentation rather than
organellar contamination. Composition, not length, is the evidence.

Usage:
    python3 raw_scaffold_composition.py [--min-scaffold 100000]
                                        [--max-organellar 50]
"""

import os
import sys
import glob
from collections import defaultdict, Counter

BASE = os.environ.get("BASE_DIR", os.getcwd())
GENOME_DIRS = [
    f"{BASE}/downloads/genome_fasta",
    f"{BASE}/multi_organelles/downloads/genome_fasta",
]
DATASETS = [f"{BASE}/combined_dataset.csv"]


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


def load_fai(sp2asm):
    asm2sp = {v: k for k, v in sp2asm.items()}
    fai = {}
    for gdir in GENOME_DIRS:
        for f in glob.glob(f"{gdir}/*.fai"):
            parts = os.path.basename(f).split("_")
            if len(parts) < 2:
                continue
            sp = asm2sp.get(parts[0] + "_" + parts[1])
            if sp is None:
                continue
            d = {}
            with open(f) as fh:
                for line in fh:
                    p = line.rstrip("\n").split("\t")
                    if len(p) >= 2:
                        d[p[0]] = int(p[1])
            fai[sp] = d
    return fai



def union_bp(intervals):
    """Total bp covered by a set of intervals, counting overlaps once.

    The unmerged BED files contain one line per BLAST HSP, and HSPs
    routinely overlap (multiple local alignments over the same region).
    Summing interval lengths therefore double-counts and can report
    coverage above 100% of the scaffold.
    """
    if not intervals:
        return 0
    intervals = sorted(intervals)
    total = 0
    cur_s, cur_e = intervals[0]
    for s, e in intervals[1:]:
        if s > cur_e:
            total += cur_e - cur_s
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    total += cur_e - cur_s
    return total


def read_bed(path):
    """Returns (count per scaffold, union bp per scaffold)."""
    n = defaultdict(int)
    iv = defaultdict(list)
    if not os.path.exists(path):
        return n, {}
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            try:
                iv[p[0]].append((int(p[1]), int(p[2])))
                n[p[0]] += 1
            except ValueError:
                continue
    bp = {sc: union_bp(v) for sc, v in iv.items()}
    return n, bp


def classify(length, pct_self, pct_other, min_scaffold, max_organellar):
    if length is None:
        return "NO_FAI"
    if pct_self > 80:
        return "ORGANELLAR_SELF"
    if pct_other > 80:
        return "CROSS_ORGANELLE"
    if length < min_scaffold:
        return "SHORT_SCAFFOLD"
    if pct_self + pct_other > max_organellar:
        return "ORGANELLE_RICH"
    return "NUCLEAR"


def main():
    argv = sys.argv[1:]

    def opt(flag, default, cast=int):
        if flag in argv:
            i = argv.index(flag)
            v = cast(argv[i + 1])
            del argv[i:i + 2]
            return v
        return default

    min_scaffold = opt("--min-scaffold", 100000)
    max_organellar = opt("--max-organellar", 50.0, float)

    sp2asm = load_assembly_map()
    fai = load_fai(sp2asm)
    print(f"species with FAI: {len(fai)}")
    print(f"thresholds: scaffold >= {min_scaffold:,} bp, "
          f"combined organellar <= {max_organellar}%\n")

    per_species = []
    totals = {"numt": Counter(), "nupt": Counter()}
    totals_bp = {"numt": Counter(), "nupt": Counter()}

    species = sorted(os.path.basename(p)
                     for p in glob.glob(f"{BASE}/species/*")
                     if os.path.isdir(p))

    for sp in species:
        d = f"{BASE}/species/{sp}"
        # Composition must come from the UNMERGED intervals. The
        # merged10000 files join intervals up to 10 kb apart, so a merged
        # interval spans the intervening nuclear sequence; reading coverage
        # from them credits that nuclear DNA as organellar and grossly
        # overstates scaffold composition.
        n_mt, bp_mt = read_bed(f"{d}/numts/numts.clean.out_AllFrags.bed")
        n_pt, bp_pt = read_bed(f"{d}/nupts/nupts.clean.out_AllFrags.bed")
        if not n_mt and not n_pt:
            continue
        lengths = fai.get(sp, {})

        row = {"species": sp}
        for label, n_own, bp_own, bp_oth in (
                ("numt", n_mt, bp_mt, bp_pt),
                ("nupt", n_pt, bp_pt, bp_mt)):
            cls_n, cls_bp = Counter(), Counter()
            for sc, cnt in n_own.items():
                L = lengths.get(sc)
                if L:
                    ps = 100 * bp_own.get(sc, 0) / L
                    po = 100 * bp_oth.get(sc, 0) / L
                else:
                    ps = po = 0.0
                v = classify(L, ps, po, min_scaffold, max_organellar)
                cls_n[v] += cnt
                cls_bp[v] += bp_own.get(sc, 0)
            totals[label].update(cls_n)
            totals_bp[label].update(cls_bp)
            tot_n = sum(cls_n.values())
            tot_bp = sum(cls_bp.values())
            susp_n = sum(cls_n[k] for k in
                         ("ORGANELLAR_SELF", "CROSS_ORGANELLE"))
            susp_bp = sum(cls_bp[k] for k in
                          ("ORGANELLAR_SELF", "CROSS_ORGANELLE"))
            row[f"{label}_n"] = tot_n
            row[f"{label}_bp"] = tot_bp
            row[f"{label}_suspect_n"] = susp_n
            row[f"{label}_suspect_bp"] = susp_bp
            row[f"{label}_pct_suspect_bp"] = (100 * susp_bp / tot_bp) if tot_bp else 0.0
        per_species.append(row)

    out = "raw_scaffold_composition.tsv"
    cols = ["species",
            "numt_n", "numt_bp", "numt_suspect_n", "numt_suspect_bp", "numt_pct_suspect_bp",
            "nupt_n", "nupt_bp", "nupt_suspect_n", "nupt_suspect_bp", "nupt_pct_suspect_bp"]
    with open(out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in per_species:
            fh.write("\t".join(
                f"{r[c]:.2f}" if isinstance(r[c], float) else str(r[c])
                for c in cols) + "\n")

    order = ["NUCLEAR", "ORGANELLE_RICH", "SHORT_SCAFFOLD",
             "CROSS_ORGANELLE", "ORGANELLAR_SELF", "NO_FAI"]
    for label, name in (("numt", "NUMTs"), ("nupt", "NUPTs")):
        tn = sum(totals[label].values())
        tb = sum(totals_bp[label].values())
        print(f"=== all {name} ===")
        print(f"  total: {tn:,} intervals, {tb:,} bp")
        for v in order:
            if totals[label].get(v):
                print(f"    {v:<18} {totals[label][v]:>10,} intervals  "
                      f"{totals_bp[label][v]:>14,} bp  "
                      f"({100*totals_bp[label][v]/tb:5.2f}% of bp)")
        susp_bp = sum(totals_bp[label][k] for k in
                      ("ORGANELLAR_SELF", "CROSS_ORGANELLE"))
        print(f"  suspect bp (organellar composition only): "
              f"{susp_bp:,} ({100*susp_bp/tb:.2f}%)")
        short_bp = totals_bp[label].get("SHORT_SCAFFOLD", 0)
        print(f"  on short scaffolds (reported, not counted as suspect): "
              f"{short_bp:,} ({100*short_bp/tb:.2f}%)\n")

    print("Species where suspect scaffolds carry the largest share of NUMT bp:")
    for r in sorted(per_species, key=lambda x: -x["numt_pct_suspect_bp"])[:15]:
        print(f"  {r['species'][:38]:<38} {r['numt_pct_suspect_bp']:6.2f}%  "
              f"({r['numt_suspect_bp']:,} of {r['numt_bp']:,} bp)")

    print("\nSpecies where suspect scaffolds carry the largest share of NUPT bp:")
    for r in sorted(per_species, key=lambda x: -x["nupt_pct_suspect_bp"])[:15]:
        print(f"  {r['species'][:38]:<38} {r['nupt_pct_suspect_bp']:6.2f}%  "
              f"({r['nupt_suspect_bp']:,} of {r['nupt_bp']:,} bp)")

    print(f"\nWrote {out}")


if __name__ == "__main__":
    main()
