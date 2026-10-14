#!/usr/bin/env python3
"""
make_transcription_table.py

Build the transcription-evidence table and summary counts from
all_junction_counts.tsv, so the manuscript numbers come from the data rather
than being transcribed by hand.

Classification
--------------
  TRANSCRIBED     >= --min-junction reads at BOTH boundaries
  NOT_TRANSCRIBED fewer than that, but flanking depth >= --min-depth
                  (the region is sequenced, so the absence is informative)
  UNTESTED        flanking depth below --min-depth on both sides
                  (no reads in the region; absence carries no information)

Excludes, with a printed reason:
  * loci whose flanking sequence is itself organellar (--exclude-flank file,
    one "Species<TAB>host_gene" per line)
  * duplicate detections of the same locus: identical species and overlapping
    interval, keeping the longer fragment

Usage
-----
  python3 make_transcription_table.py \
      --counts all_junction_counts.tsv \
      --exclude-flank organellar_flanks.txt \
      --min-junction 5 --min-depth 20 \
      --out transcription
"""

import argparse
import csv
import sys
from collections import defaultdict


def num(x, default=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", default="all_junction_counts.tsv")
    ap.add_argument("--exclude-flank", default=None,
                    help="species<TAB>host_gene per line; loci whose junction-"
                         "defining flank is itself organellar")
    ap.add_argument("--min-junction", type=int, default=5)
    ap.add_argument("--min-depth", type=float, default=20.0)
    ap.add_argument("--out", default="transcription")
    args = ap.parse_args()

    excl = set()
    if args.exclude_flank:
        try:
            with open(args.exclude_flank) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) >= 2:
                        excl.add((f[0].strip(), f[1].strip()))
        except OSError as e:
            sys.exit(f"cannot read {args.exclude_flank}: {e}")

    rows = []
    with open(args.counts, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append(r)
    if not rows:
        sys.exit("no rows read - check --counts path and that it has a header")

    # collapse duplicate detections: same species, overlapping interval.
    # A composite insertion is found by both the plastid and mitochondrial
    # search and appears twice; keep the longer fragment.
    by_species = defaultdict(list)
    for r in rows:
        by_species[r["species"]].append(r)

    kept, dupes = [], []
    for sp, rs in by_species.items():
        rs.sort(key=lambda x: (x["chrom"], int(x["start"])))
        for r in rs:
            s, e = int(r["start"]), int(r["end"])
            hit = None
            for k in kept:
                if k["species"] == sp and k["chrom"] == r["chrom"] \
                   and int(k["end"]) > s and int(k["start"]) < e:
                    hit = k
                    break
            if hit is None:
                kept.append(r)
            else:
                if (e - s) > (int(hit["end"]) - int(hit["start"])):
                    kept[kept.index(hit)] = r
                    dupes.append(hit)
                else:
                    dupes.append(r)

    out, excluded = [], []
    tally = defaultdict(int)
    for r in kept:
        key = (r["species"], r["host_gene"])
        if key in excl:
            excluded.append(r)
            continue
        j5, j3 = int(num(r["j5"])), int(num(r["j3"]))
        f5, f3 = num(r["flank5_depth"]), num(r["flank3_depth"])
        if j5 >= args.min_junction and j3 >= args.min_junction:
            cls = "TRANSCRIBED"
        elif max(f5, f3) >= args.min_depth:
            cls = "NOT_TRANSCRIBED"
        else:
            cls = "UNTESTED"
        tally[cls] += 1
        r["_class"] = cls
        out.append(r)

    # ---- summary ----
    print(f"rows read:              {len(rows)}")
    print(f"duplicate detections:   {len(dupes)}"
          + (f"  ({', '.join(sorted(set(d['species'] for d in dupes)))})" if dupes else ""))
    print(f"excluded, organellar flank: {len(excluded)}"
          + (f"  ({', '.join(sorted(set(e['species'] for e in excluded)))})" if excluded else ""))
    print(f"distinct loci assessed: {len(out)}")
    for k in ("TRANSCRIBED", "NOT_TRANSCRIBED", "UNTESTED"):
        print(f"  {k:16s} {tally[k]}")
    print(f"\nthresholds: >= {args.min_junction} junction reads at both boundaries; "
          f"flank depth >= {args.min_depth}x for an informative negative")

    # ---- full TSV ----
    cols = ["species", "kind", "donor_gene", "host_gene", "verdict", "chrom",
            "start", "end", "j5", "j3", "flank5_depth", "flank3_depth",
            "inside_depth", "_class"]
    with open(args.out + "_all.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in sorted(out, key=lambda x: (x["_class"], x["species"])):
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

    # ---- LaTeX table of transcribed loci ----
    def sp_latex(s):
        p = s.replace("_", " ")
        return "\\emph{" + p + "}"

    tr = [r for r in out if r["_class"] == "TRANSCRIBED"]
    tr.sort(key=lambda x: -(int(num(x["j5"])) + int(num(x["j3"]))))
    with open(args.out + "_table.tex", "w") as fh:
        fh.write("\\begin{table}[t]\n\\centering\n\\small\n")
        fh.write("\\caption{dcNUxT loci with junction-spanning read support. "
                 "$j_5$, $j_3$: reads spanning the upstream and downstream "
                 "insertion boundary, requiring $\\ge$10\\,bp of aligned "
                 "sequence on each side and mapping quality $\\ge$20. Flank: "
                 "mean depth over 300\\,bp either side. Inside: mean depth over "
                 "the insertion, shown for context only, since organellar "
                 "transcript reads align there.}\n")
        fh.write("\\label{tab:transcription}\n")
        fh.write("\\begin{tabular}{llrrrrr}\n\\hline\n")
        fh.write("Species & Donor & $j_5$ & $j_3$ & Flank 5$'$ & Flank 3$'$ & Inside \\\\\n\\hline\n")
        for r in tr:
            fh.write(f"{sp_latex(r['species'])} & \\emph{{{r['donor_gene']}}} & "
                     f"{int(num(r['j5']))} & {int(num(r['j3']))} & "
                     f"{num(r['flank5_depth']):.1f} & {num(r['flank3_depth']):.1f} & "
                     f"{num(r['inside_depth']):.1f} \\\\\n")
        fh.write("\\hline\n\\end{tabular}\n\\end{table}\n")

    print(f"\nWrote {args.out}_all.tsv and {args.out}_table.tex "
          f"({len(tr)} transcribed loci)")


if __name__ == "__main__":
    sys.exit(main())
