#!/usr/bin/env python3
"""
02_dedup_kaks_rows.py -- count every coding insertion once.

Some organelle gene sets hold the same gene under two or three names (Allium
sativum plastid: ycf2 + gp005 + gp025; Lycium: ycf1 + orf489; Zea: IR copies
ZemaCp066/ZemaCp110; ATPB + atpB; cemA + ycf10; mito: ATP1 + atpA, nad2 + ND2,
sdh3 + sdh3p, orf127 + rps12 ...). One nuclear locus then appears in the Ka/Ks
table once per name, with the same fragment. These rows are collapsed here:

  same species + same scaffold + different gene name + coordinates overlapping
  by >= 80% of BOTH rows (IDs are scaffold_start-end_gene[_REV]); the BLAST
  hits against two copies of one gene differ slightly at the ends, so exact
  coordinates or identical fragments are not required. Genuinely different
  neighbouring genes (psbD/psbC, atpB/atpE, ndhC/ndhK) overlap by a few dozen bp
  at most and stay separate.
  -> one row per group, kept under the best name:
       real gene name  >  locus tag / ORF / "unknown"
       conventional    >  synonym (ycf10, pafII, sdh3p, nd2, atpA in mito ...)
       normal case     >  ALL CAPS
       then alphabetical

The same rule also collapses rows of the SAME gene that overlap >= 80% at one
nuclear locus (one insertion hit twice, e.g. against two copies of the gene or
by two overlapping BLAST hits); the longest row is kept. Both kinds are
reported separately.

Rows whose gene is an ORF, hypothetical protein or locus tag are kept (they
are only preferred less when choosing the name of a group).

Works on species/<SP>/results/<SP>_{nupt,numt}_kaks.csv (01_*_gene_alignment.sh)
in place, after copying each file to <file>.pre_dedup (mtime kept). Re-running
always starts from that backup, so the script can be run more than once.
Writes duplicate_rows_dropped.tsv (every dropped ID with the ID kept).

Usage (from the working directory):
  python3 02_dedup_kaks_rows.py combined_dataset.csv            # dry run: counts only
  python3 02_dedup_kaks_rows.py combined_dataset.csv --apply    # rewrite the files
"""
import collections, os, re, shutil, sys
import pandas as pd

APPLY = "--apply" in sys.argv
SUF = ".pre_dedup"
CSV = next((a for a in sys.argv[1:] if not a.startswith("--")), "combined_dataset.csv")
AXT = {"nupt": "nupt_input.axt", "numt": "numt_input.axt"}
TAG = re.compile(r"^(orf|hypothetical|gp\d|p\d{3}$|locus|unknown$|[A-Za-z]{2,6}(cp|mp)\d{2,}$)", re.I)
ALT = {"nupt": re.compile(r"^(ycf10|ycf9|pafi|pafii|ycf5|ycf6|ycf7|clpp1)$", re.I),
       "numt": re.compile(r"^(atpa|nd\d|nd4l|sdh\d+p|cox\d+p|nad\d+p)$", re.I)}

MIN_OV = 0.8
LOC = re.compile(r"^(.+?)_(\d+)-(\d+)(?:_[^-]*)?$")   # scaffold_start-end[_LOCUSPREFIX]

def gene_of(q):
    return q.replace(".mafft", "").replace("_REV", "").split("_")[-1]

def locus_of(q):
    return "_".join(q.replace(".mafft", "").replace("_REV", "").split("_")[:-1])

def span(q):
    m = LOC.match(locus_of(q))
    return abs(int(m.group(3)) - int(m.group(2))) + 1 if m else 0

def rank(q, cls):
    name = gene_of(q)
    return (bool(TAG.match(name)), bool(ALT[cls].match(name)),
            name.isupper() and any(c.isalpha() for c in name), -span(q), name)

def axt(p):
    d = {}
    if not os.path.isfile(p):
        return d
    L = [l.rstrip("\n") for l in open(p)]
    i = 0
    while i < len(L):
        if L[i].strip() and i + 2 < len(L):
            d[L[i].split()[0].replace(".mafft", "")] = L[i + 2].replace("-", "").upper(); i += 3
        else:
            i += 1
    return d

species = pd.read_csv(CSV)["SPECIES"].dropna().astype(str).str.strip().unique()
log = []
EX = collections.defaultdict(list)
for cls in ("nupt", "numt"):
    tot = collections.Counter()
    bysp = collections.Counter()
    leftover = []
    n_found = 0
    for sp in species:
        path = f"species/{sp}/results/{sp}_{cls}_kaks.csv"
        src = path + SUF if os.path.isfile(path + SUF) else path
        if not os.path.isfile(src):
            continue
        n_found += 1
        t = pd.read_csv(src, sep="\t", dtype=str, keep_default_na=False)
        if t.empty:
            continue
        a = axt(f"species/{sp}/results/{AXT[cls]}")
        ids = t.iloc[:, 0].tolist()
        # union-find over rows: same scaffold, different gene name, reciprocal overlap >= MIN_OV
        parent = list(range(len(ids)))
        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]; i = parent[i]
            return i
        byscaf = collections.defaultdict(list)
        exact = collections.defaultdict(list)
        for i, q in enumerate(ids):
            exact[locus_of(q)].append(i)
            m = LOC.match(locus_of(q))
            if not m:
                tot["unparsed"] += 1
                if len(EX["unparsed", cls]) < 4: EX["unparsed", cls].append(f"{sp}: {q}")
            if m:
                s0, e0 = sorted((int(m.group(2)), int(m.group(3))))
                byscaf[m.group(1)].append((s0, e0, i))
        for iv in byscaf.values():
            iv.sort()
            for j, (s1, e1, i1) in enumerate(iv):
                for s2, e2, i2 in iv[j + 1:]:
                    if s2 > e1:
                        break
                    ov = min(e1, e2) - s2 + 1
                    if ov >= MIN_OV * (e1 - s1 + 1) and ov >= MIN_OV * (e2 - s2 + 1):
                        parent[find(i1)] = find(i2)
        for idx in exact.values():          # identical locus string -> same insertion (also if unparsed)
            for i in idx[1:]:
                if gene_of(ids[i]).lower() != gene_of(ids[idx[0]]).lower():
                    parent[find(i)] = find(idx[0])
        groups = collections.defaultdict(list)
        for i in range(len(ids)):
            groups[find(i)].append(i)
        drop = set()
        ks_all = pd.to_numeric(t["Ks"], errors="coerce").tolist()
        for key, idx in groups.items():
            if len(idx) < 2:
                continue
            keep = min(idx, key=lambda i: rank(ids[i], cls))
            for i in idx:
                if i != keep:
                    drop.add(i)
                    same = gene_of(ids[i]).lower() == gene_of(ids[keep]).lower()
                    tot["drop_same_gene" if same else "drop_other_name"] += 1
                    if same and ks_all[i] > 0.1: tot["drop_same_gene_dc"] += 1
                    if not same and ks_all[i] > 0.1: tot["drop_other_name_dc"] += 1
                    log.append((cls, sp, ids[i], ids[keep], gene_of(ids[i]), gene_of(ids[keep])))
        ks = pd.to_numeric(t["Ks"], errors="coerce")
        kept = [i for i in range(len(t)) if i not in drop]
        tot["rows"] += len(t); tot["rows_after"] += len(kept)
        tot["coding"] += int(ks.notna().sum()); tot["coding_after"] += int(ks.iloc[kept].notna().sum())
        tot["dc"] += int((ks > 0.1).sum()); tot["dc_after"] += int((ks.iloc[kept] > 0.1).sum())
        if drop:
            bysp[sp] = len(drop)
        # identical fragments under different genes at DIFFERENT loci (remain after this step)
        fr = collections.defaultdict(set)
        for i in kept:
            if ks.iloc[i] > 0.1 and a.get(ids[i].replace(".mafft", "")):
                fr[a[ids[i].replace(".mafft", "")]].add(gene_of(ids[i]).lower())
        leftover += [(sp, sorted(g)) for g in fr.values() if len(g) > 1]
        if APPLY and (drop or os.path.isfile(path + SUF)):   # also restores files filtered in an earlier run
            if not os.path.isfile(path + SUF):
                shutil.copy2(path, path + SUF)
            st = os.stat(src)
            t.iloc[kept].to_csv(path, sep="\t", index=False)
            os.utime(path, (st.st_atime, st.st_mtime))
    print(f"\n== {cls.upper()} {'(applied)' if APPLY else '(dry run)'}")
    print(f"  Ka/Ks tables found   {n_found} of {len(species)} species")
    if not n_found:
        print(f"  WARNING: no species/<SP>/results/<SP>_{cls}_kaks.csv; run 01_{cls}_gene_alignment.sh first")
    print(f"  all rows        {tot['rows']:>6} -> {tot['rows_after']}")
    print(f"  coding (Ks ok)  {tot['coding']:>6} -> {tot['coding_after']}")
    print(f"  dcNUxT (Ks>0.1) {tot['dc']:>6} -> {tot['dc_after']}")
    print(f"  species changed: {len(bysp)}  {dict(bysp.most_common(8))}")
    print(f"  IDs whose coordinates could not be parsed: {tot['unparsed']}  e.g. {EX['unparsed', cls]}")
    print(f"  dropped as same gene under another name: {tot['drop_other_name']} rows ({tot['drop_other_name_dc']} dc)")
    print(f"  dropped as same gene hit twice at one locus: {tot['drop_same_gene']} rows ({tot['drop_same_gene_dc']} dc)")
    print(f"  identical dc fragments left under different genes (different loci): {len(leftover)} {leftover[:5]}")

pd.DataFrame(log, columns=["class", "species", "dropped_id", "kept_id", "dropped_gene", "kept_gene"]) \
  .to_csv("duplicate_rows_dropped.tsv", sep="\t", index=False)
kept_pairs = collections.Counter((c, d, k) for c, _, _, _, d, k in log)
print("\nmost common (class, dropped -> kept):", kept_pairs.most_common(12))
print("log -> duplicate_rows_dropped.tsv")
