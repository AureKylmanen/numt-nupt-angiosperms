#!/usr/bin/env python3
"""
dcnuxt_annotation_overlap.py   (v2 - fixed)

Classify the genomic context of dcNUPTs and dcNUMTs against genome annotation.

Changes from v1
---------------
 * smart_open() sniffs magic bytes instead of trusting the .gz suffix
   (some downloads are truncated or are plain text named .gz).
 * --audit-only uses a cheap seqid-only reader, no full feature parse.
 * Unreadable/corrupt GFFs are recorded as 'gff_unreadable' in the audit
   instead of aborting the whole run.
 * GFF discovery also tries the accession with the version stripped, so
   GCA_xxx.1 on disk still matches GCA_xxx.2 in the CSV.
 * Audit summary reports dcNUxTs lost per status, not just species counts.

Pipeline
--------
 1. Read combined_dataset.csv (SPECIES,ASSEMBLY,PLASTID,MITO).
 2. Per species read the corrected KaKs tables and keep Ks > MIN_KS.
 3. Recover genome coordinates by rebuilding the exact filename
    sanitisation used by the alignment scripts from the AllFrags BED,
    then longest-prefix-matching the KaKs `Sequence` column.
 4. Deduplicate identical intervals (duplicate BLAST HSPs against
    ycf1/ycf2/ndhB inflate counts otherwise).
 5. Locate a GFF, verify seqid concordance, translate GCA<->GCF seqids
    via *_assembly_report.txt / *sequence_report.jsonl when present.
 6. Classify: CDS > UTR/exon > intronic > pseudogene > intergenic.
 7. Permutation null: shuffle each locus within its own scaffold.

Outputs
-------
  <prefix>_loci.tsv        one row per dcNUxT locus
  <prefix>_summary.tsv     per-class counts, observed vs expected
  <prefix>_gff_audit.tsv   per-species GFF status and seqid concordance

Usage
-----
  # fast audit first
  python3 dcnuxt_annotation_overlap.py \
      --base /vol/storage/5T_disk/organelles \
      --gff-dir /vol/storage/5T_disk/organelles/downloads/genome_fasta \
      --gff-dir /vol/storage/5T_disk/organelles/test_dataset/gff \
      --audit-only --out-prefix dcnuxt_overlap_v2

  # full run
  python3 dcnuxt_annotation_overlap.py \
      --base /vol/storage/5T_disk/organelles \
      --gff-dir /vol/storage/5T_disk/organelles/downloads/genome_fasta \
      --gff-dir /vol/storage/5T_disk/organelles/test_dataset/gff \
      --perms 1000 --out-prefix dcnuxt_overlap_v2
"""

import argparse
import bisect
import csv
import glob
import gzip
import os
import random
import re
import sys
import zlib
from collections import defaultdict, Counter

# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def sanitize(name):
    """Reproduce exactly the awk sanitisation in nupt/numt_gene_alignment*.sh:
         gsub(/[^A-Za-z0-9_.-]/,"_",safe)
       applied to the bedtools getfasta header 'chrom:start-end'."""
    return re.sub(r"[^A-Za-z0-9_.\-]", "_", name)


def smart_open(path):
    """Sniff magic bytes rather than trusting the extension - some files
    carry a .gz suffix but are plain text or truncated."""
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", errors="replace")
    return open(path, "rt", errors="replace")


def sniff_table(path):
    """KaKs_Calculator writes TSV even when the file is named .csv."""
    with open(path, errors="replace") as fh:
        head = fh.readline()
    return "\t" if head.count("\t") >= head.count(",") else ","


class IntervalIndex:
    """Sorted starts + running max-end, for bounded backward scan queries."""

    def __init__(self):
        self._raw = defaultdict(list)
        self._built = {}

    def add(self, seqid, start, end, payload):
        self._raw[seqid].append((start, end, payload))

    def build(self):
        for seqid, items in self._raw.items():
            items.sort(key=lambda x: x[0])
            starts = [i[0] for i in items]
            maxend, running = [], 0
            for _, e, _p in items:
                running = max(running, e)
                maxend.append(running)
            self._built[seqid] = (starts, maxend, items)
        self._raw = None

    def overlaps(self, seqid, qstart, qend):
        """Half-open [qstart, qend). Returns list of (start, end, payload)."""
        entry = self._built.get(seqid)
        if entry is None:
            return []
        starts, maxend, items = entry
        hi = bisect.bisect_left(starts, qend)      # all i < hi have start < qend
        out = []
        for i in range(hi - 1, -1, -1):
            if maxend[i] <= qstart:                # nothing further back can reach
                break
            s, e, p = items[i]
            if e > qstart and s < qend:
                out.append((s, e, p))
        return out


# --------------------------------------------------------------------------
# annotation discovery and loading
# --------------------------------------------------------------------------

GFF_PATTERNS = [
    "{acc}_genomic.gff", "{acc}_genomic.gff3",
    "{acc}*.gff3.gz", "{acc}*.gff.gz", "{acc}*.gff3", "{acc}*.gff",
    "{acc}/*.gff3.gz", "{acc}/*.gff.gz", "{acc}/*.gff3", "{acc}/*.gff",
    "{sp}*.gff3.gz", "{sp}*.gff.gz", "{sp}*.gff3", "{sp}*.gff",
    "{sp}/*.gff3.gz", "{sp}/*.gff.gz", "{sp}/*.gff3", "{sp}/*.gff",
]


def _accession_variants(assembly):
    """The CSV and the filesystem often disagree on assembly version, and
    annotation frequently exists only on the RefSeq mirror of a GenBank
    assembly (or vice versa). Try all of it."""
    base = assembly.split(".")[0]
    out = [assembly, base]
    if assembly.startswith("GCA_"):
        out += ["GCF_" + assembly[4:], "GCF_" + base[4:]]
    elif assembly.startswith("GCF_"):
        out += ["GCA_" + assembly[4:], "GCA_" + base[4:]]
    seen, uniq = set(), []
    for a in out:
        if a not in seen:
            seen.add(a)
            uniq.append(a)
    return uniq


def find_gff(gff_dirs, assembly, species):
    for d in gff_dirs:
        if not os.path.isdir(d):
            continue
        for acc in _accession_variants(assembly):
            for pat in GFF_PATTERNS:
                hits = sorted(glob.glob(os.path.join(d, pat.format(acc=acc, sp=species))))
                hits = [h for h in hits if os.path.getsize(h) > 0]
                if hits:
                    return hits[0]
    return None


def find_assembly_report(gff_dirs, assembly):
    for d in gff_dirs:
        if not os.path.isdir(d):
            continue
        for acc in _accession_variants(assembly):
            for pat in ("{acc}*assembly_report.txt", "{acc}/*assembly_report.txt",
                        "{acc}*sequence_report.jsonl", "{acc}/*sequence_report.jsonl"):
                hits = sorted(glob.glob(os.path.join(d, pat.format(acc=acc))))
                if hits:
                    return hits[0]
    return None


def load_seqid_map(report_path):
    """Every-alias -> every-alias translation.

    NCBI assembly_report.txt columns:
      Sequence-Name, Sequence-Role, Assigned-Molecule, Assigned-Molecule-Location,
      GenBank-Accn, Relationship, RefSeq-Accn, Assembly-Unit, Sequence-Length, ...
    sequence_report.jsonl uses genbankAccession / refseqAccession / chrName.
    """
    aliases = []
    try:
        if report_path.endswith(".jsonl"):
            import json
            with open(report_path, "rt", errors="replace") as fh:
                for line in fh:
                    line = line.strip()
                    if not line:
                        continue
                    try:
                        rec = json.loads(line)
                    except ValueError:
                        continue
                    group = {rec.get(k) for k in ("genbankAccession", "refseqAccession",
                                                  "chrName", "ucscStyleName")}
                    group = {g for g in group if g}
                    if len(group) > 1:
                        aliases.append(group)
        else:
            with open(report_path, "rt", errors="replace") as fh:
                for line in fh:
                    if line.startswith("#") or not line.strip():
                        continue
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 7:
                        continue
                    group = {x for x in (f[0], f[4], f[6]) if x and x != "na"}
                    if len(group) > 1:
                        aliases.append(group)
    except OSError:
        return None

    m = {}
    for group in aliases:
        for a in group:
            m.setdefault(a, set()).update(group)
    return m or None


def gff_seqids(path):
    """Column 1 + ##sequence-region only. Used by --audit-only."""
    seqids, seqlen = set(), {}
    with smart_open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                if line.startswith("##sequence-region"):
                    p = line.split()
                    if len(p) >= 4:
                        try:
                            seqlen[p[1]] = int(p[3])
                        except ValueError:
                            pass
                        seqids.add(p[1])
                continue
            i = line.find("\t")
            if i > 0:
                seqids.add(line[:i])
    return seqids, seqlen


def load_gff(path, wanted_seqids, seqid_map=None):
    """Return (indexes, seq_lengths, translated_flag).

    GFF 1-based closed -> 0-based half-open, matching BED."""
    idx = {k: IntervalIndex() for k in ("CDS", "exon", "gene", "pseudogene")}
    seqlen = {}
    translated = False

    def translate(s):
        nonlocal translated
        if s in wanted_seqids or not seqid_map:
            return s
        for cand in seqid_map.get(s, ()):
            if cand in wanted_seqids:
                translated = True
                return cand
        return s

    with smart_open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                if line.startswith("##sequence-region"):
                    p = line.split()
                    if len(p) >= 4:
                        try:
                            seqlen[translate(p[1])] = int(p[3])
                        except ValueError:
                            pass
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            ftype = f[2]
            if ftype not in ("CDS", "exon", "gene", "pseudogene"):
                continue
            seqid = translate(f[0])
            if seqid not in wanted_seqids:
                continue
            try:
                start, end = int(f[3]) - 1, int(f[4])
            except ValueError:
                continue
            attrs = f[8]

            gene = re.search(r"(?:^|;)gene=([^;]+)", attrs)
            locus = re.search(r"(?:^|;)locus_tag=([^;]+)", attrs)
            fid = re.search(r"(?:^|;)ID=([^;]+)", attrs)
            product = re.search(r"(?:^|;)product=([^;]+)", attrs)
            parent = re.search(r"(?:^|;)Parent=([^;]+)", attrs)

            payload = {
                "type": ftype,
                "label": (gene.group(1) if gene else
                          locus.group(1) if locus else
                          fid.group(1) if fid else "."),
                "product": product.group(1) if product else ".",
                "parent": parent.group(1) if parent else "",
                "pseudo": ("pseudo=true" in attrs or ftype == "pseudogene"),
            }
            idx[ftype].add(seqid, start, end, payload)

    for i in idx.values():
        i.build()
    return idx, seqlen, translated


# --------------------------------------------------------------------------
# dcNUxT coordinate recovery
# --------------------------------------------------------------------------

def load_bed_map(bed_path):
    """safe fragment name -> (chrom, start, end), matching bedtools getfasta."""
    m = {}
    if not os.path.exists(bed_path):
        return m
    with open(bed_path, errors="replace") as fh:
        for line in fh:
            f = line.split()
            if len(f) < 3:
                continue
            try:
                chrom, start, end = f[0], int(f[1]), int(f[2])
            except ValueError:
                continue
            m[sanitize(f"{chrom}:{start}-{end}")] = (chrom, start, end)
    return m


def match_fragment(seq_id, keys, bed_map):
    """seq_id looks like  <safe_frag>_<ACC>__<GENE>[_REV].mafft
    Longest-prefix match against BED-derived keys. Never regex-split:
    chrom names and gene names both contain '_' and digits."""
    stem = seq_id
    for suf in (".mafft", ".fa", ".fasta", ".axt"):
        if stem.endswith(suf):
            stem = stem[: -len(suf)]
    if stem.endswith("_REV"):
        stem = stem[:-4]

    i = bisect.bisect_right(keys, stem)
    for j in range(i - 1, max(-1, i - 200), -1):
        k = keys[j]
        if stem.startswith(k) and (len(stem) == len(k) or stem[len(k)] == "_"):
            chrom, start, end = bed_map[k]
            rest = stem[len(k):].lstrip("_")
            gene = rest.split("__")[-1] if "__" in rest else rest
            return chrom, start, end, gene
    return None


def read_kaks(path, min_ks):
    sep = sniff_table(path)
    out = []
    with open(path, newline="", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter=sep)
        cols = {c.strip().lower(): c for c in (rd.fieldnames or [])}
        c_seq, c_ka = cols.get("sequence"), cols.get("ka")
        c_ks, c_w = cols.get("ks"), cols.get("ka/ks")
        c_p = next((cols[c] for c in cols if c.startswith("p-value")), None)
        if not (c_seq and c_ks):
            return out
        for row in rd:
            try:
                ks = float(row[c_ks])
            except (TypeError, ValueError):
                continue
            if not (ks > min_ks):
                continue

            def num(c):
                try:
                    return float(row[c])
                except (TypeError, ValueError, KeyError):
                    return float("nan")

            out.append({"Sequence": (row[c_seq] or "").strip(),
                        "Ka": num(c_ka), "Ks": ks, "KaKs": num(c_w),
                        "P": num(c_p) if c_p else float("nan")})
    return out


# --------------------------------------------------------------------------
# classification
# --------------------------------------------------------------------------

def classify(idx, chrom, start, end, min_frac):
    """Precedence: CDS > UTR/exon > intronic > pseudogene > intergenic."""
    L = max(1, end - start)

    def best(hits):
        if not hits:
            return None, 0
        top = max(hits, key=lambda h: min(h[1], end) - max(h[0], start))
        return top, min(top[1], end) - max(top[0], start)

    cds, cds_bp = best(idx["CDS"].overlaps(chrom, start, end))
    exon, exon_bp = best(idx["exon"].overlaps(chrom, start, end))
    gene, gene_bp = best(idx["gene"].overlaps(chrom, start, end))
    pgene, pgene_bp = best(idx["pseudogene"].overlaps(chrom, start, end))

    def call(frac):
        if cds and cds_bp / L >= frac:
            return "CDS", cds
        if exon and exon_bp / L >= frac:
            return "UTR_exon", exon
        if gene and gene_bp / L >= frac:
            return "intronic", gene
        if pgene and pgene_bp / L >= frac:
            return "pseudogene", pgene
        return "intergenic", None

    strict, feat = call(min_frac)
    permissive, _ = call(1.0 / L)          # >= 1 bp
    return {
        "class_strict": strict,
        "class_permissive": permissive,
        "cds_bp": cds_bp, "exon_bp": exon_bp,
        "gene_bp": gene_bp, "pseudogene_bp": pgene_bp,
        "feature_label": feat[2]["label"] if feat else ".",
        "feature_product": feat[2]["product"] if feat else ".",
        "feature_pseudo": (feat[2]["pseudo"] if feat else False),
    }


def permute(idx, loci, seqlen, nperm, min_frac, rng):
    """Shuffle each locus uniformly within its own scaffold, controlling for
    the strong scaffold bias in NUxT density."""
    acc = Counter()
    usable = [(c, s, e) for (c, s, e) in loci if seqlen.get(c, 0) > (e - s) + 1]
    if not usable:
        return acc, 0
    for _ in range(nperm):
        for chrom, s, e in usable:
            L = e - s
            ns = rng.randrange(0, seqlen[chrom] - L)
            acc[classify(idx, chrom, ns, ns + L, min_frac)["class_strict"]] += 1
    for k in list(acc):
        acc[k] /= nperm
    return acc, len(usable)


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

# Corrected tables ONLY. The older _kaks.csv / _kaks_c11.csv predate the
# genetic-code and codon-frame fixes; reading them as a fallback would mix
# pre- and post-fix data silently. A species without a fixed table is
# reported as 'no_kaks_table', not substituted.
KAKS_SUFFIXES = {
    "NUPT": ["_kaks_c11_fixed.csv"],
    "NUMT": ["_numt_kaks_fixed.csv"],
}
BED_REL = {
    "NUPT": "nupts/nupts.clean.out_AllFrags.bed",
    "NUMT": "numts/numts.clean.out_AllFrags.bed",
}


def species_dirs(base, species):
    return [d for d in (os.path.join(base, "species", species),
                        os.path.join(base, "multi_organelles", "species", species))
            if os.path.isdir(d)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="/vol/storage/5T_disk/organelles")
    ap.add_argument("--csv", default=None, help="default: <base>/combined_dataset.csv")
    ap.add_argument("--gff-dir", action="append", required=True,
                    help="repeat for each directory that may hold GFFs")
    ap.add_argument("--min-ks", type=float, default=0.1)
    ap.add_argument("--min-frac", type=float, default=0.5,
                    help="min fraction of the dcNUxT inside a feature for the strict call")
    ap.add_argument("--perms", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out-prefix", default="dcnuxt_overlap_v2")
    ap.add_argument("--audit-only", action="store_true")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    csv_path = args.csv or os.path.join(args.base, "combined_dataset.csv")
    rng = random.Random(args.seed)

    rows, audit = [], []
    obs, exp = defaultdict(Counter), defaultdict(Counter)
    unmatched_total = Counter()

    with open(csv_path, newline="", errors="replace") as fh:
        for f in csv.reader(fh):
            if len(f) < 2:
                continue
            species, assembly = f[0].strip(), f[1].strip()
            # robust header skip - do not rely on the literal string "SPECIES"
            if not re.match(r"^GC[AF]_", assembly):
                continue

            sdirs = species_dirs(args.base, species)
            if not sdirs:
                audit.append([species, assembly, "no_species_dir", "", "0", "", ""])
                continue

            # ---- collect dcNUxTs of both kinds ---------------------------
            per_kind = {}
            tables_found = set()
            for kind in ("NUPT", "NUMT"):
                bed_map = {}
                for sd in sdirs:
                    bed_map.update(load_bed_map(os.path.join(sd, BED_REL[kind])))
                if not bed_map:
                    continue
                keys = sorted(bed_map)

                # Exactly ONE table per species per kind. A species can appear
                # in both species/ and multi_organelles/species/, and can have
                # both results/ and results_v2/; without this the same rows
                # were read more than once.
                hits, used = [], None
                for sd in sdirs:
                    for sub in ("results", "results_v2"):
                        for suf in KAKS_SUFFIXES[kind]:
                            p = os.path.join(sd, sub, species + suf)
                            if os.path.exists(p):
                                hits, used = read_kaks(p, args.min_ks), p
                                break
                        if used:
                            break
                    if used:
                        break
                if used:
                    tables_found.add(kind)
                if not hits:
                    continue

                seen, recs, unmatched = set(), [], 0
                for h in hits:
                    m = match_fragment(h["Sequence"], keys, bed_map)
                    if m is None:
                        unmatched += 1
                        continue
                    chrom, start, end, gene = m
                    key = (chrom, start, end, gene)   # collapse duplicate HSP intervals
                    if key in seen:
                        continue
                    seen.add(key)
                    recs.append({**h, "chrom": chrom, "start": start,
                                 "end": end, "gene": gene})
                unmatched_total[kind] += unmatched
                if recs:
                    per_kind[kind] = recs

            if not per_kind:
                status = "no_dcnuxt" if tables_found else "no_kaks_table"
                audit.append([species, assembly, status, "", "0", "",
                              "" if tables_found else "no corrected KaKs table found"])
                continue

            n_dc = sum(len(r) for r in per_kind.values())
            wanted = {r["chrom"] for recs in per_kind.values() for r in recs}

            # ---- annotation ----------------------------------------------
            gff = find_gff(args.gff_dir, assembly, species)
            if gff is None:
                audit.append([species, assembly, "no_gff", "", str(n_dc), "", ""])
                continue

            report = find_assembly_report(args.gff_dir, assembly)
            smap = load_seqid_map(report) if report else None

            try:
                if args.audit_only:
                    gff_ids, seqlen = gff_seqids(gff)
                    covered = len(wanted & gff_ids)
                    translated = False
                    if not covered and smap:
                        alias = {a for c in wanted for a in smap.get(c, ())}
                        covered = len(alias & gff_ids)
                        translated = covered > 0
                    idx = None
                else:
                    idx, seqlen, translated = load_gff(gff, wanted, smap)
                    covered = sum(1 for c in wanted
                                  if c in seqlen or any(c in v._built for v in idx.values()))
            except (OSError, EOFError, zlib.error, gzip.BadGzipFile, UnicodeError) as exc:
                audit.append([species, assembly, "gff_unreadable", os.path.basename(gff),
                              str(n_dc), "", type(exc).__name__])
                continue

            status = "ok" if covered else "SEQID_MISMATCH"
            audit.append([species, assembly, status, os.path.basename(gff), str(n_dc),
                          f"{covered}/{len(wanted)}",
                          "translated" if translated else
                          ("report_absent" if report is None else "report_unused")])

            if not covered or args.audit_only:
                continue

            if args.verbose:
                print(f"  {species}: {n_dc} dcNUxT, gff={os.path.basename(gff)}",
                      file=sys.stderr)

            # ---- classify + permute --------------------------------------
            for kind, recs in per_kind.items():
                loci = []
                for r in recs:
                    res = classify(idx, r["chrom"], r["start"], r["end"], args.min_frac)
                    obs[kind][res["class_strict"]] += 1
                    loci.append((r["chrom"], r["start"], r["end"]))
                    rows.append([
                        species, kind, r["gene"], r["chrom"], r["start"], r["end"],
                        r["end"] - r["start"],
                        f"{r['Ks']:.4f}", f"{r['KaKs']:.4f}", f"{r['P']:.4g}",
                        res["class_strict"], res["class_permissive"],
                        res["cds_bp"], res["exon_bp"], res["gene_bp"],
                        res["feature_label"], res["feature_product"],
                        "yes" if res["feature_pseudo"] else "no",
                        seqlen.get(r["chrom"], ""), os.path.basename(gff),
                    ])
                if args.perms and seqlen:
                    e, _n = permute(idx, loci, seqlen, args.perms, args.min_frac, rng)
                    for k, v in e.items():
                        exp[kind][k] += v

    # ---- write ---------------------------------------------------------
    with open(args.out_prefix + "_gff_audit.tsv", "w") as fh:
        fh.write("species\tassembly\tstatus\tgff\tn_dcnuxt\tseqid_concordance\tnote\n")
        for a in audit:
            fh.write("\t".join(map(str, a)) + "\n")

    st = Counter(a[2] for a in audit)
    dc_by_status = defaultdict(int)
    for a in audit:
        try:
            dc_by_status[a[2]] += int(a[4] or 0)
        except ValueError:
            pass

    print("GFF audit:")
    for k, v in st.most_common():
        print(f"  {k:18s} {v:4d} species   {dc_by_status[k]:6d} dcNUxTs")
    if sum(unmatched_total.values()):
        print(f"  [warn] KaKs rows with no BED coordinate match: "
              f"NUPT {unmatched_total['NUPT']}, NUMT {unmatched_total['NUMT']}")

    if args.audit_only:
        print(f"\nWrote {args.out_prefix}_gff_audit.tsv")
        return

    with open(args.out_prefix + "_loci.tsv", "w") as fh:
        fh.write("\t".join([
            "species", "kind", "gene", "chrom", "start", "end", "len",
            "Ks", "KaKs", "P", "class_strict", "class_permissive",
            "cds_bp", "exon_bp", "gene_bp", "feature", "product",
            "feature_is_pseudogene", "scaffold_len", "gff"]) + "\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")

    order = ("CDS", "UTR_exon", "intronic", "pseudogene", "intergenic")
    with open(args.out_prefix + "_summary.tsv", "w") as fh:
        fh.write("kind\tclass\tobserved\tpct\texpected_perm\tobs_over_exp\n")
        for kind in ("NUPT", "NUMT"):
            tot = sum(obs[kind].values())
            if not tot:
                continue
            for cls in order:
                o, e = obs[kind][cls], exp[kind].get(cls, 0.0)
                oe = f"{o / e:.2f}" if e > 0 else "NA"
                fh.write(f"{kind}\t{cls}\t{o}\t{100 * o / tot:.1f}\t{e:.1f}\t{oe}\n")
            fh.write(f"{kind}\tTOTAL\t{tot}\t100.0\t\t\n")

    for kind in ("NUPT", "NUMT"):
        tot = sum(obs[kind].values())
        if not tot:
            continue
        print(f"\n{kind}  n={tot}")
        for cls in order:
            o, e = obs[kind][cls], exp[kind].get(cls, 0.0)
            oe = f"   o/e {o / e:.2f}" if e > 0 else ""
            print(f"  {cls:12s} obs {o:5d} ({100 * o / tot:5.1f}%)   exp {e:7.1f}{oe}")

    print(f"\nWrote {args.out_prefix}_loci.tsv, _summary.tsv, _gff_audit.tsv")


if __name__ == "__main__":
    sys.exit(main())
