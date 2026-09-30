#!/usr/bin/env python3
"""
dcnuxt_exon_cluster.py

Follow-up on dcNUxT loci that fall inside an annotated gene: reconstruct the
host gene's structure and decide whether the overlap is credible evidence of
integration into a pre-existing nuclear gene, or an artefact.

Input : <prefix>_loci.tsv from dcnuxt_annotation_overlap.py
Output: <prefix>_exon_cluster.tsv   one row per locus inside a gene
        <prefix>_rna_targets.bed    coordinates for RNA-seq validation
        <prefix>_rna_targets.fa     locus + flanking sequence (with --flank)

What it reports per locus
-------------------------
  n_exons            exons in the host gene's longest transcript.
                     1 exon is a red flag: ab initio predictors call an
                     intact NUxT ORF a single-exon gene, so the "gene" is
                     the insertion itself, not a pre-existing gene.
  hits_exon          does the dcNUxT overlap a real exon, or only an intron?
  frac_gene_is_nuxt  fraction of the host gene span covered by the dcNUxT.
                     Near 1.0 means the gene IS the insertion (circular).
  donor_like         host product/name restates the donor organellar gene
                     (atpA -> "ATP synthase alpha", nad3 -> "NADH ...").
                     Homology-driven annotation, not independent evidence.
  verdict            CANDIDATE / CIRCULAR / SINGLE_EXON / INTRONIC_ONLY /
                     ORGANELLAR_CONTIG

Usage
-----
  python3 dcnuxt_exon_cluster.py \
      --loci dcnuxt_overlap_loci.tsv \
      --csv  /vol/storage/5T_disk/organelles/combined_dataset.csv \
      --gff-dir /vol/storage/5T_disk/organelles/downloads/genome_fasta \
      --genome-dir /vol/storage/5T_disk/organelles/downloads/genome_fasta \
      --flank 300 \
      --out-prefix dcnuxt_overlap
"""

import argparse
import csv
import glob
import gzip
import os
import re
import subprocess
import sys
from collections import defaultdict

# scaffolds below this are treated as probable unplaced organellar contigs
ORGANELLAR_SCAFFOLD_MAX = 200_000

# donor gene symbol -> words that would appear in a homology-driven annotation
DONOR_KEYWORDS = {
    "atpa": ["atp synthase", "atpase", "atp1"],
    "atp1": ["atp synthase", "atpase", "atpa"],
    "atpb": ["atp synthase", "atpase"],
    "atpf": ["atp synthase"], "atph": ["atp synthase"], "atpi": ["atp synthase"],
    "nad1": ["nadh"], "nad2": ["nadh"], "nad3": ["nadh"], "nad4": ["nadh"],
    "nad5": ["nadh"], "nad6": ["nadh"], "nad7": ["nadh"], "nad9": ["nadh"],
    "ndhf": ["nadh"], "ndhj": ["nadh"], "ndha": ["nadh"], "ndhb": ["nadh"],
    "rbcl": ["ribulose", "rubisco"],
    "psba": ["photosystem"], "psbb": ["photosystem"], "psbc": ["photosystem"],
    "psbd": ["photosystem"], "psaa": ["photosystem"], "psab": ["photosystem"],
    "rpob": ["rna polymerase"], "rpoc1": ["rna polymerase"],
    "rpoc2": ["rna polymerase"], "rpoa": ["rna polymerase"],
    "ycf1": ["ycf1"], "ycf2": ["ycf2"],
    "accd": ["acetyl-coa", "carboxyl"],
    "peta": ["cytochrome"], "petb": ["cytochrome"], "petd": ["cytochrome"],
    "cob": ["cytochrome"], "cox1": ["cytochrome"], "cox2": ["cytochrome"],
    "cox3": ["cytochrome"], "ccmfc": ["cytochrome"], "ccmfn": ["cytochrome"],
    "matr": ["maturase"], "sdh3": ["succinate"], "sdh4": ["succinate"],
    "cema": ["cema"], "clpp": ["clp protease"],
    "rpl16": ["ribosomal protein l16"], "rpl22": ["ribosomal protein l22"],
    "rpl2": ["ribosomal protein l2"], "rps4": ["ribosomal protein s4"],
    "rps12": ["ribosomal protein s12"],
}

GENERIC = ("hypothetical protein", "uncharacterized", "unknown", "predicted protein", ".")


def smart_open(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", errors="replace")
    return open(path, "rt", errors="replace")


def accession_variants(assembly):
    base = assembly.split(".")[0]
    out = [assembly, base]
    if assembly.startswith("GCA_"):
        out += ["GCF_" + assembly[4:], "GCF_" + base[4:]]
    elif assembly.startswith("GCF_"):
        out += ["GCA_" + assembly[4:], "GCA_" + base[4:]]
    seen, uniq = set(), []
    for a in out:
        if a not in seen:
            seen.add(a); uniq.append(a)
    return uniq


def find_one(dirs, assembly, patterns):
    for d in dirs:
        if not os.path.isdir(d):
            continue
        for acc in accession_variants(assembly):
            for pat in patterns:
                hits = sorted(glob.glob(os.path.join(d, pat.format(acc=acc))))
                hits = [h for h in hits if os.path.getsize(h) > 0]
                if hits:
                    return hits[0]
    return None


def load_gene_models(gff_path, wanted_seqids):
    """Build gene -> {span, exons per transcript, product} for every gene on
    the scaffolds we care about. Matching is done by coordinate, not by the
    label in the loci file: submitter-supplied GFFs frequently omit gene= and
    locus_tag= from exon lines, which leaves the label empty."""
    genes = {}
    tx2gene = {}

    with smart_open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] not in wanted_seqids:
                continue
            ftype, attrs = f[2], f[8]
            fid = re.search(r"(?:^|;)ID=([^;]+)", attrs)
            parent = re.search(r"(?:^|;)Parent=([^;]+)", attrs)
            gene = re.search(r"(?:^|;)gene=([^;]+)", attrs)
            locus = re.search(r"(?:^|;)locus_tag=([^;]+)", attrs)
            product = re.search(r"(?:^|;)product=([^;]+)", attrs)
            label = (gene.group(1) if gene else
                     locus.group(1) if locus else
                     (fid.group(1) if fid else "."))

            if ftype in ("gene", "pseudogene"):
                gid = fid.group(1) if fid else f"{f[0]}:{f[3]}-{f[4]}"
                genes[gid] = {
                    "label": label, "seqid": f[0],
                    "start": int(f[3]) - 1, "end": int(f[4]),
                    "type": ftype, "product": ".",
                    "tx": defaultdict(list),
                }
            elif ftype in ("mRNA", "transcript", "lnc_RNA", "ncRNA",
                           "primary_transcript"):
                if parent and fid and parent.group(1) in genes:
                    tx2gene[fid.group(1)] = parent.group(1)
                    if product and genes[parent.group(1)]["product"] == ".":
                        genes[parent.group(1)]["product"] = product.group(1)
            elif ftype == "CDS":
                if parent and product:
                    for p in parent.group(1).split(","):
                        g = tx2gene.get(p)
                        if g and genes[g]["product"] == ".":
                            genes[g]["product"] = product.group(1)

    # second pass for exons (Parent can precede the mRNA line)
    with smart_open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon" or f[0] not in wanted_seqids:
                continue
            parent = re.search(r"(?:^|;)Parent=([^;]+)", f[8])
            if not parent:
                continue
            for p in parent.group(1).split(","):
                g = tx2gene.get(p)
                if g:
                    genes[g]["tx"][p].append((int(f[3]) - 1, int(f[4])))

    by_seqid = defaultdict(list)
    for g in genes.values():
        by_seqid[g["seqid"]].append(g)
    return by_seqid


def donor_like(donor_gene, label, product):
    d = (donor_gene or "").lower().strip()
    text = f"{label} {product}".lower()
    if d and d in text:
        return True
    for kw in DONOR_KEYWORDS.get(d, []):
        if kw in text:
            return True
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--loci", required=True)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--gff-dir", action="append", required=True)
    ap.add_argument("--genome-dir", action="append", default=[])
    ap.add_argument("--flank", type=int, default=0,
                    help="bp of flanking sequence to extract for RNA validation")
    ap.add_argument("--out-prefix", default="dcnuxt_overlap")
    args = ap.parse_args()

    asm = {}
    with open(args.csv, newline="", errors="replace") as fh:
        for f in csv.reader(fh):
            if len(f) >= 2 and re.match(r"^GC[AF]_", f[1].strip()):
                asm[f[0].strip()] = f[1].strip()

    INSIDE = {"CDS", "UTR_exon", "intronic", "pseudogene"}
    by_species = defaultdict(list)
    with open(args.loci, newline="", errors="replace") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for row in rd:
            if row["class_strict"] in INSIDE or row["class_permissive"] in INSIDE:
                by_species[row["species"]].append(row)

    out_rows, bed_rows = [], []

    for species, rows in sorted(by_species.items()):
        assembly = asm.get(species)
        if not assembly:
            continue
        gff = find_one(args.gff_dir, assembly,
                       ["{acc}_genomic.gff", "{acc}_genomic.gff3",
                        "{acc}*.gff3.gz", "{acc}*.gff.gz", "{acc}*.gff3", "{acc}*.gff"])
        if not gff:
            continue

        wanted_seqids = {r["chrom"] for r in rows}
        by_seqid = load_gene_models(gff, wanted_seqids)

        for r in rows:
            start, end = int(r["start"]), int(r["end"])
            L = end - start
            scaf = int(r["scaffold_len"]) if r["scaffold_len"] else 0

            cand = [x for x in by_seqid.get(r["chrom"], ())
                    if x["end"] > start and x["start"] < end]
            g = max(cand, key=lambda x: x["end"] - x["start"]) if cand else None

            if g:
                best_tx = max(g["tx"].values(), key=len) if g["tx"] else []
                n_exons = len(best_tx)
                exon_bp = sum(max(0, min(e, end) - max(s, start)) for s, e in best_tx)
                gene_span = g["end"] - g["start"]
                frac = L / gene_span if gene_span else 0.0
                product = g["product"] if g["product"] != "." else r["product"]
                host = g["label"]
                gtype = g["type"]
            else:
                n_exons, exon_bp, gene_span, frac = 0, 0, 0, 0.0
                product, host, gtype = r["product"], r["feature"], "?"

            hits_exon = exon_bp > 0
            dlike = donor_like(r["gene"], host, product)
            generic = any(x in product.lower() for x in GENERIC)

            if scaf and scaf < ORGANELLAR_SCAFFOLD_MAX:
                verdict = "ORGANELLAR_CONTIG"
            elif dlike:
                verdict = "CIRCULAR"
            elif n_exons <= 1:
                verdict = "SINGLE_EXON"
            elif not hits_exon:
                verdict = "INTRONIC_ONLY"
            elif frac > 0.8:
                verdict = "CIRCULAR"
            else:
                verdict = "CANDIDATE"

            out_rows.append([
                species, r["kind"], r["gene"], r["chrom"], start, end, L,
                r["Ks"], r["KaKs"], r["P"], r["class_strict"],
                host, gtype, product, n_exons,
                "yes" if hits_exon else "no", exon_bp,
                gene_span, f"{frac:.3f}",
                "yes" if dlike else "no", "yes" if generic else "no",
                scaf, verdict,
            ])

            if verdict == "CANDIDATE":
                fs = max(0, start - args.flank)
                bed_rows.append([species, r["chrom"], fs, end + args.flank,
                                 f"{species}|{r['gene']}|{r['feature']}", start, end])

    hdr = ["species", "kind", "donor_gene", "chrom", "start", "end", "len",
           "Ks", "KaKs", "P", "class", "host_gene", "host_type", "host_product",
           "n_exons", "hits_exon", "exon_overlap_bp", "gene_span",
           "frac_gene_is_nuxt", "donor_like_annotation", "generic_product",
           "scaffold_len", "verdict"]
    with open(args.out_prefix + "_exon_cluster.tsv", "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for r in out_rows:
            fh.write("\t".join(map(str, r)) + "\n")

    with open(args.out_prefix + "_rna_targets.bed", "w") as fh:
        for sp, chrom, s, e, name, ts, te in bed_rows:
            fh.write(f"{chrom}\t{s}\t{e}\t{name}\t0\t+\n")

    # optional sequence extraction for junction-spanning read checks
    if args.flank and args.genome_dir:
        written = 0
        with open(args.out_prefix + "_rna_targets.fa", "w") as out:
            for sp, chrom, s, e, name, ts, te in bed_rows:
                a = asm.get(sp)
                gen = find_one(args.genome_dir, a,
                               ["{acc}_*_genomic.fna", "{acc}_genomic.fna",
                                "{acc}_*_genomic.fna.gz"]) if a else None
                if not gen or gen.endswith(".gz"):
                    continue
                if not os.path.exists(gen + ".fai"):
                    continue
                try:
                    seq = subprocess.run(
                        ["samtools", "faidx", gen, f"{chrom}:{s+1}-{e}"],
                        capture_output=True, text=True, check=True).stdout
                except (subprocess.CalledProcessError, FileNotFoundError):
                    continue
                body = "".join(seq.split("\n")[1:])
                out.write(f">{name} {chrom}:{s+1}-{e} insert={ts}-{te}\n{body}\n")
                written += 1
        print(f"extracted {written} flanked sequences")

    from collections import Counter
    v = Counter(r[-1] for r in out_rows)
    print(f"\n{len(out_rows)} dcNUxT loci inside an annotated gene")
    for k, n in v.most_common():
        print(f"  {k:20s} {n}")
    print(f"\nWrote {args.out_prefix}_exon_cluster.tsv, _rna_targets.bed")


if __name__ == "__main__":
    sys.exit(main())
