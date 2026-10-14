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
                     PSEUDOGENE / NO_GENE_STRUCTURE / ORGANELLAR_CONTIG

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
from collections import defaultdict, Counter

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
    cds_by_gene = defaultdict(list)

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
            elif ftype.endswith("RNA") or ftype in ("transcript", "primary_transcript",
                                                    "mRNA", "V_gene_segment"):
                if parent and fid and parent.group(1) in genes:
                    tx2gene[fid.group(1)] = parent.group(1)
                    if product and genes[parent.group(1)]["product"] == ".":
                        genes[parent.group(1)]["product"] = product.group(1)
            elif ftype == "CDS":
                if parent:
                    for p in parent.group(1).split(","):
                        g = tx2gene.get(p) or (p if p in genes else None)
                        if g:
                            cds_by_gene[g].append((int(f[3]) - 1, int(f[4])))
                            if product and genes[g]["product"] == ".":
                                genes[g]["product"] = product.group(1)

    # second pass for exons. Parent may point at a transcript OR directly at
    # the gene (e.g. "Parent=gene-ERUC_LOCUS56;gbkey=exon"), and may precede
    # the line it refers to.
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
                g = tx2gene.get(p) or (p if p in genes else None)
                if g:
                    genes[g]["tx"][p].append((int(f[3]) - 1, int(f[4])))

    # genes with no exon lines at all: use CDS segments as the exon proxy.
    # n_exons==0 then means "no structure in the GFF", which is recorded
    # separately rather than silently treated as a single-exon gene.
    for gid, g in genes.items():
        if not g["tx"] and cds_by_gene.get(gid):
            g["tx"]["_cds"] = sorted(set(cds_by_gene[gid]))
            g["exon_source"] = "CDS"
        else:
            g["exon_source"] = "exon" if g["tx"] else "none"

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
    ap.add_argument("--targets", default="CANDIDATE",
                    help="comma-separated verdicts to extract flanked sequence "
                         "for, or 'all'. Transcription is an empirical question "
                         "independent of annotation quality: pseudogene, intronic "
                         "and circular loci can all be transcribed, so widen this "
                         "rather than pre-judging on the verdict.")
    ap.add_argument("--exclude-scaffolds", default=None,
                    help="file of scaffold IDs verified as organellar contigs")
    ap.add_argument("--out-prefix", default="dcnuxt_overlap")
    args = ap.parse_args()

    ALL_VERDICTS = {"CANDIDATE","INTRONIC_ONLY","SINGLE_EXON","CIRCULAR","PSEUDOGENE","ORGANELLAR_CONTIG","NO_GENE_STRUCTURE"}
    targets = (ALL_VERDICTS if args.targets.strip().lower() == "all"
               else {t.strip().upper() for t in args.targets.split(",") if t.strip()})
    unknown = targets - ALL_VERDICTS
    if unknown:
        sys.exit(f"unknown verdict(s) in --targets: {', '.join(sorted(unknown))}\n"
                 f"choose from: {', '.join(sorted(ALL_VERDICTS))}, or 'all'")

    excluded = set()
    if args.exclude_scaffolds and os.path.exists(args.exclude_scaffolds):
        with open(args.exclude_scaffolds) as fh:
            excluded = {ln.split()[0] for ln in fh if ln.strip()}

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
        pending = []

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
            exon_src = g.get("exon_source", "none") if g else "none"

            pending.append({
                "row": r, "start": start, "end": end, "L": L, "scaf": scaf,
                "host": host, "gtype": gtype, "product": product,
                "n_exons": n_exons, "hits_exon": hits_exon, "exon_bp": exon_bp,
                "gene_span": gene_span, "frac": frac, "exon_src": exon_src,
                "dlike": dlike, "generic": generic,
            })

        # Organellar evidence is a property of the scaffold, not the locus.
        # Stylosanthes scaffold 46011 carries both a "hypothetical protein"
        # host and RPOC1_1 (plastid RNA polymerase): the plastid gene name
        # condemns the whole contig, so propagate it.
        organellar_scaf = {p["row"]["chrom"] for p in pending
                           if p["dlike"] and p["scaf"]
                           and p["scaf"] < ORGANELLAR_SCAFFOLD_MAX}

        for p in pending:
            r = p["row"]
            if r["chrom"] in excluded or r["chrom"] in organellar_scaf:
                verdict = "ORGANELLAR_CONTIG"
            elif p["gtype"] == "pseudogene":
                verdict = "PSEUDOGENE"
            elif p["dlike"]:
                verdict = "CIRCULAR"
            elif p["exon_src"] == "none":
                verdict = "NO_GENE_STRUCTURE"
            elif p["n_exons"] <= 1:
                verdict = "SINGLE_EXON"
            elif not p["hits_exon"]:
                verdict = "INTRONIC_ONLY"
            elif p["frac"] > 0.8:
                verdict = "CIRCULAR"
            else:
                verdict = "CANDIDATE"

            out_rows.append([
                species, r["kind"], r["gene"], r["chrom"], p["start"], p["end"], p["L"],
                r["Ks"], r["KaKs"], r["P"], r["class_strict"],
                p["host"], p["gtype"], p["product"], p["n_exons"],
                "yes" if p["hits_exon"] else "no", p["exon_bp"],
                p["gene_span"], f'{p["frac"]:.3f}', p["exon_src"],
                "yes" if p["dlike"] else "no", "yes" if p["generic"] else "no",
                p["scaf"], verdict,
            ])

            if verdict in targets:
                fs = max(0, p["start"] - args.flank)
                bed_rows.append([species, r["chrom"], fs, p["end"] + args.flank,
                                 f'{species}|{r["kind"]}|{r["gene"]}|{p["host"]}|{verdict}',
                                 p["start"], p["end"]])

    hdr = ["species", "kind", "donor_gene", "chrom", "start", "end", "len",
           "Ks", "KaKs", "P", "class", "host_gene", "host_type", "host_product",
           "n_exons", "hits_exon", "exon_overlap_bp", "gene_span",
           "frac_gene_is_nuxt", "exon_source", "donor_like_annotation",
           "generic_product", "scaffold_len", "verdict"]
    with open(args.out_prefix + "_exon_cluster.tsv", "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for r in out_rows:
            fh.write("\t".join(map(str, r)) + "\n")

    with open(args.out_prefix + "_rna_targets.bed", "w") as fh:
        for sp, chrom, s, e, name, ts, te in bed_rows:
            fh.write(f"{chrom}\t{s}\t{e}\t{name}\t0\t+\n")

    # optional sequence extraction for junction-spanning read checks.
    # Reads that fall wholly inside the insertion are ambiguous with the
    # organellar transcript; only reads crossing insertion<->flank are
    # informative, which is what the flanks are for.
    if args.flank and args.genome_dir:
        written, skipped = 0, Counter()
        gen_cache = {}
        with open(args.out_prefix + "_rna_targets.fa", "w") as out:
            for sp, chrom, s0, e0, name, ts, te in bed_rows:
                a = asm.get(sp)
                if a not in gen_cache:
                    gen_cache[a] = find_one(
                        args.genome_dir, a,
                        ["{acc}_*_genomic.fna", "{acc}_genomic.fna",
                         "{acc}_*_genomic.fna.gz"]) if a else None
                gen = gen_cache[a]
                if not gen:
                    skipped["no_genome_fasta"] += 1
                    continue
                if gen.endswith(".gz"):
                    skipped["genome_gzipped"] += 1
                    continue
                if not os.path.exists(gen + ".fai"):
                    try:
                        subprocess.run(["samtools", "faidx", gen],
                                       capture_output=True, check=True)
                    except (subprocess.CalledProcessError, FileNotFoundError):
                        skipped["faidx_failed"] += 1
                        continue
                try:
                    res = subprocess.run(
                        ["samtools", "faidx", gen, f"{chrom}:{s0+1}-{e0}"],
                        capture_output=True, text=True, check=True).stdout
                except (subprocess.CalledProcessError, FileNotFoundError):
                    skipped["faidx_region_failed"] += 1
                    continue
                body = "".join(res.split("\n")[1:])
                if not body:
                    skipped["empty_region"] += 1
                    continue
                # offsets of the two junctions within the extracted sequence
                j1, j2 = ts - s0, te - s0
                out.write(f">{name} {chrom}:{s0+1}-{e0} insert={ts}-{te} "
                          f"junctions={j1},{j2}\n{body}\n")
                written += 1
        print(f"extracted {written} flanked sequences")
        for k, n in skipped.most_common():
            print(f"  skipped ({k}): {n}")

    v = Counter(r[-1] for r in out_rows)
    print(f"\n{len(out_rows)} dcNUxT loci inside an annotated gene")
    for k, n in v.most_common():
        print(f"  {k:20s} {n}")
    print(f"\nWrote {args.out_prefix}_exon_cluster.tsv, _rna_targets.bed")


if __name__ == "__main__":
    sys.exit(main())
