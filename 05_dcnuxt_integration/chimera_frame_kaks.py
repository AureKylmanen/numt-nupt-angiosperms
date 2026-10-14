#!/usr/bin/env python3
"""
chimera_frame_kaks.py

For a dcNUxT that sits inside a host CDS and is shared across species,
test whether the insertion-derived part of the protein is under the same
constraint as the rest of the gene.

Why this is not the Ka/Ks already in the results table
------------------------------------------------------
The pipeline's omega compares the insertion to its ORGANELLAR donor, in the
ORGANELLAR reading frame. Once an insertion is incorporated into a nuclear
gene it is translated in a frame set by the host gene, which is generally a
different frame. Selection, if any, acts on the nuclear frame. So the
donor-based omega does not test the functional hypothesis at all.

What this does
--------------
 1. Extract the reference CDS (exons concatenated, strand-corrected).
 2. BLAST it against each comparison genome, take the best locus, extract
    the orthologous region.
 3. Use BLAST's own aligned qseq/sseq - the query is the reference CDS in
    frame, so subject bases map straight onto reference codons. No frame
    inference, which is where the earlier version went wrong.
 4. Split the alignment into insertion-derived and host-derived columns.
 5. Run KaKs_Calculator with -c 1 (standard code - the NUCLEAR frame) on
    each partition separately, for every pairwise comparison.

Reading the result
------------------
  omega(insertion) ~ omega(rest)   -> the insertion is constrained as part of
                                      the protein: evidence of co-option
  omega(insertion) >> omega(rest)  -> passively carried inside a transcribed
                                      gene; retained, but not maintained by
                                      selection on its own sequence
  omega(insertion) near 1          -> neutral

Usage
-----
  python3 chimera_frame_kaks.py \
     --ref-species Salix_udensis \
     --ref-chrom CM055022.1 \
     --cds 3001243-3001494 --cds 3000330-3000926 \
     --strand - \
     --insertion 3000332-3000653 \
     --compare Populus_tomentosa --compare Populus_trichocarpa \
     --out salix_sdh4_chimera
"""

import argparse
import csv
import glob
import os
import re
import subprocess
import sys
import tempfile

CODON_STOP = {"TAA", "TAG", "TGA"}
COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def rc(s):
    return s.translate(COMP)[::-1]


def sh(cmd, **kw):
    return subprocess.run(cmd, shell=isinstance(cmd, str), capture_output=True,
                          text=True, check=True, **kw).stdout


def genome_for(base, species):
    acc = None
    with open(os.path.join(base, "combined_dataset.csv"), errors="replace") as fh:
        for f in csv.reader(fh):
            if len(f) >= 2 and f[0].strip() == species:
                acc = f[1].strip()
                break
    if not acc:
        return None
    for pat in (f"{acc}_*_genomic.fna", f"{acc}_genomic.fna"):
        hits = sorted(glob.glob(os.path.join(base, "downloads/genome_fasta", pat)))
        if hits:
            return hits[0]
    return None


def faidx(genome, region):
    out = sh(["samtools", "faidx", genome, region])
    return "".join(out.split("\n")[1:]).upper()


def translate(seq):
    tbl = {}
    bases = "TCAG"
    aas = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    i = 0
    for a in bases:
        for b in bases:
            for c in bases:
                tbl[a + b + c] = aas[i]
                i += 1
    return "".join(tbl.get(seq[j:j + 3], "X") for j in range(0, len(seq) - 2, 3))


def build_ref_cds(genome, chrom, blocks, strand):
    """blocks: list of (start,end) 1-based inclusive. Returns (cds, mapping)
    where mapping[i] = genomic coordinate of CDS base i (0-based index)."""
    pieces, coords = [], []
    for s, e in sorted(blocks):
        pieces.append(faidx(genome, f"{chrom}:{s}-{e}"))
        coords.extend(range(s, e + 1))
    cds = "".join(pieces)
    if strand == "-":
        cds = rc(cds)
        coords = coords[::-1]
    return cds, coords


def aligned_pairs(query_fa, genome, restrict=None, min_cov=0.0, ref_len=0,
                  task="megablast", word_size=None):
    """Return {ref_cds_index0: subject_base} from the best BLASTN HSP.

    Using BLAST's own qseq/sseq removes all frame guessing: the query IS the
    reference CDS in frame, so mapping subject bases onto query positions
    yields a codon-aware alignment directly. The previous version extracted
    a flanked region and re-derived the frame by "longest stop-free run",
    which silently produced out-of-frame comparisons (Ka ~ 4, i.e. saturated,
    despite 88% nucleotide identity).
    """
    cmd = ["blastn", "-query", query_fa, "-subject", genome, "-task", task,
           "-outfmt", "6 bitscore pident length qstart qend sseqid sstart send qseq sseq",
           "-evalue", "1e-5", "-max_hsps", "1"]
    if word_size:
        cmd += ["-word_size", str(word_size)]
    out = sh(cmd)
    cands = []
    for line in out.strip().split("\n"):
        if not line:
            continue
        f = line.split("\t")
        if len(f) < 10:
            continue
        if restrict and f[5] != restrict:
            continue
        cov = (int(f[4]) - int(f[3]) + 1) / ref_len if ref_len else 1.0
        if cov < min_cov:
            continue
        cands.append((float(f[0]), f))
    if not cands:
        return None
    best = max(cands, key=lambda x: x[0])
    _, f = best
    pident, length, qs, qe, sid, qseq, sseq = (
        float(f[1]), int(f[2]), int(f[3]), int(f[4]), f[5], f[8].upper(), f[9].upper())

    mapping = {}
    qpos = qs - 1                      # 0-based index into the reference CDS
    for qc, sc in zip(qseq, sseq):
        if qc == "-":                  # insertion in subject: no ref position
            continue
        mapping[qpos] = sc
        qpos += 1
    return {"seqid": sid, "pident": pident, "length": length,
            "qstart": qs, "qend": qe, "map": mapping}


def partition(ref_cds, ref_coords, hits, ins_lo, ins_hi):
    """Build per-comparison codon strings, split by whether the reference
    codon overlaps the insertion. Only codons where the reference AND the
    comparison both have three real bases are kept."""
    ncod = len(ref_cds) // 3
    cls = []
    for c in range(ncod):
        gs = ref_coords[c * 3:c * 3 + 3]
        cls.append("INS" if any(ins_lo <= g <= ins_hi for g in gs) else "HOST")

    parts = {}
    for sp, h in hits.items():
        m = h["map"]
        parts[sp] = {"INS": ["", ""], "HOST": ["", ""]}
        for c in range(ncod):
            idx = [c * 3, c * 3 + 1, c * 3 + 2]
            rc_ = ref_cds[idx[0]:idx[0] + 3]
            sc_ = "".join(m.get(i, "-") for i in idx)
            if len(rc_) != 3 or "-" in sc_ or "N" in sc_ or "N" in rc_:
                continue
            if rc_ in CODON_STOP or sc_ in CODON_STOP:
                continue
            parts[sp][cls[c]][0] += rc_
            parts[sp][cls[c]][1] += sc_
    return parts, cls



def count_diffs(a, b):
    """Raw observed codon differences, split into synonymous and
    nonsynonymous. Ka/Ks and any ratio test built on them are meaningless
    when these counts are in the single digits, which is common at the
    divergences typical of within-genus comparisons."""
    syn = nonsyn = 0
    for i in range(0, min(len(a), len(b)) - 2, 3):
        ca, cb = a[i:i + 3], b[i:i + 3]
        if len(ca) < 3 or len(cb) < 3 or ca == cb:
            continue
        aa_a, aa_b = translate(ca), translate(cb)
        d = sum(1 for x, y in zip(ca, cb) if x != y)
        if aa_a == aa_b:
            syn += d
        else:
            nonsyn += d
    return syn, nonsyn


def run_kaks(a_name, a_seq, b_name, b_seq, tag, workdir, kaks_bin):
    """Pairwise KaKs_Calculator in the NUCLEAR frame (-c 1)."""
    A, B = a_seq, b_seq
    if len(A) < 60:
        return None, len(A) // 3
    axt = os.path.join(workdir, f"{tag}.axt")
    with open(axt, "w") as fh:
        fh.write(f"{tag}\n{A}\n{B}\n\n")
    outf = os.path.join(workdir, f"{tag}.kaks")
    try:
        subprocess.run([kaks_bin, "-i", axt, "-o", outf, "-m", "YN", "-c", "1"],
                       capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        return f"ERROR {e}", len(A) // 3
    with open(outf) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return (rows[0] if rows else None), len(A) // 3


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default="/vol/storage/5T_disk/organelles")
    ap.add_argument("--ref-species", required=True)
    ap.add_argument("--ref-chrom", required=True)
    ap.add_argument("--cds", action="append", required=True,
                    help="start-end, 1-based inclusive; repeat per CDS block")
    ap.add_argument("--strand", default="+", choices=["+", "-"])
    ap.add_argument("--insertion", required=True, help="start-end in genome coords")
    ap.add_argument("--compare", action="append", required=True)
    ap.add_argument("--blast-task", default="megablast",
                    choices=["megablast", "dc-megablast", "blastn"],
                    help="megablast (default) only seeds above ~90%% identity. "
                         "For cross-genus or cross-family comparisons use "
                         "dc-megablast, which is designed for ~70-85%% identity "
                         "and finds hits megablast misses entirely.")
    ap.add_argument("--word-size", type=int, default=None,
                    help="blastn word size; try 7-9 with --blast-task blastn "
                         "if dc-megablast still fails")
    ap.add_argument("--restrict", action="append", default=[],
                    help="SPECIES:SEQID - pin the BLAST target for that species. "
                         "Essential for haplotype-resolved or polyploid assemblies, "
                         "where the best hit is often the other haplotype or a "
                         "different NUxT copy rather than the orthologous locus.")
    ap.add_argument("--min-coverage", type=float, default=0.5,
                    help="reject HSPs covering less than this fraction of the "
                         "reference CDS (default 0.5)")
    ap.add_argument("--kaks-bin", default=os.environ.get("KAKS_BIN", "KaKs_Calculator"),
                    help="path to KaKs_Calculator (or set $KAKS_BIN); it is often "
                         "not on PATH, e.g. /vol/storage/software/KaKs_Calculator2.0/"
                         "bin/Linux/KaKs_Calculator")
    ap.add_argument("--out", default="chimera_kaks")
    args = ap.parse_args()

    for t in ("samtools", "blastn"):
        if not any(os.access(os.path.join(p, t), os.X_OK)
                   for p in os.environ["PATH"].split(":")):
            sys.exit(f"ERROR: {t} not in PATH")
    kaks = args.kaks_bin
    if not os.access(kaks, os.X_OK):
        found = [os.path.join(p, kaks) for p in os.environ["PATH"].split(":")
                 if os.access(os.path.join(p, kaks), os.X_OK)]
        if not found:
            sys.exit(f"ERROR: KaKs_Calculator not found at '{kaks}'.\n"
                     f"       pass --kaks-bin /path/to/KaKs_Calculator or set $KAKS_BIN")
        kaks = found[0]

    ref_gen = genome_for(args.base, args.ref_species)
    if not ref_gen:
        sys.exit(f"no genome for {args.ref_species}")

    blocks = [tuple(int(x) for x in b.split("-")) for b in args.cds]
    ins_lo, ins_hi = (int(x) for x in args.insertion.split("-"))

    ref_cds, ref_coords = build_ref_cds(ref_gen, args.ref_chrom, blocks, args.strand)
    prot = translate(ref_cds)
    n_ins = sum(1 for g in ref_coords if ins_lo <= g <= ins_hi)
    print(f"reference CDS {len(ref_cds)} bp ({len(ref_cds)//3} codons), "
          f"{n_ins} bp insertion-derived")
    print(f"internal stops in reference: {prot.rstrip('*').count('*')}")

    work = tempfile.mkdtemp(prefix="chimera_")
    qfa = os.path.join(work, "ref_cds.fa")
    with open(qfa, "w") as fh:
        fh.write(f">ref\n{ref_cds}\n")

    hits = {}
    for sp in args.compare:
        g = genome_for(args.base, sp)
        if not g:
            print(f"  {sp}: no genome, skipped")
            continue
        rst = None
        for r in args.restrict:
            if ":" in r and r.split(":", 1)[0] == sp:
                rst = r.split(":", 1)[1]
        h = aligned_pairs(qfa, g, restrict=rst, min_cov=args.min_coverage,
                          ref_len=len(ref_cds), task=args.blast_task,
                          word_size=args.word_size)
        if not h:
            print(f"  {sp}: no HSP passing filters"
                  + (f" on {rst}" if rst else "")
                  + f" (task={args.blast_task}, min coverage "
                    f"{args.min_coverage:.0%}), skipped"
                  + ("\n      try --blast-task dc-megablast for divergent species"
                     if args.blast_task == "megablast" else ""))
            continue
        cov = 100.0 * (h["qend"] - h["qstart"] + 1) / len(ref_cds)
        print(f"  {sp}: {h['seqid']} {h['pident']}% over {h['length']} bp "
              f"(query {h['qstart']}-{h['qend']}, {cov:.0f}% of CDS)")
        hits[sp] = h

    if not hits:
        sys.exit("no comparison sequences")

    parts, _cls = partition(ref_cds, ref_coords, hits, ins_lo, ins_hi)

    rows = []
    for sp in hits:
        for cls in ("HOST", "INS"):
            a, b = parts[sp][cls]
            res, ncod = run_kaks("ref", a, sp, b, f"{sp}_{cls}", work, kaks)
            nsyn, nnon = count_diffs(a, b)
            if isinstance(res, str) or res is None:
                rows.append([sp, cls, ncod, nsyn, nnon, "NA", "NA", "NA", "NA",
                             res or "too short"])
            else:
                note = "" if (nsyn >= 5 and nnon >= 5) else "LOW COUNTS - omega unreliable"
                rows.append([sp, cls, ncod, nsyn, nnon, res.get("Ka", "NA"),
                             res.get("Ks", "NA"), res.get("Ka/Ks", "NA"),
                             res.get("P-Value(Fisher)", "NA"), note])

    hdr = ["comparison", "region", "codons", "n_syn", "n_nonsyn",
           "Ka", "Ks", "Ka/Ks", "P", "note"]
    with open(args.out + "_frame_kaks.tsv", "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")

    w = [max(len(str(r[i])) for r in [hdr] + rows) for i in range(len(hdr))]
    print()
    print("  ".join(h.ljust(w[i]) for i, h in enumerate(hdr)))
    for r in rows:
        print("  ".join(str(c).ljust(w[i]) for i, c in enumerate(r)))
    # Internal control: HOST and INS are partitions of ONE gene in ONE
    # pairwise comparison, so synonymous divergence must be comparable.
    # A large Ks disparity means the two partitions are not aligned to the
    # same orthologous locus and omega is meaningless, however good the
    # p-value looks.
    print()
    for sp in hits:
        ks = {}
        for r in rows:
            if r[0] == sp and r[6] not in ("NA", None):
                try:
                    ks[r[1]] = float(r[6])
                except ValueError:
                    pass
        nsyn = {}
        for r in rows:
            if r[0] == sp:
                try:
                    nsyn[r[1]] = int(r[3])
                except (ValueError, TypeError):
                    pass
        if len(ks) == 2 and min(ks.values()) > 0:
            ratio = max(ks.values()) / min(ks.values())
            small = min(nsyn.values()) < 5 if len(nsyn) == 2 else True
            if small:
                verdict = ("counts too small to judge - a Ks ratio built on "
                           "<5 synonymous differences is noise, not evidence "
                           "of misalignment")
            elif ratio >= 2:
                verdict = "SUSPECT - partitions may not be orthologous"
            else:
                verdict = "OK"
            print(f"Ks control {sp}: HOST Ks={ks.get('HOST'):.4g} (n_syn="
                  f"{nsyn.get('HOST','?')})  INS Ks={ks.get('INS'):.4g} (n_syn="
                  f"{nsyn.get('INS','?')})  ratio={ratio:.1f}x\n  {verdict}")
            # raw divergence per codon is more robust than Ks at low counts
            tot = {}
            for r in rows:
                if r[0] == sp:
                    try:
                        tot[r[1]] = (int(r[3]) + int(r[4])) / int(r[2])
                    except (ValueError, TypeError, ZeroDivisionError):
                        pass
            if len(tot) == 2:
                print(f"  raw divergence/codon: HOST={tot['HOST']:.3f} "
                      f"INS={tot['INS']:.3f}  ratio={max(tot.values())/max(1e-9,min(tot.values())):.1f}x")

    print(f"\nWrote {args.out}_frame_kaks.tsv")
    print("HOST = codons from the host gene; INS = codons from the insertion.")
    print("Similar omega in both = the insertion is constrained as part of the protein.")


if __name__ == "__main__":
    sys.exit(main())
