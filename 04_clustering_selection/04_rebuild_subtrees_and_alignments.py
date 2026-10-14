#!/usr/bin/env python3
"""
04_rebuild_subtrees_and_alignments.py

Usage:
    python3 04_rebuild_subtrees_and_alignments.py [--levels 3] [--max-taxa 50]
                                               [--threads 5] [--organelle nupt|numt]
                                               [selection_dir]
"""

import os
import re
import sys
import glob
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed

import dendropy

BASE_DIR = os.environ.get("BASE_DIR", os.getcwd())


def normalise_organelle(kind):
    """Accept nupt/plastid and numt/mito; return 'plastid' or 'mito'."""
    k = kind.lower()
    if k in ("nupt", "plastid"):
        return "plastid"
    if k in ("numt", "mito"):
        return "mito"
    sys.exit(f"--organelle must be nupt or numt (got {kind!r})")


def default_selection_dir(organelle):
    return os.path.join(BASE_DIR, "selection_numt" if organelle == "mito"
                        else "selection_nupt")

# Set by --organelle (nupt default, or numt). pal2nal -codontable uses NCBI
# numbering: 11 is the plant plastid code, 1 the standard code, which is
# correct for plant mitochondria. This is a different scheme from codeml's
# icode (set to 0 elsewhere).
CLUSTER_DIR = os.path.join(BASE_DIR, "clustering_nupt", "clusters")
SINGLES_DIR = os.path.join(BASE_DIR, "clustering_nupt", "singles")
NCBI_CODON_TABLE = 11
BG_ALL_NAME = "all_plastid.fasta"

MIN_TAXA = 3


def set_organelle(kind):
    """Point the module at plastid or mitochondrial inputs."""
    global CLUSTER_DIR, SINGLES_DIR, NCBI_CODON_TABLE, BG_ALL_NAME
    if kind == "mito":
        CLUSTER_DIR = os.path.join(BASE_DIR, "clustering_numt", "clusters")
        SINGLES_DIR = os.path.join(BASE_DIR, "clustering_numt", "singles")
        NCBI_CODON_TABLE = 1
        BG_ALL_NAME = "all_mito.fasta"
    else:
        CLUSTER_DIR = os.path.join(BASE_DIR, "clustering_nupt", "clusters")
        SINGLES_DIR = os.path.join(BASE_DIR, "clustering_nupt", "singles")
        NCBI_CODON_TABLE = 11
        BG_ALL_NAME = "all_plastid.fasta"


def nupt_ids_for(name):
    for d in (CLUSTER_DIR, SINGLES_DIR):
        fa = os.path.join(d, f"{name}.fa")
        if os.path.exists(fa):
            ids = []
            with open(fa) as f:
                for line in f:
                    if line.startswith(">"):
                        ids.append(line[1:].strip())
            return ids, fa
    return None, None


def gene_for(name):
    if "_cluster_" in name:
        return re.sub(r"_cluster_.*", "", name)
    return re.sub(r"_REV$", "", name).rsplit("_", 1)[-1]


def resolve_bg(out_dir, gene, suffix):
    """Background files keep the plastid FASTA capitalisation (atpA, rbcL),
    while cluster names are lowercased, so match case-insensitively."""
    bg_dir = os.path.join(out_dir, "background")
    exact = os.path.join(bg_dir, f"{gene}_{suffix}.fa")
    if os.path.exists(exact):
        return exact
    target = f"{gene}_{suffix}.fa".lower()
    for f in os.listdir(bg_dir):
        if f.lower() == target:
            return os.path.join(bg_dir, f)
    return None


def read_fasta(path):
    recs, name, buf = {}, None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    recs[name] = "".join(buf)
                name = line[1:].strip().split()[0]
                buf = []
            elif name is not None:
                buf.append(line.strip())
    if name is not None:
        recs[name] = "".join(buf)
    return recs


def build_subtree(tree_path, nupt_ids, levels, max_taxa):
    """Resolve polytomies first, then climb. Returns (subtree, n_taxa, note)."""
    tree = dendropy.Tree.get(path=tree_path, schema="newick",
                             preserve_underscores=True)
    tree.is_rooted = True

    # Resolve BEFORE locating the MRCA so that climbing traverses
    # bifurcations rather than jumping across a wide polytomy.
    tree.resolve_polytomies(limit=2, update_bipartitions=True)
    for e in tree.postorder_edge_iter():
        if e.length is None:
            e.length = 0.0

    label_to_taxon = {t.label: t for t in tree.taxon_namespace}
    taxa = [label_to_taxon[i] for i in nupt_ids if i in label_to_taxon]
    if not taxa:
        return None, 0, "nupt_tips_absent_from_tree"

    if len(taxa) == 1:
        node = tree.find_node_with_taxon_label(taxa[0].label)
    else:
        node = tree.mrca(taxa=taxa)
    if node is None:
        return None, 0, "no_mrca"

    # One climb is mandatory: a bare tip has no background to compare
    # against. Further climbs are taken only while they stay within
    # max_taxa. If even the immediate parent exceeds the cap, the analysis
    # is not tractable and is reported rather than silently truncated.
    note = "ok"
    if node.parent_node is None:
        return None, 0, "nupt_at_root"
    node = node.parent_node
    n_now = sum(1 for lf in node.leaf_iter() if lf.taxon is not None)
    if n_now > max_taxa:
        return None, n_now, "min_clade_exceeds_cap"

    for _ in range(levels - 1):
        parent = node.parent_node
        if parent is None:
            break
        n_up = sum(1 for lf in parent.leaf_iter() if lf.taxon is not None)
        if n_up > max_taxa:
            note = "capped_at_max_taxa"
            break
        node = parent

    clade_taxa = [lf.taxon for lf in node.leaf_iter() if lf.taxon is not None]
    if len(clade_taxa) < MIN_TAXA:
        return None, len(clade_taxa), "too_few_taxa"

    subtree = tree.extract_tree_with_taxa(taxa=clade_taxa)
    subtree.resolve_polytomies(limit=2, update_bipartitions=True)
    for e in subtree.postorder_edge_iter():
        if e.length is None:
            e.length = 0.0
    return subtree, len(clade_taxa), note


def process_one(args):
    name, out_dir, levels, max_taxa, organelle = args
    set_organelle(organelle)

    tree_path = os.path.join(out_dir, "trees", f"{name}.nwk")
    if not os.path.exists(tree_path):
        return name, "no_tree", 0

    nupt_ids, nupt_fa = nupt_ids_for(name)
    if not nupt_ids:
        return name, "no_nupt_fasta", 0

    try:
        subtree, n_taxa, note = build_subtree(tree_path, nupt_ids,
                                              levels, max_taxa)
    except Exception as e:
        return name, f"tree_error:{type(e).__name__}", 0
    if subtree is None:
        return name, note, n_taxa

    subtree_out = os.path.join(out_dir, "subtrees", f"{name}.nwk")
    subtree.write(path=subtree_out, schema="newick",
                  suppress_rooting=True, unquoted_underscores=True)

    # NOTE: subtree.taxon_namespace is the shared namespace object and
    # still lists every taxon from the parent tree, so using it here left
    # the alignment unpruned (539 sequences against a 16-taxon tree, which
    # codeml rejects). Iterate the subtree's actual leaves instead.
    keep = {lf.taxon.label for lf in subtree.leaf_node_iter()
            if lf.taxon is not None}

    aln_dir = os.path.join(out_dir, "prank_alignments", name)
    prot_path = os.path.join(aln_dir, "combined_aln.fa")
    if not os.path.exists(prot_path):
        return name, "no_prot_alignment", n_taxa

    gene = gene_for(name)
    bg_nuc = resolve_bg(out_dir, gene, "nuc_filt")
    if bg_nuc is None:
        return name, "no_background_nuc", n_taxa

    prot = read_fasta(prot_path)
    nuc = read_fasta(bg_nuc)
    nuc.update(read_fasta(nupt_fa))

    common = [k for k in prot if k in keep and k in nuc]
    if len(common) < MIN_TAXA:
        return name, "too_few_after_pruning", len(common)

    pruned_prot = os.path.join(aln_dir, "pruned_prot.fa")
    pruned_nuc = os.path.join(aln_dir, "pruned_nuc.fa")
    with open(pruned_prot, "w") as f:
        for k in common:
            f.write(f">{k}\n{prot[k]}\n")
    with open(pruned_nuc, "w") as f:
        for k in common:
            f.write(f">{k}\n{nuc[k]}\n")

    phy_out = os.path.join(out_dir, "codon_alignments", f"{name}.phy")
    # -nogap deliberately omitted: it drops every column gapped in any
    # sequence, which emptied alignments entirely. codeml's cleandata
    # handles gaps instead.
    with open(phy_out, "w") as fh:
        r = subprocess.run(
            ["pal2nal.pl", pruned_prot, pruned_nuc,
             "-codontable", str(NCBI_CODON_TABLE), "-output", "paml"],
            stdout=fh, stderr=subprocess.DEVNULL)
    if r.returncode != 0 or os.path.getsize(phy_out) == 0:
        if os.path.exists(phy_out):
            os.remove(phy_out)
        return name, "pal2nal_failed", n_taxa

    with open(phy_out) as f:
        header = f.readline().split()
    if len(header) < 2 or int(header[1]) == 0:
        os.remove(phy_out)
        return name, "zero_sites", n_taxa

    return name, "ok", len(common)


def main():
    argv = list(sys.argv[1:])

    def opt(flag, default, cast=int):
        if flag in argv:
            i = argv.index(flag)
            v = cast(argv[i + 1])
            del argv[i:i + 2]
            return v
        return default

    levels = opt("--levels", 3)
    max_taxa = opt("--max-taxa", 50)
    threads = opt("--threads", 5)
    organelle = normalise_organelle(opt("--organelle", "nupt", str))
    set_organelle(organelle)
    out_dir = argv[0] if argv else default_selection_dir(organelle)
    print(f"organelle={organelle}  codontable={NCBI_CODON_TABLE}")

    for sub in ("subtrees", "codon_alignments"):
        os.makedirs(os.path.join(out_dir, sub), exist_ok=True)

    names = sorted(os.path.basename(p)[:-4]
                   for p in glob.glob(os.path.join(out_dir, "trees", "*.nwk")))
    print(f"Trees found: {len(names)}")
    print(f"levels={levels}  max_taxa={max_taxa}  threads={threads}\n")

    counts, sizes, done = {}, [], 0
    with ProcessPoolExecutor(max_workers=threads) as ex:
        futures = [ex.submit(process_one,
                             (nm, out_dir, levels, max_taxa, organelle))
                   for nm in names]
        for fut in as_completed(futures):
            name, status, n = fut.result()
            counts[status] = counts.get(status, 0) + 1
            if status == "ok":
                sizes.append(n)
            done += 1
            if done % 100 == 0:
                print(f"  {done}/{len(names)} processed")

    print("\nSummary:")
    for k in sorted(counts, key=lambda x: -counts[x]):
        print(f"  {k}: {counts[k]}")

    if sizes:
        s = sorted(sizes)
        print(f"\nSubtree size (taxa) for the {len(s)} usable analyses:")
        print(f"  min={s[0]}  median={s[len(s)//2]}  "
              f"90th={s[int(len(s)*0.9)]}  max={s[-1]}")
        print("\nNext: python3 05_run_codeml_stage.py --threads 5 --organelle "
              + ("numt" if organelle == "mito" else "nupt"))


if __name__ == "__main__":
    main()
