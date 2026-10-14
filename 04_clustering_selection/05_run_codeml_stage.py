#!/usr/bin/env python3
"""
05_run_codeml_stage.py

Labels trees, runs codeml under both models, then parses the results.

Only processes analyses that have BOTH a codon alignment and a subtree,
and skips any that already have parseable output unless --force is given.

Usage:
    python3 05_run_codeml_stage.py [--threads N] [--force]
                                   [--organelle nupt|numt] [selection_dir]

Re-running without --force skips finished codeml runs and only re-parses,
so it also replaces the former standalone parser (06_parse_codeml_results.py).
"""

import os
import re
import sys
import glob
import shutil
import subprocess
import dendropy
from concurrent.futures import ProcessPoolExecutor, as_completed

from scipy import stats

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
CLUSTER_DIR = os.path.join(BASE_DIR, "clustering_nupt", "clusters")
SINGLES_DIR = os.path.join(BASE_DIR, "clustering_nupt", "singles")


CODEML_TIMEOUT = 3600  # seconds per model; guards against non-convergence

def set_organelle(kind):
    """Point the module at plastid or mitochondrial inputs.

    codeml's icode stays 0 (universal) for both: the plant plastid and
    plant mitochondrial codes share the universal code's
    codon-to-amino-acid assignments, so synonymous and nonsynonymous site
    classification is identical.
    """
    global CLUSTER_DIR, SINGLES_DIR
    if kind == "mito":
        CLUSTER_DIR = os.path.join(BASE_DIR, "clustering_numt", "clusters")
        SINGLES_DIR = os.path.join(BASE_DIR, "clustering_numt", "singles")
    else:
        CLUSTER_DIR = os.path.join(BASE_DIR, "clustering_nupt", "clusters")
        SINGLES_DIR = os.path.join(BASE_DIR, "clustering_nupt", "singles")

OMEGA_BOUND = 10.0   # above this, dS is ~0 and omega is uninformative

CTL_TEMPLATE = """      seqfile  = {codeml_dir}/input.phy
      treefile = {codeml_dir}/tree_labeled.nwk
      outfile  = {codeml_dir}/{model}.out
      noisy    = 0
      verbose  = 0
      runmode  = 0
      seqtype  = 1
      CodonFreq = 2
      model    = 2
      NSsites  = 0
      icode    = 0
      fix_kappa = 0
      kappa    = 2
      fix_omega = {fix_omega}
      omega    = 1
      cleandata = 1
      fix_alpha = 1
      alpha    = 0
"""


def nupt_ids_for(name):
    for d in (CLUSTER_DIR, SINGLES_DIR):
        fa = os.path.join(d, f"{name}.fa")
        if os.path.exists(fa):
            ids = []
            with open(fa) as f:
                for line in f:
                    if line.startswith(">"):
                        ids.append(line[1:].strip())
            return ids
    return None


def has_lnl(path):
    if not os.path.exists(path):
        return False
    with open(path) as f:
        return "lnL" in f.read()


def label_tree(tree_path, nupt_ids):
    """Write PAML branch labels following Chen et al.'s numt_dnds.7.py.

    Clusters: every branch strictly BELOW the MRCA of the NUPT sequences is
    foreground (#1). The stem branch leading to the MRCA is deliberately
    left unlabelled, so it falls into the background class.

    That stem is the branch on which nuclear integration took place, and it
    spans time during which the sequence was still organellar. Chen et al.
    describe this as a node-based rather than branch-based last common
    ancestor approach, so that the estimate reflects selection on the
    insertion rather than on the organelle DNA. Their script marks the
    daughter clades of the MRCA and never the MRCA's own branch (version
    note: "Tree marking after unroot, exclude the marking the ancestral
    branch").

    Note their Figure 5B caption instead describes the stem as a separate
    background class with its own omega_1, which would give three ratios
    rather than two. Their codeml control files were not available to
    resolve this. The implementation here follows the script, which is
    unambiguous that the stem is not foreground.

    Singletons: the single foreground branch is #1, matching Figure 5A.

    codeml's fix_omega constrains the last rate class, so the null model
    fixes the foreground in both cases.
    """
    tree = dendropy.Tree.get(path=tree_path, schema="newick",
                             preserve_underscores=True)
    tree.is_rooted = True

    lab = {t.label: t for t in tree.taxon_namespace}
    taxa = [lab[i] for i in nupt_ids if i in lab]
    if not taxa:
        return None

    is_cluster = len(taxa) > 1
    if is_cluster:
        focal = tree.mrca(taxa=taxa)
    else:
        focal = tree.find_node_with_taxon_label(taxa[0].label)
    if focal is None:
        return None

    # Foreground = branches strictly below the MRCA. The MRCA's own branch
    # (the integration branch) is left unlabelled and joins the background.
    foreground = set()
    if is_cluster:
        for n in focal.preorder_iter():
            if n is not focal:
                foreground.add(id(n))
    else:
        foreground.add(id(focal))

    def rec(node):
        if node.is_leaf():
            s = node.taxon.label if node.taxon is not None else ""
        else:
            s = "(" + ",".join(rec(ch) for ch in node.child_nodes()) + ")"
        if id(node) in foreground:
            s += " #1"
        if node.edge.length is not None and node.parent_node is not None:
            s += ":%.9g" % node.edge.length
        return s

    return rec(tree.seed_node) + ";"


def run_one(args):
    name, out_dir, force, organelle = args
    set_organelle(organelle)
    phy = os.path.join(out_dir, "codon_alignments", f"{name}.phy")
    tree = os.path.join(out_dir, "subtrees", f"{name}.nwk")
    codeml_dir = os.path.join(out_dir, "codeml", name)

    if not (os.path.exists(phy) and os.path.exists(tree)):
        return (name, "missing_input")

    if not force and has_lnl(os.path.join(codeml_dir, "modelA.out")) \
            and has_lnl(os.path.join(codeml_dir, "modelB.out")):
        return (name, "already_done")

    ids = nupt_ids_for(name)
    if not ids:
        return (name, "no_nupt_fasta")

    os.makedirs(codeml_dir, exist_ok=True)

    try:
        nwk = label_tree(tree, ids)
    except Exception as e:
        return (name, f"labelling_failed:{type(e).__name__}")
    if nwk is None:
        return (name, "nupt_tips_absent_from_tree")
    with open(os.path.join(codeml_dir, "tree_labeled.nwk"), "w") as f:
        f.write(nwk + "\n")

    shutil.copy(phy, os.path.join(codeml_dir, "input.phy"))

    for model, fix_omega in (("modelA", 0), ("modelB", 1)):
        ctl = os.path.join(codeml_dir, f"{model}.ctl")
        with open(ctl, "w") as f:
            f.write(CTL_TEMPLATE.format(codeml_dir=codeml_dir, model=model,
                                        fix_omega=fix_omega))
        # codeml exits non-zero with "error: end of tree file." even after a
        # successful single-tree run, so exit status is not informative.
        #
        # stdin=DEVNULL matters: codeml prompts interactively when it cannot
        # parse an input, and with stdout suppressed the prompt is invisible,
        # so the process sleeps forever consuming no CPU. Feeding it EOF turns
        # that into an immediate failure instead of an indefinite hang. The
        # timeout is a second guard against non-convergence on large
        # alignments.
        try:
            with open(os.path.join(codeml_dir, f"{model}.log"), "w") as log:
                subprocess.run(["codeml", f"{model}.ctl"], cwd=codeml_dir,
                               stdin=subprocess.DEVNULL,
                               stdout=log, stderr=subprocess.STDOUT,
                               timeout=CODEML_TIMEOUT)
        except subprocess.TimeoutExpired:
            return (name, "codeml_timeout")

    ok_a = has_lnl(os.path.join(codeml_dir, "modelA.out"))
    ok_b = has_lnl(os.path.join(codeml_dir, "modelB.out"))
    if ok_a and ok_b:
        return (name, "ok")
    return (name, "no_lnl")


def parse_lnl(outfile):
    if not os.path.exists(outfile):
        return None
    with open(outfile) as f:
        for line in f:
            m = re.search(r"lnL\s*\([^)]*\)\s*:\s*(-?[\d.]+)", line)
            if m:
                return float(m.group(1))
    return None


def parse_omegas(outfile):
    if not os.path.exists(outfile):
        return []
    with open(outfile) as f:
        content = f.read()
    m = re.search(r"w \(dN/dS\) for branches:\s*(.*)", content)
    if not m:
        return []
    return [float(v) for v in
            re.findall(r"\d+\.?\d*(?:[eE][+\-]?\d+)?", m.group(1))]


def parse_all(out_dir):
    results = []
    n_no_lnl = 0
    for codeml_dir in sorted(glob.glob(os.path.join(out_dir, "codeml", "*"))):
        name = os.path.basename(codeml_dir)
        lnlA = parse_lnl(os.path.join(codeml_dir, "modelA.out"))
        lnlB = parse_lnl(os.path.join(codeml_dir, "modelB.out"))
        if lnlA is None or lnlB is None:
            n_no_lnl += 1
            continue
        omegas = parse_omegas(os.path.join(codeml_dir, "modelA.out"))
        # Two rate classes: [background, foreground]. The stem branch of a
        # cluster is unlabelled and therefore part of the background, so
        # the foreground is the last class in both clusters and singletons -
        # which is also the class fix_omega constrains in the null model.
        bg = omegas[0] if len(omegas) > 0 else None
        stem = None
        fg = omegas[-1] if len(omegas) > 1 else None
        lrt = max(0.0, 2 * (lnlA - lnlB))
        results.append({
            "name": name, "omega_bg": bg, "omega_stem": stem, "omega_fg": fg,
            "uninformative": (fg is not None and fg >= OMEGA_BOUND),
            "lnlA": lnlA, "lnlB": lnlB, "lrt": lrt,
            "pval": stats.chi2.sf(lrt, df=1),
        })

    if not results:
        print("No parseable results.")
        return

    results.sort(key=lambda r: r["pval"])
    n = len(results)
    for i, r in enumerate(results):
        r["pval_bh"] = min(r["pval"] * n / (i + 1), 1.0)

    def verdict(r):
        if r["omega_fg"] is None:
            return "no_omega"
        if r["uninformative"]:
            return "omega_uninformative"
        if r["pval_bh"] >= 0.05:
            return "not_significant"
        return "purifying" if r["omega_fg"] < 1 else "positive"

    res_dir = os.path.join(out_dir, "results")
    os.makedirs(res_dir, exist_ok=True)
    tsv = os.path.join(res_dir, "selection_results.tsv")
    with open(tsv, "w") as f:
        f.write("name\ttype\tomega_background\tomega_stem\tomega_nupt\t"
                "omega_uninformative\t"
                "lnlA\tlnlB\tLRT\tpval\tpval_BH\tverdict\n")
        for r in results:
            typ = "cluster" if "_cluster_" in r["name"] else "singleton"
            bg = f"{r['omega_bg']:.4f}" if r["omega_bg"] is not None else "NA"
            st = f"{r['omega_stem']:.4f}" if r["omega_stem"] is not None else "NA"
            fg = f"{r['omega_fg']:.4f}" if r["omega_fg"] is not None else "NA"
            f.write(f"{r['name']}\t{typ}\t{bg}\t{st}\t{fg}\t{r['uninformative']}\t"
                    f"{r['lnlA']:.4f}\t{r['lnlB']:.4f}\t{r['lrt']:.4f}\t"
                    f"{r['pval']:.6g}\t{r['pval_bh']:.6g}\t{verdict(r)}\n")

    verdicts = [verdict(r) for r in results]
    n_clu = sum(1 for r in results if "_cluster_" in r["name"])

    def tally(v, kind):
        return sum(1 for r in results
                   if verdict(r) == v
                   and (("_cluster_" in r["name"]) == (kind == "cluster")))

    print(f"\nResults: {tsv}")
    print(f"Parsed: {n}   (clusters {n_clu}, singletons {n - n_clu})")
    print(f"Discarded (missing lnL): {n_no_lnl}")
    print("\n                        all    clusters  singletons")
    for v in ["purifying", "positive", "not_significant",
              "omega_uninformative", "no_omega"]:
        print(f"  {v:<21} {verdicts.count(v):>5} "
              f"{tally(v, 'cluster'):>9} {tally(v, 'singleton'):>11}")

    good = sorted(r["omega_fg"] for r in results
                  if r["omega_fg"] is not None and not r["uninformative"])
    if good:
        print(f"\nForeground omega (informative, n={len(good)}): "
              f"min={good[0]:.4f} median={good[len(good)//2]:.4f} "
              f"max={good[-1]:.4f}")


def main():
    argv = [a for a in sys.argv[1:]]
    force = "--force" in argv
    if force:
        argv.remove("--force")
    threads = 5
    if "--threads" in argv:
        i = argv.index("--threads")
        threads = int(argv[i + 1])
        del argv[i:i + 2]
    organelle = "plastid"
    if "--organelle" in argv:
        i = argv.index("--organelle")
        organelle = normalise_organelle(argv[i + 1])
        del argv[i:i + 2]
    set_organelle(organelle)
    out_dir = argv[0] if argv else default_selection_dir(organelle)
    print(f"organelle={organelle}")

    names = sorted(
        os.path.basename(p)[:-4]
        for p in glob.glob(os.path.join(out_dir, "codon_alignments", "*.phy"))
    )
    print(f"Codon alignments found: {len(names)}")
    print(f"Subtrees found: "
          f"{len(glob.glob(os.path.join(out_dir, 'subtrees', '*.nwk')))}")
    print(f"Running codeml with {threads} workers "
          f"({'forcing rerun' if force else 'skipping completed'})...\n")

    counts = {}
    done = 0
    with ProcessPoolExecutor(max_workers=threads) as ex:
        futures = [ex.submit(run_one, (nm, out_dir, force, organelle))
                   for nm in names]
        for fut in as_completed(futures):
            nm, status = fut.result()
            counts[status] = counts.get(status, 0) + 1
            done += 1
            if done % 50 == 0:
                print(f"  {done}/{len(names)} processed")

    print("\ncodeml stage summary:")
    for k in sorted(counts):
        print(f"  {k}: {counts[k]}")

    parse_all(out_dir)


if __name__ == "__main__":
    main()
