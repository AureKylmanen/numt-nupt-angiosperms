#!/usr/bin/env python3
"""
03_build_trees_iqtree.py

ML protein trees (IQ-TREE 2, LG, -fast) for every combined background +
dcNUxT alignment. For NUMTs pass the selection_numt directory.

Usage:
    python3 03_build_trees_iqtree.py [--jobs 7] [--threads 4]
                                  [--model LG] [--force] [selection_dir]
"""

import os
import sys
import glob
import shutil
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed

BASE_DIR = "/vol/storage/5T_disk/organelles"


def build_one(args):
    name, out_dir, threads, model, force = args

    aln = os.path.join(out_dir, "prank_alignments", name, "combined_aln.fa")
    tree_out = os.path.join(out_dir, "trees", f"{name}.nwk")

    if not os.path.exists(aln) or os.path.getsize(aln) == 0:
        return name, "no_alignment"

    if not force and os.path.exists(tree_out) and os.path.getsize(tree_out) > 0:
        return name, "already_done"

    # IQ-TREE writes several files per prefix; keep them out of the way
    with tempfile.TemporaryDirectory() as tmp:
        prefix = os.path.join(tmp, "t")
        r = subprocess.run(
            ["iqtree2", "-s", aln, "-m", model, "-fast",
             "-nt", str(threads), "-pre", prefix, "-quiet"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        treefile = prefix + ".treefile"
        if r.returncode != 0 or not os.path.exists(treefile):
            return name, "iqtree_failed"
        shutil.copy(treefile, tree_out)

    return name, "ok"


def main():
    argv = list(sys.argv[1:])
    force = "--force" in argv
    if force:
        argv.remove("--force")

    def opt(flag, default, cast=int):
        if flag in argv:
            i = argv.index(flag)
            v = cast(argv[i + 1])
            del argv[i:i + 2]
            return v
        return default

    jobs = opt("--jobs", 7)
    threads = opt("--threads", 4)
    model = opt("--model", "LG", str)
    out_dir = argv[0] if argv else os.path.join(BASE_DIR, "selection")

    os.makedirs(os.path.join(out_dir, "trees"), exist_ok=True)

    names = sorted(
        os.path.basename(os.path.dirname(p))
        for p in glob.glob(os.path.join(out_dir, "prank_alignments",
                                        "*", "combined_aln.fa"))
    )
    print(f"Alignments found: {len(names)}")
    print(f"jobs={jobs} threads={threads} model={model} "
          f"({'forcing rerun' if force else 'skipping completed'})")
    print(f"Using {jobs * threads} cores total\n")

    counts, done = {}, 0
    with ProcessPoolExecutor(max_workers=jobs) as ex:
        futures = [ex.submit(build_one, (nm, out_dir, threads, model, force))
                   for nm in names]
        for fut in as_completed(futures):
            nm, status = fut.result()
            counts[status] = counts.get(status, 0) + 1
            done += 1
            if done % 25 == 0:
                print(f"  {done}/{len(names)} trees")

    print("\nSummary:")
    for k in sorted(counts, key=lambda x: -counts[x]):
        print(f"  {k}: {counts[k]}")

    # Report the worst multifurcation remaining, as a check that the
    # trees really are bifurcating now.
    try:
        import dendropy
        worst = 0
        sample = sorted(glob.glob(os.path.join(out_dir, "trees", "*.nwk")))[:50]
        for f in sample:
            t = dendropy.Tree.get(path=f, schema="newick",
                                  preserve_underscores=True)
            for n in t.preorder_node_iter():
                worst = max(worst, len(n.child_nodes()))
        print(f"\nMax children at any node across {len(sample)} sampled "
              f"trees: {worst}")
        print("(3 at the root of an unrooted tree is normal; more than that "
              "means multifurcations remain)")
    except ImportError:
        pass

    print("\nNext: python3 04_rebuild_subtrees_and_alignments.py "
          "--levels 3 --max-taxa 50 --threads 5 [--organelle mito]")


if __name__ == "__main__":
    main()
