#!/usr/bin/env python3
"""
plot_proportion_fan.py  —  Fan phylogenetic tree coloured by log10(proportion)
Dependencies:  pip install dendropy matplotlib numpy pandas

Usage:
    python plot_proportion_fan.py \
        --tree        tree.nwk \
        --counts      numt_proportion.tsv \
        --taxonomy    species_taxonomy.tsv \
        --label       "NUMT proportion of genome" \
        --rank        order \
        --min-clade   5 \
        --out         numt_fan.pdf

    # Run a second time for NUPT:
    python plot_proportion_fan.py \
        --tree        tree.nwk \
        --counts      nupt_proportion.tsv \
        --label       "NUPT proportion of genome" \
        --out         nupt_fan.pdf

    # Force both plots onto the same colour scale for direct comparison:
    python plot_proportion_fan.py --counts numt_proportion.tsv --vmin -6 --vmax -1 --out numt_fan.pdf
    python plot_proportion_fan.py --counts nupt_proportion.tsv --vmin -6 --vmax -1 --out nupt_fan.pdf
"""

import argparse, math, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import dendropy
import re

# ── colour map: deep blue → cyan → teal-green → yellow ───────────────────────
CMAP = mcolors.LinearSegmentedColormap.from_list(
    "numt",
    [(0.00, "#1a006e"),
     (0.25, "#0088cc"),
     (0.55, "#00cc88"),
     (0.80, "#88dd00"),
     (1.00, "#ffee00")],
)

# pseudo-count added before log10 to handle zeros (= effectively -9 on log10 scale)
PSEUDOCOUNT = 1e-9

# ─────────────────────────────────────────────────────────────────────────────
# Name normalisation  (spaces ↔ underscores, strip quotes)
# ─────────────────────────────────────────────────────────────────────────────

SYNONYMS = {
    "citrus_x_aurantiifolia": "citrus_aurantiifolia",
    "citrus_x_aurantium":     "citrus_aurantium",
    "pleuropterus_multiflorus": "fallopia_multiflora",
    "argentina_anserina":     "potentilla_anserina",
    "eriobotrya_japonica":    "rhaphiolepis_bibas",
}


def norm(s):
    s = str(s).strip().strip("'\"")
    s = s.replace("\xa0", " ")
    s = s.replace(" ", "_")
    s = s.replace("_x_", "_")
    s = re.sub(r"_+", "_", s)
    s = s.lower()
    s = SYNONYMS.get(s, s)
    return s


# ─────────────────────────────────────────────────────────────────────────────
# Tree helpers
# ─────────────────────────────────────────────────────────────────────────────

def get_leaves_in_order(tree):
    leaves = []
    def _walk(node):
        if node.is_leaf():
            leaves.append(node)
        else:
            for ch in node.child_nodes():
                _walk(ch)
    _walk(tree.seed_node)
    return leaves


def assign_leaf_angles(leaves):
    n = len(leaves)
    return {norm(leaf.taxon.label): 2 * math.pi * i / n
            for i, leaf in enumerate(leaves)}


def node_mean_angle(node, leaf_angles):
    labels = [norm(l.taxon.label) for l in node.leaf_iter()]
    vals = [leaf_angles[l] for l in labels if l in leaf_angles]
    if not vals:
        return 0.0
    return math.atan2(np.mean(np.sin(vals)), np.mean(np.cos(vals)))


def assign_radii(tree):
    use_bl = any(
        e.length and e.length > 0
        for e in tree.postorder_edge_iter()
        if e.length is not None
    )
    depth = {}
    def _walk(node, d):
        depth[node] = d
        for ch in node.child_nodes():
            inc = (ch.edge.length or 0) if use_bl else 1
            _walk(ch, d + inc)
    _walk(tree.seed_node, 0.0)
    max_d = max(depth.values()) or 1.0
    return {node: v / max_d for node, v in depth.items()}


def node_mean_value(node, value_map, vmin):
    """
    Average the log10(proportion) values across all leaves of this node.
    Missing leaves fall back to vmin (i.e. the lowest observed value)
    rather than 0, so they don't pull the colour toward the middle of
    the scale.
    """
    vals = [value_map.get(norm(l.taxon.label), vmin) for l in node.leaf_iter()]
    return float(np.mean(vals)) if vals else vmin


# ─────────────────────────────────────────────────────────────────────────────
# Clade arc helpers
# ─────────────────────────────────────────────────────────────────────────────

def get_clade_arcs(leaves, leaf_angles, taxonomy_df, rank, min_size):
    sp2rank  = dict(zip(taxonomy_df["species"].map(norm), taxonomy_df[rank]))
    size_map = taxonomy_df[rank].value_counts().to_dict()

    runs = []
    for leaf in leaves:
        name  = norm(leaf.taxon.label)
        label = sp2rank.get(name, "Unknown")
        angle = leaf_angles[name]
        if runs and runs[-1][0] == label:
            runs[-1][1].append(angle)
        else:
            runs.append((label, [angle]))

    arcs = []
    leaf_step = 2 * math.pi / len(leaves)

    for label, angles in runs:
        if label == "Unknown":
            continue
        if size_map.get(label, 0) < min_size:
            continue
        a0 = min(angles) - leaf_step * 0.5
        a1 = max(angles) + leaf_step * 0.5
        arcs.append(dict(label=label, start=a0, end=a1))

    return arcs


# ─────────────────────────────────────────────────────────────────────────────
# Drawing
# ─────────────────────────────────────────────────────────────────────────────

def draw_fan(tree, value_map, taxonomy_df, rank, min_clade,
             out_path, dpi, figsize, cbar_label, vmin_override, vmax_override,
             no_cbar=False):

    leaves      = get_leaves_in_order(tree)
    leaf_angles = assign_leaf_angles(leaves)
    radii       = assign_radii(tree)

    # ── colour scale: vmin/vmax on the log10 axis ─────────────────────────────
    all_vals = list(value_map.values())
    vmin = vmin_override if vmin_override is not None else min(all_vals)
    vmax = vmax_override if vmax_override is not None else max(all_vals)
    cnorm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    fig = plt.figure(figsize=figsize, facecolor="white")
    # when no colourbar, the tree can use the full vertical space
    ax_rect = [0.02, 0.02, 0.72, 0.96] if not no_cbar else [0.02, 0.02, 0.96, 0.96]
    ax  = fig.add_axes(ax_rect)
    ax.set_aspect("equal")
    ax.axis("off")

    # ── branches ──────────────────────────────────────────────────────────────
    for node in tree.preorder_node_iter():
        if node is tree.seed_node:
            continue
        parent = node.parent_node
        r_p    = radii[parent]
        r_n    = radii[node]
        a_n    = node_mean_angle(node,   leaf_angles)
        a_p    = node_mean_angle(parent, leaf_angles)
        color  = CMAP(cnorm(node_mean_value(node, value_map, vmin)))
        lw     = 1.5

        # radial segment
        ax.plot([r_p * math.cos(a_n), r_n * math.cos(a_n)],
                [r_p * math.sin(a_n), r_n * math.sin(a_n)],
                color=color, lw=lw, solid_capstyle="round")

        # connecting arc at parent radius
        if abs(a_p - a_n) > 1e-6 and r_p > 0:
            a_lo, a_hi = min(a_p, a_n), max(a_p, a_n)
            if a_hi - a_lo > math.pi:
                a_lo, a_hi = a_hi, a_lo + 2 * math.pi
            arc_a = np.linspace(a_lo, a_hi,
                                max(30, int(abs(a_hi - a_lo) * 60)))
            ax.plot(r_p * np.cos(arc_a), r_p * np.sin(arc_a),
                    color=color, lw=lw, solid_capstyle="round")

    # ── clade arc labels ───────────────────────────────────────────────────────
    if taxonomy_df is not None:
        ARC_R   = 1.08
        TICK_IN = 1.04
        LABEL_R = 1.15
        arcs = get_clade_arcs(leaves, leaf_angles, taxonomy_df, rank, min_clade)
        print(f"  Drawing {len(arcs)} clade arcs for rank '{rank}'")

        for arc in arcs:
            a0, a1 = arc["start"], arc["end"]
            if a1 < a0:
                a0, a1 = a1, a0
            span = min(a1 - a0, 2 * math.pi)
            a1   = a0 + span
            mid  = (a0 + a1) / 2

            for a_t in [a0, a1]:
                ax.plot([TICK_IN * math.cos(a_t), ARC_R * math.cos(a_t)],
                        [TICK_IN * math.sin(a_t), ARC_R * math.sin(a_t)],
                        color="black", lw=2, zorder=5)

            arc_a = np.linspace(a0, a1, max(30, int(span * 80)))
            ax.plot(ARC_R * np.cos(arc_a), ARC_R * np.sin(arc_a),
                    color="black", lw=2, solid_capstyle="butt", zorder=5)

            deg = math.degrees(mid) % 360
            ha  = "left" if deg < 90 or deg > 270 else "right"
            rot = deg if ha == "left" else deg + 180
            ax.text(LABEL_R * math.cos(mid), LABEL_R * math.sin(mid),
                    arc["label"], ha=ha, va="center",
                    rotation=rot, rotation_mode="anchor",
                    fontsize=10, fontweight="bold", zorder=6)

    # ── colour bar ────────────────────────────────────────────────────────────

    plt.savefig(out_path, dpi=dpi, bbox_inches="tight")
    print(f"Plot saved to: {out_path}")


# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tree",      required=True)
    parser.add_argument("--counts",    required=True,
                        help="TSV with columns: species, <proportion_column>")
    parser.add_argument("--taxonomy",  default=None)
    parser.add_argument("--rank",      default="order")
    parser.add_argument("--min-clade", type=int, default=5)
    parser.add_argument("--out",       default="proportion_fan.pdf")
    parser.add_argument("--dpi",       type=int, default=300)
    parser.add_argument("--figsize",   default="14,14")
    parser.add_argument("--label",     default="Proportion of genome",
                        help="Colour bar label (e.g. 'NUMT proportion of genome')")
    parser.add_argument("--vmin",      type=float, default=None,
                        help="Fix log10 colour scale minimum (e.g. -6). "
                             "Useful for comparing two plots.")
    parser.add_argument("--vmax",      type=float, default=None,
                        help="Fix log10 colour scale maximum (e.g. -1). "
                             "Useful for comparing two plots.")
    parser.add_argument("--no-cbar",   action="store_true",
                        help="Suppress the colour bar (for multi-panel figures "
                             "where only one panel needs it).")
    args = parser.parse_args()

    # ── tree ──────────────────────────────────────────────────────────────────
    print("Loading tree …")
    tree     = dendropy.Tree.get(path=args.tree, schema="newick")
    n_leaves = sum(1 for _ in tree.leaf_node_iter())
    print(f"  {n_leaves} leaves")

    # ── proportions ───────────────────────────────────────────────────────────
    print("Loading proportion values …")
    df = pd.read_csv(args.counts, sep="\t")
    df.columns = [c.strip() for c in df.columns]

    value_col = [c for c in df.columns if c.lower() != "species"][0]
    print(f"  Value column: '{value_col}'")

    df["_sp"] = df["species"].map(norm)

    # warn if values look like counts rather than proportions
    raw_max = df[value_col].astype(float).max()
    if raw_max > 1.0:
        print(f"  WARNING: max value is {raw_max:.1f} (> 1.0). "
              "Are these already proportions? Proceeding anyway.")

    # log10 transform with pseudocount to handle zeros
    df["log10_prop"] = np.log10(df[value_col].astype(float) + PSEUDOCOUNT)
    value_map = dict(zip(df["_sp"], df["log10_prop"]))

    # report the actual data range so the user can set --vmin/--vmax sensibly
    data_min = df["log10_prop"].min()
    data_max = df["log10_prop"].max()
    print(f"  log10 proportion range: {data_min:.2f} to {data_max:.2f}  "
          f"(i.e. 10^{data_min:.2f} to 10^{data_max:.2f})")

    # check tree / count overlap
    tree_names  = {norm(n.taxon.label) for n in tree.leaf_node_iter()}
    count_names = set(value_map.keys())
    matched     = tree_names & count_names
    missing     = tree_names - count_names
    print(f"  Matched {len(matched)}/{n_leaves} tips to proportion data")
    if missing:
        print(f"  WARNING: {len(missing)} tips have no data → coloured as vmin")
        print(f"  First 5 unmatched tree names  : {sorted(missing)[:5]}")
        print(f"  First 5 proportion file names : {sorted(count_names)[:5]}")

    # ── taxonomy ──────────────────────────────────────────────────────────────
    taxonomy_df = None
    if args.taxonomy:
        print("Loading taxonomy …")
        taxonomy_df = pd.read_csv(args.taxonomy, sep="\t")
        taxonomy_df["species"] = taxonomy_df["species"].map(norm)
        tax_names   = set(taxonomy_df["species"])
        missing_tax = tree_names - tax_names
        print(f"  Taxonomy overlap: {len(tree_names & tax_names)}/{n_leaves}")
        if missing_tax:
            print(f"  WARNING: {len(missing_tax)} tree tips missing taxonomy")
            print(f"  First 10: {sorted(missing_tax)[:10]}")
        if args.rank not in taxonomy_df.columns:
            sys.exit(f"ERROR: rank '{args.rank}' not in taxonomy file. "
                     f"Columns: {list(taxonomy_df.columns)}")
        print(f"  {taxonomy_df[args.rank].nunique()} unique values for rank '{args.rank}'")

    fw, fh = [float(x) for x in args.figsize.split(",")]

    # ── shared colour scale for NUMT and NUPT comparison ─────────────────────
    # Set these to the same values in both scripts so the plots are comparable.
    # Use the printed "log10 proportion range" from your first run to pick these.
    FIXED_VMIN = -5.4   # ← change to your desired minimum
    FIXED_VMAX = -0.3   # ← change to your desired maximum

    # args.vmin/args.vmax (command-line) take priority if provided;
    # otherwise the hardcoded values above are used
    vmin_to_use = args.vmin if args.vmin is not None else FIXED_VMIN
    vmax_to_use = args.vmax if args.vmax is not None else FIXED_VMAX

    draw_fan(tree, value_map, taxonomy_df, args.rank, args.min_clade,
             args.out, args.dpi, (fw, fh), args.label, vmin_to_use, vmax_to_use,
             no_cbar=args.no_cbar)

if __name__ == "__main__":
    main()
