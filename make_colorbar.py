#!/usr/bin/env python3
"""
make_colorbar.py  —  Generate a standalone colourbar PDF to go with the fan plots.

Usage:
    python make_colorbar.py \
        --vmin   -5.4 \
        --vmax   -0.3 \
        --label  "Proportion of genome (log₁₀ scale)" \
        --out    colorbar.pdf
"""

import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm

CMAP = mcolors.LinearSegmentedColormap.from_list(
    "numt",
    [(0.00, "#1a006e"),
     (0.25, "#0088cc"),
     (0.55, "#00cc88"),
     (0.80, "#88dd00"),
     (1.00, "#ffee00")],
)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--vmin",  type=float, required=True,
                        help="log10 minimum (e.g. -5.4)")
    parser.add_argument("--vmax",  type=float, required=True,
                        help="log10 maximum (e.g. -0.3)")
    parser.add_argument("--label", default="Proportion of genome (log₁₀ scale)")
    parser.add_argument("--out",   default="colorbar.pdf")
    parser.add_argument("--width", type=float, default=6,
                        help="Figure width in inches")
    parser.add_argument("--height", type=float, default=0.8,
                        help="Figure height in inches")
    args = parser.parse_args()

    cnorm = mcolors.Normalize(vmin=args.vmin, vmax=args.vmax)

    fig, ax = plt.subplots(figsize=(args.width, args.height),
                           facecolor="white")

    sm = cm.ScalarMappable(cmap=CMAP, norm=cnorm)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=ax, orientation="horizontal")
    cb.set_label(args.label, fontsize=11)
    cb.ax.tick_params(labelsize=10)

    mid_val = (args.vmin + args.vmax) / 2
    cb.set_ticks([args.vmin, mid_val, args.vmax])
    cb.set_ticklabels([
        f"10^{args.vmin:.1f}",
        f"10^{mid_val:.1f}",
        f"10^{args.vmax:.1f}",
    ])

    plt.savefig(args.out, dpi=300, bbox_inches="tight")
    print(f"Colourbar saved to: {args.out}")

if __name__ == "__main__":
    main()
