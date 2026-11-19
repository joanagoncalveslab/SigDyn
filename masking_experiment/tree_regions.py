#!/usr/bin/env python3
"""
Fig 4a (left) over the missingness levels: cells ordered by tree traversal, coloured by region.
Same as Case_studies/Figure4_LSCC/Fig4a_left_trees.ipynb (load tree, ladderize, iter_leaves,
one bar per cell coloured by its region prefix), run on the ROOTED tree of every level.

Reads   <out>/<name>/pXX/results/trees/cell_tree.rooted        (made by root_trees.py)
Writes  <out>/<name>/analysis/trees/
          distribution_by_region_pXX.png/pdf     one per level (as the notebook)
          distribution_by_region_all_levels.png/pdf   all levels stacked, p00 on top
          region_legend.pdf
          leaf_order_pXX.csv                     leaf order and region per level
          tree_regions_summary.csv               level, cells, root cell (outgroup)

  python tree_regions.py --name neck
  python tree_regions.py --name neck --levels 0 90
  python tree_regions.py --name neck --no-ladderize
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ete3 import Tree
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent

# from Fig4a_left_trees.ipynb
region_colors = {
    "neck15": "#E69F00",  # orange
    "neck17": "#56B4E9",  # sky blue
    "neck18": "#D55E00",  # reddish
    "neck16": "#CC79A7",  # purple/magenta
}
region_mapping = {
    "neck15": "Normal laryngeal mucosal epithelia (N)",
    "neck17": "Margins of cancer (R)",
    "neck18": "Carcinoma in situ (T)",
    "neck16": "Lymph nodes metastasis (L)",
}
legend_order = ["neck15", "neck17", "neck18", "neck16"]  # N, R, T, L


def leaf_order(tree_file, ladderize):
    t = Tree(str(tree_file), format=1)
    if ladderize:
        t.ladderize()   # before iter_leaves(); applied to every tree compared
    leaves = [leaf.name for leaf in t.iter_leaves()]
    # outgroup = leaf child of the root (the cell the tree was rooted on)
    root_cells = [c.name for c in t.children if c.is_leaf()]
    return leaves, ";".join(root_cells)


def draw_strip(ax, regions):
    # one column per cell, same as plt.bar(i, 1, width=1) per cell, but fast for 56k cells
    keys = list(region_colors) + ["other"]
    cmap = ListedColormap([region_colors[k] for k in region_colors] + ["black"])
    idx = np.array([keys.index(r) if r in region_colors else len(keys) - 1 for r in regions])
    ax.imshow(idx[None, :], aspect="auto", cmap=cmap, vmin=-0.5, vmax=len(keys) - 0.5,
              interpolation="nearest", extent=(-0.5, len(idx) - 0.5, 0, 1), rasterized=True)
    ax.set_xticks([])
    ax.set_yticks([])


def save(fig, stem, png_dpi=1000):
    fig.savefig(f"{stem}.png", format="png", dpi=png_dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(f"{stem}.pdf", format="pdf", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", default="neck")
    ap.add_argument("--levels", nargs="*", help="levels to plot, e.g. 0 30 (default: every rooted level)")
    ap.add_argument("--no-ladderize", action="store_true", help="keep the tree's own child order")
    ap.add_argument("--out", default=str(HERE / "levels"))
    a = ap.parse_args()

    base = Path(a.out) / a.name
    dirs = sorted(p for p in base.glob("p*") if (p / "results/trees/cell_tree.rooted").exists())
    if a.levels:
        want = {f"p{int(x.lstrip('p')):02d}" for x in a.levels}
        dirs = [d for d in dirs if d.name in want]
    if not dirs:
        raise SystemExit(f"no rooted tree under {base} - run root_trees.py first")

    outdir = base / "analysis" / "trees"
    outdir.mkdir(parents=True, exist_ok=True)

    results, summary = {}, []
    for d in dirs:
        print(f"{d.name}: reading rooted tree ...", flush=True)
        leaves, root_cell = leaf_order(d / "results/trees/cell_tree.rooted", not a.no_ladderize)
        regions = [leaf.split("_")[0] for leaf in leaves]
        results[d.name] = regions
        pd.DataFrame({"leaf": leaves, "region": regions}).to_csv(outdir / f"leaf_order_{d.name}.csv", index=False)
        counts = pd.Series(regions).value_counts()
        summary.append({"level": d.name, "cells": len(leaves), "root_cell": root_cell,
                        **{f"n_{k}": int(counts.get(k, 0)) for k in legend_order}})

        fig, ax = plt.subplots(figsize=(12, 4))
        draw_strip(ax, regions)
        ax.set_xlabel("Cells (ordered by tree traversal)")
        ax.set_title(f"Distribution of Cells by Region (Ordered by Tree) - {d.name}")
        save(fig, outdir / f"distribution_by_region_{d.name}")
        print(f"{d.name}: {len(leaves)} cells, rooted on {root_cell}")

    # all levels stacked, one row per level
    fig, axes = plt.subplots(len(results), 1, figsize=(12, 1.2 * len(results) + 0.6), squeeze=False)
    for ax, (lvl, regions) in zip(axes[:, 0], results.items()):
        draw_strip(ax, regions)
        ax.set_ylabel(lvl, rotation=0, ha="right", va="center")
    axes[-1, 0].set_xlabel("Cells (ordered by tree traversal)")
    axes[0, 0].set_title("Distribution of Cells by Region (Ordered by Tree)")
    save(fig, outdir / "distribution_by_region_all_levels")

    fig = plt.figure(figsize=(6, 2))
    plt.legend(handles=[Patch(facecolor=region_colors[k], label=region_mapping[k]) for k in legend_order],
               loc="center", frameon=False, ncol=2)
    plt.axis("off")
    fig.savefig(outdir / "region_legend.pdf", format="pdf", dpi=300, bbox_inches="tight")
    plt.close(fig)

    pd.DataFrame(summary).to_csv(outdir / "tree_regions_summary.csv", index=False)
    print(f"\nwritten to {outdir}")


if __name__ == "__main__":
    main()
