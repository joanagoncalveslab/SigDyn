#!/usr/bin/env python3
"""
SH-like branch support distribution over the missingness levels
(the 'Support distribution' histogram of tree_support/unrooted_support_breast_neck2.ipynb).
Same support values as the tree-support summary of analyze.py (tree_support_vs_missingness.csv).

Per level, from the unrooted FastTree tree, two histograms (30 bins, as the notebook):
  absolute : # internal branches
  relative : fraction of internal branches
red line = mean support; grey dashed = 0.5 / 0.7 / 0.9 (as the notebook).

Reads   <out>/<name>/pXX/results/trees/cell_tree.unrooted
Writes  <out>/<name>/analysis/support/
          support_distribution_pXX.pdf/png       absolute + relative, one per level
          support_distribution_all_levels.pdf/png   rows = levels
          support_values_pXX.csv                 support per internal branch

  python support_distribution.py --name neck
"""
import argparse
from pathlib import Path

import dendropy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
COLORS = {"breast": "#C44E52", "neck": "#4C72B0"}   # from the notebook
BINS = np.linspace(0, 1, 31)                         # 30 bins
plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 10,
                     "axes.grid": True, "grid.alpha": 0.25})


def support_of(node):
    # FastTree stores SH-like support as the internal-node label
    try:
        return float(node.label)
    except (TypeError, ValueError):
        return None


def supports(tree_file):
    """Support of every internal branch, as tree_support() in analyze.py."""
    t = dendropy.Tree.get(path=str(tree_file), schema="newick", preserve_underscores=True)
    t.encode_bipartitions()
    sup = []
    for e in t.preorder_edge_iter():
        n = e.head_node
        if n is None or n.is_leaf() or e.bipartition is None:
            continue
        s = support_of(n)
        if s is not None:
            sup.append(s)
    return np.array(sup)


def hist(ax, sup, relative, color, title):
    w = np.full(len(sup), 1 / len(sup)) if relative else None
    ax.hist(sup, bins=BINS, weights=w, color=color, alpha=0.8, edgecolor="white")
    for t_ in (0.5, 0.7, 0.9):
        ax.axvline(t_, color="grey", ls="--", lw=0.8)
    m = sup.mean()
    ax.axvline(m, color="red", lw=1.8, label=f"mean = {m:.3f}")
    ax.set_xlim(0, 1)
    ax.set_xlabel("SH-like local support")
    ax.set_ylabel("fraction of internal branches" if relative else "# internal branches")
    ax.set_title(title)
    ax.legend(loc="upper center", framealpha=1, edgecolor="none")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", default="neck")
    ap.add_argument("--out", default=str(HERE / "levels"))
    a = ap.parse_args()

    base = Path(a.out) / a.name
    dirs = sorted(p for p in base.glob("p*") if (p / "results/trees/cell_tree.unrooted").exists())
    if not dirs:
        raise SystemExit(f"no unrooted tree under {base} - run run_levels.py stage1 first")
    pct = {}
    if (base / "levels_summary.csv").exists():
        s = pd.read_csv(base / "levels_summary.csv")
        pct = {f"p{int(l):02d}": 100 * f for l, f in zip(s["level"], s["frac_missing"])}
    outdir = base / "analysis" / "support"
    outdir.mkdir(parents=True, exist_ok=True)
    color = COLORS.get(a.name, "#4C72B0")

    sups = {}
    for d in dirs:
        sup = supports(d / "results/trees/cell_tree.unrooted")
        sups[d.name] = sup
        pd.DataFrame({"support": sup}).to_csv(outdir / f"support_values_{d.name}.csv", index=False)
        label = f"{d.name} ({pct.get(d.name, float(d.name[1:])):.0f} % missing)"

        fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
        hist(ax[0], sup, False, color, f"Support distribution - {label}")
        hist(ax[1], sup, True, color, f"Support distribution (relative) - {label}")
        fig.tight_layout()
        for ext in ("pdf", "png"):
            fig.savefig(outdir / f"support_distribution_{d.name}.{ext}", dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"{d.name}: {len(sup)} internal branches, mean {sup.mean():.4f}, median {np.median(sup):.3f}")

    # all levels: rows = levels, columns = absolute / relative; same y-axis per column
    n = len(sups)
    fig, ax = plt.subplots(n, 2, figsize=(11, 2.9 * n), squeeze=False, sharex=True, sharey="col")
    for i, (lv, sup) in enumerate(sups.items()):
        label = f"{lv} ({pct.get(lv, float(lv[1:])):.0f} % missing)"
        hist(ax[i, 0], sup, False, color, label)
        hist(ax[i, 1], sup, True, color, label)
        if i < n - 1:
            ax[i, 0].set_xlabel("")
            ax[i, 1].set_xlabel("")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(outdir / f"support_distribution_all_levels.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"\nwritten to {outdir}")


if __name__ == "__main__":
    main()
