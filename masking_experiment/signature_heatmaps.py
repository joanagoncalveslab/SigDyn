#!/usr/bin/env python3
"""
Fig 4a (right) over the missingness levels: per-cell signature activities, cells in tree order.
Same as Case_studies/Figure4_LSCC/Fig4a_right_signature_heatmaps.ipynb, run on every level:
leaves of the rooted tree in order, one viridis column per signature with its own colour scale.
Cells without activities are left out, as in the notebook.
Signatures: by default the ones active at that level, i.e. SigDyn's own list
(results/signatures/Custom_COSMIC_SBS_v3.4.txt = signatures in >= 10 % of cells);
--signatures uses one fixed list for every level instead (e.g. the paper's).

The tree is ladderized, as the paper's rerooted_tree_no_germline.nwk was, so the rows follow the
same cell order as tree_regions.py (Fig 4a left) for that level.

Reads   <out>/<name>/pXX/results/trees/cell_tree.rooted
        <out>/<name>/pXX/results/signatures/Custom_COSMIC_SBS_v3.4.txt   (active signatures)
        <out>/<name>/pXX/results/signatures/per_leaf/Assignment_Solution/Activities/Assignment_Solution_Activities.txt
Writes  <out>/<name>/analysis/trees/signature_heatmaps_pXX.pdf

  python signature_heatmaps.py --name neck
  python signature_heatmaps.py --name neck --signatures SBS1 SBS11 SBS19 SBS54 SBS87 SBS5 SBS12 SBS58   # the paper's
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from ete3 import Tree

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'

HERE = Path(__file__).resolve().parent
ACTIVITIES = "results/signatures/per_leaf/Assignment_Solution/Activities/Assignment_Solution_Activities.txt"
ACTIVE = "results/signatures/Custom_COSMIC_SBS_v3.4.txt"


def active_signatures(level_dir):
    return pd.read_csv(level_dir / ACTIVE, sep="\t", nrows=0).columns[1:].tolist()


def heatmaps(df_interest, title, path):
    signatures = df_interest.columns
    fig, axes = plt.subplots(1, len(signatures), figsize=(15 * max(len(signatures), 8) / 8, 10), sharey=True, squeeze=False)
    for ax, sig in zip(axes[0], signatures):
        hm = sns.heatmap(
            df_interest[[sig]],
            yticklabels=False,
            xticklabels=True,
            cmap="viridis",
            cbar=True,
            ax=ax,
            rasterized=True,
        )
        cbar = hm.collections[0].colorbar
        cbar.ax.tick_params(labelsize=16)
        ax.set_title(sig, fontsize=16)
        ax.tick_params(axis="x", labelrotation=90, labelsize=16)
        ax.set_ylabel("")
    plt.subplots_adjust(wspace=0.6)
    fig.suptitle(title, fontsize=18)
    fig.savefig(path, format="pdf", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", default="neck")
    ap.add_argument("--signatures", nargs="+", help="fixed list for every level (default: each level's active ones)")
    ap.add_argument("--out", default=str(HERE / "levels"))
    a = ap.parse_args()

    base = Path(a.out) / a.name
    dirs = sorted(p for p in base.glob("p*")
                  if (p / "results/trees/cell_tree.rooted").exists() and (p / ACTIVITIES).exists())
    if not dirs:
        raise SystemExit(f"no level with a rooted tree and per-leaf activities under {base}")
    outdir = base / "analysis" / "trees"
    outdir.mkdir(parents=True, exist_ok=True)

    for d in dirs:
        tree = Tree(str(d / "results/trees/cell_tree.rooted"), format=1)
        tree.ladderize()
        sig = pd.read_csv(d / ACTIVITIES, sep="\t").set_index("Samples")
        sigs = a.signatures or active_signatures(d)
        missing = [s for s in sigs if s not in sig.columns]
        if missing:
            raise SystemExit(f"{d.name}: signatures not in the activities: {missing}")

        leaves = [lf.name for lf in tree.iter_leaves()]
        order = [l for l in leaves if l in sig.index]
        df_interest = sig.loc[order, sigs]

        heatmaps(df_interest, d.name, outdir / f"signature_heatmaps_{d.name}.pdf")
        print(f"{d.name}: {len(order)} of {len(leaves)} leaves have activities; {len(sigs)} signatures: {' '.join(sigs)}")

    print(f"\nwritten to {outdir}")


if __name__ == "__main__":
    main()
