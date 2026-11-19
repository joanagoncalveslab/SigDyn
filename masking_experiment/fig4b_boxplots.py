#!/usr/bin/env python3
"""
Fig 4b over the missingness levels: SBS32 activity per region, in SBS32-active cells.
Same as Case_studies/Figure4_LSCC/Fig4_b_box_plots.ipynb, run on every level's per-leaf activities.
The H32/L32 split of the DGEA (Fig 4c/d) is SBS32 > 0 / = 0 on these same activities.

Reads   <out>/<name>/pXX/results/signatures/per_leaf/Assignment_Solution/Activities/Assignment_Solution_Activities.txt
Writes  <out>/<name>/dgea/pXX/                  (next to that level's volcano and GO plots)
          SBS32_activity_boxplots.pdf          as the notebook
          SBS32_region_summary.csv             overall + per region: % active, mean (all / nonzero), n cells
          SBS32_mannwhitney.csv                Mann-Whitney U between regions (SBS32-active cells), as in the paper
        <out>/<name>/dgea/
          SBS32_activity_boxplots_all_levels.pdf    levels side by side
          SBS32_region_summary_vs_missingness.csv, SBS32_mannwhitney_vs_missingness.csv

  python fig4b_boxplots.py --name neck
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'

HERE = Path(__file__).resolve().parent
ACTIVITIES = "results/signatures/per_leaf/Assignment_Solution/Activities/Assignment_Solution_Activities.txt"

# from the notebook
region_label_map = {
    "neck15": "Normal laryngeal mucosal epithelia (N)",
    "neck17": "Margins of cancer (R)",
    "neck18": "Carcinoma in situ (T)",
    "neck16": "Lymph nodes metastasis (L)",
}
region_order_full = [
    "Normal laryngeal mucosal epithelia (N)",
    "Margins of cancer (R)",
    "Carcinoma in situ (T)",
    "Lymph nodes metastasis (L)",
]
region_order = ["N", "R", "T", "L"]


def load(path):
    activity_matrix = pd.read_csv(path, sep="\t", index_col=0)
    activity_matrix["region"] = (
        activity_matrix.index.astype(str).str.extract(r"(neck1[5-8])", expand=False).map(region_label_map)
    )
    activity_matrix["SBS32_active"] = activity_matrix["SBS32"] > 0
    return activity_matrix


def summary(activity_matrix):
    def stats(df):
        return {
            "percent_cells_SBS32_active": 100 * df["SBS32_active"].mean(),
            "average_SBS32_all_cells": df["SBS32"].mean(),
            "average_SBS32_nonzero": df.loc[df["SBS32"] > 0, "SBS32"].mean(),
            "median_SBS32_nonzero": df.loc[df["SBS32"] > 0, "SBS32"].median(),
            "n_cells": len(df),
            "n_cells_SBS32_active": int((df["SBS32"] > 0).sum()),
        }
    rows = [{"region": "overall", **stats(activity_matrix)}]
    rows += [{"region": r, **stats(activity_matrix[activity_matrix["region"] == r])} for r in region_order_full]
    return pd.DataFrame(rows)


def mann_whitney(activity_matrix):
    """Two-sided Mann-Whitney U on SBS32 of SBS32-active cells (the boxplot data), every pair of regions.
    r = |rank-biserial correlation| = |1 - 2U / (n1 n2)|. At p00 this gives the paper's
    T vs N p = 3.4e-147, r = 0.25 and L vs N r = 0.22 (p 4.9e-116 vs 5.9e-116 in the paper)."""
    act = activity_matrix.loc[activity_matrix["SBS32"] > 0]
    short = dict(zip(region_order_full, region_order))
    rows = []
    for i, ra in enumerate(region_order_full):
        for rb in region_order_full[:i]:
            a = act.loc[act["region"] == ra, "SBS32"].values
            b = act.loc[act["region"] == rb, "SBS32"].values
            U, p = mannwhitneyu(a, b, alternative="two-sided")
            rows.append({"comparison": f"{short[ra]} vs {short[rb]}", "n1": len(a), "n2": len(b),
                         "median1": np.median(a), "median2": np.median(b), "U": U, "p": p,
                         "r_rank_biserial": abs(1 - 2 * U / (len(a) * len(b)))})
    return pd.DataFrame(rows)


def boxplot(ax, activity_matrix, title):
    # Keep only SBS32-positive cells
    plot_df = activity_matrix.loc[activity_matrix["SBS32"] > 0]
    data = [plot_df.loc[plot_df["region"] == r, "SBS32"].values for r in region_order_full]

    positions = np.arange(len(data)) * 0.25  # <--- spacing factor (smaller = closer)
    ax.boxplot(
        data,
        positions=positions,
        tick_labels=region_order,
        showfliers=True,
        patch_artist=False,
        widths=0.2,
        boxprops=dict(linewidth=0.6),
        whiskerprops=dict(linewidth=0.6),
        capprops=dict(linewidth=0.6),
        medianprops=dict(linewidth=0.6),
        flierprops=dict(marker="o", markersize=0.5, alpha=0.5),
    )
    ax.set_ylabel("SBS32 activity", fontsize=15)
    ax.set_title(title, fontsize=15, pad=12)
    ax.tick_params(axis="x", rotation=35, labelsize=13)
    ax.tick_params(axis="y", labelsize=13)
    ax.set_ylim(0, 99)
    ax.set_xlim(positions[0] - 0.14, positions[-1] + 0.14)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.4)
    ax.spines["left"].set_linewidth(0.4)
    ax.grid(axis="y", linestyle="--", linewidth=0.4, alpha=0.4)
    return plot_df["SBS32"].max() if len(plot_df) else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", default="neck")
    ap.add_argument("--out", default=str(HERE / "levels"))
    a = ap.parse_args()

    base = Path(a.out) / a.name
    dirs = sorted(p for p in base.glob("p*") if (p / ACTIVITIES).exists())
    if not dirs:
        raise SystemExit(f"no per-leaf activities under {base} - run run_levels.py stage2 first")
    dg = base / "dgea"

    mats, summ, tests = {}, [], []
    for d in dirs:
        m = load(d / ACTIVITIES)
        mats[d.name] = m
        s = summary(m)
        od = dg / d.name
        od.mkdir(parents=True, exist_ok=True)
        s.to_csv(od / "SBS32_region_summary.csv", index=False)
        summ.append(s.assign(level=d.name))
        mw = mann_whitney(m)
        mw.to_csv(od / "SBS32_mannwhitney.csv", index=False)
        tests.append(mw.assign(level=d.name))

        fig, ax = plt.subplots(figsize=(8, 6))
        vmax = boxplot(ax, m, "Distribution of SBS32 activity across regions\n(SBS32-active cells)")
        plt.tight_layout()
        fig.savefig(od / "SBS32_activity_boxplots.pdf", format="pdf", bbox_inches="tight")
        plt.close(fig)
        if vmax > 99:
            print(f"  WARNING {d.name}: SBS32 up to {vmax}, above the notebook's ylim of 99")
        pct = s.set_index("region")["percent_cells_SBS32_active"].round(1)
        print(f"{d.name}: % SBS32-active overall {pct['overall']}, "
              + ", ".join(f"{k} {v}" for k, v in zip(region_order, pct[region_order_full])))
        for _, t in mw.iterrows():
            print(f"    {t['comparison']}: p = {t['p']:.2g}, r = {t['r_rank_biserial']:.2f}")

    # all levels side by side
    fig, axes = plt.subplots(1, len(mats), figsize=(4.5 * len(mats), 6), sharey=True, squeeze=False)
    for j, (ax, (lvl, m)) in enumerate(zip(axes[0], mats.items())):
        boxplot(ax, m, lvl)
        if j:
            ax.set_ylabel("")
    fig.suptitle("Distribution of SBS32 activity across regions (SBS32-active cells)", fontsize=15)
    plt.tight_layout()
    fig.savefig(dg / "SBS32_activity_boxplots_all_levels.pdf", format="pdf", bbox_inches="tight")
    plt.close(fig)

    out = pd.concat(summ)
    out = out[["level"] + [c for c in out.columns if c != "level"]]
    out.to_csv(dg / "SBS32_region_summary_vs_missingness.csv", index=False)
    tt = pd.concat(tests)
    tt[["level"] + [c for c in tt.columns if c != "level"]].to_csv(dg / "SBS32_mannwhitney_vs_missingness.csv", index=False)
    print(f"\nwritten to {dg}")


if __name__ == "__main__":
    main()
