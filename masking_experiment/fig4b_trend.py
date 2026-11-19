#!/usr/bin/env python3
"""
How the paper's Fig 4b result (T and L have higher SBS32 than N) changes over the missingness levels.
Reads the tables written by fig4b_boxplots.py (run that first).

  signed rank-biserial r of T vs N and L vs N (> 0 = region higher than N, as in the paper);
         filled marker = Mann-Whitney p < 0.05, open = not significant

Reads   <out>/<name>/dgea/SBS32_mannwhitney_vs_missingness.csv
        <out>/<name>/levels_summary.csv     (real % missing per level)
Writes  <out>/<name>/dgea/SBS32_vs_N_over_levels.pdf/png

  python fig4b_trend.py --name neck
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'

HERE = Path(__file__).resolve().parent
ALPHA = 0.05

# region colours of Fig4a_left_trees.ipynb
colors = {"N": "#E69F00", "R": "#56B4E9", "T": "#D55E00", "L": "#CC79A7"}


def style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.6)
    ax.spines["left"].set_linewidth(0.6)
    ax.grid(axis="y", linestyle="--", linewidth=0.4, alpha=0.4)
    ax.tick_params(labelsize=11)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", default="neck")
    ap.add_argument("--out", default=str(HERE / "levels"))
    a = ap.parse_args()

    base = Path(a.out) / a.name
    dg = base / "dgea"
    mw = pd.read_csv(dg / "SBS32_mannwhitney_vs_missingness.csv")
    lv = pd.read_csv(base / "levels_summary.csv")
    pct = {f"p{int(l):02d}": 100 * f for l, f in zip(lv["level"], lv["frac_missing"])}

    # signed rank-biserial: 2U/(n1 n2) - 1, > 0 when the first region is higher than N
    mw["r_signed"] = 2 * mw["U"] / (mw["n1"] * mw["n2"]) - 1
    mw["pct"] = mw["level"].map(pct)

    fig, ax1 = plt.subplots(figsize=(5.5, 4.2))

    # ---- left: effect size vs N ----
    ax1.axhline(0, color="gray", lw=0.8)
    for reg in ("T", "L"):
        d = mw[mw["comparison"] == f"{reg} vs N"].sort_values("pct")
        ax1.plot(d["pct"], d["r_signed"], color=colors[reg], lw=2, zorder=2)
        sig = d["p"] < ALPHA
        ax1.scatter(d.loc[sig, "pct"], d.loc[sig, "r_signed"], s=60, color=colors[reg],
                    edgecolor="white", linewidth=1.5, zorder=3)
        ax1.scatter(d.loc[~sig, "pct"], d.loc[~sig, "r_signed"], s=60, facecolor="white",
                    edgecolor=colors[reg], linewidth=1.5, zorder=3)
        last = d.iloc[-1]
        ax1.annotate(f"{reg} vs N", (last["pct"], last["r_signed"]), xytext=(8, 0),
                     textcoords="offset points", va="center", fontsize=11)
    ax1.text(0.01, 0.98, "region higher than N", transform=ax1.transAxes, va="top", fontsize=9, color="dimgray")
    ax1.text(0.01, 0.02, "region lower than N", transform=ax1.transAxes, va="bottom", fontsize=9, color="dimgray")
    lim = max(0.35, mw["r_signed"].abs().max() + 0.03)
    ax1.set_ylim(-lim, lim)
    ax1.set_ylabel("Effect size vs N (rank-biserial r)", fontsize=12)
    ax1.set_title("SBS32 activity vs normal (N)", fontsize=13)
    ax1.scatter([], [], s=60, color="gray", label=f"p < {ALPHA}")
    ax1.scatter([], [], s=60, facecolor="white", edgecolor="gray", linewidth=1.5, label="not significant")
    ax1.legend(frameon=False, fontsize=9, loc="upper right")

    for ax in (ax1,):
        ax.set_xlabel("% of genotype matrix missing", fontsize=12)
        ax.set_xticks(sorted(pct.values()), [f"{v:.0f}" for v in sorted(pct.values())])
        ax.set_xlim(-5, max(pct.values()) + 14)
        style(ax)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(dg / f"SBS32_vs_N_over_levels.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(mw.loc[mw["comparison"].isin(["T vs N", "L vs N"]), ["level", "comparison", "r_signed", "p"]].to_string(index=False))
    print(f"\nwritten {dg / 'SBS32_vs_N_over_levels.pdf'}")


if __name__ == "__main__":
    main()
