"""
plot_pearson.py
Pearson r between each cell's true and estimated signature share vs mutations per cell,
biopsy3 experiment. One figure per signature and clone arrangement:
  SBS1  : clone 1 (SBS5 clone), clone 2 (SBS26 clone), all cells pooled
          (trunk / SBS1-only cells alone are not shown: their true SBS1 share is always 1,
          so r is undefined for them)
  SBS5  : clone 1 cells only
  SBS26 : clone 2 cells only
Shares: estimate / total exposure of the cell (all signatures, incl. non-simulated);
truth = the cell's UNCAPPED real share (same at every cap). Cells of the 3 replicates pooled
(1 run when uncapped). SBS12 is counted as SBS26 (merged results).
Cap where the signature is never called (estimate 0 everywhere) -> r undefined, marked x.

  python scripts/plot_pearson.py   # -> <output>/figures/pearson/
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

M = os.environ.get("SIM_OUT", os.getcwd())   # pipeline output folder (set by main.sh)
RUNS = [(os.environ.get("RUN_NAME", "sbs26_large_clone"), "large", "clone 2 (SBS26) = larger clone (369 cells)")]
OUT = os.path.join(M, "figures/pearson")
os.makedirs(OUT, exist_ok=True)
INK, INK2, MUTED, GRID, NOCALL = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e0", "#d03b3b"
# line label -> (cell group, colour, linestyle); group None = all cells
LINES = {"SBS1": [("clone 1", "clone_SBS5", "#eb6834", "-"),
                  ("clone 2", "clone_SBS26", "#2a78d6", "-"),
                  ("all cells", None, INK, "--")],
         "SBS5": [("clone 1", "clone_SBS5", "#eb6834", "-")],
         "SBS26": [("clone 2", "clone_SBS26", "#2a78d6", "-")]}
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False})
# editable text in Adobe Illustrator
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'


def pearson(x, y):
    return float(np.corrcoef(x, y)[0, 1]) if len(x) > 2 and x.std() > 0 and y.std() > 0 else np.nan


summary = []
for X, tag, title in RUNS:
    pc = pd.read_csv(f"{M}/{X}/results_merged_SBS12/leaf_per_cell.csv")
    pc["tot_true"] = pc[["true_SBS1", "true_SBS5", "true_SBS26"]].sum(1)
    pc["tot_est"] = pc[["est_SBS1", "est_SBS5", "est_SBS26", "est_non_simulated"]].sum(1)
    full_max = pc.loc[pc.mask_pct == 0, "tot_true"].max()
    pc = pc[(pc.mask_pct == 0) | (pc.mask_pct < full_max)].copy()   # caps above every burden = uncapped
    caps = sorted(pc.mask_pct.unique(), key=lambda c: -(c or 10**6))
    med = int(pc.loc[pc.mask_pct == 0, "tot_true"].median())
    labels = [f"all\n(~{med})"] + [str(int(c)) for c in caps[1:]]
    x = np.arange(len(caps))

    for sig, lines in LINES.items():
        true_c = pc[f"true_{sig}"] / pc["tot_true"]
        # each cell's real share before any downsampling (same at every cap, as in the boxplot)
        u = pc.mask_pct == 0
        true = pc["cell"].map(true_c[u].set_axis(pc.loc[u, "cell"]))
        est =(pc[f"est_{sig}"] / pc["tot_est"]).fillna(0)
        fig, ax = plt.subplots(figsize=(4.6, 3.6))
        nocall = []
        for lab, grp, col, ls in lines:
            cells = (pc.true_clone == grp) if grp else pd.Series(True, index=pc.index)
            r = []
            for i, c in enumerate(caps):
                sel = (pc.mask_pct == c) & cells
                r.append(pearson(true[sel], est[sel]))
                if est[pc.mask_pct == c].max() == 0 and i not in nocall:
                    nocall.append(i)
                summary.append({"arrangement": tag, "signature": sig, "cells": lab, "cap": int(c),
                                "n_cell_runs": int(sel.sum()), "pearson_r": r[-1]})
            if grp is None:   # all cells: wide light band underneath so clone lines stay visible
                ax.plot(x, r, color="#bdbcb7", lw=7, solid_capstyle="round", zorder=1, label=lab)
            else:
                ax.plot(x, r, color=col, lw=2.2, ls=ls, marker="o", ms=5,
                        markeredgecolor="white", zorder=3, label=lab)
        if nocall:
            ax.scatter(nocall, [0] * len(nocall), marker="x", color=NOCALL, s=45, zorder=3,
                       label=f"{sig} not called")
        ax.axhline(0, color=MUTED, lw=0.8)
        ax.set_xticks(x); ax.set_xticklabels(labels)
        ax.set_ylim(-0.05, 1.05)
        ax.set_xlabel("mutations per cell (cap)"); ax.set_ylabel("Pearson r (true vs estimated)")
        ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)
        ax.legend(frameon=False, fontsize=8.5, loc="lower left")
        ax.set_title(f"{sig} share per cell\n{title}", loc="left", fontsize=10, color=INK)
        fig.tight_layout()
        for ext in ("png", "pdf"):
            fig.savefig(f"{OUT}/{tag}_{sig}_pearson.{ext}", dpi=200)
        plt.close(fig)

pd.DataFrame(summary).round(4).to_csv(f"{OUT}/pearson_values.csv", index=False)
print("wrote", sorted(f for f in os.listdir(OUT) if f.endswith(".png")))
