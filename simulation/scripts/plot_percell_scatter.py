"""
plot_percell_scatter.py
Per-cell TRUE vs ESTIMATED share, biopsy3 experiment. One figure per signature
(SBS26, SBS5, SBS1) and clone arrangement:
  * one scatter panel per cap: x = cell's true share, y = SigDyn's estimated share,
    one point per cell (uncapped panel: cells of a lineage are identical, so points get
    a small display-only jitter and clusters are labelled with their cell counts),
    every cell of every replicate (3 replicates pooled; 1 run when uncapped),
    coloured by the cell's true group, dashed = identity (perfect estimate)
  * last panel: Pearson r vs mutations per cell
      solid  = all cells (largely driven by clone vs non-clone separation)
      dotted = carrier cells only (does each clone cell get ITS OWN share right?)
Shares: estimate / total exposure of the cell (all signatures, incl. non-simulated);
truth / the cell's surviving mutations. SBS12 is counted as SBS26 (merged results).

  python scripts/plot_percell_scatter.py   # -> <output>/figures/percell_scatter/
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

M = os.environ.get("SIM_OUT", os.getcwd())   # pipeline output folder (set by main.sh)
RUNS = [(os.environ.get("RUN_NAME", "sbs26_large_clone"), "large", "SBS26 = larger clone (369 cells)")]
OUT = os.path.join(M, "figures/percell_scatter")
os.makedirs(OUT, exist_ok=True)
GROUPS = {"clone_SBS26": ("SBS26 clone", "#2a78d6"),
          "clone_SBS5": ("SBS5 clone", "#eb6834"),
          "SBS1-only": ("SBS1-only", "#8a8984")}
CARRIER = {"SBS26": "clone_SBS26", "SBS5": "clone_SBS5", "SBS1": None}
JITTER = 0.01   # uncapped panel only: identical cells are spread by up to +/- this (display only)
LIM = {"SBS26": 0.65, "SBS5": 0.65, "SBS1": 1.1}   # a bit past the largest values
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e0"
plt.rcParams.update({"font.family": "sans-serif", "font.size": 9, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42})


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
    cap_lab = {c: (f"all (~{med})" if c == 0 else str(int(c))) for c in caps}

    for sig in ["SBS26", "SBS5", "SBS1"]:
        true = pc[f"true_{sig}"] / pc["tot_true"]
        est = (pc[f"est_{sig}"] / pc["tot_est"]).fillna(0)
        grp = CARRIER[sig]
        carrier = (pc.true_clone == grp) if grp else pd.Series(True, index=pc.index)
        n = len(caps) + 1
        ncol = 4
        fig, axes = plt.subplots(int(np.ceil(n / ncol)), ncol, figsize=(13, 3.2 * np.ceil(n / ncol)))
        axes = axes.ravel()
        r_all, r_car = [], []
        for i, c in enumerate(caps):
            ax = axes[i]; lv = pc.mask_pct == c
            for g, (lab, col) in GROUPS.items():
                sel = lv & (pc.true_clone == g)
                if c == 0:
                    # uncapped: cells of one lineage are identical and sit on the same spot;
                    # draw every cell with a small random jitter (display only; r uses real values)
                    # and label each cluster with its number of cells
                    jit = np.random.default_rng(0)
                    ax.scatter(true[sel] + jit.uniform(-JITTER, JITTER, sel.sum()),
                               est[sel] + jit.uniform(-JITTER, JITTER, sel.sum()),
                               s=5, alpha=0.35, color=col, lw=0)
                    pts = pd.DataFrame({"t": true[sel].round(4), "e": est[sel].round(4)}).value_counts()
                    for (t, e), k in pts.items():
                        ax.annotate(str(k), (t, e), xytext=(7, -10), textcoords="offset points",
                                    fontsize=7, color=col)
                else:
                    ax.scatter(true[sel], est[sel], s=7, alpha=0.35, color=col, lw=0)
            lim = LIM[sig]
            ax.plot([0, lim], [0, lim], color=INK2, ls="--", lw=0.9)
            ax.set_xlim(-0.03, lim); ax.set_ylim(-0.03, lim)
            if c == 0:
                k = len(pd.DataFrame({"t": true[lv].round(4), "e": est[lv].round(4)}).drop_duplicates())
                ax.text(0.03, 0.80, f"{k} distinct values, all {int(lv.sum())} cells shown\n"
                        f"(jitter ±{JITTER}; numbers = cells per cluster)", transform=ax.transAxes,
                        ha="left", va="top", fontsize=7.5, color=INK2)
            ra, rc = pearson(true[lv], est[lv]), pearson(true[lv & carrier], est[lv & carrier])
            r_all.append(ra); r_car.append(rc)
            if est[lv].max() == 0:
                txt = f"{sig} not called\n(estimate 0 in every cell; r undefined)"
            else:
                txt = f"r all = {ra:.2f}" + ("" if grp is None else f"\nr {GROUPS[grp][0]} = {rc:.2f}")
            ax.text(0.03, 0.97, txt, transform=ax.transAxes, va="top", fontsize=8.5, color=INK)
            ax.set_title(f"cap {cap_lab[c]}  ({int(lv.sum())} cell-runs)", loc="left", fontsize=9, color=INK)
            ax.set_xlabel(f"true {sig} share"); ax.set_ylabel(f"estimated {sig} share")
            ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)
            summary.append({"arrangement": tag, "signature": sig, "cap": int(c),
                            "n_cell_runs": int(lv.sum()), "pearson_all_cells": ra,
                            "pearson_carriers": rc if grp else np.nan,
                            "mean_abs_error_carriers": float((est[lv & carrier] - true[lv & carrier]).abs().mean())})
        # r vs cap
        ax = axes[len(caps)]; x = np.arange(len(caps))
        ax.plot(x, r_all, color=INK, lw=2, marker="o", ms=4, label="all cells")
        if grp:
            ax.plot(x, r_car, color=GROUPS[grp][1], lw=2, ls=":", marker="o", ms=4,
                    label=f"{GROUPS[grp][0]} cells only")
        nc = [i for i, c in enumerate(caps) if np.isnan(r_all[i])]
        if nc:
            ax.scatter(nc, [0] * len(nc), marker="x", color="#d03b3b", s=40, zorder=3,
                       label=f"{sig} not called")
        ax.axhline(0, color=MUTED, lw=0.8)
        ax.set_xticks(x); ax.set_xticklabels([cap_lab[c].split(" ")[0] for c in caps])
        ax.set_ylim(-0.5, 1.05); ax.set_xlabel("mutations per cell (cap)"); ax.set_ylabel("Pearson r")
        ax.set_title("Pearson r (true vs estimated)", loc="left", fontsize=9, color=INK)
        ax.legend(frameon=False, fontsize=8, loc="lower left"); ax.grid(axis="y", color=GRID, lw=0.6)
        for ax in axes[len(caps) + 1:]:
            ax.axis("off")
        from matplotlib.lines import Line2D
        handles = [Line2D([], [], ls="", marker="o", ms=6, color=col, alpha=0.7, label=lab)
                   for lab, col in GROUPS.values()]
        fig.legend(handles=handles, frameon=False, fontsize=9, ncol=3, loc="upper right",
                   bbox_to_anchor=(0.99, 1.0))
        fig.suptitle(f"{sig}: true vs estimated share per cell — {title}", x=0.01, ha="left",
                     fontsize=11.5, color=INK)
        fig.tight_layout()
        for ext in ("png", "pdf"):
            fig.savefig(f"{OUT}/{tag}_{sig}_scatter.{ext}", dpi=200)
        plt.close(fig)

pd.DataFrame(summary).round(4).to_csv(f"{OUT}/percell_scatter_values.csv", index=False)
print("wrote", sorted(f for f in os.listdir(OUT) if f.endswith(".png")))
