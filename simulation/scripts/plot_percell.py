"""
plot_percell.py
Per-cell view of the biopsy3 experiment: one figure per signature (SBS26, SBS5, SBS1)
and clone arrangement, three panels vs mutations per cell:
  1  box plots of each cell's estimated share (coloured) next to its real share (grey;
     the cell's uncapped real share, so the same at every cap)
     (cells of the 3 replicates pooled; box = 25-75%, whiskers = 5-95%)
  2  % of cells called correctly (call = estimated share >= --call, default 10%):
       solid  = detected in cells that truly carry the signature
       dotted = not detected in cells that don't (not shown for SBS1: every cell carries it)
  3  % of carrier cells whose estimated share is within +/- --tol (default 0.05)
     of that cell's own true share
Estimated share = estimate / total exposure of the cell (all signatures, incl. non-simulated).
Carrier cells: SBS26 -> SBS26 clone, SBS5 -> SBS5 clone, SBS1 -> all cells.

  python scripts/plot_percell.py        # -> <output>/figures/percell/*.png + .pdf + csv
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

M = os.environ.get("SIM_OUT", os.getcwd())   # pipeline output folder (set by main.sh)
ap = argparse.ArgumentParser()
ap.add_argument("--call", type=float, default=0.10)
ap.add_argument("--tol", type=float, default=0.05)
a = ap.parse_args()
RUNS = [(os.environ.get("RUN_NAME", "sbs26_large_clone"), "large", "SBS26 = larger clone (369 cells)")]
OUT = os.path.join(M, "figures/percell")
os.makedirs(OUT, exist_ok=True)
SIGS = {"SBS26": ("clone_SBS26", "#2a78d6", 0.5),
        "SBS5": ("clone_SBS5", "#eb6834", 0.5),
        "SBS1": (None, "#4a3aa7", 1.0)}
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e0"
plt.rcParams.update({"font.family": "sans-serif", "font.size": 10, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False})
# editable text in Adobe Illustrator
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
from matplotlib.patches import Patch

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
    for sig, (grp, col, ymax) in SIGS.items():
        est = pc[f"est_{sig}"] / pc["tot_est"]
        true = pc[f"true_{sig}"] / pc["tot_true"]
        # each cell's real share before any downsampling (constant across caps)
        true0 = pc["cell"].map(true[pc.mask_pct == 0].set_axis(pc.loc[pc.mask_pct == 0, "cell"]))
        carrier = (pc.true_clone == grp) if grp else pd.Series(True, index=pc.index)
        box, tbox, true_mean, det, absn, within = [], [], [], [], [], []
        for c in caps:
            lv = pc.mask_pct == c
            box.append(est[lv & carrier].values)
            tbox.append(true0[lv & carrier & (pc.rep == pc.loc[lv, "rep"].min())].values)
            true_mean.append(true[lv & carrier].mean())
            det.append(100 * (est[lv & carrier] >= a.call).mean())
            absn.append(100 * (est[lv & ~carrier] < a.call).mean() if grp else np.nan)
            within.append(100 * ((est[lv & carrier] - true[lv & carrier]).abs() <= a.tol).mean())
            summary.append({"arrangement": tag, "signature": sig, "cap": int(c),
                            "n_cells_pooled": int((lv & carrier).sum()),
                            "est_share_median": float(np.median(est[lv & carrier])),
                            "true_share_mean": true_mean[-1],
                            "true_share_uncapped_median": float(np.median(tbox[-1])),
                            "pct_detected_in_carriers": det[-1],
                            "pct_absent_in_noncarriers": absn[-1],
                            f"pct_within_{a.tol}": within[-1]})

        fig, axes = plt.subplots(1, 3, figsize=(14, 3.9))
        # 1 box plot: estimated (filled) next to real (grey) at each cap
        ax = axes[0]
        off, w = 0.19, 0.34
        common = dict(widths=w, whis=(5, 95), showfliers=False, patch_artist=True,
                      medianprops={"color": INK, "lw": 1.2})
        ax.boxplot(box, positions=x - off, boxprops={"facecolor": col, "alpha": 0.45, "edgecolor": col},
                   whiskerprops={"color": col}, capprops={"color": col}, **common)
        ax.boxplot(tbox, positions=x + off, boxprops={"facecolor": "#cfcec9", "edgecolor": MUTED},
                   whiskerprops={"color": MUTED}, capprops={"color": MUTED}, **common)
        ax.set_xlim(-0.6, len(caps) - 0.4)
        ax.set_ylim(-0.03, ymax + (0.03 if sig == "SBS1" else 0)); ax.set_ylabel("share per cell")
        ax.legend(handles=[Patch(facecolor=col, alpha=0.45, edgecolor=col, label="estimated"),
                           Patch(facecolor="#cfcec9", edgecolor=MUTED, label="real")],
                  frameon=False, fontsize=8.5, loc="lower left", ncol=2)
        ax.set_title("1  per-cell share (box 25–75%, whiskers 5–95%)", loc="left", fontsize=9.5, color=INK)
        # 2 correctly called
        ax = axes[1]
        ax.plot(x, det, color=col, lw=2.2, marker="o", ms=5, markeredgecolor="white",
                label=f"detected (≥{a.call:.0%}) where present")
        if grp:
            ax.plot(x, absn, color=col, lw=2, ls=":", marker="o", ms=5, markeredgecolor="white",
                    label=f"absent (<{a.call:.0%}) where not present")
        ax.set_ylim(0, 105); ax.set_ylabel("% of cells")
        ax.legend(frameon=False, fontsize=8.5, loc="lower left")
        ax.set_title("2  cells called correctly", loc="left", fontsize=9.5, color=INK)
        # 3 within tolerance
        ax = axes[2]
        ax.plot(x, within, color=col, lw=2.2, marker="o", ms=5, markeredgecolor="white")
        ax.set_ylim(0, 105); ax.set_ylabel("% of cells")
        ax.set_title(f"3  cells within ±{a.tol} of their true share", loc="left", fontsize=9.5, color=INK)
        for ax in axes:
            ax.set_xticks(x); ax.set_xticklabels(labels)
            ax.set_xlabel("mutations per cell (cap)")
            ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)
        cells = {"clone_SBS26": "SBS26-clone cells", "clone_SBS5": "SBS5-clone cells", None: "all cells"}[grp]
        fig.suptitle(f"{sig} in {cells} — {title}", x=0.01, ha="left", fontsize=11.5, color=INK)
        fig.tight_layout()
        for ext in ("png", "pdf"):
            fig.savefig(f"{OUT}/{tag}_{sig}_percell.{ext}", dpi=200)
        plt.close(fig)

pd.DataFrame(summary).round(4).to_csv(f"{OUT}/percell_values.csv", index=False)
print("wrote", sorted(f for f in os.listdir(OUT) if f.endswith(".png")))
