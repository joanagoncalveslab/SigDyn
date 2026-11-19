"""
plot_simple_nonsimulated.py
Companion to plot_simple.py: the signatures SigDyn called that were NEVER simulated
(anything other than SBS1, SBS5, SBS26), each as its own line, vs mutations per cell.
One figure per clone arrangement of the biopsy3 experiment, one panel per true cell
group (SBS26-clone, SBS5-clone, SBS1-only cells). Their true share is 0 everywhere,
so every point above 0 is error.
  line   = mean share of exposure in the cell (mean of replicates; uncapped = 1 run)
  grey   = all non-simulated signatures together
Estimated share = exposure / total exposure of the cell (all signatures), as in plot_simple.py.

  python scripts/plot_simple_nonsimulated.py   # -> <output>/figures/simple/*_non_simulated.png/.pdf + csv
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

M = os.environ.get("SIM_OUT", os.getcwd())   # pipeline output folder (set by main.sh)
RUNS = [(os.environ.get("RUN_NAME", "sbs26_large_clone"), "large", "SBS26 = larger clone (369 cells)")]
OUT = os.path.join(M, "figures/simple")
SIMULATED = {"SBS1", "SBS5", "SBS26"}
GROUPS = [("clone_SBS26", "SBS26-clone cells"), ("clone_SBS5", "SBS5-clone cells"),
          ("SBS1-only", "SBS1-only cells")]
PALETTE = ["#2e9e6b", "#c9356f", "#d4a017", "#6b4fbb", "#1f9bb4", "#9a6b3f"]
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e0"
plt.rcParams.update({"font.family": "sans-serif", "font.size": 10, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False, "pdf.fonttype": 42})

rows = []
for X, tag, title in RUNS:
    pc = pd.read_csv(f"{M}/{X}/results_merged_SBS12/leaf_per_cell.csv")
    full_max = pc[pc.mask_pct == 0][["true_SBS1", "true_SBS5", "true_SBS26"]].sum(1).max()
    for (cap, rep), s in pc.groupby(["mask_pct", "rep"]):
        if cap and cap >= full_max:          # cap above every cell's burden = same as uncapped
            continue
        d = f"{M}/{X}/runs/cap{'all' if cap == 0 else f'{cap:03d}'}_r{rep}"
        E = pd.read_csv(f"{d}/results/signatures/exposures_active_signatures.txt", sep="\t", index_col=0).T
        s = s.set_index("cell")
        F = E.loc[s.index]; F = F.div(F.sum(1), axis=0)
        tot = s[["true_SBS1", "true_SBS5", "true_SBS26"]].sum(1)
        extra = [c for c in F.columns if c not in SIMULATED]
        for grp, _ in GROUPS:
            m = s.true_clone == grp
            row = {"arrangement": tag, "cell_group": grp, "cap": int(cap), "rep": int(rep),
                   "median_muts": int(tot.median()),
                   "all_non_simulated": F.loc[m, extra].sum(1).mean() if extra else 0.0}
            for sig in extra:
                row[sig] = F.loc[m, sig].mean()
            rows.append(row)
vals = pd.DataFrame(rows).fillna(0.0)
sigs = sorted(c for c in vals.columns if c.startswith("SBS"))
mean = vals.groupby(["arrangement", "cell_group", "cap"])[sigs + ["all_non_simulated", "median_muts"]].mean().reset_index()
mean.to_csv(f"{OUT}/non_simulated_values.csv", index=False)
colour = {s: PALETTE[i % len(PALETTE)] for i, s in enumerate(sigs)}

for X, tag, title in RUNS:
    d0 = mean[mean.arrangement == tag]
    shown = [s for s in sigs if d0[s].max() > 0.001]       # signatures that appear in this arrangement
    ymax = max(0.05, float(np.ceil(d0["all_non_simulated"].max() * 20) / 20) + 0.05)
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.9), sharey=True)
    for ax, (grp, glab) in zip(axes, GROUPS):
        d = d0[d0.cell_group == grp].copy()
        d["order"] = d.cap.replace(0, 10**6)
        d = d.sort_values("order", ascending=False)
        x = np.arange(len(d))
        labels = [f"all\n(~{int(d.median_muts.iloc[0])})"] + [str(c) for c in d.cap.iloc[1:]]
        ax.plot(x, d["all_non_simulated"], color=MUTED, lw=4, alpha=0.35, label="all non-simulated", zorder=1)
        for sig in shown:
            ax.plot(x, d[sig], color=colour[sig], lw=2, marker="o", ms=4.5,
                    markeredgecolor="white", markeredgewidth=1, label=sig, zorder=3)
        ax.axhline(0, color=INK2, lw=1, ls="--", label="true (0)")
        ax.set_xticks(x); ax.set_xticklabels(labels)
        ax.set_xlabel("mutations per cell (cap)")
        ax.set_ylim(-0.005, ymax)
        ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)
        ax.set_title(glab, loc="left", fontsize=10, color=INK)
    axes[0].set_ylabel("mean share of exposure in cell")
    axes[0].legend(frameon=False, fontsize=8.5, loc="upper left")
    fig.suptitle(f"Signatures called but never simulated — {title}", x=0.01, ha="left",
                 fontsize=11.5, color=INK)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{OUT}/{tag}_non_simulated.{ext}", dpi=200)
    plt.close(fig)
print("wrote", sorted(f for f in os.listdir(OUT) if "non_simulated" in f))
