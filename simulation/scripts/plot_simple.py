"""
plot_simple.py
One simple plot per signature (SBS26, SBS5, SBS1) and per clone arrangement of the
biopsy3 experiment: mean share of exposure vs mutations per cell.
  solid  = SigDyn estimate (mean of replicates; uncapped = 1 run)
  dashed = true share
Cells averaged: SBS26 -> SBS26-clone cells, SBS5 -> SBS5-clone cells, SBS1 -> all cells.
Estimated share = exposure / total exposure of the cell (all signatures).

  python scripts/plot_simple.py        # -> <output>/figures/simple/*.png + .pdf + values csv
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
os.makedirs(OUT, exist_ok=True)

SIGS = {  # signature: (cells to average, colour, y max)
    "SBS26": ("clone_SBS26", "#2a78d6", 0.35),
    "SBS5": ("clone_SBS5", "#eb6834", 0.35),
    "SBS1": (None, "#4a3aa7", 1.0),          # None = all cells
}
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
        for sig, (grp, _, _) in SIGS.items():
            m = (s.true_clone == grp) if grp else pd.Series(True, index=s.index)
            rows.append({"arrangement": tag, "signature": sig, "cap": int(cap), "rep": int(rep),
                         "median_muts": int(tot.median()),
                         "estimated": F.loc[m, sig].mean() if sig in F else 0.0,
                         "true": (s.loc[m, f"true_{sig}"] / tot[m]).mean()})
vals = pd.DataFrame(rows)
mean = vals.groupby(["arrangement", "signature", "cap"])[["estimated", "true", "median_muts"]].mean().reset_index()
mean.to_csv(f"{OUT}/simple_values.csv", index=False)

for X, tag, title in RUNS:
    for sig, (grp, col, ymax) in SIGS.items():
        d = mean[(mean.arrangement == tag) & (mean.signature == sig)].copy()
        d["order"] = d.cap.replace(0, 10**6)
        d = d.sort_values("order", ascending=False)
        x = np.arange(len(d))
        labels = [f"all\n(~{int(d.median_muts.iloc[0])})"] + [str(c) for c in d.cap.iloc[1:]]
        fig, ax = plt.subplots(figsize=(5.2, 3.8))
        ax.plot(x, d["true"], color=col, lw=1.6, ls="--", label="true", zorder=2)
        ax.plot(x, d["estimated"], color=col, lw=2.2, marker="o", ms=5,
                markeredgecolor="white", markeredgewidth=1.2, label="SigDyn estimate", zorder=3)
        ax.set_xticks(x); ax.set_xticklabels(labels)
        ax.set_xlabel("mutations per cell (cap)")
        ax.set_ylabel("mean share of exposure in cell")
        ax.set_ylim(0, ymax)
        ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)
        cells = {"clone_SBS26": "SBS26-clone cells", "clone_SBS5": "SBS5-clone cells", None: "all cells"}[grp]
        ax.set_title(f"{sig} in {cells}\n{title}", loc="left", fontsize=10.5, color=INK)
        ax.legend(frameon=False, fontsize=9, loc="lower left")
        fig.tight_layout()
        for ext in ("png", "pdf"):
            fig.savefig(f"{OUT}/{tag}_{sig}.{ext}", dpi=200)
        plt.close(fig)
print("wrote", sorted(f for f in os.listdir(OUT) if f.endswith(".png")))
