"""
plot_downsampling.py
Figures for the per-cell downsampling experiments (leaf exposures vs mutations per cell).

Fig 1  signature recovery:  mean share of each signature in its true clone's cells, per
       per-cell cap (mean of replicates, whiskers = replicate min-max), one panel per
       clone arrangement. Dashed blue / orange = true SBS26 / SBS5 share. Shaded = SBS26 not called by SigDyn.
Fig 2  accuracy:  per-cell cosine (true vs estimated exposures, all signatures counted,
       SBS12 counted as SBS26) median and 10th percentile, and share of cells assigned
       to the right group (SBS5 clone / SBS26 clone / SBS1-only).
Fig 4  SBS1: mean SBS1 share in each true group (SBS1-only, SBS5 clone, SBS26 clone),
       dashed = true share of that group.

  python scripts/plot_downsampling.py            # writes <output>/figures/*.png + .pdf
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
ap.add_argument("--runs", nargs="+", default=[
    os.environ.get("RUN_NAME", "sbs26_large_clone") + "|SBS26 = larger clone (369 cells)"])
ap.add_argument("--out", default="figures")
a = ap.parse_args()
OUT = os.path.join(M, a.out)
os.makedirs(OUT, exist_ok=True)

# reference palette (dataviz skill, categorical slots 1-4, light mode) + text tokens
C = {"SBS26": "#2a78d6", "SBS5": "#eb6834", "SBS12": "#1baf7a", "other": "#eda100", "SBS1only": "#4a3aa7"}
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e0"
plt.rcParams.update({"font.family": "sans-serif", "font.size": 10, "axes.edgecolor": MUTED,
                     "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "figure.facecolor": "white"})


def per_run_shares(X):
    """one row per run: mean shares by true group + whether SBS26 was called."""
    pc = pd.read_csv(f"{M}/{X}/results_merged_SBS12/leaf_per_cell.csv")
    met = pd.read_csv(f"{M}/{X}/results_merged_SBS12/leaf_metrics.csv")
    rows = []
    for _, r in met.iterrows():
        cap, rep = int(r.cap_per_cell), int(r.rep)
        d = f"{M}/{X}/runs/cap{'all' if cap == 0 else f'{cap:03d}'}_r{rep}"
        E = pd.read_csv(f"{d}/results/signatures/exposures_active_signatures.txt", sep="\t", index_col=0).T
        s = pc[(pc.mask_pct == cap) & (pc.rep == rep)].set_index("cell")
        F = E.loc[s.index]; F = F.div(F.sum(1), axis=0)
        g = lambda grp, c: F[s.true_clone == grp][c].mean() if c in F else 0.0
        tot = s[["true_SBS1", "true_SBS5", "true_SBS26"]].sum(1)
        sim = [c for c in ("SBS1", "SBS5", "SBS26") if c in F]
        rows.append({
            "cap": cap, "rep": rep, "muts": int(tot.median()),
            "SBS26 in SBS26 clone": g("clone_SBS26", "SBS26"),
            "SBS12 in SBS26 clone": g("clone_SBS26", "SBS12"),
            "SBS5 in SBS5 clone": g("clone_SBS5", "SBS5"),
            "SBS5 in SBS26 clone (wrong)": g("clone_SBS26", "SBS5"),
            "SBS5 in SBS1-only (false)": g("SBS1-only", "SBS5"),
            "non-simulated (all cells)": (1 - F[sim].sum(1) - (F["SBS12"] if "SBS12" in F else 0)).mean(),
            "SBS1 in SBS1-only": g("SBS1-only", "SBS1"),
            "SBS1 in SBS5 clone": g("clone_SBS5", "SBS1"),
            "SBS1 in SBS26 clone": g("clone_SBS26", "SBS1"),
            "true SBS1 in SBS1-only": (s.true_SBS1 / tot)[s.true_clone == "SBS1-only"].mean(),
            "true SBS1 in SBS5 clone": (s.true_SBS1 / tot)[s.true_clone == "clone_SBS5"].mean(),
            "true SBS1 in SBS26 clone": (s.true_SBS1 / tot)[s.true_clone == "clone_SBS26"].mean(),
            "true SBS26": (s.true_SBS26 / tot)[s.true_clone == "clone_SBS26"].mean(),
            "true SBS5": (s.true_SBS5 / tot)[s.true_clone == "clone_SBS5"].mean(),
            "SBS26 called": "SBS26" in r.active_signatures.split(";"),
            "cos median": r.cosine_median, "cos p10": r.cosine_p10,
            "clone accuracy": r.clone_assignment_accuracy,
        })
    df = pd.DataFrame(rows)
    # drop caps above every cell's burden (identical to the uncapped run)
    maxb = df.loc[df.cap == 0, "muts"].max()
    full_max = pc[pc.mask_pct == 0][["true_SBS1", "true_SBS5", "true_SBS26"]].sum(1).max()
    df = df[(df.cap == 0) | (df.cap < full_max)]
    df["x"] = np.where(df.cap == 0, df.muts, df.cap)
    return df


def summarise(df, cols):
    g = df.groupby("x")[cols]
    return g.mean(), g.min(), g.max()


def line(ax, x, mean, lo, hi, color, label, ls="-", label_at="right", dy=0):
    x = np.arange(len(mean))
    ax.plot(x, mean, color=color, lw=2, ls=ls, marker="o", ms=5,
            markeredgecolor="white", markeredgewidth=1.2, zorder=3)
    ax.errorbar(x, mean, yerr=[np.clip(mean - lo, 0, None), np.clip(hi - mean, 0, None)], fmt="none", ecolor=color,
                elinewidth=1, capsize=0, zorder=2)
    i = -1 if label_at == "left" else 0          # "left" = label after the last point
    ax.annotate(label, (x[i], mean[i]), xytext=(7 if i == -1 else 0, dy), textcoords="offset points",
                va="center", ha="left" if i == -1 else "center", fontsize=8.5, color=INK2)


def style(ax, xs, ylabel=None):
    ax.set_xticks(range(len(xs)))
    ax.set_xticklabels([f"all\n(~{int(xs[0])})"] + [str(int(v)) for v in xs[1:]])
    ax.set_xlim(-0.4, len(xs) - 1 + 2.6)            # room on the right for line labels
    ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)
    ax.set_xlabel("mutations per cell (cap)")
    if ylabel:
        ax.set_ylabel(ylabel)


def shade_not_called(ax, df, ytop):
    called = df.groupby("x")["SBS26 called"].mean()
    xs = sorted(called.index, reverse=True)
    miss = [i for i, v in enumerate(xs) if called[v] == 0]
    if miss:
        ax.axvspan(min(miss) - 0.5, max(miss) + 0.5, color="#f3f2ee", zorder=0)
        ax.text((min(miss) + max(miss)) / 2, ytop * 0.97, "SBS26 not called by SigDyn",
                ha="center", va="top", fontsize=8.5, color=MUTED)


runs = [(r.split("|")[0], r.split("|")[1]) for r in a.runs]
data = [(X, title, per_run_shares(X)) for X, title in runs]

# ---------------- Fig 1: signature shares ----------------
fig, axes = plt.subplots(1, len(data), figsize=(6.2 * len(data), 4.6), sharey=True)
for ax, (X, title, df) in zip(np.atleast_1d(axes), data):
    cols = ["SBS26 in SBS26 clone", "SBS5 in SBS5 clone", "SBS5 in SBS26 clone (wrong)",
            "non-simulated (all cells)", "true SBS26", "true SBS5"]
    mean, lo, hi = summarise(df, cols)
    xs = mean.index.values[::-1]
    mean, lo, hi = mean.loc[xs], lo.loc[xs], hi.loc[xs]
    shade_not_called(ax, df, 0.36)
    for tcol, col in (("true SBS26", C["SBS26"]), ("true SBS5", C["SBS5"])):
        ax.axhline(mean[tcol].mean(), color=col, ls="--", lw=1, alpha=0.6, zorder=1)
    ax.text(0, 0.335, f"dashed = true share\nSBS26 {mean['true SBS26'].mean():.2f} (blue), "
            f"SBS5 {mean['true SBS5'].mean():.2f} (orange)", fontsize=8, color=MUTED, ha="left", va="top")
    line(ax, xs, mean["SBS26 in SBS26 clone"].values, lo["SBS26 in SBS26 clone"].values,
         hi["SBS26 in SBS26 clone"].values, C["SBS26"], "SBS26 (SBS26 clone)", label_at="left")
    line(ax, xs, mean["SBS5 in SBS5 clone"].values, lo["SBS5 in SBS5 clone"].values,
         hi["SBS5 in SBS5 clone"].values, C["SBS5"], "SBS5 (SBS5 clone)", label_at="left")
    line(ax, xs, mean["SBS5 in SBS26 clone (wrong)"].values, lo["SBS5 in SBS26 clone (wrong)"].values,
         hi["SBS5 in SBS26 clone (wrong)"].values, C["SBS5"], "SBS5 in SBS26 clone (wrong)", ls=":", label_at="left")
    line(ax, xs, mean["non-simulated (all cells)"].values, lo["non-simulated (all cells)"].values,
         hi["non-simulated (all cells)"].values, C["other"], "non-simulated signatures", label_at="left")
    style(ax, xs, "mean share of exposure in cells" if ax is np.atleast_1d(axes)[0] else None)
    ax.set_ylim(0, 0.36)
    ax.set_title(title, loc="left", fontsize=11, color=INK)
fig.suptitle("Signature recovery as cells lose mutations\n(3 replicates, whiskers = min–max)",
             x=0.01, ha="left", fontsize=12, color=INK)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}/fig1_signature_shares.{ext}", dpi=200)

# ---------------- Fig 2: accuracy ----------------
fig, axes = plt.subplots(1, len(data), figsize=(6.2 * len(data), 4.2), sharey=True)
for ax, (X, title, df) in zip(np.atleast_1d(axes), data):
    mean, lo, hi = summarise(df, ["cos median", "cos p10", "clone accuracy"])
    xs = mean.index.values[::-1]
    mean, lo, hi = mean.loc[xs], lo.loc[xs], hi.loc[xs]
    shade_not_called(ax, df, 1.17)
    line(ax, xs, mean["cos median"].values, lo["cos median"].values, hi["cos median"].values,
         C["SBS26"], "cosine, median cell", label_at="left")
    line(ax, xs, mean["cos p10"].values, lo["cos p10"].values, hi["cos p10"].values,
         C["SBS26"], "cosine, 10th percentile cell", ls=":", label_at="left")
    line(ax, xs, mean["clone accuracy"].values, lo["clone accuracy"].values, hi["clone accuracy"].values,
         C["SBS5"], "cells in correct group", label_at="left")
    style(ax, xs, "value" if ax is np.atleast_1d(axes)[0] else None)
    ax.set_ylim(0, 1.17)
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.set_title(title, loc="left", fontsize=11, color=INK)
fig.suptitle("Per-cell accuracy\n(true vs estimated exposures; SBS12 counted as SBS26)",
             x=0.01, ha="left", fontsize=12, color=INK)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}/fig2_accuracy.{ext}", dpi=200)


# ---------------- Fig 4: SBS1 per true group ----------------
fig, axes = plt.subplots(1, len(data), figsize=(6.2 * len(data), 4.4), sharey=True)
groups = [("SBS1 in SBS1-only", "true SBS1 in SBS1-only", C["SBS1only"], "SBS1-only cells"),
          ("SBS1 in SBS5 clone", "true SBS1 in SBS5 clone", C["SBS5"], "SBS5 clone"),
          ("SBS1 in SBS26 clone", "true SBS1 in SBS26 clone", C["SBS26"], "SBS26 clone")]
for ax, (X, title, df) in zip(np.atleast_1d(axes), data):
    cols = [c for g_ in groups for c in g_[:2]]
    mean, lo, hi = summarise(df, cols)
    xs = mean.index.values[::-1]
    mean, lo, hi = mean.loc[xs], lo.loc[xs], hi.loc[xs]
    shade_not_called(ax, df, 1.12)
    for est, true, col, lab in groups:
        ax.axhline(mean[true].mean(), color=col, ls="--", lw=1, alpha=0.6, zorder=1)
    ax.text(0, 1.005, "dashed = true SBS1 share of each group", fontsize=8, color=MUTED, ha="left", va="bottom")
    for est, true, col, lab in groups:
        line(ax, xs, mean[est].values, lo[est].values, hi[est].values, col, "", label_at="left")
    # end labels, pushed apart so crossing lines don't overlap (min gap in data units)
    ends = sorted(((mean[est].values[-1], lab) for est, _, _, lab in groups), reverse=True)
    ypos, gap = [], 0.035
    for v, lab in ends:
        y = v if not ypos or ypos[-1] - v >= gap else ypos[-1] - gap
        ypos.append(y)
        ax.annotate(lab, (len(xs) - 1, v), xytext=(len(xs) - 1 + 0.12, y), textcoords="data",
                    va="center", ha="left", fontsize=8.5, color=INK2)
    style(ax, xs, "mean SBS1 share of exposure" if ax is np.atleast_1d(axes)[0] else None)
    ax.set_ylim(0.5, 1.12)
    ax.set_yticks(np.arange(0.5, 1.01, 0.1))
    ax.set_title(title, loc="left", fontsize=11, color=INK)
fig.suptitle("SBS1 recovery by true group\n(3 replicates, whiskers = min–max)",
             x=0.01, ha="left", fontsize=12, color=INK)
fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(f"{OUT}/fig4_sbs1.{ext}", dpi=200)

# table view of everything plotted
pd.concat([d.assign(run=X) for X, _, d in data]).to_csv(f"{OUT}/plotted_values.csv", index=False)
print("wrote", sorted(os.listdir(OUT)))
