#!/usr/bin/env python3
"""
Fig 3a (signatures overall) over the missingness levels.
Same as Case_studies/Figure4_LSCC/Fig3a_signatures_overall.ipynb, run on every level:
  - per cell   : results/signatures/per_leaf/...   signatures in >= 10 % of cells
  - pseudobulk : results/signatures/all_leaves/... no threshold
Mutations per Mb (activity / 3200), signatures sorted by median of nonzero values,
one dot per sample, red line = median, bottom axis = number of samples with > 0.

Reads   <out>/<name>/pXX/results/signatures/{per_leaf,all_leaves}/Assignment_Solution/Activities/Assignment_Solution_Activities.txt
Writes  <out>/<name>/analysis/signatures/
          mutations_per_mb_per_cell_pXX.pdf       one per level (as the notebook)
          mutations_per_mb_all_pXX.pdf            pseudobulk, one per level
          mutations_per_mb_per_cell_all_levels.pdf / mutations_per_mb_all_all_levels.pdf   levels side by side
          signatures_per_cell_vs_missingness.csv / signatures_pseudobulk_vs_missingness.csv

  python signatures_overall.py --name neck
  python signatures_overall.py --name breast --levels 0 86
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'

HERE = Path(__file__).resolve().parent

# ---- PARAMETERS (from the notebook) ----
GENOME_SIZE_MB = 3200  # adjust depending on WGS/WES/panel
ANALYSES = {
    # name: (SigDyn folder, MIN_CELL_FRACTION, output stem)
    "per_cell": ("per_leaf", 0.10, "mutations_per_mb_per_cell"),   # 10% threshold
    "pseudobulk": ("all_leaves", 0, "mutations_per_mb_all"),
}
ACTIVITIES = "Assignment_Solution/Activities/Assignment_Solution_Activities.txt"


def load(path, min_cell_fraction):
    df = pd.read_csv(path, sep="\t").set_index("Samples")

    # Remove signatures with all zeros
    df = df.loc[:, (df != 0).any(axis=0)]

    # Remove signatures in <10% of cells
    cell_threshold = int(len(df) * min_cell_fraction)
    df = df.loc[:, (df > 0).sum(axis=0) >= cell_threshold]

    # Normalize to mutations per megabase
    df = df / GENOME_SIZE_MB

    # Sort signatures by nonzero median
    nonzero_medians = df.apply(lambda col: col[col > 0].median(), axis=0)
    return df[nonzero_medians.sort_values().index]


def plot(ax, df, ylabel=True):
    for i, sig in enumerate(df.columns, start=1):
        values = df[sig].values
        values = np.sort(values[values > 0])
        if len(values) > 0:
            x_positions = np.linspace(i - 0.4, i + 0.4, len(values))
            ax.scatter(x_positions, values, color="black", s=1, zorder=3)
            ax.hlines(np.median(values), i - 0.4, i + 0.4, colors="red", lw=1, zorder=4)

    # Add alternating grey background
    for i in range(1, len(df.columns) + 1, 2):
        ax.axvspan(i - 0.5, i + 0.5, facecolor="grey", alpha=0.3, zorder=0)

    # Y-axis
    ax.set_yscale("log")
    if ylabel:
        ax.set_ylabel("Somatic Mutations per Megabase", fontsize=25)
    ax.tick_params(axis="y", labelsize=20)
    ax.grid(axis="y", linestyle="--", alpha=0.6, zorder=0)

    # Bottom x-axis → number of samples with >0 mutations
    counts = (df > 0).sum(axis=0).values
    ax.set_xticks(range(1, len(df.columns) + 1))
    ax.set_xticklabels(counts, rotation=40, ha="center", fontsize=20)

    # Top x-axis → signature names
    ax_top = ax.secondary_xaxis("top")
    ax_top.set_xticks(range(1, len(df.columns) + 1))
    ax_top.set_xticklabels(df.columns, rotation=40, fontsize=20)

    # Remove first and last column space
    ax.set_xlim(0.5, len(df.columns) + 0.5)


def summary_rows(level, df):
    pos = df > 0
    return [{"level": level, "signature": sig, "n_samples_gt0": int(pos[sig].sum()),
             "frac_samples_gt0": pos[sig].mean(), "median_nonzero_mut_per_mb": df[sig][pos[sig]].median()}
            for sig in df.columns]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", default="neck")
    ap.add_argument("--levels", nargs="*", help="levels to plot, e.g. 0 30 (default: every finished level)")
    ap.add_argument("--out", default=str(HERE / "levels"))
    a = ap.parse_args()

    base = Path(a.out) / a.name
    dirs = sorted(p for p in base.glob("p*")
                  if all((p / "results/signatures" / f / ACTIVITIES).exists() for f, _, _ in ANALYSES.values()))
    if a.levels:
        want = {f"p{int(x.lstrip('p')):02d}" for x in a.levels}
        dirs = [d for d in dirs if d.name in want]
    if not dirs:
        raise SystemExit(f"no finished stage 2 under {base} - run run_levels.py stage2 first")

    outdir = base / "analysis" / "signatures"
    outdir.mkdir(parents=True, exist_ok=True)

    for analysis, (folder, min_frac, stem) in ANALYSES.items():
        dfs, rows = {}, []
        for d in dirs:
            df = load(d / "results/signatures" / folder / ACTIVITIES, min_frac)
            dfs[d.name] = df
            rows += summary_rows(d.name, df)

            fig, ax = plt.subplots(figsize=(7, 7), dpi=300)
            plot(ax, df)
            plt.tight_layout()
            fig.savefig(outdir / f"{stem}_{d.name}.pdf", format="pdf")
            plt.close(fig)
            print(f"{analysis} {d.name}: {len(df)} samples, {len(df.columns)} signatures: {', '.join(df.columns)}")

        # all levels side by side, same y-axis
        n = len(dfs)
        widths = [max(len(df.columns), 1) for df in dfs.values()]
        fig, axes = plt.subplots(1, n, figsize=(sum(7 * w / 6 for w in widths) + 1, 7), dpi=300,
                                 sharey=True, squeeze=False, gridspec_kw={"width_ratios": widths})
        for j, (ax, (lvl, df)) in enumerate(zip(axes[0], dfs.items())):
            plot(ax, df, ylabel=(j == 0))
            ax.annotate(lvl, xy=(0.5, 1), xycoords="axes fraction", xytext=(0, 75),
                        textcoords="offset points", ha="center", fontsize=25)
        plt.tight_layout()
        fig.savefig(outdir / f"{stem}_all_levels.pdf", format="pdf")
        plt.close(fig)

        pd.DataFrame(rows).to_csv(outdir / f"signatures_{analysis}_vs_missingness.csv", index=False)

    print(f"\nwritten to {outdir}")


if __name__ == "__main__":
    main()
