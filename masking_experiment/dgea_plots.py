#!/usr/bin/env python3
"""
Fig 4c/d over the missingness levels, step 3: volcano + GO dotplot per level, and a comparison
of every level against the 0 % level (= the paper's DGEA).

Everything is made twice, in a subfolder per p-value:
  fdr/   significance and volcano y-axis from the FDR (BH-adjusted; label "Adjusted p-value")
  rawp/  from MAST's raw Pr(>Chisq), as the paper's Fig4_c_d_volcano_plot.ipynb did

Per level (as Fig4_c_d_volcano_plot.ipynb), in <out>/<name>/dgea/pXX/{fdr,rawp}/:
  volcano.pdf/png        significant = p < 0.01 and |log2FC| > 0.5 (red: up in H32, blue: up in L32)
  significant_genes.csv
  go_enrichment.csv      Enrichr GO_Biological_Process_2021 on the genes up in L32 (log2FC < 0), as in the paper
  go_dotplot.pdf/png     top 5 terms by adjusted p
Across levels, in <out>/<name>/dgea/{fdr,rawp}/:
  dgea_vs_missingness.csv, go_top_terms_vs_missingness.csv, summary.pdf/png

  python dgea_plots.py --name neck
  python dgea_plots.py --name neck --no-go          # skip Enrichr (needs internet); uses saved go_enrichment.csv if present

Enrichr results are saved; delete go_enrichment.csv to query again.
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
LOGFC_THRESH, PVAL_THRESH = 0.5, 0.01          # as Fig4_c_d_volcano_plot.ipynb
GENE_SET, TOP_GO = "GO_Biological_Process_2021", 5
REF = "p00"
# subfolder: (column used, volcano y label, short name)
MODES = {"fdr": ("FDR", "-log10 Adjusted p-value", "FDR"),
         "rawp": ("pval", "-log10 p-value", "p")}

plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"]})


def load_mast(path, pcol):
    df = pd.read_csv(path)
    df.columns = ["gene", "pval", "logfc", "FDR"]          # pval = Pr(>Chisq) (raw); FDR = BH-adjusted (dgea.R)
    df["-log10(p)"] = -np.log10(df[pcol].clip(lower=1e-300))
    sig = df[pcol] < PVAL_THRESH
    df["direction"] = np.where(sig & (df["logfc"] > LOGFC_THRESH), "up_H32",
                      np.where(sig & (df["logfc"] < -LOGFC_THRESH), "up_L32", "ns"))
    df["significant"] = df["direction"] != "ns"
    return df


def style(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.4)
    ax.spines["left"].set_linewidth(0.4)


def volcano(df, title, ylabel, path):
    colors = df["direction"].map({"up_H32": "red", "up_L32": "blue", "ns": "lightgray"})
    fig, ax = plt.subplots(figsize=(10, 7))
    ax.scatter(df["logfc"], df["-log10(p)"], c=colors, s=20, edgecolor="none")
    ax.axhline(-np.log10(PVAL_THRESH), linestyle="--", color="gray", lw=0.4)
    ax.axvline(-LOGFC_THRESH, linestyle="--", color="gray", lw=0.4)
    ax.axvline(LOGFC_THRESH, linestyle="--", color="gray", lw=0.4)
    # top 5 genes each way by logFC (the paper placed these labels by hand)
    ymax = df["-log10(p)"].max()
    for side, top in ((1, df.nlargest(5, "logfc")), (-1, df.nsmallest(5, "logfc"))):
        for k, (_, r) in enumerate(top.iterrows()):
            ax.annotate(r["gene"], xy=(r["logfc"], r["-log10(p)"]),
                        xytext=(r["logfc"] + side * 0.25, ymax * (0.95 - 0.08 * k)),
                        fontsize=12, ha="left" if side > 0 else "right", va="center",
                        arrowprops=dict(arrowstyle="->", color="black", lw=0.4))
    ax.legend(handles=[Line2D([0], [0], marker="o", color="w", label=l, markerfacecolor=c, markersize=10)
                       for l, c in (("Upregulated in H32", "red"), ("Upregulated in L32", "blue"),
                                    ("Not significant", "lightgray"))],
              loc="lower right", bbox_to_anchor=(1, 0.2), frameon=False, fontsize=14)
    ax.set_title(title, fontsize=18)
    ax.set_xlabel("log2 Fold change", fontsize=18)
    ax.set_ylabel(ylabel, fontsize=18)
    lo, hi = df["logfc"].min(), df["logfc"].max()
    ax.set_xlim(min(-1.7, lo - 0.6), max(1.5, hi + 0.6))
    ax.tick_params(labelsize=16)
    style(ax)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(path.with_suffix("." + ext), dpi=300)
    plt.close(fig)


def go_enrichment(genes, path, query):
    if path.exists():
        return pd.read_csv(path)
    if not query:
        return None
    if not genes:
        res = pd.DataFrame(columns=["Term", "Overlap", "P-value", "Adjusted P-value", "Genes"])
    else:
        import gseapy as gp
        res = gp.enrichr(gene_list=genes, gene_sets=GENE_SET, organism="Human", outdir=None, cutoff=0.05).results
    res.to_csv(path, index=False)
    return res


def go_dotplot(res, n_genes, title, path):
    """As the dotplot in Fig4_c_d_volcano_plot.ipynb (incl. its 'Fold Enrichment' = gene count / term size)."""
    df = res.copy()
    if df.empty:
        return
    df["Gene Count"] = df["Overlap"].str.split("/").str[0].astype(int)
    df["GO Term Total Genes"] = df["Overlap"].str.split("/").str[1].astype(int)
    df["Fold Enrichment"] = (df["Gene Count"] / n_genes) / (df["GO Term Total Genes"] / n_genes)
    df["-log10(Adjusted P-value)"] = -np.log10(df["Adjusted P-value"] + 1e-10)
    df = df.sort_values("Adjusted P-value").head(TOP_GO)
    fig, ax = plt.subplots(figsize=(10, 6))
    sca = ax.scatter(df["-log10(Adjusted P-value)"], df["Term"], s=df["Gene Count"] * 10,
                     c=df["Fold Enrichment"], cmap="viridis", edgecolor="black", alpha=0.8)
    ax.set_xlabel("-log10(Adjusted P-value)", fontsize=16)
    ax.set_ylabel("GO Term", fontsize=16)
    ax.set_title(title, fontsize=16)
    ax.tick_params(labelsize=14)
    cb = fig.colorbar(sca)
    cb.set_label("Fold Enrichment", fontsize=16)
    ax.set_ylim(len(df) - 0.5, -0.5)
    style(ax)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(path.with_suffix("." + ext), dpi=300)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", default="neck")
    ap.add_argument("--out", default=str(HERE / "levels"))
    ap.add_argument("--no-go", action="store_true", help="do not query Enrichr")
    a = ap.parse_args()
    for mode in MODES:
        run(a, mode)


def run(a, mode):
    pcol, ylabel, pname = MODES[mode]
    base = Path(a.out) / a.name
    dg = base / "dgea"
    pct = {}
    if (base / "levels_summary.csv").exists():
        s = pd.read_csv(base / "levels_summary.csv")
        pct = {f"p{int(l):02d}": 100 * f for l, f in zip(s["level"], s["frac_missing"])}
    sizes = pd.read_csv(dg / "group_sizes.csv").set_index("level") if (dg / "group_sizes.csv").exists() else None

    levels = sorted(p.name for p in dg.glob("p*") if (p / "mast_results.csv").exists())
    if not levels:
        raise SystemExit(f"no mast_results.csv under {dg}/p*/ (run dgea_mast.R first)")
    print(f"[{mode}] levels:", " ".join(levels))
    for lv in levels:
        (dg / lv / mode).mkdir(exist_ok=True)

    mast, go = {}, {}
    for lv in levels:
        df = load_mast(dg / lv / "mast_results.csv", pcol)
        od = dg / lv / mode
        mast[lv] = df
        title = f"{int(lv[1:])} % missing (H32 vs L32)"
        volcano(df, title, ylabel, od / "volcano")
        df.loc[df["significant"], ["gene", "logfc", "pval", "FDR", "direction"]] \
          .sort_values(pcol).to_csv(od / "significant_genes.csv", index=False)
        genes_l32 = df.loc[df["direction"] == "up_L32", "gene"].tolist()
        go[lv] = go_enrichment(genes_l32, od / "go_enrichment.csv", not a.no_go)
        if go[lv] is not None:
            go_dotplot(go[lv], len(genes_l32), title, od / "go_dotplot")
        print(f"[{mode}] {lv}: {(df['direction'] == 'up_H32').sum()} up in H32, {len(genes_l32)} up in L32")

    # ---- comparison with the 0 % level ----
    ref = mast.get(REF)
    rows = []
    for lv, df in mast.items():
        r = {"level": lv, "pct_missing": pct.get(lv, np.nan)}
        if sizes is not None and lv in sizes.index:
            r["n_cells_H32"], r["n_cells_L32"] = sizes.loc[lv, "group2_H32"], sizes.loc[lv, "group1_L32"]
        r["n_genes_tested"] = len(df)
        for d in ("up_H32", "up_L32"):
            s = set(df.loc[df["direction"] == d, "gene"])
            r[f"n_{d}"] = len(s)
            r[f"n_{d}_FDR05"] = int(((df["FDR"] < 0.05) & (np.sign(df["logfc"]) == (1 if d == "up_H32" else -1))
                                     & (df["logfc"].abs() > LOGFC_THRESH)).sum())
            if ref is not None:
                s0 = set(ref.loc[ref["direction"] == d, "gene"])
                r[f"recall_{d}"] = len(s & s0) / len(s0) if s0 else np.nan        # fraction of 0 % genes found again
                r[f"precision_{d}"] = len(s & s0) / len(s) if s else np.nan     # fraction of this level's genes that are 0 % genes
                r[f"jaccard_{d}"] = len(s & s0) / len(s | s0) if s | s0 else np.nan
        if ref is not None:
            m = ref.merge(df, on="gene", suffixes=("_ref", ""))
            r["spearman_logfc_all"] = spearmanr(m["logfc_ref"], m["logfc"])[0]
            ms = m[m["significant_ref"]]
            r["spearman_logfc_ref_sig"] = spearmanr(ms["logfc_ref"], ms["logfc"])[0] if len(ms) > 2 else np.nan
            r["sign_agreement_ref_sig"] = (np.sign(ms["logfc_ref"]) == np.sign(ms["logfc"])).mean() if len(ms) else np.nan
        rows.append(r)
    summ = pd.DataFrame(rows)
    sd = dg / mode
    sd.mkdir(exist_ok=True)
    summ.to_csv(sd / "dgea_vs_missingness.csv", index=False)
    print(summ.T.to_string(header=False))

    # GO: the 0 % level's top terms, followed over the levels
    go_tab = None
    if go.get(REF) is not None and not go[REF].empty:
        top_terms = go[REF].sort_values("Adjusted P-value").head(10)["Term"].tolist()
        go_tab = pd.DataFrame({lv: (g.set_index("Term")["Adjusted P-value"].reindex(top_terms)
                                    if g is not None and not g.empty else pd.Series(np.nan, index=top_terms))
                               for lv, g in go.items()})
        go_tab.index.name = "Term (top 10 at 0 %)"
        go_tab.to_csv(sd / "go_top_terms_vs_missingness.csv")

    # ---- summary figure ----
    x = summ["pct_missing"].fillna(summ["level"].str[1:].astype(float)).values
    fig, axs = plt.subplots(2, 2, figsize=(12, 9))
    ax = axs[0, 0]
    ax.plot(x, summ["n_up_H32"], "o-", color="red", label="up in H32")
    ax.plot(x, summ["n_up_L32"], "o-", color="blue", label="up in L32")
    ax.set_ylabel("significant genes"); ax.set_title(f"DE genes ({pname} < 0.01, |log2FC| > 0.5)"); ax.legend(frameon=False)
    if ref is not None:
        ax = axs[0, 1]
        for d, c in (("up_H32", "red"), ("up_L32", "blue")):
            ax.plot(x, summ[f"recall_{d}"], "o-", color=c, label=f"recall, {d.replace('_', ' in ')}")
            ax.plot(x, summ[f"precision_{d}"], "s--", color=c, mfc="white", label=f"precision, {d.replace('_', ' in ')}")
        ax.plot(x, summ["spearman_logfc_all"], "o-", color="black", label="Spearman log2FC (all genes)")
        ax.set_ylim(-0.05, 1.05); ax.set_title("agreement with 0 %"); ax.legend(frameon=False, fontsize=8)
        ax = axs[1, 0]
        for lv, c in zip(levels, plt.cm.viridis(np.linspace(0, 0.9, len(levels)))):
            if lv == REF:
                continue
            m = ref.merge(mast[lv], on="gene", suffixes=("_ref", ""))
            ax.scatter(m["logfc_ref"], m["logfc"], s=3, color=c, alpha=0.5, label=f"{int(lv[1:])} %")
        lim = ax.get_xlim()
        ax.plot(lim, lim, color="gray", lw=0.5, ls="--")
        ax.set_xlabel("log2FC at 0 %"); ax.set_ylabel("log2FC at level"); ax.set_title("log2FC vs 0 %")
        ax.legend(frameon=False, markerscale=4)
    ax = axs[1, 1]
    if go_tab is not None:
        val = -np.log10(go_tab.astype(float).clip(lower=1e-300))
        im = ax.imshow(val.values, aspect="auto", cmap="viridis")
        ax.set_xticks(range(val.shape[1]), [f"{int(c[1:])} %" for c in val.columns])
        ax.set_yticks(range(val.shape[0]), [t if len(t) < 45 else t[:42] + "..." for t in val.index], fontsize=7)
        for i in range(val.shape[0]):
            for j in range(val.shape[1]):
                p = go_tab.iat[i, j]
                ax.text(j, i, "n/a" if pd.isna(p) else ("*" if p < 0.05 else ""), ha="center", va="center",
                        color="white", fontsize=8)
        fig.colorbar(im, ax=ax, label="-log10 adj. p")
        ax.set_title("GO BP (up in L32): top 10 terms at 0 %\n(* adj. p < 0.05)")
    else:
        ax.axis("off")
    for a_ in axs.flat:
        if a_.get_xlabel() == "" and a_ is not axs[1, 1]:
            a_.set_xlabel("% missing")
        style(a_)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(sd / f"summary.{ext}", dpi=200)
    plt.close(fig)
    print("wrote", sd)


if __name__ == "__main__":
    main()
