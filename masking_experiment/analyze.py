#!/usr/bin/env python3
"""
Collect results over the missingness levels from the SigDyn outputs in <out>/<name>/pXX/results/.

  always       : SH-like branch support of results/trees/cell_tree.unrooted
                 (as tree_support/unrooted_support_breast_neck2.ipynb)
  --entropy    : region entropy on results/trees/cell_tree.rooted, K=30
                 (as entropy_neck/region_entropy_scalable.ipynb; region = cell-name prefix, e.g. neck15)
  --signature  : % cells with that signature > 0, overall and per region, from
                 results/signatures/per_leaf (as Fig4_b_box_plots.ipynb)
  --mmr DIR    : Spearman/Pearson of MLH1/MSH2/MSH6/PMS2 expression vs SBS26, per cell
                 (as single_cell_exploration.ipynb; DIR = filtered_feature_bc_matrix)

Writes CSVs + summary.png/pdf to <out>/<name>/analysis/.

  python analyze.py --name neck --entropy --signature SBS32
  python analyze.py --name breast --mmr ../filtered_feature_bc_matrix
"""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import dendropy
from ete3 import Tree

HERE = Path(__file__).resolve().parent
LABELS = {"neck15": "Normal mucosa (N)", "neck17": "Cancer margin (R)",
          "neck18": "Carcinoma in situ (T)", "neck16": "Lymph-node met. (L)"}
COLORS = {"neck15": "#E69F00", "neck17": "#56B4E9", "neck18": "#D55E00", "neck16": "#CC79A7", "overall": "black"}
K, MAX_POOL, SEED, FLOOR = 30, 200, 0, 1e-6
MMR_GENES = ["MLH1", "MSH2", "MSH6", "PMS2"]
ACTIVITIES = "results/signatures/per_leaf/Assignment_Solution/Activities/Assignment_Solution_Activities.txt"


def leaf_to_region(name):
    return name.split("_")[0]


# ---------- entropy (copied from region_entropy_scalable.ipynb) ----------
def local_neighborhoods(tree, k, max_pool=200, seed=0):
    rng = np.random.default_rng(seed)
    lvs = tree.get_leaves()
    idx = {lf: i for i, lf in enumerate(lvs)}
    n = len(lvs)
    for node in tree.traverse("postorder"):
        node.add_feature("nleaf", 1 if node.is_leaf() else sum(c.nleaf for c in node.children))
    cache, members = {}, []
    indptr = np.empty(n + 1, dtype=np.int64); indptr[0] = 0
    for i, lf in enumerate(lvs):
        a = lf
        while a.nleaf < k + 1 and a.up is not None:
            a = a.up
        pool = cache.get(id(a))
        if pool is None:
            pool = np.fromiter((idx[x] for x in a.get_leaves()), dtype=np.int32)
            if pool.size > max_pool:
                pool = rng.choice(pool, max_pool, replace=False)
            cache[id(a)] = pool
        nb = pool[pool != i]
        members.append(nb)
        indptr[i + 1] = indptr[i] + nb.size
    return np.concatenate(members), indptr, [lf.name for lf in lvs]


def nbhd_entropy_bits(members, indptr, y_int, M):
    lab = y_int[members]
    sizes = np.diff(indptr).astype(float)
    counts = np.empty((len(sizes), M))
    for c in range(M):
        counts[:, c] = np.add.reduceat((lab == c).astype(np.int64), indptr[:-1])
    p = counts / sizes[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(p > 0, -p * np.log2(p), 0.0).sum(1)


def region_entropy(tree_file):
    t = Tree(str(tree_file), format=1)
    members, indptr, order = local_neighborhoods(t, K, MAX_POOL, SEED)
    regions = np.array([leaf_to_region(n) for n in order])
    region_order = sorted(np.unique(regions))
    code = {r: i for i, r in enumerate(region_order)}
    H = nbhd_entropy_bits(members, indptr, np.array([code[r] for r in regions]), len(region_order))
    Hmax = np.log2(len(region_order))
    return {str(r): float(H[regions == r].mean() / Hmax) for r in region_order}


# ---------- tree support (copied from unrooted_support_breast_neck2.ipynb) ----------
def support_of(node):
    try:
        return float(node.label)
    except (TypeError, ValueError):
        return None


def tree_support(tree_file):
    t = dendropy.Tree.get(path=str(tree_file), schema="newick", preserve_underscores=True)
    t.encode_bipartitions()
    N = len(t.leaf_nodes())
    sup, bl, minside = [], [], []
    for e in t.preorder_edge_iter():
        n = e.head_node
        if n is None or n.is_leaf() or e.bipartition is None:
            continue
        s = support_of(n)
        if s is None:
            continue
        sup.append(s)
        bl.append(e.length if e.length is not None else np.nan)
        k = bin(e.bipartition.leafset_bitmask).count("1")
        minside.append(min(k, N - k))
    sup, bl, minside = np.array(sup), np.array(bl, float), np.array(minside)
    all_bl = np.array([e.length for e in t.preorder_edge_iter() if e.length is not None])
    m = np.isfinite(bl)
    return {"median support": np.median(sup), "mean support": sup.mean(),
            "% support < 0.5": 100 * (sup < 0.5).mean(),
            "% support >= 0.7": 100 * (sup >= 0.7).mean(),
            "% support >= 0.9": 100 * (sup >= 0.9).mean(),
            "% branches at ~0 floor": 100 * (all_bl < FLOOR).mean(),
            "corr(support, branch length)": np.corrcoef(sup[m], bl[m])[0, 1],
            "corr(support, split balance)": np.corrcoef(sup, minside)[0, 1]}


# ---------- signature presence (as Fig4_b_box_plots.ipynb) ----------
def signature_presence(act, sig):
    s = act[sig] if sig in act else pd.Series(0.0, index=act.index)
    reg = act.index.map(leaf_to_region)
    rows = [("overall", s)] + [(r, s[reg == r]) for r in sorted(reg.unique())]
    return [{"region": r, "pct_cells": 100 * (x > 0).mean(), "mean_activity": x.mean(), "n_cells": len(x)}
            for r, x in rows]


# ---------- MMR vs SBS26 (as single_cell_exploration.ipynb) ----------
def mmr_expression(matrix_dir):
    """Log-normalized expression (normalize_total 1e4 + log1p). The notebook also ran sc.pp.scale,
    which is linear per gene and does not change Spearman/Pearson."""
    import scanpy as sc
    from scipy.io import mmread
    d = Path(matrix_dir)
    X = mmread(d / "matrix.mtx").T.tocsr()
    genes = pd.read_csv(d / "features.tsv", header=None, sep="\t")[1].values
    barcodes = pd.read_csv(d / "barcodes.tsv", header=None, sep="\t")[0].values
    ad = sc.AnnData(X=X.astype(float), var=pd.DataFrame(index=genes), obs=pd.DataFrame(index=barcodes))
    ad.var_names_make_unique()
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)
    return pd.DataFrame(ad[:, MMR_GENES].X.toarray(), index=ad.obs_names, columns=MMR_GENES)


def mmr_sbs26_corr(act, expr):
    from scipy.stats import pearsonr, spearmanr
    s = act["SBS26"] if "SBS26" in act else pd.Series(0.0, index=act.index)
    bc = s.index.str.rsplit("_", n=1).str[-1]          # cell name -> barcode (= name[20:] for breast)
    keep = bc.isin(expr.index)
    s, bc = s[keep].astype(float), bc[keep]
    rows = []
    for g in MMR_GENES:
        e = expr.loc[bc, g].values
        r_s, p_s = spearmanr(e, s.values) if s.std() > 0 else (np.nan, np.nan)
        r_p, p_p = pearsonr(e, s.values) if s.std() > 0 else (np.nan, np.nan)
        rows.append({"gene": g, "spearman_r": r_s, "spearman_p": p_s, "pearson_r": r_p, "pearson_p": p_p,
                     "n_cells": int(keep.sum()), "pct_cells_SBS26": 100 * (s > 0).mean()})
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", required=True)
    ap.add_argument("--entropy", action="store_true", help="region entropy on the rooted tree")
    ap.add_argument("--signature", help="e.g. SBS32: % cells with activity > 0, per region")
    ap.add_argument("--mmr", metavar="MATRIX_DIR", help="filtered_feature_bc_matrix for MMR vs SBS26")
    ap.add_argument("--out", default=str(HERE / "levels"))
    a = ap.parse_args()

    base = Path(a.out) / a.name
    res = base / "analysis"; res.mkdir(exist_ok=True)
    levels = pd.read_csv(base / "levels_summary.csv").set_index("level")["frac_missing"]
    expr = mmr_expression(a.mmr) if a.mmr else None

    tables = {"tree_support": [], "entropy": [], f"{a.signature}_presence": [], "MMR_SBS26_corr": []}
    for lvl, frac in levels.items():
        run, info = base / f"p{lvl:02d}", {"level": lvl, "pct_missing": round(100 * frac, 1)}
        if (run / "results/trees/cell_tree.unrooted").exists():
            tables["tree_support"].append({**info, **tree_support(run / "results/trees/cell_tree.unrooted")})
        if a.entropy and (run / "results/trees/cell_tree.rooted").exists():
            for r, e in region_entropy(run / "results/trees/cell_tree.rooted").items():
                tables["entropy"].append({**info, "region": r, "label": LABELS.get(r, r), "entropy_norm": e})
        if (run / ACTIVITIES).exists() and (a.signature or a.mmr):
            act = pd.read_csv(run / ACTIVITIES, sep="\t", index_col=0)
            if a.signature:
                for row in signature_presence(act, a.signature):
                    tables[f"{a.signature}_presence"].append({**info, **row, "label": LABELS.get(row["region"], row["region"])})
            if a.mmr:
                for row in mmr_sbs26_corr(act, expr):
                    tables["MMR_SBS26_corr"].append({**info, **row})

    tables = {k: pd.DataFrame(v) for k, v in tables.items() if v}
    for k, df in tables.items():
        df.round(4).to_csv(res / f"{k}_vs_missingness.csv", index=False)
        print(f"\n== {a.name}: {k}\n{df.round(3).to_string(index=False)}")

    # ---------- one figure, one panel per analysis ----------
    fig, axes = plt.subplots(1, len(tables), figsize=(5 * len(tables), 4), squeeze=False)
    for ax, (k, df) in zip(axes[0], tables.items()):
        if k == "tree_support":
            ax.plot(df.pct_missing, df["median support"], "o-", label="median support")
            ax.plot(df.pct_missing, df["% support >= 0.7"] / 100, "s--", label="fraction ≥ 0.7")
            ax.plot(df.pct_missing, df["% branches at ~0 floor"] / 100, "^:", label="fraction at ~0 length")
            ax.set(ylim=(0, 1.05), ylabel="SH-like support / fraction", title="Tree support")
        elif k == "entropy":
            for r, d in df.groupby("region"):
                ax.plot(d.pct_missing, d.entropy_norm, "o-", color=COLORS.get(r), label=LABELS.get(r, r))
            ax.set(ylim=(0, 1.05), ylabel="normalized entropy (0 clustered, 1 mixed)", title=f"Region mixing on tree (K={K})")
        elif k.endswith("_presence"):
            for r, d in df.groupby("region"):
                ax.plot(d.pct_missing, d.pct_cells, "o-", color=COLORS.get(r), lw=3 if r == "overall" else 1.5,
                        label=LABELS.get(r, r))
            ax.set(ylim=(0, 100), ylabel=f"% cells with {a.signature} > 0", title=f"{a.signature} detection")
        else:
            for g, d in df.groupby("gene"):
                ax.plot(d.pct_missing, d.spearman_r, "o-", label=g)
            ax.axhline(0, color="grey", lw=0.8)
            ax.set(ylabel="Spearman r (expression vs SBS26)", title="MMR genes vs SBS26")
        ax.set_xlabel("% missing entries")
        ax.legend(fontsize=8, frameon=False)
    fig.suptitle(a.name)
    fig.tight_layout()
    fig.savefig(res / "summary.png", dpi=200)
    fig.savefig(res / "summary.pdf")
    print(f"\nsaved -> {res}")


if __name__ == "__main__":
    main()
