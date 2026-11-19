#!/usr/bin/env python3
"""
Fig 4c/d over the missingness levels, step 1: expression matrix + SBS32 groups per level.

As bef_dgea.ipynb (Case_studies/Figure4_LSCC), except that the groups are
SBS32 = 0 (group1, L32) vs SBS32 > 0 (group2, H32), from each level's per-leaf activities.
Cells in the 10x data but not in the tree are "other" (as in the paper).

Expression does not depend on the level, so it is built once:
  <out>/<name>/dgea/expression.h5ad        normalize_total(1e4) + log1p, genes in >= 10 % of cells
  <out>/<name>/dgea/pXX/groups.csv         cell, dge_group (group1 / group2 / other)
  <out>/<name>/dgea/group_sizes.csv

  python dgea_prepare.py --name neck --tenx /path/with/neck15..neck18
  python dgea_prepare.py --name neck --tenx ... --cells 2000 --out test     # subsample test

--tenx must contain neck15/, neck16/, neck17/, neck18/, each with filtered_feature_bc_matrix/
(matrix.mtx, features.tsv, barcodes.tsv; .gz also accepted).
"""
import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread

HERE = Path(__file__).resolve().parent
ACTIVITIES = "results/signatures/per_leaf/Assignment_Solution/Activities/Assignment_Solution_Activities.txt"
SAMPLES = [("neck15", "Normal (N)"), ("neck16", "Metastasis (M)"),
           ("neck17", "Margins (R)"), ("neck18", "Primary (P)")]
SIGNATURE = "SBS32"


def _open(path):
    """path or path.gz"""
    if path.exists():
        return path
    gz = Path(str(path) + ".gz")
    if gz.exists():
        return gz
    raise FileNotFoundError(path)


def load_sample(tenx, folder, sample_name):
    """As load_sample in bef_dgea.ipynb: cells named <folder>_<barcode>."""
    d = Path(tenx) / folder / "filtered_feature_bc_matrix"
    matrix = mmread(str(_open(d / "matrix.mtx"))).T.tocsr()
    genes = pd.read_csv(_open(d / "features.tsv"), header=None, sep="\t")[1].values
    barcodes = pd.read_csv(_open(d / "barcodes.tsv"), header=None, sep="\t")[0].values
    barcodes = [f"{folder}_{b}" for b in barcodes]
    adata = sc.AnnData(X=matrix, var=pd.DataFrame(index=genes), obs=pd.DataFrame(index=barcodes))
    adata.obs["sample"] = sample_name
    dup = adata.var.index.duplicated()
    if dup.sum() > 0:
        print(f"{folder}: removing {dup.sum()} duplicated genes.")
        adata = adata[:, ~dup].copy()
    return adata


def level_dirs(base):
    return sorted(p for p in base.glob("p*") if (p / ACTIVITIES).exists())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", default="neck")
    ap.add_argument("--tenx", required=True, help="folder containing neck15..neck18")
    ap.add_argument("--out", default=str(HERE / "levels"))
    ap.add_argument("--cells", type=int, help="test: random subsample of N tree cells (+ the 'other' cells)")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()

    base = Path(a.out) / a.name
    out = base / "dgea"
    out.mkdir(parents=True, exist_ok=True)
    levels = level_dirs(base)
    if not levels:
        raise SystemExit(f"no finished level (with {ACTIVITIES}) under {base}")
    print("levels:", " ".join(p.name for p in levels))

    acts = {p.name: pd.read_csv(p / ACTIVITIES, sep="\t", index_col=0)[SIGNATURE] for p in levels}

    # expression (once)
    adata = sc.concat([load_sample(a.tenx, f, s) for f, s in SAMPLES],
                      join="outer", label="sample", keys=["N", "M", "R", "P"])
    print(f"10x: {adata.n_obs} cells x {adata.n_vars} genes")
    tree_cells = set().union(*(s.index for s in acts.values()))
    missing = tree_cells - set(adata.obs_names)
    if missing:
        raise SystemExit(f"{len(missing)} tree cells not in the 10x data, e.g. {sorted(missing)[:3]}")

    if a.cells:
        rng = np.random.default_rng(a.seed)
        in_tree = np.array([c in tree_cells for c in adata.obs_names])
        keep = rng.choice(np.flatnonzero(in_tree), size=min(a.cells, in_tree.sum()), replace=False)
        adata = adata[np.sort(np.concatenate([keep, np.flatnonzero(~in_tree)]))].copy()
        print(f"subsample: {adata.n_obs} cells ({len(keep)} in the tree)")

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.filter_genes(adata, min_cells=int(0.1 * adata.n_obs))
    print(f"after filtering: {adata.n_obs} cells x {adata.n_vars} genes")

    # groups per level: SBS32 = 0 -> group1 (L32), > 0 -> group2 (H32), not in tree -> other
    sizes = []
    for lv, s in acts.items():
        g = pd.Series("other", index=adata.obs_names, name="dge_group")
        s = s.reindex(adata.obs_names).dropna()
        g[s.index[s == 0]] = "group1"
        g[s.index[s > 0]] = "group2"
        (out / lv).mkdir(exist_ok=True)
        g.rename_axis("cell").to_csv(out / lv / "groups.csv")
        n = g.value_counts()
        sizes.append({"level": lv, "group1_L32": n.get("group1", 0), "group2_H32": n.get("group2", 0),
                      "other": n.get("other", 0)})
        print(f"{lv}: L32 {sizes[-1]['group1_L32']}  H32 {sizes[-1]['group2_H32']}  other {sizes[-1]['other']}")
    pd.DataFrame(sizes).to_csv(out / "group_sizes.csv", index=False)

    adata.obs["dge_group"] = "other"   # placeholder; dgea_mast.R sets it per level from groups.csv
    tmp = out / "expression.h5ad.tmp"
    adata.write_h5ad(tmp)
    os.replace(tmp, out / "expression.h5ad")
    print("wrote", out / "expression.h5ad")


if __name__ == "__main__":
    main()
