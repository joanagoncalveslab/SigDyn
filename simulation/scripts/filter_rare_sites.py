"""
filter_rare_sites.py
Drop mutations seen in fewer than --min-freq of the biopsy's cells (scRNA variant
calling can't detect them). Applied AFTER signatures are assigned, so the labels
and clone placement come from the full biopsy.

Writes filtered genotypes.csv + mutation_types.csv (+ copies true_tree.nwk) to --outdir
and prints per-cell burden and clone-signature shares.

  python scripts/filter_rare_sites.py --genotypes biopsy/genotypes.csv \
      --mutation-types biopsy/mutation_types.csv --tree biopsy/true_tree.nwk \
      --min-freq 0.20 --outdir cellcoal_out
"""
import argparse
import os
import shutil

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--genotypes", required=True)
ap.add_argument("--mutation-types", required=True)
ap.add_argument("--tree", required=True)
ap.add_argument("--min-freq", type=float, default=0.20)
ap.add_argument("--outdir", required=True)
a = ap.parse_args()

g = pd.read_csv(a.genotypes, index_col=0)
m = pd.read_csv(a.mutation_types)
keep = g.columns[g.mean() >= a.min_freq]
g = g[keep]
m = m[m.site.astype(str).isin(set(keep.astype(str)))]

os.makedirs(a.outdir, exist_ok=True)
g.to_csv(f"{a.outdir}/genotypes.csv")
m.to_csv(f"{a.outdir}/mutation_types.csv", index=False)
shutil.copy(a.tree, f"{a.outdir}/true_tree.nwk")

b = g.sum(1)
print(f"kept {g.shape[1]} sites (>= {a.min_freq:.0%} of {len(g)} cells); shared by all {(g.mean() == 1).sum()}; "
      f"per-cell min/median/max {b.min()}/{int(b.median())}/{b.max()}")
print(pd.crosstab(m["clone"], m["signature"]).to_string())
gs = g.copy(); gs.columns = gs.columns.astype(str)
sig = m.assign(site=m.site.astype(str)).set_index("site").signature
t = pd.DataFrame({s: gs[sig.index[sig == s]].sum(1) for s in sorted(sig.unique())})
f = t.div(t.sum(1), axis=0)
for s in [c for c in t.columns if c != "SBS1"]:
    c = f[s][t[s] > 0]
    print(f"{s}: {len(c)} cells, share in those cells median {c.median():.2f}, "
          f"{(m.signature == s).sum()} distinct mutations")
only = (t.drop(columns=["SBS1"], errors="ignore").sum(1) == 0).sum()
print(f"SBS1-only cells: {only}")
