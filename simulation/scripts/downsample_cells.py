"""
downsample_cells.py
Per-cell downsampling: every cell keeps at most --cap of its mutations (chosen at
random, the rest set to R). Cells already at or below the cap are untouched, so
for caps below the least-mutated cell every cell has exactly --cap mutations.
Sites left with no carrier are dropped (PhylinSic tables never contain them).

  python downsample_cells.py in.txt out.txt --cap 100 --seed 1      (--cap 0 = keep all)
"""
import argparse
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("inp")
ap.add_argument("out")
ap.add_argument("--cap", type=int, required=True)
ap.add_argument("--seed", type=int, default=1)
a = ap.parse_args()

df = pd.read_csv(a.inp, sep="\t", index_col=0)
rng = np.random.default_rng(a.seed)
if a.cap > 0:
    for cell in df.columns:
        mut = np.flatnonzero(df[cell].isin(["H", "A"]).to_numpy())
        if len(mut) > a.cap:
            drop = rng.choice(mut, size=len(mut) - a.cap, replace=False)
            df.iloc[drop, df.columns.get_loc(cell)] = "R"
df = df.loc[df.isin(["H", "A"]).any(axis=1)]
df.to_csv(a.out, sep="\t")
b = df.isin(["H", "A"]).sum()
print(f"cap={a.cap or 'all'}: {len(df)} sites; per-cell min/median/max {b.min()}/{int(b.median())}/{b.max()}")
