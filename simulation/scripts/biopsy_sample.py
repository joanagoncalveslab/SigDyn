"""
biopsy_sample.py
Mimic a biopsy of a large CellCoal tumour: take ONE clade (a spatially coherent
lineage), sample --n-cells cells from it, and keep only mutations detectable in the
sample (carried by >= --min-freq of the sampled cells, i.e. scRNA variant calling
cannot see mutations present in only a handful of cells).

Mutations above the clade's MRCA become shared by every sampled cell (the biopsy's
"trunk"), which is what makes a biopsy look highly clonal.

Clade choice: among clades with --clade-min..--clade-max cells, the one whose
sample is closest to --target-sites AND --target-shared (mutations carried by every
sampled cell), each as relative error (ties -> larger clade).

Outputs (same format as run_cellcoal_binary.py)
  <outdir>/genotypes.csv   sampled cells x kept sites (0/1)
  <outdir>/true_tree.nwk   true tree pruned to the sampled cells
  optional --clade-outdir: the WHOLE chosen clade (all its cells, every site present
  in it, no frequency filter) + its tree, e.g. to assign signatures before sampling
"""
import argparse
import os

import dendropy
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--tumour", required=True, help="dir with genotypes.csv + true_tree.nwk of the full tumour")
ap.add_argument("--n-cells", type=int, default=688)
ap.add_argument("--clade-min", type=int, default=700)
ap.add_argument("--clade-max", type=int, default=2000)
ap.add_argument("--min-freq", type=float, default=0.01)
ap.add_argument("--target-sites", type=int, default=460)
ap.add_argument("--target-shared", type=int, default=101)
ap.add_argument("--seed", type=int, default=1)
ap.add_argument("--outdir", required=True)
ap.add_argument("--clade-outdir", default=None, help="also write the whole clade (unfiltered) here")
a = ap.parse_args()

geno = pd.read_csv(f"{a.tumour}/genotypes.csv", index_col=0)
tree = dendropy.Tree.get(path=f"{a.tumour}/true_tree.nwk", schema="newick",
                         preserve_underscores=True)


def take(nd, rng):
    leaves = [l.taxon.label for l in nd.leaf_nodes()]
    samp = sorted(rng.choice(leaves, a.n_cells, replace=False))
    s = geno.loc[samp]
    s = s.loc[:, s.mean() >= max(a.min_freq, 1e-12)]
    return samp, s


cands = []
for nd in tree.preorder_node_iter():
    n = len(nd.leaf_nodes())
    if a.clade_min <= n <= a.clade_max and n >= a.n_cells:
        samp, s = take(nd, np.random.default_rng(a.seed))
        err = (abs(s.shape[1] - a.target_sites) / a.target_sites
               + abs(int((s.mean() == 1).sum()) - a.target_shared) / a.target_shared)
        cands.append((err, -n, nd, samp, s))
if not cands:
    raise SystemExit(f"no clade with {a.clade_min}-{a.clade_max} cells; change the range")
for c in sorted(cands, key=lambda x: (x[0], x[1])):
    s, n = c[4], -c[1]
    b = s.sum(1)
    print(f"  clade {n:5d} cells: sites {s.shape[1]:4d}  density {s.values.mean():.2f}  "
          f"shared-by-all {(s.mean() == 1).sum():3d}  per-cell min/median/max "
          f"{b.min()}/{int(b.median())}/{b.max()}")
_, negn, node, samp, sub = sorted(cands, key=lambda x: (x[0], x[1]))[0]

pruned = dendropy.Tree.get(path=f"{a.tumour}/true_tree.nwk", schema="newick",
                           preserve_underscores=True)
pruned.retain_taxa_with_labels(samp)
pruned.is_rooted = True

os.makedirs(a.outdir, exist_ok=True)
sub.index.name = "cell"
sub.to_csv(f"{a.outdir}/genotypes.csv")
with open(f"{a.outdir}/true_tree.nwk", "w") as fh:
    fh.write(pruned.as_string(schema="newick").strip() + "\n")
b = sub.sum(1)
print(f"BIOPSY: clade of {-negn} cells -> {len(samp)} sampled, {sub.shape[1]} sites, "
      f"density {sub.values.mean():.2f}, shared-by-all {(sub.mean() == 1).sum()}, "
      f"per-cell min/median/max {b.min()}/{int(b.median())}/{b.max()}")

if a.clade_outdir:
    clade_cells = sorted(l.taxon.label for l in node.leaf_nodes())
    cg = geno.loc[clade_cells]
    cg = cg.loc[:, cg.sum() > 0]
    ct = dendropy.Tree.get(path=f"{a.tumour}/true_tree.nwk", schema="newick",
                           preserve_underscores=True)
    ct.retain_taxa_with_labels(clade_cells)
    ct.is_rooted = True
    os.makedirs(a.clade_outdir, exist_ok=True)
    cg.index.name = "cell"
    cg.to_csv(f"{a.clade_outdir}/genotypes.csv")
    with open(f"{a.clade_outdir}/true_tree.nwk", "w") as fh:
        fh.write(ct.as_string(schema="newick").strip() + "\n")
    print(f"CLADE: {len(clade_cells)} cells, {cg.shape[1]} sites (unfiltered) -> {a.clade_outdir}")
