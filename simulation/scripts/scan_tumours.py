"""
scan_tumours.py
Search CellCoal tumours (seeds x tumour sizes x n-sites) for a biopsy clade close to
688 cells whose sample looks like the real dataset: ~460 sites and a ~101-mutation trunk after the
rare-site filter, with enough mutations left below the trunk for clones.

Every tumour uses run_cellcoal_binary.py (growth 1e-4, Ne 10000); seed / cells /
sites / mutation rate vary. Expected mutations scale with n-sites x mut-rate. CellCoal
aborts when a genealogy needs more mutations than sites (infinite-sites check), so to
avoid failures use more sites with a proportionally lower rate, e.g.
--n-sites 11600 --mut-rate 1.75e-6 gives the same expected burden as 5800 / 3.5e-6.
Sampling of cells is identical to biopsy_sample.py (same rng, same seed), so the
numbers here are what biopsy_sample.py will produce for that clade.

Per candidate clade (after dropping sites seen in < --min-freq of sampled cells):
  sites, shared-by-all (trunk), per-cell min/median/max,
  the two largest sub-clades below the biopsy root and how many mutations sit
  on/below each (= what a clone starting there could carry).

Output: <outdir>/scan_results.csv (best first). Tumour folders are deleted except the
--keep-top best ones (CellCoal raw output is always removed to save disk).

Example:
  python scripts/scan_tumours.py --seeds 1-20 --tumour-cells 1500 2000 3000 \
      --n-sites 11600 14000 --mut-rate 1.75e-6 --outdir output/scan
"""
import argparse
import os
import shutil
import subprocess
import sys

import dendropy
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RUN_CC = os.path.join(HERE, "run_cellcoal_binary.py")
CC_BIN = os.environ.get("CELLCOAL", os.path.join(HERE, "..", "bin", "cellcoal-1.1.0"))

ap = argparse.ArgumentParser()
ap.add_argument("--seeds", default="1-20", help="e.g. 1-20 or 1,5,9")
ap.add_argument("--tumour-cells", type=int, nargs="+", default=[3000])
ap.add_argument("--n-sites", type=int, nargs="+", default=[11600])
ap.add_argument("--mut-rate", default="1.75e-6", help="CellCoal somatic mutation rate")
ap.add_argument("--n-cells", type=int, default=688, help="cells in the biopsy")
ap.add_argument("--clade-min", type=int, default=688)
ap.add_argument("--clade-max", type=int, default=800)
ap.add_argument("--min-freq", type=float, default=0.20)
ap.add_argument("--target-sites", type=int, default=460)
ap.add_argument("--target-shared", type=int, default=101)
ap.add_argument("--sample-seed", type=int, default=1, help="same as biopsy_sample.py --seed")
ap.add_argument("--keep-top", type=int, default=3)
ap.add_argument("--outdir", required=True)
a = ap.parse_args()

if "-" in a.seeds:
    lo, hi = map(int, a.seeds.split("-")); seeds = list(range(lo, hi + 1))
else:
    seeds = [int(s) for s in a.seeds.split(",")]
os.makedirs(a.outdir, exist_ok=True)


def branch_counts(sub, node, sampled):
    """mutations on the branch above `node` and below it (within its sampled cells)."""
    L = frozenset(l.taxon.label for l in node.leaf_nodes()) & sampled
    carriers = {c: frozenset(sub.index[sub[c] == 1]) for c in sub.columns}
    on = sum(1 for s in carriers.values() if s == L)
    below = sum(1 for s in carriers.values() if s < L)
    return len(L), on, below


rows = []
for cells in a.tumour_cells:
    for nsites in a.n_sites:
        for seed in seeds:
            tag = f"tumour_s{seed}_c{cells}_l{nsites}_u{a.mut_rate}"
            d = os.path.join(a.outdir, tag)
            if not os.path.exists(os.path.join(d, "genotypes.csv")):
                # run inside outdir with a short relative path: CellCoal truncates paths > 120 chars
                subprocess.run([sys.executable, RUN_CC, "--cellcoal", CC_BIN,
                                "--n-cells", str(cells), "--n-sites", str(nsites),
                                "--mut-rate", a.mut_rate, "--growth", "1e-4", "--Ne", "10000",
                                "--seed", str(seed), "--outdir", tag],
                               capture_output=True, text=True, cwd=a.outdir)
                shutil.rmtree(os.path.join(d, "_cellcoal_raw"), ignore_errors=True)
            if not os.path.exists(os.path.join(d, "genotypes.csv")):
                print(f"{tag}: CellCoal produced no genotypes (usually the infinite-sites check: "
                      "use more --n-sites with a lower --mut-rate) - skipped", flush=True)
                rows.append({"tumour": tag, "seed": seed, "tumour_cells": cells, "n_sites": nsites,
                             "mut_rate": a.mut_rate, "status": "cellcoal_failed"})
                continue

            geno = pd.read_csv(os.path.join(d, "genotypes.csv"), index_col=0)
            tree = dendropy.Tree.get(path=os.path.join(d, "true_tree.nwk"), schema="newick",
                                     preserve_underscores=True)
            n_cand = 0
            for nd in tree.preorder_node_iter():
                n = len(nd.leaf_nodes())
                if not (a.clade_min <= n <= a.clade_max and n >= a.n_cells):
                    continue
                n_cand += 1
                # identical to biopsy_sample.py: take()
                leaves = [l.taxon.label for l in nd.leaf_nodes()]
                samp = sorted(np.random.default_rng(a.sample_seed).choice(leaves, a.n_cells, replace=False))
                s = geno.loc[samp]
                s = s.loc[:, s.mean() >= a.min_freq]
                b = s.sum(1)
                shared = int((s.mean() == 1).sum())
                sampled = frozenset(samp)
                kids = sorted((branch_counts(s, ch, sampled) for ch in nd.child_nodes()),
                              key=lambda x: -x[0])[:2]
                err = (abs(s.shape[1] - a.target_sites) / a.target_sites
                       + abs(shared - a.target_shared) / a.target_shared)
                r = {"tumour": tag, "seed": seed, "tumour_cells": cells, "n_sites": nsites,
                     "mut_rate": a.mut_rate, "status": "ok", "clade_size": n, "sites": s.shape[1], "shared_by_all": shared,
                     "cell_min": int(b.min()), "cell_median": int(b.median()), "cell_max": int(b.max()),
                     "score": round(err, 3)}
                for i, (sz, on, below) in enumerate(kids, 1):
                    r[f"sub{i}_cells"], r[f"sub{i}_muts_on_branch"], r[f"sub{i}_muts_below"] = sz, on, below
                rows.append(r)
                print(f"{tag}: clade {n} -> sites {s.shape[1]}, shared {shared}, per-cell "
                      f"{b.min()}/{int(b.median())}/{b.max()}, sub-clades "
                      + ", ".join(f"{sz} cells ({on}+{below} muts)" for sz, on, below in kids), flush=True)
            if n_cand == 0:
                print(f"{tag}: no clade with {a.clade_min}-{a.clade_max} cells", flush=True)
                rows.append({"tumour": tag, "seed": seed, "tumour_cells": cells, "n_sites": nsites,
                             "mut_rate": a.mut_rate, "status": "no_clade_in_range"})

res = pd.DataFrame(rows)
ok = res[res.status == "ok"].sort_values("score") if "score" in res else res.iloc[0:0]
res = pd.concat([ok, res[res.status != "ok"]])
res.to_csv(os.path.join(a.outdir, "scan_results.csv"), index=False)

keep = set(ok.tumour.drop_duplicates().head(a.keep_top))
for tag in set(res.tumour) - keep:
    shutil.rmtree(os.path.join(a.outdir, tag), ignore_errors=True)

pd.set_option("display.width", 250, "display.max_columns", 30)
print(f"\nbest candidates (kept tumours: {sorted(keep)}):")
print(ok.head(15).to_string(index=False))
print(f"\nwrote {os.path.join(a.outdir, 'scan_results.csv')}")
