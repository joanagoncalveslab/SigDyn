"""
evaluate_leaves.py
Leaf-level recovery of the simulated signatures for every run in runs/.

Truth per leaf = counts of the mutations that SURVIVED masking in that cell,
split by the signature that generated them (sim/site_signature_map.csv).

Per run it reports
  * which signatures SigDyn called active (Custom_COSMIC) and the global fit
  * leaf exposures from exposures_active_signatures.txt (SigDyn's final output)
      - cosine(truth, estimate) per cell, fraction correlation / MAE per signature
      - clone detection: a cell "has" SBS5/SBS26 if its estimated fraction >= 10%
  * the unrestricted per-leaf COSMIC fit: share of exposure on signatures that
    were never simulated (spurious mass)

Writes results/leaf_metrics.csv (one row per run) and results/leaf_per_cell.csv.

  python evaluate_leaves.py                                   # masking runs (sim/, runs/p*_r*)
  python evaluate_leaves.py --sim from_scratch_460/sim \
         --runs "from_scratch_460/run" --out from_scratch_460/results
"""
import argparse
import glob
import os
import re

import numpy as np
import pandas as pd

M = os.environ.get("SIM_OUT", os.getcwd())   # pipeline output folder (set by main.sh)
ap = argparse.ArgumentParser()
ap.add_argument("--sim", default="sim", help="simulation dir, relative to the output folder")
ap.add_argument("--runs", default="runs/p*_r*", help="glob of SigDyn run dirs, relative to the output folder")
ap.add_argument("--out", default="results", help="output dir, relative to the output folder")
ap.add_argument("--merge", nargs="*", default=[], metavar="EST:SIM",
                help="count an estimated signature as a simulated one, e.g. SBS12:SBS26")
args = ap.parse_args()
SIM, OUT = os.path.join(M, args.sim), os.path.join(M, args.out)   # absolute paths also work
SIGS = ["SBS1", "SBS5", "SBS26"]
CLONE_SIGS = ["SBS5", "SBS26"]
CALL_FRAC = 0.10

sitemap = pd.read_csv(f"{SIM}/site_signature_map.csv").set_index("site_id")
geno_full = pd.read_csv(f"{SIM}/cellcoal_out/genotypes.csv", index_col=0)
geno_full.columns = geno_full.columns.astype(str)
mt = pd.read_csv(f"{SIM}/cellcoal_out/mutation_types.csv")
# true clone of each cell = the clone signature it carries (none -> SBS1-only)
cell_clone = pd.Series("SBS1-only", index=geno_full.index)
for s in CLONE_SIGS:
    sites = mt.loc[mt.signature == s, "site"].astype(str)
    cell_clone[geno_full[sites].sum(axis=1) > 0] = f"clone_{s}"


def cos(a, b):
    return (a * b).sum(1) / np.sqrt((a ** 2).sum(1) * (b ** 2).sum(1))


rows, percell = [], []
for d in sorted(glob.glob(f"{M}/{args.runs}")):
    exp_f = f"{d}/results/signatures/exposures_active_signatures.txt"
    if not os.path.exists(exp_f):
        continue
    hit = re.findall(r"p(\d+)_r(\d+)$", d)
    pct, rep = map(int, hit[0]) if hit else (0, 1)
    cap = re.findall(r"cap(\d+|all)_r(\d+)$", d)      # per-cell downsampling runs
    if cap:
        pct, rep = (0 if cap[0][0] == "all" else int(cap[0][0])), int(cap[0][1])

    G = pd.read_csv(f"{d}/results/genotypes/genotypes_no_PON.txt", sep="\t", index_col=0)
    sig = sitemap.signature.reindex(G.index)
    alt = G == "H"
    truth = pd.DataFrame({s: alt[sig == s].sum() for s in SIGS})

    E = pd.read_csv(exp_f, sep="\t", index_col=0).T
    for pair in args.merge:                      # e.g. SBS12:SBS26 -> SBS12 counted as SBS26
        src, dst = pair.split(":")
        if src in E:
            E[dst] = E.get(dst, 0) + E.pop(src)
    leaves = [c for c in E.index if c in truth.index]
    t = truth.loc[leaves]
    # keep EVERY estimated signature: non-simulated ones count as error (truth = 0)
    cols = SIGS + [c for c in E.columns if c not in SIGS]
    est_all = E.loc[leaves].reindex(columns=cols).fillna(0.0)
    t_all = t.reindex(columns=cols).fillna(0.0)
    est = est_all[SIGS]
    active = list(pd.read_csv(f"{d}/results/signatures/Custom_COSMIC_SBS_v3.4.txt",
                              sep="\t", index_col=0).columns)

    tf = t.div(t.sum(1), axis=0)
    ef = est_all.div(est_all.sum(1).replace(0, np.nan), axis=0).fillna(0)   # share of ALL exposure
    c = cos(t_all, est_all)

    # unrestricted per-leaf COSMIC fit (first SigDyn step)
    pl = pd.read_csv(f"{d}/results/signatures/per_leaf/Assignment_Solution/Activities/"
                     "Assignment_Solution_Activities.txt", sep="\t", index_col=0)
    pl = pl.loc[pl.index.intersection(leaves)]
    spurious = pl.drop(columns=[s for s in SIGS if s in pl], errors="ignore").sum(1) / pl.sum(1)

    glob_fit = pd.read_csv(f"{d}/results/signatures/all_leaves/Assignment_Solution/Activities/"
                           "Assignment_Solution_Activities.txt", sep="\t", index_col=0).iloc[0]
    glob_fit = glob_fit[glob_fit > 0]

    r = {"mask_pct": pct, "cap_per_cell": pct if cap else None, "rep": rep, "n_sites": len(G),
         "median_muts_per_leaf": int(t.sum(1).median()),
         "active_signatures": ";".join(active),
         "global_fit": ";".join(f"{k}:{int(v)}" for k, v in glob_fit.items()),
         "merged": ";".join(args.merge),
         "cosine_median": c.median(), "cosine_p10": c.quantile(0.1),
         "non_simulated_share_mean": ef[[x for x in cols if x not in SIGS]].sum(1).mean(),
         "spurious_frac_perleaf_cosmic_median": spurious.median()}
    for s in SIGS:
        r[f"{s}_frac_r"] = np.corrcoef(tf[s], ef[s])[0, 1] if ef[s].std() > 0 else np.nan
        r[f"{s}_frac_MAE"] = (tf[s] - ef[s]).abs().mean()
    for s in CLONE_SIGS:
        in_clone = cell_clone.loc[leaves] == f"clone_{s}"
        called = ef[s] >= CALL_FRAC
        r[f"{s}_sensitivity"] = (called & in_clone).sum() / in_clone.sum()
        r[f"{s}_specificity"] = (~called & ~in_clone).sum() / (~in_clone).sum()
    # leaves assigned to the right group (SBS1-only / clone_SBS5 / clone_SBS26)
    pred = pd.Series("SBS1-only", index=leaves)
    for s in CLONE_SIGS:
        pred[(ef[s] >= CALL_FRAC) & (ef[s] >= ef[CLONE_SIGS].max(1))] = f"clone_{s}"
    r["clone_assignment_accuracy"] = (pred == cell_clone.loc[leaves]).mean()
    rows.append(r)

    pc = pd.concat([t.add_prefix("true_"), est.add_prefix("est_")], axis=1)
    pc["est_non_simulated"] = est_all[[x for x in cols if x not in SIGS]].sum(1)
    pc["true_clone"] = cell_clone.loc[leaves]; pc["pred_clone"] = pred; pc["cosine"] = c
    pc["mask_pct"], pc["rep"] = pct, rep
    percell.append(pc.rename_axis("cell").reset_index())

os.makedirs(OUT, exist_ok=True)
res = pd.DataFrame(rows).sort_values(["mask_pct", "rep"])
if res["cap_per_cell"].notna().all():   # downsampling: 0 = no cap, show it first
    res = res.drop(columns="mask_pct").sort_values(
        ["cap_per_cell", "rep"], key=lambda c: c.replace(0, 10**6) if c.name == "cap_per_cell" else c,
        ascending=[False, True])
else:
    res = res.drop(columns="cap_per_cell")
res.to_csv(f"{OUT}/leaf_metrics.csv", index=False)
if percell:
    pd.concat(percell).to_csv(f"{OUT}/leaf_per_cell.csv", index=False)
pd.set_option("display.width", 250, "display.max_columns", 40)
print(res.drop(columns=["global_fit"]).round(3).to_string(index=False))
