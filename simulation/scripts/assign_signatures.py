"""
assign_signatures.py  (Step 2, additive model)
Type each mutation in genotypes.csv by the signature that generated it (no
mutations are ever removed), given
per-clone signature ACTIVITY.

  * SBS1 active EVERYWHERE at a constant weight (background clock).
  * Each clone = a disjoint subtree with one onset signature added ON TOP of SBS1
    at a given rate (mixture weight): SBS40c high in one clone, SBS26 very high in
    a parallel clone. Disjoint -> SBS40c and SBS26 never co-occur in a cell.

A mutation's region = the deepest onset subtree containing the MRCA of its
carrying cells (else base). In that region draw the generating signature from the
additive mixture {SBS1: base_weight, onset: rate}, then draw an SBS96 channel.
"rate" = exposure weight (fraction of a clone's mutations from that signature),
not extra mutation count.
"""
from __future__ import annotations
import argparse, hashlib
import numpy as np, pandas as pd, dendropy

SUBS = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
BASES = ["A", "C", "G", "T"]
CHANNELS = [f"{a}[{s}]{b}" for s in SUBS for a in BASES for b in BASES]

def _norm(v):
    v = np.asarray(v, float); return v / v.sum()

def toy_signature(name):
    if name in ("SBS1", "SBS1-like"):
        v = np.full(96, 0.001)
        for i, ch in enumerate(CHANNELS):
            if ch[2:5] == "C>T" and ch.endswith("G"):
                v[i] += 0.25
        return _norm(v)
    seed = int(hashlib.md5(name.encode()).hexdigest(), 16) % (2**32)
    return _norm(np.random.default_rng(seed).dirichlet(np.ones(96) * 0.3))

def load_signatures(path, names):
    if path is None:
        return {n: toy_signature(n) for n in names}
    df = pd.read_csv(path, sep=None, engine="python")
    df = df.set_index(df.columns[0]).loc[CHANNELS]
    return {n: _norm(df[n].to_numpy()) for n in names}

def pick_clone_subtrees(tree, k, min_frac, max_frac, rng):
    """Cut the tree into DISJOINT clades no bigger than max_frac, then take the
    k largest that are >= min_frac. Guarantees parallel (disjoint) clones and
    lets 'earlier' just mean a larger max_frac."""
    N = len(tree.leaf_nodes())
    frontier = list(tree.seed_node.child_nodes())
    changed = True
    while changed:
        changed = False
        nxt = []
        for nd in frontier:
            if len(nd.leaf_nodes()) / N > max_frac and nd.child_nodes():
                nxt.extend(nd.child_nodes()); changed = True
            else:
                nxt.append(nd)
        frontier = nxt
    cand = sorted([nd for nd in frontier if len(nd.leaf_nodes()) / N >= min_frac],
                  key=lambda nd: -len(nd.leaf_nodes()))
    return cand[:k]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genotypes", default="cellcoal_out/genotypes.csv")
    ap.add_argument("--tree", default="cellcoal_out/true_tree.nwk")
    ap.add_argument("--base-signature", default="SBS1")
    ap.add_argument("--base-weight", type=float, default=1.0)
    ap.add_argument("--clone-signatures", nargs="*", default=["SBS26", "SBS5"])
    ap.add_argument("--clone-rates", nargs="*", type=float, default=[8.0, 20.0])
    ap.add_argument("--min-frac", type=float, default=0.2)
    ap.add_argument("--max-frac", type=float, default=0.6)
    ap.add_argument("--signatures-file", default=None)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="cellcoal_out/mutation_types.csv")
    ap.add_argument("--root-candidates-out", default="cellcoal_out/root_candidates.csv",
                    help="list of naturally SBS1-only cells (ancestral region), by burden")
    args = ap.parse_args()

    if not args.clone_signatures:
        raise SystemExit("ERROR: no --clone-signatures given, so every "
            "mutation is typed as the base signature. Pass names, e.g. "
            "--clone-signatures SBS40c SBS26 --clone-rates 8 20")
    if len(args.clone_rates) != len(args.clone_signatures):
        raise SystemExit("ERROR: --clone-rates must match --clone-signatures in count")
    rng = np.random.default_rng(args.seed)
    geno = pd.read_csv(args.genotypes, index_col=0)
    tns = dendropy.TaxonNamespace()
    tree = dendropy.Tree.get(path=args.tree, schema="newick",
                             preserve_underscores=True, taxon_namespace=tns)
    tree.is_rooted = True
    taxon = {t.label: t for t in tns}

    names = [args.base_signature] + list(args.clone_signatures)
    sigs = load_signatures(args.signatures_file, names)

    onset_nodes = pick_clone_subtrees(tree, len(args.clone_signatures),
                                      args.min_frac, args.max_frac, rng)
    if len(onset_nodes) < len(args.clone_signatures):
        print(f"WARNING: only placed {len(onset_nodes)} of "
              f"{len(args.clone_signatures)} clones — widen --max-frac or "
              f"lower --min-frac (tree split too uneven).")
    node_region = {nd: "root" for nd in tree.preorder_node_iter()}
    region_mix = {"root": {args.base_signature: args.base_weight}}
    clones = []
    for (sig, rate), onode in zip(zip(args.clone_signatures, args.clone_rates), onset_nodes):
        label = f"clone_{sig}"
        clones.append((label, sig, rate, onode))
        region_mix[label] = {args.base_signature: args.base_weight, sig: rate}
        for nd in onode.preorder_iter():
            node_region[nd] = label

    def draw(region):
        mix = region_mix[region]; comps = list(mix)
        p = _norm([mix[c] for c in comps])
        comp = comps[rng.choice(len(comps), p=p)]
        return comp, CHANNELS[rng.choice(96, p=sigs[comp])]

    records = []
    for site in geno.columns:
        carriers = geno.index[geno[site] == 1].tolist()
        if not carriers:
            continue
        mrca = (tree.find_node_with_taxon_label(carriers[0]) if len(carriers) == 1
                else tree.mrca(taxa=[taxon[c] for c in carriers]))
        region = node_region.get(mrca, "root")
        comp, ch = draw(region)
        records.append({"site": site, "n_cells": len(carriers), "clone": region,
                        "signature": comp, "sbs96": ch, "sub_type": ch[2:5],
                        "context": f"{ch[0]}{ch[2]}{ch[-1]}"})

    out = pd.DataFrame(records)
    out.to_csv(args.out, index=False)

    # --- suggest a root: least-mutated NATURALLY SBS1-only cell (no edits) ---
    onset_sites = [str(s) for s in out.loc[out["signature"].isin(args.clone_signatures), "site"]
                   if str(s) in geno.columns]
    burden = geno.sum(axis=1)
    carries_onset = geno[onset_sites].sum(axis=1) > 0 if onset_sites else pd.Series(False, index=geno.index)
    sbs1_only = geno.index[~carries_onset].tolist()

    print("\nclones (parallel, disjoint subtrees):")
    for label, sig, rate, onode in clones:
        print(f"  {label}: {len(onode.leaf_nodes())} cells, mix = "
              f"{args.base_signature}:{args.base_weight} + {sig}:{rate}")
    print("\nmutations by clone x generating signature:")
    print(pd.crosstab(out["clone"], out["signature"]).to_string())

    if not sbs1_only:
        print("\nWARNING: no naturally SBS1-only cells — the clones cover the whole "
              "tree. Lower --max-frac so an ancestral SBS1-only region remains.")
    else:
        cand = burden[sbs1_only].sort_values()
        pd.DataFrame({"cell": cand.index, "burden": cand.values}).to_csv(
            args.root_candidates_out, index=False)
        gmin = burden.idxmin()
        print(f"\n{len(sbs1_only)} cells are naturally SBS1-only (ancestral region, NO edits).")
        print(f"suggested root = least-mutated SBS1-only cell: {cand.index[0]} "
              f"({int(cand.iloc[0])} muts)  ->  {args.root_candidates_out}")
        print(f"global least-mutated cell {gmin} ({int(burden.min())} muts) "
              f"is SBS1-only: {gmin in set(sbs1_only)}")

    onset_sigs = list(args.clone_signatures)
    cell_has = {s: set() for s in onset_sigs}
    for _, rr in out.iterrows():
        if rr["signature"] in onset_sigs:
            cell_has[rr["signature"]].update(geno.index[geno[rr["site"]] == 1])
    overlap = set.intersection(*cell_has.values()) if len(onset_sigs) > 1 else set()
    print(f"\ncells carrying BOTH onset signatures: {len(overlap)} (should be 0)")
    print(f"wrote {args.out}  ({len(out)} mutations)")

if __name__ == "__main__":
    main()
