"""
assign_signatures_onsets.py
Same additive model as assign_signatures.py (in this folder; imported,
not modified), but YOU choose where each clone starts instead of the script taking
the largest clades. Cells outside every clone stay SBS1-only.

  * SBS1 (base) active everywhere at --base-weight
  * clone i starts at the node with exactly --onset-sizes[i] cells and adds
    --clone-signatures[i] at --clone-rates[i] on top of SBS1
  * a mutation's region = the clone whose subtree contains the MRCA of its carriers
    (else base); signature drawn from that region's mixture, then an SBS96 channel

Draw order and RNG use are identical to assign_signatures.py, so with the same clone nodes
the output is identical (checked against assign_signatures.py for the 434/254 clones).

Pick node sizes from inspect_biopsy_tree.py.

Labelling a LARGER set of cells than the biopsy (e.g. the whole clade, before
sampling 688 of its cells): pass that set's genotypes/tree and
--onset-ref-tree <biopsy>/true_tree.nwk. --onset-sizes then refer to nodes of the
biopsy tree (the map printed by inspect_biopsy_tree.py); each clone starts at the
node of the larger tree that is the MRCA of that biopsy node's cells.

  python scripts/assign_signatures_onsets.py --genotypes .../genotypes.csv \
      --tree .../true_tree.nwk --clone-signatures SBS5 SBS26 --clone-rates 1 1 \
      --onset-sizes 434 178 --signatures-file COSMIC_v3.4_SBS_GRCh38.txt \
      --out .../mutation_types.csv
"""
import argparse
import importlib.util
import os

import dendropy
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("assign_signatures",
                                              os.path.join(HERE, "assign_signatures.py"))
yours = importlib.util.module_from_spec(spec)
spec.loader.exec_module(yours)
CHANNELS, _norm, load_signatures = yours.CHANNELS, yours._norm, yours.load_signatures

ap = argparse.ArgumentParser()
ap.add_argument("--genotypes", required=True)
ap.add_argument("--tree", required=True)
ap.add_argument("--base-signature", default="SBS1")
ap.add_argument("--base-weight", type=float, default=1.0)
ap.add_argument("--clone-signatures", nargs="+", required=True)
ap.add_argument("--clone-rates", nargs="+", type=float, required=True)
ap.add_argument("--onset-sizes", nargs="+", type=int, required=True,
                help="cell count of the node where each clone starts (must be unique in the tree)")
ap.add_argument("--onset-ref-tree", default=None,
                help="biopsy tree: --onset-sizes are looked up here and mapped to --tree via MRCA")
ap.add_argument("--signatures-file", default=None)
ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--out", required=True)
args = ap.parse_args()
if not (len(args.clone_signatures) == len(args.clone_rates) == len(args.onset_sizes)):
    raise SystemExit("ERROR: --clone-signatures, --clone-rates and --onset-sizes must have the same length")

rng = np.random.default_rng(args.seed)
geno = pd.read_csv(args.genotypes, index_col=0)
tns = dendropy.TaxonNamespace()
tree = dendropy.Tree.get(path=args.tree, schema="newick", preserve_underscores=True, taxon_namespace=tns)
tree.is_rooted = True
taxon = {t.label: t for t in tns}

names = [args.base_signature] + list(args.clone_signatures)
sigs = load_signatures(args.signatures_file, names)

# --- the only difference from assign_signatures.py: onset nodes chosen by size ---
ref = (dendropy.Tree.get(path=args.onset_ref_tree, schema="newick", preserve_underscores=True)
       if args.onset_ref_tree else tree)
onset_nodes = []
for size in args.onset_sizes:
    hits = [nd for nd in ref.preorder_node_iter() if len(nd.leaf_nodes()) == size]
    if len(hits) != 1:
        raise SystemExit(f"ERROR: {len(hits)} nodes have exactly {size} cells; need exactly 1 "
                         "(see inspect_biopsy_tree.py)")
    if ref is tree:
        onset_nodes.append(hits[0])
    else:   # biopsy node -> tumour node = MRCA of its cells in the tumour tree
        cells = [l.taxon.label for l in hits[0].leaf_nodes()]
        onset_nodes.append(tree.mrca(taxa=[taxon[c] for c in cells]))
leafsets = [set(nd.leaf_nodes()) for nd in onset_nodes]
for i in range(len(leafsets)):
    for j in range(i + 1, len(leafsets)):
        if leafsets[i] & leafsets[j]:
            raise SystemExit(f"ERROR: clones starting at {args.onset_sizes[i]} and "
                             f"{args.onset_sizes[j]} cells overlap; choose disjoint nodes")

# --- from here on identical to assign_signatures.py ---
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
os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
out.to_csv(args.out, index=False)

in_clone = set()
biopsy = {l.taxon.label for l in ref.leaf_nodes()} if ref is not tree else None
print("clones (chosen start nodes, disjoint):")
for label, sig, rate, onode in clones:
    cells = {l.taxon.label for l in onode.leaf_nodes()}
    in_clone |= cells
    extra = f" ({len(cells & biopsy)} of them in the biopsy)" if biopsy else ""
    print(f"  {label}: {len(cells)} cells{extra}, mix = {args.base_signature}:{args.base_weight} + {sig}:{rate}")
print(f"  SBS1-only (outside every clone): {len(geno) - len(in_clone)} cells"
      + (f" ({len(biopsy - in_clone)} of them in the biopsy)" if biopsy else ""))
print("\nmutations by clone x generating signature (all labelled mutations, before biopsy/filter):")
print(pd.crosstab(out["clone"], out["signature"]).to_string())
print(f"\nwrote {args.out}  ({len(out)} mutations)")
