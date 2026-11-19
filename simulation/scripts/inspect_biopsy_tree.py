"""
inspect_biopsy_tree.py
Print the biopsy's true tree as a node map, to choose where clones start.

For each node (>= --min-cells cells) it shows, after the rare-site filter
(sites carried by >= --min-freq of cells, i.e. what SigDyn will actually see):
  * mutations on the branch above the node
  * mutations further below it, inside the node's subtree
  * their sum = mutations a clone starting at this node would carry
A clone that starts at a node also takes the branch above it. Cells outside every
clone stay SBS1-only, so choosing deeper nodes creates an SBS1-only group but gives
the clones fewer mutations.

Use the printed cell counts with assign_signatures_onsets.py --onset-sizes; counts
marked (NOT UNIQUE) can't be used that way.

  python scripts/inspect_biopsy_tree.py --dir output/<name>/sim/biopsy
"""
import argparse
from collections import Counter

import dendropy
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--dir", required=True, help="folder with genotypes.csv + true_tree.nwk")
ap.add_argument("--min-freq", type=float, default=0.20)
ap.add_argument("--min-cells", type=int, default=30)
ap.add_argument("--depth", type=int, default=4)
a = ap.parse_args()

g = pd.read_csv(f"{a.dir}/genotypes.csv", index_col=0)
N = len(g)
gf = g.loc[:, g.mean() >= a.min_freq]
t = dendropy.Tree.get(path=f"{a.dir}/true_tree.nwk", schema="newick", preserve_underscores=True)

carriers_all = Counter(frozenset(g.index[g[c] == 1]) for c in g.columns)
carriers_f = Counter(frozenset(gf.index[gf[c] == 1]) for c in gf.columns)
size_count = Counter(len(nd.leaf_nodes()) for nd in t.preorder_node_iter())


def counts(L, carriers):
    on = carriers.get(L, 0)
    below = sum(v for s, v in carriers.items() if s < L)
    return on, below


print(f"{N} cells; {g.shape[1]} sites before filter, {gf.shape[1]} after (>= {a.min_freq:.0%} of cells)")
print("node: cells (% of biopsy) | after filter: on branch + below = clone total  [before filter: total]\n")


def show(nd, depth):
    L = frozenset(l.taxon.label for l in nd.leaf_nodes())
    if len(L) < a.min_cells:
        return
    on, below = counts(L, carriers_f)
    on0, below0 = counts(L, carriers_all)
    uniq = "" if size_count[len(L)] == 1 else "  (NOT UNIQUE)"
    label = "ROOT (trunk)" if nd is t.seed_node else f"{len(L)} cells ({len(L) / N:.0%})"
    print("    " * depth + f"{label}: {on} + {below} = {on + below}  [{on0 + below0}]{uniq}")
    if depth < a.depth:
        for ch in sorted(nd.child_nodes(), key=lambda c: -len(c.leaf_nodes())):
            show(ch, depth + 1)


show(t.seed_node, 0)
