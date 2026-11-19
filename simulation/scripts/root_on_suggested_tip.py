"""
root_on_suggested_tip.py
Automates SigDyn's manual FigTree step: root results/trees/cell_tree.unrooted on
the tip SigDyn suggested (results/trees/root_cell.txt) and save it as
results/trees/cell_tree.rooted (branch lengths kept, no internal labels).
Run from inside a run directory.
"""
from ete3 import Tree

root = open("results/trees/root_cell.txt").read().strip()
t = Tree("results/trees/cell_tree.unrooted", format=0)
t.set_outgroup(t & root)
t.write(outfile="results/trees/cell_tree.rooted", format=5)
print(f"rooted on {root}")
