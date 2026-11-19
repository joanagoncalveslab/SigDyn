#!/usr/bin/env python3
"""
Root the unrooted FastTree tree of one level on a cell YOU choose
(replaces the manual FigTree step of SigDyn).

Reads   <out>/<name>/pXX/results/trees/cell_tree.unrooted
Writes  <out>/<name>/pXX/results/trees/cell_tree.rooted
in the same form as a FigTree export (leaf names + branch lengths, no support labels),
which is what SigDyn's build_ancestral_profiles.py expects.

  python root_trees.py --name neck --show                              # suggestion of every level
  python root_trees.py --name neck --level 45 --cell neck15_XXXX-1
  python root_trees.py --name neck --level 45 --cell suggested         # the level's root_cell.txt
  python root_trees.py --name neck --level 45 --cell neck15_XXXX-1 --force   # redo
"""
import argparse
from pathlib import Path

from ete3 import Tree

HERE = Path(__file__).resolve().parent


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--name", required=True)
    ap.add_argument("--cell", help="root cell name, or 'suggested' to use the level's root_cell.txt")
    ap.add_argument("--level", help="the level to root, e.g. 45 or p45")
    ap.add_argument("--show", action="store_true", help="only print the suggested root cells")
    ap.add_argument("--force", action="store_true", help="overwrite an existing cell_tree.rooted")
    ap.add_argument("--out", default=str(HERE / "levels"))
    a = ap.parse_args()

    out = Path(a.out) / a.name
    dirs = sorted(p for p in out.glob("p*") if (p / "results/trees/cell_tree.unrooted").exists())
    if a.cell and not a.level:
        ap.error("rooting needs --level (one level at a time)")
    if a.level:
        dirs = [d for d in dirs if d.name == f"p{int(a.level.lstrip('p')):02d}"]
    if not dirs:
        raise SystemExit(f"no unrooted tree for this level under {out} - run run_levels.py stage1 first")

    for d in dirs:
        trees = d / "results/trees"
        suggested = (trees / "root_cell.txt").read_text().strip() if (trees / "root_cell.txt").exists() else None
        if a.show or not a.cell:
            print(f"{d.name}: suggested root (max RR) = {suggested}")
            continue
        rooted = trees / "cell_tree.rooted"
        if rooted.exists() and not a.force:
            print(f"{d.name}: already rooted (use --force to redo; then rerun stage2)")
            continue
        cell = suggested if a.cell == "suggested" else a.cell
        t = Tree(str(trees / "cell_tree.unrooted"), format=0)   # FastTree support read as support
        hits = t.search_nodes(name=cell)
        if not hits:
            raise SystemExit(f"{d.name}: cell {cell!r} not in the tree")
        t.set_outgroup(hits[0])
        t.write(format=5, outfile=str(rooted))                  # leaf names + branch lengths
        print(f"{d.name}: rooted on {cell} -> {rooted}")

    if not a.cell and not a.show:
        print("\nChoose a cell: python root_trees.py --name <name> --level <XX> --cell <CELL>")


if __name__ == "__main__":
    main()
