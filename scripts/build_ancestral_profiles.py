#!/usr/bin/env python3
"""
Build mutation profiles for all nodes (leaf + internal) of a rooted cell phylogeny.

For each node, compute mutations shared across all descendant leaves
under the infinite-sites assumption.

Inputs:
- results/trees/cell_tree.rooted  (manual)
- filtered genotype matrix (after PON removal)

Outputs:
- genotype matrix for all nodes (leaf + internal)
- rooted tree with explicit node names
"""

import sys

log = open(snakemake.log[0], "w")
sys.stdout = log
sys.stderr = log


from ete3 import Tree
import pandas as pd
import random
import string
import os

########################################
# Snakemake I/O
########################################

TREE_FILE = "results/trees/cell_tree.rooted"
GENOTYPES_FILE = snakemake.input["genotypes"]

OUT_GENOTYPES = snakemake.output["genotypes_all_nodes"]
OUT_TREE = snakemake.output["tree_named"]

########################################
# Validate inputs
########################################

if not os.path.exists(TREE_FILE):
    raise FileNotFoundError(
        f"Rooted tree not found at {TREE_FILE}. "
        "Create it manually (FigTree) and save as results/trees/cell_tree.rooted."
    )

mutations = pd.read_csv(GENOTYPES_FILE, sep="\t")
if "Site" not in mutations.columns:
    # tolerate first-column-is-site
    mutations = mutations.rename(columns={mutations.columns[0]: "Site"})
mutations.set_index("Site", inplace=True)

# Genotype columns are cell barcodes
genotype_cells = set(mutations.columns)

########################################
# Load tree
########################################

tree = Tree(TREE_FILE, format=1)

########################################
# Optional leaf name normalization
########################################
# If your tree leaves are not exactly the genotype column names,
# you can normalize them here.
#
# Common pattern: tree leaf name contains a long prefix and ends with barcode.
# Example: "..._CELLBARCODE" -> want "CELLBARCODE"
#
# You can control behavior by setting in config.yaml (optional):
#   tree_leaf_name:
#     split_on: "_"     # take last token after splitting on this string
#
split_on = None
try:
    split_on = snakemake.config.get("tree_leaf_name", {}).get("split_on", None)
except Exception:
    split_on = None

def normalize_leaf_name(name: str) -> str:
    if split_on and split_on in name:
        return name.split(split_on)[-1]
    return name

# Build mapping from original leaf names to normalized names
leaf_name_map = {}
for leaf in tree.iter_leaves():
    leaf_name_map[leaf.name] = normalize_leaf_name(leaf.name)

# Check mismatch
normalized_leaves = {leaf_name_map[n] for n in leaf_name_map}
missing = sorted(list(normalized_leaves - genotype_cells))

if missing:
    # show only first few to avoid huge spam
    preview = "\n".join(missing[:20])
    raise KeyError(
        "Leaf names in rooted tree do not match genotype matrix columns.\n"
        f"Example missing (first 20):\n{preview}\n\n"
        "Fix options:\n"
        "  - Ensure FigTree saves leaf names unchanged, OR\n"
        "  - Add config.yaml:\n"
        "      tree_leaf_name:\n"
        "        split_on: \"_\"\n"
        "    (or another delimiter to extract the barcode)\n"
    )

########################################
# Assign unique names to internal nodes
########################################

existing_names = set(node.name for node in tree.traverse() if node.name)

def generate_unique_name(existing, length=7):
    while True:
        name = "".join(random.choices(string.ascii_letters, k=length))
        if name not in existing:
            return name

for node in tree.traverse("preorder"):
    if not node.name:
        new_name = generate_unique_name(existing_names)
        node.name = new_name
        existing_names.add(new_name)

########################################
# Export tree with node names
########################################

tree.write(outfile=OUT_TREE, format=1)

########################################
# Compute shared mutations per node
########################################

def leaf_mutations_for_cell(cell_name: str):
    """Return set of mutation sites for a cell (H or A)."""
    s = mutations[cell_name]
    return set(s.index[s.isin(["H", "A"])])

def get_mutations(node):
    """Mutations present in ALL descendant leaves (intersection)."""
    if node.is_leaf():
        cell = leaf_name_map[node.name]  # normalized
        return leaf_mutations_for_cell(cell)

    child_sets = [get_mutations(child) for child in node.children]
    return set.intersection(*child_sets) if child_sets else set()

dict_nodes_mutations = {}
for node in tree.traverse("preorder"):
    dict_nodes_mutations[node.name] = get_mutations(node)

########################################
# Convert node → mutations to genotype matrix
########################################

all_mutations = sorted({m for muts in dict_nodes_mutations.values() for m in muts})

df_nodes = pd.DataFrame(index=all_mutations, columns=dict_nodes_mutations.keys())
for node, muts in dict_nodes_mutations.items():
    df_nodes[node] = ["H" if mut in muts else "R" for mut in df_nodes.index]

df_nodes = df_nodes.reset_index().rename(columns={"index": "Site"})

########################################
# Write output
########################################

os.makedirs(os.path.dirname(OUT_GENOTYPES), exist_ok=True)
df_nodes.to_csv(OUT_GENOTYPES, sep="\t", index=False)
