"""
run_cellcoal_binary.py
======================
Step 1 only: run CellCoal in clean binary mode and emit exactly two files —
the ground-truth tree and the genotype matrix. No signatures, no noise; those
come in later steps.

What it does:
  * writes a minimal CellCoal parameter file (binary alphabet, mutations
    accumulated along the coalescent genealogy at rate proportional to branch
    length; ADO / genotyping error / reads all OFF, so genotypes are the TRUE
    ones),
  * runs the CellCoal binary,
  * parses CellCoal's `snv_gen` into a cells x sites binary matrix
    (00 -> 0 reference, 01/10 -> 1 mutation present) and copies the `trees`
    Newick out.

Requires a compiled CellCoal binary (github.com/dapogon/cellcoal -> `make`);
pass its path with --cellcoal.

Outputs:  <outdir>/genotypes.csv   (rows = cells, cols = SNV site positions, 0/1)
          <outdir>/true_tree.nwk
"""

from __future__ import annotations
import argparse
import os
import subprocess

import pandas as pd
import dendropy


def write_params(path, n_cells, n_sites, mut_rate, growth, Ne, seed, results_dir):
    params = f"""[CellCoal binary: clean true tree + genotypes, no noise]

[COALESCENT]
	[replicates] n1
	[cells] s{n_cells}
	[sites] l{n_sites}
	[effective population size] e{Ne}
	[exponential growth rate] g{growth}

[GENEALOGY MODIFIERS]
	[root branch length ratio] k0.1

[MUTATION MODEL]
	[alphabet binary:0 DNA:1] b0
	[somatic mutation rate] u{mut_rate}

[OUTPUT]
	[print SNV genotypes] 1
	[print trees] 6
	[results folder] o{results_dir}
	[seed] #{seed}
"""
    with open(path, "w") as fh:
        fh.write(params)


def parse_snv_gen(path):
    """CellCoal binary snv_gen -> cells x sites DataFrame of {0,1}."""
    with open(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip() != ""]
    # lines: [rep], [#1], "<ntax> <nsites>", "<site positions>", then taxa rows
    positions = lines[3].split()
    rows, labels = [], []
    for ln in lines[4:]:
        parts = ln.split()
        labels.append(parts[0])
        genos = parts[1:]
        # diploid binary: '00' ref -> 0 ; '01'/'10' het mutation -> 1
        rows.append([1 if "1" in g else 0 for g in genos])
    df = pd.DataFrame(rows, index=labels, columns=positions[: len(rows[0])])
    df.index.name = "cell"
    return df


def read_tree(path):
    with open(path) as fh:
        for ln in fh:
            ln = ln.strip()
            if ln.startswith("(") and ln.endswith(";"):
                return ln
    raise ValueError(f"no Newick tree found in {path}")


def main():
    ap = argparse.ArgumentParser(description="Run CellCoal (binary) -> tree + genotypes.")
    ap.add_argument("--cellcoal", required=True, help="path to compiled cellcoal binary")
    ap.add_argument("--n-cells", type=int, default=100)
    ap.add_argument("--n-sites", type=int, default=1000)
    ap.add_argument("--mut-rate", default="2e-6")
    ap.add_argument("--growth", default="1e-4")
    ap.add_argument("--Ne", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--outdir", default="cellcoal_out")
    ap.add_argument("--keep-outgroup", action="store_true",
                    help="keep the healthy 'outgcell' outgroup (dropped by default)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    results_dir = os.path.join(args.outdir, "_cellcoal_raw")
    params_path = os.path.join(args.outdir, "params_binary")
    write_params(params_path, args.n_cells, args.n_sites, args.mut_rate,
                 args.growth, args.Ne, args.seed, results_dir)

    #subprocess.run([os.path.abspath(args.cellcoal), f"-F{params_path}"],
    #               check=True, capture_output=True, text=True)
    

    cmd = [os.path.abspath(args.cellcoal), f"-F{params_path}"]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True
    )

    print(result.stdout)

    tree_file = os.path.join(results_dir, "trees")
    snv_file = os.path.join(results_dir, "snv_gen")

    if not os.path.exists(tree_file):
        print(result.stderr)
        raise RuntimeError(f"Missing tree file: {tree_file}")

    if not os.path.exists(snv_file):
        print(result.stderr)
        raise RuntimeError(f"Missing genotype file: {snv_file}")


    geno = parse_snv_gen(os.path.join(results_dir, "snv_gen"))
    newick = read_tree(os.path.join(results_dir, "trees"))

    if not args.keep_outgroup:
        # CellCoal's default outgroup tip is 'outgcell' (a healthy reference)
        tns = dendropy.TaxonNamespace()
        t = dendropy.Tree.get(data=newick, schema="newick",
                              preserve_underscores=True, taxon_namespace=tns)
        t.prune_taxa_with_labels(["outgcell"])
        newick = t.as_string(schema="newick").strip()
        geno = geno.drop(index=[i for i in geno.index if i == "outgcell"])
        # sites that were mutated only in the outgroup are now all-zero: drop them
        geno = geno.loc[:, geno.sum(axis=0) > 0]

    geno.to_csv(os.path.join(args.outdir, "genotypes.csv"))
    with open(os.path.join(args.outdir, "true_tree.nwk"), "w") as fh:
        fh.write(newick + "\n")

    print(f"genotypes: {geno.shape[0]} cells x {geno.shape[1]} sites"
          f"  (density {geno.to_numpy().mean():.3f})")
    print(f"wrote {args.outdir}/genotypes.csv and {args.outdir}/true_tree.nwk")


if __name__ == "__main__":
    main()
