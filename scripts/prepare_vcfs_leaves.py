#!/usr/bin/env python3
import os
import pandas as pd

# -------------------------
# Snakemake I/O
# -------------------------
GENOTYPES_FILE = snakemake.input["genotypes"]

GLOBAL_DIR = snakemake.output["global_dir"]     # e.g. results/vcf_sig_leaves
PER_CELL_DIR = snakemake.output["per_cell_dir"] # e.g. results/vcf_per_cell_leaves

GLOBAL_INPUT_DIR = os.path.join(GLOBAL_DIR, "input")
PER_CELL_INPUT_DIR = os.path.join(PER_CELL_DIR, "input")

os.makedirs(GLOBAL_INPUT_DIR, exist_ok=True)
os.makedirs(PER_CELL_INPUT_DIR, exist_ok=True)

GLOBAL_VCF = os.path.join(GLOBAL_INPUT_DIR, "vcf_sig.vcf")

# -------------------------
# LOAD DATA
# -------------------------
df = pd.read_csv(GENOTYPES_FILE, sep="\t")

# Ensure first column is Site
if "Site" not in df.columns:
    df = df.rename(columns={df.columns[0]: "Site"})

mutation_col = "Site"
nodes = [c for c in df.columns if c != mutation_col]

if len(nodes) == 0:
    raise ValueError("No sample columns found in genotype matrix (expected columns besides 'Site').")

# -------------------------
# GLOBAL VCF (all mutations present anywhere)
# -------------------------
mut_info_vcf = []
for _, row in df.iterrows():
    mutation = row[mutation_col]
    if any(val in ["A", "H"] for val in row[nodes]):
        mut_info_vcf.append(mutation)

with open(GLOBAL_VCF, "w") as f:
    f.write("#CHROM\tPOS\tFILTER\tREF\tALT\n")
    for mut in mut_info_vcf:
        chrom, pos, ref, alt = mut.split("_")
        chrom = chrom.replace("chr", "")
        f.write(f"{chrom}\t{pos}\tPASS\t{ref}\t{alt}\n")

# -------------------------
# PER-CELL / PER-NODE VCFs
# -------------------------
dict_mut_vcf = {node: [] for node in nodes}

for _, row in df.iterrows():
    mutation = row[mutation_col]
    for node in nodes:
        if row[node] in ["A", "H"]:
            dict_mut_vcf[node].append(mutation)

for node, mutations in dict_mut_vcf.items():
    out_vcf = os.path.join(PER_CELL_INPUT_DIR, f"{node}.vcf")
    with open(out_vcf, "w") as f:
        f.write("#CHROM\tPOS\tFILTER\tREF\tALT\n")
        for mut in mutations:
            chrom, pos, ref, alt = mut.split("_")
            chrom = chrom.replace("chr", "")
            f.write(f"{chrom}\t{pos}\tPASS\t{ref}\t{alt}\n")

print(f"Generated global VCF: {GLOBAL_VCF}")
print(f"Generated {len(nodes)} per-node VCF files in {PER_CELL_INPUT_DIR}")
