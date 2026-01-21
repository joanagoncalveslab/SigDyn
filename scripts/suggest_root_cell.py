#!/usr/bin/env python3
import os
import pandas as pd

GENOTYPES_FILE = snakemake.input[0]
OUT_FILE = snakemake.output[0]

df = pd.read_csv(GENOTYPES_FILE, sep="\t")

# Ensure first column is Site; tolerate weird headers
if "Site" not in df.columns:
    df = df.rename(columns={df.columns[0]: "Site"})

# Cell columns are everything except Site
cell_cols = [c for c in df.columns if c != "Site"]
if len(cell_cols) == 0:
    raise ValueError("No cell columns found in genotype matrix (expected columns besides 'Site').")

# Count RR (R) per cell
r_counts = (df[cell_cols] == "R").sum(axis=0)

# Choose cell with maximum RR count
root_cell = r_counts.idxmax()
root_rr = int(r_counts.loc[root_cell])

# Write output
os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)
with open(OUT_FILE, "w") as f:
    f.write(f"{root_cell}\n")

print(f"Suggested root tip (max RR): {root_cell} (RR count = {root_rr})")
print(f"Wrote root suggestion to: {OUT_FILE}")
