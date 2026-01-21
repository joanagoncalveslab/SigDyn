# =========================
# ASSIGN ACTIVE SIGNATURES USING NNLS
# =========================

import pandas as pd
import numpy as np
from scipy.optimize import nnls
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen

# -------------------------
# Inputs
# -------------------------
VCF_DIR = snakemake.input.vcf_dir
CUSTOM_COSMIC = snakemake.input.custom_cosmic
OUT_FILE = snakemake.output[0]
GENOME = snakemake.params.genome

# -------------------------
# GENERATE MUTATION MATRICES
# -------------------------
print("Generating trinucleotide matrices...")

output = matGen.SigProfilerMatrixGeneratorFunc(
    project="mutation_contexts_active",
    reference_genome=GENOME,
    path_to_input_files=VCF_DIR,
    exome=False
)

mutation_matrix = output["96"]

# -------------------------
# LOAD CUSTOM COSMIC SIGNATURES
# -------------------------
custom_signatures = pd.read_csv(CUSTOM_COSMIC, sep="\t", index_col=0)

# Align contexts
mutation_matrix = mutation_matrix.loc[custom_signatures.index]

# -------------------------
# NNLS FITTING
# -------------------------
results = {}

A = custom_signatures.values

for sample in mutation_matrix.columns:
    y = mutation_matrix[sample].values
    exposures, residual = nnls(A, y)
    results[sample] = exposures

# -------------------------
# OUTPUT
# -------------------------
exposure_df = pd.DataFrame(
    results,
    index=custom_signatures.columns
)

exposure_df.to_csv(OUT_FILE, sep="\t")

print(f"Active-signature exposures written to: {OUT_FILE}")
