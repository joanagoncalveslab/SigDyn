# =========================
# ASSIGN MUTATIONAL SIGNATURES
# =========================


import sys
log = open(snakemake.log[0], "w")
sys.stdout = log
sys.stderr = log


import os
import SigProfilerAssignment as spa
from SigProfilerAssignment import Analyzer as Analyze

# -------------------------
# Snakemake parameters
# -------------------------
VCF_ALL_DIR = snakemake.input.global_vcf_dir
VCF_PER_CELL_DIR = snakemake.input.per_cell_vcf_dir

OUT_ALL = snakemake.output.global_out
OUT_PER_CELL = snakemake.output.per_cell_out

GENOME_BUILD = snakemake.params.genome_build
CONTEXT = snakemake.params.context
COSMIC_VERSION = snakemake.params.cosmic_version

os.makedirs(OUT_ALL, exist_ok=True)
os.makedirs(OUT_PER_CELL, exist_ok=True)

# -------------------------
# GLOBAL SIGNATURE ASSIGNMENT
# -------------------------
print("Running SigProfilerAssignment on all nodes together...")

Analyze.cosmic_fit(
    samples=VCF_ALL_DIR,
    output=OUT_ALL,
    input_type="vcf",
    context_type=CONTEXT,
    genome_build=GENOME_BUILD,
    cosmic_version=COSMIC_VERSION
)

# -------------------------
# PER-CELL / PER-NODE ASSIGNMENT
# -------------------------
print("Running SigProfilerAssignment per node/cell...")

def run_per_cell_assignment():
    Analyze.cosmic_fit(
        samples=VCF_PER_CELL_DIR,
        output=OUT_PER_CELL,
        input_type="vcf",
        context_type=CONTEXT,
        genome_build=GENOME_BUILD,
        cosmic_version=COSMIC_VERSION
    )

# SigProfilerAssignment sometimes throws a permission error on first run
# (known behavior when creating reference files).
# We retry once automatically.

try:
    run_per_cell_assignment()
except PermissionError as e:
    print("PermissionError encountered. Retrying once...")
    run_per_cell_assignment()

print("Signature assignment completed.")
