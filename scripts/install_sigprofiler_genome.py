#!/usr/bin/env python3
import os
from SigProfilerMatrixGenerator import install as geninstall

# Install GRCh38 reference genome for SigProfiler (downloads files)
geninstall.install("GRCh38")

# Create sentinel file
out = snakemake.output[0]
os.makedirs(os.path.dirname(out), exist_ok=True)
with open(out, "w") as f:
    f.write("GRCh38 installed\n")

print("SigProfiler genome GRCh38 installed.")
