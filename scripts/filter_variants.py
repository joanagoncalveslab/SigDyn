
import sys
sys.stderr = open(snakemake.log[0], "w")
sys.stdout = open(snakemake.log[0], "a")


"""
Filter germline variants from a PhylinSic genotype matrix using a Panel of Normals (PON).

This script removes variant sites that are present in the PON,
indicating likely germline variants or recurrent sequencing artifacts.

Input:
  - Genotype matrix from PhylinSic (txt)
    Rows: variant sites (chrom_pos)
    Columns: cells
    First column must be named "Site"

Output:
  - Filtered genotype matrix (txt)
"""

import pandas as pd
import yaml
from collections import defaultdict

########################################
# Snakemake I/O
########################################

genotype_file = snakemake.input[0]
output_file = snakemake.output[0]

########################################
# Load config
########################################

with open("config/config.yaml") as f:
    config = yaml.safe_load(f)

pon_vcf = config["filtering"]["pon_vcf"]
use_pon = config["filtering"]["use_pon"]

########################################
# Helper: load PON variants
########################################

def load_pon_variants(vcf_file):
    """
    Load PON variants into a dictionary:
    {chromosome: set(positions)}
    """
    pon = defaultdict(set)

    with open(vcf_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chrom, pos = fields[0], fields[1]
            pon[chrom].add(pos)

    return pon

########################################
# Main logic
########################################

# Load genotype matrix
genotypes = pd.read_csv(genotype_file, sep="\t")
genotypes = genotypes.set_index("Site")

if use_pon:
    pon_dict = load_pon_variants(pon_vcf)

    sites_to_remove = []

    for site in genotypes.index:
        chrom, pos = site.split("_", 1)

        if chrom in pon_dict and pos in pon_dict[chrom]:
            sites_to_remove.append(site)

    genotypes_filtered = genotypes.drop(index=sites_to_remove)

else:
    genotypes_filtered = genotypes

########################################
# Write output
########################################

genotypes_filtered.reset_index().to_csv(
    output_file,
    sep="\t",
    index=False
)
