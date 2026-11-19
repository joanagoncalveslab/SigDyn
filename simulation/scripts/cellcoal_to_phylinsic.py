"""
cellcoal_to_phylinsic.py
========================
CellCoal binary genotypes + per-site SBS96 typing  ->  PhylinSic-style table for
SigDyn, with every mutation placed at a REAL GRCh38 locus whose trinucleotide
context matches its assigned SBS96 channel, so SigProfiler re-derives the intended
spectrum.

Inputs
  genotypes.csv       cells x sites, {0,1}, columns = integer site IDs
  mutation_types.csv  per site: site,n_cells,clone,signature,sbs96,sub_type,context

Output
  <out.txt>           sites x cells, {R,H}, tab-sep, Site = chr_pos_ref_alt
  --mapping-out CSV   site_id <-> original site, clone, signature, sbs96 (optional)

Genome (for meaningful positions)
  --chrom chr22       single chromosome to place all mutations on. SBS96 is
                      chromosome-independent, so one small chromosome is sufficient
                      and fast. If --reference is omitted, this chromosome is
                      downloaded from UCSC hg38 (~11 MB).
  --reference FASTA   local GRCh38 FASTA that contains --chrom (e.g. the one
                      SigProfiler installed). Skips the download.
  --chr-name NAME     chromosome label written into the site IDs (default: --chrom).
                      Match this to what your VCF/SigProfiler step expects
                      (e.g. '22' vs 'chr22'), or mutations get dropped as not-found.

CONSISTENCY: the FASTA sequence must match the GRCh38 SigProfiler uses (any standard
primary assembly does) and --chr-name must match SigProfiler's expected naming.

Usage
  # download chr22 automatically:
  python cellcoal_to_phylinsic.py genotypes.csv mutation_types.csv out.txt \
         --chrom chr22 --mapping-out map.csv
  # or use a local genome you already have:
  python cellcoal_to_phylinsic.py genotypes.csv mutation_types.csv out.txt \
         --reference /path/GRCh38.fa --chrom chr22 --chr-name 22 --mapping-out map.csv
"""
from __future__ import annotations
import argparse
import gzip
import os
import random
import re
import shutil
import urllib.request
from collections import defaultdict

import pandas as pd
from pyfaidx import Fasta

UCSC = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/{c}.fa.gz"


def load_inputs(geno_csv, types_csv):
    g = pd.read_csv(geno_csv, index_col=0)
    g.columns = g.columns.astype(int)
    m = pd.read_csv(types_csv, index_col=0)
    missing = set(g.columns) - set(m.index)
    if missing:
        raise ValueError(f"{len(missing)} genotype sites absent from mutation_types.csv")
    return g, m.loc[list(g.columns)]


def fetch_chrom(chrom, dest_dir):
    """Download one UCSC hg38 chromosome FASTA (runs on the user's machine)."""
    os.makedirs(dest_dir, exist_ok=True)
    fa = os.path.join(dest_dir, f"{chrom}.fa")
    if not os.path.exists(fa):
        gz = fa + ".gz"
        print(f"downloading {UCSC.format(c=chrom)} ...")
        urllib.request.urlretrieve(UCSC.format(c=chrom), gz)
        with gzip.open(gz, "rb") as fi, open(fa, "wb") as fo:
            shutil.copyfileobj(fi, fo)
        os.remove(gz)
    return fa


def context_pools(seq, contexts_needed, cap):
    """trinucleotide -> list of 1-based +strand centre positions (overlapping search)."""
    return {tri: [mm.start() + 2 for mm in re.finditer(f"(?={tri})", seq)][:cap]
            for tri in contexts_needed}


def place_on_genome(m, fasta_path, chrom, chrom_label=None, seed=1):
    rng = random.Random(seed)
    fa = Fasta(fasta_path, sequence_always_upper=True)
    if chrom not in fa:
        raise ValueError(f"{chrom!r} not in {fasta_path} (have {list(fa.keys())[:3]}...)")
    label = chrom_label or chrom
    seq = str(fa[chrom])

    need = defaultdict(int)
    for _, r in m.iterrows():
        need[r["context"]] += 1
    pools = context_pools(seq, set(need), cap=max(need.values()) * 200)

    ids, used = {}, set()
    for site, r in m.iterrows():
        tri = r["context"]
        ref, alt = r["sub_type"].split(">")
        cand = [p for p in pools[tri] if p not in used]
        if not cand:
            raise RuntimeError(f"ran out of {tri} loci on {chrom}; pick a larger chromosome")
        pos = rng.choice(cand)
        used.add(pos)
        ids[site] = f"{label}_{pos}_{ref}_{alt}"
    return ids


def convert(geno_csv, types_csv, out_txt, reference=None, chrom="chr22",
            chr_name=None, cache_dir="resources/ref_cache", seed=1,
            mapping_out=None, het_symbol="H"):
    g, m = load_inputs(geno_csv, types_csv)

    if reference is None:                      # no local genome -> fetch this one chrom
        reference = fetch_chrom(chrom, cache_dir)

    site_id = place_on_genome(m, reference, chrom, chrom_label=chr_name, seed=seed)

    mat = g.T.replace({0: "R", 1: het_symbol})
    mat.index = [site_id[s] for s in mat.index]
    mat.index.name = "Site"
    mat.to_csv(out_txt, sep="\t")

    if mapping_out:
        mp = m[["clone", "signature", "sbs96", "sub_type", "context"]].copy()
        mp.insert(0, "site_id", [site_id[s] for s in mp.index])
        mp.to_csv(mapping_out)
    return mat


if __name__ == "__main__":
    if "snakemake" in globals():
        convert(snakemake.input["genotypes"], snakemake.input["mutation_types"],   # noqa: F821
                snakemake.output[0],
                reference=snakemake.params.get("reference", None),                  # noqa: F821
                chrom=snakemake.params.get("chrom", "chr22"),                       # noqa: F821
                chr_name=snakemake.params.get("chr_name", None),                    # noqa: F821
                mapping_out=snakemake.output.get("mapping", None))                  # noqa: F821
    else:
        ap = argparse.ArgumentParser()
        ap.add_argument("genotypes")
        ap.add_argument("mutation_types")
        ap.add_argument("out_txt")
        ap.add_argument("--reference", default=None,
                        help="local GRCh38 FASTA containing --chrom; if omitted --chrom is downloaded")
        ap.add_argument("--chrom", default="chr22",
                        help="single chromosome to place all mutations on (default chr22)")
        ap.add_argument("--chr-name", default=None,
                        help="chrom label in site IDs (match SigProfiler, e.g. '22' vs 'chr22')")
        ap.add_argument("--mapping-out", default=None)
        ap.add_argument("--seed", type=int, default=1)
        a = ap.parse_args()
        convert(a.genotypes, a.mutation_types, a.out_txt, reference=a.reference,
                chrom=a.chrom, chr_name=a.chr_name, seed=a.seed, mapping_out=a.mapping_out)
