# Missingness levels for SigDyn

Runs **your SigDyn snakefile, unchanged** on genotype matrices with increasing missingness
(0, 15, 45, 65 % and the true missingness), then collects tree support, region entropy,
signature presence and MMR–SBS26 correlations over the levels.

## Setup

Copy this folder into your SigDyn folder, so it sits next to the `snakefile`:

```
SigDyn-main/
  snakefile, scripts/, envs/, config/, resources/, .snakemake/conda/   (yours, untouched)
  filtered_feature_bc_matrix/          # breast 10x matrix (already there; used by --mmr)
  missingness/
    make_levels.py  run_levels.py  root_trees.py  analyze.py  README.md
    inputs/                            # <- you put the inputs here (see below)
      breast/
        genotypes_no_PON.txt
        ref_count.all.txt
        alt_count.all.txt
      neck/
        genotypes_no_PON.txt
        ref_count.all.txt
        alt_count.all.txt
    levels/                            # created by the scripts: one folder per level + analysis/
```

You need:
- `snakemake` on your PATH (the one you use for SigDyn).
- A Python with pandas, numpy, ete3, dendropy, matplotlib and scanpy, such as your `snakemake` conda env:
  `/Users/saracosta/miniforge3/envs/snakemake/bin/python`.

Run all commands from `SigDyn-main/missingness/`.

**If this is a new or copied SigDyn folder:** conda envs cannot be copied, so snakemake builds fresh ones
in `SigDyn-main/.snakemake/conda/` on the first run. The new sigprofiler env has no GRCh38 genome yet,
but a copied `resources/sigprofiler_genomes/GRCh38.installed` makes SigDyn skip its install rule.
SigProfiler then tries to download the genome by itself during stage 2, and with little disk space
this leaves a corrupt copy ("MD5 Mismatch ... genome GRCh38 has not been installed"). Either delete
`resources/sigprofiler_genomes/GRCh38.installed` before the first run (SigDyn's own rule then installs
it; needs ~3 GB free), or copy a working genome into the new env:

```bash
G=lib/python3.10/site-packages/SigProfilerMatrixGenerator/references/chromosomes/tsb/GRCh38
cp -Rp <old SigDyn>/.snakemake/conda/<old sigprofiler env>/$G <new SigDyn>/.snakemake/conda/<new sigprofiler env>/$G
```

## Inputs

For each dataset:
- **genotypes**: imputed PhylinSic genotypes, **already without PON sites** (every entry H/R/A);
- **ref_count.all.txt** and **alt_count.all.txt**: from the **same** PhylinSic run. An entry with
  no coverage (empty in both files) is "truly missing".

**Where to put them:** `missingness/inputs/<name>/`, using exactly these file names
(`<name>` is what you pass as `--name`):

| put this file ... | ... at | comes from (paper data) |
|---|---|---|
| breast genotypes (460 sites × 688 cells) | `inputs/breast/genotypes_no_PON.txt` | `missing_values/3_09/genotypes_no_PON.txt` |
| breast ref counts | `inputs/breast/ref_count.all.txt` | `missing_values/3_09/ref_count.all.txt` |
| breast alt counts | `inputs/breast/alt_count.all.txt` | `missing_values/3_09/alt_count.all.txt` |
| neck genotypes (731 sites × 56,145 cells) | `inputs/neck/genotypes_no_PON.txt` | `filtered_genotype.txt` (rename it) |
| neck ref counts | `inputs/neck/ref_count.all.txt` | `neck_2025_6/output/08_genotype/ref_count.all.txt` |
| neck alt counts | `inputs/neck/alt_count.all.txt` | `neck_2025_6/output/08_genotype/alt_count.all.txt` |

(`missing_values/neck_paper/` holds copies of the three neck files. `neck/genotypes.txt` is a
different 1000-site run, and is **not** the paper's neck data.)

The count files may contain more sites or cells than the genotypes; only the genotype sites and
cells are used. Stage 1 stops with an error if a site or cell is missing from the counts.

To keep the files somewhere else, pass the paths instead: `--geno`, `--ref`, `--alt`.

## How the levels are made (`make_levels.py`)

- Real (covered) entries always keep their genotype.
- Level X % blanks a random subset of the truly-missing entries so that X % of the whole matrix
  is empty. The non-blanked entries keep their imputed letter.
- Levels are **nested** (15 ⊂ 45 ⊂ 65 ⊂ top), from one random draw (`--seed`, default 0).
- The top level blanks **all** truly-missing entries and is named by its real percentage
  (breast `p86`, which equals `3_09/genotypes_masked_no_PON.txt`; neck `p90`).

Each level becomes a SigDyn working directory:

```
levels/<name>/levels_summary.csv           # real % missing per level
levels/<name>/pXX/                          # snakemake --directory for this level
    config, resources                       # symlinks to your SigDyn config/ and resources/
    results/genotypes/genotypes_no_PON.txt  # the level's genotypes (SigDyn input)
    results/...                             # normal SigDyn outputs
    snakemake_stage1_unrooted.log, snakemake_stage2_full.log
```

## Steps

Every run is **one dataset** (`--name`) and, for SigDyn, **one level** (`--level`). Nothing runs
in parallel: run the levels one after another, in any order, whenever you want.

```bash
PY=/Users/saracosta/miniforge3/envs/snakemake/bin/python
```

**0. Build the levels** (once per dataset, takes a few minutes):

```bash
$PY run_levels.py levels --name breast      # reads inputs/breast/
$PY run_levels.py levels --name neck        # reads inputs/neck/
```

This writes the genotype files of all levels together. They must come from the same random
draw to be nested; no SigDyn is run here. It prints the level names (e.g. `0 15 45 65 86`).

Then, for **each level** (example: neck, 45 %), steps 1 to 3:

**1. SigDyn stage 1** (fasta, FastTree, root suggestion):

```bash
$PY run_levels.py stage1 --name neck --level 45
```

It ends by printing SigDyn's suggested root cell (max RR) for this level.

**2. Root the tree** (replaces the manual FigTree step). You choose the cell:

```bash
$PY root_trees.py --name neck --show                                 # suggestion of every level done so far
$PY root_trees.py --name neck --level 45 --cell <CELL>
$PY root_trees.py --name neck --level 45 --cell suggested            # this level's own suggestion
$PY root_trees.py --name neck --level 45 --cell <CELL> --force       # redo
```

This writes `results/trees/cell_tree.rooted` like a FigTree export: leaf names and branch
lengths, no support labels. If you re-root after stage 2, delete that level's stage-2 outputs
(or the whole level folder and redo stage 1) and run stage 2 again.

**3. SigDyn stage 2** (ancestral profiles, signatures):

```bash
$PY run_levels.py stage2 --name neck --level 45
```

`--level` accepts `45` or `p45`. Outputs: `levels/neck/p45/results/...` (normal SigDyn outputs) and
the logs `levels/neck/p45/snakemake_stage1_unrooted.log` / `snakemake_stage2_full.log`.

**4. Analyses**, per dataset, whenever you like; it uses every level finished so far (from the SigDyn outputs; writes `levels/<name>/analysis/*.csv` and `summary.png/pdf`):

```bash
$PY analyze.py --name neck --entropy --signature SBS32
$PY analyze.py --name breast --mmr ../filtered_feature_bc_matrix
```

| analysis | reads | same as |
|---|---|---|
| tree support (always) | `results/trees/cell_tree.unrooted` | `tree_support/unrooted_support_breast_neck2.ipynb` |
| `--entropy` (K=30) | `results/trees/cell_tree.rooted` | `entropy_neck/region_entropy_scalable.ipynb` |
| `--signature SBS32` | `results/signatures/per_leaf/.../Assignment_Solution_Activities.txt` | `Fig4_b_box_plots.ipynb` |
| `--mmr DIR` | per-leaf activities + 10x matrix | `single_cell_exploration.ipynb` (Spearman/Pearson, MLH1/MSH2/MSH6/PMS2 vs SBS26) |

## How SigDyn is called

For the chosen level, `run_levels.py` runs:

```
snakemake -s ../snakefile --directory levels/<name>/pXX --use-conda --conda-frontend conda \
          --conda-prefix ../.snakemake/conda --cores 4 --rerun-triggers mtime \
          --allowed-rules <every rule except filter_variants and phylinsic_external> \
          stage1_unrooted | stage2_full
```

- **The snakefile and scripts are not modified.** The PON rules are left out with `--allowed-rules`,
  because the genotypes are given already PON-filtered at `results/genotypes/genotypes_no_PON.txt`.
- `--conda-prefix` reuses your existing SigDyn conda envs, so nothing is reinstalled.
- `--directory` keeps each level in its own folder; your `SigDyn-main/results/` is never touched.
- `--cores N` sets snakemake's cores (default 4).

## Testing on a subsample

Add `--cells N` and a separate `--out` to the levels step, then pass the same `--out` to every other step:

```bash
$PY run_levels.py levels --name neck --cells 400 --out test
$PY run_levels.py stage1 --name neck --level 45 --out test
$PY root_trees.py --name neck --level 45 --out test --cell suggested
$PY run_levels.py stage2 --name neck --level 45 --out test
$PY analyze.py --name neck --out test --entropy --signature SBS32
```

Tested this way with 120 breast cells and 400 neck cells (both under 5 minutes). At 0 % the test
reproduced the paper's SBS32 pattern (46.5 % of cells; N 65 %, R 45 %, T 30 %, L 80 %).

## One SigDyn folder for everything

Use **one** SigDyn folder for all datasets and levels; no copies are needed. Each level reads
from it (snakefile, `scripts/`, `envs/`, the conda envs in `.snakemake/conda/`, and `config/` and
`resources/` through symlinks) and writes only into its own `missingness/levels/<name>/pXX/`.
Your `SigDyn-main/results/` is never touched.

**Nothing has to be cleaned or swapped between runs.** Every level has its own folder, and the
`levels` step already put that level's input in it (`pXX/results/genotypes/genotypes_no_PON.txt`).
So you do not empty `SigDyn-main/results/` or replace the input before the next level; just run the
next `--level`.

Do not edit `config/`, `resources/`, the snakefile or `scripts/` while a set of levels is in progress,
because every level uses them as they are at the time it runs. Make a separate SigDyn folder only if
you want to run a modified SigDyn next to this one.

## Outputs saved per level

Everything SigDyn normally writes, inside `levels/<name>/pXX/`:

| path | from | used by analyze.py |
|---|---|---|
| `results/genotypes/genotypes_no_PON.txt` | level input (`levels` step) | |
| `results/trees/mutations.fa`, `cell_tree.unrooted`, `root_cell.txt` | stage 1 | tree support (unrooted) |
| `results/trees/cell_tree.rooted` | `root_trees.py` | entropy |
| `results/trees/.rooting_done`, `cell_tree.rooted.named` | stage 2 | |
| `results/genotypes/genotypes_all_nodes.txt` | stage 2 (ancestral profiles) | |
| `results/vcf_sig_leaves/`, `vcf_per_cell_leaves/` | stage 2 (VCFs + SigProfiler matrices, leaves) | |
| `results/signatures/all_leaves/`, `signatures/per_leaf/` | stage 2 (SigProfilerAssignment) | `per_leaf/.../Assignment_Solution_Activities.txt` → SBS32, MMR |
| `results/signatures/Custom_COSMIC_SBS_v3.4.txt` | stage 2 (active signatures, ≥10 % of cells) | |
| `results/vcf_sig_all/`, `vcf_per_cell_all/` | stage 2 (VCFs, all nodes) | |
| *(neck only)* the 4 `results/vcf_*` folders are **deleted** after stage 2 succeeds; `VCFS_DELETED` lists them | `run_levels.py` | |
| `results/signatures/exposures_active_signatures.txt` | stage 2 (final exposures, all nodes) | |
| `results/logs/`, `snakemake_*.log`, `.snakemake/` | logs and snakemake metadata | |

Plus, per dataset: `levels/<name>/levels_summary.csv` (real % missing per level) and
`levels/<name>/analysis/` (CSV tables + `summary.png/pdf`).

## Disk space

| | per level | 5 levels |
|---|---|---|
| breast | about 30 MB | about 150 MB |
| neck | about 6 GB | about 30 GB |

For neck, about 5.5 GB of each level is the per-cell VCF folders. So **for neck, `run_levels.py` deletes
`results/vcf_sig_leaves/`, `vcf_per_cell_leaves/`, `vcf_sig_all/` and `vcf_per_cell_all/` once stage 2
has finished successfully**. That leaves about 0.5 GB per level. The analyses don't use these
folders, and every other output is kept. The deleted folders are listed in `pXX/VCFS_DELETED`. Breast keeps
everything. To change which datasets are cleaned, edit `DELETE_VCFS_FOR = {"neck"}` at the top
of `run_levels.py`.

Because snakemake would see those folders as missing and redo all of stage 2, `stage2` refuses
to run again on a cleaned level. To really redo it, delete that level's `results/signatures/` and
`VCFS_DELETED` first.

Each neck level still needs the full ~6 GB **while** stage 2 runs. Check the free space first (`df -h ~`).
If it's short, put the levels on another disk with `--out /path/on/big/disk`, added to **every**
command (levels, stage1, root_trees, stage2, analyze).

## Runtime (full data)

- breast: about a minute per level and stage.
- neck, per level: stage 1 about 1 h (FastTree on 56k cells). Stage 2 is several hours, because it
  also builds and assigns signatures for about 112k ancestral nodes.
- Steps already done are skipped by snakemake, so after a crash just rerun the same command.

## DGEA over the levels (Fig 4c/d, neck)

Repeats the paper's DGEA (`Case_studies/Figure4_LSCC`: `bef_dgea.ipynb` → `dgea.R` → `Fig4_c_d_volcano_plot.ipynb`)
for every finished neck level, then compares each level with 0 % (= the paper's run).

- Groups: **L32 = SBS32 = 0** (group1) vs **H32 = SBS32 > 0** (group2), from each level's per-leaf activities.
  At 0 % this is exactly the paper's split (30,474 / 25,671). Cells in the 10x data but not in the tree are
  "other", and stay in the model as in the paper.
- Expression (normalize_total 1e4 + log1p, genes in ≥ 10 % of cells) is the same for every level, so it is built once.
- MAST: `find_de_MAST_RE_simple` copied unchanged from `dgea.R`.
- Volcano/GO as the notebook: significant = p < 0.01 and |log2FC| > 0.5; Enrichr GO_Biological_Process_2021
  on the genes up in L32 (needs internet). Everything after MAST is made **twice**, in subfolders:
  `fdr/` uses the BH-adjusted FDR (volcano y-axis "Adjusted p-value"); `rawp/` uses MAST's raw Pr(>Chisq),
  as the paper's Fig 4c did. The summary also reports FDR < 0.05 counts.

Needs the neck 10x data: a folder with `neck15/ neck16/ neck17/ neck18/`, each with `filtered_feature_bc_matrix/`.
Python: the `dgea` env (`/Users/saracosta/miniforge3/envs/dgea/bin/python`, has scanpy + gseapy). R: `Rscript` with MAST and zellkonverter.

```bash
PY=/Users/saracosta/miniforge3/envs/dgea/bin/python
$PY dgea_prepare.py --name neck --tenx /path/to/neck_10x     # once; a few minutes
Rscript dgea_mast.R --level 0                                 # one level per call; repeat for 30, 60, 90
$PY dgea_plots.py --name neck                                 # whenever; uses every level with mast_results.csv
```

Outputs, in `levels/neck/dgea/`:

| file | content |
|---|---|
| `expression.h5ad`, `group_sizes.csv` | shared expression matrix; L32/H32/other cells per level |
| `pXX/groups.csv`, `pXX/mast_results.csv` | the level's groups; MAST output (as `mast_results_26_sep.csv`) |
| `pXX/{fdr,rawp}/volcano.pdf`, `significant_genes.csv`, `go_enrichment.csv`, `go_dotplot.pdf` | Fig 4c and 4d for that level |
| `{fdr,rawp}/dgea_vs_missingness.csv` | per level: # DE genes each way (at that folder's p < 0.01, and at FDR < 0.05), recall / precision / Jaccard of the 0 % genes, Spearman of log2FC with 0 % |
| `{fdr,rawp}/go_top_terms_vs_missingness.csv` | adjusted p of the 10 top GO terms at 0 %, at every level |
| `{fdr,rawp}/summary.pdf/png` | the above as figures |

Test on a subsample: add `--cells 2000 --out test` to `dgea_prepare.py`, and `--out test` to the other two
(the `test` folder needs the levels' activity files, e.g. from the SigDyn test run).
`dgea_mast.R --cores N` runs MAST's per-gene fits on N cores (same result, more memory).
