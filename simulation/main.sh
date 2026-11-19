#!/usr/bin/env bash
# Biopsy3 simulation, start to end.
#
#   bash main.sh            # everything: simulate -> label -> sigdyn -> evaluate -> plots
#   bash main.sh <step>     # one step only: simulate | label | sigdyn | evaluate | plots
#   bash main.sh scan       # optional, not part of "all": search CellCoal tumours for a
#                           # good biopsy clade (how SEED / TUMOUR_CELLS / CLADE_SIZE were picked)
#
# Order: CellCoal tumour -> clade (all cells) -> label clade with signatures -> sample
#        BIOPSY_CELLS cells -> drop sites in < MIN_FREQ of them -> PhylinSic format
#        -> SigDyn at several per-cell caps x replicates -> leaf evaluation -> figures
#
# Clones: the larger clone (ONSETS[0], 369 biopsy cells) gets SBS26, the smaller one
# (ONSETS[1], 232 biopsy cells) gets SBS5; cells outside both clones are SBS1-only.
#
# Any setting can be overridden per call, e.g.  REPS="1" CAPS="100 50" bash main.sh
set -euo pipefail

# ---------------- settings ----------------
HERE=$(cd "$(dirname "$0")" && pwd)
OUT="${OUT:-$HERE/output}"                        # everything is written here
NAME="${NAME:-sbs26_large_clone}"                # run folder: $OUT/$NAME

SEED="${SEED:-3}"                                 # CellCoal seed
TUMOUR_CELLS="${TUMOUR_CELLS:-2000}"              # CellCoal --n-cells
NSITES="${NSITES:-11600}"                         # CellCoal --n-sites
MUT_RATE="${MUT_RATE:-1.75e-6}"                   # CellCoal --mut-rate
CLADE_SIZE="${CLADE_SIZE:-751}"                   # biopsy clade, exact size; ALL its cells are labelled
BIOPSY_CELLS="${BIOPSY_CELLS:-688}"
MIN_FREQ="${MIN_FREQ:-0.20}"                      # drop sites seen in < 20% of biopsy cells (after labelling)
CLONE_SIGS="${CLONE_SIGS:-SBS26 SBS5}"            # signature of each clone, in ONSETS order
CLONE_RATES="${CLONE_RATES:-2 2}"
ONSETS="${ONSETS:-369 232}"                       # clone start nodes, as cell counts in the biopsy tree map
CAPS="${CAPS:-250 200 150 100 75 50 25}"          # per-cell caps (plus one uncapped run)
REPS="${REPS:-1 2 3}"
CORES="${CORES:-8}"

# ---------------- tools / external inputs ----------------
# Requirements: conda, snakemake (>= 8), and a python with pandas, numpy, dendropy, pyfaidx, matplotlib.
PY="${PY:-python}"
SNAKEMAKE="${SNAKEMAKE:-snakemake}"
# CellCoal 1.1.0 binary (github.com/dapogon/cellcoal). The bundled one was compiled for macOS
# (arm64); on another system compile CellCoal and point CELLCOAL at it.
CELLCOAL="${CELLCOAL:-$HERE/bin/cellcoal-1.1.0}"
COSMIC="${COSMIC:-$HERE/resources/COSMIC_v3.4_SBS_GRCh38.txt}"
# GRCh38 chr22 FASTA used to give each mutation a real trinucleotide context.
# Downloaded once from UCSC hg38 into resources/ if it is not there yet.
REFERENCE="${REFERENCE:-$HERE/resources/chr22.fa}"
# conda envs for SigDyn, built by snakemake --use-conda on the first run
CONDA_PREFIX_DIR="${CONDA_PREFIX_DIR:-$HERE/.snakemake/conda}"
# python of SigDyn's "base" env (has ete3), used for automatic rooting; found automatically if unset
BASE_PY="${BASE_PY:-}"

S=$HERE/scripts
SIGDYN=$HERE/SigDyn-main
export SIM_OUT=$OUT RUN_NAME=$NAME               # evaluate_leaves.py and plot_*.py resolve paths here
mkdir -p "$OUT"
# ------------------------------------------

# 1) CellCoal tumour -> whole clade + 688-cell biopsy (unfiltered) -> biopsy tree map
simulate () {
  local SIM=$OUT/$NAME/sim
  mkdir -p "$SIM"; cd "$SIM"
  $PY $S/run_cellcoal_binary.py --cellcoal "$CELLCOAL" \
    --n-cells $TUMOUR_CELLS --n-sites $NSITES --mut-rate $MUT_RATE --growth 1e-4 --Ne 10000 \
    --seed $SEED --outdir tumour
  [ -s tumour/genotypes.csv ] || { echo "CellCoal produced no genotypes (raise NSITES, lower MUT_RATE)"; exit 1; }
  rm -rf tumour/_cellcoal_raw
  $PY $S/biopsy_sample.py --tumour tumour --n-cells $BIOPSY_CELLS --min-freq 0 \
    --clade-min $CLADE_SIZE --clade-max $CLADE_SIZE --seed 1 --outdir biopsy --clade-outdir clade
  $PY $S/inspect_biopsy_tree.py --dir biopsy --min-freq $MIN_FREQ | tee biopsy/tree_map.txt
}

# 2) label every mutation of the clade, keep the biopsy's sites, rare-site filter, PhylinSic format
label () {
  local SIM=$OUT/$NAME/sim
  [ -f "$SIM/biopsy/genotypes.csv" ] || { echo "run the simulate step first"; exit 1; }
  [ -z "$(ls -A "$OUT/$NAME/runs" 2>/dev/null)" ] || echo "WARNING: $OUT/$NAME/runs already has SigDyn runs made from previous labels"
  cd "$SIM"
  $PY $S/assign_signatures_onsets.py \
    --genotypes clade/genotypes.csv --tree clade/true_tree.nwk \
    --onset-ref-tree biopsy/true_tree.nwk \
    --clone-signatures $CLONE_SIGS --clone-rates $CLONE_RATES --onset-sizes $ONSETS \
    --signatures-file "$COSMIC" --seed 0 --out clade/mutation_types.csv
  $PY $S/filter_rare_sites.py --genotypes biopsy/genotypes.csv \
    --mutation-types clade/mutation_types.csv --tree biopsy/true_tree.nwk \
    --min-freq $MIN_FREQ --outdir cellcoal_out
  ensure_reference
  $PY $S/cellcoal_to_phylinsic.py \
    cellcoal_out/genotypes.csv cellcoal_out/mutation_types.csv phylinsic_genotypes.txt \
    --reference "$REFERENCE" --chrom chr22 --mapping-out site_signature_map.csv
}
ensure_reference () {
  [ -s "$REFERENCE" ] && return
  echo "downloading GRCh38 chr22 from UCSC -> $REFERENCE"
  mkdir -p "$(dirname "$REFERENCE")"
  curl -fL https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz | gunzip > "$REFERENCE.tmp"
  mv "$REFERENCE.tmp" "$REFERENCE"
}

# 3) SigDyn on every per-cell cap x replicate (finished runs are skipped)
run_cap () {   # $1 = run name, $2 = cap (0 = all), $3 = replicate
  local B=$OUT/$1 CAP=$2 REP=$3 TAG
  TAG=$([ "$CAP" = 0 ] && echo all || printf %03d "$CAP")
  local D=$B/runs/cap${TAG}_r$REP
  if [ -f "$D/results/signatures/exposures_active_signatures.txt" ]; then echo "done: $D"; return; fi
  mkdir -p "$D/config" "$D/resources/sigprofiler_genomes" "$D/results/genotypes"
  sed 's/use_pon: true/use_pon: false/' "$SIGDYN/config/config.yaml" > "$D/config/config.yaml"   # simulated data: no germline
  cp "$COSMIC" "$D/resources/"
  touch "$D/resources/sigprofiler_genomes/GRCh38.installed"
  $PY $S/downsample_cells.py "$B/sim/phylinsic_genotypes.txt" \
    "$D/results/genotypes/phylinsic_genotypes.txt" --cap "$CAP" --seed "$REP"
  local SM="$SNAKEMAKE -s $SIGDYN/snakefile --directory $D --use-conda --conda-frontend conda --conda-prefix $CONDA_PREFIX_DIR --cores $CORES"
  $SM stage1_unrooted || return 1
  find_base_py || return 1
  (cd "$D" && "$BASE_PY" "$S/root_on_suggested_tip.py")   # replaces the manual FigTree rooting
  $SM stage2_full
  find "$D/results/signatures" -name "Decomposed_MutationType_Probabilities.txt" -exec gzip -f {} \;
  echo "finished $D"
}
find_base_py () {   # SigDyn's conda env that has ete3 (exists once stage1 has run)
  [ -n "$BASE_PY" ] && return
  for p in "$CONDA_PREFIX_DIR"/*/bin/python; do
    if [ -x "$p" ] && "$p" -c "import ete3" 2>/dev/null; then BASE_PY=$p; return; fi
  done
  echo "no conda env with ete3 found in $CONDA_PREFIX_DIR (set BASE_PY)"; return 1
}
sigdyn_one () {   # $1 = run name
  local B=$OUT/$1
  [ -f "$B/sim/phylinsic_genotypes.txt" ] || { echo "run the label step first"; exit 1; }
  mkdir -p "$B/runs" "$B/logs"
  run_cap "$1" 0 1 > "$B/logs/capall_r1.log" 2>&1 && echo "ok $1 cap all" || echo "FAILED $1 cap all (see $B/logs)"
  for rep in $REPS; do
    for cap in $CAPS; do
      run_cap "$1" $cap $rep > "$B/logs/cap$(printf %03d $cap)_r$rep.log" 2>&1 \
        && echo "ok $1 cap$cap r$rep" || echo "FAILED $1 cap$cap r$rep (see $B/logs)"
    done
  done
}
sigdyn () { sigdyn_one "$NAME"; }

# 4) leaf-level scoring of every run (plain, and with SBS12 counted as SBS26)
evaluate () {
  cd "$OUT"
  $PY $S/evaluate_leaves.py --sim $NAME/sim --runs "$NAME/runs/cap*_r*" --out $NAME/results
  $PY $S/evaluate_leaves.py --sim $NAME/sim --runs "$NAME/runs/cap*_r*" --out $NAME/results_merged_SBS12 --merge SBS12:SBS26
}

# 5) figures -> $OUT/figures
plots () {
  cd "$OUT"
  $PY $S/plot_downsampling.py --out figures \
    --runs "$NAME|SBS26 = larger clone (369 cells)"
  for p in plot_simple plot_simple_nonsimulated plot_percell plot_percell_scatter plot_pearson; do
    $PY $S/$p.py
  done
}

# optional: search tumours (seeds x sizes) for a biopsy clade like the real data
scan () {
  CELLCOAL="$CELLCOAL" $PY $S/scan_tumours.py --seeds 1-20 --tumour-cells 1500 2000 3000 \
    --n-sites 11600 --mut-rate $MUT_RATE --outdir "$OUT/scan"
}

STEP=${1:-all}
case $STEP in
  all)      simulate; label; sigdyn; evaluate; plots ;;
  simulate) simulate ;;
  label)    label ;;
  sigdyn)   sigdyn ;;
  evaluate) evaluate ;;
  plots)    plots ;;
  scan)     scan ;;
  *) echo "usage: bash main.sh [all|simulate|label|sigdyn|evaluate|plots|scan]"; exit 1 ;;
esac
