#!/usr/bin/env python3
"""
Run YOUR SigDyn snakefile on one missingness level at a time (snakefile and scripts untouched).

Each level gets its own SigDyn working directory  <out>/<name>/pXX/  (snakemake --directory),
with config/ and resources/ linked from the SigDyn folder and the level's genotypes placed at
results/genotypes/genotypes_no_PON.txt. The PON rules (filter_variants, phylinsic_external)
are excluded with --allowed-rules, because the input genotypes are already PON-filtered.

  levels:  build the genotype files of ALL levels of a dataset (one random draw -> nested levels)
  stage1:  one level: snakemake stage1_unrooted   (fasta, FastTree, root suggestion)
  -> root that level's tree with root_trees.py (replaces the manual FigTree step)
  stage2:  one level: snakemake stage2_full       (ancestral profiles, signatures)

Usage (from SigDyn-main/missingness/):
  python run_levels.py levels --name neck                 # inputs from inputs/neck/ (see README)
  python run_levels.py stage1 --name neck --level 45
  python root_trees.py --name neck --level 45 --cell <CELL>
  python run_levels.py stage2 --name neck --level 45

--fasttree BINARY (stage1): build the tree with that FastTree binary instead of the snakefile's conda
FastTree (e.g. the single-precision FastTree 2.1.11 of the paper). The snakefile still makes
mutations.fa and root_cell.txt; only the FastTree command (the same "FastTree -nt mutations.fa")
is run with the given binary. The version used is written to results/trees/FASTTREE_USED.

--reuse-leaves LEVEL_DIR (stage2): the per-leaf signature results do not depend on the tree, so copy
them (results/signatures/all_leaves, per_leaf, Custom_COSMIC_SBS_v3.4.txt) from LEVEL_DIR, a finished
run of the same level genotypes, and leave out the leaf rules (prepare_vcfs_leaves,
assign_signatures_leaves, make_custom_cosmic_signatures) with --allowed-rules. Only the tree-dependent
part (ancestral profiles, all-node VCFs, active-signature exposures) is run.
"""
import argparse
import re
import shutil
import subprocess
import time
from pathlib import Path

from make_levels import DEFAULT_LEVELS, make_levels

HERE = Path(__file__).resolve().parent
EXCLUDED_RULES = {"filter_variants", "phylinsic_external"}
# datasets whose results/vcf_* folders are deleted after a successful stage 2 (disk space:
# ~5.5 GB per neck level; the analyses do not use them)
DELETE_VCFS_FOR = {"neck"}
LEAF_RULES = {"prepare_vcfs_leaves", "assign_signatures_leaves", "make_custom_cosmic_signatures"}
LEAF_OUTPUTS = ["results/signatures/all_leaves", "results/signatures/per_leaf",
                "results/signatures/Custom_COSMIC_SBS_v3.4.txt"]


def snakemake_cmd(sigdyn, workdir, target, cores, exclude=()):
    snakefile = sigdyn / "snakefile"
    rules = re.findall(r"^rule (\w+):", snakefile.read_text(), flags=re.M)
    allowed = [r for r in rules if r not in EXCLUDED_RULES | set(exclude)]
    return ["snakemake", "-s", str(snakefile), "--directory", str(workdir),
            "--use-conda", "--conda-frontend", "conda",
            "--conda-prefix", str(sigdyn / ".snakemake/conda"),   # reuse SigDyn's conda envs
            "--cores", str(cores), "--rerun-triggers", "mtime",
            "--allowed-rules", *allowed, "-p", target]


def setup_workdir(sigdyn, workdir):
    for d in ["config", "resources"]:
        link = workdir / d
        if not link.exists():
            link.symlink_to(sigdyn / d)


def run(sigdyn, workdir, target, cores, exclude=()):
    log = workdir / f"snakemake_{Path(target).name}.log"
    print(f"[{time.strftime('%H:%M:%S')}] {workdir.name}: {target} (log: {log})", flush=True)
    with open(log, "w") as fh:
        rc = subprocess.run(snakemake_cmd(sigdyn, workdir, target, cores, exclude), stdout=fh, stderr=subprocess.STDOUT).returncode
    print(f"[{time.strftime('%H:%M:%S')}] {workdir.name}: {target} {'done' if rc == 0 else f'FAILED (exit {rc}), see log'}", flush=True)
    return rc


def run_fasttree(sigdyn, w, binary, cores):
    """snakefile up to mutations.fa, then `<binary> -nt mutations.fa > cell_tree.unrooted` (as rule fasttree)."""
    trees = w / "results/trees"
    tree = trees / "cell_tree.unrooted"
    version = subprocess.run([str(binary), "-expert"], capture_output=True, text=True).stderr.splitlines()[0]
    if tree.exists():
        used = (trees / "FASTTREE_USED").read_text().strip() if (trees / "FASTTREE_USED").exists() else "unknown"
        if used != version:
            raise SystemExit(f"{tree} exists, made with: {used}. Delete it to rebuild with {binary}.")
        print(f"{w.name}: {tree.name} already made with {version}")
        return 0
    rc = run(sigdyn, w, "results/trees/mutations.fa", cores)
    if rc:
        return rc
    print(f"[{time.strftime('%H:%M:%S')}] {w.name}: {binary} -nt mutations.fa  ({version})", flush=True)
    tmp = trees / "cell_tree.unrooted.tmp"
    with open(tmp, "w") as out, open(w / "fasttree.log", "w") as log:
        rc = subprocess.run([str(binary), "-nt", str(trees / "mutations.fa")], stdout=out, stderr=log).returncode
    if rc:
        print(f"{w.name}: FastTree FAILED (exit {rc}), see {w / 'fasttree.log'}")
        return rc
    tmp.rename(tree)                                  # newer than mutations.fa -> snakemake keeps it
    (trees / "FASTTREE_USED").write_text(f"{version}\n{binary}\n")
    print(f"[{time.strftime('%H:%M:%S')}] {w.name}: FastTree done", flush=True)
    return 0


def level_dir(out, level):
    d = out / f"p{int(str(level).lstrip('p')):02d}"
    if not (d / "results/genotypes/genotypes_no_PON.txt").exists():
        have = sorted(p.name for p in out.glob("p*")) or "none - run the 'levels' step first"
        raise SystemExit(f"no level {d.name} in {out} (available: {have})")
    return d


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step", choices=["levels", "stage1", "stage2"])
    ap.add_argument("--name", required=True, help="dataset name, e.g. breast or neck")
    ap.add_argument("--level", help="[stage1/stage2] the level to run, e.g. 45 or p45")
    ap.add_argument("--geno", help="[levels] imputed genotypes, already without PON sites "
                                   "(default: inputs/<name>/genotypes_no_PON.txt)")
    ap.add_argument("--ref", help="[levels] ref counts of the same PhylinSic run (default: inputs/<name>/ref_count.all.txt)")
    ap.add_argument("--alt", help="[levels] alt counts of the same PhylinSic run (default: inputs/<name>/alt_count.all.txt)")
    ap.add_argument("--levels", type=int, nargs="+", default=DEFAULT_LEVELS, help="[levels] besides the top level")
    ap.add_argument("--cells", type=int, default=None, help="[levels] random subsample of cells (testing)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--sigdyn", default=str(HERE.parent), help="SigDyn folder (default: parent of this folder)")
    ap.add_argument("--out", default=str(HERE / "levels"))
    ap.add_argument("--cores", type=int, default=4, help="cores for snakemake")
    ap.add_argument("--reuse-leaves", metavar="LEVEL_DIR",
                    help="[stage2] copy the tree-independent per-leaf signature results from this finished level")
    ap.add_argument("--fasttree", help="[stage1] FastTree binary to use instead of the snakefile's conda FastTree")
    a = ap.parse_args()

    out = Path(a.out).resolve() / a.name

    if a.step == "levels":
        inp = HERE / "inputs" / a.name
        a.geno = a.geno or inp / "genotypes_no_PON.txt"
        a.ref = a.ref or inp / "ref_count.all.txt"
        a.alt = a.alt or inp / "alt_count.all.txt"
        for f in (a.geno, a.ref, a.alt):
            if not Path(f).exists():
                ap.error(f"input not found: {f} (put it in {inp}/ or pass --geno/--ref/--alt)")
        if (out / "levels_summary.csv").exists():
            raise SystemExit(f"levels already exist in {out} (delete that folder to rebuild them)")
        made = make_levels(a.geno, a.ref, a.alt, out, a.levels, a.cells, a.seed)
        print(f"\nNext, one level at a time: python run_levels.py stage1 --name {a.name} --level <{'|'.join(map(str, made))}>")
        return

    if not a.level:
        ap.error(f"{a.step} needs --level (e.g. --level 45)")
    sigdyn = Path(a.sigdyn).resolve()
    assert (sigdyn / "snakefile").exists(), f"no snakefile in {sigdyn}"
    w = level_dir(out, a.level)
    setup_workdir(sigdyn, w)

    if a.step == "stage1":
        rc = 0
        if a.fasttree:
            rc = run_fasttree(sigdyn, w, Path(a.fasttree).resolve(), a.cores)
        if rc == 0:
            rc = run(sigdyn, w, "stage1_unrooted", a.cores)
        f = w / "results/trees/root_cell.txt"
        if f.exists():
            print(f"\nSuggested root cell (max RR): {f.read_text().strip()}")
            print(f"Next: python root_trees.py --name {a.name} --level {w.name} --cell <CELL>")
    else:
        if not (w / "results/trees/cell_tree.rooted").exists():
            raise SystemExit(f"{w.name} is not rooted yet - run root_trees.py first")
        if (w / "VCFS_DELETED").exists():   # snakemake would otherwise redo all of stage 2
            raise SystemExit(f"{w.name}: stage 2 already finished (its vcf folders were deleted). "
                             f"To redo it, delete {w}/results/signatures and VCFS_DELETED first.")
        exclude = ()
        if a.reuse_leaves:
            src = Path(a.reuse_leaves).resolve()
            g1, g2 = (src / "results/genotypes/genotypes_no_PON.txt").read_bytes(), (w / "results/genotypes/genotypes_no_PON.txt").read_bytes()
            if g1 != g2:
                raise SystemExit(f"--reuse-leaves: {src} has different genotypes than {w.name}")
            for rel in LEAF_OUTPUTS:
                if not (w / rel).exists():
                    (w / rel).parent.mkdir(parents=True, exist_ok=True)
                    subprocess.run(["cp", "-cpR", str(src / rel), str(w / rel)], check=True)   # -p keeps mtimes
            (w / "LEAVES_REUSED_FROM").write_text(f"{src}\n")
            print(f"{w.name}: per-leaf signatures copied from {src}")
            exclude = LEAF_RULES
        rc = run(sigdyn, w, "stage2_full", a.cores, exclude)
        if rc == 0 and a.name in DELETE_VCFS_FOR:
            vcfs = sorted((w / "results").glob("vcf_*"))
            for d in vcfs:
                shutil.rmtree(d)
            (w / "VCFS_DELETED").write_text("\n".join(d.name for d in vcfs) + "\n")
            print(f"{w.name}: deleted {', '.join(d.name for d in vcfs)} (DELETE_VCFS_FOR={sorted(DELETE_VCFS_FOR)})")
    raise SystemExit(rc)


if __name__ == "__main__":
    main()
