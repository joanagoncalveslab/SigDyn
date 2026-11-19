#!/usr/bin/env python3
"""
Build genotype matrices with increasing missingness.

Same idea as make_masked_genotypes.ipynb:
  - start from the imputed, PON-filtered genotypes (every entry has H/R/A);
  - an entry is "truly missing" when it is empty in BOTH ref_count and alt_count;
  - real (covered) entries are NEVER blanked.

Each truly-missing entry gets one random number u; level X blanks the truly-missing
entries with u < q_X, where q_X makes X % of the whole matrix blank. The same u is used
for every level, so levels are nested (15 % ⊂ 45 % ⊂ 65 % ⊂ top). The top level blanks
ALL truly-missing entries and is named by its real percentage (e.g. p86, p90).

Each level is written where SigDyn expects its (PON-filtered) input:
  <out>/pXX/results/genotypes/genotypes_no_PON.txt
plus <out>/levels_summary.csv.

Usually called by run_levels.py; can be run alone:
  python make_levels.py --geno G.txt --ref ref_count.all.txt --alt alt_count.all.txt --out levels/neck
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_LEVELS = [0, 15, 45, 65]   # + the top level (all truly-missing) is always added


def make_levels(geno_file, ref_file, alt_file, out, levels=DEFAULT_LEVELS, cells=None, seed=0):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)

    geno = pd.read_csv(geno_file, sep="\t", index_col=0, dtype=str)
    if cells:   # random subsample of cells, for testing
        cols = rng.choice(geno.columns, size=min(cells, geno.shape[1]), replace=False)
        geno = geno[sorted(cols)]
    assert geno.notna().values.all(), "starting genotypes should have no empty entries"

    ref = pd.read_csv(ref_file, sep="\t", index_col=0)
    alt = pd.read_csv(alt_file, sep="\t", index_col=0)
    for name, m in [("ref", ref), ("alt", alt)]:
        assert geno.index.isin(m.index).all(), f"{(~geno.index.isin(m.index)).sum()} sites missing from {name} counts"
        assert geno.columns.isin(m.columns).all(), f"{(~geno.columns.isin(m.columns)).sum()} cells missing from {name} counts"
    ref = ref.reindex(index=geno.index, columns=geno.columns)
    alt = alt.reindex(index=geno.index, columns=geno.columns)
    n_disc = int((ref.isna() ^ alt.isna()).values.sum())
    missing = (ref.isna() & alt.isna()).values
    del ref, alt

    n_total, n_miss = missing.size, int(missing.sum())
    top = round(100 * n_miss / n_total)
    print(f"{geno.shape[0]} sites x {geno.shape[1]} cells | truly missing {100 * n_miss / n_total:.2f} % "
          f"| ref/alt discordant entries: {n_disc}")

    u = rng.random(missing.shape, dtype=np.float32)   # one random number per entry
    values = geno.values
    rows = []
    for lvl in sorted(set(levels)) + [None]:
        if lvl is None:
            lvl, q = top, 1.0
        elif lvl >= top:
            print(f"  skip p{lvl:02d}: not below the true missingness ({top} %)")
            continue
        else:
            q = lvl / 100 * n_total / n_miss
        blank = missing & (u < q)
        assert not (blank & ~missing).any()          # real entries untouched
        path = out / f"p{lvl:02d}" / "results" / "genotypes" / "genotypes_no_PON.txt"
        path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(np.where(blank, "", values), index=geno.index, columns=geno.columns).to_csv(path, sep="\t")
        frac = blank.sum() / n_total
        rows.append(dict(level=lvl, frac_missing=round(float(frac), 4), sites=geno.shape[0], cells=geno.shape[1]))
        print(f"  p{lvl:02d}: {100 * frac:.2f} % blank -> {path}")

    pd.DataFrame(rows).to_csv(out / "levels_summary.csv", index=False)
    return [r["level"] for r in rows]


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--geno", required=True, help="imputed genotypes, already without PON sites")
    ap.add_argument("--ref", required=True, help="ref_count.all.txt of the same PhylinSic run")
    ap.add_argument("--alt", required=True, help="alt_count.all.txt of the same PhylinSic run")
    ap.add_argument("--out", required=True)
    ap.add_argument("--levels", type=int, nargs="+", default=DEFAULT_LEVELS)
    ap.add_argument("--cells", type=int, default=None, help="random subsample of cells (testing)")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    make_levels(a.geno, a.ref, a.alt, a.out, a.levels, a.cells, a.seed)
