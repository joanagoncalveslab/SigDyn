# =========================
# BUILD CUSTOM COSMIC SIGNATURE MATRIX
# =========================

import sys
log = open(snakemake.log[0], "w")
sys.stdout = log
sys.stderr = log

#!/usr/bin/env python3
import os
import glob
import pandas as pd

EXPOSURES_DIR = snakemake.input["exposures_dir"]
OUT_FILE = snakemake.output[0]
MIN_FRACTION = float(snakemake.params["min_fraction"])

# Prefer explicit input in Snakefile; fallback to default path
COSMIC_FILE = snakemake.input.get("cosmic", "resources/COSMIC_SBS_v3.4.txt")

# ---- expected location in your setup ----
expected = os.path.join(
    EXPOSURES_DIR,
    "Assignment_Solution",
    "Activities",
    "Assignment_Solution_Activities.txt"
)

# Some versions omit ".txt" in the filename you see in Finder vs actual file.
# We'll accept both.
expected_noext = os.path.join(
    EXPOSURES_DIR,
    "Assignment_Solution",
    "Activities",
    "Assignment_Solution_Activities"
)

# Pick the correct exposures file
if os.path.exists(expected):
    EXPOSURE_FILE = expected
elif os.path.exists(expected_noext):
    EXPOSURE_FILE = expected_noext
else:
    # fallback: find any Activities file under exposures dir
    candidates = glob.glob(os.path.join(EXPOSURES_DIR, "**", "*Activities*.txt"), recursive=True)
    if not candidates:
        candidates = glob.glob(os.path.join(EXPOSURES_DIR, "**", "*Activities*"), recursive=True)

    if not candidates:
        raise FileNotFoundError(
            f"Could not find an Activities file under '{EXPOSURES_DIR}'.\n"
            f"Tried:\n  {expected}\n  {expected_noext}\n"
        )

    EXPOSURE_FILE = sorted(candidates, key=len)[0]

print("Using exposure file:", EXPOSURE_FILE)

# Load exposures + COSMIC
exposures = pd.read_csv(EXPOSURE_FILE, sep="\t", index_col=0)
cosmic = pd.read_csv(COSMIC_FILE, sep="\t", index_col=0)

# Orientation fix: ensure rows are signatures, columns are samples
idx_overlap = len(set(exposures.index).intersection(set(cosmic.columns)))
col_overlap = len(set(exposures.columns).intersection(set(cosmic.columns)))
if col_overlap > idx_overlap:
    exposures = exposures.T
    print("Transposed exposures so rows are signatures.")

# Active signatures: exposure > 0 in >= MIN_FRACTION of samples
active_fraction = (exposures > 0).sum(axis=1) / exposures.shape[1]
active_signatures = active_fraction[active_fraction >= MIN_FRACTION].index.tolist()

if not active_signatures:
    raise ValueError(
        f"No active signatures found with min_fraction={MIN_FRACTION}.\n"
        f"Exposure file used: {EXPOSURE_FILE}"
    )

print(f"Active signatures ({len(active_signatures)}): {active_signatures}")

# Subset COSMIC
missing = [s for s in active_signatures if s not in cosmic.columns]
if missing:
    raise KeyError(
        "Active signatures not found in COSMIC reference:\n"
        + "\n".join(missing)
        + f"\nCOSMIC file: {COSMIC_FILE}"
    )

custom_cosmic = cosmic[active_signatures]

os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)
custom_cosmic.to_csv(OUT_FILE, sep="\t")
print("Wrote:", OUT_FILE)
