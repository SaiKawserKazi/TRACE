"""Paths, resample counts and seeds for TRACE, in one place.

ROOT is the repository root (or TRACE_ROOT if set). Each count has its own name and the value used in the paper as default;
an environment variable with the same name overrides it (for quick test runs). Every analysis writes the counts it used
into its output JSON."""
import os
from pathlib import Path

ROOT = Path(os.environ.get("TRACE_ROOT", Path(__file__).resolve().parents[1]))
RAW, PROC, RES, FIG = ROOT / "data" / "raw", ROOT / "data" / "proc", ROOT / "results", ROOT / "figures"
for _p in (PROC, RES, FIG):
    _p.mkdir(parents=True, exist_ok=True)


def _int(name, default):
    return int(os.environ.get(name, default))


# confirmatory (protocol v1.0)
N_H5_BOOT_PRESPECIFIED = _int("N_H5_BOOT_PRESPECIFIED", 200)      # H5a donor bootstraps
N_H3_NULLS_PRESPECIFIED = 10                                       # H3a degree-preserving randomizations (seeds 0-9)
N_H1B_BOOT = _int("N_H1B_BOOT", 2000)
N_LABEL_PERM = _int("N_LABEL_PERM", 10000)
CV_FOLDS, CV_SEED = 5, 0

# exploratory
N_COHORT_BOOT = _int("N_COHORT_BOOT", 2000)       # donor bootstraps per cohort specification
N_SPLIT_HALVES = _int("N_SPLIT_HALVES", 500)
N_COHORT_PERM = _int("N_COHORT_PERM", 1000)       # permuted-age null (pre-specified model)
N_SPEC_PERM = _int("N_SPEC_PERM", 1000)           # permuted-age null per discovery specification; permuted-sex null
N_POWER_SUBSETS = _int("N_POWER_SUBSETS", 1000)
N_PLANTED = _int("N_PLANTED", 200)
N_H5_BOOT_EXTENDED = _int("N_H5_BOOT_EXTENDED", 2000)
N_NETWORK_NULLS = _int("N_NETWORK_NULLS", 500)    # seeds 1000-1499
N_NETWORK_SEEDS = _int("N_NETWORK_SEEDS", 10)
N_BLOCK_BOOT = _int("N_BLOCK_BOOT", 2000)
N_RANDOM_SIGNATURES = _int("N_RANDOM_SIGNATURES", 100)
N_RELIABILITY_SPLITS = _int("N_RELIABILITY_SPLITS", 20)
N_CI_BOOT = _int("N_CI_BOOT", 2000)
N_PROC = _int("N_PROC", 4)
