"""Reproduce every result and figure of TRACE from the inputs listed in DATA_MANIFEST.csv.

    python run_all.py              # full pipeline: download and check, preprocess, confirmatory, exploratory, figures
    python run_all.py --tests      # unit tests only (seconds; no downloads needed)
    python run_all.py --from X     # start at the step whose name starts with X (for example --from x19)

Resample counts and seeds are in src/config.py (the paper's values are the defaults). Each step prints its wall time.
The manuscript text itself is not built here; every number it reports is in results/ (see results/README.md)."""
import os, subprocess, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
E18 = {"TRACE_DISCOVERY": "gtex_fibroblasts", "TRACE_RESULTS": "results/e18_gtex_discovery"}
STEPS = [
    # name, script and arguments, extra environment
    ("tests", ["tests/test_core.py"], {}),
    ("download", ["src/download.py"], {}),
    ("orthologs", ["src/orthologs.py"], {}),
    ("perturb-cache", ["src/perturb.py", "K562", "K562gw", "RPE1"], {}),
    ("pseudobulk-hepg2", ["src/pseudobulk_sc.py", "hepg2"], {}),
    ("pseudobulk-jurkat", ["src/pseudobulk_sc.py", "jurkat"], {}),
    # confirmatory: H1-H5 as pre-specified (protocol v1.0)
    ("confirmatory", ["src/trace_analysis.py"], {}),
    # exploratory, in the order of the amendment log (PROTOCOL_AMENDMENTS.md)
    ("x01-cohorts-split-halves", ["src/rev_cohorts.py"], {}),
    ("x02-score-definitions", ["src/rev_scores.py"], {}),
    ("x03-network-evidence", ["src/rev_network_evidence.py"], {}),
    ("x06-h5-bootstraps", ["src/rev_h5_boot.py"], {}),
    ("x07-network-500-nulls", ["src/rev_network_null.py"], {}),
    ("x08-x10-cohort-sensitivity", ["src/rev3_cohort.py"], {}),
    ("x10-bias-corrected-intervals", ["src/rev3_cohort_bc.py"], {}),
    ("x11-x16-scores", ["src/rev3_scores.py"], {}),
    ("x13-depmap", ["src/rev3_depmap_overlap.py"], {}),
    ("x12-reliability", ["src/rev3_reliability.py"], {}),
    ("x15-x22-network-in-fold", ["src/rev3_network.py"], {}),
    ("x14-h5-multi-line", ["src/rev3_h5_multiline.py"], {}),
    ("x18-gtex-discovery", ["src/trace_analysis.py"], E18),
    ("x19-x21-cohort-nulls", ["src/x_cohort_nulls.py"], {}),
    ("x23-gtex-signature", ["src/x_gtex_signature.py"], {}),
    ("x25-cellage-method", ["src/x_cellage_method.py"], {}),
    ("x25-v6-v7", ["src/x_v6v7_effects.py"], {}),
    ("x25-exploratory-tests", ["src/x_exploratory_tests.py"], {}),
    # figures (final numbering) and their source data
    ("figures-main", ["src/make_figures.py"], {}),
    ("figures-supplement", ["src/make_figures_supp.py"], {}),
]

if __name__ == "__main__":
    steps = STEPS[:1] if "--tests" in sys.argv else STEPS
    if "--from" in sys.argv:
        start = sys.argv[sys.argv.index("--from") + 1]; i = next(i for i, s in enumerate(steps) if s[0].startswith(start)); steps = steps[i:]
    elif "--tests" not in sys.argv:
        for p in (ROOT / "data" / "proc").glob("ctx_*.pkl"): p.unlink()  # never reuse a stale perturbation cache in a full run
    t_all = time.time()
    for name, cmd, env in steps:
        t0 = time.time(); print(f">> {name}: {' '.join(cmd)}", flush=True)
        subprocess.run([sys.executable, str(ROOT / cmd[0])] + cmd[1:], check=True, cwd=ROOT, env={**os.environ, **env})
        print(f"   {name} done in {(time.time() - t0) / 60:.1f} min", flush=True)
    print(f"all {len(steps)} steps completed in {(time.time() - t_all) / 60:.1f} min")
