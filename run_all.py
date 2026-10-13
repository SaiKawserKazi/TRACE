"""Reproduce every result, figure and the manuscript from the raw inputs listed in DATA_MANIFEST.csv.

    python run_all.py            # full pipeline (about 90 min on 8 CPU cores / 8 GB RAM, after downloads)
    python run_all.py --tests    # unit and smoke tests only (seconds; no downloads needed)

Steps (all deterministic; seeds fixed in code):
  1  tests/test_core.py                    unit tests
  2  src/orthologs.py                      human orthologs (Alliance)                  -> data/raw/ortho_*.tsv
  3  src/pseudobulk_sc.py hepg2|jurkat     HepG2 / Jurkat pseudobulk                   -> data/raw/*_normalized_bulk.h5ad
  4  src/trace_analysis.py                 pre-registered analyses H1-H5               -> results/results.json
  5  src/rev_cohorts.py                    review round: cohort sensitivity            -> results/rev_cohorts.json
  6  src/rev_scores.py                     review round: score-definition sensitivity  -> results/rev_scores.json
  7  src/rev_network_evidence.py           review round: network CIs, evidence audit   -> results/rev_network_evidence.json
  8  src/figures_v1.py                     graphical abstract, Figs 1-5, S1-S5        -> figures/
  9  src/build_v1.py                       manuscript V1.0, Cell Press format (.docx/.pdf) -> manuscript/
"""
import subprocess, sys, os

ROOT = os.path.dirname(os.path.abspath(__file__))
PY = [sys.executable, "-I"]
STEPS = [["tests/test_core.py"], ["src/orthologs.py"], ["src/pseudobulk_sc.py", "hepg2"], ["src/pseudobulk_sc.py", "jurkat"],
         ["src/trace_analysis.py"], ["src/rev_cohorts.py"], ["src/rev_scores.py"], ["src/rev_network_evidence.py"],
         ["src/figures_v1.py"], ["src/build_v1.py"]]

if __name__ == "__main__":
    steps = STEPS[:1] if "--tests" in sys.argv else STEPS
    for s in steps:
        print(">>", " ".join(s), flush=True)
        subprocess.run(PY + [os.path.join(ROOT, s[0])] + s[1:], check=True, cwd=ROOT)
    print("all steps completed")

