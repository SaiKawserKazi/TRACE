# TRACE — Transferable Rejuvenation Assessment across Contexts and Evidence

Code, frozen protocol and outputs for the manuscript *"Transcriptomic rejuvenation scores from perturbation
atlases are dominated by proliferation and do not replicate across aging cohorts: a pre-registered,
network-constrained benchmark"*.

## Layout

| Path | Content |
|---|---|
| `RESEARCH_PROTOCOL.md` | Pre-registered protocol v1.0 (git tag `protocol-v1.0`, committed before any scoring) |
| `src/download.ps1`, `src/pardl.py` | Data download (all public) |
| `src/orthologs.py` | Human orthologs from the Alliance of Genome Resources |
| `src/aging.py` | Aim 1: donor-aware aging signatures (Fleischer, GTEx) |
| `src/perturb.py`, `src/pseudobulk_sc.py` | Perturb-seq pseudobulk shifts (Replogle release; HepG2/Jurkat pseudobulked here) |
| `src/trace_analysis.py` | All hypothesis tests (H1–H5), network model, prioritisation → `results/results.json` |
| `src/figures.py` | Figures 1–5, S1–S2 and per-panel source data (`figures/source_data/`) |
| `src/build_manuscript.py` | Builds `manuscript/TRACE_manuscript.pdf` from `results/results.json` |

## Reproduce

```
powershell -File src/download.ps1          # small/medium files
python src/pardl.py <url> <out> <bytes>     # HepG2 / Jurkat single-cell files (GEO GSE264667)
python src/orthologs.py
python src/pseudobulk_sc.py hepg2
python src/pseudobulk_sc.py jurkat
python src/trace_analysis.py                # TRACE_NBOOT=200 by default
python src/figures.py
python src/build_manuscript.py
```

Requirements: Python 3.13, numpy, scipy, pandas, statsmodels, networkx, anndata, h5py, matplotlib, pymupdf;
Microsoft Edge (headless) for HTML→PDF. 8 GB RAM is sufficient.
