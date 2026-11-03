# TRACE: Transferable Rejuvenation Assessment across Contexts and Evidence

Code, pre-specified protocol, results and figures for the paper *"Rejuvenation rankings from CRISPRi screens change with the
aging cohort and track cell growth programs"*. TRACE is a set of tests, fixed in part before analysis, for rankings that score gene knockdowns by how well they reverse an aging
expression signature: it swaps the aging cohort, compares score definitions with random signatures and a proliferation-only
baseline, compares lines and replicate screens, and tests rankings against held-out labels.

## Install

Python 3.13.7 on Linux, macOS or Windows.

```bash
pip install -r requirements.txt
```

or `conda env create -f environment.yml`.

## Run

```bash
python run_all.py --tests     # unit tests on synthetic data (seconds; no downloads)
python run_all.py             # everything: download and check inputs, preprocess, confirmatory H1-H5, exploratory analyses, report, figures
python run_all.py --from x19  # resume at a step
```

The summary is written to `results/trace_report.csv`: one row per check (cohort transfer with its power check, programs inside the
signature, random signatures, replicates, labels against a proliferation-only baseline, knockdown quality), with its value, the reference
it is judged against and the result (pass, fail, flag or info), following Box 1 of the paper.

Inputs (about 16 GB) are listed with URL, version, access date and SHA-256 in `DATA_MANIFEST.csv`; `src/download.py` fetches and
checks them. Resample counts and seeds are in `src/config.py` with the paper's values as defaults. Set `TRACE_ROOT` to run from
another folder. Each step prints its wall time; the analyses were run on one 8-core, 8 GB workstation.

To check another aging cohort, give `src/trace_analysis.py` two CSV files (genes × samples expression and a sample table with
`age` and `sex`): `TRACE_DISCOVERY=table TRACE_EXPR=expr.csv TRACE_META=meta.csv TRACE_RESULTS=results/my_run python src/trace_analysis.py`.

## Layout

| Path | Content |
|---|---|
| `RESEARCH_PROTOCOL.md` | Protocol v1.0, frozen at tag `protocol-v1.0` (commit `1eb8681`) before any knockdown was scored |
| `PROTOCOL_FREEZE_RECORD.json` | Hashes and commit time of the frozen protocol |
| `PROTOCOL_AMENDMENTS.md` | Deviations (D1–D11) and dated exploratory analyses (X01–X33) |
| `DATA_MANIFEST.csv` | Every input file with URL, version, access date, size and SHA-256 |
| `src/config.py` | Paths, resample counts and seeds |
| `src/download.py` | Cross-platform download and checksum check |
| `src/aging.py`, `src/fleischer_meta.py` | Aging cohorts, age models and signatures |
| `src/perturb.py`, `src/pseudobulk_sc.py`, `src/orthologs.py` | Perturb-seq pseudobulk shifts (HepG2 and Jurkat pseudobulked here) and orthologs |
| `src/trace_analysis.py` | Confirmatory tests H1–H5 → `results/results.json` |
| `src/rev_*.py`, `src/rev3_*.py`, `src/x_*.py` | Exploratory analyses, in the order of the amendment log |
| `src/make_figures.py`, `src/make_figures_supp.py` | Figures 1–5, S1–S10 and the graphical abstract, with per-panel source data |
| `src/figstyle.py`, `src/schematics.py`, `src/icons.py` | Figure style, schematic panels, licensed icons (`figures/icons/CREDITS.md`) |
| `results/` | Machine-readable results; data dictionary in `results/README.md` |
| `figures/` | Final figures (PDF and PNG) and `figures/source_data/` |
| `tests/test_core.py` | Unit tests |

## Figures and tables → scripts

| Item | Script (function) | Main inputs |
|---|---|---|
| Figure 1A, 2A, graphical abstract | `src/schematics.py` via `make_figures.py` | – |
| Figure 1B–1H | `make_figures.py` (`fig1`) | `rev_cohort_table.csv`, `rev3_point_tstats.pkl`, `rev_cohorts.json`, `x_cohort_*`, `rev3_cohort_sensitivity*.csv`, `signature_hallmark_enrichment.csv`, `rev3_cohort_age_prediction.csv` |
| Figure 2B–2F | `make_figures.py` (`fig2`) | `scores_*.csv`, `rev3_scores_prolif_dependence.csv`, `rev_scores.json`, `rev3_scores_random_signatures.csv` |
| Figure 3A–3E | `make_figures.py` (`fig3`) | `results.json`, `rev_scores.json`, `rev3_reliability_*.csv`, `rev_h5_boot.csv`, `rev3_h5_multiline.csv` |
| Figure 4A–4E | `make_figures.py` (`fig4`) | `rev_network*.{csv,json}`, `rev3_network_*.csv` |
| Figure 5A–5D | `make_figures.py` (`fig5`) | `rev3_scores_evidence.csv`, `rev3_scores_paired_auroc.csv`, `rev3_depmap_label_fisher.csv`, `rev3_scores_label_essential_overlap.csv` |
| Figures S1–S10 | `make_figures_supp.py` (`S1` … `S10`) | see the function for each figure |
| Table 1 | – | `results.json` |
| Supplementary tables | – | see `results/README.md`; `x_exploratory_tests.csv` lists every exploratory test |

## License and citation

Code: MIT (`LICENSE`). Results, figures and source data: CC BY 4.0. Third-party input data keep their own terms (OncoKB, GTEx
and the Perturb-seq releases are not redistributed here). Icons: Servier Medical Art (CC BY 3.0) and Bioicons (CC0, CC BY 4.0).
Please cite the paper (see `CITATION.cff`). Code archive: Zenodo, https://doi.org/10.5281/zenodo.23763631 (all versions).
