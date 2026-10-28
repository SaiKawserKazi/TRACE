# Results: data dictionary

Every file below is written by the script named in the second column (run order in `run_all.py`). Correlations are Spearman ρ
unless stated; AUROC = area under the ROC curve; screens are `K562` (essential-gene), `K562gw` (genome-wide), `RPE1`, `HepG2`, `Jurkat`.
Confirmatory results are in `results.json` only; everything else is exploratory (see `PROTOCOL_AMENDMENTS.md`, IDs in the third column).

| File | Written by | Log ID | Content |
|---|---|---|---|
| `results.json` | trace_analysis.py | H1–H5 | All pre-specified tests (statistics, P, Holm P, decisions), aging summary, per-screen H4 AUROCs |
| `age_t_fleischer.csv`, `age_t_gtex_fibroblast.csv`, `age_t_gtex_skin.csv` | trace_analysis.py | H5 | Per-gene age t-statistic, P and q (pre-specified model) per cohort |
| `signature_hallmark_enrichment.csv` | trace_analysis.py | – | Hallmark hypergeometric enrichment of the up and down signature arms |
| `scores_<screen>.csv` | trace_analysis.py | H1–H4 | Per-knockdown raw score R, R_adj (V1), program shifts π, σ, ι, μ |
| `robust_scores.csv` | trace_analysis.py | H4, H5 | Per-gene R_rob, number of lines, sign consistency |
| `h3_<screen>.csv` | trace_analysis.py | H3 | Held-out predictions of every H3 model per knockdown |
| `h5_boot_rho.csv` | trace_analysis.py | H5a | ρ of each of the 200 pre-specified donor bootstraps |
| `prioritisation.csv` | trace_analysis.py | – | Pareto objectives and front membership (genes in ≥ 2 lines) |
| `e18_gtex_discovery/` | trace_analysis.py with `TRACE_DISCOVERY=gtex_fibroblasts` | X18 | The same files for the GTEx-derived signature (here `age_t_fleischer.csv` holds the GTEx discovery t-statistics) |
| `rev_cohorts.json`, `rev_cohort_table.csv` | rev_cohorts.py | X01 | Split-half references (500 halvings) and the cohort description table |
| `rev_cohort_sensitivity.csv`, `rev_cohort_bootstrap.csv.gz` | rev_cohorts.py | X01, X05 | Superseded by `rev3_cohort_sensitivity.csv` (skin G0 contained ischemic time, deviation D8); not used by figures or text |
| `rev_scores.json`, `rev_scores_*.csv` | rev_scores.py | X02 | Score definitions V0–V5: cross-line ρ, ρ and top-5% overlap between definitions, candidates |
| `rev_network.csv`, `rev_network_evidence.json`, `rev_evidence_*.csv` | rev_network_evidence.py | X03, X04 | Pre-specified H3 with CIs; label coverage and sensitivity subsets |
| `rev_h5_boot.{json,csv}`, `rev_h5_top5_support.csv` | rev_h5_boot.py | X06 | 2,000 donor bootstraps of R_rob |
| `rev_network_null.json`, `rev_network_null_dist.csv` | rev_network_null.py | X07 | 500 degree-preserving randomized networks |
| `rev3_cohort.json`, `rev3_cohort_*.csv` | rev3_cohort.py | X08–X10, X17 | Cohort sensitivity (F0–F9 × G0–G3) with 2,000 donor bootstraps; permutations; power-matched subsets; cell-cycle and pre-named genes; age prediction |
| `rev3_cohort_sensitivity_bc.csv` | rev3_cohort_bc.py | X10 | Bias-corrected intervals for the cohort sensitivity table |
| `rev3_point_tstats.pkl`, `rev3_point_signatures.pkl` | rev3_cohort.py | X10 | Point-estimate t-statistics and signatures per specification (pandas pickles) |
| `rev3_scores.json`, `rev3_scores_*.csv`, `rev3_E5_inputs.pkl` | rev3_scores.py | X11, X13, X14, X16 | Random signatures; evidence AUROCs; paired DeLong tests; logistic models; proliferation dependence; Pareto filter; senescence-gene ranks; up-minus-down score |
| `rev3_depmap.json`, `rev3_depmap_label_fisher.csv` | rev3_depmap_overlap.py | X13 | DepMap 24Q2 common essentials and Fisher tests; also adds DepMap columns to `rev3_scores_label_essential_overlap.csv` |
| `rev3_reliability*.{json,csv}` | rev3_reliability.py | X12 | Split-cell and test–retest reliability, disattenuated ρ, size-matched non-targeting groups, held-out proliferation |
| `rev3_network*.{json,csv}` | rev3_network.py | X15, X22 | In-fold network analysis (10 seeds), baselines, block bootstrap, λ choice |
| `rev3_h5_multiline.{json,csv}` | rev3_h5_multiline.py | X14 | H5 bootstraps on all genes and on genes in ≥ 2 lines |
| `x_cohort_*.{json,csv}` | x_cohort_nulls.py | X19–X21 | Permuted-age null per specification; permuted-sex null; planted effect; hidden factors |
| `x_gtex_signature*.{json,csv}` | x_gtex_signature.py | X23 | Both signatures side by side: H1–H5, programs, evidence, paired tests, composition |
| `x_cellage_method.csv` | x_cellage_method.py | X25 | CellAge AUROC by manipulation method |
| `x_v6v7_effects.json` | x_v6v7_effects.py | X25 | V6 and V7 effects on top-5% lists and the Pareto front |
| `x_exploratory_tests.csv` | x_exploratory_tests.py | X25 | Every exploratory test with P value and family size |

Figure source data (one CSV per panel, final numbering) are in `figures/source_data/`.
