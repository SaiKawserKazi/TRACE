# Deviations and exploratory analysis log

Protocol: `RESEARCH_PROTOCOL.md` v1.0, frozen at commit `1eb8681` (tag `protocol-v1.0`, 2026-10-12 16:08 +0600, self-attested
local time). SHA-256 `01a573574e8225e5ba68ddd2757b81dcd765083274d9b3290b8a1bb0ca879fa7`. Unchanged since.

Confirmatory = H1–H5 exactly as written in the protocol. Everything else is exploratory.
Dates below come from the authors' records; the repository history shows when each change was committed.

## A. Deviations from the protocol

| ID | Date | Before or after the related result | Change | Reason | Possible influence |
|---|---|---|---|---|---|
| D1 | 2026-10-12 | Before any H4 result | Orthologs from the Alliance of Genome Resources instead of Ensembl BioMart | The Ensembl BioMart service was unavailable on 12 October 2026 | Changes the mapped label set; not tuned to results |
| D2 | 2026-10-12 | Before results | The K562 essential-gene and genome-wide screens are averaged into one line for R_rob, sign consistency and the Pareto ranking | Protocol silent | Defines R_rob |
| D3 | 2026-10-12 | Before results | H1 uses the K562 essential-gene screen as "K562"; H1b uses the 1,784 knockdowns measured in all screens rather than the 2,053 shared by the two K562 screens | Protocol silent on the first point; the second keeps one knockdown set for within- and between-line agreement | None expected |
| D4 | 2026-10-12 | Before those screens' results | HepG2 and Jurkat pseudobulked by us; their non-targeting controls are 100 random groups of non-targeting cells | No processed release with the same format | Normalization may differ slightly from the published release |
| D5 | 2026-10-12 | Before results | Ridge penalty λ = 1 for the H3 models | Protocol said "ridge" without λ | Untuned; may understate the combined model |
| D6 | 2026-10-12 | After the first H4 results | A file-size check (files under 1,000 bytes were skipped) silently dropped the yeast ortholog table in the first H4b run; in the corrected run the table passed the check and H4b was rerun. The check now raises an error instead of skipping | Bug | H4b AUROC 0.52 (Holm P = 0.18) before, 0.51 (Holm P = 0.27) after; conclusion unchanged |
| D7 | 2026-10-12 | After the pre-specified results | Early exploratory additions: early sensitivity run on ages 18 to 80 (replaced by F1, ages 20 to 79), adjusted within-K562 agreement, calibrated line-transfer model, single-screen evidence tests, on-target annotation | – | Labeled exploratory |
| D8 | 2026-10-26 | After Amendment 2 results | In Amendments 1 and 2 the GTEx-skin G0 model also contained ischemic time by mistake. G0 is now age + sex for both GTEx tissues; ischemic time enters the skin models from G1 | Bug | Skin G0 values changed (0.21 to 0.15; fibroblasts versus skin 0.02 to 0.06). No pre-specified result changed |
| D9 | 2026-10-26 | – | Exclusion flow written out: 143 GEO samples = 133 healthy donors + 10 progeria samples; 22 donors under 18 excluded, leaving 111 donors aged 19–96 years (GEO age field) | Clarity | None |
| D10 | 2026-10-26 | – | Loader fix: the fibroblast loader read sample columns [8:-1] and dropped the last column (one donor aged 17). All sample columns are now read | Bug | None: the donor is under 18 and excluded by the pre-specified age filter |
| D11 | 2026-10-28 | – | Identity genes: log-normalized gene mean of a line minus the mean of the other lines, among genes above that line's median expression, top 100 (protocol: z-score of the line's mean relative to the other lines on control expression). K562 and RPE1 use the gene means distributed with the Replogle release; HepG2 and Jurkat use non-targeting cells only | The code has used this rule since the first analysis; recorded here | Defines the identity program ι used by V1, V3 and V7 |

## B. Exploratory analyses (added after the pre-specified results were known)

### Amendment 1 (13 October 2026)

| ID | Analysis | Outputs |
|---|---|---|
| X01 | GEO metadata (platform, repository, biopsy site); cohort sensitivity F0–F7 × G0–G2 with donor-bootstrap CIs and split-half references | results/rev_cohorts.json, rev_cohort_sensitivity.csv |
| X02 | Score definitions V2–V5 | results/rev_scores.json |
| X03 | Network baselines and CIs; calibrated line transfer | results/rev_network.csv, rev_network_evidence.json |
| X04 | Evidence coverage, label permutations, label subsets and proliferation-only baseline | results/rev_evidence_*.csv |

### Amendment 2 (13 October 2026)

| ID | Analysis | Outputs |
|---|---|---|
| X05 | 2,000 donor bootstraps for every cohort specification (500 split halves) | results/rev_cohort_sensitivity.csv |
| X06 | 2,000 donor bootstraps for H5 (the 200 bootstraps in the protocol give the main result) | results/rev_h5_boot.json |
| X07 | 500 degree-preserving randomized networks for H3a (the 10 nulls in the protocol give the main result); 2,000 bootstrap resamples for all CIs; 10,000 label permutations | results/rev_network_null.json |

### Amendment 3 (26 to 27 October 2026)

| ID | Analysis | Outputs |
|---|---|---|
| X08 | Sex as a positive control (all genes, autosomes, X/Y, top-200 genes); sign check of pre-named age genes and canonical sex genes | results/rev3_cohort_sensitivity.csv, rev3_cohort_known_genes.csv |
| X09 | Permuted-age null (1,000 permutations, main model); power-matched 74-donor subsets; age per repository; variance inflation factor of age | results/rev3_cohort.json, rev3_cohort_permutations.csv |
| X10 | F8 (without AG donors aged 80 or older), F9 and G3 (five hidden factors); 2,000 donor bootstraps for all specifications; bias-corrected intervals (27 October) | results/rev3_cohort_sensitivity*.csv |
| X11 | 100 random expression-matched signatures through H1, H2 and H4 | results/rev3_scores_random_signatures.csv |
| X12 | Split-cell reliability (HepG2, Jurkat; 20 splits); test–retest reliability (two K562 screens); size-matched non-targeting groups; disattenuated agreement; proliferation measured in held-out cells | results/rev3_reliability*.csv |
| X13 | Paired DeLong tests, incremental logistic regression, minimum detectable AUROC; GenAge mouse-only and without yeast; label overlap with Hart CEGv2 and (27 October) DepMap 24Q2 common essentials | results/rev3_scores_*.csv, rev3_depmap*.{json,csv} |
| X14 | R_rob for single-line genes (equals R_adj); H4 and H5 on genes measured in at least two lines | results/rev3_h5_multiline.*, rev3_scores_evidence.csv |
| X15 | In-fold network analysis: in-fold residualization, inner-CV ridge penalty, 10 CV seeds, common subset per screen; nearest-neighbor, gene-covariate, physical-network, calibrated-network and proliferation baselines; block bootstrap over response clusters | results/rev3_network*.csv |
| X16 | Proliferation deciles, distance correlation, mutual information; definitions V6 (spline) and V7 (no magnitude); filtered Pareto front; ranks of senescence genes; up-minus-down score | results/rev3_scores_*.csv |
| X17 | Cross-cohort elastic-net age prediction | results/rev3_cohort_age_prediction.csv |
| X18 | TRACE run end to end with a GTEx-derived discovery signature | results/e18_gtex_discovery/ |

### Amendment 4 (28 October 2026)

| ID | Analysis | Outputs |
|---|---|---|
| X19 | Permuted-age null for every discovery specification (F0–F9) against GTEx G0, both tissues; empirical two-sided P for genome-wide and top-200 agreement | results/x_cohort_perm_p.csv, x_cohort_perm_by_spec.csv |
| X20 | Permuted-sex null for the sex positive control; planted shared age-like effect (500 genes, two effect sizes, 200 replicates each) on age-permuted data | results/x_cohort_nulls.json, x_cohort_sex_perm.csv, x_cohort_planted.csv |
| X21 | Discovery hidden factors: association with age, repository and platform; Hallmark enrichment of the factor loadings | results/x_cohort_hidden_factors.csv, x_cohort_hidden_factor_hallmark.csv |
| X22 | Network baselines rerun with a wider ridge grid (10^-6 to 10^3) and a nearest-neighbor predictor that uses the gene's response in the other lines only; the own-profile predictor kept as an upper reference | results/rev3_network*.csv |
| X23 | GTEx-derived signature (X18): proliferation-only baseline, paired AUROC tests, oxidative-phosphorylation and stress program tests, Hallmark composition | results/x_gtex_signature*.{json,csv} |
| X24 | Random-signature null compared with the real signature without proliferation genes (V2 input) | results/rev3_scores.json, figures/source_data/Fig2E.csv |
| X25 | CellAge labels by manipulation method (CellAge web table); V6 and V7 effects on the top-5% lists and the Pareto front; table of all exploratory tests with family sizes | results/x_cellage_method.csv, x_v6v7_effects.json, x_exploratory_tests.csv |

### Amendment 5 (29 October 2026)

| ID | Analysis | Outputs |
|---|---|---|
| X26 | Random-signature null: every raw-score statistic compared with the real signature without proliferation genes (V2 input), with one-sided empirical P | results/x_random_signature_v2.csv |
| X27 | Agreement among the k strongest discovery age genes, k = 50 to 1,000 (the choice k = 200 was made after the first results) | results/x_topk_sweep.csv |
| X28 | Permuted-age null for the agreement between GTEx fibroblasts and GTEx skin | results/x_gtex_pair_null.json |
| X29 | Per-family Holm and Benjamini–Hochberg values for every exploratory test; report table with one row per check | results/x_exploratory_tests.csv, results/trace_report.csv |

### Amendment 6 (3 November 2026)

| ID | Analysis | Outputs |
|---|---|---|
| X30 | Permuted-age null (1,000 permutations) for agreement among the k strongest discovery genes at every k (50 to 1,000) and for the sign agreement of signature genes | results/x_cohort_topk_null.csv, x_cohort_checks_v2.json |
| X31 | Hidden factors computed after a sex-only model (F9b), their correlation with age, and top-200 agreement with its own permuted-age null | results/x_cohort_checks_v2.json |
| X32 | Planted shared effect added to the real (not permuted) data of both cohorts (planted part alone with an absolute t of about 4 and 6; 200 replicates each) | results/x_cohort_planted_real.csv, x_cohort_checks_v2.json |
| X33 | Variance inflation factor of age in the GTEx batch models; Box 1 checks applied to this study (pass, fail or flag per check); duplicate rows removed from the exploratory-test table | results/x_cohort_checks_v2.json, trace_report.csv, x_exploratory_tests.csv |
