# Protocol amendments and deviations (public log)

The protocol `RESEARCH_PROTOCOL.md` was frozen at commit `6ab7550` (tag `protocol-v1.0`, local time 2026-10-12 16:08:26 +0600)
and has not been modified since. Its SHA-256 and git blob ID are in `deposit/PREREGISTRATION_RECORD.json`.
Commit times are recorded by local version control and are self-attested. The public archive verifies the content of the
protocol, but it cannot independently prove the freeze time.

**Confirmatory** = the pre-specified test in protocol v1.0, run as specified (H1–H5, Table 1 of the manuscript).
**Exploratory** = everything else, including every analysis listed under the amendments below.

| ID | Change | When (relative to results) | Possible influence |
|---|---|---|---|
| D1 | Orthologs from the Alliance of Genome Resources instead of Ensembl BioMart (service returned HTTP 405) | 2026-10-12, before any H4 result | Changes the mapped label set; not tuned to results |
| D2 | The K562 essential-gene and genome-wide screens are averaged into one line for R_rob and the Pareto ranking (protocol silent) | Before results | Defines R_rob |
| D3 | H1 uses the K562 essential-gene screen as "K562" (protocol silent) | Before results | None expected |
| D4 | HepG2/Jurkat were pseudobulked by us; their non-targeting controls are 100 random groups of non-targeting cells | Before those screens' results | Normalisation may differ slightly from the published release |
| D5 | Ridge penalty λ = 1 for the H3 models (protocol said "ridge" without λ) | Before results | Untuned; may understate the combined model |
| D6 | Bug fix: an ortholog-file size threshold excluded yeast from H4b | 2026-10-12, **after** first H4 results | H4b AUROC 0.52 (P_Holm = 0.184) before → 0.51 (0.270) after; conclusion unchanged |
| D7 | First-draft exploratory additions: 18–80-year sensitivity, adjusted within-K562 agreement, calibrated transfer-only model, single-screen evidence tests, on-target annotation | 2026-10-12, after pre-specified results | Labelled exploratory |
| AMENDMENT-R1 | Review round 1: GEO metadata; cohort sensitivity (8 × 3 specifications) with donor-bootstrap CIs and split-half references; score definitions V2–V5; network baselines and CIs; evidence coverage, permutation, subsets and proliferation-only baseline; reproducibility package | 2026-10-13, after all pre-specified results | All exploratory |
| AMENDMENT-R2 | Review round 2: 2,000 donor bootstraps for cohort sensitivity (500 split halves); 2,000 donor bootstraps for H5 (the pre-specified 200 is kept as confirmatory); 500 degree-preserving randomised networks for H3a (pre-specified 10 kept as confirmatory); 2,000 bootstrap resamples for all CIs and 10,000 label permutations; GTEx-skin comparison relabelled exploratory (it was never part of H5) | 2026-10-13, after all pre-specified results | All exploratory. Confirmatory outcomes kept as specified; H3a is fragile (2/5 screens against the mean of 500 nulls versus 3/5 with the pre-specified 10) |

## AMENDMENT-R3 (v0.1 internal review by the senior author; logged 2026-10-26 before the analyses were run)

| ID | Change | When | Possible influence |
|---|---|---|---|
| D8 | Bug fix: in AMENDMENT-R1/R2 the GTEx-skin "G0 (age + sex)" model also contained ischemic time. G0 is now age + sex for both GTEx tissues; ischemic time enters skin models from G1. | 2026-10-26, after R2 results | Changes the reported skin G0 values (0.21 → 0.15; fibroblasts vs skin 0.02 → 0.06). Pre-specified analyses (H1–H5) were not affected. |
| D9 | Exclusion flow made explicit: 143 GEO samples = 133 healthy donors + 10 progeria samples; 22 donors under 18 excluded, leaving 111 donors aged 19–96 years (GEO age field). No change to the analysis. | 2026-10-26 | None |
| E1 | Sex as a positive control (same models and bootstraps; all genes, autosomes, X/Y; top-200 genes) | 2026-10-26 | Exploratory |
| E2 | Permuted-age null (1,000 permutations, age shuffled within cohort, covariates fixed) | 2026-10-26 | Exploratory |
| E3 | Sign check of pre-named age genes and canonical sex genes | 2026-10-26 | Exploratory |
| E4 | Hidden-factor adjustment (5 residual principal components); discovery without the 31 NIA donors aged 80+; 1,000 power-matched 74-donor subsets; VIF of age; age per repository | 2026-10-26 | Exploratory |
| E5 | 100 random expression-matched signatures through H1, H2, H4 | 2026-10-26 | Exploratory |
| E6 | Split-cell reliability (HepG2, Jurkat; 20 splits); test-retest reliability for K562 (two screens); size-matched non-targeting groups; disattenuated agreement | 2026-10-26 | Exploratory; RPE1 has no replicate or single-cell data in hand |
| E7 | Paired DeLong tests, incremental logistic regression, minimum detectable AUROC | 2026-10-26 | Exploratory |
| E8/E9 | CellAge Build 3 convention documented (labels are gene-function based); overlaps with core-essential genes (Hart CEGv2; DepMap download is behind a CAPTCHA); GenAge mammal-only and without yeast | 2026-10-26 | Exploratory |
| E10 | Rrob for single-line genes stated (equals Radj); H4 and H5 on genes measured in ≥2 lines | 2026-10-26 | Exploratory |
| E11–E12 | In-fold residualization, inner-CV ridge penalty, 10 CV seeds, common subset per screen; kNN, gene-covariate, physical-network, calibrated-network and proliferation baselines; block bootstrap over gene clusters | 2026-10-26 | Exploratory |
| E13–E16 | Proliferation deciles, distance correlation, mutual information, V6 spline and V7 variants; filtered Pareto front; ranks of senescence genes; up-minus-down score | 2026-10-26 | Exploratory |
| E17 | Cross-cohort elastic-net age prediction | 2026-10-26 | Exploratory |
| E18 | TRACE run end to end with a GTEx-derived discovery signature | 2026-10-26 | Exploratory |
| D10 | Loader fix: load_fleischer read sample columns [8:-1], which dropped the last column (sample 99, age 17). Now all sample columns are read. A generic two-CSV cohort loader (load_table) was added for reuse. | 2026-10-26 | None: the dropped sample is under 18 and excluded by the pre-specified age filter; all results are unchanged |
| Bootstrap count | Revised cohort sensitivity rerun with 2,000 donor bootstraps (as in AMENDMENT-R2), 1,000 permutations and 1,000 power-matched subsets; the count is also stored in rev3_cohort.json (n_boot) | 2026-10-26 | Exploratory |
| E8/E9 update | DepMap 24Q2 CRISPR-inferred common essentials obtained from the official Figshare+ release (doi:10.25452/figshare.plus.25880521.v1); overlaps and Fisher tests added next to Hart CEGv2 | 2026-10-27 | Exploratory |
| Table S3 intervals | Bias-corrected (BC) bootstrap intervals computed from the saved 2,000 donor bootstraps and reported next to percentile intervals; BCa acceleration not computed (needs a jackknife refit per donor) | 2026-10-27 | Exploratory |
| E6 update | RPE1 split-cell reliability (20 splits) from the Replogle RPE1 normalized single-cell release (Figshare 20029387, file 35775554; cells already z-scored per gem group, so the standardization step is skipped); RPE1 pairs added to the disattenuated table | 2026-10-27, before running | Exploratory |
