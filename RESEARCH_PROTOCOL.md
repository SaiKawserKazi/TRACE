# TRACE research protocol (v1.0)

**TRACE** — Transferable Rejuvenation Assessment across Contexts and Evidence.

This protocol is fixed before any perturbation response is scored against an aging
signature. It is committed to git and tagged `protocol-v1.0`. Any later deviation is
listed in the "Deviations" section of the manuscript with a reason.

## 1. Question

Transcriptomic "rejuvenation" screens rank interventions by how strongly they reverse an
aging gene-expression signature. We ask: (i) do these rankings reproduce across
independent donor cohorts, across cell contexts and across experiments; (ii) how much of
the apparent rejuvenation is explained by proliferation, identity loss or stress
programs rather than by aging biology; (iii) does biological-network information add
measurable predictive value; and (iv) do robust candidates agree with independent
experimental longevity and senescence evidence?

## 2. Data (all public)

* Aging discovery cohort: primary dermal fibroblasts, Fleischer et al. 2018 (GEO
  GSE113957). Include healthy donors aged >= 18 years; exclude Hutchinson-Gilford
  progeria samples and any sample without numeric age.
* Aging replication cohorts: GTEx v10 cultured fibroblasts, and GTEx v10 skin (not sun
  exposed). Age is available as 10-year bins; bin midpoints are used.
* Perturbation atlases (CRISPRi Perturb-seq pseudobulk): K562 essential-gene screen,
  K562 genome-wide screen, RPE1 essential-gene screen (Replogle et al. 2022); HepG2 and
  Jurkat essential-gene screens (Nadig et al. 2025), pseudobulked by us with the
  Replogle normalisation (per-gem-group z-score against non-targeting controls, then
  mean per perturbation). If HepG2/Jurkat cannot be processed with available memory or
  bandwidth this is reported and the analysis runs on the remaining contexts.
* Network: STRING v12 human, combined score >= 700.
* Gene sets: MSigDB Hallmark (via Enrichr).
* Evidence never used for fitting or selection: CellAge (genes whose experimental
  manipulation induces or inhibits cellular senescence), GenAge model-organism genes
  (pro-/anti-longevity) mapped to human orthologs, OncoKB cancer gene list
  (tumour-suppressor annotation for safety flags).

## 3. Aging signature (Aim 1)

* Expression: log2(FPKM/TPM + 1); genes expressed (mean log2 >= 1) in the cohort.
* Model per gene: expression ~ age + sex (OLS). Age t-statistic `a_g`.
* Discovery signature `S`: Fleischer genes with BH-FDR < 0.05, capped at the 250 most
  up- and 250 most down-regulated genes by t-statistic (fixed a priori).
* Replication: Spearman correlation of age t-statistics between Fleischer and each GTEx
  cohort over shared expressed genes; sign concordance of `S` genes in GTEx.
* Donor robustness: 200 donor-level bootstrap resamples of the Fleischer cohort; the
  signature and all downstream perturbation scores are recomputed per resample.

## 4. Rejuvenation score

For perturbation `p` in context `c` with pseudobulk z-shift vector `d_pc`:

    R_pc = - SpearmanCorr( d_pc[g], a_g ),   g in S ∩ measured genes

`R > 0` means the knockdown moves expression opposite to the aging direction. The
perturbation's own target gene is removed from `S` when scoring it. Non-targeting
controls (when present as pseudobulk entries) are the negative control.

Variants (fixed a priori):
* `R_raw`: as above.
* `R_adj` (confound-adjusted): genes in Hallmark E2F targets, G2M checkpoint, Mitotic
  spindle, and MYC targets V1/V2 are removed from `S`, and `R` is residualised (across
  perturbations within a context) on the perturbation's proliferation-program shift,
  identity-program shift and total response magnitude (mean |z|).

Program shifts (per perturbation, per context): mean z-shift over
* proliferation: Hallmark E2F targets ∪ G2M checkpoint;
* stress/inflammation: Hallmark TNFA signalling via NFKB ∪ Inflammatory response ∪ p53
  pathway;
* identity: the 100 genes with the highest control expression specificity for that
  context (z-score of the context's mean relative to the other contexts, computed on
  control/non-targeting expression only; for single-context files, the top-100 most
  expressed measured genes).

## 5. Hypotheses and decision rules

All tests two-sided unless stated; Holm correction within each family; alpha = 0.05.

* **H1 (cross-context reproducibility).** For shared perturbations, Spearman correlation
  of `R_raw` between every pair of distinct cell lines is > 0 (one-sided, Holm over
  pairs). *Supported* if all pairs are significant.
* **H1b (context vs. experiment).** Concordance of `R_raw` between two experiments in the
  same line (K562 essential vs. K562 genome-wide, shared perturbations) exceeds the mean
  between-line concordance. Tested by bootstrap over perturbations (2,000 resamples),
  one-sided. *Supported* if the 95% one-sided bootstrap interval of the difference
  excludes 0.
* **H2 (confounding).** Across perturbations, `R_raw` is correlated with the
  proliferation-program shift in each context (Spearman, Holm over contexts).
  *Supported* if significant in all contexts. Secondary (descriptive): overlap of top-5%
  candidates between `R_raw` and `R_adj`.
* **H3 (network value).** Task: predict a held-out perturbation's `R_adj` in a target
  context. Predictors (5-fold CV over perturbations, seed 0):
  (a) Transfer: mean `R_adj` of the same perturbation in other contexts;
  (b) Network: mean `R_adj` of the target gene's STRING neighbours in the training folds
  of the target context (weighted by STRING score);
  (c) Transfer + network: ridge regression on (a) and (b) fitted in training folds;
  null: (b) with degree-preserving randomised network (10 randomisations).
  H3a: network (b) beats the randomised-network null (paired Wilcoxon on absolute
  error over perturbations, one-sided, Holm over contexts). H3b: (c) beats (a)
  (same test). *Supported* if significant in a majority of contexts.
* **H4 (independent evidence).** Robust score `R_rob` = mean `R_adj` over contexts in
  which the perturbation was measured, divided by (1 + SD across contexts) — fixed
  a priori. H4a: CellAge genes whose activity *induces* senescence have higher `R_rob`
  than genes whose activity *inhibits* senescence (one-sided Mann-Whitney).
  H4b: human orthologs of GenAge *anti-longevity* genes have higher `R_rob` than
  orthologs of *pro-longevity* genes (one-sided Mann-Whitney). Holm over H4a/H4b.
  Exploratory: same tests with single-context `R_raw`.
* **H5 (donor robustness).** Across 200 donor bootstraps, the median Spearman
  correlation between bootstrap and full-cohort `R_rob` rankings is >= 0.8
  (descriptive threshold, fixed a priori), and the ranking from the GTEx-fibroblast
  signature correlates positively with the Fleischer ranking (one-sided Spearman).

## 6. Multi-objective prioritisation

Objectives per perturbation (maximise): `R_rob`; cross-context sign consistency
(fraction of contexts with `R_adj > 0`); identity preservation (−|identity shift|);
low proliferation gain (−max(0, proliferation shift)); low stress induction
(−max(0, stress shift)). Report the Pareto-non-dominated set among perturbations
measured in >= 2 contexts, with bootstrap support (fraction of donor bootstraps in
which the perturbation is in the top 5% by `R_rob`). Safety annotations: OncoKB
tumour suppressor (knockdown risk) and common-essential status (essential in all
screened lines) are reported as flags, never as evidence of safety. Missing evidence
is labelled "unknown".

## 7. Leakage controls

* The aging signature uses only the aging cohorts; no perturbation data.
* CellAge, GenAge and OncoKB are never used in fitting, scoring or selection.
* All thresholds above are fixed in this document.
* The perturbation's own target is excluded from its score.
* Network predictor uses only training-fold perturbations of the target context.

## 8. Claims

No claim of lifespan extension or therapeutic efficacy. Candidates are hypotheses for
experimental testing.
