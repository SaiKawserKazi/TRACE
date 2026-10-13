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
