"""TRACE report: one row per check in Box 1 of the paper, with its value, the threshold it is judged against and the
result, for the discovery signature and for the GTEx-derived signature. This is the summary a user reads after
`python run_all.py`; it is written to results/trace_report.csv.

Results are "pass" (the ranking can be trusted on this point), "fail", "flag" (a warning that needs reporting, not a
failure), "info" (reported without a rule) or "not run" (the input for this check was not computed for this signature).
"Corrected P" means Holm-adjusted P for tests fixed before the study and Benjamini-Hochberg q for the others, within each
check. The stability and label rules follow tests fixed before this study (protocol v1.0: H5a, H4); the other thresholds
are suggestions set after this study and have not been tested on other screens. A ranking that fails a check should not be
used to choose targets; its top candidates are hypotheses only."""
import json
import pandas as pd
from statsmodels.stats.multitest import multipletests
from config import RES

J = lambda n: json.load(open(RES / n, encoding="utf-8"))
CSV = lambda n: pd.read_csv(RES / n)
SCR = {"K562": "K562 ess.", "K562gw": "K562 gw", "RPE1": "RPE1", "HepG2": "HepG2", "Jurkat": "Jurkat"}
NR = ("–", "not run")


def main():
    R = J("results.json"); E = J("e18_gtex_discovery/results.json"); NU = J("x_cohort_nulls.json"); CK = J("x_cohort_checks_v2.json")
    H5M = J("rev3_h5_multiline.json"); rows = []

    def add(item, check, value, threshold, result, gx=NR):
        rows.append({"item": item, "check": check, "value": value, "threshold": threshold, "result": result, "value_gtex": gx[0], "result_gtex": gx[1]})
    pf_ = lambda ok: "pass" if ok else "fail"
    # exploratory tests are corrected within their Table S5 families (BH q), so the report reads q values from there
    XT = CSV("x_exploratory_tests.csv")
    xq = lambda fam, test: float(XT[(XT.Family == fam) & (XT.Test == test)]["BH q (within family)"].iloc[0])

    # 1. Cohorts. Discovery: the 200 strongest discovery genes against GTEx fibroblasts; GTEx-derived: the 200 strongest
    # GTEx-fibroblast genes against the discovery cohort (same permutations). Sex control and planted effect use the same
    # cohorts and pipeline for both signatures.
    pp = CSV("x_cohort_perm_p.csv"); S = CSV("rev3_cohort_sensitivity.csv")
    p0 = pp[pp.spec.str.startswith("F0") & (pp.tissue == "fibroblasts")].iloc[0]
    g12 = [S[S.discovery_spec.str.startswith("F0") & S.replication_spec.str.startswith(f"fibroblasts: {g}")].top_rho.iloc[0] for g in ("G1", "G2")]
    rv = J("x_reverse_top200.json")
    add("Cohorts", "Age effects, 200 strongest genes of the signature's cohort, versus the other fibroblast cohort (Spearman ρ; main model; for the discovery signature G1 and G2 for information)",
        f"{p0.top_rho:.2f} (G1 {g12[0]:.2f}, G2 {g12[1]:.2f}); null 97.5th percentile {p0.null_top_hi:.2f}", "above the permuted-age 97.5th percentile under the main model",
        pf_(p0.top_rho > p0.null_top_hi), (f"{rv['rho']:.2f}; null 97.5th percentile {rv['null_hi']:.2f}", pf_(rv["rho"] > rv["null_hi"])))
    k0 = str(min(CK["plant_sizes"])); pl = CK["planted_real"][k0]["fibroblasts"]; tpl = CK["planted_real"][k0]["median_abs_t_planted_only"]
    mt = J("x_topk_f9_check.json")["F0_top200_median_abs_t"]; okp = tpl <= mt and pl["frac_top_rho_beyond_null"] >= 0.8
    add("Cohorts", "Planted shared effect added to the real data: median |t| of the planted part alone (against the median |t| of the observed top 200); runs with top-200 ρ above the null",
        f"|t| {tpl:.1f} (observed {mt:.1f}); {pl['frac_top_rho_beyond_null']:.2f}", "planted |t| no larger than the observed median |t|, recovered in at least 0.80 of runs",
        pf_(okp), NR)   # planted in the discovery direction only; not run for the GTEx-derived signature
    sxp = [NU["sex_perm"][f"sex_top_rho_{t}"] for t in ("fibroblasts", "skin")]
    sxq = [xq("Sex positive control (permuted-sex null)", f"Sex effects, 200 most sex-biased genes, discovery vs GTEx {t}") for t in ("fibroblasts", "skin")]
    for t, sx, q in zip(("fibroblasts", "skin"), sxp, sxq):
        add("Cohorts", f"Positive control: sex effects, 200 most sex-biased genes, versus GTEx {t} (permuted-sex P; BH q within the sex-control family, Table S5)",
            f"ρ = {sx['observed']:.2f}; P = {sx['p']:.3f}; q = {q:.2g}", "corrected P < 0.05", pf_(q < 0.05), ("same cohorts and pipeline", pf_(q < 0.05)))
    # 2. Programs
    rs = CSV("x_random_signature_v2.csv").set_index("statistic"); h2 = [R["H2"]["contexts"][c]["rho_raw_prolif"] for c in SCR]
    e2 = [E["H2"]["contexts"][c]["rho_raw_prolif"] for c in SCR]
    spec = all(rs.loc[f"H2_rho_{c}", "frac_null_ge_real"] < 0.05 for c in SCR)
    add("Programs", "Raw score versus proliferation shift (Spearman ρ, range over screens); real signature versus 100 random signatures (discovery only)",
        f"{min(h2):.2f} to {max(h2):.2f}; above all random signatures in every screen" if spec else f"{min(h2):.2f} to {max(h2):.2f}",
        "flag if |ρ| ≥ 0.3 in any screen or the real signature exceeds 95% of random signatures", "flag" if max(map(abs, h2)) >= 0.3 or spec else "pass",
        (f"{min(e2):.2f} to {max(e2):.2f}", "flag" if max(map(abs, e2)) >= 0.3 else "pass"))
    # 3. Replicates: reliability of the adjusted score in every screen and adjusted between-line agreement above the
    # random-signature null; within- and between-line agreement reported side by side
    hb, eb = R["H1b"], E["H1b"]
    add("Replicates", "Raw score: within-line (two K562 screens) and mean between-line agreement, side by side (one-sided 95% bound of the difference)",
        f"{hb['rho_within_K562']:.2f} and {hb['rho_between_mean']:.2f} ({hb['diff_ci_lower_one_sided95']:.2f})", "reported", "info",
        (f"{eb['rho_within_K562']:.2f} and {eb['rho_between_mean']:.2f} ({eb['diff_ci_lower_one_sided95']:.2f})", "info"))
    rel = CSV("rev3_reliability_summary.csv"); r1 = rel[rel.variant == "V1 pre-specified"].reliability_SB
    add("Replicates", "Adjusted score (V1): reliability (range over screens; split cells in HepG2 and Jurkat, test–retest in K562)", f"{r1.min():.2f} to {r1.max():.2f}", "≥ 0.5 in every screen",
        pf_(r1.min() >= 0.5), (f"K562 test–retest {eb['rho_adj_within_allshared']:.2f} (one screen only); split cells not run", pf_(eb["rho_adj_within_allshared"] >= 0.5)))
    ra = rs.loc["H1_median_rho_adj"]
    add("Replicates", "Adjusted score: median between-line agreement against 100 random signatures (real signature without proliferation genes; one-sided empirical P)",
        f"{ra.real_V2_input:.2f} (null maximum {ra.null_max:.2f}; P = {ra.empirical_p_one_sided:.3f})", "above the random-signature null (P < 0.05)", pf_(ra.empirical_p_one_sided < 0.05))
    # 4. Labels
    ev = CSV("rev3_scores_evidence.csv"); gp = CSV("x_gtex_signature_paired_auroc.csv")
    # the GTEx-derived label tests are exploratory, so Box 1's rule uses BH q within their Table S5 family
    pro = lambda lab: ev[(ev.labels == lab) & (ev.gene_set == "all genes") & (ev.score == "proliferation shift alone")].auroc.iloc[0]
    for key, lab in [("H4a_cellage", "CellAge"), ("H4b_genage", "GenAge")]:
        d, e = R["H4"][key], E["H4"][key]
        g = gp[(gp.labels == lab) & (gp.gene_set == "all genes") & (gp.a == "R_rob (V1)") & (gp.b == "proliferation shift alone")].iloc[0]
        qe = xq("Held-out labels, GTEx-derived signature", f"{lab}, all genes, R_rob (V1)")
        qd = xq("Paired AUROC differences (DeLong), GTEx-derived signature", f"{lab}, all genes: R_rob (V1) minus proliferation shift alone")
        ok_e = qe < 0.05 and e["auroc"] > pro(lab) and qd < 0.05
        add("Labels", f"{lab}: Rrob (V1) AUROC, all genes, and proliferation shift alone (discovery: Holm P, fixed in advance; GTEx-derived: BH q, exploratory, with paired DeLong q)",
            f"{d['auroc']:.2f} (P = {d['p_holm']:.2f}); baseline {pro(lab):.2f}",
            "AUROC above 0.5 (corrected P < 0.05) and above the proliferation-only AUROC (paired DeLong, corrected P < 0.05)", pf_(d["p_holm"] < 0.05 and d["auroc"] > pro(lab)),
            (f"{e['auroc']:.2f} (q = {qe:.2g}); baseline {pro(lab):.2f} (paired q = {qd:.2g})", pf_(ok_e)))
    # 5. Power (the label sets and genes are the same for both signatures)
    mda = lambda lab: ev[(ev.labels == lab) & (ev.gene_set == "all genes") & (ev.score == "R_rob (V1, pre-specified)")].min_detectable_auroc_80pct.iloc[0]
    for lab, nm in [("CellAge", "CellAge"), ("GenAge", "GenAge"), ("GenAge mammal-only", "GenAge, mouse genes only")]:
        res = "info" if mda(lab) <= 0.60 else "uninformative"
        add("Power", f"Minimum detectable AUROC at 80% power, {nm}", f"{mda(lab):.2f}", "a label test above 0.60 is reported as uninformative", res, (f"{mda(lab):.2f}", res))
    # 6. Effect sizes
    boot, eboot = H5M["all"]["median_top5"], E["H5"]["median_top5_overlap"]
    add("Effect sizes", "Top-5% overlap between the ranking and the ranking from the other fibroblast cohort",
        f"{R['H5']['gtex_fibroblast']['top5_overlap']:.2f} (bootstraps {boot:.2f})", "at least half the overlap between donor bootstraps of one cohort",
        pf_(R["H5"]["gtex_fibroblast"]["top5_overlap"] >= boot / 2),
        (f"{R['H5']['gtex_fibroblast']['top5_overlap']:.2f} (bootstraps {eboot:.2f})", pf_(R["H5"]["gtex_fibroblast"]["top5_overlap"] >= eboot / 2)))
    add("Effect sizes", "Top-5% overlap between the ranking and the ranking from GTEx skin",
        f"{R['H5']['gtex_skin']['top5_overlap']:.2f}", "at least half the overlap between donor bootstraps of one cohort", pf_(R["H5"]["gtex_skin"]["top5_overlap"] >= boot / 2),
        (f"{E['H5']['gtex_skin']['top5_overlap']:.2f}", pf_(E["H5"]["gtex_skin"]["top5_overlap"] >= eboot / 2)))
    # 7. Knockdown quality
    s3 = J("rev3_scores.json")["E14_front"]; pf = CSV("rev3_scores_pareto_filtered.csv"); weak = int(pf[pf.pareto].flag_weak_or_unmeasured_knockdown.sum())
    add("Knockdown quality", "Pareto-front genes with weak or unmeasured on-target knockdown before filtering; front after filtering",
        f"{weak} of {s3['unfiltered']['n']}; {s3['filtered']['n']} genes after filtering", "report the share, then rank only after filtering", "info",
        (f"{E['pareto']['n_front_weak_kd']} of {E['pareto']['n_front']}; filtering not run", "info"))
    # 8. Stability
    add("Stability", "Median donor-bootstrap ρ of the ranking (all genes; genes in ≥ 2 lines)", f"{R['H5']['median_rho']:.3f}; {H5M['multi']['median_rho']:.3f}", "≥ 0.8",
        f"{pf_(R['H5']['median_rho'] >= 0.8)} (all genes); {pf_(H5M['multi']['median_rho'] >= 0.8)} (genes in ≥ 2 lines)",
        (f"{E['H5']['median_rho']:.3f}; not run", f"{pf_(E['H5']['median_rho'] >= 0.8)} (all genes)"))
    out = pd.DataFrame(rows); out.to_csv(RES / "trace_report.csv", index=False)
    print(out.drop(columns=["check", "threshold"]).to_string(index=False))


if __name__ == "__main__":
    main()
