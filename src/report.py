"""TRACE report: one row per check in Box 1 of the paper, with its value, the threshold it is judged against and the
result. This is the summary a user reads after `python run_all.py`; it is written to results/trace_report.csv.

Results are "pass" (the ranking can be trusted on this point), "fail", "flag" (a warning that needs reporting, not a
failure) or "info" (reported without a rule). The thresholds are suggestions set after this study; they have not been
tested on other screens."""
import json
import pandas as pd
from config import RES

J = lambda n: json.load(open(RES / n, encoding="utf-8"))
CSV = lambda n: pd.read_csv(RES / n)
SCR = {"K562": "K562 ess.", "K562gw": "K562 gw", "RPE1": "RPE1", "HepG2": "HepG2", "Jurkat": "Jurkat"}


def main():
    R = J("results.json"); NU = J("x_cohort_nulls.json"); CK = J("x_cohort_checks_v2.json"); H5M = J("rev3_h5_multiline.json"); rows = []

    def add(item, check, value, threshold, result):
        rows.append({"item": item, "check": check, "value": value, "threshold": threshold, "result": result})
    pf_ = lambda ok: "pass" if ok else "fail"

    # 1. Cohorts
    pp = CSV("x_cohort_perm_p.csv"); S = CSV("rev3_cohort_sensitivity.csv")
    p0 = pp[pp.spec.str.startswith("F0") & (pp.tissue == "fibroblasts")].iloc[0]
    g12 = [S[S.discovery_spec.str.startswith("F0") & S.replication_spec.str.startswith(f"fibroblasts: {g}")].top_rho.iloc[0] for g in ("G1", "G2")]
    add("Cohorts", "Age effects, 200 strongest discovery genes, versus GTEx fibroblasts (Spearman ρ; main model, then G1 and G2)",
        f"{p0.top_rho:.2f} (G1 {g12[0]:.2f}, G2 {g12[1]:.2f})", f"above the permuted-age 97.5th percentile ({p0.null_top_hi:.2f}) and positive under every replication model",
        pf_(p0.top_rho > p0.null_top_hi and min(g12) > 0))
    pl = CK["planted_real"][str(min(CK["plant_sizes"]))]["fibroblasts"]
    add("Cohorts", f"Planted shared effect added to the real data (planted-only |t| about {CK['planted_real'][str(min(CK['plant_sizes']))]['median_abs_t_planted_only']:.0f}): runs with top-200 ρ above the null",
        f"{pl['frac_top_rho_beyond_null']:.2f}", "at least 0.80", pf_(pl["frac_top_rho_beyond_null"] >= 0.8))
    for t in ("fibroblasts", "skin"):
        sx = NU["sex_perm"][f"sex_top_rho_{t}"]
        add("Cohorts", f"Positive control: sex effects, 200 most sex-biased genes, versus GTEx {t} (permuted-sex P)", f"ρ = {sx['observed']:.2f}; P = {sx['p']:.3f}", "P < 0.05", pf_(sx["p"] < 0.05))
    # 2. Programs
    rs = CSV("x_random_signature_v2.csv").set_index("statistic"); h2 = [R["H2"]["contexts"][c]["rho_raw_prolif"] for c in SCR]
    spec = all(rs.loc[f"H2_rho_{c}", "frac_null_ge_real"] < 0.05 for c in SCR)
    add("Programs", "Raw score versus proliferation shift (Spearman ρ, range over screens); real signature versus 100 random signatures",
        f"{min(h2):.2f} to {max(h2):.2f}; above all random signatures in every screen" if spec else f"{min(h2):.2f} to {max(h2):.2f}",
        "flag if |ρ| ≥ 0.3 in any screen or the real signature exceeds 95% of random signatures", "flag" if max(map(abs, h2)) >= 0.3 or spec else "pass")
    # 3. Replicates
    hb = R["H1b"]; rel = CSV("rev3_reliability_summary.csv")
    rv = {v: rel[rel.variant == v].reliability_SB for v in ("V0 raw", "V1 pre-specified")}
    add("Replicates", "Raw score: within-line minus between-line agreement (one-sided 95% bound); score reliability (range over screens)",
        f"{hb['diff_ci_lower_one_sided95']:.2f}; {rv['V0 raw'].min():.2f} to {rv['V0 raw'].max():.2f}", "bound > 0 and reliability ≥ 0.5 in every screen",
        pf_(hb["diff_ci_lower_one_sided95"] > 0 and rv["V0 raw"].min() >= 0.5))
    add("Replicates", "Adjusted score (V1): reliability (range over screens)", f"{rv['V1 pre-specified'].min():.2f} to {rv['V1 pre-specified'].max():.2f}", "≥ 0.5 in every screen",
        pf_(rv["V1 pre-specified"].min() >= 0.5))
    # 4. Labels
    ev = CSV("rev3_scores_evidence.csv"); H4 = R["H4"]
    pro = lambda lab: ev[(ev.labels == lab) & (ev.gene_set == "all genes") & (ev.score == "proliferation shift alone")].auroc.iloc[0]
    for key, lab in [("H4a_cellage", "CellAge"), ("H4b_genage", "GenAge")]:
        add("Labels", f"{lab}: Rrob (V1) AUROC, all genes (Holm P), and proliferation shift alone", f"{H4[key]['auroc']:.2f} (P = {H4[key]['p_holm']:.2f}); baseline {pro(lab):.2f}",
            "AUROC above 0.5 (corrected P < 0.05) and above the proliferation-only AUROC (paired corrected P < 0.05)", pf_(H4[key]["p_holm"] < 0.05 and H4[key]["auroc"] > pro(lab)))
    # 5. Power
    mda = lambda lab: ev[(ev.labels == lab) & (ev.gene_set == "all genes") & (ev.score == "R_rob (V1, pre-specified)")].min_detectable_auroc_80pct.iloc[0]
    for lab, nm in [("CellAge", "CellAge"), ("GenAge", "GenAge"), ("GenAge mammal-only", "GenAge, mouse genes only")]:
        add("Power", f"Minimum detectable AUROC at 80% power, {nm}", f"{mda(lab):.2f}", "a label test above 0.60 is reported as uninformative", "info" if mda(lab) <= 0.60 else "uninformative")
    # 6. Effect sizes
    boot = H5M["all"]["median_top5"]
    for t, key in [("GTEx fibroblasts", "gtex_fibroblast"), ("GTEx skin", "gtex_skin")]:
        add("Effect sizes", f"Top-5% overlap between the discovery ranking and the ranking from {t}", f"{R['H5'][key]['top5_overlap']:.2f}",
            f"at least half the overlap between donor bootstraps ({boot:.2f} / 2 = {boot / 2:.2f})", pf_(R["H5"][key]["top5_overlap"] >= boot / 2))
    # 7. Knockdown quality
    s3 = J("rev3_scores.json")["E14_front"]; pf = CSV("rev3_scores_pareto_filtered.csv"); weak = int(pf[pf.pareto].flag_weak_or_unmeasured_knockdown.sum())
    add("Knockdown quality", "Pareto-front genes with weak or unmeasured on-target knockdown before filtering; front after filtering",
        f"{weak} of {s3['unfiltered']['n']}; {s3['filtered']['n']} genes after filtering", "report the share, then rank only after filtering", "info")
    # 8. Stability
    add("Stability", "Median donor-bootstrap ρ of the ranking (all genes; genes in ≥ 2 lines)", f"{R['H5']['median_rho']:.3f}; {H5M['multi']['median_rho']:.3f}", "≥ 0.8",
        f"{pf_(R['H5']['median_rho'] >= 0.8)} (all genes); {pf_(H5M['multi']['median_rho'] >= 0.8)} (genes in ≥ 2 lines)")
    out = pd.DataFrame(rows); out.to_csv(RES / "trace_report.csv", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
