"""TRACE report: one row per audit check, with its value, the reference it is judged against, the decision rule and
pass/fail. This is the summary a user reads after `python run_all.py`; it is written to results/trace_report.csv.

The rules are the ones in Box 1 of the paper. A "fail" means the ranking should not be trusted on that point without
further evidence; "info" rows are reported without a rule."""
import json
import pandas as pd
from config import RES

J = lambda n: json.load(open(RES / n, encoding="utf-8"))
CSV = lambda n: pd.read_csv(RES / n)
SCR = {"K562": "K562 ess.", "K562gw": "K562 gw", "RPE1": "RPE1", "HepG2": "HepG2", "Jurkat": "Jurkat"}


def main():
    R = J("results.json"); NU = J("x_cohort_nulls.json"); rows = []
    add = lambda check, value, reference, rule, ok: rows.append({"check": check, "value": value, "reference": reference, "rule": rule,
                                                                "result": "info" if ok is None else ("pass" if ok else "fail")})
    # 1. cohort transfer among the strongest age genes, with its power check
    pp = CSV("x_cohort_perm_p.csv"); p0 = pp[pp.spec.str.startswith("F0") & (pp.tissue == "fibroblasts")].iloc[0]
    add("Age effects, 200 strongest genes, discovery vs second cohort (Spearman ρ)", round(p0.top_rho, 3),
        f"permuted-age 95% range {p0.null_top_lo:.2f} to {p0.null_top_hi:.2f}", "ρ above the 97.5th percentile of the permuted-age null", p0.top_rho > p0.null_top_hi)
    pl = NU["planted"]["0.2"]["fibroblasts"]
    add("Power check: planted shared effect recovered among the 200 strongest genes (fraction of runs)", pl["frac_top_rho_beyond_null"],
        f"median ρ {pl['median_top_rho']:.2f}", "recovered in at least 80% of runs", pl["frac_top_rho_beyond_null"] >= 0.8)
    add("Power check: planted shared effect recovered genome-wide (fraction of runs)", pl["frac_rho_all_beyond_null"], "",
        "genome-wide agreement is uninformative when this is below 80%", None)
    sx = NU["sex_perm"]["sex_top_rho_fibroblasts"]
    add("Positive control: sex effects, 200 most sex-biased genes (Spearman ρ)", round(sx["observed"], 3), f"permuted-sex P = {sx['p']:.3f}", "P < 0.05", sx["p"] < 0.05)
    # 2. ranking transfer between cohorts
    g = R["H5"]["gtex_fibroblast"]
    add("Ranking from a second-cohort signature: top-5% overlap", round(g["top5_overlap"], 3), "5% expected by chance", "overlap clearly above chance (at least 2 × 5%)", g["top5_overlap"] >= 0.10)
    add("Ranking from a second-cohort signature: Spearman ρ", round(g["rho"], 3), "", "reported", None)
    # 3. programs inside the signature
    for c in SCR:
        h = R["H2"]["contexts"][c]
        add(f"Raw score vs proliferation shift, {SCR[c]} (Spearman ρ)", round(h["rho_raw_prolif"], 3), "", "reported per screen; |ρ| ≥ 0.3 means the ranking mostly follows growth", None)
    rs = CSV("x_random_signature_v2.csv").set_index("statistic")
    add("Real signature vs 100 random signatures: adjusted cross-line agreement", round(rs.loc["H1_median_rho_adj", "real_V2_input"], 3),
        f"null maximum {rs.loc['H1_median_rho_adj', 'null_max']:.2f}", "above the 95th percentile of the random signatures", rs.loc["H1_median_rho_adj", "frac_null_ge_real"] < 0.05)
    # 4. replicates
    rel = CSV("rev3_reliability_summary.csv"); v1 = rel[rel.variant == "V1 pre-specified"]
    for _, r in v1.iterrows():
        add(f"Reliability of the adjusted score, {SCR[r.screen]}", round(r.reliability_SB, 3), "", "reported", None)
    # 5. labels against a proliferation-only baseline
    pa = CSV("rev3_scores_paired_auroc.csv")
    for lab in ["CellAge", "GenAge"]:
        x = pa[(pa.labels == lab) & (pa.gene_set == "all genes") & (pa.a == "mean raw score") & (pa.b == "proliferation shift alone")].iloc[0]
        add(f"{lab}: raw score minus proliferation-only baseline (AUROC difference)", round(x.diff_a_minus_b, 3), f"paired DeLong P = {x.p_two_sided:.2g}",
            "score beats the baseline (difference > 0, P < 0.05)", bool(x.diff_a_minus_b > 0 and x.p_two_sided < 0.05))
    # 6. knockdown quality before multi-objective ranking
    s3 = J("rev3_scores.json")["E14_front"]; pf = CSV("rev3_scores_pareto_filtered.csv")
    weak = int(pf[pf.pareto].flag_weak_or_unmeasured_knockdown.sum())
    add("Pareto front genes with weak or unmeasured on-target knockdown", f"{weak} of {s3['unfiltered']['n']}", "", "none after filtering on knockdown and response size", weak == 0)
    out = pd.DataFrame(rows); out.to_csv(RES / "trace_report.csv", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
