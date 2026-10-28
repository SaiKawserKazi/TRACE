"""Catalog of every exploratory test with its P value and the size of its family (Table: all exploratory tests).

Reads the exploratory result files and writes results/x_exploratory_tests.csv. AUROC tests without a stored P value get a
two-sided Mann-Whitney normal approximation from the AUROC and the class sizes. No P value here is corrected unless the
column p_holm is filled; the family size is given so that readers can apply their own correction."""
import json
import numpy as np, pandas as pd
from scipy import stats
from config import RES

J = lambda n: json.load(open(RES / n, encoding="utf-8"))
CSV = lambda n: pd.read_csv(RES / n)
rows = []


def add(family, test, stat_name, stat, p, n=None, p_holm=np.nan, source=""):
    rows.append({"family": family, "test": test, "statistic": stat_name, "value": stat, "p": p, "p_holm": p_holm, "n": n, "source": source})


def p_auroc(a, n1, n2):
    se = np.sqrt((n1 + n2 + 1) / (12 * n1 * n2)); return float(2 * stats.norm.sf(abs(a - 0.5) / se))


# cohort agreement: permuted-age null per discovery specification (G0), genome-wide and top-200
for _, r in CSV("x_cohort_perm_p.csv").iterrows():
    add("Cohort agreement, permuted-age null", f"{r.spec} vs GTEx {r.tissue}, all genes", "Spearman rho", r.rho_all, r.p_rho_all, source="x_cohort_perm_p.csv")
    add("Cohort agreement, permuted-age null", f"{r.spec} vs GTEx {r.tissue}, 200 strongest genes", "Spearman rho", r.top_rho, r.p_top_rho, source="x_cohort_perm_p.csv")
# sex control, permuted-sex null
for k, v in J("x_cohort_nulls.json")["sex_perm"].items():
    if isinstance(v, dict): add("Sex positive control, permuted-sex null", k, "Spearman rho", v["observed"], v["p"], source="x_cohort_nulls.json")
# network: 500 nulls
for c, v in J("rev_network_null.json")["screens"].items():
    add("Network vs 500 randomized networks, rank", f"{c}", "Spearman rho", v["rho_network"], v["empirical_p_rho"], v["n_eval"], source="rev_network_null.json")
    add("Network vs mean of 500 randomized networks, absolute error", f"{c}", "MAE", v["mae_network"], v["h3a_wilcoxon_vs_mean_null_p"], v["n_eval"],
        v["h3a_wilcoxon_vs_mean_null_p_holm"], "rev_network_null.json")
# held-out labels: every score, label subset and gene set
for _, r in CSV("rev3_scores_evidence.csv").iterrows():
    add("Held-out labels, discovery signature", f"{r.labels}, {r.gene_set}, {r.score}", "AUROC", r.auroc, p_auroc(r.auroc, r.n_pos, r.n_neg), f"{r.n_pos}/{r.n_neg}", source="rev3_scores_evidence.csv")
for _, r in CSV("x_gtex_signature_evidence.csv").iterrows():
    add("Held-out labels, GTEx-derived signature", f"{r.labels}, {r.gene_set}, {r.score}", "AUROC", r.auroc, p_auroc(r.auroc, r.n_pos, r.n_neg), f"{r.n_pos}/{r.n_neg}", source="x_gtex_signature_evidence.csv")
for _, r in CSV("x_cellage_method.csv").iterrows():
    add("Held-out labels, CellAge by manipulation", f"{r.signature}, {r.subset}, {r.score}", "AUROC", r.auroc, p_auroc(r.auroc, r.n_pos, r.n_neg), f"{r.n_pos}/{r.n_neg}", source="x_cellage_method.csv")
H4 = J("results.json")["H4"]
for key, lab in [("H4a_cellage", "CellAge"), ("H4b_genage", "GenAge")]:
    for c in ["K562", "K562gw", "RPE1", "HepG2", "Jurkat"]:
        for typ in ["raw", "adj"]:
            e = H4[key][f"{typ}_{c}"]
            add("Held-out labels, single screens", f"{lab}, {c}, {'raw score' if typ == 'raw' else 'R_adj (V1)'}", "AUROC", e["auroc"],
                e.get("p", p_auroc(e["auroc"], e["n_pos"], e["n_neg"])), f"{e['n_pos']}/{e['n_neg']}", source="results.json")
# paired AUROC (DeLong) and incremental logistic models
for _, r in CSV("rev3_scores_paired_auroc.csv").iterrows():
    add("Paired AUROC differences (DeLong), discovery signature", f"{r.labels}, {r.gene_set}: {r.a} minus {r.b}", "AUROC difference", r.diff_a_minus_b, r.p_two_sided, f"{r.n_pos}/{r.n_neg}", source="rev3_scores_paired_auroc.csv")
for _, r in CSV("x_gtex_signature_paired_auroc.csv").iterrows():
    add("Paired AUROC differences (DeLong), GTEx-derived signature", f"{r.labels}, {r.gene_set}: {r.a} minus {r.b}", "AUROC difference", r.diff_a_minus_b, r.p_two_sided, f"{r.n_pos}/{r.n_neg}", source="x_gtex_signature_paired_auroc.csv")
for _, r in CSV("rev3_scores_logistic.csv").iterrows():
    add("Incremental logistic model (raw score added to proliferation)", f"{r.labels}, {r.gene_set}", "coefficient", r.coef_raw, r.lr_p_raw_added, r.n, source="rev3_scores_logistic.csv")
for _, r in CSV("rev3_depmap_label_fisher.csv").iterrows():
    add("Label classes and essential genes (Fisher)", f"{r.labels}, {r.reference}: {r.positive_class} vs {r.negative_class}", "odds ratio", r.odds_ratio, r.fisher_p_two_sided, source="rev3_depmap_label_fisher.csv")

T = pd.DataFrame(rows)
T["family_size"] = T.groupby("family").test.transform("size")
T = T[["family", "family_size", "test", "statistic", "value", "n", "p", "p_holm", "source"]]
T.to_csv(RES / "x_exploratory_tests.csv", index=False)
print(len(T), "tests in", T.family.nunique(), "families")
print(T.groupby("family").size().to_string())
