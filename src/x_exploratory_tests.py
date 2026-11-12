"""Catalog of every exploratory test with its P value, the size of its family, and Holm and Benjamini-Hochberg values
within the family (Table: all exploratory tests).

Pre-specified tests (H1-H5) are not listed. Each test appears once (the CellAge-by-manipulation family lists the method subsets only): proliferation-only baselines, which do not depend on
the signature, are listed under the discovery signature only. AUROC tests without a stored P value get a two-sided
Mann-Whitney normal approximation from the AUROC and the class sizes. Supplementary-table references are written as ⟨key⟩ and numbered by the manuscript builder. Writes results/x_exploratory_tests.csv."""
import json
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from config import RES

J = lambda n: json.load(open(RES / n, encoding="utf-8"))
CSV = lambda n: pd.read_csv(RES / n)
SCR = {"K562": "K562 ess.", "K562gw": "K562 gw", "RPE1": "RPE1", "HepG2": "HepG2", "Jurkat": "Jurkat"}
LAB = {"GenAge mammal-only": "GenAge, mouse genes only", "GenAge without yeast": "GenAge without yeast"}
SCORE = {"R_rob (V1, pre-specified)": "R_rob (V1)", "mean raw score": "raw score", "proliferation shift alone": "proliferation shift alone"}
rows = []


def add(family, test, stat_name, stat, p, n="", where=""):
    rows.append({"Family": family, "Test": test, "Statistic": stat_name, "Value": stat, "n": n, "P": p, "Shown in": where})


def hallmark_name(s):
    """Hallmark set names in the text style (E2F targets, G2M checkpoint, TNF-α signaling via NF-κB)."""
    s = " ".join(str(s).split()).replace("Pperoxisome", "Peroxisome").replace("G2-M", "G2M").replace("TNF-alpha", "TNF-α").replace("NF-kB", "NF-κB").replace(
        "Wnt-beta Catenin", "Wnt/β-catenin").replace("TGF-beta", "TGF-β")
    w = s.split(" ")
    return " ".join([w[0][0].upper() + w[0][1:]] + [x if (any(ch.isdigit() for ch in x) or sum(ch.isupper() for ch in x) > 1) else x.lower() for x in w[1:]])


def p_auroc(a, n1, n2):
    se = np.sqrt((n1 + n2 + 1) / (12 * n1 * n2)); return float(2 * stats.norm.sf(abs(a - 0.5) / se))


def score_name(s):
    s = SCORE.get(s, s)
    return s.replace("R_rob (V0 raw)", "R_rob (V0, raw)").replace("(V1 pre-specified)", "(V1)").replace("(V2 genes removed)", "(V2)").replace("(V3 resid. (all 3))", "(V3)").replace(
        "(V4 resid. (prolif.))", "(V4)").replace("(V5 stratified)", "(V5)").replace("(V6 spline resid.)", "(V6)").replace("(V7 resid. (prolif. + identity))", "(V7)")


SPEC = lambda s: s.replace("F8 without NIA donors aged 80+", "F8 without AG donors aged 80 or older").replace(" (age 18+; age + sex)", "").replace(" (GTEx range)", "")

# cohort agreement: permuted-age null per discovery specification (G0)
COH = "Cohort agreement (permuted-age null)"
for _, r in CSV("x_cohort_perm_p.csv").iterrows():
    add(COH, f"{SPEC(r.spec)} vs GTEx {r.tissue}, all genes", "Spearman ρ", r.rho_all, r.p_rho_all, where="Figure 1E; ⟨cohortsens⟩")
    add(COH, f"{SPEC(r.spec)} vs GTEx {r.tissue}, 200 strongest genes", "Spearman ρ", r.top_rho, r.p_top_rho,
        where=("Figures 1D and 1F; ⟨cohortsens⟩" if r.spec.startswith("F0") and r.tissue == "fibroblasts" else "Figure 1F; ⟨cohortsens⟩"))
# same permuted-age run (x_cohort_checks_v2.py): F9b and the sign agreement of signature genes
CK = J("x_cohort_checks_v2.json")
for t in ("fibroblasts", "skin"):
    v = CK["nulls"][t]["F9b_top_rho"]
    add(COH, f"F9b factors after a sex-only model vs GTEx {t}, 200 strongest genes", "Spearman ρ", v["obs"], v["p"], where="⟨addstats⟩")
for t in ("fibroblasts", "skin"):
    v = CK["nulls"][t]["sig_sign_agree"]
    add(COH, f"F0 pre-specified vs GTEx {t}, signature genes with the same age direction", "fraction", v["obs"], v["p"], where="Figure S2; ⟨addstats⟩")
# G2 model with age permuted within nucleic-acid isolation batch (x_g2_null.py)
G2N = J("x_g2_null.json")
for t in ("fibroblasts", "skin"):
    add(COH, f"F0 pre-specified vs GTEx {t}, 200 strongest genes, G2 model (age permuted within nucleic-acid isolation batch)", "Spearman ρ", G2N[t]["top_rho"], G2N[t]["p"], where="⟨addstats⟩ (1)")
# the 200 strongest GTEx-fibroblast genes against the discovery cohort (x_reverse_top200.py; Box 1, GTEx-derived signature)
RV = J("x_reverse_top200.json")
add(COH, "200 strongest GTEx-fibroblast genes vs discovery (F0)", "Spearman ρ", RV["rho"], RV["p"], where="⟨boxcheck⟩; ⟨addstats⟩ (1)")
# top-k sweep: the F0 top-200 test repeated at other k (k = 200 is listed above), same permuted-age run
for t in ("fibroblasts", "skin"):
    for k in CK["ks"]:
        if k == 200: continue
        v = CK["nulls"][t][f"k{k}"]
        add("Top-k sweep (F0 test at other k, permuted-age null)", f"F0 pre-specified vs GTEx {t}, {k:,} strongest genes", "Spearman ρ", v["rho"], v["p"], where="Figure S1E; ⟨addstats⟩")
# sex control
SEXN = {"sex_rho_all": "all genes", "sex_top_rho": "200 most sex-biased genes", "sex_top_autosomal_rho": "200 most sex-biased autosomal genes"}
for k, v in J("x_cohort_nulls.json")["sex_perm"].items():
    if isinstance(v, dict):
        base, tissue = k.rsplit("_", 1)
        add("Sex positive control (permuted-sex null)", f"Sex effects, {SEXN[base]}, discovery vs GTEx {tissue}", "Spearman ρ", v["observed"], v["p"], where="Figure 1D; ⟨sexctrl⟩")
# random signatures (real side without proliferation genes)
RSN = {"H1_median_rho_raw": "median cross-line ρ, raw score", "H1_median_rho_adj": "median cross-line ρ, R_adj", "H4_CellAge_auroc_raw": "CellAge AUROC, raw R_rob",
       "H4_CellAge_auroc_adj": "CellAge AUROC, V1 R_rob", "H4_GenAge_auroc_raw": "GenAge AUROC, raw R_rob", "H4_GenAge_auroc_adj": "GenAge AUROC, V1 R_rob"}
for _, r in CSV("x_random_signature_v2.csv").iterrows():
    nm = RSN.get(r.statistic) or f"ρ(score, proliferation), {SCR[r.statistic.split('_')[-1]]}"
    add("Real signature vs 100 random signatures (one-sided empirical P)", nm, r.statistic.split("_")[0].replace("H2", "ρ").replace("H1", "ρ").replace("H4", "AUROC"),
        r.real_V2_input, r.empirical_p_one_sided, 100, "Figures 2E and S6")
# top-5% retention under V1 (binomial test against 5%)
for _, r in pd.read_csv(RES.parent / "figures" / "source_data" / "Fig2F.csv").iterrows():
    kept = int(round(r.k * r.kept_pct / 100))
    add("Raw top 5% kept under V1 (binomial, against 5%)", SCR[r.screen], "fraction kept", r.kept_pct / 100, stats.binomtest(kept, int(r.k), 0.05, "greater").pvalue, int(r.k), "Figure 2F")
# network: 500 nulls
for c, v in J("rev_network_null.json")["screens"].items():
    add("Network vs 500 randomized networks, rank (one-sided empirical P)", SCR[c], "Spearman ρ", v["rho_network"], v["empirical_p_rho"], v["n_eval"], "Figure 4B")
    add("Network vs mean of 500 randomized networks, absolute error (one-sided Wilcoxon)", SCR[c], "MAE", v["mae_network"], v["h3a_wilcoxon_vs_mean_null_p"], v["n_eval"], "Figure 4A")
# held-out labels, discovery signature (pre-specified H4a/H4b rows excluded)
for _, r in CSV("rev3_scores_evidence.csv").iterrows():
    if r.score == "R_rob (V1, pre-specified)" and r.gene_set == "all genes" and r.labels in ("CellAge", "GenAge"): continue
    add("Held-out labels, discovery signature", f"{LAB.get(r.labels, r.labels)}, {r.gene_set.replace('>=', '≥')}, {score_name(r.score)}", "AUROC", r.auroc,
        p_auroc(r.auroc, r.n_pos, r.n_neg), f"{r.n_pos}/{r.n_neg}", "Figures 5A and 5B; ⟨evidence⟩")
# held-out labels, GTEx-derived signature (proliferation baseline already listed above)
for _, r in CSV("x_gtex_signature_evidence.csv").iterrows():
    if r.score == "proliferation shift alone": continue
    add("Held-out labels, GTEx-derived signature", f"{LAB.get(r.labels, r.labels)}, {r.gene_set.replace('>=', '≥')}, {score_name(r.score)}", "AUROC", r.auroc,
        p_auroc(r.auroc, r.n_pos, r.n_neg), f"{r.n_pos}/{r.n_neg}", "Figure S7B; ⟨evidence⟩")
# CellAge by manipulation: the method subsets only (the all-genes rows are listed in the families above or are H4a);
# proliferation listed once
seen = set()
for _, r in CSV("x_cellage_method.csv").iterrows():
    if r.subset == "all CellAge genes": continue
    if r.score == "proliferation shift alone":
        if r.subset in seen: continue
        seen.add(r.subset); sig = "same for both signatures"
    else: sig = r.signature.replace("GTEx-fibroblast signature", "GTEx-derived signature").replace("Discovery signature", "discovery signature")
    add("Held-out labels, CellAge by manipulation", f"{r.subset}, {score_name(r.score)} ({sig})", "AUROC", r.auroc, p_auroc(r.auroc, r.n_pos, r.n_neg), f"{r.n_pos}/{r.n_neg}", "⟨evidence⟩")
H4 = J("results.json")["H4"]
for key, lab in [("H4a_cellage", "CellAge"), ("H4b_genage", "GenAge")]:
    for c in SCR:
        for typ in ["raw", "adj"]:
            e = H4[key][f"{typ}_{c}"]
            add("Held-out labels, single screens", f"{lab}, {SCR[c]}, {'raw score' if typ == 'raw' else 'R_adj (V1)'}", "AUROC", e["auroc"],
                e.get("p", p_auroc(e["auroc"], e["n_pos"], e["n_neg"])), f"{e['n_pos']}/{e['n_neg']}", "Figure S9")
# paired AUROC and logistic models
for _, r in CSV("rev3_scores_paired_auroc.csv").iterrows():
    add("Paired AUROC differences (DeLong), discovery signature", f"{LAB.get(r.labels, r.labels)}, {r.gene_set.replace('>=', '≥')}: {score_name(r.a)} minus {score_name(r.b)}",
        "AUROC difference", r.diff_a_minus_b, r.p_two_sided, f"{r.n_pos}/{r.n_neg}", "Figure 5C; ⟨evidence⟩")
for _, r in CSV("x_gtex_signature_paired_auroc.csv").iterrows():
    add("Paired AUROC differences (DeLong), GTEx-derived signature", f"{LAB.get(r.labels, r.labels)}, {r.gene_set.replace('>=', '≥')}: {score_name(r.a)} minus {score_name(r.b)}",
        "AUROC difference", r.diff_a_minus_b, r.p_two_sided, f"{r.n_pos}/{r.n_neg}", "⟨evidence⟩")
for _, r in CSV("rev3_scores_logistic.csv").iterrows():
    add("Raw score added to proliferation (logistic, likelihood ratio)", f"{LAB.get(r.labels, r.labels)}, {r.gene_set.replace('>=', '≥')}", "coefficient", r.coef_raw, r.lr_p_raw_added, r.n, "⟨evidence⟩")
CLS = {"Induces": "inducers", "Inhibits": "inhibitors", "Anti": "anti-longevity", "Pro": "pro-longevity"}
REF = {"DepMap_common_essential": "DepMap common essentials", "CEGv2": "Hart CEGv2 core essentials"}
for _, r in CSV("rev3_depmap_label_fisher.csv").iterrows():
    add("Label classes and essential genes (Fisher)", f"{r.labels}, {REF.get(r.reference, r.reference)}: {CLS.get(r.positive_class, r.positive_class)} vs {CLS.get(r.negative_class, r.negative_class)}",
        "odds ratio", r.odds_ratio, r.fisher_p_two_sided, where="Figure 5D")
# hidden factors of the discovery cohort: Hallmark enrichment of the loadings
for _, r in CSV("x_cohort_hidden_factor_hallmark.csv").iterrows():
    add("Hidden-factor loadings, Hallmark enrichment (hypergeometric)", f"Factor {r.factor}, {r.direction} loadings: {hallmark_name(r.set)}", "overlap", r.overlap, r.p, f"{r.overlap}/{r['size']}", "⟨addstats⟩")
# GTEx fibroblasts vs skin
gp = RES / "x_gtex_pair_null.json"
if gp.exists():
    g = J("x_gtex_pair_null.json")
    add("GTEx fibroblasts vs GTEx skin (permuted-age null)", "age effects, all shared genes", "Spearman ρ", g["observed_rho"], g["p_two_sided"], g["n_genes"], "⟨addstats⟩ (1)")

T = pd.DataFrame(rows)
T["Tests in family"] = T.groupby("Family").Test.transform("size")
T["Holm P (within family)"] = T.groupby("Family").P.transform(lambda p: multipletests(p.fillna(1), method="holm")[1])
T["BH q (within family)"] = T.groupby("Family").P.transform(lambda p: multipletests(p.fillna(1), method="fdr_bh")[1])
T = T[["Family", "Tests in family", "Test", "Statistic", "Value", "n", "P", "Holm P (within family)", "BH q (within family)", "Shown in"]]
T.to_csv(RES / "x_exploratory_tests.csv", index=False)
print(len(T), "tests in", T.Family.nunique(), "families")
print(T.groupby("Family").agg(n=("P", "size"), p05=("P", lambda p: (p < 0.05).sum()), holm05=("Holm P (within family)", lambda p: (p < 0.05).sum()),
                              bh05=("BH q (within family)", lambda p: (p < 0.05).sum())).to_string())
