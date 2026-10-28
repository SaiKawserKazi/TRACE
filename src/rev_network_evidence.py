"""Exploratory (Amendment 1): network confidence intervals and evidence audit.

Network: bootstrap (over held-out knockdowns) 95% CIs for every predictor's MAE and Spearman rho, paired
differences against the transfer-only and constant baselines, relative error change, and a check that the
randomised networks preserve every node degree.
Evidence: label coverage and mapping table, AUROC with bootstrap CIs, label-permutation p values, and
label-subset sensitivity (CellAge non-cancer / replicative evidence; GenAge per organism).
Outputs: results/rev_network.csv, results/rev_evidence_coverage.csv, results/rev_evidence_sensitivity.csv,
results/rev_network_evidence.json
"""
import os, sys, json, io, zipfile, time
import numpy as np, pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
import trace_analysis as T

RES, RAW = T.RES, T.RAW
rng = np.random.default_rng(3)
NB = 2000


def log(*a): print(time.strftime("%H:%M:%S"), *a, flush=True)


def boot_ci(f, n, nb=NB):
    vals = [f(rng.integers(0, n, n)) for _ in range(nb)]
    return float(np.quantile(vals, 0.025)), float(np.quantile(vals, 0.975))


def network_part(out):
    rows = []
    for c in ["K562", "K562gw", "RPE1", "HepG2", "Jurkat"]:
        f = os.path.join(RES, f"h3_{c}.csv")
        if not os.path.exists(f): continue
        D = pd.read_csv(f, index_col=0)
        for task, cols in [("H3a (network available)", ["network", "null", "mean"]),
                           ("H3b (transfer and network available)", ["transfer", "transfer_cal", "combined", "mean", "network"])]:
            d = D.dropna(subset=cols)
            y = d.y.values; n = len(d)
            for col in cols:
                pred = d[col].values; err = np.abs(pred - y)
                r = {"screen": c, "task": task, "predictor": col, "n": n, "mae": err.mean(),
                     "rho": stats.spearmanr(pred, y)[0]}
                r["mae_lo"], r["mae_hi"] = boot_ci(lambda i: err[i].mean(), n)
                r["rho_lo"], r["rho_hi"] = boot_ci(lambda i: stats.spearmanr(pred[i], y[i])[0], n, 2000)
                ref = "null" if task.startswith("H3a") else "transfer_cal"
                if col != ref:
                    e0 = np.abs(d[ref].values - y); dd = err - e0
                    r[f"mae_diff_vs_{ref}"] = dd.mean(); r["rel_change_pct"] = 100 * dd.mean() / e0.mean()
                    r["diff_lo"], r["diff_hi"] = boot_ci(lambda i: dd[i].mean(), n)
                    r["wilcoxon_p_two_sided"] = stats.wilcoxon(err, e0).pvalue
                    r["frac_better"] = float((err < e0).mean())
                rows.append(r)
    N = pd.DataFrame(rows)
    N.to_csv(os.path.join(RES, "rev_network.csv"), index=False)
    out["network_table"] = "results/rev_network.csv"
    sel = N[(N.task.str.startswith("H3b")) & (N.predictor == "combined")]
    out["combined_vs_cal_rel_change_pct"] = dict(zip(sel.screen, sel.rel_change_pct.round(2)))
    out["combined_vs_cal_diff_ci"] = {s: [float(a), float(b)] for s, a, b in zip(sel.screen, sel.diff_lo, sel.diff_hi)}
    out["combined_rho_ci"] = {s: [float(a), float(b), float(c)] for s, a, b, c in zip(sel.screen, sel.rho, sel.rho_lo, sel.rho_hi)}
    sel = N[(N.task.str.startswith("H3b")) & (N.predictor == "transfer_cal")]
    out["transfer_cal_rho_ci"] = {s: [float(a), float(b), float(c)] for s, a, b, c in zip(sel.screen, sel.rho, sel.rho_lo, sel.rho_hi)}
    sel = N[(N.task.str.startswith("H3a")) & (N.predictor == "network")]
    out["network_vs_null_rel_change_pct"] = dict(zip(sel.screen, sel.rel_change_pct.round(2)))
    out["network_rho_ci"] = {s: [float(a), float(b), float(c)] for s, a, b, c in zip(sel.screen, sel.rho, sel.rho_lo, sel.rho_hi)}
    sel = N[(N.task.str.startswith("H3a")) & (N.predictor == "mean")]
    out["mean_vs_null_rel_change_pct"] = dict(zip(sel.screen, sel.rel_change_pct.round(2)))
    # degree preservation of the randomised networks (same construction and seeds as the main analysis)
    G = T.load_string()
    deg = dict(G.degree())
    chk = []
    for s in range(10):
        H = T.randomised(G, s)
        same = all(H.degree(n) == d for n, d in deg.items())
        E0 = set(map(frozenset, G.edges())); E1 = set(map(frozenset, H.edges()))
        chk.append({"seed": s, "degrees_identical": same, "edges": H.number_of_edges(),
                    "frac_edges_rewired": 1 - len(E0 & E1) / len(E0)})
    out["null_checks"] = chk
    log("network done", chk[0])


def evidence_part(out):
    rob = pd.read_csv(os.path.join(RES, "robust_scores.csv"), index_col=0)
    R = rob.R_rob
    screens = {c: pd.read_csv(os.path.join(RES, f"scores_{c}.csv"), index_col=0) for c in ["K562", "K562gw", "RPE1", "HepG2", "Jurkat"]}
    z = zipfile.ZipFile(os.path.join(RAW, "cellage.zip"))
    cdf = pd.read_csv(io.BytesIO(z.read([x for x in z.namelist() if x.endswith(".tsv")][0])), sep="\t")
    ca = T.cellage()
    gl, n_org = T.genage_human_labels()
    # GenAge mapping detail
    zg = zipfile.ZipFile(os.path.join(RAW, "genage_models.zip"))
    gdf = pd.read_csv(io.BytesIO(zg.read("genage_models.csv")))
    org = {"Mus musculus": "mmusculus", "Caenorhabditis elegans": "celegans", "Drosophila melanogaster": "dmelanogaster", "Saccharomyces cerevisiae": "scerevisiae"}
    cov_rows = []
    cov_rows.append({"source": "CellAge", "subset": "all records", "records": len(cdf), "genes": cdf["Gene symbol"].nunique(),
                     "induces": int((ca == "Induces").sum()), "inhibits": int((ca == "Inhibits").sum()),
                     "conflicting_or_unclear_dropped": int((~ca.isin(["Induces", "Inhibits"])).sum()),
                     "evaluable_pos": int(R.reindex(ca[ca == "Induces"].index).notna().sum()),
                     "evaluable_neg": int(R.reindex(ca[ca == "Inhibits"].index).notna().sum())})
    votes = []
    for o, key in org.items():
        sub = gdf[gdf.organism == o]
        f = os.path.join(RAW, f"ortho_{key}.tsv")
        orth = pd.read_csv(f, sep="\t").dropna(); orth.columns = ["human", "model", "type"]
        lab = sub[sub["longevity influence"].isin(["Pro-Longevity", "Anti-Longevity"])]
        mm = lab.merge(orth, left_on="symbol", right_on="model")
        votes.append(mm.assign(org=o))
        cov_rows.append({"source": "GenAge", "subset": o, "records": len(sub), "genes": sub.symbol.nunique(),
                         "induces": int((lab["longevity influence"] == "Anti-Longevity").sum()),
                         "inhibits": int((lab["longevity influence"] == "Pro-Longevity").sum()),
                         "mapped_model_genes": int(mm.symbol.nunique()), "human_orthologs": int(mm.human.nunique())})
    Vt = pd.concat(votes)
    tab = Vt.groupby(["human", "longevity influence"]).size().unstack(fill_value=0)
    ties = int((tab.get("Anti-Longevity", 0) == tab.get("Pro-Longevity", 0)).sum())
    conflicting = int(((tab.get("Anti-Longevity", 0) > 0) & (tab.get("Pro-Longevity", 0) > 0)).sum())
    cov_rows.append({"source": "GenAge", "subset": "human consensus", "genes": len(tab), "induces": int((gl == "Anti").sum()),
                     "inhibits": int((gl == "Pro").sum()), "conflicting_or_unclear_dropped": ties,
                     "conflicting_resolved_by_majority": conflicting - ties,
                     "evaluable_pos": int(R.reindex(gl[gl == "Anti"].index).notna().sum()),
                     "evaluable_neg": int(R.reindex(gl[gl == "Pro"].index).notna().sum())})
    for c, s in screens.items():
        cov_rows.append({"source": "coverage", "subset": c, "genes": len(s),
                         "evaluable_pos": int(s.index.isin(ca[ca == "Induces"].index).sum()),
                         "evaluable_neg": int(s.index.isin(ca[ca == "Inhibits"].index).sum()),
                         "mapped_model_genes": int(s.index.isin(gl[gl == "Anti"].index).sum()),
                         "human_orthologs": int(s.index.isin(gl[gl == "Pro"].index).sum())})
    pd.DataFrame(cov_rows).to_csv(os.path.join(RES, "rev_evidence_coverage.csv"), index=False)

    def test(score, pos, neg, label):
        x = score.reindex(pos).dropna().values; y = score.reindex(neg).dropna().values
        if len(x) < 5 or len(y) < 5: return None
        a = stats.mannwhitneyu(x, y, alternative="greater").statistic / (len(x) * len(y))
        bs = [stats.mannwhitneyu(rng.choice(x, len(x)), rng.choice(y, len(y))).statistic / (len(x) * len(y)) for _ in range(2000)]
        allv = np.r_[x, y]; perm = []
        for _ in range(10000):
            p_ = rng.permutation(allv); perm.append(stats.mannwhitneyu(p_[:len(x)], p_[len(x):]).statistic / (len(x) * len(y)))
        return {"analysis": label, "n_pos": len(x), "n_neg": len(y), "auroc": a, "lo": float(np.quantile(bs, 0.025)),
                "hi": float(np.quantile(bs, 0.975)), "perm_p_one_sided": float((np.sum(np.array(perm) >= a) + 1) / (len(perm) + 1))}
    rows = []
    ind = ca[ca == "Induces"].index; inh = ca[ca == "Inhibits"].index
    rows.append(test(R, ind, inh, "CellAge, all (pre-specified H4a)"))
    nc = cdf[cdf["Cancer Cell"] == "No"].groupby("Gene symbol")["Senescence Effect"].agg(lambda s: s.mode()[0] if s.nunique() == 1 else "Mixed")
    rows.append(test(R, nc[nc == "Induces"].index, nc[nc == "Inhibits"].index, "CellAge, non-cancer-cell evidence only"))
    rep = cdf[cdf["Type of senescence"].astype(str).str.contains("Replicative")].groupby("Gene symbol")["Senescence Effect"].agg(lambda s: s.mode()[0] if s.nunique() == 1 else "Mixed")
    rows.append(test(R, rep[rep == "Induces"].index, rep[rep == "Inhibits"].index, "CellAge, replicative senescence only"))
    rows.append(test(R, gl[gl == "Anti"].index, gl[gl == "Pro"].index, "GenAge, human consensus (pre-specified H4b)"))
    for o in org:
        sub = Vt[Vt.org == o].groupby("human")["longevity influence"].agg(lambda s: s.mode()[0] if s.nunique() == 1 else "Mixed")
        rows.append(test(R, sub[sub == "Anti-Longevity"].index, sub[sub == "Pro-Longevity"].index, f"GenAge, {o} only"))
    # proliferation shift alone as the "score" (confounding check): mean over screens, higher = more proliferative
    pro = pd.concat([s.prolif for s in screens.values()], axis=1).mean(axis=1)
    rows.append(test(pro, ind, inh, "CellAge, proliferation shift alone as score"))
    rows.append(test(pro, gl[gl == "Anti"].index, gl[gl == "Pro"].index, "GenAge, proliferation shift alone as score"))
    raw = pd.concat([s.R_raw for s in screens.values()], axis=1).mean(axis=1)
    rows.append(test(raw, ind, inh, "CellAge, mean raw score (no adjustment)"))
    rows.append(test(raw, gl[gl == "Anti"].index, gl[gl == "Pro"].index, "GenAge, mean raw score (no adjustment)"))
    S = pd.DataFrame([r for r in rows if r])
    S.to_csv(os.path.join(RES, "rev_evidence_sensitivity.csv"), index=False)
    out["evidence"] = S.round(4).to_dict("records")
    out["coverage"] = pd.DataFrame(cov_rows).to_dict("records")
    log("evidence done")


if __name__ == "__main__":
    out = {}
    evidence_part(out)
    network_part(out)
    json.dump(out, open(os.path.join(RES, "rev_network_evidence.json"), "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    print("done")


