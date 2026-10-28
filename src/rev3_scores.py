"""Exploratory (Amendment 3): score-level controls.

E5  random expression-matched signature null (100 signatures through H1, H2 and H4)
E7  paired AUROC tests (DeLong), incremental logistic regression, paired drops V0 -> V1..V7, minimum detectable AUROC
E8  CellAge label convention, release, overlap with core-essential genes
E9  GenAge mammal-only and no-yeast tests; overlap with core-essential genes
E10 Rrob for single-line genes; H4 on all genes and on genes measured in >= 2 lines
E13 proliferation deciles, distance correlation, mutual information; V6 (spline residual) and V7 (no magnitude term)
E14 Pareto front before and after filtering (on-target knockdown, response magnitude, cell count)
E15 ranks of known senescence genes under V0..V7
E16 up-minus-down score form; per-arm coverage
Outputs: results/rev3_scores_*.csv, results/rev3_scores.json
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import config as _cfg
import os, sys, json, time, io, zipfile
import multiprocessing as mp
import numpy as np, pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
import aging, perturb
import trace_analysis as T

RES, RAW = T.RES, T.RAW
CTX = ["K562", "K562gw", "RPE1", "HepG2", "Jurkat"]; LINES = ["K562", "RPE1", "HepG2", "Jurkat"]
NSIG = _cfg.N_RANDOM_SIGNATURES; NPROC = _cfg.N_PROC
VAR = ["V0 raw", "V1 pre-specified", "V2 genes removed", "V3 resid. (all 3)", "V4 resid. (prolif.)", "V5 stratified",
       "V6 spline resid.", "V7 resid. (prolif. + identity)"]
SEN = ["CDKN1A", "CDKN2A", "TP53", "RB1", "LMNB1"]


def log(*a): print(time.strftime("%H:%M:%S"), *a, flush=True)


# ------------------------------------------------------------------ scores
def spline_basis(x, df=5):
    from patsy import dmatrix
    return pd.DataFrame(np.asarray(dmatrix(f"bs(x, df={df}) - 1", {"x": np.asarray(x)})), index=x.index)


def variants(ctx, sig_t, progs):
    sig_rm = sig_t[[g for g in sig_t.index if g not in T.PROLIF_REMOVE]]
    out = {v: {} for v in VAR}
    for c, (shift, nt, gm, nc) in ctx.items():
        P = progs[c]; raw = T.score(shift, sig_t); rm = T.score(shift, sig_rm)
        out["V0 raw"][c] = raw
        out["V1 pre-specified"][c] = T.residualise(rm, P[["prolif", "identity", "magnitude"]])
        out["V2 genes removed"][c] = rm
        out["V3 resid. (all 3)"][c] = T.residualise(raw, P[["prolif", "identity", "magnitude"]])
        out["V4 resid. (prolif.)"][c] = T.residualise(raw, P[["prolif"]])
        q = pd.qcut(P.prolif.rank(method="first"), 5, labels=False)
        out["V5 stratified"][c] = raw.groupby(q.values).rank(pct=True) - 0.5
        out["V6 spline resid."][c] = T.residualise(rm, pd.concat([spline_basis(P.prolif), P[["identity", "magnitude"]]], axis=1))
        out["V7 resid. (prolif. + identity)"][c] = T.residualise(rm, P[["prolif", "identity"]])
    return out


def robust_of(vd): return T.robust({c: pd.DataFrame({"R_adj": s}) for c, s in vd.items()})


# ------------------------------------------------------------------ AUROC machinery
def auc(x, y): return stats.mannwhitneyu(x, y).statistic / (len(x) * len(y))


def delong(pos_a, neg_a, pos_b, neg_b):
    """Paired DeLong test for two scores on the same positives/negatives. Returns (aucA, aucB, diff, se, lo, hi, p two-sided)."""
    def comps(pos, neg):
        m, n = len(pos), len(neg)
        allv = np.r_[pos, neg]; r = stats.rankdata(allv)
        rp = stats.rankdata(pos); rn = stats.rankdata(neg)
        a = (r[:m].sum() - m * (m + 1) / 2) / (m * n)
        v10 = (r[:m] - rp) / n          # structural components for positives
        v01 = 1 - (r[m:] - rn) / m      # for negatives
        return a, v10, v01
    aA, v10A, v01A = comps(pos_a, neg_a); aB, v10B, v01B = comps(pos_b, neg_b)
    s10 = np.cov(np.vstack([v10A, v10B])); s01 = np.cov(np.vstack([v01A, v01B]))
    m, n = len(pos_a), len(neg_a)
    var = (s10[0, 0] + s10[1, 1] - 2 * s10[0, 1]) / m + (s01[0, 0] + s01[1, 1] - 2 * s01[0, 1]) / n
    se = np.sqrt(max(var, 1e-12)); d = aA - aB
    return aA, aB, d, se, d - 1.96 * se, d + 1.96 * se, 2 * stats.norm.sf(abs(d) / se)


def min_detectable_auc(n1, n2, alpha=0.05, power=0.8):
    z = stats.norm.ppf(1 - alpha) + stats.norm.ppf(power)
    for A in np.arange(0.5005, 0.9, 0.0005):
        q1, q2 = A / (2 - A), 2 * A * A / (1 + A)
        se = np.sqrt((A * (1 - A) + (n1 - 1) * (q1 - A * A) + (n2 - 1) * (q2 - A * A)) / (n1 * n2))
        if (A - 0.5) / se >= z: return float(round(A, 3))
    return np.nan


def labels():
    ca = T.cellage(); gl, _ = T.genage_human_labels()
    z = zipfile.ZipFile(os.path.join(RAW, "genage_models.zip")); gdf = pd.read_csv(io.BytesIO(z.read("genage_models.csv")))
    org = {"Mus musculus": "mmusculus", "Caenorhabditis elegans": "celegans", "Drosophila melanogaster": "dmelanogaster", "Saccharomyces cerevisiae": "scerevisiae"}
    votes = []
    for o, key in org.items():
        orth = pd.read_csv(os.path.join(RAW, f"ortho_{key}.tsv"), sep="\t").dropna(); orth.columns = ["human", "model", "type"]
        sub = gdf[(gdf.organism == o) & gdf["longevity influence"].isin(["Pro-Longevity", "Anti-Longevity"])]
        votes.append(sub.merge(orth, left_on="symbol", right_on="model")[["human", "longevity influence", "organism"]])
    V = pd.concat(votes)

    def consensus(df):
        tab = df.groupby(["human", "longevity influence"]).size().unstack(fill_value=0)
        for c in ["Pro-Longevity", "Anti-Longevity"]:
            if c not in tab: tab[c] = 0
        lab = pd.Series(np.where(tab["Anti-Longevity"] > tab["Pro-Longevity"], "Anti", np.where(tab["Pro-Longevity"] > tab["Anti-Longevity"], "Pro", "Tie")), index=tab.index)
        return lab[lab != "Tie"]
    return {"CellAge": (ca, "Induces", "Inhibits"), "GenAge": (gl, "Anti", "Pro"),
            "GenAge mammal-only": (consensus(V[V.organism == "Mus musculus"]), "Anti", "Pro"),
            "GenAge without yeast": (consensus(V[V.organism != "Saccharomyces cerevisiae"]), "Anti", "Pro")}


def evid(score, lab, pos, neg):
    x = score.reindex(lab[lab == pos].index).dropna(); y = score.reindex(lab[lab == neg].index).dropna()
    return x, y


def dcor(x, y, n=1500, reps=5, seed=0):
    r = np.random.default_rng(seed); vals = []
    x, y = np.asarray(x), np.asarray(y)
    for _ in range(reps):
        i = r.choice(len(x), min(n, len(x)), replace=False); a, b = x[i], y[i]
        A = np.abs(a[:, None] - a[None, :]); B = np.abs(b[:, None] - b[None, :])
        A = A - A.mean(0) - A.mean(1)[:, None] + A.mean(); B = B - B.mean(0) - B.mean(1)[:, None] + B.mean()
        vals.append(np.sqrt(max((A * B).mean(), 0) / np.sqrt((A * A).mean() * (B * B).mean())))
    return float(np.mean(vals))


# ------------------------------------------------------------------ X11 random signatures (parallel)
_W = {}


def _init():
    ctx = {c: perturb.load_context(c) for c in CTX}
    ident = T.identity_genes(ctx)
    _W.update(ctx=ctx, progs={c: T.program_shifts(v[0], ident[c]) for c, v in ctx.items()}, labs=labels())
    _W["sig"] = pd.read_pickle(os.path.join(RES, "rev3_E5_inputs.pkl"))


def summarize(per_raw, per_adj, progs, labs):
    out = {}
    pairs = [(a, b) for i, a in enumerate(LINES) for b in LINES[i + 1:]]
    out["H1_median_rho_raw"] = float(np.median([stats.spearmanr(*[per_raw[x].reindex(per_raw[a].index.intersection(per_raw[b].index)) for x in (a, b)])[0] for a, b in pairs]))
    out["H1_median_rho_adj"] = float(np.median([stats.spearmanr(*[per_adj[x].reindex(per_adj[a].index.intersection(per_adj[b].index)) for x in (a, b)])[0] for a, b in pairs]))
    for c in CTX: out[f"H2_rho_{c}"] = float(stats.spearmanr(per_raw[c], progs[c].prolif)[0])
    rraw, radj = robust_of(per_raw).R_rob, robust_of(per_adj).R_rob
    for nm in ["CellAge", "GenAge"]:
        lab, pos, neg = labs[nm]
        out[f"H4_{nm}_auroc_raw"] = auc(*evid(rraw, lab, pos, neg)); out[f"H4_{nm}_auroc_adj"] = auc(*evid(radj, lab, pos, neg))
    return out


def _sig_task(i):
    d = _W["sig"]; r = np.random.default_rng(60_000 + i)
    real, pool_by_bin, gbin = d["real"], d["pool"], d["bin"]
    genes, n_exact = [], 0
    for g in real.index:
        b = gbin[g]; bm, bv = divmod(b, 10)
        for w in range(10):  # widen to neighboring mean/variance deciles only if the exact bin is exhausted
            bins = [i * 10 + j for i in range(max(0, bm - w), min(9, bm + w) + 1) for j in range(max(0, bv - w), min(9, bv + w) + 1)]
            cand = [x for bb in bins for x in pool_by_bin.get(bb, []) if x not in genes]
            if cand: break
        n_exact += (w == 0)
        genes.append(cand[r.integers(len(cand))])
    sig = pd.Series(real.values, index=genes)
    sig_rm = sig  # random genes are not proliferation genes by construction of the pool
    per_raw, per_adj = {}, {}
    for c, (shift, nt, gm, nc) in _W["ctx"].items():
        P = _W["progs"][c]; per_raw[c] = T.score(shift, sig)
        per_adj[c] = T.residualise(T.score(shift, sig_rm), P[["prolif", "identity", "magnitude"]])
    return {"sig": i, "frac_exact_bin": n_exact / len(real), **summarize(per_raw, per_adj, _W["progs"], _W["labs"])}


def main():
    e, m = aging.load_fleischer(); ts = aging.age_tstats(e, m); sig_f = aging.signature(ts)
    ctx = {c: perturb.load_context(c) for c in CTX}
    ident = T.identity_genes(ctx); progs = {c: T.program_shifts(v[0], ident[c]) for c, v in ctx.items()}
    labs = labels(); out = {}
    Vd = variants(ctx, sig_f.t, progs)
    ROB = {v: robust_of(Vd[v]) for v in VAR}
    rob1 = ROB["V1 pre-specified"]
    # ---------- X14 Rrob definition and gene sets
    n_lines = rob1.n_lines
    out["E10"] = {"n_genes_total": int(len(rob1)), "n_single_line": int((n_lines == 1).sum()), "n_multi_line": int((n_lines >= 2).sum()),
                  "single_line_equals_radj": bool(np.allclose(rob1[n_lines == 1].R_rob, rob1[n_lines == 1].filter(like="Radj_").max(axis=1, skipna=True)))}
    single = rob1[n_lines == 1].filter(like="Radj_").notna().idxmax(axis=1).str.replace("Radj_", "").value_counts().to_dict()
    out["E10"]["single_line_by_line"] = single
    # ---------- H4 / evidence on all and multi-line genes (+ X13 paired tests)
    raw_mean = pd.concat([Vd["V0 raw"][c] for c in CTX], axis=1).mean(axis=1)
    pro_mean = pd.concat([progs[c].prolif for c in CTX], axis=1).mean(axis=1)
    rows, paired, logit = [], [], []
    multi = rob1.index[n_lines >= 2]
    for nm, (lab, pos, neg) in labs.items():
        for gs_name, gs in [("all genes", rob1.index), ("genes in >= 2 lines", multi)]:
            for sname, s in [("R_rob (V1, pre-specified)", rob1.R_rob), ("mean raw score", raw_mean), ("proliferation shift alone", pro_mean)] + \
                            [(f"R_rob ({v})", ROB[v].R_rob) for v in VAR if v not in ("V1 pre-specified",)]:
                x, y = evid(s.reindex(gs), lab, pos, neg)
                if len(x) >= 5 and len(y) >= 5:
                    rows.append({"labels": nm, "gene_set": gs_name, "score": sname, "n_pos": len(x), "n_neg": len(y), "auroc": auc(x, y),
                                 "min_detectable_auroc_80pct": min_detectable_auc(len(x), len(y))})
            # X13 paired DeLong: raw vs proliferation alone; V0 vs each other definition
            common = rob1.index.intersection(gs)
            for a_name, a_s, b_name, b_s in [("mean raw score", raw_mean, "proliferation shift alone", pro_mean)] + \
                                             [("R_rob (V0 raw)", ROB["V0 raw"].R_rob, f"R_rob ({v})", ROB[v].R_rob) for v in VAR[1:]]:
                both = a_s.reindex(common).dropna().index.intersection(b_s.reindex(common).dropna().index)
                xa, ya = evid(a_s.reindex(both), lab, pos, neg); xb, yb = evid(b_s.reindex(both), lab, pos, neg)
                if len(xa) >= 5 and len(ya) >= 5:
                    aA, aB, d, se, lo, hi, p = delong(xa.values, ya.values, xb.values, yb.values)
                    paired.append({"labels": nm, "gene_set": gs_name, "a": a_name, "b": b_name, "n_pos": len(xa), "n_neg": len(ya),
                                   "auroc_a": aA, "auroc_b": aB, "diff_a_minus_b": d, "se": se, "lo": lo, "hi": hi, "p_two_sided": p})
            # incremental logistic regression: label ~ proliferation + raw score (z-scored)
            import statsmodels.api as sm
            dfl = pd.DataFrame({"pro": pro_mean, "raw": raw_mean}).reindex(common).dropna()
            yv = lab.reindex(dfl.index); keep = yv.isin([pos, neg]); dfl, yv = dfl[keep], (yv[keep] == pos).astype(int)
            if yv.sum() >= 5:
                Xz = sm.add_constant((dfl - dfl.mean()) / dfl.std())
                fit_full = sm.Logit(yv, Xz).fit(disp=0); fit_pro = sm.Logit(yv, Xz[["const", "pro"]]).fit(disp=0)
                lr = 2 * (fit_full.llf - fit_pro.llf)
                logit.append({"labels": nm, "gene_set": gs_name, "n": int(len(yv)), "coef_raw": fit_full.params["raw"],
                              "coef_raw_lo": fit_full.conf_int().loc["raw", 0], "coef_raw_hi": fit_full.conf_int().loc["raw", 1],
                              "coef_prolif": fit_full.params["pro"], "lr_p_raw_added": float(stats.chi2.sf(lr, 1))})
    pd.DataFrame(rows).to_csv(os.path.join(RES, "rev3_scores_evidence.csv"), index=False)
    pd.DataFrame(paired).to_csv(os.path.join(RES, "rev3_scores_paired_auroc.csv"), index=False)
    pd.DataFrame(logit).to_csv(os.path.join(RES, "rev3_scores_logistic.csv"), index=False)
    log("evidence done")
    # ---------- X13 / X13 overlaps with core-essential genes (Hart CEGv2) and the essential-gene library
    ceg = set(pd.read_csv(os.path.join(RAW, "hart_CEGv2.txt"), sep="\t").iloc[:, 0])
    esslib = set(ctx["K562"][0].index)
    ov = []
    for nm, (lab, pos, neg) in labs.items():
        for cls in (pos, neg):
            g = set(lab[lab == cls].index)
            ov.append({"labels": nm, "class": cls, "n": len(g), "in_CEGv2": len(g & ceg), "frac_CEGv2": len(g & ceg) / max(1, len(g)),
                       "in_essential_library": len(g & esslib), "frac_essential_library": len(g & esslib) / max(1, len(g))})
    pd.DataFrame(ov).to_csv(os.path.join(RES, "rev3_scores_label_essential_overlap.csv"), index=False)
    out["CEGv2_size"] = len(ceg)
    # ---------- X16 deciles, distance correlation, mutual information
    from sklearn.feature_selection import mutual_info_regression
    dec, dep = [], []
    for c in CTX:
        P = progs[c]; qd = pd.qcut(P.prolif.rank(method="first"), 10, labels=False)
        for v in VAR:
            s = Vd[v][c].reindex(P.index)
            g = s.groupby(qd.values).agg(["mean", "std", "count"])
            for d_, r in g.iterrows(): dec.append({"screen": c, "variant": v, "decile": int(d_) + 1, "mean": r["mean"], "sd": r["std"], "n": int(r["count"])})
            ok = s.notna()
            dep.append({"screen": c, "variant": v, "spearman": float(stats.spearmanr(s[ok], P.prolif[ok])[0]), "distance_corr": dcor(s[ok], P.prolif[ok]),
                        "mutual_info": float(mutual_info_regression(P.prolif[ok].values.reshape(-1, 1), s[ok].values, random_state=0)[0])})
    pd.DataFrame(dec).to_csv(os.path.join(RES, "rev3_scores_prolif_deciles.csv"), index=False)
    pd.DataFrame(dep).to_csv(os.path.join(RES, "rev3_scores_prolif_dependence.csv"), index=False)
    log("E13 done")
    # ---------- top-5% movement for V6, V7 and Pareto fronts
    M = pd.DataFrame({v: ROB[v].R_rob for v in VAR}); Mm = M.loc[multi]
    k = max(1, int(0.05 * len(Mm)))
    top = {v: set(Mm[v].nlargest(k).index) for v in VAR}
    out["top5_overlap_with_V1_multi"] = {v: len(top[v] & top["V1 pre-specified"]) / k for v in VAR}
    out["k_top5_multi"] = k
    # ---------- X16 Pareto before / after filtering, for V1, V6, V7
    P0 = pd.read_csv(os.path.join(RES, "prioritisation.csv"), index_col=0)
    ncells = {}
    for c in CTX:
        nc = ctx[c][3]
        if nc is None:  # HepG2/Jurkat pseudobulk: read cell counts written by pseudobulk_sc
            import anndata as ad
            a = ad.read_h5ad(os.path.join(RAW, perturb.CONTEXTS[c]), backed="r")
            obs = a.obs; nc = pd.Series(obs.num_cells_filtered.values, index=[i.split("_")[1] for i in obs.index]).groupby(level=0).sum()
        ncells[c] = nc
    cells_mean = pd.concat([ncells[c].rename(c) for c in CTX], axis=1).mean(axis=1)
    nt_mag95 = {c: float(ctx[c][1].abs().mean(axis=1).quantile(0.95)) for c in CTX}
    out["nt_magnitude_95th"] = nt_mag95
    mag_ok = pd.concat([(progs[c].magnitude > nt_mag95[c]).rename(c) for c in CTX], axis=1).mean(axis=1) >= 0.5
    filt = (P0.ontarget_z < -0.5) & mag_ok.reindex(P0.index).fillna(False) & (cells_mean.reindex(P0.index) >= 25)
    out["E14_filter"] = {"n_before": int(len(P0)), "n_after": int(filt.sum()), "rule": "on-target z < -0.5; magnitude above the screen's non-targeting 95th percentile in at least half of the measured screens; >= 25 cells (mean over screens)"}
    objs = ["R_rob", "sign_consistency", "obj_identity", "obj_prolif", "obj_stress"]
    fronts = {}
    for lbl, sub in [("unfiltered", P0), ("filtered", P0[filt])]:
        f = T.pareto(sub[sub.R_rob > 0], objs); fronts[lbl] = list(f)
    out["E14_front"] = {k_: {"n": len(v_), "genes": v_} for k_, v_ in fronts.items()}
    P0.assign(passes_filter=filt, pareto_filtered=P0.index.isin(fronts["filtered"]), cells_mean=cells_mean.reindex(P0.index)).to_csv(os.path.join(RES, "rev3_scores_pareto_filtered.csv"))
    # ---------- X16 known senescence genes
    rk = []
    for g in SEN:
        row = {"gene": g, "measured_lines": int(rob1.n_lines.get(g, 0)) if g in rob1.index else 0}
        for v in VAR:
            s = ROB[v].R_rob.dropna()
            row[v] = int(s.rank(ascending=False)[g]) if g in s.index else None
        row["n_ranked"] = int(ROB["V1 pre-specified"].R_rob.notna().sum())
        rk.append(row)
    pd.DataFrame(rk).to_csv(os.path.join(RES, "rev3_scores_senescence_gene_ranks.csv"), index=False)
    # ---------- X16 up-minus-down score and arm coverage
    ud = []
    up = sig_f[sig_f.t > 0].index; dn = sig_f[sig_f.t < 0].index
    for c in CTX:
        sh = ctx[c][0]; u = up.intersection(sh.columns); d_ = dn.intersection(sh.columns)
        s_ud = -(sh[u].mean(axis=1) - sh[d_].mean(axis=1))
        ud.append({"screen": c, "n_up_measured": len(u), "n_down_measured": len(d_), "rho_with_R": float(stats.spearmanr(s_ud, Vd["V0 raw"][c])[0]),
                   "rho_ud_with_prolif": float(stats.spearmanr(s_ud, progs[c].prolif)[0])})
    pd.DataFrame(ud).to_csv(os.path.join(RES, "rev3_scores_updown.csv"), index=False)
    log("E14-E16 done")
    # ---------- real values of the X11 statistics, and X11 null
    out["E5_real"] = summarize(Vd["V0 raw"], Vd["V1 pre-specified"], progs, labs)
    gmean = e.loc[e.mean(axis=1) >= 1].mean(axis=1); gvar = e.loc[gmean.index].var(axis=1)
    measured = set.intersection(*[set(ctx[c][0].columns) for c in LINES])
    cand = gmean.index[gmean.index.isin(measured) & ~gmean.index.isin(T.PROLIF_REMOVE) & ~gmean.index.isin(sig_f.index)]
    qm = pd.qcut(gmean.rank(method="first"), 10, labels=False); qv = pd.qcut(gvar.rank(method="first"), 10, labels=False)
    gbin = (qm * 10 + qv)
    pool = {b: list(cand[gbin.reindex(cand).values == b]) for b in range(100)}
    real = sig_f.t[[g for g in sig_f.index if g not in T.PROLIF_REMOVE and g in gbin.index]]
    pd.to_pickle({"real": real, "pool": pool, "bin": gbin.to_dict()}, os.path.join(RES, "rev3_E5_inputs.pkl"))
    out["E5_real_noprolif_genes"] = int(len(real))
    # run in-process, reusing the loaded screens (parallel workers would each hold a full copy; 8 GB RAM)
    _W.update(ctx=ctx, progs=progs, labs=labs, sig=pd.read_pickle(os.path.join(RES, "rev3_E5_inputs.pkl")))
    nul = []
    for i in range(NSIG):
        nul.append(_sig_task(i))
        if i % 10 == 0: log("random signature", i)
    N = pd.DataFrame(nul); N.to_csv(os.path.join(RES, "rev3_scores_random_signatures.csv"), index=False)
    out["E5_null"] = {c: {"median": float(N[c].median()), "q025": float(N[c].quantile(0.025)), "q975": float(N[c].quantile(0.975)),
                          "frac_null_ge_real": float((N[c] >= out["E5_real"][c]).mean())} for c in out["E5_real"]}
    out["E5_frac_exact_bin_mean"] = float(N.frac_exact_bin.mean())
    json.dump(out, open(os.path.join(RES, "rev3_scores.json"), "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    log("done")


if __name__ == "__main__":
    main()
