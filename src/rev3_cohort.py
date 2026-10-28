"""Exploratory (Amendment 3): cohort analyses.

Analyses:
  * bug fix: the GTEx-skin G0 model now contains age + sex only (Amendments 1 and 2 also included ischemic time in G0; deviation D8)
  * E1  sex as a positive control for cohort transfer (all genes, autosomes only, X/Y only)
  * E2  permuted-age null (age shuffled within cohort, covariates fixed)
  * E3  sign check on a pre-named list of age genes
  * E4  hidden-factor adjustment (residual principal components), discovery without the 31 NIA donors aged 80+,
        power-matched random 74-donor subsets, VIF of age in the repository-adjusted model, age range per repository
  * E15 cohort transfer without Hallmark cell-cycle genes; age t of E2F/G2-M genes per cohort
  * E17 cross-cohort age prediction (elastic net)
Outputs: results/rev3_cohort_*.csv and results/rev3_cohort.json
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import config as _cfg
import os, sys, json, time, gzip
import multiprocessing as mp
import numpy as np, pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
import aging
from fleischer_meta import load_meta

ROOT = str(_cfg.ROOT); RAW, RES = os.path.join(ROOT, "data", "raw"), os.path.join(ROOT, "results")
NB = _cfg.N_COHORT_BOOT; NPERM = _cfg.N_COHORT_PERM
NSUB = _cfg.N_POWER_SUBSETS; NPROC = _cfg.N_PROC; KPC = 5
KNOWN = ["LMNB1", "CDKN2A", "CDKN1A", "COL1A1", "COL1A2", "COL3A1", "MKI67", "TOP2A", "IL6", "MMP1"]
SEXGENES = ["XIST", "RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "USP9Y", "ZFY", "KDM5C", "KDM6A"]
TOPK = 200


def log(*a): print(time.strftime("%H:%M:%S"), *a, flush=True)


def hallmark():
    sets = {}
    for l in open(os.path.join(RAW, "hallmark.txt"), encoding="utf-8"):
        f = l.rstrip("\n").split("\t"); sets[f[0]] = set(g for g in f[2:] if g)
    return sets


HM = hallmark()
CELLCYCLE = HM["E2F Targets"] | HM["G2-M Checkpoint"] | HM["Mitotic Spindle"]


def chrom_map():
    df = pd.read_csv(os.path.join(RAW, "GSE113957_fpkm.txt.gz"), sep="\t", usecols=["chr", "Annotation/Divergence"])
    df["gene"] = df["Annotation/Divergence"].str.split("|").str[0]
    return df.groupby("gene")["chr"].first().str.replace("chr", "")


# ------------------------------------------------------------------ cohorts
def fleischer():
    e, m = aging.load_fleischer()
    m = m.join(load_meta().reindex(m.index))
    m["site2"] = np.where(m.site.str.lower().str.contains("arm"), "arm", "other")
    return e[e.mean(axis=1) >= 0.5].copy(), m


def gtex(fname):
    e, m = aging.load_gtex(fname)
    s = pd.read_csv(os.path.join(RAW, "gtex_sample.txt"), sep="\t", low_memory=False, index_col=0).reindex(m.index)
    subj = pd.read_csv(os.path.join(RAW, "gtex_subject.txt"), sep="\t", index_col=0)
    m = m.assign(center=s.SMCENTER.values, batch_type=s.SMGEBTCHT.values, na_batch=s.SMNABTCH.values, rin=s.SMRIN.values,
                 ischemic=s.SMTSISCH.values, hardy=subj.reindex(m.donor).DTHHRDY.fillna(-1).astype(int).astype(str).values)
    return e[e.mean(axis=1) >= 0.5].copy(), m


FL = {  # name -> (mask, covariates, min_mean, hidden factors)
    "F0 pre-specified (age 18+; age + sex)": (lambda m: m.age >= 18, ["sex"], 1.0, 0),
    "F1 ages 20-79 (GTEx range)": (lambda m: (m.age >= 20) & (m.age < 80), ["sex"], 1.0, 0),
    "F2 + platform, repository, site": (lambda m: m.age >= 18, ["sex", "platform", "repository", "site2"], 1.0, 0),
    "F3 ages 20-79 + platform, repository, site": (lambda m: (m.age >= 20) & (m.age < 80), ["sex", "platform", "repository", "site2"], 1.0, 0),
    "F4 NextSeq only": (lambda m: (m.age >= 18) & (m.platform == "NextSeq 500"), ["sex"], 1.0, 0),
    "F5 arm biopsies only": (lambda m: (m.age >= 18) & (m.site2 == "arm"), ["sex"], 1.0, 0),
    "F6 expression filter 0.5": (lambda m: m.age >= 18, ["sex"], 0.5, 0),
    "F7 expression filter 2": (lambda m: m.age >= 18, ["sex"], 2.0, 0),
    "F8 without NIA donors aged 80+": (lambda m: (m.age >= 18) & ~((m.repository == "AG") & (m.age >= 80)), ["sex"], 1.0, 0),
    "F9 + 5 hidden factors": (lambda m: m.age >= 18, ["sex"], 1.0, KPC),
}
GT = {  # corrected: G0 is age + sex for both tissues; ischemic time enters skin models from G1 on
    "G0 pre-specified (age + sex)": (["sex"], [], 0),
    "G1 + center, batch type, RIN, Hardy": (["sex", "center", "batch_type", "rin", "hardy"], ["ischemic"], 0),
    "G2 G1 + nucleic-acid batch": (["sex", "center", "batch_type", "rin", "hardy", "na_batch"], ["ischemic"], 0),
    "G3 age + sex + 5 hidden factors": (["sex"], [], KPC),
}


def design(meta, covars, age_col="age"):
    cols, names = [np.ones(len(meta)), meta[age_col].values.astype(float)], ["const", "age"]
    for c in covars:
        v = meta[c]
        if v.dtype.kind in "biuf":
            cols.append(v.fillna(v.median()).values.astype(float)); names.append(c)
        else:
            d = pd.get_dummies(v.astype(str), drop_first=True)
            for k in d.columns:
                x = d[k].values.astype(float)
                if x.sum() > 0: cols.append(x); names.append(f"{c}:{k}")
    X = np.column_stack(cols)
    q, r = np.linalg.qr(X); keep = np.abs(np.diag(r)) > 1e-8 * np.abs(np.diag(r)).max(); keep[:2] = True
    return X[:, keep], [n for n, k in zip(names, keep) if k]


def fit(e, m, covars, min_mean=1.0, k_hidden=0, idx=None, age_col="age"):
    """OLS per gene; returns t for age and for sex (male indicator) plus BH q for age."""
    if idx is not None:
        e = e.iloc[:, idx]; m = m.iloc[idx]
    e = e[e.mean(axis=1) >= min_mean]
    X, names = design(m, covars, age_col)
    Y = e.values.T
    if k_hidden:
        Yz = (Y - Y.mean(0)) / (Y.std(0) + 1e-9)
        B0 = np.linalg.lstsq(X, Yz, rcond=None)[0]
        U, s, Vt = np.linalg.svd(Yz - X @ B0, full_matrices=False)
        X = np.column_stack([X, U[:, :k_hidden]]); names += [f"hf{i}" for i in range(k_hidden)]
    XtXi = np.linalg.pinv(X.T @ X); B = XtXi @ X.T @ Y
    dof = X.shape[0] - X.shape[1]; s2 = ((Y - X @ B) ** 2).sum(0) / dof
    out = pd.DataFrame(index=e.index)
    for nm in ["age"] + [n for n in names if n.startswith("sex:")]:
        j = names.index(nm); se = np.sqrt(s2 * XtXi[j, j])
        out["t" if nm == "age" else "t_sex"] = B[j] / np.where(se > 0, se, np.nan)
    out["p"] = 2 * stats.t.sf(np.abs(out.t), dof)
    out = out.dropna(subset=["t"]); out["q"] = stats.false_discovery_control(out.p.values)
    return out


def sig_from(ts):
    s = aging.signature(ts)
    return s if (s.t > 0).sum() >= 50 and (s.t < 0).sum() >= 50 else aging.signature_topk(ts)


def compare(tf, tg, sig_point, chrom):
    sh = tf.index.intersection(tg.index)
    nocc = sh.difference(pd.Index(sorted(CELLCYCLE)))
    sig_b = sig_from(tf); sb = sig_b.index.intersection(tg.index); sp = sig_point.index.intersection(tg.index)
    out = {"rho_all": stats.spearmanr(tf.t[sh], tg.t[sh])[0],
           "rho_no_cellcycle": stats.spearmanr(tf.t[nocc], tg.t[nocc])[0],
           "sig_sign_agree": float((np.sign(sig_b.t[sb]) == np.sign(tg.t[sb])).mean()),
           "sig_sign_agree_fixed": float((np.sign(sig_point.t[sp]) == np.sign(tg.t[sp])).mean())}
    top = tf.t[sh].abs().nlargest(TOPK).index
    out["top_sign_agree"] = float((np.sign(tf.t[top]) == np.sign(tg.t[top])).mean())
    out["top_rho"] = stats.spearmanr(tf.t[top], tg.t[top])[0]
    if "t_sex" in tf and "t_sex" in tg:
        c = chrom.reindex(sh)
        ta = tf.t_sex[sh[~c.isin(["X", "Y"]).values & c.notna().values]].abs().nlargest(TOPK).index
        tall = tf.t_sex[sh].abs().nlargest(TOPK).index
        out["sex_top_sign_agree"] = float((np.sign(tf.t_sex[tall]) == np.sign(tg.t_sex[tall])).mean())
        out["sex_top_rho"] = stats.spearmanr(tf.t_sex[tall], tg.t_sex[tall])[0]
        out["sex_top_autosomal_sign_agree"] = float((np.sign(tf.t_sex[ta]) == np.sign(tg.t_sex[ta])).mean())
        out["sex_top_autosomal_rho"] = stats.spearmanr(tf.t_sex[ta], tg.t_sex[ta])[0]
        auto = sh[~c.isin(["X", "Y"]).values & c.notna().values]; xy = sh[c.isin(["X", "Y"]).values]
        out.update(sex_rho_all=stats.spearmanr(tf.t_sex[sh], tg.t_sex[sh])[0],
                   sex_rho_autosomes=stats.spearmanr(tf.t_sex[auto], tg.t_sex[auto])[0],
                   sex_rho_xy=stats.spearmanr(tf.t_sex[xy], tg.t_sex[xy])[0] if len(xy) > 5 else np.nan, n_xy=len(xy))
    return out


# ------------------------------------------------------------------ parallel workers
_W = {}


def _init():
    _W["f"] = fleischer(); _W["g"] = {"fibroblasts": gtex("gtex_fibroblasts_tpm.gct.gz"), "skin": gtex("gtex_skin_ns_tpm.gct.gz")}
    _W["chrom"] = chrom_map(); _W["sig_point"] = pd.read_pickle(os.path.join(RES, "rev3_point_signatures.pkl"))


def _fit_all(rng_, perm=False):
    ef, mf = _W["f"]
    F = {}
    for name, (mask, cov, mm, kh) in FL.items():
        k = mask(mf).values; e_, m_ = ef.loc[:, k], mf[k].copy()
        idx = None if perm else rng_.integers(0, k.sum(), k.sum())
        if perm: m_["age"] = rng_.permutation(m_.age.values)
        F[name] = fit(e_, m_, cov, mm, kh, idx)
    G = {}
    for coh, (e_, m_) in _W["g"].items():
        m2 = m_.copy()
        idx = None if perm else rng_.integers(0, len(m_), len(m_))
        if perm: m2["age"] = rng_.permutation(m2.age.values)
        for name, (cov, skin_extra, kh) in GT.items():
            G[f"{coh}: {name}"] = fit(e_, m2, cov + (skin_extra if coh == "skin" else []), 1.0, kh, idx)
    return F, G


def _boot(b):
    F, G = _fit_all(np.random.default_rng(30_000 + b))
    return [dict(b=b, fn=fn, gn=gn, **compare(tf, tg, _W["sig_point"][fn], _W["chrom"])) for fn, tf in F.items() for gn, tg in G.items()]


def _perm(b):
    """Permuted-age null for the pre-specified models only (age shuffled within each cohort, covariates fixed)."""
    r = np.random.default_rng(40_000 + b); ef, mf = _W["f"]
    m_ = mf[mf.age >= 18].copy(); e_ = ef.loc[:, m_.index]; m_["age"] = r.permutation(m_.age.values)
    tf = fit(e_, m_, ["sex"])
    out = {"b": b, "fdr05_discovery": int((tf.q < 0.05).sum())}
    for coh, (eg, mg) in _W["g"].items():
        m2 = mg.copy(); m2["age"] = r.permutation(m2.age.values); tg = fit(eg, m2, ["sex"])
        sh = tf.index.intersection(tg.index)
        out[f"rho_{coh}"] = stats.spearmanr(tf.t[sh], tg.t[sh])[0]; out[f"fdr05_{coh}"] = int((tg.q < 0.05).sum())
        top = tf.t[sh].abs().nlargest(TOPK).index
        out[f"top_rho_{coh}"] = stats.spearmanr(tf.t[top], tg.t[top])[0]
    return out


def _subset(b):
    r = np.random.default_rng(50_000 + b); ef, mf = _W["f"]
    m_ = mf[mf.age >= 18]; idx = np.sort(r.choice(len(m_), 74, replace=False))
    return int((fit(ef.loc[:, m_.index], m_, ["sex"], idx=idx).q < 0.05).sum())


def vif_age(meta, covars):
    X, names = design(meta, covars)
    a = X[:, 1]; Z = np.delete(X, 1, axis=1)
    r2 = 1 - ((a - Z @ np.linalg.lstsq(Z, a, rcond=None)[0]) ** 2).sum() / ((a - a.mean()) ** 2).sum()
    return float(1 / (1 - r2))


def age_predict(ef, mf, cohorts):
    """X17: elastic net on z-scored expression of shared high-variance genes; cross-cohort and within-cohort CV error."""
    from sklearn.linear_model import ElasticNetCV, ElasticNet
    from sklearn.model_selection import cross_val_predict, KFold
    data = {"Discovery": (ef.loc[:, mf.age >= 18], mf[mf.age >= 18].age.values)}
    for k, (e_, m_) in cohorts.items(): data[k] = (e_, m_.age.values)
    genes = None
    for e_, _ in data.values():
        e1 = e_[e_.mean(axis=1) >= 1]; genes = e1.index if genes is None else genes.intersection(e1.index)
    var = data["Discovery"][0].loc[genes].var(axis=1); genes = var.nlargest(2000).index
    Z = {k: ((e_.loc[genes].T - e_.loc[genes].mean(axis=1)) / (e_.loc[genes].std(axis=1) + 1e-9)).values for k, (e_, a) in data.items()}
    out = []
    for tr in data:
        a_tr = data[tr][1]
        model = ElasticNetCV(l1_ratio=0.5, cv=5, alphas=np.logspace(-1, 1.5, 10), max_iter=3000, tol=1e-3, selection="random", random_state=0).fit(Z[tr], a_tr)
        cvp = cross_val_predict(ElasticNet(alpha=model.alpha_, l1_ratio=0.5, max_iter=3000, tol=1e-3, selection="random", random_state=0), Z[tr], a_tr,
                                cv=KFold(5, shuffle=True, random_state=0))
        out.append({"train": tr, "test": tr + " (5-fold CV)", "n_test": len(a_tr), "mae_years": float(np.abs(cvp - a_tr).mean()),
                    "pearson_r": float(np.corrcoef(cvp, a_tr)[0, 1]), "n_genes_used": int((model.coef_ != 0).sum())})
        for te in data:
            if te == tr: continue
            p = model.predict(Z[te]); a_te = data[te][1]
            out.append({"train": tr, "test": te, "n_test": len(a_te), "mae_years": float(np.abs(p - a_te).mean()),
                        "pearson_r": float(np.corrcoef(p, a_te)[0, 1]), "n_genes_used": int((model.coef_ != 0).sum())})
    return pd.DataFrame(out)


def main():
    ef, mf = fleischer(); cohorts = {"fibroblasts": gtex("gtex_fibroblasts_tpm.gct.gz"), "skin": gtex("gtex_skin_ns_tpm.gct.gz")}
    chrom = chrom_map()
    # point estimates
    F, G = {}, {}
    for name, (mask, cov, mm, kh) in FL.items():
        k = mask(mf).values; F[name] = (fit(ef.loc[:, k], mf[k], cov, mm, kh), int(k.sum()))
    for coh, (e_, m_) in cohorts.items():
        for name, (cov, skin_extra, kh) in GT.items():
            G[f"{coh}: {name}"] = (fit(e_, m_, cov + (skin_extra if coh == "skin" else []), 1.0, kh), len(m_))
    pd.to_pickle({fn: sig_from(v[0]) for fn, v in F.items()}, os.path.join(RES, "rev3_point_signatures.pkl"))
    pd.to_pickle({**{fn: v[0] for fn, v in F.items()}, **{gn: v[0] for gn, v in G.items()}}, os.path.join(RES, "rev3_point_tstats.pkl"))
    log("point estimates", {k: int((v[0].q < 0.05).sum()) for k, v in {**F, **G}.items()})
    out = {"n_boot": NB, "n_perm": NPERM, "n_subsets": NSUB, "hidden_factors": KPC}
    # X08 known genes
    rows = []
    for g in KNOWN + SEXGENES:
        row = {"gene": g, "set": "age" if g in KNOWN else "sex"}
        for lab, ts in [("Discovery", F["F0 pre-specified (age 18+; age + sex)"][0]), ("GTEx fibroblasts", G["fibroblasts: G0 pre-specified (age + sex)"][0]),
                        ("GTEx skin", G["skin: G0 pre-specified (age + sex)"][0])]:
            row[f"t_{lab}"] = float(ts.t.get(g, np.nan)); row[f"q_{lab}"] = float(ts.q.get(g, np.nan)); row[f"t_sex_{lab}"] = float(ts.t_sex.get(g, np.nan))
        rows.append(row)
    pd.DataFrame(rows).to_csv(os.path.join(RES, "rev3_cohort_known_genes.csv"), index=False)
    # X16 E2F/G2-M age t per cohort
    cc = []
    for lab, ts in [("Discovery", F["F0 pre-specified (age 18+; age + sex)"][0]), ("GTEx fibroblasts", G["fibroblasts: G0 pre-specified (age + sex)"][0]),
                    ("GTEx skin", G["skin: G0 pre-specified (age + sex)"][0])]:
        for setname in ["E2F Targets", "G2-M Checkpoint"]:
            gs = ts.index.intersection(pd.Index(sorted(HM[setname])))
            cc.extend({"cohort": lab, "set": setname, "gene": g, "t": float(ts.t[g])} for g in gs)
    CC = pd.DataFrame(cc); CC.to_csv(os.path.join(RES, "rev3_cohort_cellcycle_t.csv"), index=False)
    out["cellcycle_median_t"] = CC.groupby(["cohort", "set"]).t.median().round(3).reset_index().to_dict("records")
    sigF0 = F["F0 pre-specified (age 18+; age + sex)"][0]; sig = sig_from(sigF0)
    out["signature_frac_cellcycle"] = {"down": float(np.mean([g in CELLCYCLE for g in sig[sig.t < 0].index])), "up": float(np.mean([g in CELLCYCLE for g in sig[sig.t > 0].index]))}
    # X09 VIF and age per repository
    m18 = mf[mf.age >= 18]
    out["vif_age_F2"] = vif_age(m18, ["sex", "platform", "repository", "site2"])
    out["age_by_repository"] = m18.groupby("repository").age.describe()[["count", "min", "50%", "max"]].round(1).reset_index().to_dict("records")
    out["nia80_n"] = int(((m18.repository == "AG") & (m18.age >= 80)).sum()); out["n80_total"] = int((m18.age >= 80).sum())
    # X17 age prediction
    AP = age_predict(ef, mf, cohorts); AP.to_csv(os.path.join(RES, "rev3_cohort_age_prediction.csv"), index=False)
    log("age prediction done")
    with mp.Pool(NPROC, initializer=_init) as pool:
        subs = list(pool.imap_unordered(_subset, range(NSUB), chunksize=8)); log("power-matched subsets done")
        perms = list(pool.imap_unordered(_perm, range(NPERM), chunksize=4)); log("permutations done")
        pd.DataFrame(perms).to_csv(os.path.join(RES, "rev3_cohort_permutations.csv"), index=False)
        pd.Series(subs).to_csv(os.path.join(RES, "rev3_cohort_power_subsets.csv"), index=False)
        boot = []
        for i, res in enumerate(pool.imap_unordered(_boot, range(NB), chunksize=2)):
            boot.extend(res)
            if i % 25 == 0: log("bootstrap", i)
    BT = pd.DataFrame(boot); BT.to_csv(os.path.join(RES, "rev3_cohort_bootstrap.csv.gz"), index=False)
    PM = pd.DataFrame(perms); PM.to_csv(os.path.join(RES, "rev3_cohort_permutations.csv"), index=False)
    out["power_matched_fdr05"] = {"median": float(np.median(subs)), "q025": float(np.quantile(subs, 0.025)), "q975": float(np.quantile(subs, 0.975)),
                                  "frac_le_observed_20_79": float(np.mean(np.array(subs) <= int((F["F1 ages 20-79 (GTEx range)"][0].q < 0.05).sum())))}
    pd.Series(subs).to_csv(os.path.join(RES, "rev3_cohort_power_subsets.csv"), index=False)
    out["perm"] = {c: {"median": float(PM[c].median()), "q025": float(PM[c].quantile(0.025)), "q975": float(PM[c].quantile(0.975)), "max": float(PM[c].max())}
                   for c in PM.columns if c != "b"}
    rows = []
    for fn, (tf, nf) in F.items():
        for gn, (tg, ng) in G.items():
            pt = compare(tf, tg, sig_from(tf), chrom)
            B = BT[(BT.fn == fn) & (BT.gn == gn)]
            row = {"discovery_spec": fn, "n_discovery": nf, "fdr05_discovery": int((tf.q < 0.05).sum()), "replication_spec": gn,
                   "n_replication": ng, "fdr05_replication": int((tg.q < 0.05).sum()), **pt}
            for k in [c for c in pt if not c.startswith("n_")]:
                row[f"{k}_boot_mean"] = float(B[k].mean()); row[f"{k}_lo"] = float(B[k].quantile(0.025)); row[f"{k}_hi"] = float(B[k].quantile(0.975))
            rows.append(row)
    S = pd.DataFrame(rows); S.to_csv(os.path.join(RES, "rev3_cohort_sensitivity.csv"), index=False)
    tfib, tskin = G["fibroblasts: G0 pre-specified (age + sex)"][0], G["skin: G0 pre-specified (age + sex)"][0]
    sh = tfib.index.intersection(tskin.index)
    out["gtex_fibro_vs_skin_rho_G0"] = stats.spearmanr(tfib.t[sh], tskin.t[sh])[0]
    out["gtex_fibro_vs_skin_sex_rho_G0"] = stats.spearmanr(tfib.t_sex[sh], tskin.t_sex[sh])[0]
    out["fdr05"] = {**{k: int((v[0].q < 0.05).sum()) for k, v in F.items()}, **{k: int((v[0].q < 0.05).sum()) for k, v in G.items()}}
    out["n"] = {**{k: v[1] for k, v in F.items()}, **{k: v[1] for k, v in G.items()}}
    json.dump(out, open(os.path.join(RES, "rev3_cohort.json"), "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    log("done")


if __name__ == "__main__":
    main()
