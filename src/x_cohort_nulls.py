"""Exploratory (Amendment 4, logged 2026-10-28 before running): nulls and checks for the cohort comparison.

* Permuted-age null for every discovery specification F0-F9 against GTEx G0 (fibroblasts, skin): genome-wide rho and
  rho over the 200 strongest discovery genes; empirical two-sided P per specification and tissue.
* Permuted-sex null for the sex positive control (F0 x G0): sex labels shuffled within each cohort.
* Planted-effect check: with age permuted (no real shared age effect), a shared synthetic age effect of the size of the
  discovery age effects is added to 500 random genes in both cohorts; how often do the genome-wide and top-200 metrics
  exceed their permuted-age 97.5th percentile?
* Hidden factors of the discovery cohort (F9): association of factor 1 with age, repository and platform, and Hallmark
  enrichment of the genes that load on it.
Outputs: results/x_cohort_nulls.json, results/x_cohort_perm_by_spec.csv, results/x_cohort_sex_perm.csv,
results/x_cohort_planted.csv, results/x_cohort_hidden_factor_hallmark.csv
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import config as _cfg
import os, sys, json, time
import multiprocessing as mp
import numpy as np, pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
import rev3_cohort as C

RES = C.RES
NPERM = _cfg.N_SPEC_PERM; NPLANT = _cfg.N_PLANTED; NPROC = _cfg.N_PROC
TOPK, NGENE_PLANT = 200, 500


def log(*a): print(time.strftime("%H:%M:%S"), *a, flush=True)


_W = {}


def _init():
    _W["f"] = C.fleischer(); _W["g"] = {"fibroblasts": C.gtex("gtex_fibroblasts_tpm.gct.gz"), "skin": C.gtex("gtex_skin_ns_tpm.gct.gz")}
    _W["chrom"] = C.chrom_map()


def metrics(tf, tg, col="t"):
    sh = tf.index.intersection(tg.index)
    top = tf[col][sh].abs().nlargest(TOPK).index
    return stats.spearmanr(tf[col][sh], tg[col][sh])[0], stats.spearmanr(tf[col][top], tg[col][top])[0]


def _perm(b):
    r = np.random.default_rng(60_000 + b); ef, mf = _W["f"]; out = []
    G = {}
    for coh, (eg, mg) in _W["g"].items():
        m2 = mg.copy(); m2["age"] = r.permutation(m2.age.values); G[coh] = C.fit(eg, m2, ["sex"])
    for name, (mask, cov, mm, kh) in C.FL.items():
        k = mask(mf).values; m_ = mf[k].copy(); m_["age"] = r.permutation(m_.age.values)
        tf = C.fit(ef.loc[:, k], m_, cov, mm, kh)
        for coh, tg in G.items():
            ra, rt = metrics(tf, tg); out.append({"b": b, "spec": name, "tissue": coh, "rho_all": ra, "top_rho": rt, "fdr05": int((tf.q < 0.05).sum())})
    return out


def _sexperm(b):
    r = np.random.default_rng(61_000 + b); ef, mf = _W["f"]
    m_ = mf[mf.age >= 18].copy(); m_["sex"] = r.permutation(m_.sex.values); tf = C.fit(ef.loc[:, m_.index], m_, ["sex"])
    chrom = _W["chrom"]; out = {"b": b}
    for coh, (eg, mg) in _W["g"].items():
        m2 = mg.copy(); m2["sex"] = r.permutation(m2.sex.values); tg = C.fit(eg, m2, ["sex"])
        sh = tf.index.intersection(tg.index); c = chrom.reindex(sh); auto = sh[~c.isin(["X", "Y"]).values & c.notna().values]
        top = tf.t_sex[sh].abs().nlargest(TOPK).index; ta = tf.t_sex[auto].abs().nlargest(TOPK).index
        out.update({f"sex_rho_all_{coh}": stats.spearmanr(tf.t_sex[sh], tg.t_sex[sh])[0], f"sex_top_rho_{coh}": stats.spearmanr(tf.t_sex[top], tg.t_sex[top])[0],
                    f"sex_top_autosomal_rho_{coh}": stats.spearmanr(tf.t_sex[ta], tg.t_sex[ta])[0]})
    return out


def _plant(b, k_effect):
    """Permute age in both cohorts, then add k * SD(gene) * z(age) with a shared random sign to 500 shared random genes."""
    r = np.random.default_rng(62_000 + b); ef, mf = _W["f"]
    m_ = mf[mf.age >= 18].copy(); m_["age"] = r.permutation(m_.age.values); e_ = ef.loc[:, m_.index]
    fib, skin = _W["g"]["fibroblasts"], _W["g"]["skin"]
    pool = e_.index[e_.mean(axis=1) >= 1].intersection(fib[0].index[fib[0].mean(axis=1) >= 1]).intersection(skin[0].index[skin[0].mean(axis=1) >= 1])
    genes = r.choice(np.array(pool), NGENE_PLANT, replace=False); sign = r.choice([-1.0, 1.0], NGENE_PLANT)
    def planted(e, ages):
        z = (ages - ages.mean()) / ages.std(); e2 = e.copy()
        sd = e2.loc[genes].std(axis=1).values[:, None]
        e2.loc[genes] = e2.loc[genes].values + k_effect * sd * sign[:, None] * z[None, :]
        return e2
    tf = C.fit(planted(e_, m_.age.values.astype(float)), m_, ["sex"]); out = {"b": b, "k": k_effect}
    for coh, (eg, mg) in _W["g"].items():
        m2 = mg.copy(); m2["age"] = r.permutation(m2.age.values)
        tg = C.fit(planted(eg, m2.age.values.astype(float)), m2, ["sex"])
        ra, rt = metrics(tf, tg); out.update({f"rho_all_{coh}": ra, f"top_rho_{coh}": rt})
        sh = tf.index.intersection(tg.index); top = tf.t[sh].abs().nlargest(TOPK).index
        out[f"top_frac_planted_{coh}"] = float(np.isin(top, genes).mean())
    out["median_abs_t_planted_discovery"] = float(tf.t.reindex(genes).abs().median())
    return out


def _plant_task(args): return _plant(*args)


def hidden_factor_check():
    ef, mf = C.fleischer(); m_ = mf[mf.age >= 18]; e_ = ef.loc[:, m_.index]; e_ = e_[e_.mean(axis=1) >= 1.0]
    X, _ = C.design(m_, ["sex"]); Y = e_.values.T; Yz = (Y - Y.mean(0)) / (Y.std(0) + 1e-9)
    R = Yz - X @ np.linalg.lstsq(X, Yz, rcond=None)[0]; U, s, Vt = np.linalg.svd(R, full_matrices=False)
    rows, enr = [], []
    for k in range(C.KPC):
        f = U[:, k]
        rows.append({"factor": k + 1, "var_explained": float(s[k] ** 2 / (s ** 2).sum()), "rho_age": stats.spearmanr(f, m_.age)[0],
                     "auroc_repository_AG": float(stats.mannwhitneyu(f[(m_.repository == "AG").values], f[(m_.repository != "AG").values]).statistic
                                                  / ((m_.repository == "AG").sum() * (m_.repository != "AG").sum())),
                     "auroc_platform_NextSeq": float(stats.mannwhitneyu(f[(m_.platform == "NextSeq 500").values], f[(m_.platform != "NextSeq 500").values]).statistic
                                                     / ((m_.platform == "NextSeq 500").sum() * (m_.platform != "NextSeq 500").sum()))})
    load = pd.Series(np.corrcoef(U[:, 0], Yz.T)[0, 1:], index=e_.index)
    bg = set(e_.index); N = len(bg)
    for direction, genes in [("positive", set(load.nlargest(200).index)), ("negative", set(load.nsmallest(200).index))]:
        for name, gs in C.HM.items():
            gs = gs & bg; k = len(gs & genes)
            if len(gs) < 10: continue
            enr.append({"factor": 1, "direction": direction, "set": name, "overlap": k, "size": len(gs), "p": stats.hypergeom.sf(k - 1, N, len(gs), len(genes))})
    E = pd.DataFrame(enr); E["q"] = stats.false_discovery_control(E.p.values); E = E.sort_values("p")
    return pd.DataFrame(rows), E


def emp_p(obs, null):
    null = np.asarray(null); null = null[np.isfinite(null)]
    return float((1 + np.sum(np.abs(null) >= abs(obs))) / (1 + len(null)))


def main():
    out = {"n_perm": NPERM, "n_plant": NPLANT, "topk": TOPK, "n_genes_planted": NGENE_PLANT}
    HF, HE = hidden_factor_check(); HF.to_csv(os.path.join(RES, "x_cohort_hidden_factors.csv"), index=False)
    HE.to_csv(os.path.join(RES, "x_cohort_hidden_factor_hallmark.csv"), index=False); log("hidden factors done")
    with mp.Pool(NPROC, initializer=_init) as pool:
        P = pd.DataFrame([row for res in pool.imap_unordered(_perm, range(NPERM), chunksize=4) for row in res]); log("age permutations done")
        P.to_csv(os.path.join(RES, "x_cohort_perm_by_spec.csv"), index=False)
        SX = pd.DataFrame(list(pool.imap_unordered(_sexperm, range(NPERM), chunksize=4))); log("sex permutations done")
        SX.to_csv(os.path.join(RES, "x_cohort_sex_perm.csv"), index=False)
        tasks = [(b, k) for k in (0.20, 0.41) for b in range(NPLANT)]
        PL = pd.DataFrame(list(pool.imap_unordered(_plant_task, tasks, chunksize=2))); log("planted done")
        PL.to_csv(os.path.join(RES, "x_cohort_planted.csv"), index=False)
    S = pd.read_csv(os.path.join(RES, "rev3_cohort_sensitivity.csv"))
    G0 = {"fibroblasts": "fibroblasts: G0 pre-specified (age + sex)", "skin": "skin: G0 pre-specified (age + sex)"}
    rows = []
    for spec in C.FL:
        for coh, gn in G0.items():
            r = S[(S.discovery_spec == spec) & (S.replication_spec == gn)].iloc[0]; nul = P[(P.spec == spec) & (P.tissue == coh)]
            rows.append({"spec": spec, "tissue": coh, "rho_all": r.rho_all, "p_rho_all": emp_p(r.rho_all, nul.rho_all),
                         "null_rho_all_lo": nul.rho_all.quantile(0.025), "null_rho_all_hi": nul.rho_all.quantile(0.975),
                         "top_rho": r.top_rho, "p_top_rho": emp_p(r.top_rho, nul.top_rho),
                         "null_top_lo": nul.top_rho.quantile(0.025), "null_top_hi": nul.top_rho.quantile(0.975)})
    SP = pd.DataFrame(rows); SP.to_csv(os.path.join(RES, "x_cohort_perm_p.csv"), index=False)
    f0 = S[S.discovery_spec.str.startswith("F0")]
    sexp = {}
    for coh, gn in G0.items():
        r = f0[f0.replication_spec == gn].iloc[0]
        for k in ["sex_rho_all", "sex_top_rho", "sex_top_autosomal_rho"]:
            sexp[f"{k}_{coh}"] = {"observed": r[k], "p": emp_p(r[k], SX[f"{k}_{coh}"]), "null_lo": SX[f"{k}_{coh}"].quantile(0.025), "null_hi": SX[f"{k}_{coh}"].quantile(0.975)}
        sexp[f"sex_rho_autosomes_all_genes_{coh}"] = float(r.sex_rho_autosomes)
    out["sex_perm"] = sexp
    f0n = P[P.spec.str.startswith("F0")]
    rec = {}
    for k_, g in PL.groupby("k"):
        rec[str(k_)] = {"median_abs_t_planted_discovery": float(g.median_abs_t_planted_discovery.median())}
        for coh in G0:
            nn = f0n[f0n.tissue == coh]
            rec[str(k_)][coh] = {"frac_rho_all_beyond_null": float((g[f"rho_all_{coh}"] > nn.rho_all.quantile(0.975)).mean()),
                                 "frac_top_rho_beyond_null": float((g[f"top_rho_{coh}"] > nn.top_rho.quantile(0.975)).mean()),
                                 "median_rho_all": float(g[f"rho_all_{coh}"].median()), "median_top_rho": float(g[f"top_rho_{coh}"].median()),
                                 "median_top_frac_planted": float(g[f"top_frac_planted_{coh}"].median())}
    out["planted"] = rec
    PM = pd.read_csv(os.path.join(RES, "rev3_cohort_permutations.csv"))
    out["perm_fdr05_discovery_quantiles"] = {q: float(PM.fdr05_discovery.quantile(float(q))) for q in ["0.5", "0.95", "0.99"]} | {"max": int(PM.fdr05_discovery.max()),
                                                                                                                                 "frac_zero": float((PM.fdr05_discovery == 0).mean())}
    json.dump(out, open(os.path.join(RES, "x_cohort_nulls.json"), "w"), indent=1, default=float)
    log("done"); print(SP.round(3).to_string()); print(json.dumps(out, indent=1, default=float)[:3000])


if __name__ == "__main__":
    main()
