"""Exploratory (Amendment 6, 3 November 2026): further checks of the gene-level cohort comparison.

* X30 Permuted-age null per k: agreement among the k strongest discovery genes (k = 50 to 1,000) for F0 against GTEx G0,
  both tissues, with the permuted-age 95% range at each k; and a permuted-age null for the sign agreement of the
  signature genes (age permuted in both cohorts; the signature is rebuilt from the permuted discovery fit, as in the
  observed analysis).
* X31 Hidden factors computed without protecting age: the first five principal components of the residuals after a
  sex-only model (F9b). Their correlation with age, and top-200 agreement under F9b with its own permuted-age null.
* X32 Planted shared effect added to the real (not permuted) data of both cohorts at two sizes (planted-only |t| about
  4 and 6 in the discovery cohort): how often planted genes enter the discovery top 200, and how often top-200 agreement
  turns positive or exceeds the permuted-age 97.5th percentile.
* X33 Variance inflation factor of age in the GTEx models with batch covariates (G1, G2; and age + sex + nucleic-acid
  isolation batch).
Outputs: results/x_cohort_checks_v2.json, x_cohort_topk_null.csv, x_cohort_planted_real.csv
"""
import os, sys, json, time
import multiprocessing as mp
import numpy as np, pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as _cfg
import rev3_cohort as C

RES = C.RES
NPERM = _cfg.N_SPEC_PERM; NPLANT = _cfg.N_PLANTED; NPROC = _cfg.N_PROC
KS = [50, 100, 200, 300, 500, 1000]
NGENE_PLANT, TOPK = 500, 200
PLANT_SIZES = (0.38, 0.57)          # multiples of the gene SD per SD of age; give planted-only |t| of about 4 and 6


def log(*a): print(time.strftime("%H:%M:%S"), *a, flush=True)


def sex_only_factors(e, m, k=C.KPC):
    """Principal components of z-scored expression after regressing out sex only (age not protected)."""
    X, _ = C.design(m, ["sex"]); X = X[:, [0] + list(range(2, X.shape[1]))]      # drop the age column
    Y = e.values.T; Yz = (Y - Y.mean(0)) / (Y.std(0) + 1e-9)
    U = np.linalg.svd(Yz - X @ np.linalg.lstsq(X, Yz, rcond=None)[0], full_matrices=False)[0]
    return U[:, :k]


def fit_extra(e, m, covars, F):
    """C.fit with extra numeric columns F (samples x k) added to the design."""
    m2 = m.copy(); cols = []
    for j in range(F.shape[1]): m2[f"sv{j}"] = F[:, j]; cols.append(f"sv{j}")
    return C.fit(e, m2, covars + cols)


def topk_stats(tf, tg):
    sh = tf.index.intersection(tg.index); order = tf.t[sh].abs().sort_values(ascending=False).index; out = {}
    for k in KS:
        g = order[:k]; out[f"rho_k{k}"] = stats.spearmanr(tf.t[g], tg.t[g])[0]; out[f"sign_k{k}"] = float((np.sign(tf.t[g]) == np.sign(tg.t[g])).mean())
    sig = C.sig_from(tf); s = sig.index.intersection(tg.index)
    out["sig_sign_agree"] = float((np.sign(sig.t[s]) == np.sign(tg.t[s])).mean())
    return out


_W = {}


def _init():
    ef, mf = C.fleischer(); m18 = mf[mf.age >= 18].copy(); e18 = ef.loc[:, m18.index]; e18 = e18[e18.mean(axis=1) >= 1.0]
    _W["f"] = (e18, m18); _W["sv"] = sex_only_factors(e18, m18)
    _W["g"] = {"fibroblasts": C.gtex("gtex_fibroblasts_tpm.gct.gz"), "skin": C.gtex("gtex_skin_ns_tpm.gct.gz")}


def _perm(b):
    r = np.random.default_rng(70_000 + b); e, m = _W["f"]
    m_ = m.copy(); m_["age"] = r.permutation(m_.age.values)
    tf = C.fit(e, m_, ["sex"]); tf9 = fit_extra(e, m_, ["sex"], _W["sv"]); out = []
    for coh, (eg, mg) in _W["g"].items():
        m2 = mg.copy(); m2["age"] = r.permutation(m2.age.values); tg = C.fit(eg, m2, ["sex"])
        sh = tf9.index.intersection(tg.index); top = tf9.t[sh].abs().nlargest(TOPK).index
        out.append({"b": b, "tissue": coh, **topk_stats(tf, tg), "F9b_top_rho": stats.spearmanr(tf9.t[top], tg.t[top])[0]})
    return out


def _plant(args):
    b, k_eff = args
    r = np.random.default_rng(71_000 + b); e, m = _W["f"]; fib, skin = _W["g"]["fibroblasts"], _W["g"]["skin"]
    pool = e.index.intersection(fib[0].index[fib[0].mean(axis=1) >= 1]).intersection(skin[0].index[skin[0].mean(axis=1) >= 1])
    genes = r.choice(np.array(pool), NGENE_PLANT, replace=False); sign = r.choice([-1.0, 1.0], NGENE_PLANT)

    def planted(e_, ages):
        z = (ages - ages.mean()) / ages.std(); e2 = e_.copy(); sd = e2.loc[genes].std(axis=1).values[:, None]
        e2.loc[genes] = e2.loc[genes].values + k_eff * sd * sign[:, None] * z[None, :]
        return e2
    tf = C.fit(planted(e, m.age.values.astype(float)), m, ["sex"])
    # size of the planted part alone: the same addition on age-permuted data
    mp_ = m.copy(); mp_["age"] = r.permutation(m.age.values); tp = C.fit(planted(e, mp_.age.values.astype(float)), mp_, ["sex"])
    out = {"b": b, "k": k_eff, "median_abs_t_planted_only": float(tp.t.reindex(genes).abs().median()),
           "median_abs_t_planted_genes_real": float(tf.t.reindex(genes).abs().median())}
    for coh, (eg, mg) in _W["g"].items():
        tg = C.fit(planted(eg, mg.age.values.astype(float)), mg, ["sex"])
        sh = tf.index.intersection(tg.index); top = tf.t[sh].abs().nlargest(TOPK).index
        out.update({f"top_rho_{coh}": stats.spearmanr(tf.t[top], tg.t[top])[0], f"top_frac_planted_{coh}": float(np.isin(top, genes).mean())})
    return out


def main():
    _init(); e, m = _W["f"]; sv = _W["sv"]
    out = {"n_perm": NPERM, "n_plant": NPLANT, "ks": KS, "plant_sizes": PLANT_SIZES, "n_genes_planted": NGENE_PLANT}
    # observed
    tf = C.fit(e, m, ["sex"]); tf9 = fit_extra(e, m, ["sex"], sv); obs = {}
    G0 = {coh: C.fit(eg, mg, ["sex"]) for coh, (eg, mg) in _W["g"].items()}
    for coh, tg in G0.items():
        sh = tf9.index.intersection(tg.index); top = tf9.t[sh].abs().nlargest(TOPK).index
        obs[coh] = {**topk_stats(tf, tg), "F9b_top_rho": stats.spearmanr(tf9.t[top], tg.t[top])[0]}
    out["observed"] = obs
    out["F9b_factors"] = [{"factor": j + 1, "rho_age": float(stats.spearmanr(sv[:, j], m.age)[0]),
                           "auroc_repository_AG": float(stats.mannwhitneyu(sv[(m.repository == "AG").values, j], sv[(m.repository != "AG").values, j]).statistic
                                                        / ((m.repository == "AG").sum() * (m.repository != "AG").sum()))} for j in range(sv.shape[1])]
    out["F9b_n_fdr05"] = int((tf9.q < 0.05).sum())
    # X33 VIF of age in GTEx batch models
    vif = {}
    for coh, (eg, mg) in _W["g"].items():
        extra = ["ischemic"] if coh == "skin" else []
        vif[coh] = {"G1": C.vif_age(mg, ["sex", "center", "batch_type", "rin", "hardy"] + extra),
                    "G2": C.vif_age(mg, ["sex", "center", "batch_type", "rin", "hardy", "na_batch"] + extra),
                    "age_sex_na_batch": C.vif_age(mg, ["sex", "na_batch"]), "n_na_batches": int(mg.na_batch.nunique())}
    out["vif_age_gtex"] = vif; log("observed and VIF done", json.dumps(vif))
    with mp.Pool(NPROC, initializer=_init) as pool:
        P = pd.DataFrame([row for res in pool.imap_unordered(_perm, range(NPERM), chunksize=4) for row in res]); log("permutations done")
        PL = pd.DataFrame(list(pool.imap_unordered(_plant, [(b, k) for k in PLANT_SIZES for b in range(NPLANT)], chunksize=2))); log("planted done")
    P.to_csv(os.path.join(RES, "x_cohort_topk_null.csv"), index=False); PL.to_csv(os.path.join(RES, "x_cohort_planted_real.csv"), index=False)
    emp = lambda o, n: float((1 + np.sum(np.abs(n - np.nanmedian(n)) >= abs(o - np.nanmedian(n)))) / (1 + np.isfinite(n).sum()))
    summ = {}
    for coh in G0:
        n = P[P.tissue == coh]; d = {}
        for k in KS:
            d[f"k{k}"] = {"rho": obs[coh][f"rho_k{k}"], "null_lo": float(n[f"rho_k{k}"].quantile(0.025)), "null_hi": float(n[f"rho_k{k}"].quantile(0.975)),
                          "p": emp(obs[coh][f"rho_k{k}"], n[f"rho_k{k}"].values)}
        d["sig_sign_agree"] = {"obs": obs[coh]["sig_sign_agree"], "null_lo": float(n.sig_sign_agree.quantile(0.025)),
                               "null_hi": float(n.sig_sign_agree.quantile(0.975)), "p": emp(obs[coh]["sig_sign_agree"], n.sig_sign_agree.values)}
        d["F9b_top_rho"] = {"obs": obs[coh]["F9b_top_rho"], "null_lo": float(n.F9b_top_rho.quantile(0.025)), "null_hi": float(n.F9b_top_rho.quantile(0.975)),
                            "p": emp(obs[coh]["F9b_top_rho"], n.F9b_top_rho.values)}
        summ[coh] = d
    out["nulls"] = summ
    pl = {}
    for k_, g in PL.groupby("k"):
        pl[str(k_)] = {"median_abs_t_planted_only": float(g.median_abs_t_planted_only.median()),
                       "median_abs_t_planted_genes_real": float(g.median_abs_t_planted_genes_real.median())}
        for coh in G0:
            hi = summ[coh]["k200"]["null_hi"]
            pl[str(k_)][coh] = {"median_top_frac_planted": float(g[f"top_frac_planted_{coh}"].median()), "median_top_rho": float(g[f"top_rho_{coh}"].median()),
                                "frac_top_rho_positive": float((g[f"top_rho_{coh}"] > 0).mean()), "frac_top_rho_beyond_null": float((g[f"top_rho_{coh}"] > hi).mean())}
    out["planted_real"] = pl
    json.dump(out, open(os.path.join(RES, "x_cohort_checks_v2.json"), "w"), indent=1, default=float)
    log("done"); print(json.dumps(out, indent=1, default=float)[:6000])


if __name__ == "__main__":
    main()
