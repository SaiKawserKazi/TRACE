"""Exploratory (Amendment 7, 6 November 2026): permuted-age null for the GTEx batch model G2.

Under G2 (age + sex + center, batch type, RIN, Hardy scale and nucleic-acid isolation batch; skin also ischemic time),
most GTEx age signal is absorbed by batch, so the top-200 agreement with the discovery cohort cannot be judged against
the G0 null. Here age is permuted within nucleic-acid isolation batch in GTEx (so batch and age stay unlinked under the
null, as the G2 model assumes) and across donors in the discovery cohort (F0), 1,000 times, and the G2 model is refitted.
Two-sided empirical P about zero, as for the other correlation nulls. Also records the number of FDR < 0.05 GTEx age
genes under the null, for comparison with the observed G2 counts.
Outputs: results/x_g2_null.json, results/x_g2_null.csv
"""
import os, sys, json, time
import multiprocessing as mp
import numpy as np, pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as _cfg
import rev3_cohort as C

RES = C.RES
NPERM = _cfg.N_SPEC_PERM; NPROC = _cfg.N_PROC; TOPK = 200
G2 = ["sex", "center", "batch_type", "rin", "hardy", "na_batch"]


def log(*a): print(time.strftime("%H:%M:%S"), *a, flush=True)


_W = {}


def _init():
    ef, mf = C.fleischer(); m18 = mf[mf.age >= 18].copy(); e18 = ef.loc[:, m18.index]
    _W["f"] = (e18, m18); _W["g"] = {"fibroblasts": C.gtex("gtex_fibroblasts_tpm.gct.gz"), "skin": C.gtex("gtex_skin_ns_tpm.gct.gz")}


def cov(coh): return G2 + (["ischemic"] if coh == "skin" else [])


def within(r, age, groups):
    a = age.copy()
    for g in np.unique(groups):
        i = np.where(groups == g)[0]; a[i] = r.permutation(a[i])
    return a


def top_rho(tf, tg):
    sh = tf.index.intersection(tg.index); top = tf.t[sh].abs().nlargest(TOPK).index
    return stats.spearmanr(tf.t[top], tg.t[top])[0]


def _perm(b):
    r = np.random.default_rng(80_000 + b); e, m = _W["f"]
    m_ = m.copy(); m_["age"] = r.permutation(m_.age.values); tf = C.fit(e, m_, ["sex"]); out = []
    for coh, (eg, mg) in _W["g"].items():
        m2 = mg.copy(); m2["age"] = within(r, m2.age.values.copy(), m2.na_batch.astype(str).values)
        tg = C.fit(eg, m2, cov(coh)); out.append({"b": b, "tissue": coh, "top_rho": top_rho(tf, tg), "fdr05_gtex": int((tg.q < 0.05).sum())})
    return out


def main():
    _init(); e, m = _W["f"]; tf = C.fit(e, m, ["sex"]); obs = {}
    for coh, (eg, mg) in _W["g"].items():
        tg = C.fit(eg, mg, cov(coh)); obs[coh] = {"top_rho": float(top_rho(tf, tg)), "fdr05_gtex": int((tg.q < 0.05).sum()),
                                                 "n_na_batches": int(mg.na_batch.nunique()), "single_donor_batches": int((mg.na_batch.value_counts() == 1).sum())}
    log("observed", obs)
    with mp.Pool(NPROC, initializer=_init) as pool:
        P = pd.DataFrame([row for res in pool.imap_unordered(_perm, range(NPERM), chunksize=4) for row in res])
    P.to_csv(os.path.join(RES, "x_g2_null.csv"), index=False)
    out = {"n_perm": NPERM, "topk": TOPK, "design": "age permuted across donors (discovery, F0) and within nucleic-acid isolation batch (GTEx); G2 refitted"}
    for coh in obs:
        n = P[P.tissue == coh].top_rho.values; n = n[np.isfinite(n)]; o = obs[coh]["top_rho"]
        out[coh] = {**obs[coh], "null_lo": float(np.quantile(n, 0.025)), "null_hi": float(np.quantile(n, 0.975)),
                    "p": float((1 + np.sum(np.abs(n) >= abs(o))) / (1 + len(n))),
                    "null_fdr05_gtex_median": float(P[P.tissue == coh].fdr05_gtex.median()), "null_fdr05_gtex_q95": float(P[P.tissue == coh].fdr05_gtex.quantile(0.95)),
                    "null_fdr05_gtex_max": int(P[P.tissue == coh].fdr05_gtex.max()), "frac_null_fdr05_ge_obs": float((P[P.tissue == coh].fdr05_gtex >= obs[coh]["fdr05_gtex"]).mean())}
    json.dump(out, open(os.path.join(RES, "x_g2_null.json"), "w"), indent=1); log("done"); print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
