"""Exploratory (Amendment 7, 6 November 2026): the cohort check of Box 1 for the GTEx-derived signature.

Agreement of age effects among the 200 strongest GTEx-fibroblast age genes (main model G0, age + sex) with the discovery
cohort (F0), against the permuted-age null. The permutations are the same 1,000 age labelings as the F0 test
(x_cohort_nulls.py: seed 60,000 + b; draw order GTEx fibroblasts, GTEx skin, discovery). Two-sided empirical P about zero.
Outputs: results/x_reverse_top200.json, results/x_reverse_top200_null.csv
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
_W = {}


def _init():
    ef, mf = C.fleischer(); m18 = mf[mf.age >= 18].copy(); _W["f"] = (ef.loc[:, m18.index], m18)
    _W["fib"] = C.gtex("gtex_fibroblasts_tpm.gct.gz"); _W["n_skin"] = len(C.gtex("gtex_skin_ns_tpm.gct.gz")[1])


def rev_rho(tg, tf):
    sh = tg.index.intersection(tf.index); top = tg.t[sh].abs().nlargest(TOPK).index
    return stats.spearmanr(tg.t[top], tf.t[top])[0]


def _perm(b):
    r = np.random.default_rng(60_000 + b); e, m = _W["f"]; eg, mg = _W["fib"]
    m2 = mg.copy(); m2["age"] = r.permutation(m2.age.values); tg = C.fit(eg, m2, ["sex"])
    r.permutation(_W["n_skin"])                                       # the skin draw, kept so the discovery permutation matches
    m_ = m.copy(); m_["age"] = r.permutation(m_.age.values); tf = C.fit(e, m_, ["sex"])
    return {"b": b, "rho": rev_rho(tg, tf), "f0_top_rho": rev_rho(tf, tg)}


def main():
    _init(); e, m = _W["f"]; eg, mg = _W["fib"]
    obs = float(rev_rho(C.fit(eg, mg, ["sex"]), C.fit(e, m, ["sex"])))
    with mp.Pool(NPROC, initializer=_init) as pool:
        P = pd.DataFrame(list(pool.imap_unordered(_perm, range(NPERM), chunksize=4)))
    P.to_csv(os.path.join(RES, "x_reverse_top200_null.csv"), index=False)
    n = P.rho.values[np.isfinite(P.rho.values)]
    # check: the F0 direction of the same permutations must reproduce the main null
    ref = pd.read_csv(os.path.join(RES, "x_cohort_perm_by_spec.csv")); ref = ref[ref.spec.str.startswith("F0") & (ref.tissue == "fibroblasts")].set_index("b").top_rho
    same = bool(np.allclose(P.set_index("b").f0_top_rho.sort_index(), ref.sort_index(), atol=1e-9))
    out = {"n_perm": NPERM, "topk": TOPK, "rho": obs, "null_lo": float(np.quantile(n, 0.025)), "null_hi": float(np.quantile(n, 0.975)),
           "p": float((1 + np.sum(np.abs(n) >= abs(obs))) / (1 + len(n))), "same_permutations_as_F0_test": same}
    json.dump(out, open(os.path.join(RES, "x_reverse_top200.json"), "w"), indent=1); print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
