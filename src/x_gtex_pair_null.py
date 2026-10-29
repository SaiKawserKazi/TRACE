"""Permuted-age null for the agreement of age effects between GTEx fibroblasts and GTEx skin (exploratory, Amendment 5).

Pre-specified model (expression ~ age + sex) in each tissue; age is shuffled within each tissue 1,000 times with sex
fixed, and the genome-wide Spearman correlation of age t-statistics between the two tissues is recomputed.
Writes results/x_gtex_pair_null.json."""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy import stats
import aging
from config import RES, N_SPEC_PERM


def tvals(Y, age, male):
    X = np.column_stack([np.ones(len(age)), age, male]); XtXi = np.linalg.inv(X.T @ X)
    B = XtXi @ X.T @ Y; r = Y - X @ B; s2 = (r ** 2).sum(0) / (X.shape[0] - 3)
    return B[1] / np.sqrt(s2 * XtXi[1, 1])


def main():
    D = {}
    for k, fn in [("fibroblasts", "gtex_fibroblasts_tpm.gct.gz"), ("skin", "gtex_skin_ns_tpm.gct.gz")]:
        e, m = aging.load_gtex(fn); e = e[e.mean(axis=1) >= 1.0]; D[k] = (e, m)
    genes = D["fibroblasts"][0].index.intersection(D["skin"][0].index)
    Y = {k: D[k][0].loc[genes].values.T for k in D}; age = {k: D[k][1].age.values.astype(float) for k in D}
    male = {k: (D[k][1].sex.values == "M").astype(float) for k in D}
    obs = stats.spearmanr(tvals(Y["fibroblasts"], age["fibroblasts"], male["fibroblasts"]), tvals(Y["skin"], age["skin"], male["skin"]), nan_policy="omit")[0]
    r = np.random.default_rng(64_000); null = []
    for _ in range(N_SPEC_PERM):
        tf = tvals(Y["fibroblasts"], r.permutation(age["fibroblasts"]), male["fibroblasts"]); ts = tvals(Y["skin"], r.permutation(age["skin"]), male["skin"])
        null.append(stats.spearmanr(tf, ts, nan_policy="omit")[0])
    null = np.array(null)
    out = {"n_genes": int(len(genes)), "n_perm": int(len(null)), "observed_rho": float(obs), "null_q025": float(np.quantile(null, 0.025)),
           "null_q975": float(np.quantile(null, 0.975)), "p_two_sided": float((1 + (np.abs(null) >= abs(obs)).sum()) / (1 + len(null)))}
    json.dump(out, open(RES / "x_gtex_pair_null.json", "w"), indent=1); print(out)


if __name__ == "__main__":
    main()
