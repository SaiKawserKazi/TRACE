"""Exploratory (Amendment 2): H5 donor bootstrap extended from the pre-specified 200 to 2,000 replicates.

The pre-specified 200-replicate result (results.json) is kept as the confirmatory test; this run tightens the
uncertainty estimate. Each replicate resamples discovery donors with replacement, re-derives the signature and
recomputes R_rob in every screen exactly as in trace_analysis. Outputs results/rev_h5_boot.json and .csv.
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import config as _cfg
import os, sys, json, time
import multiprocessing as mp
import numpy as np, pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
import aging, perturb
import trace_analysis as T

RES = T.RES
NB = _cfg.N_H5_BOOT_EXTENDED
NPROC = _cfg.N_PROC
_W = {}


def log(*a): print(time.strftime("%H:%M:%S"), *a, flush=True)


def _init():
    e, m = aging.load_fleischer()
    ctx = {c: perturb.load_context(c) for c in perturb.CONTEXTS}
    ident = T.identity_genes(ctx)
    progs = {c: T.program_shifts(v[0], ident[c]) for c, v in ctx.items()}
    rob = pd.read_csv(os.path.join(RES, "robust_scores.csv"), index_col=0).R_rob
    _W.update(e=e, m=m, ctx=ctx, progs=progs, rob=rob)


def _task(b):
    r = np.random.default_rng(50_000 + b)
    e, m, rob = _W["e"], _W["m"], _W["rob"]
    ix = r.integers(0, len(m), len(m))
    mb = m.iloc[ix].copy(); eb = e.iloc[:, ix]; eb.columns = mb.index = [f"s{i}" for i in range(len(ix))]
    sb = aging.signature(aging.age_tstats(eb, mb))
    _, rb = T.all_scores(_W["ctx"], sb.t, _W["progs"])
    sh = rob.index.intersection(rb.index)
    k = max(1, int(0.05 * rob.notna().sum()))
    top_full = set(rob.nlargest(k).index); tb = set(rb.R_rob.nlargest(k).index)
    return b, stats.spearmanr(rob[sh], rb.R_rob[sh])[0], len(tb & top_full) / k, sorted(tb)


if __name__ == "__main__":
    rows, cnt = [], {}
    with mp.Pool(NPROC, initializer=_init) as pool:
        for i, (b, rho, ov, tb) in enumerate(pool.imap_unordered(_task, range(NB), chunksize=2)):
            rows.append((b, rho, ov))
            for g in tb: cnt[g] = cnt.get(g, 0) + 1
            if i % 100 == 0: log("boot", i)
    D = pd.DataFrame(rows, columns=["b", "rho", "top5_overlap"]).sort_values("b")
    D.to_csv(os.path.join(RES, "rev_h5_boot.csv"), index=False)
    pd.Series(cnt).div(NB).rename("bootstrap_top5_support").to_csv(os.path.join(RES, "rev_h5_top5_support.csv"))
    # Monte-Carlo precision of the median: bootstrap of the bootstrap distribution
    rr = np.random.default_rng(1); med = [np.median(rr.choice(D.rho, len(D))) for _ in range(2000)]
    out = {"n_boot": NB, "median_rho": float(D.rho.median()), "q05_rho": float(D.rho.quantile(0.05)), "q95_rho": float(D.rho.quantile(0.95)),
           "q025_rho": float(D.rho.quantile(0.025)), "q975_rho": float(D.rho.quantile(0.975)),
           "median_rho_mc_ci": [float(np.quantile(med, 0.025)), float(np.quantile(med, 0.975))],
           "frac_ge_0.8": float((D.rho >= 0.8).mean()), "median_top5_overlap": float(D.top5_overlap.median()),
           "stability_criterion_met": bool(D.rho.median() >= 0.8)}
    json.dump(out, open(os.path.join(RES, "rev_h5_boot.json"), "w"), indent=1)
    log("done", out)
