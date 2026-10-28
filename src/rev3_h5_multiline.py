"""Exploratory (Amendment 3): H5 donor-bootstrap stability of Rrob on all genes
and on the genes measured in at least two cell lines (Rrob equals Radj for genes measured in one line).
Same design as rev_h5_boot.py (2,000 donor bootstraps, seeds 50,000 + b). Outputs results/rev3_h5_multiline.json/.csv
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

RES = T.RES; NB = _cfg.N_H5_BOOT_EXTENDED; NPROC = _cfg.N_PROC
_W = {}


def _init():
    e, m = aging.load_fleischer(); ctx = {c: perturb.load_context(c) for c in perturb.CONTEXTS}
    ident = T.identity_genes(ctx); progs = {c: T.program_shifts(v[0], ident[c]) for c, v in ctx.items()}
    rob = pd.read_csv(os.path.join(RES, "robust_scores.csv"), index_col=0)
    _W.update(e=e, m=m, ctx=ctx, progs=progs, rob=rob.R_rob, multi=rob.index[rob.n_lines >= 2])


def _task(b):
    r = np.random.default_rng(50_000 + b); e, m, rob, multi = _W["e"], _W["m"], _W["rob"], _W["multi"]
    ix = r.integers(0, len(m), len(m)); mb = m.iloc[ix].copy(); eb = e.iloc[:, ix]; eb.columns = mb.index = [f"s{i}" for i in range(len(ix))]
    _, rb = T.all_scores(_W["ctx"], aging.signature(aging.age_tstats(eb, mb)).t, _W["progs"])
    out = {"b": b}
    for lab, gs in [("all", rob.index), ("multi", multi)]:
        sh = gs.intersection(rb.index); k = max(1, int(0.05 * len(sh)))
        out[f"rho_{lab}"] = stats.spearmanr(rob[sh], rb.R_rob[sh])[0]
        out[f"top5_{lab}"] = len(set(rob[sh].nlargest(k).index) & set(rb.R_rob[sh].nlargest(k).index)) / k
    return out


if __name__ == "__main__":
    with mp.Pool(NPROC, initializer=_init) as pool:
        rows = []
        for i, r_ in enumerate(pool.imap_unordered(_task, range(NB), chunksize=2)):
            rows.append(r_)
            if i % 200 == 0: print(time.strftime("%H:%M:%S"), "boot", i, flush=True)
    D = pd.DataFrame(rows).sort_values("b"); D.to_csv(os.path.join(RES, "rev3_h5_multiline.csv"), index=False)
    out = {"n_boot": NB}
    for lab in ["all", "multi"]:
        out[lab] = {"median_rho": float(D[f"rho_{lab}"].median()), "q025": float(D[f"rho_{lab}"].quantile(0.025)), "q975": float(D[f"rho_{lab}"].quantile(0.975)),
                    "frac_ge_0.8": float((D[f"rho_{lab}"] >= 0.8).mean()), "median_top5": float(D[f"top5_{lab}"].median()),
                    "top5_q025": float(D[f"top5_{lab}"].quantile(0.025)), "top5_q975": float(D[f"top5_{lab}"].quantile(0.975))}
    json.dump(out, open(os.path.join(RES, "rev3_h5_multiline.json"), "w"), indent=1)
    print("done", out)
