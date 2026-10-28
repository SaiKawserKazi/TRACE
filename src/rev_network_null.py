"""Exploratory (Amendment 2): 500 degree-preserving randomised STRING networks for the H3a network test.

For every randomisation and screen we recompute the network predictor exactly as in the pre-specified analysis
(same five folds, same training-fold rule) and record its Spearman rho and MAE against held-out R_adj. Outputs:
results/rev_network_null_dist.csv (one row per randomisation x screen), results/rev_network_null.json (empirical
P values, null quantiles, and H3a re-tested against the mean prediction over all randomisations).
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import config as _cfg
import os, sys, json, time
import multiprocessing as mp
import numpy as np, pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
import trace_analysis as T

RES = T.RES
CTX = ["K562", "K562gw", "RPE1", "HepG2", "Jurkat"]
NNULL = _cfg.N_NETWORK_NULLS
NPROC = _cfg.N_PROC
_W = {}


def log(*a): print(time.strftime("%H:%M:%S"), *a, flush=True)


def _folds(y, G):
    idx = np.array([p for p in y.index if p in G])
    folds = np.random.default_rng(0).permutation(len(idx)) % 5  # identical to trace_analysis.h3_network
    return idx, folds


def predict(G, y, idx, folds):
    out = {}
    for f in range(5):
        tr, te = idx[folds != f], idx[folds == f]
        out.update(T.neighbour_mean(G, y[tr].to_dict(), te).to_dict())
    return pd.Series(out, dtype=float)


def _init():
    G = T.load_string()
    ys = {c: pd.read_csv(os.path.join(RES, f"scores_{c}.csv"), index_col=0).R_adj for c in CTX}
    _W.update(G=G, ys=ys, fl={c: _folds(ys[c], G) for c in CTX}, deg=dict(G.degree()))


def _task(seed):
    G = _W["G"]; H = T.randomised(G, 1000 + seed)
    degs_ok = all(H.degree(n) == d for n, d in _W["deg"].items())
    rows, preds = [], {}
    for c in CTX:
        y = _W["ys"][c]; idx, folds = _W["fl"][c]
        pr = predict(H, y, idx, folds)
        preds[c] = pr
        rows.append((seed, c, degs_ok, len(pr)))
    return rows, preds


def main():
    G = T.load_string(); ys = {c: pd.read_csv(os.path.join(RES, f"scores_{c}.csv"), index_col=0).R_adj for c in CTX}
    real = {c: predict(G, ys[c], *_folds(ys[c], G)) for c in CTX}
    log("real network predictions done")
    psum = {c: pd.Series(0.0, index=real[c].index) for c in CTX}; pcnt = {c: pd.Series(0, index=real[c].index) for c in CTX}
    dist = []
    with mp.Pool(NPROC, initializer=_init) as pool:
        for i, (rows, preds) in enumerate(pool.imap_unordered(_task, range(NNULL))):
            for (seed, c, ok, n), pr in zip(rows, [preds[c_] for c_ in CTX]):
                y = ys[c]
                common = real[c].index.intersection(pr.index)
                rho = stats.spearmanr(pr[common], y[common])[0]
                mae = float((pr[common] - y[common]).abs().mean())
                dist.append({"seed": seed, "screen": c, "degrees_preserved": ok, "n_eval": len(common), "rho": rho, "mae": mae})
                pr = pr.reindex(psum[c].index)
                psum[c] = psum[c].add(pr.fillna(0)); pcnt[c] = pcnt[c].add(pr.notna().astype(int))
            if i % 25 == 0: log("null", i)
    D = pd.DataFrame(dist); D.to_csv(os.path.join(RES, "rev_network_null_dist.csv"), index=False)
    out = {"n_null": NNULL, "all_degrees_preserved": bool(D.degrees_preserved.all()), "screens": {}}
    pvals = []
    for c in CTX:
        y = ys[c]; r = real[c]; d = D[D.screen == c]
        rho_real = stats.spearmanr(r, y[r.index])[0]; mae_real = float((r - y[r.index]).abs().mean())
        mean_null = (psum[c] / pcnt[c].replace(0, np.nan)).reindex(r.index)
        ok = mean_null.notna()
        e_r, e_n = (r[ok] - y[r.index][ok]).abs(), (mean_null[ok] - y[r.index][ok]).abs()
        w = stats.wilcoxon(e_r, e_n, alternative="less").pvalue
        pvals.append(w)
        out["screens"][c] = {"n_eval": int(ok.sum()), "rho_network": rho_real, "mae_network": mae_real,
                             "null_rho_median": float(d.rho.median()), "null_rho_q025": float(d.rho.quantile(0.025)), "null_rho_q975": float(d.rho.quantile(0.975)),
                             "null_rho_max": float(d.rho.max()), "empirical_p_rho": float((np.sum(d.rho >= rho_real) + 1) / (len(d) + 1)),
                             "null_mae_median": float(d.mae.median()), "null_mae_q025": float(d.mae.quantile(0.025)), "null_mae_q975": float(d.mae.quantile(0.975)),
                             "empirical_p_mae": float((np.sum(d.mae <= mae_real) + 1) / (len(d) + 1)),
                             "h3a_wilcoxon_vs_mean_null_p": w, "mae_mean_null": float(e_n.mean())}
    for c, q in zip(CTX, T.holm(pvals)): out["screens"][c]["h3a_wilcoxon_vs_mean_null_p_holm"] = q
    out["h3a_supported_500"] = sum(1 for c in CTX if out["screens"][c]["h3a_wilcoxon_vs_mean_null_p_holm"] < 0.05) > len(CTX) / 2
    json.dump(out, open(os.path.join(RES, "rev_network_null.json"), "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    log("done", {c: (round(v["rho_network"], 3), round(v["null_rho_q975"], 3), v["empirical_p_rho"]) for c, v in out["screens"].items()})


if __name__ == "__main__":
    main()
