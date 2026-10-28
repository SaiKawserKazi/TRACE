"""Exploratory (Amendment 3): network prediction, redone.

E11 in-fold fitting: the target (V1 adjusted score) is residualized on proliferation, identity and magnitude with
    coefficients fitted on the training fold only; ridge penalties tuned by inner 5-fold CV; 10 CV seeds; one common
    subset of knockdowns per screen for every predictor. Line transfer uses only other cell lines (for K562 essential it
    never includes the K562 genome-wide screen, and vice versa).
E12 baselines: constant, line transfer, calibrated transfer, STRING network (combined score >= 700), calibrated network,
    STRING physical network, transfer + network, k-nearest neighbors in the screen's own response space (non-signature
    genes only), gene covariates, and proliferation shift as a predictor of the raw score. R^2 is reported next to the
    reliability ceiling. Block bootstrap over response clusters for the network increment.
Outputs: results/rev3_network_*.csv, results/rev3_network.json
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import config as _cfg
import os, sys, json, time
import numpy as np, pandas as pd, networkx as nx
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
import aging, perturb
import trace_analysis as T

RAW, RES = T.RAW, T.RES
CTX = ["K562", "K562gw", "RPE1", "HepG2", "Jurkat"]; LINE = T.LINE_OF
NSEED = _cfg.N_NETWORK_SEEDS; NBLOCK = _cfg.N_BLOCK_BOOT; KNN = 20
LAMBDAS = np.logspace(-6, 3, 19)  # widened (Amendment 4): the old lower edge 1e-3 was chosen almost always


def log(*a): print(time.strftime("%H:%M:%S"), *a, flush=True)


def ridge_fit(X, y, lam):
    A = np.column_stack([np.ones(len(X)), X]); P = np.eye(A.shape[1]) * lam; P[0, 0] = 0
    return np.linalg.solve(A.T @ A + P, A.T @ y)


def ridge_pred(b, X): return np.column_stack([np.ones(len(X)), X]) @ b


def ridge_cv(X, y, rng):
    """Inner 5-fold CV over LAMBDAS (features standardized on the training data)."""
    mu, sd = X.mean(0), X.std(0) + 1e-12; Xs = (X - mu) / sd
    f = rng.permutation(len(y)) % 5; err = []
    for lam in LAMBDAS:
        e = 0.0
        for k in range(5):
            tr, te = f != k, f == k
            e += np.abs(ridge_pred(ridge_fit(Xs[tr], y[tr], lam), Xs[te]) - y[te]).sum()
        err.append(e)
    lam = LAMBDAS[int(np.argmin(err))]
    b = ridge_fit(Xs, y, lam)
    return lambda Xn: ridge_pred(b, (Xn - mu) / sd), lam


def load_physical(min_score=700):
    info = pd.read_csv(os.path.join(RAW, "string_info.txt.gz"), sep="\t", usecols=[0, 1]); name = dict(zip(info.iloc[:, 0], info.iloc[:, 1]))
    L = pd.read_csv(os.path.join(RAW, "string_physical_links.txt.gz"), sep=" "); L = L[L.combined_score >= min_score]
    L["a"] = L.protein1.map(name); L["b"] = L.protein2.map(name); L = L[L.a < L.b]
    G = nx.Graph(); G.add_weighted_edges_from(zip(L.a, L.b, L.combined_score / 1000.0)); return G


def main():
    e, m = aging.load_fleischer(); sig_f = aging.signature(aging.age_tstats(e, m))
    sig_rm = sig_f.t[[g for g in sig_f.index if g not in T.PROLIF_REMOVE]]
    ctx = {c: perturb.load_context(c) for c in CTX}
    ident = T.identity_genes(ctx); progs = {c: T.program_shifts(v[0], ident[c]) for c, v in ctx.items()}
    rm = {c: T.score(ctx[c][0], sig_rm) for c in CTX}; raw = {c: T.score(ctx[c][0], sig_f.t) for c in CTX}
    radj_full = {c: T.residualise(rm[c], progs[c][["prolif", "identity", "magnitude"]]) for c in CTX}
    G = T.load_string(); GP = load_physical(); log("networks", G.number_of_edges(), GP.number_of_edges())
    deg = dict(G.degree())
    ceg = set(pd.read_csv(os.path.join(RAW, "hart_CEGv2.txt"), sep="\t").iloc[:, 0])
    rel = pd.read_csv(os.path.join(RES, "rev3_reliability_summary.csv")) if os.path.exists(os.path.join(RES, "rev3_reliability_summary.csv")) else None
    rows, seeds_rows, blocks = [], [], []
    for c in CTX:
        shift = ctx[c][0]; P = progs[c]
        other = [o for o in CTX if LINE[o] != LINE[c]]
        transfer = pd.concat([radj_full[o] for o in other], axis=1).mean(axis=1)
        nonsig = [g for g in shift.columns if g not in sig_f.index]
        resp = shift[nonsig]; resp = resp.sub(resp.mean(axis=1), axis=0); resp = resp.div(np.sqrt((resp ** 2).sum(axis=1)) + 1e-12, axis=0)
        ontarget = pd.Series({g: shift.at[g, g] for g in shift.index if g in shift.columns})
        gmean = ctx[c][2]
        cov = pd.DataFrame({"expr": np.log1p(gmean.reindex(shift.index).fillna(0)), "core_essential": shift.index.isin(ceg).astype(float),
                            "ontarget": ontarget.reindex(shift.index).fillna(0), "magnitude": P.magnitude, "degree": pd.Series(deg).reindex(shift.index).fillna(0)},
                           index=shift.index)
        common = [p for p in shift.index if p in G and pd.notna(transfer.get(p, np.nan))]
        common = np.array(common)
        # fair kNN (Amendment 4): the held-out gene's response profile in the OTHER cell lines (never its profile in this screen)
        ors = []
        for o in other:
            so = ctx[o][0]; cols = [g for g in nonsig if g in so.columns]
            ro = so[cols].reindex(common); ro = ro.sub(ro.mean(axis=1), axis=0); ors.append(ro)
        cols_o = sorted(set.intersection(*[set(r_.columns) for r_ in ors]))
        OP = pd.concat([r_[cols_o] for r_ in ors]).groupby(level=0).mean().reindex(common)
        OP = OP.div(np.sqrt((OP ** 2).sum(axis=1)) + 1e-12, axis=0).fillna(0.0)
        RT = resp.reindex(index=common, columns=cols_o).fillna(0.0); RT = RT.sub(RT.mean(axis=1), axis=0)
        RT = RT.div(np.sqrt((RT ** 2).sum(axis=1)) + 1e-12, axis=0)
        # response clusters for the block bootstrap
        from sklearn.cluster import KMeans
        cl = pd.Series(KMeans(100, n_init=3, random_state=0).fit_predict(resp.loc[common].values), index=common)
        for seed in range(NSEED):
            rng = np.random.default_rng(seed); folds = rng.permutation(len(common)) % 5
            pred = {k: pd.Series(np.nan, index=common) for k in ["constant", "transfer", "transfer_calibrated", "network", "network_calibrated",
                                                                   "physical_network", "transfer_plus_network", "knn_response", "knn_other_lines", "gene_covariates"]}
            ytrue = pd.Series(np.nan, index=common); praw = pd.Series(np.nan, index=common); rawtrue = raw[c].reindex(common)
            lams = []
            for k in range(5):
                tr, te = common[folds != k], common[folds == k]
                Xtr = P.loc[tr, ["prolif", "identity", "magnitude"]].values; Xte = P.loc[te, ["prolif", "identity", "magnitude"]].values
                bres = np.linalg.lstsq(np.column_stack([np.ones(len(tr)), Xtr]), rm[c][tr].values, rcond=None)[0]
                ytr = rm[c][tr].values - np.column_stack([np.ones(len(tr)), Xtr]) @ bres
                yte = rm[c][te].values - np.column_stack([np.ones(len(te)), Xte]) @ bres
                ytrue[te] = yte; ytr_s = pd.Series(ytr, index=tr)
                pred["constant"][te] = ytr.mean()
                pred["transfer"][te] = transfer[te].values
                f_, lam = ridge_cv(transfer[tr].values[:, None], ytr, rng); pred["transfer_calibrated"][te] = f_(transfer[te].values[:, None]); lams.append(lam)
                nb_tr = T.neighbour_mean(G, ytr_s.to_dict(), tr).reindex(tr); nb_te = T.neighbour_mean(G, ytr_s.to_dict(), te).reindex(te)
                nb_tr_f, nb_te_f = nb_tr.fillna(ytr.mean()), nb_te.fillna(ytr.mean())
                pred["network"][te] = nb_te_f.values
                f_, _ = ridge_cv(nb_tr_f.values[:, None], ytr, rng); pred["network_calibrated"][te] = f_(nb_te_f.values[:, None])
                ph_te = T.neighbour_mean(GP, ytr_s.to_dict(), te).reindex(te).fillna(ytr.mean()); pred["physical_network"][te] = ph_te.values
                Xa = np.column_stack([transfer[tr].values, nb_tr_f.values]); Xb = np.column_stack([transfer[te].values, nb_te_f.values])
                f_, _ = ridge_cv(Xa, ytr, rng); pred["transfer_plus_network"][te] = f_(Xb)
                sim = resp.loc[te].values @ resp.loc[tr].values.T
                nn = np.argsort(-sim, axis=1)[:, :KNN]; pred["knn_response"][te] = ytr[nn].mean(axis=1)
                sim_o = OP.loc[te].values @ RT.loc[tr].values.T
                nn_o = np.argsort(-sim_o, axis=1)[:, :KNN]; pred["knn_other_lines"][te] = ytr[nn_o].mean(axis=1)
                f_, _ = ridge_cv(cov.loc[tr].values, ytr, rng); pred["gene_covariates"][te] = f_(cov.loc[te].values)
                f_, _ = ridge_cv(P.loc[tr, ["prolif"]].values, rawtrue[tr].values, rng); praw[te] = f_(P.loc[te, ["prolif"]].values)
            for nm, p in pred.items():
                err = (p - ytrue).abs(); sse = ((p - ytrue) ** 2).sum(); sst = ((ytrue - ytrue.mean()) ** 2).sum()
                seeds_rows.append({"screen": c, "seed": seed, "predictor": nm, "n": len(common), "mae": err.mean(), "rho": stats.spearmanr(p, ytrue)[0], "r2": 1 - sse / sst})
            sse = ((praw - rawtrue) ** 2).sum(); sst = ((rawtrue - rawtrue.mean()) ** 2).sum()
            seeds_rows.append({"screen": c, "seed": seed, "predictor": "prolif_predicts_raw_score", "n": len(common), "mae": (praw - rawtrue).abs().mean(),
                               "rho": stats.spearmanr(praw, rawtrue)[0], "r2": 1 - sse / sst})
            if seed == 0:
                d = ((pred["transfer_plus_network"] - ytrue).abs() - (pred["transfer_calibrated"] - ytrue).abs())
                base = (pred["transfer_calibrated"] - ytrue).abs().mean()
                ucl = cl.unique(); r2 = np.random.default_rng(99); bs = []
                for _ in range(NBLOCK):
                    pick = r2.choice(ucl, len(ucl)); idx = np.concatenate([cl.index[cl == u] for u in pick]); bs.append(d[idx].mean())
                blocks.append({"screen": c, "n": len(common), "n_clusters": len(ucl), "delta_mae": d.mean(), "rel_change_pct": 100 * d.mean() / base,
                               "block_lo": float(np.quantile(bs, 0.025)), "block_hi": float(np.quantile(bs, 0.975)),
                               "lambda_transfer_median": float(np.median(lams)), "frac_lambda_at_lower_edge": float(np.mean(np.isclose(lams, LAMBDAS[0]))),
                               "frac_lambda_at_upper_edge": float(np.mean(np.isclose(lams, LAMBDAS[-1])))})
        log(c, "done", len(common))
    S = pd.DataFrame(seeds_rows); S.to_csv(os.path.join(RES, "rev3_network_seeds.csv"), index=False)
    A = S.groupby(["screen", "predictor"]).agg(n=("n", "first"), mae=("mae", "mean"), mae_sd=("mae", "std"), rho=("rho", "mean"), rho_sd=("rho", "std"),
                                                r2=("r2", "mean"), r2_sd=("r2", "std")).reset_index()
    if rel is not None:
        rmap = {(r.screen, r.variant): r.reliability_SB for r in rel.itertuples()}
        A["reliability_ceiling_r2"] = [rmap.get((s.replace("K562gw", "K562"), "V1 pre-specified"), np.nan) for s in A.screen]
    A.to_csv(os.path.join(RES, "rev3_network_summary.csv"), index=False)
    B = pd.DataFrame(blocks); B.to_csv(os.path.join(RES, "rev3_network_block_bootstrap.csv"), index=False)
    out = {"n_seeds": NSEED, "knn_k": KNN, "lambdas": LAMBDAS.tolist(), "block_bootstrap_resamples": NBLOCK,
           "n_screens_block_ci_below_zero": int((B.block_hi < 0).sum()), "common_n": B.set_index("screen").n.to_dict(),
           "physical_edges": GP.number_of_edges()}
    json.dump(out, open(os.path.join(RES, "rev3_network.json"), "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    log("done", out)


if __name__ == "__main__":
    main()
