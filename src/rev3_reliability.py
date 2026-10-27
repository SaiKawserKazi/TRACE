"""v0.1 revision (exploratory, added in revision; AMENDMENT-R3), E6 and part of E13: score reliability.

* Split-cell reliability for HepG2 and Jurkat (raw single-cell data in hand): for 20 random splits, the cells of each
  knockdown are split in two halves; pseudobulk shifts, program shifts and scores V0..V7 are recomputed per half and
  correlated between halves (Spearman), then Spearman-Brown corrected to full-sample reliability.
  Program shifts and magnitude are computed on the scored gene subset (signature, proliferation, stress, identity genes).
* Size-matched non-targeting pseudo-replicates (cell count = median knockdown cell count) for Figure S1.
* Held-out proliferation (E13): correlation of R from half B with proliferation from half A versus from half B.
* K562 test-retest reliability: Spearman of each definition between the essential-gene and genome-wide K562 screens.
* Disattenuated between-line agreement for pairs whose reliabilities are known (K562, HepG2, Jurkat).
Outputs: results/rev3_reliability.json, results/rev3_reliability_*.csv
"""
import os, sys, json, time
import numpy as np, pandas as pd, h5py
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
import aging, perturb
import trace_analysis as T
from rev3_scores import spline_basis

RAW, RES = T.RAW, T.RES
NSPLIT = int(os.environ.get("TRACE_NSPLIT", 20)); BLOCK = 4000
VAR = ["V0 raw", "V1 pre-specified", "V2 genes removed", "V3 resid. (all 3)", "V4 resid. (prolif.)", "V5 stratified",
       "V6 spline resid.", "V7 resid. (prolif. + identity)"]


def log(*a): print(time.strftime("%H:%M:%S"), *a, flush=True)


def read_col(g, k):
    x = g[k]
    if "__categories" in g and k in g["__categories"]:
        return g["__categories"][k][:].astype(str)[x[:]]
    if isinstance(x, h5py.Group): return x["categories"][:].astype(str)[x["codes"][:]]
    v = x[:]; return v.astype(str) if v.dtype.kind in "SO" else v


def scores_from(shift, sig_t, ident):
    """V0..V7 for one pseudobulk matrix (genes restricted to the scored subset)."""
    P = pd.DataFrame({"prolif": shift[[g for g in T.PROLIF if g in shift.columns]].mean(axis=1),
                      "stress": shift[[g for g in T.STRESS if g in shift.columns]].mean(axis=1),
                      "identity": shift[[g for g in ident if g in shift.columns]].mean(axis=1),
                      "magnitude": shift.abs().mean(axis=1)})
    sig_rm = sig_t[[g for g in sig_t.index if g not in T.PROLIF_REMOVE]]
    raw = T.score(shift, sig_t); rm = T.score(shift, sig_rm)
    q = pd.qcut(P.prolif.rank(method="first"), 5, labels=False)
    out = {"V0 raw": raw, "V1 pre-specified": T.residualise(rm, P[["prolif", "identity", "magnitude"]]), "V2 genes removed": rm,
           "V3 resid. (all 3)": T.residualise(raw, P[["prolif", "identity", "magnitude"]]), "V4 resid. (prolif.)": T.residualise(raw, P[["prolif"]]),
           "V5 stratified": raw.groupby(q.values).rank(pct=True) - 0.5,
           "V6 spline resid.": T.residualise(rm, pd.concat([spline_basis(P.prolif), P[["identity", "magnitude"]]], axis=1)),
           "V7 resid. (prolif. + identity)": T.residualise(rm, P[["prolif", "identity"]])}
    return out, P


def split_line(name, sig_t, ident, keep_genes, rng, fname=None, normalized=False):
    """normalized=True reads cells already z-scored against non-targeting cells of their gem group (Replogle release),
    so the per-gem-group standardization step is skipped; everything after it is identical."""
    f = h5py.File(os.path.join(RAW, fname or f"{name}_raw_singlecell.h5ad"), "r")
    obs, var, X = f["obs"], f["var"], f["X"]
    tgt = read_col(obs, "gene").astype(str); gem = read_col(obs, "gem_group").astype(str)
    umi = read_col(obs, "UMI_count").astype(float) if "UMI_count" in obs else np.ones(len(tgt))
    gname = read_col(var, "gene_name") if "gene_name" in var else var[var.attrs.get("_index", "_index")][:].astype(str)
    gsel = np.array([i for i, g in enumerate(gname) if g in keep_genes]); genes = gname[gsel]
    n = len(tgt); is_nt = tgt == "non-targeting"; med = np.median(umi[is_nt])
    gems = sorted(set(gem)); gi = {g: i for i, g in enumerate(gems)}; G = len(gsel)
    S1 = np.zeros((len(gems), G)); S2 = np.zeros((len(gems), G)); N = np.zeros(len(gems))
    def blocks():
        for s in range(0, n, BLOCK):
            e = min(n, s + BLOCK)
            D = X[s:e, :][:, gsel].astype(np.float64)
            if not normalized: D = D / np.maximum(umi[s:e, None], 1) * med
            yield s, e, D
    if not normalized:
        for s, e, D in blocks():
            for j in np.where(is_nt[s:e])[0]:
                k = gi[gem[s + j]]; S1[k] += D[j]; S2[k] += D[j] ** 2; N[k] += 1
        mu = S1 / np.maximum(N[:, None], 1); sd = np.sqrt(np.maximum(S2 / np.maximum(N[:, None], 1) - mu ** 2, 0)); sd[sd == 0] = np.nan
    else:
        mu = np.zeros((len(gems), G)); sd = np.ones((len(gems), G))
    perts = sorted(set(tgt[~is_nt])); pi = {p: i for i, p in enumerate(perts)}
    counts = pd.Series(tgt[~is_nt]).value_counts()
    medcells = int(counts[counts >= 20].median())
    # split assignment per cell for NSPLIT splits (balanced within each knockdown)
    half = np.zeros((NSPLIT, n), dtype=np.int8)
    idx_by = pd.Series(np.arange(n)).groupby(tgt).apply(list)
    for p, ids in idx_by.items():
        ids = np.array(ids)
        for s_ in range(NSPLIT):
            perm = rng.permutation(len(ids)); half[s_, ids[perm[: len(ids) // 2]]] = 1
    # size-matched non-targeting pseudo-replicates
    nt_ids = np.where(is_nt)[0]; NNT = 200
    nt_groups = [rng.choice(nt_ids, medcells, replace=False) for _ in range(NNT)]
    nt_member = {}
    for gidx, ids in enumerate(nt_groups):
        for c in ids: nt_member.setdefault(c, []).append(gidx)
    Z = np.zeros((NSPLIT, 2, len(perts), G), dtype=np.float32); C = np.zeros((NSPLIT, 2, len(perts)), dtype=np.int32)
    ZF = np.zeros((len(perts), G)); CF = np.zeros(len(perts)); ZN = np.zeros((NNT, G)); CN = np.zeros(NNT)
    for s, e, D in blocks():
        k = np.array([gi[g] for g in gem[s:e]])
        Zb = np.nan_to_num((D - mu[k]) / sd[k], nan=0.0)
        for j in range(e - s):
            c_ = s + j; t = tgt[c_]
            if t == "non-targeting":
                for gidx in nt_member.get(c_, []): ZN[gidx] += Zb[j]; CN[gidx] += 1
                continue
            p = pi[t]; ZF[p] += Zb[j]; CF[p] += 1
            for s_ in range(NSPLIT):
                h = half[s_, c_]; Z[s_, h, p] += Zb[j]; C[s_, h, p] += 1
        if s % 40000 == 0: log(name, "cells", s)
    keep = CF >= 20
    P_index = np.array(perts)[keep]
    full = pd.DataFrame(ZF[keep] / CF[keep, None], index=P_index, columns=genes)
    full = full.T.groupby(level=0).mean().T
    nt = pd.DataFrame(ZN / CN[:, None], columns=genes).T.groupby(level=0).mean().T
    rel, heldout = [], []
    for s_ in range(NSPLIT):
        halves = []
        for h in (0, 1):
            m = pd.DataFrame(Z[s_, h][keep] / np.maximum(C[s_, h][keep], 1)[:, None], index=P_index, columns=genes).T.groupby(level=0).mean().T
            halves.append(scores_from(m, sig_t, ident))
        for v in VAR:
            a, b = halves[0][0][v], halves[1][0][v]
            rel.append({"screen": name, "split": s_, "variant": v, "r_half": stats.spearmanr(a, b.reindex(a.index))[0]})
        rB = halves[1][0]["V0 raw"]
        heldout.append({"screen": name, "split": s_, "rho_R_B_prolif_A": stats.spearmanr(rB, halves[0][1].prolif.reindex(rB.index))[0],
                        "rho_R_B_prolif_B": stats.spearmanr(rB, halves[1][1].prolif.reindex(rB.index))[0]})
    full_scores, full_P = scores_from(full, sig_t, ident)
    nt_scores = T.score(nt, sig_t)
    sizes = {"median_knockdown_cells": medcells, "nt_groups": NNT, "nt_sd_R": float(nt_scores.std()), "kd_sd_R": float(full_scores["V0 raw"].std()),
             "nt_q95_R": float(nt_scores.quantile(0.95))}
    return pd.DataFrame(rel), pd.DataFrame(heldout), sizes, nt_scores


def main():
    e, m = aging.load_fleischer(); sig_f = aging.signature(aging.age_tstats(e, m))
    ctx = {c: perturb.load_context(c) for c in perturb.CONTEXTS}
    ident = T.identity_genes(ctx)
    rng = np.random.default_rng(70_000)
    out, rels, helds, ntd = {}, [], [], {}
    for name, c in [("hepg2", "HepG2"), ("jurkat", "Jurkat")]:
        keep = set(sig_f.index) | T.PROLIF_REMOVE | T.PROLIF | T.STRESS | set(ident[c])
        r, h, sizes, nts = split_line(name, sig_f.t, ident[c], keep, rng)
        rels.append(r); helds.append(h); out[c] = sizes; ntd[c] = nts
        log(name, "done", sizes)
    R = pd.concat(rels); H = pd.concat(helds)
    R.to_csv(os.path.join(RES, "rev3_reliability_splits.csv"), index=False); H.to_csv(os.path.join(RES, "rev3_reliability_heldout_prolif.csv"), index=False)
    pd.concat({k: v for k, v in ntd.items()}).rename("R_raw").to_csv(os.path.join(RES, "rev3_reliability_nt_sizematched.csv"))
    agg = R.groupby(["screen", "variant"]).r_half.agg(["mean", "std"]).reset_index()
    agg["reliability_SB"] = 2 * agg["mean"] / (1 + agg["mean"])
    # K562 test-retest per definition, and cross-line agreement per definition
    from rev3_scores import variants
    progs = {c: T.program_shifts(v[0], ident[c]) for c, v in ctx.items()}
    Vd = variants(ctx, sig_f.t, progs)
    tr, dis = [], []
    for v in VAR:
        a, b = Vd[v]["K562"], Vd[v]["K562gw"]; sh = a.index.intersection(b.index)
        tr.append({"screen": "k562", "variant": v, "mean": stats.spearmanr(a[sh], b[sh])[0], "std": np.nan,
                   "reliability_SB": stats.spearmanr(a[sh], b[sh])[0]})
    REL = pd.concat([agg, pd.DataFrame(tr)], ignore_index=True)
    REL["screen"] = REL.screen.map({"hepg2": "HepG2", "jurkat": "Jurkat", "k562": "K562"})
    REL.to_csv(os.path.join(RES, "rev3_reliability_summary.csv"), index=False)
    relmap = {(r.screen, r.variant): r.reliability_SB for r in REL.itertuples()}
    for v in VAR:
        for a_, b_ in [("K562", "HepG2"), ("K562", "Jurkat"), ("HepG2", "Jurkat")]:
            sa, sb = Vd[v][a_], Vd[v][b_]; sh = sa.index.intersection(sb.index); rho = stats.spearmanr(sa[sh], sb[sh])[0]
            ra, rb = relmap[(a_, v)], relmap[(b_, v)]
            dis.append({"variant": v, "pair": f"{a_}-{b_}", "n": len(sh), "rho": rho, "rel_a": ra, "rel_b": rb,
                        "rho_disattenuated": rho / np.sqrt(ra * rb) if ra > 0 and rb > 0 else np.nan})
    DIS = pd.DataFrame(dis); DIS.to_csv(os.path.join(RES, "rev3_reliability_disattenuated.csv"), index=False)
    out["heldout_prolif"] = H.groupby("screen")[["rho_R_B_prolif_A", "rho_R_B_prolif_B"]].mean().round(3).to_dict("index")
    out["n_splits"] = NSPLIT
    json.dump(out, open(os.path.join(RES, "rev3_reliability.json"), "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    log("done")


if __name__ == "__main__":
    main()
