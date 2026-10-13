"""Review round (AMENDMENT-R1), item P0.2: alternative score definitions (exploratory).

V0 raw            R = -Spearman(d[S], a[S])
V1 protocol       proliferation genes removed from S, then residualised on proliferation, identity and magnitude shifts
V2 removal only   proliferation genes removed from S, no residualisation
V3 resid only     full S, residualised on proliferation, identity and magnitude
V4 prolif resid   full S, residualised on proliferation shift only
V5 stratified     full S; percentile rank of R within quintiles of proliferation shift (per screen), centred at 0
For each variant: association with proliferation, top-k overlap, cross-line agreement, evidence AUROCs with CIs,
and the cohort-swap test. Outputs results/rev_scores.json and CSV tables.
"""
import os, sys, json, time
import numpy as np, pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
import aging, perturb
import trace_analysis as T

RES = T.RES
rng = np.random.default_rng(2)
VARS = ["V0 raw", "V1 protocol", "V2 removal only", "V3 resid only", "V4 prolif resid", "V5 stratified"]


def log(*a): print(time.strftime("%H:%M:%S"), *a, flush=True)


def variants(ctx, sig_t, progs):
    sig_rm = sig_t[[g for g in sig_t.index if g not in T.PROLIF_REMOVE]]
    out = {v: {} for v in VARS}
    for c, (shift, nt, gm, nc) in ctx.items():
        P = progs[c]
        raw = T.score(shift, sig_t); rm = T.score(shift, sig_rm)
        cov = P[["prolif", "identity", "magnitude"]]
        out["V0 raw"][c] = raw
        out["V1 protocol"][c] = T.residualise(rm, cov)
        out["V2 removal only"][c] = rm
        out["V3 resid only"][c] = T.residualise(raw, cov)
        out["V4 prolif resid"][c] = T.residualise(raw, P[["prolif"]])
        q = pd.qcut(P.prolif.rank(method="first"), 5, labels=False)
        out["V5 stratified"][c] = raw.groupby(q.values).rank(pct=True) - 0.5
    return out


def robust_of(vd):
    return T.robust({c: pd.DataFrame({"R_adj": s}) for c, s in vd.items()})


def topk(s, frac=0.05):
    s = s.dropna(); k = max(1, int(frac * len(s)))
    return set(s.nlargest(k).index)


def ov(a, b, frac=0.05):
    sh = a.dropna().index.intersection(b.dropna().index)
    A, B = topk(a[sh], frac), topk(b[sh], frac)
    return len(A & B) / max(1, len(A))


def auc_ci(x, y, nb=1000):
    x, y = np.asarray(x), np.asarray(y)
    if len(x) < 3 or len(y) < 3: return [np.nan] * 4
    u = stats.mannwhitneyu(x, y, alternative="greater")
    a = u.statistic / (len(x) * len(y))
    bs = []
    for _ in range(nb):
        xb, yb = rng.choice(x, len(x)), rng.choice(y, len(y))
        bs.append(stats.mannwhitneyu(xb, yb).statistic / (len(x) * len(y)))
    return [a, float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975)), u.pvalue]


def main():
    e, m = aging.load_fleischer()
    sig_f = aging.signature(aging.age_tstats(e, m))
    sigs = {"fleischer": sig_f.t}
    for key, fn in [("gtex_fibroblast", "gtex_fibroblasts_tpm.gct.gz"), ("gtex_skin", "gtex_skin_ns_tpm.gct.gz")]:
        ge, gmeta = aging.load_gtex(fn); sigs[key] = aging.signature(aging.age_tstats(ge, gmeta)).t
    ctx = {c: perturb.load_context(c) for c in perturb.CONTEXTS}
    ident = T.identity_genes(ctx)
    progs = {c: T.program_shifts(v[0], ident[c]) for c, v in ctx.items()}
    V = {k: variants(ctx, s, progs) for k, s in sigs.items()}
    log("variants computed")
    lines = ["K562", "RPE1", "HepG2", "Jurkat"]
    pairs = [(a, b) for i, a in enumerate(lines) for b in lines[i + 1:]]
    ca = T.cellage(); gl, _ = T.genage_human_labels()
    res = {"variants": VARS, "per_variant": {}}
    rob = {k: {v: robust_of(V[k][v]) for v in VARS} for k in V}
    rows_ctx = []
    for v in VARS:
        d = {}
        for c in ctx:
            s = V["fleischer"][v][c]; P = progs[c]
            rows_ctx.append({"variant": v, "screen": c, "rho_prolif": stats.spearmanr(s, P.prolif)[0],
                             "rho_raw": stats.spearmanr(s, V["fleischer"]["V0 raw"][c])[0],
                             "top5_overlap_with_raw": ov(s, V["fleischer"]["V0 raw"][c]),
                             "top5_overlap_with_V1": ov(s, V["fleischer"]["V1 protocol"][c])})
        cl = {f"{a}-{b}": stats.spearmanr(*[V["fleischer"][v][x].reindex(V["fleischer"][v][a].index.intersection(V["fleischer"][v][b].index)) for x in (a, b)])[0] for a, b in pairs}
        sh = V["fleischer"][v]["K562"].index.intersection(V["fleischer"][v]["K562gw"].index)
        d["cross_line_rho"] = cl
        d["cross_line_median"] = float(np.median(list(cl.values())))
        d["within_K562_rho"] = stats.spearmanr(V["fleischer"][v]["K562"][sh], V["fleischer"][v]["K562gw"][sh])[0]
        R = rob["fleischer"][v].R_rob
        for name, lab, pos, neg in [("cellage", ca, "Induces", "Inhibits"), ("genage", gl, "Anti", "Pro")]:
            x = R.reindex(lab[lab == pos].index).dropna(); y = R.reindex(lab[lab == neg].index).dropna()
            a, lo, hi, p = auc_ci(x, y)
            d[f"{name}_auroc"], d[f"{name}_lo"], d[f"{name}_hi"], d[f"{name}_p"] = a, lo, hi, p
            d[f"{name}_n"] = [len(x), len(y)]
        for key in ("gtex_fibroblast", "gtex_skin"):
            Rg = rob[key][v].R_rob; s2 = R.index.intersection(Rg.index)
            d[f"swap_{key}_rho"] = stats.spearmanr(R[s2], Rg[s2])[0]
            d[f"swap_{key}_top5"] = ov(R, Rg)
        res["per_variant"][v] = d
        log(v, {k: (round(x, 3) if isinstance(x, float) else x) for k, x in d.items() if not isinstance(x, dict)})
    pd.DataFrame(rows_ctx).to_csv(os.path.join(RES, "rev_scores_per_screen.csv"), index=False)
    # agreement between variants (robust scores, Fleischer signature)
    M = pd.DataFrame({v: rob["fleischer"][v].R_rob for v in VARS})
    rho_m = M.corr(method="spearman")
    ovm = pd.DataFrame([[ov(M[a], M[b]) for b in VARS] for a in VARS], index=VARS, columns=VARS)
    rho_m.to_csv(os.path.join(RES, "rev_scores_variant_rho.csv")); ovm.to_csv(os.path.join(RES, "rev_scores_variant_top5.csv"))
    res["variant_rho_matrix"] = rho_m.round(3).to_dict(); res["variant_top5_matrix"] = ovm.round(3).to_dict()
    # examples: robust raw top candidates that drop or stay after protocol adjustment (lines >= 2)
    lines_mean_prolif = pd.concat([progs[c].prolif for c in ctx], axis=1).mean(axis=1)
    E = pd.DataFrame({"R_rob_raw": M["V0 raw"], "R_rob_protocol": M["V1 protocol"], "R_rob_stratified": M["V5 stratified"],
                      "n_lines": rob["fleischer"]["V1 protocol"].n_lines, "prolif_shift_mean": lines_mean_prolif}).dropna(subset=["R_rob_raw"])
    E = E[E.n_lines >= 2]
    E["rank_raw"] = E.R_rob_raw.rank(ascending=False); E["rank_protocol"] = E.R_rob_protocol.rank(ascending=False)
    E["rank_stratified"] = E.R_rob_stratified.rank(ascending=False)
    top = E.nsmallest(40, "rank_raw")
    k = int(0.05 * len(E))
    drop = top[top.rank_protocol > k].nsmallest(8, "rank_raw"); keep = top[top.rank_protocol <= k].nsmallest(8, "rank_raw")
    E.to_csv(os.path.join(RES, "rev_scores_candidates.csv"))
    res["examples_drop"] = drop.round(4).reset_index().rename(columns={"index": "gene"}).to_dict("records")
    res["examples_keep"] = keep.round(4).reset_index().rename(columns={"index": "gene"}).to_dict("records")
    res["n_candidates_ge2"] = len(E); res["k_top5"] = k
    res["top_raw_prolif_mean"] = float(E.nsmallest(k, "rank_raw").prolif_shift_mean.mean())
    res["top_protocol_prolif_mean"] = float(E.nsmallest(k, "rank_protocol").prolif_shift_mean.mean())
    res["top_stratified_prolif_mean"] = float(E.nsmallest(k, "rank_stratified").prolif_shift_mean.mean())
    res["all_prolif_mean"] = float(E.prolif_shift_mean.mean())
    json.dump(res, open(os.path.join(RES, "rev_scores.json"), "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    log("done")


if __name__ == "__main__":
    main()
