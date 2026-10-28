"""Exploratory (Amendment 4, X25): effect of the V6 (spline) and V7 (no magnitude) definitions on the top-5% list and the
Pareto front, next to V1, on genes measured in >= 2 lines. Output: results/x_v6v7_effects.json"""
import os, sys, json
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import aging, perturb
import trace_analysis as T
import rev3_scores as S

RES = T.RES; CTX = ["K562", "K562gw", "RPE1", "HepG2", "Jurkat"]


def main():
    e, m = aging.load_fleischer(); sig = aging.signature(aging.age_tstats(e, m))
    ctx = {c: perturb.load_context(c) for c in CTX}; ident = T.identity_genes(ctx); progs = {c: T.program_shifts(v[0], ident[c]) for c, v in ctx.items()}
    Vd = S.variants(ctx, sig.t, progs)
    P0 = pd.read_csv(os.path.join(RES, "rev3_scores_pareto_filtered.csv"), index_col=0)
    objs = ["R_rob", "sign_consistency", "obj_identity", "obj_prolif", "obj_stress"]
    out, fronts, tops = {}, {}, {}
    for v in ["V1 pre-specified", "V6 spline resid.", "V7 resid. (prolif. + identity)"]:
        rb = S.robust_of(Vd[v]); multi = rb.index[rb.n_lines >= 2].intersection(P0.index)
        P = P0.loc[multi].copy(); P["R_rob"] = rb.R_rob.reindex(multi); P["sign_consistency"] = rb.sign_consistency.reindex(multi)
        fronts[v] = set(T.pareto(P[P.R_rob > 0], objs)); k = max(1, int(0.05 * len(multi))); tops[v] = set(P.R_rob.nlargest(k).index)
        out[v] = {"n_genes": int(len(multi)), "n_front": len(fronts[v]), "front": sorted(fronts[v]), "k_top5": k}
    for v in ["V6 spline resid.", "V7 resid. (prolif. + identity)"]:
        out[v]["front_shared_with_V1"] = len(fronts[v] & fronts["V1 pre-specified"]); out[v]["top5_shared_with_V1"] = len(tops[v] & tops["V1 pre-specified"])
    json.dump(out, open(os.path.join(RES, "x_v6v7_effects.json"), "w"), indent=1)
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk != "front"} for k, v in out.items()}, indent=1))


if __name__ == "__main__":
    main()
