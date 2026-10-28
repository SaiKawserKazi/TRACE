"""Exploratory (Amendment 4, X23): the full run with a GTEx-fibroblast-derived signature as discovery input, reported
next to the discovery signature.
* H1-H5 side by side (from results/results.json and results/e18_gtex_discovery/results.json)
* signature composition: cell-cycle share per arm, top Hallmark sets per arm
* held-out labels for the GTEx signature: R_rob (V1), mean raw score, proliferation shift alone; all genes and genes in
  >= 2 lines; paired DeLong raw vs proliferation; incremental logistic model
* program tests (same logic as H2): Spearman of each signature's raw score with proliferation, stress, oxidative
  phosphorylation, response magnitude and identity shifts, per screen, Holm over the five screens
Outputs: results/x_gtex_signature_*.csv, results/x_gtex_signature.json
"""
import os, sys, json
import numpy as np, pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
import aging, perturb
import trace_analysis as T
import rev3_scores as S

RES = T.RES; E18 = os.path.join(RES, "e18_gtex_discovery")
CTX = ["K562", "K562gw", "RPE1", "HepG2", "Jurkat"]


def holm(p):
    p = np.asarray(p); o = np.argsort(p); m = len(p); out = np.empty(m); run = 0
    for i, k in enumerate(o): run = max(run, (m - i) * p[k]); out[k] = min(1, run)
    return out


def main():
    R0, R1 = json.load(open(os.path.join(RES, "results.json"))), json.load(open(os.path.join(E18, "results.json")))
    rows = []
    for lab, R in [("Discovery signature", R0), ("GTEx-fibroblast signature", R1)]:
        h1 = R["H1"]["pairs"]; h2 = R["H2"]["contexts"]; h3 = R["H3"]["contexts"]
        rows += [{"signature": lab, "test": "H1a raw cross-line rho", "value": f"{min(v['rho_raw'] for v in h1.values()):.2f}-{max(v['rho_raw'] for v in h1.values()):.2f}", "supported": R["H1"]["supported"]},
                 {"signature": lab, "test": "H1b within-K562 vs between-line", "value": f"{R['H1b']['rho_within_K562']:.2f} vs {R['H1b']['rho_between_mean']:.2f}", "supported": R["H1b"]["supported"]},
                 {"signature": lab, "test": "adjusted cross-line rho", "value": f"{min(v['rho_adj'] for v in h1.values()):.2f}-{max(v['rho_adj'] for v in h1.values()):.2f}", "supported": None},
                 {"signature": lab, "test": "H2 rho(R, proliferation)", "value": f"{min(v['rho_raw_prolif'] for v in h2.values()):.2f} to {max(v['rho_raw_prolif'] for v in h2.values()):.2f}", "supported": R["H2"]["supported"]},
                 {"signature": lab, "test": "H3a network vs 10 nulls (screens significant)", "value": f"{sum(1 for c in h3 if h3[c]['h3a_p_holm'] < 0.05)}/5", "supported": R["H3"]["H3a_supported"]},
                 {"signature": lab, "test": "H3b transfer + network vs line transfer", "value": f"{sum(1 for c in h3 if h3[c]['h3b_p_holm'] < 0.05)}/5", "supported": R["H3"]["H3b_supported"]},
                 {"signature": lab, "test": "H4a CellAge AUROC (Holm P)", "value": f"{R['H4']['H4a_cellage']['auroc']:.3f} ({R['H4']['H4a_cellage']['p_holm']:.4f})", "supported": R["H4"]["supported_H4a"]},
                 {"signature": lab, "test": "H4b GenAge AUROC (Holm P)", "value": f"{R['H4']['H4b_genage']['auroc']:.3f} ({R['H4']['H4b_genage']['p_holm']:.4f})", "supported": R["H4"]["supported_H4b"]},
                 {"signature": lab, "test": "H5a donor-bootstrap median rho (threshold 0.8)", "value": f"{R['H5']['median_rho']:.3f}", "supported": R["H5"]["supported_stability"]}]
        for k, v in R["H5"].items():
            if isinstance(v, dict) and "rho" in v:
                rows.append({"signature": lab, "test": f"H5b ranking vs {k.replace('_', ' ')} signature", "value": f"{v['rho']:.2f} (top-5% kept {100 * v['top5_overlap']:.0f}%)",
                             "supported": (v["p"] < 0.05) if k == "gtex_fibroblast" else None})
    pd.DataFrame(rows).to_csv(os.path.join(RES, "x_gtex_signature_h1_h5.csv"), index=False)

    # signature composition
    HM = S.T.HM if hasattr(S.T, "HM") else T.HM
    cc = HM["E2F Targets"] | HM["G2-M Checkpoint"] | HM["Mitotic Spindle"]
    e, m = aging.load_gtex("gtex_fibroblasts_tpm.gct.gz"); sig_g = aging.signature(aging.age_tstats(e, m))
    e0, m0 = aging.load_fleischer(); sig_d = aging.signature(aging.age_tstats(e0, m0))
    comp = {}
    for lab, sg in [("discovery", sig_d), ("gtex_fibroblast", sig_g)]:
        comp[lab] = {"n_up": int((sg.t > 0).sum()), "n_down": int((sg.t < 0).sum()), "frac_cellcycle_all": float(np.mean([g in cc for g in sg.index])),
                     "frac_cellcycle_up": float(np.mean([g in cc for g in sg[sg.t > 0].index])), "frac_cellcycle_down": float(np.mean([g in cc for g in sg[sg.t < 0].index]))}
    EN = pd.read_csv(os.path.join(E18, "signature_hallmark_enrichment.csv"))
    comp["gtex_top_sets"] = {d: EN[EN.dir == d].nsmallest(5, "p")[["set", "overlap", "size", "p", "q"]].to_dict("records") for d in ("up", "down")}

    # program tests (H2 logic) for both signatures
    ctx = {c: perturb.load_context(c) for c in CTX}
    ident = T.identity_genes(ctx); progs = {c: T.program_shifts(v[0], ident[c]) for c, v in ctx.items()}
    ox = HM["Oxidative Phosphorylation"]
    prog_rows = []
    for lab, folder in [("Discovery signature", RES), ("GTEx-fibroblast signature", E18)]:
        for prog in ["prolif", "stress", "oxphos", "magnitude", "identity"]:
            res = []
            for c in CTX:
                sc = pd.read_csv(os.path.join(folder, f"scores_{c}.csv"), index_col=0)
                if prog == "oxphos":
                    sh = ctx[c][0]; x = sh[[g for g in ox if g in sh.columns]].mean(axis=1)
                else:
                    x = progs[c][prog]
                common = sc.index.intersection(x.index); rho, p = stats.spearmanr(sc.R_raw[common], x[common]); res.append((c, rho, p, len(common)))
            ph = holm([r[2] for r in res])
            prog_rows += [{"signature": lab, "program": prog, "screen": c, "rho": rho, "p": p, "p_holm": h, "n": n} for (c, rho, p, n), h in zip(res, ph)]
    PR = pd.DataFrame(prog_rows); PR.to_csv(os.path.join(RES, "x_gtex_signature_programs.csv"), index=False)

    # held-out labels for the GTEx signature
    labs = S.labels()
    rob = pd.read_csv(os.path.join(E18, "robust_scores.csv"), index_col=0)
    sc = {c: pd.read_csv(os.path.join(E18, f"scores_{c}.csv"), index_col=0) for c in CTX}
    raw_mean = pd.concat([sc[c].R_raw for c in CTX], axis=1).mean(axis=1)
    pro_mean = pd.concat([progs[c].prolif for c in CTX], axis=1).mean(axis=1)
    multi = rob.index[rob.n_lines >= 2]; ev, paired = [], []
    for nm, (lab, pos, neg) in labs.items():
        for gs_name, gs in [("all genes", rob.index), ("genes in >= 2 lines", multi)]:
            for sname, s in [("R_rob (V1)", rob.R_rob), ("mean raw score", raw_mean), ("proliferation shift alone", pro_mean)]:
                x, y = S.evid(s.reindex(gs), lab, pos, neg)
                if len(x) >= 5 and len(y) >= 5:
                    ev.append({"labels": nm, "gene_set": gs_name, "score": sname, "n_pos": len(x), "n_neg": len(y), "auroc": S.auc(x, y),
                               "min_detectable_auroc_80pct": S.min_detectable_auc(len(x), len(y))})
            both = gs.intersection(raw_mean.dropna().index).intersection(pro_mean.dropna().index)
            for a_name, a_s in [("mean raw score", raw_mean), ("R_rob (V1)", rob.R_rob)]:
                bb = both.intersection(a_s.dropna().index)
                xa, ya = S.evid(a_s.reindex(bb), lab, pos, neg); xb, yb = S.evid(pro_mean.reindex(bb), lab, pos, neg)
                if len(xa) >= 5 and len(ya) >= 5:
                    aA, aB, d, se, lo, hi, p = S.delong(xa.values, ya.values, xb.values, yb.values)
                    paired.append({"labels": nm, "gene_set": gs_name, "a": a_name, "b": "proliferation shift alone", "n_pos": len(xa), "n_neg": len(ya),
                                   "auroc_a": aA, "auroc_b": aB, "diff_a_minus_b": d, "lo": lo, "hi": hi, "p_two_sided": p})
    pd.DataFrame(ev).to_csv(os.path.join(RES, "x_gtex_signature_evidence.csv"), index=False)
    pd.DataFrame(paired).to_csv(os.path.join(RES, "x_gtex_signature_paired_auroc.csv"), index=False)
    json.dump({"composition": comp}, open(os.path.join(RES, "x_gtex_signature.json"), "w"), indent=1, default=float)
    print(pd.DataFrame(rows).to_string()); print(json.dumps(comp, indent=1, default=float)[:1500])
    print(PR.pivot_table(index=["signature", "program"], columns="screen", values="rho").round(2).to_string())
    print(pd.DataFrame(ev)[lambda d: d.labels.isin(["CellAge", "GenAge"])].round(3).to_string()); print(pd.DataFrame(paired).round(3).to_string())


if __name__ == "__main__":
    main()
