"""TRACE main analysis. Implements RESEARCH_PROTOCOL.md v1.0 sections 3-7.

Every number reported in the manuscript is written to results/results.json.
"""
import os, json, zipfile, io, sys, time
import numpy as np, pandas as pd
from scipy import stats
import networkx as nx

sys.path.insert(0, os.path.dirname(__file__))
import aging, perturb

ROOT = r"D:\Study\RejuvBench"
RAW, PROC = (os.path.join(ROOT, d) for d in ("data/raw", "data/proc"))
# TRACE_RESULTS / TRACE_DISCOVERY let the same pipeline run on another input (E18) without code changes
RES = os.path.join(ROOT, os.environ.get("TRACE_RESULTS", "results"))
DISCOVERY = os.environ.get("TRACE_DISCOVERY", "fleischer")
os.makedirs(RES, exist_ok=True)
RNG = np.random.default_rng(0)
LINE_OF = {"K562": "K562", "K562gw": "K562", "RPE1": "RPE1", "HepG2": "HepG2", "Jurkat": "Jurkat"}
R = {}  # results


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


# ---------------------------------------------------------------- gene sets
def hallmark():
    sets = {}
    for l in open(os.path.join(RAW, "hallmark.txt"), encoding="utf-8"):
        f = l.rstrip("\n").split("\t")
        sets[f[0]] = set(g for g in f[2:] if g)
    return sets


HM = hallmark()
PROLIF_REMOVE = set().union(*[HM[k] for k in ["E2F Targets", "G2-M Checkpoint", "Mitotic Spindle",
                                               "Myc Targets V1", "Myc Targets V2"]])
PROLIF = HM["E2F Targets"] | HM["G2-M Checkpoint"]
STRESS = HM["TNF-alpha Signaling via NF-kB"] | HM["Inflammatory Response"] | HM["p53 Pathway"]


# ---------------------------------------------------------------- scoring
def rank_rows(M):
    return stats.rankdata(M, axis=1)


def score(shift, sig_t):
    """R = -Spearman(d_p[S], a[S]) with the perturbation's own target excluded."""
    genes = [g for g in sig_t.index if g in shift.columns]
    a = sig_t[genes].values
    D = shift[genes].values
    ra = stats.rankdata(a); ra = (ra - ra.mean()) / ra.std()
    rD = rank_rows(D); rD = (rD - rD.mean(1, keepdims=True)) / rD.std(1, keepdims=True)
    r = -(rD @ ra) / len(genes)
    out = pd.Series(r, index=shift.index)
    gi = {g: i for i, g in enumerate(genes)}
    for p in shift.index:
        if p in gi:
            keep = np.ones(len(genes), bool); keep[gi[p]] = False
            out[p] = -stats.spearmanr(D[shift.index.get_loc(p), keep], a[keep])[0]
    return out


def program_shifts(shift, ident_genes):
    def m(gs):
        gs = [g for g in gs if g in shift.columns]
        return shift[gs].mean(axis=1)
    return pd.DataFrame({"prolif": m(PROLIF), "stress": m(STRESS), "identity": m(ident_genes),
                         "magnitude": shift.abs().mean(axis=1)})


def residualise(y, X):
    A = np.column_stack([np.ones(len(X)), X.values])
    b, *_ = np.linalg.lstsq(A, y.values, rcond=None)
    return pd.Series(y.values - A @ b, index=y.index)


def all_scores(ctx, sig_t, progs):
    """Returns dict context -> DataFrame(R_raw, R_adj) and line-level R_adj table, R_rob."""
    sig_adj = sig_t[[g for g in sig_t.index if g not in PROLIF_REMOVE]]
    per = {}
    for c, (shift, nt, gm, nc) in ctx.items():
        raw = score(shift, sig_t)
        adj0 = score(shift, sig_adj)
        adj = residualise(adj0, progs[c][["prolif", "identity", "magnitude"]])
        per[c] = pd.DataFrame({"R_raw": raw, "R_adj": adj})
    return per, robust(per)


def robust(per):
    lines = {}
    for c, df in per.items():
        lines.setdefault(LINE_OF[c], []).append(df.R_adj)
    L = pd.DataFrame({l: pd.concat(v, axis=1).mean(axis=1) for l, v in lines.items()})
    n = L.notna().sum(axis=1)
    rob = L.mean(axis=1) / (1 + L.std(axis=1, ddof=0).fillna(0))
    cons = (L > 0).sum(axis=1) / n
    return pd.DataFrame({"R_rob": rob, "n_lines": n, "sign_consistency": cons}).join(L.add_prefix("Radj_"))


# ---------------------------------------------------------------- identity genes
def identity_genes(ctx, k=100):
    lines = {}
    for c, (s, nt, gm, nc) in ctx.items():
        if LINE_OF[c] == c:
            lines[c] = gm
    G = pd.DataFrame(lines).dropna()
    G = np.log(G / G.sum() * 1e4 + 0.01)
    out = {}
    for l in G.columns:
        others = [o for o in G.columns if o != l]
        if not others:
            out[l] = G[l].nlargest(k).index.tolist(); continue
        spec = G[l] - G[others].mean(axis=1)
        expressed = G[l] > G[l].median()
        out[l] = spec[expressed].nlargest(k).index.tolist()
    return {c: out[LINE_OF[c]] for c in ctx}


# ---------------------------------------------------------------- network
def load_string(min_score=700):
    info = pd.read_csv(os.path.join(RAW, "string_info.txt.gz"), sep="\t", usecols=[0, 1])
    name = dict(zip(info.iloc[:, 0], info.iloc[:, 1]))
    L = pd.read_csv(os.path.join(RAW, "string_links.txt.gz"), sep=" ")
    L = L[L.combined_score >= min_score]
    L["a"] = L.protein1.map(name); L["b"] = L.protein2.map(name)
    L = L[L.a < L.b]  # undirected, each pair once
    G = nx.Graph()
    G.add_weighted_edges_from(zip(L.a, L.b, L.combined_score / 1000.0))
    return G


def neighbour_mean(G, values, nodes):
    out = {}
    for n in nodes:
        if n not in G: continue
        num = den = 0.0
        for m, d in G[n].items():
            v = values.get(m)
            if v is not None:
                num += d["weight"] * v; den += d["weight"]
        if den > 0: out[n] = num / den
    return pd.Series(out, dtype=float)


def randomised(G, seed):
    H = G.copy()
    nx.double_edge_swap(H, nswap=2 * H.number_of_edges(), max_tries=20 * H.number_of_edges(), seed=seed)
    w = np.array([d["weight"] for _, _, d in G.edges(data=True)])
    np.random.default_rng(seed).shuffle(w)
    for (u, v), x in zip(H.edges(), w):
        H[u][v]["weight"] = x
    return H


def h3_network(per, G, nulls):
    res = {}
    for c, df in per.items():
        y = df.R_adj
        other = [o for o in per if LINE_OF[o] != LINE_OF[c]]
        transfer = pd.concat([per[o].R_adj for o in other], axis=1).mean(axis=1) if other else pd.Series(dtype=float)
        idx = np.array([p for p in y.index if p in G])
        rng = np.random.default_rng(0)
        folds = rng.permutation(len(idx)) % 5
        rows = []
        for f in range(5):
            tr, te = idx[folds != f], idx[folds == f]
            yv = y[tr].to_dict()
            b_te = neighbour_mean(G, yv, te)
            b_tr = neighbour_mean(G, yv, tr)
            nb = [neighbour_mean(H, yv, te) for H in nulls]
            # ridge on (transfer, network) fitted in training fold
            Xtr = pd.DataFrame({"a": transfer.reindex(tr), "b": b_tr.reindex(tr)}).dropna()
            beta = beta_a = None
            if len(Xtr) > 20:
                A = np.column_stack([np.ones(len(Xtr)), Xtr.values])
                lam = 1.0
                beta = np.linalg.solve(A.T @ A + lam * np.diag([0, 1, 1]), A.T @ y[Xtr.index].values)
                A1 = A[:, :2]; beta_a = np.linalg.solve(A1.T @ A1 + lam * np.diag([0, 1]), A1.T @ y[Xtr.index].values)
            for p in te:
                a_ = transfer.get(p, np.nan); b_ = b_te.get(p, np.nan)
                n_ = np.nanmean([s.get(p, np.nan) for s in nb]) if nb else np.nan
                c_ = beta @ [1, a_, b_] if (beta is not None and np.isfinite(a_) and np.isfinite(b_)) else np.nan
                ca_ = beta_a @ [1, a_] if (beta_a is not None and np.isfinite(a_)) else np.nan
                rows.append((p, y[p], a_, b_, n_, c_, y[tr].mean(), ca_))
        T = pd.DataFrame(rows, columns=["pert", "y", "transfer", "network", "null", "combined", "mean", "transfer_cal"]).set_index("pert")
        T.to_csv(os.path.join(RES, f"h3_{c}.csv"))
        out = {"n_network": int(T.network.notna().sum())}
        d = T.dropna(subset=["network", "null"])
        if len(d) > 10:
            e_b, e_n = (d.network - d.y).abs(), (d.null - d.y).abs()
            out.update(h3a_n=len(d), h3a_mae_network=e_b.mean(), h3a_mae_null=e_n.mean(),
                       h3a_mae_mean=(d["mean"] - d.y).abs().mean(),
                       h3a_rho_network=stats.spearmanr(d.network, d.y)[0],
                       h3a_rho_null=stats.spearmanr(d.null, d.y)[0],
                       h3a_p=stats.wilcoxon(e_b, e_n, alternative="less").pvalue)
        d = T.dropna(subset=["transfer", "combined"])
        if len(d) > 10:
            e_a, e_c = (d.transfer - d.y).abs(), (d.combined - d.y).abs()
            out.update(h3b_n=len(d), h3b_mae_transfer=e_a.mean(), h3b_mae_combined=e_c.mean(),
                       h3b_rho_transfer=stats.spearmanr(d.transfer, d.y)[0],
                       h3b_rho_combined=stats.spearmanr(d.combined, d.y)[0],
                       h3b_p=stats.wilcoxon(e_c, e_a, alternative="less").pvalue,
                       h3b_mae_transfer_cal=(d.transfer_cal - d.y).abs().mean(),
                       h3b_p_vs_cal=stats.wilcoxon(e_c, (d.transfer_cal - d.y).abs(), alternative="less").pvalue)
        res[c] = out
    for key in ("h3a_p", "h3b_p"):
        cs = [c for c in res if key in res[c]]
        adj = holm([res[c][key] for c in cs])
        for c, q in zip(cs, adj): res[c][key + "_holm"] = q
    return res


def holm(p):
    p = np.asarray(p, float); m = len(p)
    if m == 0: return []
    o = np.argsort(p); adj = np.empty(m)
    run = 0
    for i, j in enumerate(o):
        run = max(run, min(1, (m - i) * p[j])); adj[j] = run
    return adj.tolist()


# ---------------------------------------------------------------- evidence
def cellage():
    z = zipfile.ZipFile(os.path.join(RAW, "cellage.zip"))
    n = [x for x in z.namelist() if x.endswith(".tsv")][0]
    df = pd.read_csv(io.BytesIO(z.read(n)), sep="\t")
    eff = df.groupby("Gene symbol")["Senescence Effect"].agg(lambda s: s.mode()[0] if s.nunique() == 1 else "Mixed")
    return eff


def genage_human_labels():
    z = zipfile.ZipFile(os.path.join(RAW, "genage_models.zip"))
    df = pd.read_csv(io.BytesIO(z.read("genage_models.csv")))
    org = {"Mus musculus": "mmusculus", "Caenorhabditis elegans": "celegans",
           "Drosophila melanogaster": "dmelanogaster", "Saccharomyces cerevisiae": "scerevisiae"}
    votes = []
    for o, key in org.items():
        f = os.path.join(RAW, f"ortho_{key}.tsv")
        if not os.path.exists(f) or os.path.getsize(f) < 1000: continue
        orth = pd.read_csv(f, sep="\t").dropna()
        orth.columns = ["human", "model", "type"]
        orth = orth[orth.type.isin(["best"])]
        sub = df[(df.organism == o) & df["longevity influence"].isin(["Pro-Longevity", "Anti-Longevity"])]
        m = sub.merge(orth, left_on="symbol", right_on="model")
        votes.append(m[["human", "longevity influence", "organism"]])
    if not votes: return pd.Series(dtype=object), 0
    V = pd.concat(votes)
    tab = V.groupby(["human", "longevity influence"]).size().unstack(fill_value=0)
    for col in ["Pro-Longevity", "Anti-Longevity"]:
        if col not in tab: tab[col] = 0
    lab = pd.Series(np.where(tab["Anti-Longevity"] > tab["Pro-Longevity"], "Anti",
                             np.where(tab["Pro-Longevity"] > tab["Anti-Longevity"], "Pro", "Tie")), index=tab.index)
    return lab[lab != "Tie"], len(votes)


def mw_greater(x, y):
    x, y = np.asarray(x), np.asarray(y)
    if len(x) < 3 or len(y) < 3: return np.nan, np.nan
    u = stats.mannwhitneyu(x, y, alternative="greater")
    return u.pvalue, u.statistic / (len(x) * len(y))  # p, AUROC


def pareto(df, cols):
    X = df[cols].values
    nd = np.ones(len(X), bool)
    for i in range(len(X)):
        dom = np.all(X >= X[i], 1) & np.any(X > X[i], 1)
        if dom.any(): nd[i] = False
    return df.index[nd]


# ================================================================ main
def main():
    # ---------------- Aim 1
    if DISCOVERY == "fleischer": e, m = aging.load_fleischer()
    elif DISCOVERY == "table": e, m = aging.load_table(os.environ["TRACE_EXPR"], os.environ["TRACE_META"])  # any cohort as two CSV files
    else: e, m = aging.load_gtex(f"{DISCOVERY}_tpm.gct.gz")
    ts_f = aging.age_tstats(e, m)
    sig_f = aging.signature(ts_f)
    R["aging"] = {"fleischer_n_donors": len(m), "fleischer_age_min": float(m.age.min()),
                  "fleischer_age_max": float(m.age.max()), "fleischer_n_female": int((m.sex == "F").sum()),
                  "fleischer_n_genes": len(ts_f), "fleischer_n_fdr05": int((ts_f.q < 0.05).sum()),
                  "sig_up": int((sig_f.t > 0).sum()), "sig_down": int((sig_f.t < 0).sum())}
    ts_f.to_csv(os.path.join(RES, "age_t_fleischer.csv"))
    k80 = (m.age <= 80).values
    ts80 = aging.age_tstats(e.loc[:, k80], m[k80])
    R["explore"] = {"fleischer_18_80_n": int(k80.sum()), "fleischer_18_80_fdr05": int((ts80.q < 0.05).sum())}
    gt = {}
    for key, fn in [("gtex_fibroblast", "gtex_fibroblasts_tpm.gct.gz"), ("gtex_skin", "gtex_skin_ns_tpm.gct.gz")]:
        if not os.path.exists(os.path.join(RAW, fn)): log("missing", fn); continue
        if fn.startswith(DISCOVERY): continue  # the discovery cohort is not its own replication
        ge, gmeta = aging.load_gtex(fn)
        tsg = aging.age_tstats(ge, gmeta)
        tsg.to_csv(os.path.join(RES, f"age_t_{key}.csv"))
        gt[key] = tsg
        sh = ts_f.index.intersection(tsg.index)
        ss = sig_f.index.intersection(tsg.index)
        R["aging"][key] = {"n_donors": len(gmeta), "n_genes": len(tsg), "n_fdr05": int((tsg.q < 0.05).sum()),
                           "rho_t_all": stats.spearmanr(ts_f.t[sh], tsg.t[sh])[0], "n_shared": len(sh),
                           "sig_sign_concord": float((np.sign(sig_f.t[ss]) == np.sign(tsg.t[ss])).mean()),
                           "sig_sign_binom_p": stats.binomtest(int((np.sign(sig_f.t[ss]) == np.sign(tsg.t[ss])).sum()), len(ss)).pvalue,
                           "sig_n": len(ss),
                           "sig_rho_t": stats.spearmanr(sig_f.t[ss], tsg.t[ss])[0]}
        log(key, R["aging"][key])
    # signature composition
    R["aging"]["sig_frac_prolif"] = float(np.mean([g in PROLIF_REMOVE for g in sig_f.index]))
    R["aging"]["sig_down_frac_prolif"] = float(np.mean([g in PROLIF_REMOVE for g in sig_f[sig_f.t < 0].index]))
    R["aging"]["sig_up_frac_prolif"] = float(np.mean([g in PROLIF_REMOVE for g in sig_f[sig_f.t > 0].index]))
    enr = {}
    for k, gs in HM.items():
        for d, sub in (("up", sig_f[sig_f.t > 0]), ("down", sig_f[sig_f.t < 0])):
            a_ = len(set(sub.index) & gs); bg = set(ts_f.index)
            K = len(gs & bg)
            p = stats.hypergeom.sf(a_ - 1, len(bg), K, len(sub))
            enr[(k, d)] = (a_, K, p)
    E = pd.DataFrame([(k, d, *v) for (k, d), v in enr.items()], columns=["set", "dir", "overlap", "size", "p"])
    E["q"] = stats.false_discovery_control(E.p.values)
    E.sort_values("p").to_csv(os.path.join(RES, "signature_hallmark_enrichment.csv"), index=False)

    # ---------------- contexts
    ctx = {}
    for c in perturb.CONTEXTS:
        if os.path.exists(os.path.join(RAW, perturb.CONTEXTS[c])) or os.path.exists(os.path.join(PROC, f"ctx_{c}.pkl")):
            ctx[c] = perturb.load_context(c)
    R["contexts"] = {c: {"n_perturbed_genes": int(v[0].shape[0]), "n_measured_genes": int(v[0].shape[1]),
                         "n_nt": int(v[1].shape[0]),
                         "sig_genes_measured": int(len(sig_f.index.intersection(v[0].columns)))} for c, v in ctx.items()}
    log("contexts", R["contexts"])
    ident = identity_genes(ctx)
    progs = {c: program_shifts(v[0], ident[c]) for c, v in ctx.items()}
    per, rob = all_scores(ctx, sig_f.t, progs)
    for c in ctx:
        per[c] = per[c].join(progs[c])
        per[c].to_csv(os.path.join(RES, f"scores_{c}.csv"))
    # negative controls: non-targeting pseudobulks
    R["negctrl"] = {}
    for c, (shift, nt, gm, nc) in ctx.items():
        nts = score(nt, sig_f.t)
        thr = np.quantile(nts, 0.95)
        R["negctrl"][c] = {"nt_mean": float(nts.mean()), "nt_sd": float(nts.std()), "nt_q95": float(thr),
                           "pert_sd": float(per[c].R_raw.std()),
                           "frac_pert_above_nt95": float((per[c].R_raw > thr).mean()),
                           "frac_pert_below_nt05": float((per[c].R_raw < np.quantile(nts, 0.05)).mean()),
                           "mw_p_abs": stats.mannwhitneyu(per[c].R_raw.abs(), nts.abs(), alternative="greater").pvalue}

    # ---------------- H1
    line_ctx = [c for c in ["K562", "RPE1", "HepG2", "Jurkat"] if c in per]
    pairs = [(a_, b_) for i, a_ in enumerate(line_ctx) for b_ in line_ctx[i + 1:]]
    h1 = {}
    for a_, b_ in pairs:
        sh = per[a_].index.intersection(per[b_].index)
        r_, p_ = stats.spearmanr(per[a_].R_raw[sh], per[b_].R_raw[sh], alternative="greater")
        ra, _ = stats.spearmanr(per[a_].R_adj[sh], per[b_].R_adj[sh])
        k = max(1, int(0.05 * len(sh)))
        ta = set(per[a_].R_raw[sh].nlargest(k).index); tb = set(per[b_].R_raw[sh].nlargest(k).index)
        h1[f"{a_}-{b_}"] = {"n": len(sh), "rho_raw": r_, "p": p_, "rho_adj": ra,
                            "top5_overlap": len(ta & tb) / k, "top5_expected": k / len(sh)}
    for (kk, q) in zip(list(h1), holm([v["p"] for v in h1.values()])): h1[kk]["p_holm"] = q
    R["H1"] = {"pairs": h1, "supported": bool(len(h1) > 0 and all(v["p_holm"] < 0.05 for v in h1.values()))}
    # H1b
    if "K562gw" in per and len(line_ctx) >= 2:
        common = per["K562"].index.intersection(per["K562gw"].index)
        for c in line_ctx: common = common.intersection(per[c].index)
        V = pd.DataFrame({c: per[c].R_raw[common] for c in line_ctx + ["K562gw"]})
        Vr = V.rank()
        def stat(ix):
            C = Vr.iloc[ix].corr().values
            cols = list(V.columns)
            within = C[cols.index("K562"), cols.index("K562gw")]
            between = np.mean([C[cols.index(a_), cols.index(b_)] for a_, b_ in pairs])
            return within, between
        w0, b0 = stat(np.arange(len(V)))
        bs = np.array([np.subtract(*stat(RNG.integers(0, len(V), len(V)))) for _ in range(2000)])
        R["H1b"] = {"n_common": len(V), "rho_within_K562": w0, "rho_between_mean": b0, "diff": w0 - b0,
                    "diff_ci_lower_one_sided95": float(np.quantile(bs, 0.05)),
                    "p_boot": float((bs <= 0).mean()),
                    "supported": bool(np.quantile(bs, 0.05) > 0)}
        sh = per["K562"].index.intersection(per["K562gw"].index)
        R["H1b"]["rho_within_allshared"] = stats.spearmanr(per["K562"].R_raw[sh], per["K562gw"].R_raw[sh])[0]
        R["H1b"]["n_within_allshared"] = len(sh)
        # same comparison for magnitude (positive control for reproducibility)
        R["H1b"]["rho_adj_within_allshared"] = stats.spearmanr(per["K562"].R_adj[sh], per["K562gw"].R_adj[sh])[0]
        R["H1b"]["prolif_rho_within"] = stats.spearmanr(per["K562"].prolif[sh], per["K562gw"].prolif[sh])[0]
        R["H1b"]["magnitude_rho_within"] = stats.spearmanr(per["K562"].magnitude[sh], per["K562gw"].magnitude[sh])[0]
    log("H1", R["H1"], R.get("H1b"))

    # ---------------- H2
    h2 = {}
    for c in per:
        r_, p_ = stats.spearmanr(per[c].R_raw, per[c].prolif)
        r2, p2 = stats.spearmanr(per[c].R_adj, per[c].prolif)
        k = max(1, int(0.05 * len(per[c])))
        top_raw = set(per[c].R_raw.nlargest(k).index); top_adj = set(per[c].R_adj.nlargest(k).index)
        h2[c] = {"rho_raw_prolif": r_, "p": p_, "rho_adj_prolif": r2,
                 "rho_raw_identity": stats.spearmanr(per[c].R_raw, per[c].identity)[0],
                 "rho_raw_stress": stats.spearmanr(per[c].R_raw, per[c].stress)[0],
                 "rho_raw_magnitude": stats.spearmanr(per[c].R_raw, per[c].magnitude)[0],
                 "rho_raw_adj": stats.spearmanr(per[c].R_raw, per[c].R_adj)[0],
                 "top5_overlap_raw_adj": len(top_raw & top_adj) / k, "k": k}
    for (kk, q) in zip(list(h2), holm([v["p"] for v in h2.values()])): h2[kk]["p_holm"] = q
    R["H2"] = {"contexts": h2, "supported": bool(all(v["p_holm"] < 0.05 for v in h2.values()))}
    log("H2", R["H2"])

    # adjusted-score cross-context concordance (descriptive)
    rob.to_csv(os.path.join(RES, "robust_scores.csv"))

    # ---------------- H3
    if os.path.exists(os.path.join(RAW, "string_links.txt.gz")):
        G = load_string()
        log("STRING", G.number_of_nodes(), G.number_of_edges())
        nulls = [randomised(G, s) for s in range(10)]
        log("randomised")
        h3 = h3_network(per, G, nulls)
        R["H3"] = {"contexts": h3, "string_nodes": G.number_of_nodes(), "string_edges": G.number_of_edges()}
        sig_a = [c for c in h3 if h3[c].get("h3a_p_holm", 1) < 0.05]
        sig_b = [c for c in h3 if h3[c].get("h3b_p_holm", 1) < 0.05]
        R["H3"]["H3a_supported"] = len(sig_a) > len(h3) / 2
        R["H3"]["H3b_supported"] = len(sig_b) > len(h3) / 2
        log("H3", R["H3"])

    # ---------------- H4
    ca = cellage()
    gl, n_org = genage_human_labels()
    R["H4"] = {"n_organisms_mapped": n_org}
    tests = []
    for name, lab, pos, neg in [("H4a_cellage", ca, "Induces", "Inhibits"), ("H4b_genage", gl, "Anti", "Pro")]:
        x = rob.R_rob.reindex(lab[lab == pos].index).dropna()
        y = rob.R_rob.reindex(lab[lab == neg].index).dropna()
        p_, auc = mw_greater(x, y)
        d = {"n_pos": len(x), "n_neg": len(y), "auroc": auc, "p": p_,
             "median_pos": float(x.median()) if len(x) else None, "median_neg": float(y.median()) if len(y) else None}
        # exploratory: single context raw
        for c in per:
            xx = per[c].R_raw.reindex(lab[lab == pos].index).dropna(); yy = per[c].R_raw.reindex(lab[lab == neg].index).dropna()
            pp, aa = mw_greater(xx, yy)
            d[f"raw_{c}"] = {"n_pos": len(xx), "n_neg": len(yy), "auroc": aa, "p": pp}
            xx = per[c].R_adj.reindex(lab[lab == pos].index).dropna(); yy = per[c].R_adj.reindex(lab[lab == neg].index).dropna()
            pp, aa = mw_greater(xx, yy)
            d[f"adj_{c}"] = {"n_pos": len(xx), "n_neg": len(yy), "auroc": aa, "p": pp}
        R["H4"][name] = d; tests.append(name)
    ps = [R["H4"][t]["p"] for t in tests]
    for t, q in zip(tests, holm([1 if not np.isfinite(p) else p for p in ps])): R["H4"][t]["p_holm"] = q
    R["H4"]["supported_H4a"] = R["H4"]["H4a_cellage"]["p_holm"] < 0.05
    R["H4"]["supported_H4b"] = R["H4"]["H4b_genage"]["p_holm"] < 0.05
    log("H4", R["H4"])

    # ---------------- H5 donor bootstrap and cross-cohort signature
    boots, tops = [], []
    k = max(1, int(0.05 * rob.R_rob.notna().sum()))
    top_full = set(rob.R_rob.nlargest(k).index)
    topcount = pd.Series(0, index=rob.index, dtype=float)
    NB = int(os.environ.get("TRACE_NBOOT", 200))
    for b in range(NB):
        ix = RNG.integers(0, len(m), len(m))
        mb = m.iloc[ix].copy(); eb = e.iloc[:, ix]
        eb.columns = [f"s{i}" for i in range(len(ix))]; mb.index = eb.columns
        tsb = aging.age_tstats(eb, mb)
        sb = aging.signature(tsb)
        _, rb = all_scores(ctx, sb.t, progs)
        sh = rob.index.intersection(rb.index)
        boots.append(stats.spearmanr(rob.R_rob[sh], rb.R_rob[sh])[0])
        tb = set(rb.R_rob.nlargest(k).index); tops.append(len(tb & top_full) / k)
        topcount[list(tb & set(topcount.index))] += 1
        if b % 20 == 0: log("boot", b, boots[-1])
    rob["bootstrap_top5_support"] = topcount / NB
    R["H5"] = {"n_boot": NB, "median_rho": float(np.median(boots)), "q05_rho": float(np.quantile(boots, 0.05)),
               "q95_rho": float(np.quantile(boots, 0.95)), "median_top5_overlap": float(np.median(tops)),
               "threshold": 0.8, "supported_stability": bool(np.median(boots) >= 0.8)}
    pd.Series(boots).to_csv(os.path.join(RES, "h5_boot_rho.csv"))
    for key, tsg in gt.items():
        sg = aging.signature(tsg)
        _, rg = all_scores(ctx, sg.t, progs)
        sh = rob.index.intersection(rg.index)
        r_, p_ = stats.spearmanr(rob.R_rob[sh], rg.R_rob[sh], alternative="greater")
        R["H5"][key] = {"rho": r_, "p": p_, "n": len(sh), "sig_up": int((sg.t > 0).sum()), "sig_down": int((sg.t < 0).sum()),
                        "sig_overlap_with_fleischer": len(set(sg.index) & set(sig_f.index))}
        kk = max(1, int(0.05 * len(sh)))
        R["H5"][key]["top5_overlap"] = len(set(rg.R_rob[sh].nlargest(kk).index) & set(rob.R_rob[sh].nlargest(kk).index)) / kk
        rob[f"R_rob_{key}"] = rg.R_rob
        for c in line_ctx:
            pc, _ = all_scores({c: ctx[c]}, sg.t, {c: progs[c]})
            R["H5"][key][f"rho_raw_{c}"] = stats.spearmanr(per[c].R_raw, pc[c].R_raw.reindex(per[c].index))[0]
    log("H5", R["H5"])

    # ---------------- multi-objective prioritisation
    P = rob[rob.n_lines >= 2].copy()
    lines_present = [l for l in ["K562", "RPE1", "HepG2", "Jurkat"] if f"Radj_{l}" in rob]
    def line_mean(col):
        cols = []
        for l in lines_present:
            cs = [c for c in per if LINE_OF[c] == l]
            cols.append(pd.concat([per[c][col] for c in cs], axis=1).mean(axis=1).rename(l))
        return pd.concat(cols, axis=1).mean(axis=1)
    P["identity_shift"] = line_mean("identity").reindex(P.index)
    P["prolif_shift"] = line_mean("prolif").reindex(P.index)
    P["stress_shift"] = line_mean("stress").reindex(P.index)
    P["magnitude"] = line_mean("magnitude").reindex(P.index)
    ont = {}
    for c, (shift, nt, gm, nc) in ctx.items():
        for g in P.index:
            if g in shift.index and g in shift.columns:
                ont.setdefault(g, []).append(float(shift.at[g, g]))
    P["ontarget_z"] = pd.Series({g: np.mean(v) for g, v in ont.items()}).reindex(P.index)
    # annotation only (not an objective): target not measured or knockdown weaker than z = -0.5
    P["flag_weak_or_unmeasured_knockdown"] = ~(P.ontarget_z < -0.5)
    P["obj_identity"] = -P.identity_shift.abs()
    P["obj_prolif"] = -P.prolif_shift.clip(lower=0)
    P["obj_stress"] = -P.stress_shift.clip(lower=0)
    objs = ["R_rob", "sign_consistency", "obj_identity", "obj_prolif", "obj_stress"]
    front = pareto(P[P.R_rob > 0], objs)
    P["pareto"] = P.index.isin(front)
    onc = pd.read_csv(os.path.join(RAW, "oncokb_cancer_genes.tsv"), sep="\t")
    tsg_set = set(onc[onc["Gene Type"].astype(str).str.contains("TSG")]["Hugo Symbol"])
    ess = set(ctx["K562"][0].index) if "K562" in ctx else set()
    P["flag_tumour_suppressor"] = P.index.isin(tsg_set)
    P["flag_essential_library"] = P.index.isin(ess)
    P["cellage"] = ca.reindex(P.index).fillna("unknown")
    P["genage"] = gl.reindex(P.index).fillna("unknown")
    P = P.join(rob[["bootstrap_top5_support"]], rsuffix="_")
    P.sort_values("R_rob", ascending=False).to_csv(os.path.join(RES, "prioritisation.csv"))
    F = P[P.pareto].sort_values("R_rob", ascending=False)
    R["pareto"] = {"n_candidates_ge2_lines": len(P), "n_positive": int((P.R_rob > 0).sum()), "n_front": len(F),
                   "front_top": F.head(25)[["R_rob", "sign_consistency", "identity_shift", "prolif_shift",
                                            "stress_shift", "magnitude", "ontarget_z", "flag_weak_or_unmeasured_knockdown", "bootstrap_top5_support", "flag_tumour_suppressor",
                                            "cellage", "genage"]].round(3).reset_index().to_dict("records")}
    # top robust by R_raw vs R_adj overlap across lines (what disappears)
    rraw = robust({c: per[c].assign(R_adj=per[c].R_raw) for c in per})
    k = max(1, int(0.05 * len(P)))
    t_raw = set(rraw.R_rob.reindex(P.index).nlargest(k).index); t_adj = set(P.R_rob.nlargest(k).index)
    R["pareto"]["n_front_weak_kd"] = int(F.flag_weak_or_unmeasured_knockdown.sum())
    R["pareto"]["front_magnitude_median"] = float(F.magnitude.median()); R["pareto"]["all_magnitude_median"] = float(P.magnitude.median())
    R["pareto"]["top5_raw_vs_adj_overlap"] = len(t_raw & t_adj) / k
    R["pareto"]["top5_raw_prolif_mean"] = float(P.prolif_shift[list(t_raw)].mean())
    R["pareto"]["top5_adj_prolif_mean"] = float(P.prolif_shift[list(t_adj)].mean())
    R["pareto"]["all_prolif_mean"] = float(P.prolif_shift.mean())
    rob.to_csv(os.path.join(RES, "robust_scores.csv"))

    json.dump(R, open(os.path.join(RES, "results.json"), "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    log("done")


if __name__ == "__main__":
    main()








