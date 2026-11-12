"""Exploratory (Amendment 8, 6 November 2026): two checks on Box 1.

A. Like-for-like cohort transfer. The discovery signature comes mostly from donors aged 80 or older, an age stratum that
   GTEx lacks. Here the ranking is rebuilt from discovery donors aged 20-79 only (the GTEx age range; signature rule as
   in the cohort sensitivity analyses: FDR < 0.05, or the top 250 per direction if fewer than 50 pass) and compared with
   the GTEx-fibroblast ranking and with the main discovery ranking. Reference: top-5% overlap between donor bootstraps
   of the same 20-79 subset. No P value (genes are not independent), as for the skin comparison.

B. Positive example for Box 1. GTEx fibroblasts are split into two random halves, treated as two cohorts, and the
   cohort (top-200 agreement against a permuted-age null), effect-size (top-5% overlap against donor bootstraps) and
   stability (median donor-bootstrap rho >= 0.8) checks are applied. The halves share platform, center and processing,
   so this shows that the checks can pass, not that they would pass for two truly independent cohorts.
   One primary split (seed 0) gets the full checks; 20 splits give the spread of the point estimates.
Outputs: results/x_box1_positive.json, results/x_box1_positive_splits.csv
"""
import os, sys, json, time
import multiprocessing as mp
import numpy as np, pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as _cfg
import aging, perturb
import trace_analysis as T
import rev3_cohort as C

RES = T.RES; E18 = os.path.join(RES, "e18_gtex_discovery")
NPROC = _cfg.N_PROC; TOPK = 200
N_SPLITS = _cfg._int("N_BOX1_SPLITS", 20); N_PERM = _cfg._int("N_BOX1_PERM", 200); N_BOOT = _cfg._int("N_BOX1_BOOT", 100)
_W = {}


def log(*a): print(time.strftime("%H:%M:%S"), *a, flush=True)


def _init():
    ctx = {c: perturb.load_context(c) for c in perturb.CONTEXTS}; ident = T.identity_genes(ctx)
    _W["ctx"] = ctx; _W["progs"] = {c: T.program_shifts(v[0], ident[c]) for c, v in ctx.items()}
    e0, m0 = aging.load_fleischer(); k = (m0.age >= 20) & (m0.age < 80); _W["f1"] = (e0.loc[:, m0.index[k]], m0[k])
    _W["g"] = aging.load_gtex("gtex_fibroblasts_tpm.gct.gz")


def ranking(e, m):
    ts = aging.age_tstats(e, m); s = C.sig_from(ts)
    return T.all_scores(_W["ctx"], s.t, _W["progs"])[1], ts, s


def k5(rob): return max(1, int(0.05 * rob.notna().sum()))


def compare(a, b):
    """Spearman over shared genes and top-5% overlap (k from the first ranking), all genes and genes in >= 2 lines."""
    out = {}
    for lab, gs in [("all", a.index), ("multi", a.index[a.n_lines >= 2])]:
        sh = gs.intersection(b.index); x, y = a.R_rob[sh], b.R_rob[sh]; ok = x.notna() & y.notna(); x, y = x[ok], y[ok]
        k = k5(x); out[lab] = {"rho": float(stats.spearmanr(x, y)[0]), "top5_overlap": len(set(x.nlargest(k).index) & set(y.nlargest(k).index)) / k, "n": int(len(x)), "k": int(k)}
    return out


def top_rho(ta, tb):
    sh = ta.index.intersection(tb.index); top = ta.t[sh].abs().nlargest(TOPK).index
    return float(stats.spearmanr(ta.t[top], tb.t[top])[0])


def halves(seed):
    e, m = _W["g"]; r = np.random.default_rng(seed); ix = r.permutation(len(m)); h = len(m) // 2
    return (e.iloc[:, ix[:h]], m.iloc[ix[:h]]), (e.iloc[:, ix[h:]], m.iloc[ix[h:]])


# ---- workers
def _f1_boot(b):
    e, m = _W["f1"]; r = np.random.default_rng(92_000 + b); ix = r.integers(0, len(m), len(m))
    mb = m.iloc[ix].copy(); eb = e.iloc[:, ix]; eb.columns = mb.index = [f"s{i}" for i in range(len(ix))]
    return ranking(eb, mb)[0]


def _split(s):
    (ea, ma), (eb, mb) = halves(s); ra, ta, _ = ranking(ea, ma); rb, tb, _ = ranking(eb, mb); c = compare(ra, rb)
    return {"split": s, "top200_rho": top_rho(ta, tb), "rank_rho": c["all"]["rho"], "top5_overlap": c["all"]["top5_overlap"],
            "rank_rho_multi": c["multi"]["rho"], "top5_overlap_multi": c["multi"]["top5_overlap"]}


def _perm(b):
    (ea, ma), (eb, mb) = halves(0); r = np.random.default_rng(90_000 + b)
    ma = ma.copy(); mb = mb.copy(); ma["age"] = r.permutation(ma.age.values); mb["age"] = r.permutation(mb.age.values)
    return top_rho(aging.age_tstats(ea, ma), aging.age_tstats(eb, mb))


def _boot_a(b):
    (ea, ma), _ = halves(0); r = np.random.default_rng(91_000 + b); ix = r.integers(0, len(ma), len(ma))
    mb = ma.iloc[ix].copy(); eb = ea.iloc[:, ix]; eb.columns = mb.index = [f"s{i}" for i in range(len(ix))]
    return ranking(eb, mb)[0]


def main():
    _init(); out = {"n_splits": N_SPLITS, "n_perm": N_PERM, "n_boot": N_BOOT}
    # A. like-for-like
    e1, m1 = _W["f1"]; r1, ts1, s1 = ranking(e1, m1)
    disc = pd.read_csv(os.path.join(RES, "robust_scores.csv"), index_col=0); gtex = pd.read_csv(os.path.join(E18, "robust_scores.csv"), index_col=0)
    with mp.Pool(NPROC, initializer=_init) as pool:
        fb = list(pool.imap(_f1_boot, range(N_BOOT), chunksize=2)); log("F1 boots done")
        sp = pd.DataFrame(pool.map(_split, range(N_SPLITS))); log("splits done")
        nul = np.array(pool.map(_perm, range(N_PERM), chunksize=4)); log("perms done")
        ab = list(pool.imap(_boot_a, range(N_BOOT), chunksize=2)); log("half-A boots done")
    fbc = [compare(r1, b)["all"] for b in fb]
    out["like_for_like"] = {"n_donors": int(len(m1)), "n_fdr05": int((ts1.q < 0.05).sum()), "n_fdr05_up": int(((ts1.q < 0.05) & (ts1.t > 0)).sum()),
                            "n_fdr05_down": int(((ts1.q < 0.05) & (ts1.t < 0)).sum()), "fallback_top250": bool(len(s1) == 500 and (s1.q >= 0.05).any()),
                            "vs_gtex": compare(r1, gtex), "vs_main_discovery": compare(r1, disc),
                            "boot_median_rho": float(np.median([c["rho"] for c in fbc])), "boot_median_top5_overlap": float(np.median([c["top5_overlap"] for c in fbc]))}
    # B. positive example
    sp.to_csv(os.path.join(RES, "x_box1_positive_splits.csv"), index=False)
    (ea, ma), (eb, mb) = halves(0); ra, ta, _ = ranking(ea, ma); rb, tb, _ = ranking(eb, mb); c = compare(ra, rb)
    abc = [compare(ra, b)["all"] for b in ab]
    obs = top_rho(ta, tb); hi = float(np.quantile(nul, 0.975)); bo = float(np.median([x["top5_overlap"] for x in abc])); br = float(np.median([x["rho"] for x in abc]))
    out["split_gtex"] = {"n_donors_half": [int(len(ma)), int(len(mb))], "primary_split": 0,
                         "cohort": {"top200_rho": obs, "null_hi": hi, "null_lo": float(np.quantile(nul, 0.025)), "result": "pass" if obs > hi else "fail"},
                         "effect_size": {"top5_overlap": c["all"]["top5_overlap"], "boot_top5_overlap": bo, "result": "pass" if c["all"]["top5_overlap"] >= bo / 2 else "fail",
                                         "rank_rho": c["all"]["rho"], "rank_rho_multi": c["multi"]["rho"], "top5_overlap_multi": c["multi"]["top5_overlap"]},
                         "stability": {"boot_median_rho": br, "boot_q025": float(np.quantile([x["rho"] for x in abc], 0.025)), "boot_q975": float(np.quantile([x["rho"] for x in abc], 0.975)),
                                       "result": "pass" if br >= 0.8 else "fail"},
                         "splits": {k: [float(sp[k].min()), float(sp[k].median()), float(sp[k].max())] for k in ["top200_rho", "rank_rho", "top5_overlap", "rank_rho_multi", "top5_overlap_multi"]}}
    json.dump(out, open(os.path.join(RES, "x_box1_positive.json"), "w"), indent=1); log("done"); print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
