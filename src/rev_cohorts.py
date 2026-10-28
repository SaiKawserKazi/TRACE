"""Exploratory (Amendment 1): cohort comparability and sensitivity of the age-effect replication.

All analyses here are exploratory and were added after the pre-specified results were known.
Outputs: results/rev_cohort_table.csv, results/rev_cohort_sensitivity.csv, results/rev_cohorts.json
"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import config as _cfg
import os, sys, json, time
import numpy as np, pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
import aging
from fleischer_meta import load_meta

ROOT = str(_cfg.ROOT)
RAW, RES = os.path.join(ROOT, "data", "raw"), os.path.join(ROOT, "results")
NB = _cfg.N_COHORT_BOOT
rng = np.random.default_rng(1)


def log(*a): print(time.strftime("%H:%M:%S"), *a, flush=True)


# ------------------------------------------------------------------ cohorts with metadata
def fleischer():
    e, m = aging.load_fleischer()
    M = load_meta().reindex(m.index)
    m = m.join(M)
    m["site2"] = np.where(m.site.str.lower().str.contains("arm"), "arm", "other")
    return e, m


def gtex(fname):
    e, m = aging.load_gtex(fname)
    s = pd.read_csv(os.path.join(RAW, "gtex_sample.txt"), sep="\t", low_memory=False, index_col=0)
    s = s.reindex(m.index)
    subj = pd.read_csv(os.path.join(RAW, "gtex_subject.txt"), sep="\t", index_col=0)
    m = m.assign(center=s.SMCENTER.values, batch_type=s.SMGEBTCHT.values, na_batch=s.SMNABTCH.values,
                 rin=s.SMRIN.values, ischemic=s.SMTSISCH.values,
                 hardy=subj.reindex(m.donor).DTHHRDY.fillna(-1).astype(int).astype(str).values)
    return e[e.mean(axis=1) >= 0.5], m


def describe(name, e, m, extra):
    bins = pd.cut(m.age, [17, 29, 39, 49, 59, 69, 79, 120], labels=["18-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"])
    d = {"cohort": name, "n_donors": len(m), "age_min": m.age.min(), "age_median": m.age.median(), "age_max": m.age.max(),
         "n_female": int((m.sex == "F").sum()), "genes_mean_log2_ge1": int((e.mean(axis=1) >= 1).sum())}
    d.update({f"age_{k}": int(v) for k, v in bins.value_counts().sort_index().items()})
    d.update(extra)
    return d


# ------------------------------------------------------------------ specifications
FL_SPECS = {  # name -> (sample mask function, covariates, min_mean)
    "F0 protocol (>=18 y; age+sex)": (lambda m: m.age >= 18, ["sex"], 1.0),
    "F1 age 20-79 (GTEx range)": (lambda m: (m.age >= 20) & (m.age < 80), ["sex"], 1.0),
    "F2 + platform, repository, site": (lambda m: m.age >= 18, ["sex", "platform", "repository", "site2"], 1.0),
    "F3 age 20-79 + platform, repository, site": (lambda m: (m.age >= 20) & (m.age < 80), ["sex", "platform", "repository", "site2"], 1.0),
    "F4 NextSeq only": (lambda m: (m.age >= 18) & (m.platform == "NextSeq 500"), ["sex"], 1.0),
    "F5 arm biopsies only": (lambda m: (m.age >= 18) & (m.site2 == "arm"), ["sex"], 1.0),
    "F6 gene filter mean log2 >= 0.5": (lambda m: m.age >= 18, ["sex"], 0.5),
    "F7 gene filter mean log2 >= 2": (lambda m: m.age >= 18, ["sex"], 2.0),
}
GT_SPECS = {
    "G0 protocol (age+sex)": ["sex"],
    "G1 + centre, batch type, RIN, Hardy": ["sex", "center", "batch_type", "rin", "hardy"],
    "G2 G1 + nucleic-acid batch": ["sex", "center", "batch_type", "rin", "hardy", "na_batch"],
}
SKIN_EXTRA = ["ischemic"]


def fit(e, m, covars, min_mean, idx=None):
    if idx is not None:
        e = e.iloc[:, idx].copy(); m = m.iloc[idx].copy()
        e.columns = m.index = [f"s{i}" for i in range(len(idx))]
    return aging.age_tstats_cov(e, m, covars, min_mean)[["t", "q"]]


def sig_from(ts):
    s = aging.signature(ts)
    return s if (s.t > 0).sum() >= 50 and (s.t < 0).sum() >= 50 else aging.signature_topk(ts)


def compare(tf, tg, sig):
    sh = tf.index.intersection(tg.index); ss = sig.index.intersection(tg.index)
    return {"rho_all": stats.spearmanr(tf.t[sh], tg.t[sh])[0], "n_shared": len(sh),
            "sig_sign_agree": float((np.sign(sig.t[ss]) == np.sign(tg.t[ss])).mean()), "sig_n": len(ss),
            "rho_sig": stats.spearmanr(sig.t[ss], tg.t[ss])[0]}


import multiprocessing as mp
NPROC = _cfg.N_PROC
NSPLIT = _cfg.N_SPLIT_HALVES
_W = {}


def _init():
    _W["f"] = fleischer()
    _W["g"] = {"fibroblast": gtex("gtex_fibroblasts_tpm.gct.gz"), "skin": gtex("gtex_skin_ns_tpm.gct.gz")}


def _boot_task(b):
    """One donor-level bootstrap replicate (seeded by b) of every specification; returns all 48 comparisons."""
    r = np.random.default_rng(10_000 + b)
    ef, mf = _W["f"]
    Fb = {}
    for name, (mask, cov, mm) in FL_SPECS.items():
        k = mask(mf).values; e_, m_ = ef.loc[:, k], mf[k]
        Fb[name] = fit(e_, m_, cov, mm, r.integers(0, k.sum(), k.sum()))
    Gb = {}
    for coh, (e_, m_) in _W["g"].items():
        idx = r.integers(0, len(m_), len(m_))  # same donor resample for all models of a cohort
        for name, cov in GT_SPECS.items():
            Gb[f"{coh}: {name}"] = fit(e_, m_, cov + (SKIN_EXTRA if coh == "skin" else []), 1.0, idx)
    out = []
    for fn, tfb in Fb.items():
        sig = sig_from(tfb)
        for gn, tgb in Gb.items():
            c = compare(tfb, tgb, sig); out.append((b, fn, gn, c["rho_all"], c["sig_sign_agree"]))
    return out


def _half_task(s):
    r = np.random.default_rng(20_000 + s)
    res = []
    for lab, (e_, m_) in {"Fleischer split-half": _W["f"], "GTEx fibroblast split-half": _W["g"]["fibroblast"],
                          "GTEx skin split-half": _W["g"]["skin"]}.items():
        perm = r.permutation(len(m_)); a, c = perm[: len(m_) // 2], perm[len(m_) // 2:]
        ta, tc = fit(e_, m_, ["sex"], 1.0, a), fit(e_, m_, ["sex"], 1.0, c)
        sh = ta.index.intersection(tc.index); res.append((s, lab, stats.spearmanr(ta.t[sh], tc.t[sh])[0]))
    return res


def main():
    ef, mf = fleischer()
    eg, mg = gtex("gtex_fibroblasts_tpm.gct.gz")
    es, ms = gtex("gtex_skin_ns_tpm.gct.gz")
    tab = [describe("Fleischer (discovery)", ef, mf, {
        "tissue": "primary dermal fibroblasts (cultured)", "source": "Coriell repositories " + ", ".join(f"{k}: {v}" for k, v in mf.repository.value_counts().items()),
        "platform": ", ".join(f"{k}: {v}" for k, v in mf.platform.value_counts().items()),
        "biopsy_site": ", ".join(f"{k}: {v}" for k, v in mf.site2.value_counts().items()),
        "age_resolution": "years", "expression": "FPKM (summed per gene), log2(x+1)",
        "covariates_available": "sex, platform, repository, biopsy site, ethnicity", "passage": "not reported"}),
        describe("GTEx v10 cultured fibroblasts", eg, mg, {
            "tissue": "fibroblasts cultured from lower-leg skin biopsies", "source": "GTEx post-mortem/organ donors",
            "platform": ", ".join(f"{k}: {v}" for k, v in mg.center.value_counts().items()) + " (centre)",
            "biopsy_site": "lower leg", "age_resolution": "10-year bins (midpoints used)", "expression": "TPM (summed per symbol), log2(x+1)",
            "covariates_available": "sex, centre, batch type, nucleic-acid batch, RIN, Hardy scale", "passage": "not reported"}),
        describe("GTEx v10 skin, not sun-exposed", es, ms, {
            "tissue": "suprapubic skin (tissue, not cultured)", "source": "GTEx post-mortem/organ donors",
            "platform": ", ".join(f"{k}: {v}" for k, v in ms.center.value_counts().items()) + " (centre)",
            "biopsy_site": "suprapubic", "age_resolution": "10-year bins (midpoints used)", "expression": "TPM (summed per symbol), log2(x+1)",
            "covariates_available": "sex, centre, batch type, nucleic-acid batch, RIN, ischaemic time, Hardy scale", "passage": "not applicable"})]
    pd.DataFrame(tab).to_csv(os.path.join(RES, "rev_cohort_table.csv"), index=False)
    log("cohort table done")

    # point estimates for every specification
    F = {}
    for name, (mask, cov, mm) in FL_SPECS.items():
        k = mask(mf).values
        F[name] = (fit(ef.loc[:, k], mf[k], cov, mm), int(k.sum()))
        log(name, "n", k.sum(), "FDR<0.05", int((F[name][0].q < 0.05).sum()))
    G = {}
    for coh, (e_, m_) in {"fibroblast": (eg, mg), "skin": (es, ms)}.items():
        for name, cov in GT_SPECS.items():
            cov = cov + (SKIN_EXTRA if coh == "skin" else [])
            G[f"{coh}: {name}"] = (fit(e_, m_, cov, 1.0), len(m_))
            log(f"{coh}: {name}", "FDR<0.05", int((G[f'{coh}: {name}'][0].q < 0.05).sum()))
    # donor-level bootstraps and split halves, streamed in parallel (each task re-fits every specification once)
    with mp.Pool(NPROC, initializer=_init) as pool:
        boot = []
        for i, res in enumerate(pool.imap_unordered(_boot_task, range(NB), chunksize=4)):
            boot.extend(res)
            if i % 50 == 0: log("bootstrap", i)
        halves = list(pool.imap_unordered(_half_task, range(NSPLIT), chunksize=4))
    BT = pd.DataFrame(boot, columns=["b", "fn", "gn", "rho_all", "sig_sign_agree"])
    BT.to_csv(os.path.join(RES, "rev_cohort_bootstrap.csv.gz"), index=False)
    HV = pd.DataFrame([r for h in halves for r in h], columns=["split", "cohort", "rho"])
    ref = {lab: {"median": float(g.rho.median()), "q025": float(g.rho.quantile(0.025)), "q975": float(g.rho.quantile(0.975)), "n_splits": len(g)}
           for lab, g in HV.groupby("cohort")}
    log("split halves", ref)
    rows = []
    for fn, (tf, nf) in F.items():
        sig = sig_from(tf)
        for gn, (tg, ng) in G.items():
            pt = compare(tf, tg, sig)
            B = BT[(BT.fn == fn) & (BT.gn == gn)]
            rows.append({"discovery_spec": fn, "n_discovery": nf, "fdr05_discovery": int((tf.q < 0.05).sum()),
                         "signature_rule": "FDR<0.05 top 250/250" if (tf.q < 0.05).sum() >= 100 else "top 250/250 by t (few FDR genes)",
                         "replication_spec": gn, "n_replication": ng, "fdr05_replication": int((tg.q < 0.05).sum()),
                         **pt, "rho_all_lo": B.rho_all.quantile(0.025), "rho_all_hi": B.rho_all.quantile(0.975),
                         "sig_sign_agree_lo": B.sig_sign_agree.quantile(0.025), "sig_sign_agree_hi": B.sig_sign_agree.quantile(0.975)})
    S = pd.DataFrame(rows)
    S.to_csv(os.path.join(RES, "rev_cohort_sensitivity.csv"), index=False)
    # GTEx fibroblast vs skin
    tfib, tskin = G["fibroblast: G0 protocol (age+sex)"][0], G["skin: G0 protocol (age+sex)"][0]
    sh = tfib.index.intersection(tskin.index)
    out = {"n_boot": NB, "split_half": ref,
           "gtex_fibro_vs_skin_rho": stats.spearmanr(tfib.t[sh], tskin.t[sh])[0],
           "fdr05": {**{k: int((v[0].q < 0.05).sum()) for k, v in F.items()}, **{k: int((v[0].q < 0.05).sum()) for k, v in G.items()}},
           "n": {**{k: v[1] for k, v in F.items()}, **{k: v[1] for k, v in G.items()}},
           "fleischer_crosstab_age_repository": pd.crosstab(pd.cut(mf.age, [17, 39, 59, 79, 120]).astype(str), mf.repository).to_dict(),
           "fleischer_crosstab_age_platform": pd.crosstab(pd.cut(mf.age, [17, 39, 59, 79, 120]).astype(str), mf.platform).to_dict(),
           "fleischer_age_repository_test": stats.mannwhitneyu(mf.age[mf.repository == "AG"], mf.age[mf.repository == "GM"]).pvalue,
           "rho_all_range_over_specs": [float(S.rho_all.min()), float(S.rho_all.max())],
           "fibro_rho_range": [float(S[S.replication_spec.str.startswith("fibro")].rho_all.min()), float(S[S.replication_spec.str.startswith("fibro")].rho_all.max())],
           "skin_rho_range": [float(S[S.replication_spec.str.startswith("skin")].rho_all.min()), float(S[S.replication_spec.str.startswith("skin")].rho_all.max())],
           "fibro_sign_range": [float(S[S.replication_spec.str.startswith("fibro")].sig_sign_agree.min()), float(S[S.replication_spec.str.startswith("fibro")].sig_sign_agree.max())],
           "skin_sign_range": [float(S[S.replication_spec.str.startswith("skin")].sig_sign_agree.min()), float(S[S.replication_spec.str.startswith("skin")].sig_sign_agree.max())]}
    # keyed views for the manuscript
    def get(fn, gn): return S[(S.discovery_spec == fn) & (S.replication_spec == gn)].iloc[0].to_dict()
    out["key"] = {
        "F0_G0_fibro": get("F0 protocol (>=18 y; age+sex)", "fibroblast: G0 protocol (age+sex)"),
        "F0_G0_skin": get("F0 protocol (>=18 y; age+sex)", "skin: G0 protocol (age+sex)"),
        "F3_G2_fibro": get("F3 age 20-79 + platform, repository, site", "fibroblast: G2 G1 + nucleic-acid batch"),
        "F3_G2_skin": get("F3 age 20-79 + platform, repository, site", "skin: G2 G1 + nucleic-acid batch"),
        "F1_G0_fibro": get("F1 age 20-79 (GTEx range)", "fibroblast: G0 protocol (age+sex)"),
        "F2_G1_fibro": get("F2 + platform, repository, site", "fibroblast: G1 + centre, batch type, RIN, Hardy"),
    }
    json.dump(out, open(os.path.join(RES, "rev_cohorts.json"), "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    log("done")


if __name__ == "__main__":
    main()



