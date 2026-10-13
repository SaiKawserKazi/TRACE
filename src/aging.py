"""Aim 1: donor-aware aging signatures (protocol section 3)."""
import gzip, re, os
import numpy as np, pandas as pd
from scipy import stats

RAW = r"D:\Study\RejuvBench\data\raw"
OUT = r"D:\Study\RejuvBench\data\proc"
os.makedirs(OUT, exist_ok=True)


def parse_fleischer_sample(name):
    parts = name.split("_")
    m = re.match(r"(\d+)", parts[1])
    age = float(m.group(1)) if m else np.nan
    if "mos" in parts[1]:
        mo = re.search(r"(\d+)mos", parts[1])
        age += float(mo.group(1)) / 12 if mo else 0
    sex_tok = parts[2].lower()
    hgps = "hgps" in sex_tok
    sex = "F" if sex_tok.startswith("f") else ("M" if sex_tok.startswith("m") else None)
    return age, sex, hgps


def load_fleischer():
    df = pd.read_csv(os.path.join(RAW, "GSE113957_fpkm.txt.gz"), sep="\t")
    df["gene"] = df["Annotation/Divergence"].str.split("|").str[0]
    samples = df.columns[8:-1]
    expr = df.groupby("gene")[list(samples)].sum()
    meta = pd.DataFrame([parse_fleischer_sample(s) for s in samples], index=samples,
                        columns=["age", "sex", "hgps"])
    keep = (~meta.hgps) & (meta.age >= 18) & meta.sex.notna()
    meta = meta[keep]
    expr = np.log2(expr[meta.index] + 1)
    return expr, meta


def load_gtex(fname):
    expr = pd.read_csv(os.path.join(RAW, fname), sep="\t", skiprows=2, index_col=0)
    sym = expr.pop("Description")
    expr = np.log2(expr.groupby(sym.values).sum() + 1)
    subj = pd.read_csv(os.path.join(RAW, "gtex_subject.txt"), sep="\t", index_col=0)
    sid = pd.Index(["-".join(c.split("-")[:2]) for c in expr.columns])
    lo = subj.loc[sid, "AGE"].str.split("-").str[0].astype(float).values
    meta = pd.DataFrame({"age": lo + 4.5,
                         "sex": np.where(subj.loc[sid, "SEX"].values == 1, "M", "F"),
                         "donor": sid}, index=expr.columns)
    # one sample per donor
    first = ~meta.donor.duplicated()
    return expr.loc[:, first.values], meta[first.values]


def age_tstats(expr, meta, min_mean=1.0):
    """OLS expression ~ age + sex per gene; returns t-statistic, p, BH q for age."""
    expr = expr[expr.mean(axis=1) >= min_mean]
    X = np.column_stack([np.ones(len(meta)), meta.age.values, (meta.sex.values == "M").astype(float)])
    Y = expr.values.T  # samples x genes
    XtXi = np.linalg.inv(X.T @ X)
    B = XtXi @ X.T @ Y
    resid = Y - X @ B
    dof = X.shape[0] - X.shape[1]
    s2 = (resid ** 2).sum(0) / dof
    se = np.sqrt(s2 * XtXi[1, 1])
    t = B[1] / np.where(se > 0, se, np.nan)
    p = 2 * stats.t.sf(np.abs(t), dof)
    out = pd.DataFrame({"t": t, "p": p, "beta_per_decade": B[1] * 10}, index=expr.index)
    out = out.dropna()
    out["q"] = stats.false_discovery_control(out.p.values)
    return out


def age_tstats_cov(expr, meta, covars=(), min_mean=1.0, age_col="age"):
    """OLS expression ~ age + covariates per gene (review-round sensitivity analyses). Categorical covariates are
    dummy-coded (first level dropped); numeric covariates are entered linearly after median imputation."""
    expr = expr[expr.mean(axis=1) >= min_mean]
    cols = [np.ones(len(meta)), meta[age_col].values.astype(float)]
    for c in covars:
        v = meta[c]
        if v.dtype.kind in "biuf":
            cols.append(v.fillna(v.median()).values.astype(float))
        else:
            d = pd.get_dummies(v.astype(str), drop_first=True).values.astype(float)
            d = d[:, d.sum(0) > 0]
            cols.extend(d.T)
    X = np.column_stack(cols)
    rank = np.linalg.matrix_rank(X)
    if rank < X.shape[1]:  # drop collinear dummies
        q, r = np.linalg.qr(X)
        keep = np.abs(np.diag(r)) > 1e-8 * np.abs(np.diag(r)).max()
        keep[:2] = True
        X = X[:, keep]
    Y = expr.values.T
    XtXi = np.linalg.pinv(X.T @ X)
    B = XtXi @ X.T @ Y
    resid = Y - X @ B
    dof = X.shape[0] - X.shape[1]
    s2 = (resid ** 2).sum(0) / dof
    se = np.sqrt(s2 * XtXi[1, 1])
    t = B[1] / np.where(se > 0, se, np.nan)
    p = 2 * stats.t.sf(np.abs(t), dof)
    out = pd.DataFrame({"t": t, "p": p, "beta_per_decade": B[1] * 10}, index=expr.index).dropna()
    out["q"] = stats.false_discovery_control(out.p.values)
    out.attrs["n_params"] = X.shape[1]; out.attrs["dof"] = dof
    return out


def signature_topk(ts, cap=250):
    """Top-k signature without an FDR filter (sensitivity analysis for cohorts with few FDR genes)."""
    return pd.concat([ts.nlargest(cap, "t"), ts.nsmallest(cap, "t")])


def signature(ts, q=0.05, cap=250):
    sig = ts[ts.q < q]
    up = sig[sig.t > 0].nlargest(cap, "t")
    dn = sig[sig.t < 0].nsmallest(cap, "t")
    return pd.concat([up, dn])


if __name__ == "__main__":
    e, m = load_fleischer()
    print("Fleischer donors", len(m), "age range", m.age.min(), m.age.max(), m.sex.value_counts().to_dict())
    ts = age_tstats(e, m)
    ts.to_csv(os.path.join(OUT, "age_t_fleischer.csv"))
    s = signature(ts)
    s.to_csv(os.path.join(OUT, "signature_fleischer.csv"))
    print("expressed genes", len(ts), "FDR<0.05", (ts.q < 0.05).sum(), "signature up/down",
          (s.t > 0).sum(), (s.t < 0).sum())
    e.to_pickle(os.path.join(OUT, "fleischer_expr.pkl")); m.to_pickle(os.path.join(OUT, "fleischer_meta.pkl"))
