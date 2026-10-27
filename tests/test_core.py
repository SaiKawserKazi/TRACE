"""Unit and smoke tests for TRACE core functions (synthetic data only; no downloads needed).
Run: python tests/test_core.py   (or: python -m pytest tests)"""
import os, sys
import numpy as np, pandas as pd, networkx as nx
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import trace_analysis as T
import aging

rng = np.random.default_rng(0)


def test_score_signs_and_target_exclusion():
    genes = [f"G{i}" for i in range(60)]
    a = pd.Series(np.linspace(-3, 3, 60), index=genes)
    shift = pd.DataFrame([-a.values, a.values, rng.normal(size=60)], index=["REV", "AGE", "G0"], columns=genes)
    r = T.score(shift, a)
    assert abs(r["REV"] - 1) < 1e-9          # exact reversal of the aging axis -> R = +1
    assert abs(r["AGE"] + 1) < 1e-9          # moving along the aging axis -> R = -1
    exp = -stats.spearmanr(shift.loc["G0", genes[1:]], a[genes[1:]])[0]
    assert abs(r["G0"] - exp) < 1e-9         # own target excluded from the signature


def test_residualise_removes_linear_covariate():
    x = pd.DataFrame({"p": rng.normal(size=500)})
    y = pd.Series(2 * x.p.values + rng.normal(scale=0.1, size=500))
    res = T.residualise(y, x)
    assert abs(np.corrcoef(res, x.p)[0, 1]) < 1e-8


def test_holm_matches_statsmodels():
    from statsmodels.stats.multitest import multipletests
    p = rng.uniform(size=20) ** 3
    assert np.allclose(T.holm(p), multipletests(p, method="holm")[1])


def test_pareto_front():
    df = pd.DataFrame({"a": [1, 2, 3, 3], "b": [3, 2, 1, 0]}, index=list("wxyz"))
    assert set(T.pareto(df, ["a", "b"])) == {"w", "x", "y"}


def test_randomised_network_preserves_degrees():
    G = nx.gnm_random_graph(300, 1500, seed=1)
    nx.set_edge_attributes(G, 1.0, "weight")
    H = T.randomised(G, 0)
    assert all(G.degree(n) == H.degree(n) for n in G)
    assert len(set(map(frozenset, G.edges())) & set(map(frozenset, H.edges()))) < 0.5 * G.number_of_edges()


def test_robust_score_definition():
    per = {"K562": pd.DataFrame({"R_adj": [0.2, -0.1]}, index=["A", "B"]),
           "K562gw": pd.DataFrame({"R_adj": [0.4, np.nan]}, index=["A", "B"]),
           "RPE1": pd.DataFrame({"R_adj": [0.1, 0.3]}, index=["A", "B"])}
    r = T.robust(per)
    k562 = 0.3  # mean of the two K562 screens
    exp = np.mean([k562, 0.1]) / (1 + np.std([k562, 0.1]))
    assert abs(r.loc["A", "R_rob"] - exp) < 1e-12 and r.loc["A", "n_lines"] == 2


def test_age_model_recovers_planted_effect():
    n, g = 120, 300
    age = rng.uniform(20, 90, n); sex = rng.choice(["M", "F"], n); batch = rng.choice(["b1", "b2", "b3"], n)
    Y = rng.normal(5, 1, (g, n)); Y[:20] += 0.05 * age; Y[20:40] += 2.0 * (batch == "b2")
    e = pd.DataFrame(Y, index=[f"g{i}" for i in range(g)], columns=[f"s{i}" for i in range(n)])
    m = pd.DataFrame({"age": age, "sex": sex, "batch": batch}, index=e.columns)
    t = aging.age_tstats_cov(e, m, ["sex", "batch"])
    assert (t.q.iloc[:20] < 0.05).all() and (t.q.iloc[20:] < 0.05).mean() < 0.1
    t0 = aging.age_tstats(e, m)
    assert np.allclose(t0.t, aging.age_tstats_cov(e, m, ["sex"]).t)


def test_generic_table_loader():
    import tempfile
    n = 40; age = rng.uniform(20, 90, n)
    raw = pd.DataFrame(rng.poisson(20, (30, n)), index=[f"g{i % 25}" for i in range(30)], columns=[f"s{i}" for i in range(n)])
    meta = pd.DataFrame({"sample": raw.columns, "age": age, "sex": rng.choice(["M", "F"], n)}); meta.loc[0, "sex"] = None
    d = tempfile.mkdtemp(); raw.to_csv(os.path.join(d, "e.csv")); meta.to_csv(os.path.join(d, "m.csv"), index=False)
    e, m = aging.load_table(os.path.join(d, "e.csv"), os.path.join(d, "m.csv"))
    assert e.shape == (25, n - 1) and list(e.columns) == list(m.index) and np.allclose(e.loc["g0"], np.log2(raw.loc["g0"].sum()[m.index] + 1))


if __name__ == "__main__":
    fs = [v for k, v in dict(globals()).items() if k.startswith("test_")]
    for f in fs:
        f(); print("PASS", f.__name__)
    print(f"{len(fs)} tests passed")
