"""Cross-cohort agreement among the k strongest discovery age genes, k = 50 to 1,000 (exploratory, Amendment 5).

Uses the point-estimate t-statistics of the pre-specified models (F0 discovery; G0 GTEx fibroblasts and skin) and ranks
discovery genes by |t| among genes shared with the GTEx tissue. Writes results/x_topk_sweep.csv."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
from scipy import stats
from config import RES

KS = [50, 100, 200, 300, 500, 1000]


def main():
    T = pd.read_pickle(RES / "rev3_point_tstats.pkl")
    f0 = next(k for k in T if k.startswith("F0")); rows = []
    for tissue in ["fibroblasts", "skin"]:
        g0 = next(k for k in T if k.startswith(f"{tissue}: G0")); a, b = T[f0], T[g0]
        sh = a.index.intersection(b.index); top = a.t[sh].abs().sort_values(ascending=False).index
        for k in KS:
            g = top[:k]; r = stats.spearmanr(a.t[g], b.t[g])[0]
            rows.append({"tissue": tissue, "k": k, "spearman_rho": r, "sign_agreement": float(((a.t[g] > 0) == (b.t[g] > 0)).mean())})
    out = pd.DataFrame(rows); out.to_csv(RES / "x_topk_sweep.csv", index=False); print(out.round(3).to_string())


if __name__ == "__main__":
    main()
