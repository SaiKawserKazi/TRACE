"""Cross-cohort agreement among the k strongest discovery age genes, k = 50 to 1,000 (exploratory, Amendment 5).

Uses the point-estimate t-statistics of the pre-specified models (F0 discovery; G0 GTEx fibroblasts and skin) and ranks
discovery genes by |t| among genes shared with the GTEx tissue. Writes results/x_topk_sweep.csv and x_topk_f9_check.json
(F9 t-statistics on the F0 top 200; median |t| of the F0 top 200)."""
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
    # F9 check: hidden factors change standard errors, not effect sizes. Agreement of the F9 t-statistics on the F0 top 200
    # (same genes) separates a change of gene set from a change of the t-statistics themselves. Also the median |t| of the
    # real F0 top 200, the reference size for the planted-effect rule in Box 1.
    f9 = next(k for k in T if k.startswith("F9")); g0 = next(k for k in T if k.startswith("fibroblasts: G0")); a, c, b = T[f0], T[f9], T[g0]
    sh = a.index.intersection(b.index).intersection(c.index); top0 = a.t[sh].abs().nlargest(200).index; top9 = c.t[sh].abs().nlargest(200).index
    chk = {"F9_t_on_F0_top200_rho": float(stats.spearmanr(c.t[top0], b.t[top0])[0]), "F0_t_on_F0_top200_rho": float(stats.spearmanr(a.t[top0], b.t[top0])[0]),
           "F9_own_top200_rho": float(stats.spearmanr(c.t[top9], b.t[top9])[0]), "F9_F0_top200_shared_genes": int(len(top0.intersection(top9))),
           "F0_top200_median_abs_t": float(a.t[top0].abs().median()), "F0_top200_min_abs_t": float(a.t[top0].abs().min())}
    import json; json.dump(chk, open(RES / "x_topk_f9_check.json", "w"), indent=1); print(chk)


if __name__ == "__main__":
    main()
