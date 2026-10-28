"""Bias-corrected (BC) bootstrap intervals for Table S3 from the saved donor bootstraps (exploratory).
BC adjusts percentile intervals for the median bias of the bootstrap distribution (Efron 1987); the acceleration term of
BCa would need a jackknife refit of every model per donor and is not computed. Output: results/rev3_cohort_sensitivity_bc.csv"""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
import config as _cfg
import os
import numpy as np, pandas as pd
from scipy import stats

RES = str(_cfg.RES)
S = pd.read_csv(os.path.join(RES, "rev3_cohort_sensitivity.csv"))
BT = pd.read_csv(os.path.join(RES, "rev3_cohort_bootstrap.csv.gz"))
METRICS = ["rho_all", "rho_no_cellcycle", "sig_sign_agree_fixed", "top_rho", "sex_top_rho", "sex_top_autosomal_rho", "sex_rho_all"]
za = stats.norm.ppf([0.025, 0.975])
rows = []
for _, r in S.iterrows():
    B = BT[(BT.fn == r.discovery_spec) & (BT.gn == r.replication_spec)]
    out = {"discovery_spec": r.discovery_spec, "replication_spec": r.replication_spec, "n_boot": len(B)}
    for m in METRICS:
        b = B[m].dropna().values; th = r[m]
        p = np.clip((b < th).mean() + 0.5 * (b == th).mean(), 1 / (len(b) + 1), len(b) / (len(b) + 1))
        z0 = stats.norm.ppf(p); q = stats.norm.cdf(2 * z0 + za)
        out.update({m: th, f"{m}_boot_mean": b.mean(), f"{m}_pct_lo": np.quantile(b, 0.025), f"{m}_pct_hi": np.quantile(b, 0.975),
                    f"{m}_bc_lo": np.quantile(b, q[0]), f"{m}_bc_hi": np.quantile(b, q[1]), f"{m}_z0": z0,
                    f"{m}_estimate_outside_pct_ci": not (np.quantile(b, 0.025) <= th <= np.quantile(b, 0.975))})
    rows.append(out)
D = pd.DataFrame(rows); D.to_csv(os.path.join(RES, "rev3_cohort_sensitivity_bc.csv"), index=False)
for m in METRICS:
    print(f"{m}: estimate outside percentile CI in {int(D[f'{m}_estimate_outside_pct_ci'].sum())} of {len(D)} rows; |z0| median {D[f'{m}_z0'].abs().median():.2f}")
k = D[D.discovery_spec.str.startswith("F0") & D.replication_spec.str.contains("G0")]
print(k[["replication_spec", "top_rho", "top_rho_pct_lo", "top_rho_pct_hi", "top_rho_bc_lo", "top_rho_bc_hi", "sex_top_rho", "sex_top_rho_bc_lo", "sex_top_rho_bc_hi",
         "rho_all", "rho_all_bc_lo", "rho_all_bc_hi"]].round(3).T.to_string())
