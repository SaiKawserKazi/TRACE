"""Random-signature null, real side without proliferation genes (exploratory, Amendment 5).

The 100 random signatures are drawn from genes outside the proliferation sets, so the fair real comparator for every
raw-score statistic is the real signature with its proliferation genes removed (V2), not the full signature (V0).
This recomputes the real H1, H2 and H4 summary statistics from V2 (raw) and V1 (adjusted) with the same function used
for the random signatures, and writes empirical P values against the existing null:
results/x_random_signature_v2.csv."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import aging, perturb
import trace_analysis as T
import rev3_scores as S
from config import RES


def main():
    e, m = aging.load_fleischer(); ts = aging.age_tstats(e, m); sig_f = aging.signature(ts)
    ctx = {c: perturb.load_context(c) for c in S.CTX}
    ident = T.identity_genes(ctx); progs = {c: T.program_shifts(v[0], ident[c]) for c, v in ctx.items()}
    labs = S.labels()
    Vd = S.variants(ctx, sig_f.t, progs)
    real_v2 = S.summarize(Vd["V2 genes removed"], Vd["V1 pre-specified"], progs, labs)
    real_v0 = S.summarize(Vd["V0 raw"], Vd["V1 pre-specified"], progs, labs)
    N = pd.read_csv(RES / "rev3_scores_random_signatures.csv")
    rows = []
    for k, v in real_v2.items():
        null = N[k]
        side = "upper"   # expected direction for every statistic: more agreement, more coupling, higher AUROC
        p = (1 + (null >= v).sum()) / (1 + len(null))
        rows.append({"statistic": k, "real_V2_input": v, "real_full_signature_V0": real_v0[k], "null_median": null.median(),
                     "null_q025": null.quantile(0.025), "null_q975": null.quantile(0.975), "null_max": null.max(), "null_min": null.min(),
                     "frac_null_ge_real": float((null >= v).mean()), "empirical_p_one_sided": float(p), "tail": side, "n_random": len(null)})
    out = pd.DataFrame(rows); out.to_csv(RES / "x_random_signature_v2.csv", index=False)
    print(out.round(3).to_string())


if __name__ == "__main__":
    main()
