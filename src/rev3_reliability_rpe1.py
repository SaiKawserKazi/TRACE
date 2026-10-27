"""E6 for RPE1 (exploratory, added in revision; logged 2026-10-27 before running): split-cell reliability from the Replogle
RPE1 normalized single-cell release (Figshare 20029387, file 35775554), same 20-split design as HepG2 and Jurkat. Appends RPE1
rows to the reliability outputs and adds the RPE1 pairs to the disattenuated table."""
import os, json
import numpy as np, pandas as pd
from scipy import stats

import aging, perturb
import trace_analysis as T
import rev3_reliability as RL
from rev3_scores import variants

RES = RL.RES


def main():
    e, m = aging.load_fleischer(); sig_f = aging.signature(aging.age_tstats(e, m))
    ctx = {c: perturb.load_context(c) for c in perturb.CONTEXTS}
    ident = T.identity_genes(ctx)
    keep = set(sig_f.index) | T.PROLIF_REMOVE | T.PROLIF | T.STRESS | set(ident["RPE1"])
    r, h, sizes, nts = RL.split_line("RPE1", sig_f.t, ident["RPE1"], keep, np.random.default_rng(71_000),
                                     fname="rpe1_normalized_singlecell_01.h5ad", normalized=True)
    RL.log("rpe1 done", sizes)
    sp = pd.read_csv(os.path.join(RES, "rev3_reliability_splits.csv")); sp = pd.concat([sp[sp.screen != "RPE1"], r]); sp.to_csv(os.path.join(RES, "rev3_reliability_splits.csv"), index=False)
    hp = pd.read_csv(os.path.join(RES, "rev3_reliability_heldout_prolif.csv")); hp = pd.concat([hp[hp.screen != "RPE1"], h]); hp.to_csv(os.path.join(RES, "rev3_reliability_heldout_prolif.csv"), index=False)
    agg = r.groupby(["screen", "variant"]).r_half.agg(["mean", "std"]).reset_index(); agg["reliability_SB"] = 2 * agg["mean"] / (1 + agg["mean"])
    REL = pd.read_csv(os.path.join(RES, "rev3_reliability_summary.csv")); REL = pd.concat([REL[REL.screen != "RPE1"], agg], ignore_index=True)
    REL.to_csv(os.path.join(RES, "rev3_reliability_summary.csv"), index=False)
    relmap = {(x.screen, x.variant): x.reliability_SB for x in REL.itertuples()}
    progs = {c: T.program_shifts(v[0], ident[c]) for c, v in ctx.items()}
    Vd = variants(ctx, sig_f.t, progs); dis = []
    for v in RL.VAR:
        for a_, b_ in [("K562", "RPE1"), ("K562", "HepG2"), ("K562", "Jurkat"), ("RPE1", "HepG2"), ("RPE1", "Jurkat"), ("HepG2", "Jurkat")]:
            sa, sb = Vd[v][a_], Vd[v][b_]; sh = sa.index.intersection(sb.index); rho = stats.spearmanr(sa[sh], sb[sh])[0]
            ra, rb = relmap[(a_, v)], relmap[(b_, v)]
            dis.append({"variant": v, "pair": f"{a_}-{b_}", "n": len(sh), "rho": rho, "rel_a": ra, "rel_b": rb,
                        "rho_disattenuated": rho / np.sqrt(ra * rb) if ra > 0 and rb > 0 else np.nan})
    pd.DataFrame(dis).to_csv(os.path.join(RES, "rev3_reliability_disattenuated.csv"), index=False)
    J = json.load(open(os.path.join(RES, "rev3_reliability.json"))); J["RPE1"] = sizes
    J["heldout_prolif"]["RPE1"] = h[["rho_R_B_prolif_A", "rho_R_B_prolif_B"]].mean().round(3).to_dict()
    J["rpe1_source"] = "Replogle et al. 2022, Figshare 20029387, rpe1_normalized_singlecell_01.h5ad (md5 2c36a053960f3fae157adacdbccd4485)"
    json.dump(J, open(os.path.join(RES, "rev3_reliability.json"), "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    print(agg.round(3).to_string())


if __name__ == "__main__":
    main()
