"""TRACE supplementary figures S1-S10 in final numbering (order of first citation). Built at 174 mm (or 114 mm), panel letters
on every multi-panel figure, true minus signs, source data per panel in figures/source_data/."""
import json, sys
import numpy as np, pandas as pd
from scipy import stats
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from figstyle import *  # noqa: F401,F403
from make_figures import J, CSV, R, H2, SC, SENS, BC, F0, G, scores, bcrow, src, PRED, hanley_se

MINUS = lambda v, d=1: f"{v:.{d}f}".replace("-", "−")


# S1 cohort controls
def S1():
    fig = plt.figure(figsize=(W2, 120 * MM)); gs = fig.add_gridspec(2, 3, hspace=0.75, wspace=0.6, left=0.08, right=0.98, top=0.93, bottom=0.1)
    PM, SUB, C3 = CSV("rev3_cohort_permutations.csv"), CSV("rev3_cohort_power_subsets.csv").iloc[:, 0], J("rev3_cohort.json"); fd = C3["fdr05"]
    ax = a1 = fig.add_subplot(gs[0, 0]); v = PM.fdr05_discovery
    bins = np.r_[0, 0.5, np.logspace(0, 4, 30)]; ax.hist(v.clip(lower=0) + 0.0, bins=bins, color=LIGHT)
    ax.set_xscale("symlog", linthresh=1); ax.set_xticks([0, 10, 100, 1000, 10000]); ax.set_xticklabels(["0", "10", "100", "1,000", "10,000"])
    obs = fd[F0]; ax.axvline(obs, color=INK, lw=1)
    ax.text(obs * 0.8, ax.get_ylim()[1] * 0.92, f"observed\n{obs:,}", ha="right", va="top", fontsize=6.5)
    ax.set_xlabel("FDR < 0.05 age genes, Discovery"); ax.set_ylabel("Permutations"); title(ax, f"Permuted age (n = {len(PM):,})")
    ax = a2 = fig.add_subplot(gs[0, 1])
    for col, lab in [("rho_fibroblasts", "GTEx fibroblasts"), ("rho_skin", "GTEx skin")]:
        ax.hist(PM[col], bins=40, histtype="step", color=COHORT_C[lab], lw=1.0, label=lab)
    ax.set_xlabel("Cross-cohort ρ, permuted age"); ax.set_ylabel("Permutations"); ax.legend(loc="upper left", fontsize=6.5); title(ax, "Null for agreement")
    ax = a3 = fig.add_subplot(gs[0, 2]); obs20 = fd["F1 ages 20-79 (GTEx range)"]
    ax.hist(SUB, bins=40, color=LIGHT); ax.axvline(obs20, color=INK, lw=1)
    ax.text(obs20 + 150, ax.get_ylim()[1] * 0.95, f"observed, ages 20–79: {obs20}", ha="left", va="top", fontsize=6.5)
    ax.set_xlabel("FDR < 0.05 genes, 74 random donors"); ax.set_ylabel("Subsets"); title(ax, f"Power-matched (n = {len(SUB):,})")
    ax = a4 = fig.add_subplot(gs[1, 0:2]); CC = CSV("rev3_cohort_cellcycle_t.csv"); pos = 0; ticks = []
    for coh in ["Discovery", "GTEx fibroblasts", "GTEx skin"]:
        for st in ["E2F Targets", "G2-M Checkpoint"]:
            vv = CC[(CC.cohort == coh) & (CC.set == st)].t.values
            ax.boxplot(vv, positions=[pos], widths=0.55, showfliers=False, medianprops=dict(color=INK, lw=1), boxprops=dict(lw=0.6), whiskerprops=dict(lw=0.6), capprops=dict(lw=0))
            ax.scatter(pos + np.random.default_rng(int(pos * 10)).uniform(-0.18, 0.18, len(vv)), vv, s=1.2, color=COHORT_C[coh], lw=0, alpha=0.6, zorder=3)
            ticks.append((pos, f"{coh}\n{'E2F targets' if st.startswith('E2F') else 'G2M checkpoint'}\nn = {len(vv)}")); pos += 1
        pos += 0.5
    ax.axhline(0, **REF); ax.set_xticks([t for t, _ in ticks]); ax.set_xticklabels([l for _, l in ticks]); ax.set_ylabel("Age t-statistic")
    title(ax, "Cell-cycle gene age effects")
    ax = a5 = fig.add_subplot(gs[1, 2]); K = CSV("rev3_cohort_known_genes.csv"); K = K[K.set == "age"].set_index("gene")
    M = K[["t_Discovery", "t_GTEx fibroblasts", "t_GTEx skin"]]
    im = ax.imshow(M.values, cmap=CORR_CMAP, vmin=-8, vmax=8, aspect="auto")
    for (i, j), v in np.ndenumerate(M.values):
        if np.isfinite(v): ax.text(j, i, MINUS(v), ha="center", va="center", fontsize=6.5, color="white" if abs(v) > 5 else INK)
        else:
            ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fc="#d9d9d9", ec="white")); ax.text(j, i, "n.e.", ha="center", va="center", fontsize=6.5, color=GRAY)
    ax.set_xticks(range(3)); ax.set_xticklabels(["Disc.", "Fibro.", "Skin"]); ax.set_yticks(range(len(M))); ax.set_yticklabels([f"$\\it{{{g}}}$" for g in M.index])
    ax.spines[:].set_visible(False); cbar(fig, im, ax, "Age t", ticks=[-8, 0, 8], frac=0.08, ticklabels=["−8", "0", "8"]); title(ax, "Pre-named genes")
    src(PM, "FigS1A", index=False); src(SUB.to_frame("fdr05"), "FigS1C", index=False); src(CC, "FigS1D", index=False); src(M, "FigS1E")
    letters(fig, [[(a1, "A"), (a2, "B"), (a3, "C")], [(a4, "D"), (a5, "E")]]); save(fig, "FigS1")


# S2 sign agreement across all specifications (bias-corrected intervals)
def S2():
    specs = list(dict.fromkeys(SENS.discovery_spec))
    fig, axs = plt.subplots(1, 2, figsize=(W2, 112 * MM), sharey=True, gridspec_kw={"wspace": 0.08, "left": 0.2, "right": 0.98, "top": 0.8, "bottom": 0.11})
    shades = {"GTEx fibroblasts": ["#e8c3d8", "#CC79A7", "#8c4f74", "#4a2a3d"], "GTEx skin": ["#f6d58f", "#E69F00", "#9c6c00", "#4f3700"]}
    rows = []
    for ax, coh, key in [(axs[0], "GTEx fibroblasts", "fibroblasts"), (axs[1], "GTEx skin", "skin")]:
        reps = [r for r in dict.fromkeys(SENS.replication_spec) if r.startswith(key)]
        for k, rp in enumerate(reps):
            d = BC[BC.replication_spec == rp].set_index("discovery_spec").reindex(specs); yy = np.arange(len(specs)) + (k - 1.5) * 0.19
            ax.errorbar(d.sig_sign_agree_fixed, yy, xerr=EB(d.sig_sign_agree_fixed, d.sig_sign_agree_fixed_bc_lo, d.sig_sign_agree_fixed_bc_hi), fmt="o", ms=2.6,
                        color=shades[coh][k], elinewidth=0.75, label=rp.split(": ", 1)[1].replace("G2 G1", "G2: G1").replace("G1 +", "G1:").replace("G0 pre-specified (age + sex)", "G0: age + sex").replace("G3 age + sex + 5 hidden factors", "G3: + 5 hidden factors"))
            rows += [{"tissue": coh, "replication": rp, "discovery": f_, "sign_agree": v, "bc_lo": l, "bc_hi": h} for f_, v, l, h in zip(specs, d.sig_sign_agree_fixed, d.sig_sign_agree_fixed_bc_lo, d.sig_sign_agree_fixed_bc_hi)]
        ax.axvline(0.5, **REF); ax.set_xlim(0.1, 0.85); title(ax, coh)
        ax.set_xlabel("Fraction of signature genes with the\nsame age direction (95% CI)")
        ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.06), ncol=1, fontsize=6.5, frameon=False, handletextpad=0.2, labelspacing=0.25)
    axs[0].set_yticks(range(len(specs))); axs[0].set_yticklabels([FSHORT[s[:2]] for s in specs]); axs[0].invert_yaxis()
    src(pd.DataFrame(rows), "FigS2", index=False)
    letters(fig, [[(axs[0], "A"), (axs[1], "B")]]); save(fig, "FigS2")


# S3 score distributions
def S3():
    fig, axs = plt.subplots(1, 5, figsize=(W2, 66 * MM), sharey=True, gridspec_kw={"wspace": 0.12, "left": 0.07, "right": 0.99, "top": 0.74, "bottom": 0.27})
    ntm = CSV("rev3_reliability_nt_sizematched.csv"); ntm.columns = ["screen", "i", "R"]
    nc = R["negctrl"]; bins = np.linspace(-0.75, 0.45, 50); rows = []
    for ax, c in zip(axs, SC):
        s = scores(c); ax.hist(s.R_raw, bins=bins, color=SCREEN_C[c], density=True)
        xx = np.linspace(-0.75, 0.45, 200); ax.plot(xx, stats.norm.pdf(xx, nc[c]["nt_mean"], nc[c]["nt_sd"]), color=INK, lw=0.9)
        if c in ("HepG2", "Jurkat"): ax.hist(ntm[ntm.screen == c].R.values, bins=bins, histtype="step", density=True, color=GRAY, lw=0.9)
        ax.set_title(f"{SCREEN_L[c]} (n = {len(s):,})\nSD, knockdowns: {fmt(nc[c]['pert_sd'])}\nSD, non-targeting: {fmt(nc[c]['nt_sd'])}", fontsize=6.5, loc="left")
        ax.set_xticks([-0.6, -0.3, 0, 0.3]); ax.set_xticklabels(["−0.6", "−0.3", "0", "0.3"])
        if c == "RPE1": ax.set_xlabel("Raw score R")
        rows.append({"screen": c, "n": len(s), "sd_knockdowns": nc[c]["pert_sd"], "sd_nt": nc[c]["nt_sd"]})
    axs[0].set_ylabel("Density")
    fig.legend([Rectangle((0, 0), 1, 1, fc=LIGHT), Line2D([], [], color=INK), Line2D([], [], color=GRAY)],
               ["knockdowns (screen color)", "non-targeting, normal fit", "non-targeting, size-matched (HepG2, Jurkat)"], loc="lower center", ncol=3, fontsize=6.5, frameon=False)
    src(pd.DataFrame(rows), "FigS3", index=False); letters(fig, [[(a, L) for a, L in zip(axs, "ABCDE")]]); save(fig, "FigS3")


# S4 raw score versus stress
def S4():
    fig, axs = plt.subplots(1, 5, figsize=(W2, 58 * MM), sharey=True, gridspec_kw={"wspace": 0.12, "left": 0.07, "right": 0.99, "top": 0.8, "bottom": 0.2})
    for ax, c in zip(axs, SC):
        s = scores(c); xmax = s.stress.quantile(0.995); s = s[s.stress <= xmax]
        ax.scatter(s.stress, s.R_raw, s=1, c=SCREEN_C[c], lw=0, alpha=0.35, rasterized=True)
        q = pd.qcut(s.stress.rank(method="first"), 10, labels=False); g = s.groupby(q.values).agg(x=("stress", "median"), y=("R_raw", "median"))
        ax.plot(g.x, g.y, color=INK, lw=1.0, marker="o", ms=1.6)
        ax.set_title(f"{SCREEN_L[c]}\nρ = {fmt(H2[c]['rho_raw_stress'])}, n = {len(s):,}", fontsize=6.5, loc="left"); ax.set_ylim(-0.7, 0.45)
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f"{v:.1f}".replace("-", "−")))
        src(s[["stress", "R_raw"]], f"FigS4_{c}")
    fig.supxlabel("Stress/inflammation shift (mean z; x clipped at the 99.5th percentile)", fontsize=6.5, y=0.04)
    axs[0].set_ylabel("Raw score R"); letters(fig, [[(a, L) for a, L in zip(axs, "ABCDE")]]); save(fig, "FigS4")


# S5 proliferation dependence by decile (V0 on its own row)
def S5():
    dec = CSV("rev3_scores_prolif_deciles.csv"); dep = CSV("rev3_scores_prolif_dependence.csv")
    fig, axs = plt.subplots(2, 5, figsize=(W2, 92 * MM), sharex=True, gridspec_kw={"hspace": 0.62, "wspace": 0.35, "left": 0.08, "right": 0.99, "top": 0.92, "bottom": 0.17},
                            height_ratios=[0.8, 1])
    sty = {"V1 pre-specified": ("-", "o", INK), "V6 spline resid.": ("--", "s", GRAY), "V7 resid. (prolif. + identity)": (":", "D", GRAY)}
    for j, c in enumerate(SC):
        d0 = dec[(dec.screen == c) & (dec.variant == "V0 raw")]
        axs[0, j].plot(d0.decile, d0["mean"], color=INK, marker="o", ms=1.8); axs[0, j].fill_between(d0.decile, d0["mean"] - d0["sd"], d0["mean"] + d0["sd"], color="#eeeeee", lw=0)
        axs[0, j].set_title(SCREEN_L[c], fontsize=6.5, loc="left"); axs[0, j].axhline(0, **REF)
        for v, (ls, mk, col) in sty.items():
            d = dec[(dec.screen == c) & (dec.variant == v)]
            axs[1, j].plot(d.decile, d["mean"], color=col, ls=ls, marker=mk, ms=1.8, label=VLABEL[v], mfc="white" if v != "V1 pre-specified" else col)
        axs[1, j].axhline(0, **REF)
        dc = {v: dep[(dep.screen == c) & (dep.variant == v)].distance_corr.iloc[0] for v in sty}
        vs = [f"{VSHORT[v]} {fmt(dc[v])}" for v in sty]
        axs[1, j].set_title("Distance correlation\n" +", ".join(vs), fontsize=6.5, loc="left")
        axs[1, j].set_xticks([1, 5, 10])
    axs[0, 0].set_ylabel("V0 raw: mean (± SD)"); axs[1, 0].set_ylabel("Adjusted: decile mean"); axs[1, 2].set_xlabel("Proliferation decile")
    h, l = axs[1, 0].get_legend_handles_labels(); fig.legend(h, l, loc="lower center", ncol=3, fontsize=6.5, frameon=False)
    src(dec, "FigS5", index=False); letters(fig, [[(axs[0, 0], "A")], [(axs[1, 0], "B")]]); save(fig, "FigS5")


# S6 random-signature nulls
def S6():
    N = CSV("rev3_scores_random_signatures.csv"); real = J("rev3_scores.json")["E5_real"]; dep = CSV("rev3_scores_prolif_dependence.csv")
    panels = [(f"H2_rho_{c}", f"{SCREEN_L[c]}: ρ(score, proliferation)", dep[(dep.screen == c) & (dep.variant == "V2 genes removed")].spearman.iloc[0], "ρ") for c in SC] + \
             [("H1_median_rho_raw", "Cross-line ρ, raw", real["H1_median_rho_raw"], "median ρ"), ("H1_median_rho_adj", f"Cross-line ρ, {RADJ}", real["H1_median_rho_adj"], "median ρ"),
              ("H4_CellAge_auroc_raw", f"CellAge, raw {RROB}", real["H4_CellAge_auroc_raw"], "AUROC"), ("H4_CellAge_auroc_adj", "CellAge, V1", real["H4_CellAge_auroc_adj"], "AUROC"),
              ("H4_GenAge_auroc_raw", f"GenAge, raw {RROB}", real["H4_GenAge_auroc_raw"], "AUROC"), ("H4_GenAge_auroc_adj", "GenAge, V1", real["H4_GenAge_auroc_adj"], "AUROC")]
    fig, axs = plt.subplots(3, 4, figsize=(W2, 120 * MM), gridspec_kw={"hspace": 0.9, "wspace": 0.45, "left": 0.07, "right": 0.99, "top": 0.94, "bottom": 0.08})
    rows = []
    for ax, (k, lab, obs, xl), L in zip(axs.ravel(), panels, "ABCDEFGHIJK"):
        ax.hist(N[k], bins=25, color=LIGHT); ax.axvline(obs, color=INK, lw=1.4)
        ax.set_title(lab, fontsize=6.5, loc="left"); ax.set_xlabel(xl); ax.set_ylabel("Signatures")
        ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True, nbins=4))
        rows.append({"statistic": k, "observed": obs, "null_median": N[k].median(), "frac_null_ge_obs": float((N[k] >= obs).mean())})
    axs.ravel()[-1].axis("off")
    axs.ravel()[-1].legend([Rectangle((0, 0), 1, 1, fc=LIGHT), Line2D([], [], color=INK, lw=1.4)], [f"{len(N)} random signatures", "real signature\n(H2 panels: V2, without\nproliferation genes)"],
                           loc="center", fontsize=6.5)
    src(pd.DataFrame(rows), "FigS6", index=False)
    letters(fig, [[(axs[r, c], "ABCDEFGHIJK"[r * 4 + c]) for c in range(4) if r * 4 + c < 11] for r in range(3)]); save(fig, "FigS6")


# S7 GTEx-derived signature: full comparison
def S7():
    T1 = CSV("x_gtex_signature_h1_h5.csv"); PR = CSV("x_gtex_signature_programs.csv"); EV = CSV("x_gtex_signature_evidence.csv"); EVd = CSV("rev3_scores_evidence.csv")
    fig = plt.figure(figsize=(W2, 150 * MM)); gt = fig.add_gridspec(1, 2, wspace=0.75, left=0.2, right=0.98, top=0.87, bottom=0.6); gs = [gt[0], None, gt[1]]
    G = [None, fig.add_gridspec(1, 1, left=0.04, right=0.98, top=0.44, bottom=0.1)[0]]
    R0, R1 = R, J("e18_gtex_discovery/results.json")
    rho_rows = [("Cross-line ρ, raw (median)", np.median([v["rho_raw"] for v in R0["H1"]["pairs"].values()]), np.median([v["rho_raw"] for v in R1["H1"]["pairs"].values()])),
                (f"Cross-line ρ, {RADJ} (median)", np.median([v["rho_adj"] for v in R0["H1"]["pairs"].values()]), np.median([v["rho_adj"] for v in R1["H1"]["pairs"].values()])),
                ("Within K562 (H1b)", R0["H1b"]["rho_within_K562"], R1["H1b"]["rho_within_K562"]),
                ("Donor bootstrap median (H5a)", R0["H5"]["median_rho"], R1["H5"]["median_rho"])]
    ax = a = fig.add_subplot(gs[0]); y = np.arange(len(rho_rows))
    ax.scatter([r[1] for r in rho_rows], y - 0.12, marker="o", color=COHORT_C["Discovery"], s=16, label="Discovery signature", zorder=3)
    ax.scatter([r[2] for r in rho_rows], y + 0.12, marker="s", color=COHORT_C["GTEx fibroblasts"], s=16, label="GTEx-fibroblast signature", zorder=3)
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rho_rows]); ax.invert_yaxis(); ax.set_xlim(0, 1); ax.set_xlabel("Spearman ρ"); gridx(ax)
    ax.legend(loc="lower left", bbox_to_anchor=(-0.05, 1.03), fontsize=6.5, frameon=False, ncol=1); title(ax, " ")
    ax = b = fig.add_subplot(G[1].subgridspec(1, 3, width_ratios=[0.25, 1, 0.08])[1]); progs = ["prolif", "stress", "oxphos", "magnitude", "identity"]
    M = PR.pivot_table(index=["program"], columns=["signature", "screen"], values="rho")
    mat = np.array([[PR[(PR.signature == sg) & (PR.program == p) & (PR.screen == c)].rho.iloc[0] for sg in ["Discovery signature", "GTEx-fibroblast signature"] for c in SC] for p in progs])
    im = ax.imshow(mat, cmap=CORR_CMAP, vmin=-1, vmax=1, aspect="auto")
    for (i, j), v in np.ndenumerate(mat): ax.text(j, i, fmt(v, 1), ha="center", va="center", fontsize=6.5, color="white" if abs(v) > 0.55 else INK)
    ax.axvline(4.5, color="white", lw=2.5)
    ax.set_xticks(range(10)); ax.set_xticklabels([SCREEN_L[c] for c in SC] * 2, rotation=45, ha="right", rotation_mode="anchor"); ax.set_yticks(range(5))
    ax.set_yticklabels(["proliferation", "stress", "oxidative\nphosphorylation", "response size", "identity"]); ax.spines[:].set_visible(False)
    ax.text(2, -0.75, "Discovery signature", ha="center", fontsize=6.5); ax.text(7, -0.75, "GTEx-fibroblast signature", ha="center", fontsize=6.5)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02); cb.set_label("Spearman ρ with score"); cb.outline.set_visible(False)
    title(ax, " ")
    ax = c_ = fig.add_subplot(gs[2]); rows = []
    for k, (lab, gsn) in enumerate([("CellAge", "all genes"), ("CellAge", "genes in >= 2 lines"), ("GenAge", "all genes"), ("GenAge", "genes in >= 2 lines")]):
        for j, (sname_d, sname_g, mk) in enumerate([("R_rob (V1, pre-specified)", "R_rob (V1)", "o"), ("proliferation shift alone", "proliferation shift alone", "D")]):
            d = EVd[(EVd.labels == lab) & (EVd.gene_set == gsn) & (EVd.score == sname_d)].iloc[0]; gg = EV[(EV.labels == lab) & (EV.gene_set == gsn) & (EV.score == sname_g)].iloc[0]
            yy = k + (j - 0.5) * 0.4
            for off, row, col in [(-0.08, d, COHORT_C["Discovery"]), (0.08, gg, COHORT_C["GTEx fibroblasts"])]:
                se = hanley_se(row.auroc, int(row.n_pos), int(row.n_neg))
                ax.errorbar(row.auroc, yy + off, xerr=1.96 * se, fmt=mk, ms=3, color=col, mfc=col if j == 0 else "white", elinewidth=0.7)
            rows.append({"labels": lab, "gene_set": gsn, "score": sname_d, "discovery": d.auroc, "gtex": gg.auroc})
    ax.axvline(0.5, **REF); ax.set_yticks(range(4)); ax.set_yticklabels(["CellAge, all", "CellAge, ≥2 lines", "GenAge, all", "GenAge, ≥2 lines"]); ax.invert_yaxis()
    ax.set_xlabel("AUROC (95% CI)"); ax.set_xlim(0.3, 0.8)
    ax.legend([Line2D([], [], marker="o", ls="none", color=INK), Line2D([], [], marker="D", ls="none", mfc="white", mec=INK),
               Line2D([], [], color=COHORT_C["Discovery"]), Line2D([], [], color=COHORT_C["GTEx fibroblasts"])],
              [f"{RROB} (V1)", "proliferation only", "discovery signature", "GTEx-fibroblast signature"],
              loc="lower center", bbox_to_anchor=(0.5, 1.03), ncol=2, fontsize=6.5, frameon=False); title(ax, " ")
    src(PR, "FigS7C", index=False); src(pd.DataFrame(rows), "FigS7B", index=False); src(pd.DataFrame(rho_rows, columns=["statistic", "discovery", "gtex"]), "FigS7A", index=False)
    letters(fig, [[(a, "A"), (c_, "B")], [(b, "C")]]); save(fig, "FigS7")


# S8 prediction error of every baseline, relative to the constant predictor
def S8():
    NS = CSV("rev3_network_summary.csv"); x = np.arange(5)
    fig, ax = plt.subplots(figsize=(W2, 72 * MM)); fig.subplots_adjust(left=0.17, right=0.98, top=0.8, bottom=0.12)
    base = NS[NS.predictor == "constant"].set_index("screen").reindex(SC).mae; rows = []
    for k, (p, lab, mk, fc, ec) in enumerate(PRED):
        if p == "constant" or p not in set(NS.predictor): continue
        d = NS[NS.predictor == p].set_index("screen").reindex(SC); v = 100 * (d.mae / base - 1); e = 100 * d.mae_sd / base
        ax.errorbar(v, x + (k - 4) * 0.075, xerr=e, fmt=mk, ms=3.2, mfc=fc, mec=ec, ecolor=ec, elinewidth=0.6, label=lab)
        rows += [{"screen": c, "predictor": p, "mae_rel_constant_pct": vv} for c, vv in zip(SC, v)]
    ax.axvline(0, **REF); ax.set_yticks(x); ax.set_yticklabels([f"{SCREEN_L[c]} (n = {int(NS[NS.screen == c].n.iloc[0]):,})" for c in SC]); ax.invert_yaxis()
    ax.set_xlabel("Mean absolute error relative to the constant predictor (%; mean ± SD over 10 seeds)"); gridx(ax)
    ax.legend(loc="lower center", bbox_to_anchor=(0.45, 1.02), ncol=4, fontsize=6.5, frameon=False, columnspacing=0.6, handletextpad=0.2)
    src(pd.DataFrame(rows), "FigS8", index=False); save(fig, "FigS8")


# S9 single-screen evidence tests
def S9():
    d = R["H4"]; fig, axs = plt.subplots(1, 2, figsize=(W15, 66 * MM), sharey=True, gridspec_kw={"wspace": 0.45, "left": 0.24, "right": 0.97, "top": 0.8, "bottom": 0.17})
    rows = []
    for ax, key, lab in [(axs[0], "H4a_cellage", "CellAge"), (axs[1], "H4b_genage", "GenAge")]:
        for i, c in enumerate(SC):
            for typ, mk, off in [("raw", "s", -0.15), ("adj", "o", 0.15)]:
                e = d[key][f"{typ}_{c}"]; A_, n1, n2 = e["auroc"], e["n_pos"], e["n_neg"]; se = hanley_se(A_, n1, n2)
                ax.errorbar(A_, i + off, xerr=1.96 * se, fmt=mk, ms=3, mfc="white" if typ == "raw" else INK, mec=INK, ecolor=GRAY, elinewidth=0.6)
                rows.append({"labels": lab, "screen": c, "score": typ, "auroc": A_, "lo": A_ - 1.96 * se, "hi": A_ + 1.96 * se, "n_pos": n1, "n_neg": n2})
            ax.text(-0.03, i, f"{n1}/{n2}", transform=ax.get_yaxis_transform(), fontsize=6.5, va="center", ha="right", color=GRAY)
        ax.text(-0.03, -0.75, "n (+/−)", transform=ax.get_yaxis_transform(), fontsize=6.5, va="center", ha="right", color=GRAY)
        ax.axvline(0.5, **REF); ax.set_xlim(0.3, 0.75); ax.set_xlabel("AUROC (95% CI)"); title(ax, lab)
    axs[0].set_yticks(range(5)); axs[0].set_yticklabels([SCREEN_L[c] for c in SC]); axs[0].invert_yaxis(); axs[0].tick_params(axis="y", pad=34)
    fig.legend([Line2D([], [], marker="s", ls="none", mfc="white", mec=INK), Line2D([], [], marker="o", ls="none", color=INK)], ["raw R", f"{RADJ} (V1)"],
               loc="upper center", ncol=2, fontsize=6.5, frameon=False)
    src(pd.DataFrame(rows), "FigS9", index=False); letters(fig, [[(axs[0], "A"), (axs[1], "B")]]); save(fig, "FigS9")


# S10 Pareto front before and after filtering
def S10():
    P = CSV("rev3_scores_pareto_filtered.csv", index_col=0)
    fig, axs = plt.subplots(1, 2, figsize=(W2, 76 * MM), sharey=True, sharex=True, gridspec_kw={"wspace": 0.08, "left": 0.08, "right": 0.98, "top": 0.86, "bottom": 0.16})
    for ax, (lab, pcol, sub) in zip(axs, [("Unfiltered", "pareto", P), ("Filtered", "pareto_filtered", P[P.passes_filter])]):
        ax.scatter(sub.prolif_shift, sub.R_rob, s=2, color=LIGHT, lw=0, rasterized=True)
        F = sub[sub[pcol].astype(bool)]; weak = F.flag_weak_or_unmeasured_knockdown.astype(bool)
        ax.scatter(F.prolif_shift[~weak], F.R_rob[~weak], s=11, color=INK, lw=0); ax.scatter(F.prolif_shift[weak], F.R_rob[weak], s=11, facecolors="white", edgecolors=INK, lw=0.7)
        ax.set_title(f"{lab}: {len(sub):,} genes, {len(F)} on the front", fontsize=6.5, loc="left")
        ax.axhline(0, color=LIGHT, lw=0.4, zorder=0); ax.axvline(0, color=LIGHT, lw=0.4, zorder=0); ax.set_xlim(-0.33, 0.16)
        if "ANLN" in F.index:
            ax.annotate("$\\it{ANLN}$", (F.prolif_shift["ANLN"], F.R_rob["ANLN"]), xytext=(F.prolif_shift["ANLN"] + 0.06, F.R_rob["ANLN"] - 0.02), fontsize=6.5,
                        arrowprops=dict(arrowstyle="-", lw=0.5, color=INK), va="center")
        ax.set_xlabel("Proliferation shift (mean z over lines)")
    axs[0].set_ylabel(RROB)
    fig.legend([Line2D([], [], marker="o", ls="none", color=INK), Line2D([], [], marker="o", ls="none", mfc="white", mec=INK)], ["Pareto-optimal", "Pareto-optimal, weak or unmeasured on-target knockdown"],
               loc="upper center", ncol=2, fontsize=6.5, frameon=False)
    src(P, "FigS10"); letters(fig, [[(axs[0], "A"), (axs[1], "B")]]); save(fig, "FigS10")


if __name__ == "__main__":
    for name in (sys.argv[1:] or ["S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10"]):
        globals()[name](); print("ok", name, flush=True)
