"""TRACE figures, final numbering: Figures 1-5, the graphical abstract and Figures S1-S10 (make_figures_supp.py).
Built at 174 mm (main) with figstyle; every data panel writes figures/source_data/<Figure><panel>.csv. Reads results/ only."""
import json, sys
import numpy as np, pandas as pd
from scipy import stats
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from figstyle import *  # noqa: F401,F403
import figstyle as FS
import schematics

J = lambda f: json.load(open(RES / f, encoding="utf-8"))
CSV = lambda f, **k: pd.read_csv(RES / f, **k)
R = J("results.json"); H2 = R["H2"]["contexts"]; SC = SCREENS
SENS, BC = CSV("rev3_cohort_sensitivity.csv"), CSV("rev3_cohort_sensitivity_bc.csv")
F0 = "F0 pre-specified (age 18+; age + sex)"
G = lambda coh, k: next(g for g in SENS.replication_spec.unique() if g.startswith(coh) and g.split(": ")[1].startswith(k))
scores = lambda c: pd.read_csv(RES / f"scores_{c}.csv", index_col=0)
OLDV = {"V0 raw": "V0 raw", "V1 pre-specified": "V1 protocol", "V2 genes removed": "V2 removal only", "V3 resid. (all 3)": "V3 resid only",
        "V4 resid. (prolif.)": "V4 prolif resid", "V5 stratified": "V5 stratified"}


def bcrow(fn, gn): return BC[(BC.discovery_spec == fn) & (BC.replication_spec == gn)].iloc[0]


def srow(fn, gn): return SENS[(SENS.discovery_spec == fn) & (SENS.replication_spec == gn)].iloc[0]


def fisher_ci(r, n):
    z = np.arctanh(r); se = 1 / np.sqrt(n - 3); return np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)


def src(df, name, **k): df.to_csv(SRC / f"{name}.csv", **k)


def mm_ax(fig, x0, top, w, h):
    """Axes placed in millimetres from the top-left corner of the figure (x0, top, width, height)."""
    W, H = fig.get_size_inches() / MM
    return fig.add_axes([x0 / W, 1 - (top + h) / H, w / W, h / H])


# ================================================================== FIGURE 1
def spec_heat(ax, metric, p, pcol, tcol):
    """Discovery specifications (rows) x replication models G0-G2 for both GTEx tissues (columns); outline = G0 cell outside
    the permuted-age null of its own specification."""
    specs = list(dict.fromkeys(SENS.discovery_spec)); colsG = [G("fibroblasts", k) for k in ("G0", "G1", "G2")] + [G("skin", k) for k in ("G0", "G1", "G2")]
    Mf = np.array([[srow(f_, g_)[metric] for g_ in colsG] for f_ in specs]); xpos = [0, 1, 2, 3.4, 4.4, 5.4]
    for jj, xp in enumerate(xpos):
        for ii in range(len(specs)):
            v = Mf[ii, jj]; c = plt.get_cmap(CORR_CMAP)((np.clip(v, -CORR_LIM, CORR_LIM) + CORR_LIM) / (2 * CORR_LIM))
            ax.add_patch(Rectangle((xp - 0.5, ii - 0.5), 1, 1, fc=c, ec="white", lw=0.6))
            ax.text(xp, ii, fmt(v), ha="center", va="center", fontsize=6.5, color=text_on(c))
    for ii, f_ in enumerate(specs):
        for jj, coh in [(0, "fibroblasts"), (3, "skin")]:
            r_ = p[(p.spec == f_) & (p.tissue == coh)]
            if len(r_) and r_[pcol].iloc[0] < 0.05: ax.add_patch(Rectangle((xpos[jj] - 0.45, ii - 0.45), 0.9, 0.9, fill=False, ec="black", lw=1.1))
    for x0, coh in [(1, "fibroblasts"), (4.4, "skin")]:
        ax.text(x0, -0.85, f"GTEx {coh}", ha="center", va="center", fontsize=6.5, color=COHORT_C[f"GTEx {coh}"], fontweight="bold")
    ax.set_xlim(-0.6, 6.0); ax.set_ylim(len(specs) - 0.5, -1.25)
    ax.set_xticks(xpos); ax.set_xticklabels(["G0", "G1", "G2"] * 2); ax.set_yticks(range(len(specs)))
    ax.spines[:].set_visible(False); ax.tick_params(length=0)
    return specs, colsG, Mf


def fig1():
    Hf = 179
    fig = plt.figure(figsize=(W2, Hf * MM))
    A = fig.add_axes([0, 1 - 43 / Hf, 1, 42 / Hf]); schematics.fig1a(A)
    TS = pd.read_pickle(RES / "rev3_point_tstats.pkl"); p = CSV("x_cohort_perm_p.csv")
    # B hexbin, sequential blue scale (log)
    ax = B = mm_ax(fig, 15, 53, 36, 30)
    tf, tg = TS[F0], TS[G("fibroblasts", "G0")]; sh = tf.index.intersection(tg.index)
    hb = ax.hexbin(tf.t[sh], tg.t[sh], gridsize=34, cmap="Blues", bins="log", mincnt=1, linewidths=0, rasterized=True)
    ax.axhline(0, color=LIGHT, lw=0.4, zorder=0); ax.axvline(0, color=LIGHT, lw=0.4, zorder=0)
    pf = p[(p.spec == F0) & (p.tissue == "fibroblasts")].iloc[0]
    ax.set_xlabel("Age t, discovery"); ax.set_ylabel("Age t, GTEx fibroblasts")
    cb = cbar(fig, hb, ax, "Genes per hexagon", frac=0.07); cb.set_ticks([1, 3, 10, 30]); cb.set_ticklabels(["1", "3", "10", "30"])
    title(ax, f"Discovery versus GTEx fibroblasts\nρ = {fmt(pf.rho_all)}, permutation P = {pf.p_rho_all:.2f}")
    src(pd.DataFrame({"t_discovery": tf.t[sh], "t_gtex_fibroblasts": tg.t[sh]}), "Fig1B")
    # C agreement matrix
    ax = Cm = mm_ax(fig, 82, 53, 27, 30)
    RC = J("rev_cohorts.json")["split_half"]; ts = [TS[F0], TS[G("fibroblasts", "G0")], TS[G("skin", "G0")]]
    M = np.full((3, 3), np.nan)
    for i in range(3):
        for j in range(i):
            s2 = ts[i].index.intersection(ts[j].index); M[i, j] = stats.spearmanr(ts[i].t[s2], ts[j].t[s2])[0]
    for i, k in enumerate(["Fleischer split-half", "GTEx fibroblast split-half", "GTEx skin split-half"]): M[i, i] = RC[k]["median"]
    im = ax.imshow(M, cmap=CORR_CMAP, vmin=-CORR_LIM, vmax=CORR_LIM, aspect="auto")
    for (i, j), v in np.ndenumerate(M):
        if np.isfinite(v): heat_text(ax, j, i, v, CORR_CMAP, -CORR_LIM, CORR_LIM)
    for i in range(3): ax.add_patch(Rectangle((i - 0.5, i - 0.5), 1, 1, fill=False, ec="black", lw=0.8))
    nm3 = ["Disc.", "GTEx\nfibro.", "GTEx\nskin"]
    ax.set_xticks(range(3)); ax.set_xticklabels(nm3); ax.set_yticks(range(3)); ax.set_yticklabels(nm3); ax.spines[:].set_visible(False)
    cbar(fig, im, ax, "", ticks=[-0.6, 0, 0.6], frac=0.06, ticklabels=["−0.6", "0", "0.6"], extend="both")
    title(ax, "Age effects,\nall genes")
    src(pd.DataFrame(M, index=["Discovery", "GTEx fibroblasts", "GTEx skin"], columns=["Discovery", "GTEx fibroblasts", "GTEx skin"]), "Fig1C")
    # D 200 strongest genes: age and sex, permuted-age band, planted effect (black triangle); key at the right
    ax = Dp = mm_ax(fig, 122, 53, 24, 30)
    NU = J("x_cohort_nulls.json"); PB = CSV("x_cohort_perm_by_spec.csv"); f0n = PB[PB.spec == F0]; rows = []
    for i, (coh, gk) in enumerate([("fibroblasts", G("fibroblasts", "G0")), ("skin", G("skin", "G0"))]):
        b = bcrow(F0, gk); nn = f0n[f0n.tissue == coh]
        ax.add_patch(Rectangle((nn.top_rho.quantile(0.025), i - 0.36), nn.top_rho.quantile(0.975) - nn.top_rho.quantile(0.025), 0.72, fc="#e6e6e6", ec="none", zorder=0))
        for k, (key, lab, mk, fc) in enumerate([("top_rho", "age", "o", INK), ("sex_top_rho", "sex", "s", "white"), ("sex_top_autosomal_rho", "sex, autosomal", "D", "white")]):
            y = i + (k - 1) * 0.22
            ax.errorbar(b[key], y, xerr=EB([b[key]], [b[f"{key}_bc_lo"]], [b[f"{key}_bc_hi"]]), fmt=mk, ms=3.2, mfc=fc, mec=INK, ecolor=INK, elinewidth=0.75, label=lab if i == 0 else None)
            rows.append({"tissue": coh, "metric": key, "estimate": b[key], "bc_lo": b[f"{key}_bc_lo"], "bc_hi": b[f"{key}_bc_hi"]})
    PLN = NU["planted"]["0.2"]
    for i, coh in enumerate(["fibroblasts", "skin"]):
        ax.scatter(PLN[coh]["median_top_rho"], i - 0.22, marker="v", s=16, color=INK, zorder=4, label="planted effect" if i == 0 else None)
    ax.axvline(0, **REF); ax.set_yticks([0, 1]); ax.set_yticklabels(["GTEx\nfibro.", "GTEx\nskin"]); ax.invert_yaxis()
    ax.set_xlim(-0.9, 1.0); ax.set_xticks([-0.5, 0, 0.5, 1]); ax.set_xticklabels(["−0.5", "0", "0.5", "1"]); ax.set_xlabel("Spearman ρ (95% CI)")
    hh, ll = ax.get_legend_handles_labels(); hh.append(Rectangle((0, 0), 1, 1, fc="#e6e6e6")); ll.append("permuted age")
    ax.legend(hh, ll, loc="center left", bbox_to_anchor=(1.0, 0.5), ncol=1, fontsize=6.5, handletextpad=0.3, labelspacing=0.3, frameon=False); gridx(ax)
    title(ax, "200 strongest\ndiscovery genes")
    src(pd.DataFrame(rows), "Fig1D", index=False)
    # E and F: every discovery specification against G0-G2, genome-wide (E) and among the 200 strongest genes (F)
    ax = E = mm_ax(fig, 33, 99, 50, 33)
    specs, colsG, Mf = spec_heat(ax, "rho_all", p, "p_rho_all", None); ax.set_yticklabels([FSHORT[s[:2]] for s in specs])
    title(ax, "All genes"); src(pd.DataFrame(Mf, index=specs, columns=[c.replace(": ", " ") for c in colsG]), "Fig1E")
    ax = Fh = mm_ax(fig, 92, 99, 50, 33)
    specs, colsG, Mt = spec_heat(ax, "top_rho", p, "p_top_rho", None); ax.set_yticklabels([])
    sm = plt.cm.ScalarMappable(cmap=CORR_CMAP, norm=plt.Normalize(-CORR_LIM, CORR_LIM))
    cbar(fig, sm, ax, "Spearman ρ", ticks=[-0.6, 0, 0.6], frac=0.035, pad=0.04, ticklabels=["−0.6", "0", "0.6"], extend="both")
    title(ax, "200 strongest genes"); src(pd.DataFrame(Mt, index=specs, columns=[c.replace(": ", " ") for c in colsG]), "Fig1F")
    # G Hallmark enrichment, k/n inside the axes
    ax = Gx = mm_ax(fig, 34, 144, 70, 24)
    En = CSV("signature_hallmark_enrichment.csv")
    top_d = En[En.dir == "down"].nsmallest(4, "p").sort_values("p"); top_u = En[En.dir == "up"].nsmallest(4, "p").sort_values("p")
    ren = {"UV Response Dn": "UV response down", "heme Metabolism": "Heme metabolism", "TGF-beta Signaling": "TGF-β signaling", "G2-M Checkpoint": "G2M checkpoint", "E2F Targets": "E2F targets", "Mitotic Spindle": "Mitotic spindle", "p53 Pathway": "p53 pathway"}
    yy = list(range(4)) + list(range(5, 9)); Et = pd.concat([top_d, top_u]); vals = -np.log10(Et.p.values)
    ax.barh(yy, vals, color=[DOWN] * 4 + [UP] * 4, height=0.66)
    xk = vals.max() * 1.03
    for y_, (_, r_) in zip(yy, Et.iterrows()): ax.text(xk, y_, f"{int(r_.overlap)}/{int(r_['size'])}", va="center", fontsize=6.5, color=GRAY)
    ax.text(xk, -1.2, "k/n", va="center", fontsize=6.5, color=GRAY)
    fdr_p = Et[Et.q < 0.05].p.max(); ax.axvline(-np.log10(fdr_p), **REF)
    ax.set_yticks(yy); ax.set_yticklabels([ren.get(s, s) for s in Et.set]); ax.invert_yaxis(); ax.set_xlim(0, xk * 1.18)
    ax.set_xlabel("−log$_{10}$ P (one-sided hypergeometric)")
    ax.legend([Rectangle((0, 0), 1, 1, fc=DOWN), Rectangle((0, 0), 1, 1, fc=UP)], ["down with age", "up with age"], loc="lower center", bbox_to_anchor=(0.62, 0.02), fontsize=6.5)
    title(ax, "Hallmark sets in the age signature")
    src(Et, "Fig1G", index=False)
    # H age prediction, Fisher CIs; key at the right
    ax = H = mm_ax(fig, 127, 144, 22, 24)
    AP = CSV("rev3_cohort_age_prediction.csv"); trn = ["Discovery", "fibroblasts", "skin"]; nmap = {"Discovery": 111, "fibroblasts": 652, "skin": 651}
    tc = {"Discovery": COHORT_C["Discovery"], "fibroblasts": COHORT_C["GTEx fibroblasts"], "skin": COHORT_C["GTEx skin"]}
    for i, t_ in enumerate(trn):
        for k, (_, r_) in enumerate(AP[AP.train == t_].iterrows()):
            within = "CV" in r_.test; te = t_ if within else r_.test; lo, hi = fisher_ci(r_.pearson_r, r_.n_test)
            x = i + (-0.25 if within else (0.05 if k == 1 else 0.25))
            ax.errorbar(x, r_.pearson_r, yerr=EB([r_.pearson_r], [lo], [hi]), fmt="o", ms=3.0, mfc=tc[te] if within else "white", mec=tc[te], ecolor=tc[te], elinewidth=0.75)
    ax.axhline(0, **REF); ax.set_xticks(range(3)); ax.set_xticklabels(["Disc.", "GTEx\nfibro.", "GTEx\nskin"])
    ax.set_ylim(-0.3, 1.0); ax.set_ylabel("Pearson r (95% CI)"); ax.set_xlabel("Training cohort")
    ax.legend([Line2D([], [], marker="o", ls="none", mfc=INK, mec=INK), Line2D([], [], marker="o", ls="none", mfc="white", mec=INK)] +
              [Line2D([], [], marker="s", ls="none", color=tc[k_], ms=3) for k_ in trn],
              ["within (CV)", "other", "Disc.", "GTEx fibro.", "GTEx skin"], loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=6.5, ncol=1,
              handletextpad=0.2, labelspacing=0.3, frameon=False, title="test cohort", title_fontsize=6.5); gridy(ax)
    title(ax, "Age predictor transfer")
    src(AP, "Fig1H", index=False)
    letters(fig, [[(A, "A", 0.003)], [(B, "B", 0.003), (Cm, "C"), (Dp, "D")], [(E, "E", 0.003), (Fh, "F")], [(Gx, "G", 0.003), (H, "H")]])
    save(fig, "Fig1")


# ================================================================== FIGURE 2
def fig2():
    Hf = 180
    fig = plt.figure(figsize=(W2, Hf * MM))
    A = fig.add_axes([0, 1 - 41 / Hf, 1, 40 / Hf]); schematics.fig2a(A)
    gs = fig.add_gridspec(3, 1, height_ratios=[44, 52, 40], hspace=0.6, left=0.145, right=0.985, top=1 - 58 / Hf, bottom=0.055)
    gb = gs[0].subgridspec(1, 5, wspace=0.12); axes_b = []
    for i, c in enumerate(SC):
        ax = fig.add_subplot(gb[i]); axes_b.append(ax); s = scores(c)
        ax.scatter(s.prolif, s.R_raw, s=1.0, c=SCREEN_C[c], lw=0, alpha=0.15 if c == "K562gw" else 0.35, rasterized=True)
        q = pd.qcut(s.prolif.rank(method="first"), 10, labels=False); g = s.groupby(q.values).agg(x=("prolif", "median"), y=("R_raw", "median"))
        ax.plot(g.x, g.y, color=INK, lw=1.0, marker="o", ms=1.8)
        ax.axhline(0, color=LIGHT, lw=0.4, zorder=0); ax.axvline(0, color=LIGHT, lw=0.4, zorder=0)
        ax.set_ylim(-0.70, 0.45); ax.set_xlim(-0.5, 0.15); ax.set_xticks([-0.4, -0.2, 0]); ax.set_xticklabels(["−0.4", "−0.2", "0"])
        ax.set_title(f"ρ = {fmt(H2[c]['rho_raw_prolif'])}, n = {len(s):,}", fontsize=6.5)
        tn = ax.text(0.5, 1.2, SCREEN_L[c], transform=ax.transAxes, ha="center", va="center", fontsize=6.5)
        fig.canvas.draw(); bb = tn.get_window_extent().transformed(ax.transAxes.inverted())
        ax.scatter([bb.x0 - 0.07], [1.2], transform=ax.transAxes, s=12, color=SCREEN_C[c], clip_on=False)
        if i == 0: ax.set_ylabel("Raw score R")
        else: ax.set_yticklabels([])
        if i == 2: ax.set_xlabel("Proliferation shift (mean z)")
        src(s[["prolif", "R_raw"]], f"Fig2B_{c}")
    axes_b[2].text(0.5, 1.36, "Raw score versus proliferation", transform=axes_b[2].transAxes, ha="center", va="bottom", fontsize=7)
    gc = gs[1].subgridspec(1, 3, width_ratios=[1.45, 0.95, 0.16], wspace=0.55); ge = gs[2].subgridspec(1, 2, wspace=0.35)
    # C correlation with proliferation by definition
    ax = Cx = fig.add_subplot(gc[0])
    dep = CSV("rev3_scores_prolif_dependence.csv"); M = dep.pivot(index="variant", columns="screen", values="spearman").loc[VARS, SC]
    im = ax.imshow(M.values, cmap=CORR_CMAP, vmin=-CORR_LIM, vmax=CORR_LIM, aspect="auto")
    for (i, j), v in np.ndenumerate(M.values): heat_text(ax, j, i, v, CORR_CMAP, -CORR_LIM, CORR_LIM)
    nsc = {c: len(scores(c)) for c in SC}
    ax.set_xticks(range(5)); ax.set_xticklabels([f"{SCREEN_L[c]}\nn = {nsc[c]:,}" for c in SC], fontsize=6.5)
    ax.set_yticks(range(len(VARS))); ax.set_yticklabels([VLABEL[v] for v in VARS]); ax.get_yticklabels()[1].set_fontweight("bold")
    ax.add_patch(Rectangle((-0.5, 0.5), 5, 1, fill=False, ec="black", lw=0.9)); ax.spines[:].set_visible(False)
    cbar(fig, im, ax, "Spearman ρ with\nproliferation shift", ticks=[-0.6, 0, 0.6], frac=0.045, ticklabels=["−0.6", "0", "0.6"], extend="both")
    title(ax, "Score versus proliferation")
    src(M, "Fig2C")
    # D top-5% overlap (V1-V5 rows, V0-V4 columns: lower triangle without the empty row and column)
    ax = D = fig.add_subplot(gc[1]); RS = J("rev_scores.json")
    V6 = VARS[:6]; O = pd.DataFrame(RS["variant_top5_matrix"]).loc[[OLDV[v] for v in V6], [OLDV[v] for v in V6]].values
    sub = O[1:, :5].copy(); sub[np.triu_indices(5, 1)] = np.nan
    im = ax.imshow(sub, cmap=COUNT_CMAP, vmin=0, vmax=1, aspect="auto")
    for (i, j), v in np.ndenumerate(sub):
        if np.isfinite(v): heat_text(ax, j, i, v, COUNT_CMAP, 0, 1, txt=f"{v:.2f}")
    ax.set_xticks(range(5)); ax.set_xticklabels([VSHORT[v] for v in V6[:5]]); ax.set_yticks(range(5)); ax.set_yticklabels([VSHORT[v] for v in V6[1:]])
    ax.spines[:].set_visible(False)
    cbar(fig, im, ax, "Fraction of\ntop 5% shared", ticks=[0, 0.5, 1], frac=0.06)
    ax.set_xlabel(f"Definition (top 5% = {RS['k_top5']} of 2,393\ngenes in ≥ 2 lines)")
    title(ax, "Shared top 5%", "exp")
    src(pd.DataFrame(sub, index=V6[1:], columns=V6[:5]), "Fig2D")
    # E random-signature null against the real signature without proliferation genes (V2)
    ax = E = fig.add_subplot(ge[0]); N = CSV("rev3_scores_random_signatures.csv")
    real_v2 = {c: dep[(dep.screen == c) & (dep.variant == "V2 genes removed")].spearman.iloc[0] for c in SC}
    for i, c in enumerate(SC):
        v = N[f"H2_rho_{c}"].values
        ax.scatter(np.full(len(v), i) + np.random.default_rng(i).uniform(-0.2, 0.2, len(v)), v, s=1.4, color=LIGHT, lw=0, rasterized=True)
        ax.plot([i - 0.3, i + 0.3], [real_v2[c]] * 2, color=ACCENT, lw=2.2, solid_capstyle="butt")
    ax.set_xticks(range(5)); ax.set_xticklabels([SCREEN_L[c] for c in SC], fontsize=6.5); ax.set_ylim(-0.4, 1)
    ax.set_ylabel("ρ(score, proliferation)"); gridy(ax)
    ax.legend([Line2D([], [], color=ACCENT, lw=2.2), Line2D([], [], marker="o", ls="none", color=LIGHT, ms=2)], ["real signature (V2 input)", f"{len(N)} random signatures"],
              loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=6.5, handlelength=1.0, frameon=False, columnspacing=0.8)
    ax.set_title("Random signatures", pad=15)
    src(N.assign(**{f"real_V2_{c}": real_v2[c] for c in SC}), "Fig2E", index=False)
    # F retention of raw top 5% under V1
    ax = Fr = fig.add_subplot(ge[1])
    from statsmodels.stats.proportion import proportion_confint
    vals, los, his, ks = [], [], [], []
    for c in SC:
        k = H2[c]["k"]; v = H2[c]["top5_overlap_raw_adj"]; lo, hi = proportion_confint(int(round(v * k)), k, method="wilson")
        vals.append(100 * v); los.append(100 * lo); his.append(100 * hi); ks.append(k)
    for i, c in enumerate(SC):
        ax.errorbar(i, vals[i], yerr=EB([vals[i]], [los[i]], [his[i]]), fmt="o", ms=3.6, color=SCREEN_C[c], elinewidth=0.75)
    ax.axhline(5, **REF); ax.text(4.4, 6.5, "chance", fontsize=6.5, ha="right", color=GRAY)
    ax.set_xticks(range(5)); ax.set_xticklabels([f"{SCREEN_L[c]}\nk = {k}" for c, k in zip(SC, ks)], fontsize=6.5); ax.set_xlim(-0.5, 4.5)
    ax.set_ylim(0, 60); ax.set_ylabel("Raw top 5% kept by V1 (%)"); gridy(ax)
    title(ax, "Kept after adjustment", "exp")
    src(pd.DataFrame({"screen": SC, "k": ks, "kept_pct": vals, "lo": los, "hi": his}), "Fig2F", index=False)
    letters(fig, [[(A, "A", 0.003)], [(axes_b[0], "B", 0.003)], [(Cx, "C", 0.003), (D, "D")], [(E, "E", 0.003), (Fr, "F")]])
    save(fig, "Fig2")


# ================================================================== FIGURE 3
def fig3():
    fig = plt.figure(figsize=(W2, 132 * MM))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.3, 1, 1], hspace=0.45, wspace=0.7, left=0.17, right=0.985, top=0.9, bottom=0.09)
    S = {c: scores(c) for c in SC}
    dis = CSV("rev3_reliability_disattenuated.csv"); dv1 = dis[dis.variant == "V1 pre-specified"].set_index("pair").rho_disattenuated
    # A dumbbell
    ax = A = fig.add_subplot(gs[0, 0]); rows = []
    for i_, a in enumerate(SC):
        for b in SC[i_ + 1:]:
            sh = S[a].index.intersection(S[b].index)
            rows.append(dict(a=a, b=b, n=len(sh), raw=stats.spearmanr(S[a].R_raw[sh], S[b].R_raw[sh])[0], adj=stats.spearmanr(S[a].R_adj[sh], S[b].R_adj[sh])[0],
                             dis=dv1.get(f"{a}-{b}", np.nan)))
    Dd = pd.DataFrame(rows); same = (Dd.a == "K562") & (Dd.b == "K562gw")
    Dd = pd.concat([Dd[same], Dd[~same].sort_values("raw", ascending=False)]).reset_index(drop=True); y = np.arange(len(Dd))
    ax.axhspan(-0.5, 0.5, color="#ececec", zorder=0)
    for yi, r in zip(y, Dd.itertuples()): ax.plot([r.adj, r.raw], [yi, yi], color=LIGHT, lw=1.4, zorder=1, solid_capstyle="round")
    ax.scatter(Dd.raw, y, s=16, facecolor="white", edgecolor=INK, lw=0.8, zorder=3, label="raw R")
    ax.scatter(Dd.adj, y, s=16, color=ACCENT, zorder=3, label=f"{RADJ} (V1)")
    ax.scatter(Dd.dis, y, marker="D", s=11, facecolor="white", edgecolor=ACCENT, lw=0.8, zorder=4, label=f"{RADJ}, reliability-corrected")
    ax.set_yticks(y); ax.set_yticklabels([f"{SCREEN_L[r.a]} – {SCREEN_L[r.b]}" for r in Dd.itertuples()]); ax.invert_yaxis()
    yt = ax.get_yaxis_transform()
    for yi, r in zip(y, Dd.itertuples()): ax.text(1.03, yi, f"{r.n:,}", transform=yt, fontsize=6.5, va="center", ha="left", color=GRAY)
    ax.text(1.03, -0.9, "n", transform=yt, fontsize=6.5, va="center", ha="left", color=GRAY)
    ax.set_xlim(0, 1); ax.set_xlabel("Spearman ρ between screens"); gridx(ax)
    ax.legend(loc="lower center", bbox_to_anchor=(0.45, 1.0), ncol=2, fontsize=6.5, handletextpad=0.1, columnspacing=0.5, frameon=False)
    ax.set_title("Agreement between screens", pad=26)
    src(Dd, "Fig3A", index=False)
    # B dot plot per definition
    ax = B = fig.add_subplot(gs[0, 1:]); PV = J("rev_scores.json")["per_variant"]; V6 = VARS[:6]; rowsb = []
    dm = dis.groupby("variant").rho_disattenuated.median()
    for xi, v in enumerate(V6):
        pv = list(PV[OLDV[v]]["cross_line_rho"].values()); med = PV[OLDV[v]]["cross_line_median"]; wk = PV[OLDV[v]]["within_K562_rho"]
        ax.scatter(np.full(len(pv), xi - 0.12) + np.linspace(-0.06, 0.06, len(pv)), pv, s=7, color="#cbc9e2", zorder=3)
        ax.plot([xi - 0.28, xi + 0.04], [med] * 2, color=ACCENT, lw=1.6, zorder=4)
        ax.scatter(xi + 0.2, wk, marker="s", s=18, facecolor="white", edgecolor=SCREEN_C["K562"], lw=0.9, zorder=5)
        if np.isfinite(dm.get(v, np.nan)): ax.plot([xi - 0.28, xi + 0.04], [dm[v]] * 2, color=ACCENT, lw=0.9, ls=":", zorder=4)
        rowsb.append({"definition": v, "pairs": pv, "median": med, "within_K562": wk, "disattenuated_median": dm.get(v, np.nan)})
    ax.set_xticks(range(6)); ax.set_xticklabels([VSHORT[v] for v in V6]); ax.get_xticklabels()[1].set_fontweight("bold")
    ax.set_ylim(0, 0.8); ax.set_ylabel("Spearman ρ"); ax.set_xlabel("Score definition (n = 1,852–2,270 knockdowns per pair)"); gridy(ax)
    ax.legend([Line2D([], [], marker="o", ls="none", color="#cbc9e2", ms=2.5), Line2D([], [], color=ACCENT, lw=1.6), Line2D([], [], color=ACCENT, lw=0.9, ls=":"),
               Line2D([], [], marker="s", ls="none", mfc="white", mec=SCREEN_C["K562"], ms=3.5)],
              ["six line pairs", "median", "median, reliability-corrected", "within K562 (two screens)"], ncol=2, loc="upper right", fontsize=6.5)
    title(ax, "Line versus experiment agreement")
    src(pd.DataFrame(rowsb), "Fig3B", index=False)
    # C reliability heatmap
    ax = Cx = fig.add_subplot(gs[1, 0]); rel = CSV("rev3_reliability_summary.csv")
    Mr = rel.pivot(index="variant", columns="screen", values="reliability_SB").reindex(index=VARS, columns=["HepG2", "Jurkat", "K562"])
    im = ax.imshow(Mr.values, cmap=COUNT_CMAP, vmin=0, vmax=1, aspect="auto")
    for (i, j), v in np.ndenumerate(Mr.values): heat_text(ax, j, i, v, COUNT_CMAP, 0, 1, txt=f"{v:.2f}")
    ax.set_xticks(range(3)); ax.set_xticklabels(["HepG2\nsplit cells", "Jurkat\nsplit cells", "K562\ntest–retest"]); ax.set_yticks(range(len(VARS))); ax.set_yticklabels([VLABEL[v] for v in VARS])
    ax.spines[:].set_visible(False); cbar(fig, im, ax, "Reliability", ticks=[0, 0.5, 1], frac=0.06)
    title(ax, "Score reliability", "exp")
    src(Mr, "Fig3C")
    # D donor bootstraps: three rows sharing x
    ax = D = fig.add_subplot(gs[1, 1])
    b2000 = CSV("rev_h5_boot.csv").rho; b200 = CSV("h5_boot_rho.csv", index_col=0).iloc[:, 0]; bm = CSV("rev3_h5_multiline.csv").rho_multi
    sets = [("all genes,\n200 bootstraps", b200), ("all genes,\n2,000", b2000), ("genes in ≥2\nlines, 2,000", bm)]
    vp = ax.violinplot([s_.values for _, s_ in sets], positions=[0, 1, 2], orientation="horizontal", widths=0.8, showextrema=False)
    for k, b_ in enumerate(vp["bodies"]): b_.set_facecolor([ACCENT, "#cbc9e2", "#cbc9e2"][k]); b_.set_edgecolor("none"); b_.set_alpha(1)
    for k, (_, s_) in enumerate(sets): ax.plot([s_.median()] * 2, [k - 0.32, k + 0.32], color=INK, lw=1.2)
    ax.axvline(0.8, color=INK, ls="--", lw=0.6); ax.text(0.81, -0.45, "threshold 0.8", fontsize=6.5, ha="left", va="center")
    ax.set_yticks(range(3)); ax.set_yticklabels([s for s, _ in sets]); ax.invert_yaxis(); ax.set_xlim(0.45, 0.97)
    ax.set_xlabel(f"ρ(bootstrap {RROB}, full {RROB})")
    title(ax, "Donor bootstraps", "pre"); ax.set_ylim(2.6, -0.75)
    src(pd.DataFrame({"rho_all_2000": b2000.values}), "Fig3D", index=False)
    # E cohort swap per definition, Fisher CIs over genes
    ax = E = fig.add_subplot(gs[1, 2]); n = 9870; rowse = []
    for xi, v in enumerate(V6):
        for k, (key, coh) in enumerate([("swap_gtex_fibroblast_rho", "GTEx fibroblasts"), ("swap_gtex_skin_rho", "GTEx skin")]):
            r_ = PV[OLDV[v]][key]; lo, hi = fisher_ci(r_, n); x = xi + (k - 0.5) * 0.3
            pre = (v == "V1 pre-specified" and k == 0)
            ax.errorbar(x, r_, yerr=EB([r_], [lo], [hi]), fmt="o", ms=3.4, color=COHORT_C[coh], mfc=COHORT_C[coh] if pre else "white", mec=COHORT_C[coh], elinewidth=0.75)
            rowse.append({"definition": v, "cohort": coh, "rho": r_, "lo": lo, "hi": hi, "n": n})
    ax.axhline(0, **REF); ax.set_xticks(range(6)); ax.set_xticklabels([VSHORT[v] for v in V6]); ax.set_ylim(-0.25, 0.62)
    ax.set_ylabel(f"ρ with Discovery-signature {RROB}"); gridy(ax)
    ax.legend([Line2D([], [], marker="o", ls="none", mfc="white", mec=COHORT_C[c]) for c in ("GTEx fibroblasts", "GTEx skin")],
              ["GTEx fibroblasts", "GTEx skin"], loc="upper left", fontsize=6.5, labelspacing=0.2, borderpad=0.2)
    ax.set_xlabel("Definition (n = 9,870 genes)")
    title(ax, "Cohort swap", "pre")
    src(pd.DataFrame(rowse), "Fig3E", index=False)
    letters(fig, [[(A, "A", 0.003), (B, "B")], [(Cx, "C", 0.003), (D, "D"), (E, "E")]])
    save(fig, "Fig3")


# ================================================================== FIGURE 4
PRED = [("constant", "constant", "o", "white", "#999999"), ("gene_covariates", "gene covariates", "^", "white", GRAY), ("network", "STRING", "s", "white", INK),
        ("physical_network", "STRING physical", "s", "#bdbdbd", "#7f7f7f"), ("knn_other_lines", "nearest neighbors, other lines", "v", "white", INK),
        ("transfer_calibrated", "calibrated line transfer", "D", "#7f7f7f", "#7f7f7f"), ("transfer_plus_network", "line transfer + network", "D", ACCENT, ACCENT),
        ("knn_response", "nearest neighbors, own profile", "x", INK, INK)]


def fig4():
    fig = plt.figure(figsize=(W2, 142 * MM))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1], hspace=0.85, wspace=0.42, left=0.08, right=0.985, top=0.88, bottom=0.1)
    NET = CSV("rev_network.csv"); ND = CSV("rev_network_null_dist.csv"); NJ = J("rev_network_null.json"); h3c = R["H3"]["contexts"]; x = np.arange(5)
    n3a = {c: int(h3c[c]["h3a_n"]) for c in SC}
    xl = [SCREEN_L[c].replace(" ", "\n") + f"\n{n3a[c]:,}" for c in SC]
    ax = A = fig.add_subplot(gs[0, 0])
    d10 = [100 * (h3c[c]["h3a_mae_network"] - h3c[c]["h3a_mae_null"]) / h3c[c]["h3a_mae_null"] for c in SC]
    d500 = [100 * (NJ["screens"][c]["mae_network"] - NJ["screens"][c]["mae_mean_null"]) / NJ["screens"][c]["mae_mean_null"] for c in SC]
    ax.bar(x - 0.2, d10, 0.38, color=ACCENT, label="10 nulls"); ax.bar(x + 0.2, d500, 0.38, color="white", ec=ACCENT, hatch="////", label="500 nulls")
    ax.axhline(0, color=INK, lw=0.5); ax.set_xticks(x); ax.set_xticklabels(xl, fontsize=6.5); ax.set_ylim(-6.5, 9)
    ax.set_ylabel("Network vs null mean,\nchange in MAE (%)"); ax.legend(fontsize=6.5, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False, columnspacing=0.6, handletextpad=0.3); gridy(ax)
    ax.set_title("Absolute error", pad=15)
    src(pd.DataFrame({"screen": SC, "n": [n3a[c] for c in SC], "dmae_10": d10, "dmae_500": d500}), "Fig4A", index=False)
    ax = B = fig.add_subplot(gs[0, 1])
    vp = ax.violinplot([ND[ND.screen == c].rho.values for c in SC], positions=x, widths=0.7, showextrema=False)
    for b_ in vp["bodies"]: b_.set_facecolor(LIGHT); b_.set_alpha(1); b_.set_edgecolor("none")
    dn = NET[(NET.task.str.startswith("H3a")) & (NET.predictor == "network")].set_index("screen").reindex(SC)
    ax.errorbar(x, dn.rho, yerr=EB(dn.rho, dn.rho_lo, dn.rho_hi), fmt="D", ms=3.2, color=ACCENT, ecolor=ACCENT, elinewidth=0.75)
    ax.set_xticks(x); ax.set_xticklabels(xl, fontsize=6.5); ax.set_ylim(-0.12, 0.45); ax.axhline(0, color=INK, lw=0.5)
    ax.set_ylabel("Spearman ρ, predicted vs observed"); gridy(ax)
    ax.legend([Line2D([], [], marker="D", ls="none", color=ACCENT), Rectangle((0, 0), 1, 1, fc=LIGHT)], ["STRING network", "500 randomized"], loc="lower center",
              bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=6.5, frameon=False, columnspacing=0.6, handletextpad=0.3)
    ax.set_title("Rank agreement", pad=15)
    src(ND, "Fig4B", index=False)
    NS = CSV("rev3_network_summary.csv"); BB = CSV("rev3_network_block_bootstrap.csv").set_index("screen").reindex(SC)
    nn = {c: int(BB.n[c]) for c in SC}; xl2 = [SCREEN_L[c].replace(" ", "\n") + f"\n{nn[c]:,}" for c in SC]
    ax = C = fig.add_subplot(gs[0, 2])
    for k, (p, lab, col) in enumerate([("transfer_calibrated", "calibrated line transfer", "#9e9e9e"), ("transfer_plus_network", "line transfer + network", ACCENT)]):
        d = NS[NS.predictor == p].set_index("screen").reindex(SC)
        ax.bar(x + (k - 0.5) * 0.38, d.rho, 0.38, color=col, label=lab, yerr=d.rho_sd, error_kw=dict(elinewidth=0.6, ecolor=GRAY))
    ax.set_xticks(x); ax.set_xticklabels(xl2, fontsize=6.5); ax.set_ylim(0, 0.5); ax.set_ylabel("Spearman ρ, held-out"); ax.legend(fontsize=6.5, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=1, frameon=False, handletextpad=0.3, labelspacing=0.2); gridy(ax)
    ax.set_title("Adding the network", pad=24)
    src(NS, "Fig4C", index=False)
    ax = Dx = fig.add_subplot(gs[1, 0])
    base = BB.delta_mae / (BB.rel_change_pct / 100); lo, hi = 100 * BB.block_lo / base, 100 * BB.block_hi / base
    for xi, c in zip(x, SC):
        sig = BB.block_hi[c] < 0
        ax.bar(xi, BB.rel_change_pct[c], 0.6, color=ACCENT if sig else "white", ec=ACCENT, yerr=EB([BB.rel_change_pct[c]], [lo[c]], [hi[c]]), error_kw=dict(elinewidth=0.6, ecolor=GRAY))
    ax.axhline(0, color=INK, lw=0.5); ax.set_xticks(x); ax.set_xticklabels(xl2, fontsize=6.5)
    ax.set_ylabel("Change in MAE vs calibrated\nline transfer (%)"); gridy(ax)
    ax.legend([Rectangle((0, 0), 1, 1, fc=ACCENT), Rectangle((0, 0), 1, 1, fc="white", ec=ACCENT)], ["CI excludes 0", "CI includes 0"], loc="lower center",
              bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=6.5, frameon=False, columnspacing=0.6, handletextpad=0.3)
    ax.set_title("Network increment", pad=15)
    src(BB.reset_index(), "Fig4D", index=False)
    ax = E = fig.add_subplot(gs[1, 1:])
    for k, (p, lab, mk, fc, ec) in enumerate(PRED):
        if p == "constant" or p not in set(NS.predictor): continue
        d = NS[NS.predictor == p].set_index("screen").reindex(SC)
        ax.scatter(d.r2, x + (k - 4) * 0.13, marker=mk, s=16, facecolor=fc, edgecolor=ec, lw=0.8, label=lab, zorder=3)
    ceil = NS[NS.predictor == "constant"].set_index("screen").reindex(SC).get("reliability_ceiling_r2")
    if ceil is not None:
        for xi, v in zip(x, ceil):
            if np.isfinite(v): ax.plot([v, v], [xi - 0.35, xi + 0.35], color=INK, lw=1.8, zorder=2)
    ax.axvline(0, **REF); ax.set_yticks(x); ax.set_yticklabels(xl2, fontsize=6.5); ax.invert_yaxis(); ax.set_xlabel(f"Held-out R² of {RADJ}")
    hs, ls = ax.get_legend_handles_labels(); hs.append(Line2D([], [], color=INK, lw=1.8)); ls.append("reliability ceiling")
    ax.legend(hs, ls, ncol=3, loc="lower center", bbox_to_anchor=(0.45, 1.0), fontsize=6.5, columnspacing=0.8, handletextpad=0.5, handlelength=1.6, frameon=False)
    gridx(ax); ax.set_title("Held-out R²", pad=38)
    src(NS, "Fig4E", index=False)
    letters(fig, [[(A, "A"), (B, "B"), (C, "C")], [(Dx, "D"), (E, "E")]])
    save(fig, "Fig4")


# ================================================================== FIGURE 5
def hanley_se(A_, n1, n2):
    q1, q2 = A_ / (2 - A_), 2 * A_ * A_ / (1 + A_)
    return np.sqrt((A_ * (1 - A_) + (n1 - 1) * (q1 - A_ * A_) + (n2 - 1) * (q2 - A_ * A_)) / (n1 * n2))


def fig5():
    fig = plt.figure(figsize=(W2, 140 * MM))
    gs = fig.add_gridspec(3, 2, width_ratios=[1.12, 1], height_ratios=[1, 0.85, 0.8], hspace=0.95, wspace=1.42, left=0.15, right=0.9, top=0.9, bottom=0.1)
    EVr, PA = CSV("rev3_scores_evidence.csv"), CSV("rev3_scores_paired_auroc.csv")
    ev = lambda lab, gsn, sc_: EVr[(EVr.labels == lab) & (EVr.gene_set == gsn) & (EVr.score == sc_)].iloc[0]
    SCOL = {"R_rob (V1, pre-specified)": ACCENT, "mean raw score": INK, "proliferation shift alone": "#4d4d4d"}
    ax = A = fig.add_subplot(gs[:, 0])
    groups = [("CellAge", "CellAge: inducers vs inhibitors"), ("GenAge", "GenAge: anti- vs pro-longevity"),
              ("GenAge without yeast", "GenAge without yeast"), ("GenAge mammal-only", "GenAge, mouse genes only")]
    scs = [("R_rob (V1, pre-specified)", f"{RROB} (V1)", "o"), ("mean raw score", "raw score", "s"), ("proliferation shift alone", "proliferation only", "D")]
    y, ticks, rows, yt = 0.0, [], [], ax.get_yaxis_transform()
    for gi, (lab, head) in enumerate(groups):
        r1, r2 = ev(lab, "all genes", scs[0][0]), ev(lab, "genes in >= 2 lines", scs[0][0])
        ax.text(0.5, y, head, transform=yt, fontsize=6.5, fontweight="bold", ha="center", va="center")
        ax.text(1.03, y, f"{int(r1.n_pos)}/{int(r1.n_neg)} ({r1.min_detectable_auroc_80pct:.2f})", transform=yt, fontsize=6.5, va="center", color=INK)
        ax.text(1.03, y - 0.36, f"{int(r2.n_pos)}/{int(r2.n_neg)} ({r2.min_detectable_auroc_80pct:.2f})", transform=yt, fontsize=6.5, va="center", color=GRAY)
        y -= 0.62; ytop = y + 0.3
        for sc_, short, mk in scs:
            for k, (gsn, off) in enumerate([("all genes", 0.13), ("genes in >= 2 lines", -0.13)]):
                r = ev(lab, gsn, sc_); A_, n1, n2 = r.auroc, int(r.n_pos), int(r.n_neg); se = hanley_se(A_, n1, n2)
                pre = sc_.startswith("R_rob") and gsn == "all genes" and lab in ("CellAge", "GenAge")
                col = SCOL[sc_] if k == 0 else "#a6a6a6"
                ax.errorbar(A_, y + off, xerr=1.96 * se, fmt=mk, ms=3.2, mfc=col if pre else "white", mec=col, ecolor=col, elinewidth=0.75, zorder=3)
                rows.append(dict(labels=lab, score=sc_, gene_set=gsn, auroc=A_, lo=A_ - 1.96 * se, hi=A_ + 1.96 * se, n_pos=n1, n_neg=n2, mda=r.min_detectable_auroc_80pct))
            ticks.append((y, short)); y -= 0.75
        ax.plot([0.5, 0.5], [ytop, y + 0.45], **REF)                          # chance line per block, stopping below its header
        y -= 0.1
        if gi < len(groups) - 1: ax.axhline(y + 0.3, color=LIGHT, lw=0.5)
        y -= 0.3
    ax.set_yticks([t for t, _ in ticks]); ax.set_yticklabels([l for _, l in ticks])
    for tl, (_, short) in zip(ax.get_yticklabels(), ticks): tl.set_color({f"{RROB} (V1)": ACCENT, "raw score": INK, "proliferation only": "#4d4d4d"}[short])
    ax.set_xlim(0.1, 0.9); ax.set_ylim(y + 0.45, 0.4); ax.set_xlabel("AUROC (95% CI)")
    ax.text(1.03, 1.0, "n (MDA)", transform=ax.transAxes, fontsize=6.5, va="bottom", fontweight="bold")
    ax.legend([Line2D([], [], marker="o", ls="none", mfc="white", mec=INK), Line2D([], [], marker="o", ls="none", mfc="white", mec="#a6a6a6")], ["all genes", "genes in ≥2 lines"],
              loc="upper center", bbox_to_anchor=(0.5, -0.06), ncol=2, fontsize=6.5, frameon=False, columnspacing=0.6, handletextpad=0.2)
    title(ax, "Held-out labels")
    src(pd.DataFrame(rows), "Fig5A", index=False)
    ax = Bx = fig.add_subplot(gs[0, 1]); rowsb = []
    for k, (lab, mk) in enumerate([("CellAge", "o"), ("GenAge", "s")]):
        for xi, v in enumerate(VARS):
            nm = "R_rob (V1, pre-specified)" if v == "V1 pre-specified" else f"R_rob ({v})"
            r = ev(lab, "all genes", nm); se = hanley_se(r.auroc, int(r.n_pos), int(r.n_neg)); col = INK if v == "V0 raw" else ACCENT
            ax.errorbar(xi + (k - 0.5) * 0.3, r.auroc, yerr=1.96 * se, fmt=mk, ms=3, color=col, mfc=col if lab == "CellAge" else "white", elinewidth=0.75)
            rowsb.append(dict(labels=lab, definition=v, auroc=r.auroc, se=se, n_pos=int(r.n_pos), n_neg=int(r.n_neg)))
    ax.axhline(0.5, **REF); ax.set_xticks(range(len(VARS))); ax.set_xticklabels([VSHORT[v] for v in VARS]); ax.get_xticklabels()[1].set_fontweight("bold")
    ax.set_ylabel("AUROC (95% CI)"); ax.set_ylim(0.36, 0.66); gridy(ax)
    ax.legend([Line2D([], [], marker="o", ls="none", color=ACCENT), Line2D([], [], marker="s", ls="none", mfc="white", mec=ACCENT)], ["CellAge", "GenAge"],
              loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=6.5, frameon=False, columnspacing=0.8)
    ax.set_xlabel("Definition (n = 257/394; 573/376)")
    title(ax, " "); ax.text(0.5, 1.2, "By definition", transform=ax.transAxes, ha="center", va="bottom", fontsize=7)
    src(pd.DataFrame(rowsb), "Fig5B", index=False)
    ax = Cx = fig.add_subplot(gs[1, 1])
    comps = [("CellAge", "mean raw score", "proliferation shift alone", "raw − proliferation"), ("CellAge", "R_rob (V0 raw)", "R_rob (V1 pre-specified)", "V0 − V1"),
             ("CellAge", "R_rob (V0 raw)", "R_rob (V6 spline resid.)", "V0 − V6"), ("GenAge", "mean raw score", "proliferation shift alone", "raw − proliferation"),
             ("GenAge", "R_rob (V0 raw)", "R_rob (V1 pre-specified)", "V0 − V1"), ("GenAge", "R_rob (V0 raw)", "R_rob (V6 spline resid.)", "V0 − V6")]
    rowsc = []
    for yi, (lab, a, b, short) in enumerate(comps):
        r = PA[(PA.labels == lab) & (PA.gene_set == "all genes") & (PA.a == a) & (PA.b == b)].iloc[0]
        mk, fc = ("o", ACCENT) if lab == "CellAge" else ("s", "white")             # same code as (B)
        ax.errorbar(r.diff_a_minus_b, yi, xerr=[[r.diff_a_minus_b - r.lo], [r.hi - r.diff_a_minus_b]], fmt=mk, ms=3, color=ACCENT, mfc=fc, mec=ACCENT, elinewidth=0.75)
        ptxt = "< 0.001" if r.p_two_sided < 0.001 else f"{float(f'{r.p_two_sided:.2g}'):g}"
        ax.text(1.04, yi, ptxt, transform=ax.get_yaxis_transform(), fontsize=6.5, va="center", ha="left", color=INK)
        rowsc.append(dict(labels=lab, comparison=short, diff=r.diff_a_minus_b, lo=r.lo, hi=r.hi, p=r.p_two_sided, n_pos=int(r.n_pos), n_neg=int(r.n_neg)))
    ax.text(1.04, -0.75, "P", transform=ax.get_yaxis_transform(), fontsize=6.5, va="center", ha="left", fontweight="bold")
    ax.axvline(0, **REF); ax.set_yticks(range(len(comps))); ax.set_yticklabels([f"{c[0]}: {c[3]}" for c in comps]); ax.invert_yaxis()
    ax.set_xlabel("Paired AUROC difference (95% CI)"); gridx(ax)
    title(ax, "Paired tests")
    src(pd.DataFrame(rowsc), "Fig5C", index=False)
    ax = Dx = fig.add_subplot(gs[2, 1])
    OV, FI = CSV("rev3_scores_label_essential_overlap.csv"), CSV("rev3_depmap_label_fisher.csv")
    sel = [("CellAge", "Induces", "CellAge inducers"), ("CellAge", "Inhibits", "CellAge inhibitors"), ("GenAge", "Anti", "GenAge anti-longevity"), ("GenAge", "Pro", "GenAge pro-longevity")]
    for i_, (a, b, l) in enumerate(sel):
        r = OV[(OV.labels == a) & (OV["class"] == b)].iloc[0]
        ax.barh(i_, 100 * r.frac_DepMap_common_essential, color="#9ecae1" if i_ % 2 == 0 else "#2171b5", height=0.6)
        ax.plot([100 * r.frac_CEGv2] * 2, [i_ - 0.3, i_ + 0.3], color=INK, lw=1.4, ls=":")
        ax.text(100 * r.frac_DepMap_common_essential + 0.6, i_, f"{100 * r.frac_DepMap_common_essential:.0f}% of {int(r.n)}", va="center", fontsize=6.5)
    yt = ax.get_yaxis_transform()
    for k, lab in enumerate(["CellAge", "GenAge"]):
        f_ = FI[(FI["labels"] == lab) & (FI.reference == "DepMap_common_essential")].iloc[0]
        ptxt = "P < 0.001" if f_.fisher_p_two_sided < 0.001 else f"P = {f_.fisher_p_two_sided:.2f}"
        ax.plot([1.03, 1.05, 1.05, 1.03], [2 * k - 0.1, 2 * k - 0.1, 2 * k + 1.1, 2 * k + 1.1], transform=yt, color=GRAY, lw=0.6, clip_on=False)
        ax.text(1.08, 2 * k + 0.5, f"OR {f_.odds_ratio:.2f}\n{ptxt}", transform=yt, fontsize=6.5, va="center", ha="left")
    ax.set_yticks(range(4)); ax.set_yticklabels([s[2] for s in sel]); ax.invert_yaxis(); ax.set_xlim(0, 40)
    ax.set_xlabel("DepMap common-essential genes (%)")
    ax.legend([Rectangle((0, 0), 1, 1, fc="#9ecae1"), Rectangle((0, 0), 1, 1, fc="#2171b5"), Line2D([], [], color=INK, ls=":", lw=1.4)],
              ["positive class", "negative class", "Hart CEGv2"], loc="upper center", bbox_to_anchor=(0.45, -0.32), ncol=3, fontsize=6.5,
              handlelength=1.0, columnspacing=0.6, frameon=False)
    title(ax, "Essential genes")
    src(OV, "Fig5D", index=False)
    bx = 0.6
    letters(fig, [[(A, "A", 0.003), (Bx, "B", bx)], [(Cx, "C", bx)], [(Dx, "D", bx)]])
    save(fig, "Fig5")


if __name__ == "__main__":
    todo = sys.argv[1:] or ["fig1", "fig2", "fig3", "fig4", "fig5", "ga"]
    for name in todo:
        if name == "ga": schematics.graphical_abstract(FIG / "Graphical_abstract")
        else: globals()[name]()
        print("ok", name, flush=True)
