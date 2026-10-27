"""Shared figure standards for TRACE v0.1 (Cell Press): Arial, 6-8 pt at final size, 0.5-1 pt lines, RGB, vector PDF,
fixed Okabe-Ito palette with one meaning per colour, aligned 9 pt panel letters, pre-specified/exploratory tags."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Ellipse, Rectangle, Polygon, Circle

ROOT = r"D:\Study\RejuvBench"
FIG = os.path.join(ROOT, "figures_v01"); SRC = os.path.join(FIG, "source_data")
os.makedirs(SRC, exist_ok=True)
MM = 1 / 25.4
W1, W15, W2 = 85 * MM, 114 * MM, 174 * MM   # Cell Press column widths

plt.rcParams.update({
    "font.family": "Arial", "font.size": 6.5, "axes.titlesize": 7, "axes.labelsize": 6.5, "xtick.labelsize": 6,
    "ytick.labelsize": 6, "legend.fontsize": 6, "legend.title_fontsize": 6, "axes.linewidth": 0.5, "xtick.major.width": 0.5,
    "ytick.major.width": 0.5, "xtick.major.size": 2, "ytick.major.size": 2, "lines.linewidth": 0.8, "patch.linewidth": 0.5,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "pdf.fonttype": 42, "ps.fonttype": 42,
    "savefig.dpi": 600, "figure.dpi": 150, "axes.titlepad": 3, "axes.labelpad": 2, "mathtext.fontset": "custom",
    "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold", "axes.unicode_minus": True,
    "axes.edgecolor": "#222222", "xtick.color": "#222222", "ytick.color": "#222222", "text.color": "#111111",
    "axes.labelcolor": "#111111", "hatch.linewidth": 0.5})

INK, INK2, GRID = "#111111", "#4d4d4d", "#e3e3e3"
# screens (only ever used for screens)
SCREEN_C = {"K562": "#0072B2", "K562gw": "#56B4E9", "RPE1": "#555555", "HepG2": "#D55E00", "Jurkat": "#009E73"}
SCREEN_L = {"K562": "K562 ess.", "K562gw": "K562 gw", "RPE1": "RPE1", "HepG2": "HepG2", "Jurkat": "Jurkat"}
SCREENS = ["K562", "K562gw", "RPE1", "HepG2", "Jurkat"]
# cohorts (only ever used for cohorts); shades for covariate level G0 / G1 / G2
COHORT_C = {"Discovery": "#000000", "GTEx fibroblasts": "#CC79A7", "GTEx skin": "#E69F00"}
COHORT_SHADES = {"GTEx fibroblasts": ["#E7BCD5", "#CC79A7", "#7A3A62"], "GTEx skin": ["#F7D392", "#E69F00", "#8A5F00"]}
# predictors (Figure 4, S4): greys for baselines, one accent for anything containing the network
PRED_C = {"constant": "#D9D9D9", "null": "#BDBDBD", "transfer": "#8C8C8C", "transfer_calibrated": "#4D4D4D",
          "network": "#E69F00", "transfer_plus_network": "#B87D00", "knn_response": "#9ECAE1", "gene_covariates": "#6BAED6",
          "physical_network": "#F3C26B", "network_calibrated": "#E69F00"}
UP, DOWN = "#B2182B", "#2166AC"
VARS = ["V0 raw", "V1 pre-specified", "V2 genes removed", "V3 resid. (all 3)", "V4 resid. (prolif.)", "V5 stratified",
        "V6 spline resid.", "V7 resid. (prolif. + identity)"]
VSHORT = {v: v.split(" ", 1)[0] for v in VARS}
CORR_CMAP, OVERLAP_CMAP = "RdBu_r", "Greys"
FSHORT = {"F0": "F0 pre-specified (age + sex)", "F1": "F1 ages 20–79", "F2": "F2 + platform/repository/site", "F3": "F3 = F1 + F2",
          "F4": "F4 NextSeq only", "F5": "F5 arm biopsies only", "F6": "F6 expression filter 0.5", "F7": "F7 expression filter 2",
          "F8": "F8 without NIA donors ≥80", "F9": "F9 + 5 hidden factors"}


def fmt(x, d=2):
    if x is None or (isinstance(x, float) and not np.isfinite(x)): return "n/a"
    s = f"{x:.{d}f}"
    if float(s) == 0: s = s.lstrip("-")
    return s.replace("-", "\u2212")


def EB(v, lo, hi):
    v, lo, hi = np.asarray(v, float), np.asarray(lo, float), np.asarray(hi, float)
    return [np.clip(v - lo, 0, None), np.clip(hi - v, 0, None)]


def gridy(ax): ax.yaxis.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True)


def tag(ax, kind, text=None, loc=(1.0, 1.0)):
    """Small corner tag: kind = 'pre' (filled black) or 'exp' (open grey)."""
    t = text or ("Pre-specified" if kind == "pre" else "Exploratory")
    if kind == "pre":
        bbox = dict(boxstyle="round,pad=0.18,rounding_size=0.15", fc="#111111", ec="#111111", lw=0.4); col = "white"
    else:
        bbox = dict(boxstyle="round,pad=0.18,rounding_size=0.15", fc="white", ec="#8c8c8c", lw=0.5); col = "#4d4d4d"
    # sit on its own line above the panel title (title lines counted), so it never covers the title
    nlines = max(len(ax.get_title(loc="left").split("\n")) if ax.get_title(loc="left") else 0,
                 len(ax.get_title().split("\n")) if ax.get_title() else 0)
    ax.annotate(t, xy=(1.0, 1.0), xycoords="axes fraction", xytext=(0, 4 + 8.5 * nlines), textcoords="offset points",
                ha="right", va="bottom", fontsize=5.5, color=col, bbox=bbox, zorder=20, annotation_clip=False)


def letters(fig, axes_letters, dx=0.0, dy=0.004):
    """Place 9 pt bold letters at the top-left of each panel's tight bounding box; panels in one row share a baseline."""
    fig.canvas.draw(); r = fig.canvas.get_renderer(); inv = fig.transFigure.inverted()
    items = []
    for ax, L in axes_letters:
        bb = ax.get_tightbbox(r).transformed(inv); items.append((ax, L, bb))
    rows = []
    for it in sorted(items, key=lambda t: -t[2].y1):
        for row in rows:
            if abs(row[0][2].y1 - it[2].y1) < 0.035: row.append(it); break
        else: rows.append([it])
    for row in rows:
        y = max(t[2].y1 for t in row)
        for ax, L, bb in row:
            fig.text(bb.x0 + dx, y + dy, L, fontsize=9, fontweight="bold", ha="left", va="bottom")


def save(fig, name, width=None):
    """Vector PDF (fonts embedded) plus an RGB 600 dpi PNG for the Word file; no transparency."""
    from PIL import Image
    fig.savefig(os.path.join(FIG, name + ".pdf"), bbox_inches="tight", pad_inches=0.02, facecolor="white")
    p = os.path.join(FIG, name + ".png")
    fig.savefig(p, bbox_inches="tight", pad_inches=0.02, facecolor="white", dpi=600)
    Image.open(p).convert("RGB").save(p, dpi=(600, 600))
    plt.close(fig)


def cbar(fig, im, ax, label, ticks=None, frac=0.05, pad=0.03):
    cb = fig.colorbar(im, ax=ax, fraction=frac, pad=pad)
    cb.ax.tick_params(labelsize=5.5, width=0.4, length=1.5); cb.outline.set_linewidth(0.4); cb.set_label(label, fontsize=6)
    if ticks is not None: cb.set_ticks(ticks)
    return cb


def lowtri(M):
    M = np.array(M, float); M[np.triu_indices_from(M, 1)] = np.nan; return M


# ---- icons for schematics and the graphical abstract (axes-fraction coordinates, aspect-corrected)
def ratio(ax):
    p = ax.get_position(); W, H = ax.figure.get_size_inches(); return (p.width * W) / (p.height * H)


def rbox(ax, x, y, w, h, fc, ec, lw=0.8, r=0.02):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0.004,rounding_size={r}", fc=fc, ec=ec, lw=lw, zorder=1, clip_on=False))


def cell(ax, x, y, w, fc, ec="#7a7a7a", ar=0.55):
    h = w * ar * ratio(ax)
    ax.add_patch(Ellipse((x, y), w, h, fc=fc, ec=ec, lw=0.6, zorder=3, clip_on=False))
    ax.add_patch(Ellipse((x, y), w * 0.45, h * 0.45, fc="white", ec=ec, lw=0.5, zorder=3.1, clip_on=False))


def cluster(ax, cx, cy, fc, ec="#7a7a7a", s=0.02, n=5):
    r = ratio(ax)
    for dx, dy in [(-0.9, 0.35), (0.0, 0.55), (0.9, 0.3), (-0.45, -0.15), (0.45, -0.25), (1.35, -0.05), (-1.3, -0.3)][:n]:
        cell(ax, cx + dx * s * 1.05, cy + dy * s * 0.9 * r, s * 1.15, fc, ec)


def person(ax, x, y, s, col):
    r = ratio(ax)
    ax.add_patch(Ellipse((x, y + s * 0.95 * r), s * 0.62, s * 0.62 * r, fc=col, ec="none", zorder=4, clip_on=False))
    ax.add_patch(FancyBboxPatch((x - s * 0.42, y - s * 0.55 * r), s * 0.84, s * 1.05 * r, boxstyle=f"round,pad=0,rounding_size={s*0.3}",
                                mutation_aspect=r, fc=col, ec="none", zorder=4, clip_on=False))


def strip(ax, x, y, w, h, vals, cmap="RdBu_r"):
    cm = plt.get_cmap(cmap); n = len(vals); bw = w / n
    for i, v in enumerate(vals):
        ax.add_patch(Rectangle((x + i * bw, y), bw * 0.92, h, fc=cm(0.5 + 0.5 * np.clip(v, -1, 1)), ec="none", zorder=3, clip_on=False))


def arrow(ax, x0, y0, x1, y1, col=INK2, lw=0.9, ms=7):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), xycoords="axes fraction",
                arrowprops=dict(arrowstyle="-|>", color=col, lw=lw, mutation_scale=ms, shrinkA=0, shrinkB=0), zorder=2)


def T(ax, x, y, s, **k):
    k.setdefault("ha", "center"); k.setdefault("va", "center"); k.setdefault("zorder", 6)
    ax.text(x, y, s, transform=ax.transAxes, **k)


def canvas(ax):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")


# ---- colour-vision check (Machado et al. 2009 deuteranopia matrix, severity 1)
def deutan_png(src_png, dst_png):
    from PIL import Image
    M = np.array([[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413], [-0.011820, 0.042940, 0.968881]])
    a = np.asarray(Image.open(src_png).convert("RGB"), dtype=float) / 255
    lin = np.where(a <= 0.04045, a / 12.92, ((a + 0.055) / 1.055) ** 2.4)
    sim = np.clip(lin @ M.T, 0, 1)
    srgb = np.where(sim <= 0.0031308, sim * 12.92, 1.055 * sim ** (1 / 2.4) - 0.055)
    Image.fromarray((srgb * 255).astype(np.uint8)).save(dst_png)
