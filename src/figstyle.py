"""Figure standards for TRACE: built at final width (174 / 114 / 85 mm), Arial, text 6.5-7 pt, subscripts >= 6 pt,
9 pt panel letters on a fixed grid, vector PDF with embedded fonts plus RGB 600 dpi PNG.
Colors: screen colors only for screens, cohort colors only for cohorts, red/blue only for age direction; everything else in
black and grays, told apart by line style, marker shape and fill."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib._mathtext as _mt
from matplotlib.patches import Rectangle, FancyBboxPatch
from config import FIG, RES

_mt.SHRINK_FACTOR = 0.93            # subscripts at >= 6 pt for 6.5-7 pt labels
SRC = FIG / "source_data"; SRC.mkdir(parents=True, exist_ok=True)
MM = 1 / 25.4
W1, W15, W2 = 85 * MM, 114 * MM, 174 * MM

plt.rcParams.update({
    "font.family": "Arial", "font.size": 6.5, "axes.titlesize": 7, "axes.labelsize": 7, "xtick.labelsize": 6.5,
    "ytick.labelsize": 6.5, "legend.fontsize": 6.5, "legend.title_fontsize": 6.5, "axes.linewidth": 0.5, "xtick.major.width": 0.5,
    "ytick.major.width": 0.5, "xtick.major.size": 2, "ytick.major.size": 2, "lines.linewidth": 1.0, "patch.linewidth": 0.5,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": True, "legend.framealpha": 1.0, "legend.edgecolor": "none",
    "legend.borderpad": 0.25, "legend.handlelength": 1.2, "legend.handletextpad": 0.4, "legend.borderaxespad": 0.3,
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 600, "figure.dpi": 150, "axes.titlepad": 3, "axes.labelpad": 2,
    "axes.titlelocation": "center", "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold", "axes.unicode_minus": True, "axes.edgecolor": "#222222", "xtick.color": "#222222",
    "ytick.color": "#222222", "text.color": "#111111", "axes.labelcolor": "#111111", "hatch.linewidth": 0.5,
    "errorbar.capsize": 0, "lines.markersize": 3})

INK, GRAY, LIGHT, GRID = "#111111", "#666666", "#bdbdbd", "#eeeeee"
REF = dict(color="#666666", lw=0.5, ls="--")                       # reference lines (chance, zero, thresholds)
SCREEN_C = {"K562": "#0072B2", "K562gw": "#56B4E9", "RPE1": "#555555", "HepG2": "#D55E00", "Jurkat": "#009E73"}
SCREEN_L = {"K562": "K562 ess.", "K562gw": "K562 gw", "RPE1": "RPE1", "HepG2": "HepG2", "Jurkat": "Jurkat"}
SCREENS = ["K562", "K562gw", "RPE1", "HepG2", "Jurkat"]
COHORT_C = {"Discovery": "#000000", "GTEx fibroblasts": "#CC79A7", "GTEx skin": "#E69F00"}
UP, DOWN = "#B2182B", "#2166AC"
CORR_CMAP, COUNT_CMAP = "RdBu_r", "Blues"
ACCENT = "#6A51A3"                                              # one accent for "real" values against gray nulls
CORR_LIM = 0.6                                                  # one diverging scale for every correlation heatmap
VARS = ["V0 raw", "V1 pre-specified", "V2 genes removed", "V3 resid. (all 3)", "V4 resid. (prolif.)", "V5 stratified",
        "V6 spline resid.", "V7 resid. (prolif. + identity)"]
VSHORT = {v: v.split(" ", 1)[0] for v in VARS}
VLABEL = {"V0 raw": "V0 raw", "V1 pre-specified": "V1 adjusted", "V2 genes removed": "V2 genes removed", "V3 resid. (all 3)": "V3 resid. (π, ι, μ)",
          "V4 resid. (prolif.)": "V4 resid. (π)", "V5 stratified": "V5 stratified", "V6 spline resid.": "V6 spline (π)", "V7 resid. (prolif. + identity)": "V7 resid. (π, ι)"}
FSHORT = {"F0": "F0 main", "F1": "F1 ages 20–79", "F2": "F2 + platform, repository, site", "F3": "F3 = F1 + F2",
          "F4": "F4 NextSeq only", "F5": "F5 arm biopsies only", "F6": "F6 expression filter 0.5", "F7": "F7 expression filter 2",
          "F8": "F8 without AG donors ≥80", "F9": "F9 + 5 hidden factors"}
RADJ, RROB = "R$_\\mathregular{adj}$", "R$_\\mathregular{rob}$"


def fmt(x, d=2):
    if x is None or (isinstance(x, float) and not np.isfinite(x)): return "n/a"
    s = f"{x:.{d}f}"
    if float(s) == 0: s = s.lstrip("-")
    return s.replace("-", "−")


def text_on(rgba):
    """One rule for text on heatmap cells: white on cells darker than 50% luminance, black otherwise."""
    r, g, b = rgba[:3]; return "white" if 0.2126 * r + 0.7152 * g + 0.0722 * b < 0.5 else INK


def heat_text(ax, x, y, v, cmap, vmin, vmax, txt=None, size=6.5):
    c = plt.get_cmap(cmap)((np.clip(v, vmin, vmax) - vmin) / (vmax - vmin))
    ax.text(x, y, txt if txt is not None else fmt(v), ha="center", va="center", fontsize=size, color=text_on(c))


def EB(v, lo, hi):
    v, lo, hi = (np.asarray(a, float) for a in (v, lo, hi))
    return [np.clip(v - lo, 0, None), np.clip(hi - v, 0, None)]


def gridx(ax): ax.xaxis.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True)


def gridy(ax): ax.yaxis.grid(True, color=GRID, lw=0.4); ax.set_axisbelow(True)


def title(ax, text, kind=None, tag=None, pill=False):
    """Plain panel title centered over the plot. The status (pre-specified or exploratory) is stated in the legend, so the
    in-figure tag is off unless pill=True."""
    t = ax.set_title(text, loc="center")
    if kind and pill:
        lab = tag or ("Pre-specified" if kind == "pre" else "Exploratory")
        box = dict(boxstyle="round,pad=0.16,rounding_size=0.12", fc="#111111" if kind == "pre" else "white", ec="#111111" if kind == "pre" else "#8c8c8c", lw=0.5)
        ax.annotate(lab, xy=(1.0, 1.0), xycoords="axes fraction", xytext=(0, 3 + 0.5), textcoords="offset points", ha="right", va="bottom",
                    fontsize=6, color="white" if kind == "pre" else "#4d4d4d", bbox=box, annotation_clip=False)
    return t


def letters(fig, rows, top_mm=1.6):
    """rows: list of lists of (ax, letter). Letters sit at the left edge of each panel's tight box; in each row they share one
    baseline (the highest panel top in that row); panels in the same grid column share the x of the leftmost box in that column."""
    fig.canvas.draw(); r = fig.canvas.get_renderer(); inv = fig.transFigure.inverted()
    H = fig.get_size_inches()[1] / MM
    placed = []
    fixed = {}
    for row in rows:
        bbs = []
        for item in row:
            ax, L = item[0], item[1]; bbs.append((ax, L, ax.get_tightbbox(r).transformed(inv)))
            if len(item) > 2: fixed[L] = item[2]
        y = max(bb.y1 for _, _, bb in bbs)
        for ax, L, bb in bbs: placed.append([ax, L, bb.x0, y])
    # align x within columns: panels whose left edges differ by < 4 mm share the smallest x
    W = fig.get_size_inches()[0] / MM
    for p in placed:
        same = [q[2] for q in placed if abs(q[2] - p[2]) * W < 4]
        p[2] = fixed.get(p[1], min(same))
    for ax, L, x, y in placed:
        fig.text(max(x, 0.002), min(y + 0.4 / H, 1 - top_mm / H), L, fontsize=9, fontweight="bold", ha="left", va="top" if y + 0.4 / H > 1 - top_mm / H else "bottom")


def save(fig, name):
    from PIL import Image
    fig.savefig(FIG / f"{name}.pdf", facecolor="white")
    p = FIG / f"{name}.png"; fig.savefig(p, facecolor="white", dpi=600)
    Image.open(p).convert("RGB").save(p, dpi=(600, 600)); plt.close(fig)


def cbar(fig, im, ax, label, ticks=None, frac=0.05, pad=0.03, ticklabels=None, extend="neither"):
    cb = fig.colorbar(im, ax=ax, fraction=frac, pad=pad, extend=extend)
    cb.ax.tick_params(labelsize=6.5, width=0.4, length=1.5); cb.outline.set_linewidth(0.4); cb.set_label(label, fontsize=6.5)
    if ticks is not None: cb.set_ticks(ticks)
    if ticklabels is not None: cb.set_ticklabels(ticklabels)
    return cb


def screen_glyph(ax, x, y, screen, size=3.2, transform=None):
    """Small cell-line glyph (licensed icon, tinted in the screen color) centered at (x, y) in axes fraction."""
    from icons import icon
    from matplotlib.offsetbox import OffsetImage, AnnotationBbox
    name = {"K562": "myeloid_cell", "K562gw": "myeloid_cell", "RPE1": "epithelial_cell", "HepG2": "liver", "Jurkat": "lymphocyte"}[screen]
    im = icon(name, 160, tint=SCREEN_C[screen])
    ab = AnnotationBbox(OffsetImage(im, zoom=size * MM * 72 / im.shape[1]), (x, y), xycoords=transform or ax.transAxes,
                        frameon=False, box_alignment=(0.5, 0.5), annotation_clip=False)
    ax.add_artist(ab)


def no_overlap_check(fig):
    """Return pairs of text artists whose rendered boxes overlap (used in the build to fail loudly)."""
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    texts = [t for t in fig.findobj(matplotlib.text.Text) if t.get_visible() and t.get_text().strip()]
    boxes = [(t, t.get_window_extent(r)) for t in texts]
    bad = []
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i][1], boxes[j][1]
            if a.overlaps(b) and a.width > 1 and b.width > 1:
                bad.append((boxes[i][0].get_text()[:30], boxes[j][0].get_text()[:30]))
    return bad
