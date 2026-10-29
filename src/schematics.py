"""Schematic panels (Figure 1A, Figure 2A) and the graphical abstract, drawn in millimetres with licensed icons
(figures/icons/CREDITS.md). Arrows keep >= 2 mm from boxes; text >= 6.5 pt in figures and >= 22 pt on the 6-inch abstract."""
import numpy as np
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Circle, Wedge
import matplotlib.pyplot as plt
from icons import place
from figstyle import INK, GRAY, LIGHT, UP, DOWN, SCREEN_C, COHORT_C, MM

BOXEC, BOXFC = "#9e9e9e", "#fafafa"


def mm_axes(ax, w, h, x0=0):
    ax.set_xlim(x0, w); ax.set_ylim(0, h); ax.set_aspect("equal"); ax.axis("off")


def box(ax, x, y, w, h, fc=BOXFC, ec=BOXEC, lw=0.6, r=1.2):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}", fc=fc, ec=ec, lw=lw, zorder=1))


def arrow(ax, x0, y0, x1, y1, col=GRAY, lw=0.9, ms=7, rad=0.0):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=ms, color=col, lw=lw, shrinkA=0, shrinkB=0,
                                 connectionstyle=f"arc3,rad={rad}", zorder=4))


def T(ax, x, y, s, size=6.5, **k):
    k.setdefault("ha", "center"); k.setdefault("va", "center"); k.setdefault("color", INK); k.setdefault("zorder", 6)
    ax.text(x, y, s, fontsize=size, **k)


def updown_bars(ax, x, y, w, h, n=8):
    """Icon of a signature: genes up (red) and down (blue) with age."""
    bw = w / n
    hs = np.r_[np.linspace(1.0, 0.45, n // 2), np.linspace(0.45, 1.0, n // 2)]
    for i, hh in enumerate(hs):
        up = i < n // 2
        ax.add_patch(Rectangle((x + i * bw + bw * 0.12, y if up else y - hh * h / 2), bw * 0.76, hh * h / 2, fc=UP if up else DOWN, ec="none", zorder=5))
    ax.plot([x, x + w], [y, y], color=INK, lw=0.6, zorder=6)


def score_bar(ax, x, y, w, h, left="looks older", right="looks younger", size=6.5):
    cm = plt.get_cmap("Greys")
    for i in range(60):
        ax.add_patch(Rectangle((x + i * w / 60, y), w / 60 + 0.01, h, fc=cm(0.75 - 0.6 * i / 59), ec="none", zorder=5))
    ax.add_patch(Rectangle((x, y), w, h, fc="none", ec=INK, lw=0.5, zorder=6))
    ax.plot([x + w / 2] * 2, [y - 0.8, y + h + 0.8], color=INK, lw=0.6, zorder=6)
    T(ax, x, y - 2.2, left, size, ha="left"); T(ax, x + w, y - 2.2, right, size, ha="right")


def cc_wheel(ax, x, y, r, paused=True):
    """Cell-cycle wheel; paused = two bars in the middle."""
    for k, (a0, a1, c) in enumerate([(90, 200, "#9e9e9e"), (200, 290, "#bdbdbd"), (290, 360, "#7f7f7f"), (0, 90, "#cfcfcf")]):
        ax.add_patch(Wedge((x, y), r, a0, a1, width=r * 0.38, fc=c, ec="white", lw=0.4, zorder=5))
    if paused:
        for dx in (-0.32 * r, 0.32 * r): ax.add_patch(Rectangle((x + dx - 0.12 * r, y - 0.32 * r), 0.24 * r, 0.64 * r, fc=INK, ec="none", zorder=6))


CELL_ICON = {"K562": "myeloid_cell", "RPE1": "cuboid_cell_2", "HepG2": "hepatocyte_like", "Jurkat": "lymphocyte"}   # suspension, epithelial, hepatocyte-like, T cell


def cell_row(ax, x0, y, dx, size, labels=True):
    """The four cell lines in their screen colors, each label 1 mm below its icon."""
    for i, c in enumerate(["K562", "RPE1", "HepG2", "Jurkat"]):
        h = place(ax, CELL_ICON[c], x0 + i * dx, y, size, tint=SCREEN_C[c])
        if labels: T(ax, x0 + i * dx, y - h / 2 - 1.0 - 1.1, c, 6.5, color=INK)


def dcas9_krab(ax, x, y, w):
    """dCas9 (no cut marks) on DNA with its KRAB repressor domain; returns the icon height."""
    h = place(ax, "dcas9", x, y, w)
    kx, ky = x + 0.30 * w, y + 0.36 * h
    ax.add_patch(FancyBboxPatch((kx - 3.2, ky - 1.3), 6.4, 2.6, boxstyle="round,pad=0,rounding_size=0.8", fc="#f4a582", ec="#b2182b", lw=0.5, zorder=7))
    T(ax, kx, ky, "KRAB", 6.5, zorder=8)
    return h


def down_bars(ax, x, y, n=4, w=11, h=4.5):
    """Genes down with age: blue bars hanging below a baseline, as in the age-signature icon."""
    bw = w / n
    for i in range(n):
        hh = h * (1 - 0.15 * i); ax.add_patch(Rectangle((x + i * bw + 0.12 * bw, y - hh), 0.76 * bw, hh, fc=DOWN, ec="none", zorder=5))
    ax.plot([x, x + w], [y, y], color=INK, lw=0.6, zorder=6)


# ------------------------------------------------------------------ Figure 1A: study design
def fig1a(ax, w=174, h=46):
    mm_axes(ax, w, h, x0=-4)
    A = dict(col=GRAY, lw=0.9, ms=7)                                        # one arrow style
    # inputs: cohorts (top) and screens (bottom)
    box(ax, 1, 24, 62, 21.5); T(ax, 32, 43.2, "Aging cohorts", 7, fontweight="bold")
    for cx, (nm, n, ic, coh) in zip([11.5, 32, 52.5], [("Discovery", "n = 111", "fibroblast", "Discovery"), ("GTEx fibroblasts", "n = 652", "fibroblast", "GTEx fibroblasts"),
                                                      ("GTEx skin", "n = 651", "skin", "GTEx skin")]):
        if ic == "fibroblast": place(ax, "fibroblast", cx, 38.4, 8.5); place(ax, "dish", cx, 35.0, 5.6)
        else: place(ax, "skin", cx, 36.6, 7.0)
        tw = 0.59 * len(nm) * 6.5 / 2.835 / 1.0
        ax.add_patch(Circle((cx - tw / 2 - 1.6, 30.2), 0.85, fc=COHORT_C[coh], ec="none", zorder=6))
        T(ax, cx, 30.2, nm, 6.5); T(ax, cx, 26.9, n, 6.5, color=GRAY)
    box(ax, 1, 0.5, 62, 21.5); T(ax, 32, 19.7, "CRISPRi Perturb-seq screens", 7, fontweight="bold")
    dcas9_krab(ax, 13.5, 12.0, 15.5); T(ax, 13.5, 4.2, "dCas9-KRAB\nsilences one gene", 6.5, linespacing=1.05)
    cell_row(ax, 32.5, 14.0, 7.6, 4.6)
    T(ax, 44, 4.4, "5 screens, 4 cell lines\n18,746 knockdown profiles", 6.5, linespacing=1.05)
    # age signature
    arrow(ax, 65, 34.75, 73, 34.75, **A)
    box(ax, 75, 24, 24, 21.5); T(ax, 87, 43.2, "Age signature", 7, fontweight="bold")
    updown_bars(ax, 80, 34.5, 14, 10.5)
    T(ax, 87, 26.8, "250 up, 250 down", 6.5)
    # reversal score
    arrow(ax, 101, 34.75, 105, 34.75, **A)
    arrow(ax, 65, 11.25, 105, 11.25, **A); T(ax, 85, 13.6, "expression shift per knockdown", 6.5, color=GRAY)
    box(ax, 107, 0.5, 26, 45); T(ax, 120, 43.2, "Reversal score", 7, fontweight="bold")
    score_bar(ax, 110, 24, 20, 3.6, "older", "younger")
    T(ax, 120, 15.5, "one score per\nknockdown", 6.5, color=GRAY, linespacing=1.05)
    # five tests
    arrow(ax, 135, 23, 139, 23, **A)
    box(ax, 141, 0.5, 32, 45); T(ax, 157, 43.2, "Five tests", 7, fontweight="bold")
    for i, lab in enumerate(["Cell lines", "Proliferation", "Network", "Labels", "Cohorts"]):
        y = 36.6 - i * 7.6; ix = 147.5
        if lab == "Cell lines":
            for k, c in enumerate(["K562", "RPE1", "HepG2", "Jurkat"]): ax.add_patch(Circle((143.9 + k * 2.4, y), 1.05, fc=SCREEN_C[c], ec="white", lw=0.3, zorder=6))
        elif lab == "Proliferation": cc_wheel(ax, ix, y, 2.5)
        elif lab == "Network": place(ax, "network", ix, y, 6.0)
        elif lab == "Labels": place(ax, "worm", ix - 2.0, y, 3.6); place(ax, "mouse", ix + 2.2, y, 3.8)
        else: place(ax, "old_man", ix - 1.4, y, 2.4); place(ax, "old_woman", ix + 1.4, y, 2.4)
        T(ax, 153, y, lab, 6.5, ha="left")


# ------------------------------------------------------------------ Figure 2A: two lanes that meet (working model)
def fig2a(ax, w=174, h=46):
    mm_axes(ax, w, h, x0=-2)
    A = dict(col=GRAY, lw=0.9, ms=7)
    T(ax, 4, 43.8, "Aging", 7, ha="left", fontweight="bold"); T(ax, 4, 20.8, "Knockdown", 7, ha="left", fontweight="bold")
    yt, yb = 33.5, 11.5
    # aging lane: young spindle fibroblasts with one dividing -> old, large, flat cell -> paused cycle -> cell-cycle genes down
    place(ax, "fibroblast_5", 12.5, yt + 2.8, 9.0); place(ax, "fibroblast_5", 15.0, yt - 1.8, 9.0)
    place(ax, "fibroblast_5", 22.0, yt + 1.4, 6.0); place(ax, "fibroblast_5", 23.5, yt - 1.2, 6.0)   # a pair of daughter cells
    T(ax, 17.5, yt - 8.6, "young: cycling\nfibroblasts", 6.5, linespacing=1.05)
    arrow(ax, 29, yt, 35, yt, **A)
    place(ax, "fibroblast_4", 45, yt, 15.5); T(ax, 45, yt - 8.6, "old: large, flat,\nfewer cycling", 6.5, linespacing=1.05)
    arrow(ax, 55, yt, 61, yt, **A)
    cc_wheel(ax, 66.5, yt, 3.5); T(ax, 66.5, yt - 8.6, "cell cycle\npaused", 6.5, linespacing=1.05)
    arrow(ax, 72, yt, 78, yt, **A)
    down_bars(ax, 79.5, yt + 2.5); T(ax, 85, yt - 8.6, "cell-cycle\ngenes down", 6.5, linespacing=1.05)
    # knockdown lane: dCas9-KRAB switches an essential gene off -> slower-growing cell -> paused cycle -> the same genes down
    dcas9_krab(ax, 15, yb + 1.2, 15.5); T(ax, 15, yb - 8.4, "dCas9-KRAB silences\nan essential gene", 6.5, linespacing=1.05)
    arrow(ax, 29, yb, 35, yb, **A)
    place(ax, "lymphocyte", 45, yb, 6.5, tint="#8c8c8c"); T(ax, 45, yb - 8.4, "slower-growing\ncells", 6.5, linespacing=1.05)
    arrow(ax, 55, yb, 61, yb, **A)
    cc_wheel(ax, 66.5, yb, 3.5); T(ax, 66.5, yb - 8.4, "cell cycle\npaused", 6.5, linespacing=1.05)
    arrow(ax, 72, yb, 78, yb, **A)
    down_bars(ax, 79.5, yb + 2.5); T(ax, 85, yb - 8.4, "same genes\ndown", 6.5, linespacing=1.05)
    # the lanes meet
    arrow(ax, 93, yt, 101, 26.5, rad=-0.2, **A); arrow(ax, 93, yb, 101, 18.5, rad=0.2, **A)
    box(ax, 103, 17.5, 23, 10, fc="white", ec=INK, lw=0.6); T(ax, 114.5, 22.5, "same genes,\nsame direction", 6.5, linespacing=1.05)
    arrow(ax, 128, 22.5, 132, 22.5, **A)
    x0, wb = 134, 38
    score_bar(ax, x0, 21, wb, 3.6, "", "")
    T(ax, x0 + wb / 2, 27.0, "0", 6.5)
    T(ax, x0 + wb / 4, 17.4, "R < 0\nlooks older", 6.5, linespacing=1.05); T(ax, x0 + 3 * wb / 4, 17.4, "R > 0\nlooks younger", 6.5, linespacing=1.05)
    ax.annotate("", xy=(x0 + 6, 25.2), xytext=(x0 + 6, 31.2), arrowprops=dict(arrowstyle="-|>", color=GRAY, lw=0.9, mutation_scale=7))
    T(ax, x0 + 6, 34.4, "most essential-gene\nknockdowns", 6.5, linespacing=1.05)
    T(ax, x0 + wb / 2, 10.6, "R > 0 only if the knockdown\nkeeps cells growing", 6.5, linespacing=1.05)
    T(ax, x0 + wb / 2, 3.0, "Holds for a signature with a cell-\ncycle arm; the GTEx-derived\nsignature has none", 6.5, color=GRAY, linespacing=1.05, style="italic")


# ------------------------------------------------------------------ graphical abstract (square, about 20 words)
def graphical_abstract(path):
    from PIL import Image
    S = 152.4  # 6 inches in mm
    fig = plt.figure(figsize=(6, 6)); ax = fig.add_axes([0, 0, 1, 1]); mm_axes(ax, S, S)
    big, mid = 23, 21
    A = dict(col=INK, lw=2.2, ms=18)
    # row 1: inputs
    box(ax, 4, 104, 68, 44, fc="white"); T(ax, 38, 140, "Aging cohorts", big, fontweight="bold")
    for i, ic in enumerate(["fibroblast", "fibroblast", "skin"]):
        place(ax, ic, 17 + i * 21, 121, 18 if ic == "fibroblast" else 15)
    for i, coh in enumerate(["Discovery", "GTEx fibroblasts", "GTEx skin"]):
        ax.add_patch(Circle((17 + i * 21, 108.5), 2.4, fc=COHORT_C[coh], ec="none", zorder=6))
    box(ax, 80.4, 104, 68, 44, fc="white"); T(ax, 114.4, 140, "CRISPRi screens", big, fontweight="bold")
    dcas9_krab(ax, 98, 121, 28)
    for i, c in enumerate(["K562", "RPE1", "HepG2", "Jurkat"]):
        place(ax, CELL_ICON[c], 124 + (i % 2) * 14, 127 - (i // 2) * 14, 10.5, tint=SCREEN_C[c])
    # row 2: score
    arrow(ax, 38, 101, 66, 87, **A); arrow(ax, 114.4, 101, 86.4, 87, **A)
    T(ax, 76.2, 80, "Reversal score", big, fontweight="bold")
    score_bar(ax, 36, 66, 80.4, 7, "", "", size=big)
    T(ax, 36, 61, "older", mid, ha="left", color=GRAY); T(ax, 116.4, 61, "younger", mid, ha="right", color=GRAY)
    for x in (25.4, 76.2, 127): arrow(ax, 76.2, 55, x, 48, **A)
    # row 3: findings
    for x, lab in [(4, "Cohort\ndecides\nranking"), (54.8, "Growth\ndrives\nscore"), (105.6, "Labels\nmatch\ngrowth")]:
        box(ax, x, 3, 42.8, 42, fc="#f5f5f5"); T(ax, x + 21.4, 13.5, lab, mid, fontweight="bold", linespacing=1.0)
    place(ax, "old_man", 16, 35, 7); place(ax, "old_woman", 34, 35, 7.5)                       # two cohorts, crossed arrows
    arrow(ax, 20, 38, 30, 29, col=INK, lw=1.6, ms=12); arrow(ax, 30, 38, 20, 29, col=INK, lw=1.6, ms=12)
    place(ax, "mitosis", 65, 34, 13); arrow(ax, 73, 34, 78, 34, col=INK, lw=1.6, ms=12)          # dividing cell slows
    cc_wheel(ax, 84, 34, 4.2); arrow(ax, 89.5, 34, 93, 34, col=INK, lw=1.6, ms=12)
    score_bar(ax, 94, 32.5, 4.5, 3, "", "")
    place(ax, "worm", 117, 35, 8); place(ax, "mouse", 128, 35, 8.5)                              # labels track the growth baseline
    T(ax, 137.5, 35, "=", mid, fontweight="bold"); cc_wheel(ax, 143.5, 35, 3.6, paused=False)
    fig.savefig(str(path) + ".pdf", facecolor="white"); fig.savefig(str(path) + ".png", dpi=300, facecolor="white")
    Image.open(str(path) + ".png").convert("RGB").save(str(path) + ".png", dpi=(300, 300))
    Image.open(str(path) + ".png").save(str(path) + ".tif", dpi=(300, 300), compression="tiff_lzw")
    Image.open(str(path) + ".png").resize((591, 591), Image.LANCZOS).save(str(path) + "_5cm_preview.png")
    plt.close(fig)
