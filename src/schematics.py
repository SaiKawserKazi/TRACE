"""Schematic panels (Figure 1A, Figure 2A) and the graphical abstract, drawn in millimetres with licensed icons
(figures/icons/CREDITS.md) and simple shapes drawn here (cobblestone epithelium, worm outline). Arrows keep >= 1.5 mm from
boxes; text >= 6.5 pt in figures and >= 16 pt on the 6-inch abstract (>= 5 pt when it is printed 5 cm wide)."""
import numpy as np
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle, Circle, Wedge, RegularPolygon
import matplotlib.pyplot as plt
from icons import place
from figstyle import INK, GRAY, LIGHT, UP, DOWN, SCREEN_C, COHORT_C, MM

BOXEC, BOXFC = "#9e9e9e", "#fafafa"
SCORE_CMAP = "BrBG"     # brown (looks older) to teal (looks younger); does not clash with the blue-red correlation scale


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


def score_bar(ax, x, y, w, h, left="looks older", right="looks younger", size=6.5, tick=True):
    cm = plt.get_cmap(SCORE_CMAP)
    for i in range(60):
        ax.add_patch(Rectangle((x + i * w / 60, y), w / 60 + 0.01, h, fc=cm(0.12 + 0.76 * i / 59), ec="none", zorder=5))
    ax.add_patch(Rectangle((x, y), w, h, fc="none", ec=INK, lw=0.5, zorder=6))
    if tick: ax.plot([x + w / 2] * 2, [y - 0.8, y + h + 0.8], color=INK, lw=0.6, zorder=6)
    if left: T(ax, x, y - 2.2, left, size, ha="left")
    if right: T(ax, x + w, y - 2.2, right, size, ha="right")


def cc_wheel(ax, x, y, r, paused=True):
    """Cell-cycle wheel; paused = two bars in the middle."""
    for k, (a0, a1, c) in enumerate([(90, 200, "#9e9e9e"), (200, 290, "#bdbdbd"), (290, 360, "#7f7f7f"), (0, 90, "#cfcfcf")]):
        ax.add_patch(Wedge((x, y), r, a0, a1, width=r * 0.38, fc=c, ec="white", lw=0.4, zorder=5))
    if paused:
        for dx in (-0.32 * r, 0.32 * r): ax.add_patch(Rectangle((x + dx - 0.12 * r, y - 0.32 * r), 0.24 * r, 0.64 * r, fc=INK, ec="none", zorder=6))


def cobblestone(ax, x, y, w, color):
    """A small sheet of four joined hexagonal epithelial cells (RPE1), w wide; returns the height."""
    r = w / 4.0; dx, dy = np.sqrt(3) * r, 1.5 * r
    for cx, cy in [(-dx / 2, dy / 2), (dx / 2, dy / 2), (0, -dy / 2), (dx, -dy / 2)]:
        ax.add_patch(RegularPolygon((x + cx - dx / 4, y + cy), 6, radius=r, orientation=0, fc=color, ec="white", lw=0.5, alpha=0.9, zorder=5))
        ax.add_patch(Circle((x + cx - dx / 4, y + cy), 0.28 * r, fc="#3b3b3b", ec="none", zorder=6))
    return 2.0 * dy + r


def worm(ax, x, y, w, lw=1.6, color="#5a5a5a"):
    """Flat C. elegans outline: a sinuous rod with a rounded head."""
    t = np.linspace(0, 1, 60); xs = x - w / 2 + w * t; ys = y + 0.16 * w * np.sin(2 * np.pi * 1.2 * t)
    ax.plot(xs, ys, color=color, lw=lw, solid_capstyle="round", zorder=6)
    ax.add_patch(Circle((xs[-1], ys[-1]), 0.045 * w, fc=color, ec="none", zorder=6))


CELL_ICON = {"K562": "myeloid_cell", "HepG2": "hepatocyte_like", "Jurkat": "lymphocyte"}   # suspension, hepatocyte-like, T cell; RPE1 is drawn


def cell_glyph(ax, c, x, y, size):
    """One cell-line glyph in its screen color; returns its height."""
    if c == "RPE1": return cobblestone(ax, x, y, size * 1.15, SCREEN_C[c])
    return place(ax, CELL_ICON[c], x, y, size, tint=SCREEN_C[c])


def cell_row(ax, x0, y, dx, size, labels=True):
    """The four cell lines in their screen colors, each label 1 mm below its icon."""
    for i, c in enumerate(["K562", "RPE1", "HepG2", "Jurkat"]):
        h = cell_glyph(ax, c, x0 + i * dx, y, size)
        if labels: T(ax, x0 + i * dx, y - max(h, size) / 2 - 1.0 - 1.1, c, 6.5, color=INK)


def dcas9_krab(ax, x, y, w, size=6.5):
    """dCas9 (no cut marks) on DNA with its KRAB repressor domain; returns the icon height."""
    h = place(ax, "dcas9", x, y, w)
    kx, ky = x + 0.30 * w, y + 0.40 * h
    bw, bh = 7.8 * size / 6.5, 3.0 * size / 6.5      # the word is about 5.6 mm wide at 6.5 pt: about 1 mm of padding on each side
    ax.add_patch(FancyBboxPatch((kx - bw / 2, ky - bh / 2), bw, bh, boxstyle="round,pad=0,rounding_size=0.8", fc="#f4a582", ec="#b2182b", lw=0.5, zorder=7))
    T(ax, kx, ky, "KRAB", size, zorder=8)
    return h


def down_bars(ax, x, y, n=4, w=11, h=4.5):
    """Genes down with age: blue bars hanging below a baseline, as in the age-signature icon."""
    bw = w / n
    for i in range(n):
        hh = h * (1 - 0.15 * i); ax.add_patch(Rectangle((x + i * bw + 0.12 * bw, y - hh), 0.76 * bw, hh, fc=DOWN, ec="none", zorder=5))
    ax.plot([x, x + w], [y, y], color=INK, lw=0.6, zorder=6)


# ------------------------------------------------------------------ Figure 1A: study design
def fig1a(ax, w=174, h=42):
    mm_axes(ax, w, h, x0=-4)
    A = dict(col=GRAY, lw=0.9, ms=7)                                        # one arrow style
    # inputs: cohorts (top) and screens (bottom)
    box(ax, 1, 21.8, 62, 19.7); T(ax, 32, 39.6, "Aging cohorts", 7, fontweight="bold")
    for cx, (nm, n, ic, coh) in zip([11.5, 32, 52.5], [("Discovery", "n = 111", "fibroblast", "Discovery"), ("GTEx fibroblasts", "n = 652", "fibroblast", "GTEx fibroblasts"),
                                                      ("GTEx skin", "n = 651", "skin", "GTEx skin")]):
        if ic == "fibroblast": place(ax, "fibroblast", cx, 35.2, 8.0); place(ax, "dish", cx, 32.6, 5.4)
        else: place(ax, "skin", cx, 34.0, 6.4)
        ax.add_patch(Circle((cx, 29.2), 0.85, fc=COHORT_C[coh], ec="none", zorder=6))     # cohort dot under its icon
        T(ax, cx, 26.9, nm, 6.5); T(ax, cx, 24.1, n, 6.5, color=GRAY)
    box(ax, 1, 0.5, 62, 19.8); T(ax, 32, 18.4, "CRISPRi Perturb-seq screens", 7, fontweight="bold")
    dcas9_krab(ax, 13.5, 11.8, 15.0); T(ax, 13.5, 3.9, "dCas9-KRAB\nsilences one gene", 6.5, linespacing=1.05)
    for i, c in enumerate(["K562", "RPE1", "HepG2", "Jurkat"]):
        cell_glyph(ax, c, 32.5 + i * 7.6, 13.4, 4.4); T(ax, 32.5 + i * 7.6, 8.9, c, 6.5)
    T(ax, 44, 3.9, "5 screens, 4 cell lines\n18,746 knockdown profiles", 6.5, linespacing=1.05)
    # age signature
    arrow(ax, 64.8, 31.6, 72.8, 31.6, **A)
    box(ax, 74.5, 21.8, 24.5, 19.7); T(ax, 86.75, 39.6, "Age signature", 7, fontweight="bold")
    updown_bars(ax, 79.75, 31.4, 14, 9.5)
    T(ax, 86.75, 24.4, "250 up, 250 down", 6.5)
    # reversal score
    arrow(ax, 100.8, 31.6, 104.8, 31.6, **A)
    arrow(ax, 64.8, 10.4, 104.8, 10.4, **A); T(ax, 84.8, 12.5, "expression shift per knockdown", 6.5, color=GRAY)
    box(ax, 106.5, 0.5, 27, 41.0); T(ax, 120, 39.6, "Reversal score", 7, fontweight="bold")
    score_bar(ax, 110, 28.0, 20, 3.6, "older", "younger")
    T(ax, 120, 15.0, "one score per\nknockdown", 6.5, color=GRAY, linespacing=1.05)
    # five tests
    arrow(ax, 135.3, 21.0, 139.3, 21.0, **A)
    box(ax, 141, 0.5, 32, 41.0); T(ax, 157, 39.6, "Five tests", 7, fontweight="bold")
    for i, lab in enumerate(["Cell lines", "Proliferation", "Network", "Labels", "Cohorts"]):
        y = 33.4 - i * 7.0; ix = 147.2
        if lab == "Cell lines":
            for k, c in enumerate(["K562", "RPE1", "HepG2", "Jurkat"]): ax.add_patch(Circle((143.5 + k * 2.45, y), 1.05, fc=SCREEN_C[c], ec="white", lw=0.3, zorder=6))
        elif lab == "Proliferation": cc_wheel(ax, ix, y, 2.4)
        elif lab == "Network": place(ax, "network", ix, y, 5.6)
        elif lab == "Labels": worm(ax, ix - 2.2, y, 3.6, lw=1.0); place(ax, "mouse", ix + 2.2, y, 3.6)
        else: place(ax, "old_man", ix - 1.4, y, 2.0); place(ax, "old_woman", ix + 1.4, y, 2.0)
        T(ax, 152.6, y, lab, 6.5, ha="left")


# ------------------------------------------------------------------ Figure 2A: two lanes that meet (working model)
def fig2a(ax, w=174, h=40):
    mm_axes(ax, w, h, x0=-2)
    A = dict(col=GRAY, lw=0.9, ms=7)
    yt, yb = 30.0, 9.5
    T(ax, 4.5, yt + 7.6, "Aging", 7, ha="left", fontweight="bold"); T(ax, 4.5, yb + 7.6, "Knockdown", 7, ha="left", fontweight="bold")
    ax.plot([4.5, 92], [20.2, 20.2], color=LIGHT, lw=0.5, zorder=2)                  # divider between the lanes
    ly_t, ly_b = yt - 7.0, yb - 7.0                                                 # one label baseline per lane
    # aging lane: young spindle fibroblasts with one dividing -> old, large, flat cell -> paused cycle -> cell-cycle genes down
    place(ax, "fibroblast_5", 12.5, yt + 2.4, 8.5); place(ax, "fibroblast_5", 15.0, yt - 1.6, 8.5)
    place(ax, "fibroblast_5", 22.0, yt + 1.2, 5.6); place(ax, "fibroblast_5", 23.5, yt - 1.0, 5.6)   # a pair of daughter cells
    T(ax, 17.5, ly_t, "young: cycling\nfibroblasts", 6.5, linespacing=1.05)
    arrow(ax, 29, yt, 35, yt, **A)
    place(ax, "fibroblast_4", 45, yt, 14.5); T(ax, 45, ly_t, "old: large, flat,\nfewer cycling", 6.5, linespacing=1.05)
    arrow(ax, 55, yt, 61, yt, **A)
    cc_wheel(ax, 66.5, yt, 3.4); T(ax, 66.5, ly_t, "cell cycle\npaused", 6.5, linespacing=1.05)
    arrow(ax, 72, yt, 78, yt, **A)
    down_bars(ax, 79.5, yt + 2.5); T(ax, 85, ly_t, "cell-cycle\ngenes down", 6.5, linespacing=1.05)
    # knockdown lane: dCas9-KRAB switches an essential gene off -> slower-growing cell -> paused cycle -> the same genes down
    dcas9_krab(ax, 15, yb + 1.2, 15.0); T(ax, 15, ly_b, "dCas9-KRAB silences\nan essential gene", 6.5, linespacing=1.05)
    arrow(ax, 29, yb, 35, yb, **A)
    place(ax, "lymphocyte", 45, yb, 6.2, tint="#8c8c8c"); T(ax, 45, ly_b, "slower-growing\ncells", 6.5, linespacing=1.05)
    arrow(ax, 55, yb, 61, yb, **A)
    cc_wheel(ax, 66.5, yb, 3.4); T(ax, 66.5, ly_b, "cell cycle\npaused", 6.5, linespacing=1.05)
    arrow(ax, 72, yb, 78, yb, **A)
    down_bars(ax, 79.5, yb + 2.5); T(ax, 85, ly_b, "same genes\ndown", 6.5, linespacing=1.05)
    # the lanes meet
    arrow(ax, 93, yt, 101, 24.0, rad=-0.2, **A); arrow(ax, 93, yb, 101, 16.0, rad=0.2, **A)
    box(ax, 103, 15.0, 23, 10, fc="white", ec=INK, lw=0.6); T(ax, 114.5, 20.0, "same genes,\nsame direction", 6.5, linespacing=1.05)
    arrow(ax, 127.5, 20.0, 131.5, 20.0, **A)
    x0, wb = 134, 38
    score_bar(ax, x0, 19.0, wb, 3.6, "", "")
    T(ax, x0 + wb / 2, 25.0, "0", 6.5)
    T(ax, x0 + wb / 4, 15.0, "R < 0\nlooks older", 6.5, linespacing=1.05); T(ax, x0 + 3 * wb / 4, 15.0, "R > 0\nlooks younger", 6.5, linespacing=1.05)
    ax.annotate("", xy=(x0 + 6, 23.2), xytext=(x0 + 6, 29.4), arrowprops=dict(arrowstyle="-|>", color=GRAY, lw=0.9, mutation_scale=7))
    T(ax, x0 + 6, 32.6, "most essential-gene\nknockdowns", 6.5, linespacing=1.05)
    T(ax, x0 + wb / 2, 8.0, "R > 0 only if the knockdown\nkeeps cells growing", 6.5, linespacing=1.05)


# ------------------------------------------------------------------ graphical abstract (square, about 20 words)
def graphical_abstract(path):
    from PIL import Image
    S = 152.4  # 6 inches in mm
    fig = plt.figure(figsize=(6, 6)); ax = fig.add_axes([0, 0, 1, 1]); mm_axes(ax, S, S)
    big, mid = 22, 16
    A = dict(col=INK, lw=2.2, ms=18)
    # row 1: inputs
    box(ax, 4, 104, 68, 44, fc="white"); T(ax, 38, 139.5, "Aging cohorts", big, fontweight="bold")
    for i, ic in enumerate(["fibroblast", "fibroblast", "skin"]):
        place(ax, ic, 17 + i * 21, 118, 17 if ic == "fibroblast" else 14)
    box(ax, 80.4, 104, 68, 44, fc="white"); T(ax, 114.4, 139.5, "CRISPRi screens", 20, fontweight="bold")
    dcas9_krab(ax, 98, 119, 26, size=13)
    for i, c in enumerate(["K562", "RPE1", "HepG2", "Jurkat"]):
        cell_glyph(ax, c, 125 + (i % 2) * 13.5, 125 - (i // 2) * 13.5, 9.5)
    # row 2: score
    arrow(ax, 38, 101, 66, 87, **A); arrow(ax, 114.4, 101, 86.4, 87, **A)
    T(ax, 76.2, 80, "Reversal score", big, fontweight="bold")
    score_bar(ax, 36, 66, 80.4, 7, "", "", size=big)
    T(ax, 36, 61, "older", mid, ha="left", color=GRAY); T(ax, 116.4, 61, "younger", mid, ha="right", color=GRAY)
    for x in (27.0, 76.2, 125.4): arrow(ax, 76.2, 55, x, 49, **A)
    # row 3: findings (labels 16 pt with normal line spacing, 3 mm from the box border)
    bw = 46.0
    for x, lab in [(3.2, "Ranking\nchanges\nwith cohort"), (53.2, "Score tracks\ncell-cycle\ngenes"), (103.2, "Labels follow\nproliferation")]:
        box(ax, x, 3, bw, 43, fc="#f5f5f5"); T(ax, x + bw / 2, 6.0, lab, mid, fontweight="bold", linespacing=1.1, va="bottom")
    place(ax, "old_man", 17, 36.0, 5.6); place(ax, "old_woman", 35, 36.0, 6.0)                   # two cohorts, crossed arrows
    arrow(ax, 21, 38.5, 31, 30.5, col=INK, lw=1.6, ms=12); arrow(ax, 31, 38.5, 21, 30.5, col=INK, lw=1.6, ms=12)
    place(ax, "fibroblast_5", 66, 37.0, 13); place(ax, "fibroblast_5", 67, 31.5, 13)              # a dividing cell slows
    arrow(ax, 75.5, 34.5, 81.5, 34.5, col=INK, lw=1.6, ms=12); cc_wheel(ax, 88.5, 34.5, 5.2)
    worm(ax, 117, 35.5, 9, lw=2.4); place(ax, "mouse", 129.5, 35.5, 9)                           # labels follow the proliferation baseline
    T(ax, 138.5, 35.5, "=", mid, fontweight="bold"); cc_wheel(ax, 144.0, 35.5, 3.6, paused=False)
    fig.savefig(str(path) + ".pdf", facecolor="white"); fig.savefig(str(path) + ".png", dpi=300, facecolor="white")
    Image.open(str(path) + ".png").convert("RGB").save(str(path) + ".png", dpi=(300, 300))
    Image.open(str(path) + ".png").save(str(path) + ".tif", dpi=(300, 300), compression="tiff_lzw")
    Image.open(str(path) + ".png").resize((591, 591), Image.LANCZOS).save(str(path) + "_5cm_preview.png")
    plt.close(fig)
