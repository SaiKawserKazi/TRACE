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


def cell_glyphs(ax, x, y, size):
    for i, (c, ic) in enumerate([("K562", "myeloid_cell"), ("RPE1", "fibroblast_flat"), ("HepG2", "liver"), ("Jurkat", "lymphocyte")]):
        cx = x + (i % 2) * size * 1.9; cy = y - (i // 2) * size * 1.45
        place(ax, ic, cx, cy, size, tint=SCREEN_C[c])
        T(ax, cx, cy - size * 0.62, {"K562": "K562", "RPE1": "RPE1", "HepG2": "HepG2", "Jurkat": "Jurkat"}[c], 6.5, color=SCREEN_C[c])


# ------------------------------------------------------------------ Figure 1A: study design
def fig1a(ax, w=174, h=46):
    mm_axes(ax, w, h, x0=-4.5)
    # inputs: cohorts (top) and screens (bottom)
    box(ax, 0.5, 24.5, 64, 21); T(ax, 3, 43, "Aging cohorts", 7, ha="left", fontweight="bold")
    rows = [("Discovery fibroblasts", "n = 111 donors", "fibroblast", "Discovery"), ("GTEx fibroblasts", "n = 652 donors", "fibroblast", "GTEx fibroblasts"),
            ("GTEx skin", "n = 651 donors", "skin", "GTEx skin")]
    for i, (nm, n, ic, coh) in enumerate(rows):
        cx = 8 + i * 20.5
        place(ax, ic, cx + 2, 36.0 if ic == "fibroblast" else 34.0, 10 if ic == "fibroblast" else 8.0, tint=None)
        if ic == "fibroblast": place(ax, "dish", cx + 2, 32.0, 6.0)
        ax.add_patch(Circle((cx - 5.6, 26.8), 0.9, fc=COHORT_C[coh], ec="none", zorder=6))
        T(ax, cx + 2, 26.8, nm.replace(" fibroblasts", "\nfibroblasts"), 6.5, linespacing=1.0)
        T(ax, cx + 2, 40.4, n, 6.5, color=GRAY)
    box(ax, 0.5, 0.5, 64, 21.5); T(ax, 3, 19.5, "CRISPRi Perturb-seq screens", 7, ha="left", fontweight="bold")
    place(ax, "cas9", 14, 10.5, 20); T(ax, 14, 3.4, "dCas9-KRAB silences\none gene per cell", 6.5, linespacing=1.05)
    cell_glyphs(ax, 31, 15.2, 4.6)
    T(ax, 55.5, 10.5, "5 screens\n4 cell lines\n18,746\nknockdown\nprofiles", 6.5, linespacing=1.1)
    # age signature
    arrow(ax, 67, 35, 74.5, 35)
    box(ax, 76.5, 24.5, 26, 21); T(ax, 89.5, 43, "Age signature", 7, fontweight="bold")
    updown_bars(ax, 81, 34, 17, 12)
    T(ax, 89.5, 27.2, "250 up, 250 down", 6.5)
    # reversal score
    arrow(ax, 104.5, 35, 108, 28.5); arrow(ax, 67, 11, 108, 20.5)
    box(ax, 110, 13.5, 26, 21); T(ax, 123, 32, "Reversal score", 7, fontweight="bold")
    score_bar(ax, 113, 22, 20, 3.4, "older", "younger")
    T(ax, 123, 16.3, "per knockdown", 6.5, color=GRAY)
    # five tests
    arrow(ax, 138, 24, 141.5, 24)
    box(ax, 143.5, 0.5, 30, 45); T(ax, 158.5, 43, "Five tests", 7, fontweight="bold")
    tests = [("network", "Cell lines"), ("mitosis", "Proliferation"), ("network", "Network"), ("fibroblast_flat", "Labels"), ("old_man", "Cohorts")]
    for i, (ic, lab) in enumerate(tests):
        y = 36.5 - i * 7.6
        if lab == "Cell lines":
            for k, c in enumerate(["K562", "HepG2", "Jurkat"]): place(ax, "lymphocyte", 147.5 + k * 2.4, y, 2.6, tint=SCREEN_C[c])
        elif lab == "Labels":
            place(ax, "fibroblast_flat", 149, y, 7)
        elif lab == "Cohorts":
            place(ax, "old_man", 148, y, 2.6); place(ax, "old_woman", 151, y, 2.6)
        elif lab == "Proliferation":
            place(ax, "mitosis", 149.5, y, 6)
        else:
            place(ax, "network", 149.5, y, 6.5)
        T(ax, 155, y, lab, 6.5, ha="left")


# ------------------------------------------------------------------ Figure 2A: two lanes that meet
def fig2a(ax, w=174, h=46):
    mm_axes(ax, w, h)
    T(ax, 5, 43.5, "Aging", 7, ha="left", fontweight="bold"); T(ax, 5, 20.5, "Knockdown", 7, ha="left", fontweight="bold")
    yt, yb = 34, 12
    # top lane
    place(ax, "mitosis", 22, yt + 1.5, 11); T(ax, 22, yt - 8.5, "young: cycling\nfibroblasts", 6.5, linespacing=1.05)
    arrow(ax, 31, yt, 37, yt)
    place(ax, "fibroblast_flat", 47, yt, 15); T(ax, 47, yt - 8, "old: flat, fewer\ncycling cells", 6.5, linespacing=1.05)
    arrow(ax, 57, yt, 62, yt)
    cc_wheel(ax, 68, yt, 3.6); T(ax, 68, yt - 8, "cell cycle\nslowed", 6.5, linespacing=1.05)
    arrow(ax, 74, yt, 79, yt)
    [ax.add_patch(Rectangle((81 + i * 2.25, yt - 4), 1.6, 4 - i * 0.6, fc=DOWN, ec="none", zorder=5)) for i in range(4)]
    T(ax, 85.5, yt - 8, "cell-cycle\ngenes down", 6.5, linespacing=1.05)
    # bottom lane
    place(ax, "cas9", 22, yb + 1.5, 15); T(ax, 22, yb - 5.2, "dCas9-KRAB silences\nan essential gene", 6.5, linespacing=1.05)
    arrow(ax, 31, yb, 37, yb)
    place(ax, "myeloid_cell", 47, yb + 1, 7, tint="#7f7f7f"); T(ax, 47, yb - 5.2, "slower-growing\ncells", 6.5, linespacing=1.05)
    arrow(ax, 57, yb, 62, yb)
    cc_wheel(ax, 68, yb + 1, 3.6); T(ax, 68, yb - 5.2, "cell cycle\nslowed", 6.5, linespacing=1.05)
    arrow(ax, 74, yb, 79, yb)
    [ax.add_patch(Rectangle((81 + i * 2.25, yb - 2.5), 1.6, 4 - i * 0.6, fc=DOWN, ec="none", zorder=5)) for i in range(4)]
    T(ax, 85.5, yb - 5.2, "same genes\ndown", 6.5, linespacing=1.05)
    # lanes meet
    arrow(ax, 92, yt - 1, 102, 25.5, rad=-0.25); arrow(ax, 92, yb + 1, 102, 20.5, rad=0.25)
    box(ax, 103.5, 18, 24, 10, fc="white", ec=INK, lw=0.6)
    T(ax, 115.5, 23, "same genes,\nsame direction", 6.5, linespacing=1.05)
    arrow(ax, 129.5, 23, 134, 23)
    score_bar(ax, 136, 21.5, 36, 3.6, "looks older", "looks younger")
    ax.annotate("", xy=(141, 25.6), xytext=(141, 30.5), arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.8, mutation_scale=6))
    T(ax, 141, 33, "most essential-gene\nknockdowns", 6.5, linespacing=1.05)
    T(ax, 154, 13.4, "R > 0 only if the knockdown\nkeeps cells growing", 6.5, linespacing=1.05)
    T(ax, 154, 4.6, "Holds for a signature with a\ncell-cycle arm; the GTEx-derived\nsignature has none", 6.5, color=GRAY, linespacing=1.05, style="italic")


# ------------------------------------------------------------------ graphical abstract (square, ~15 words)
def graphical_abstract(path):
    from PIL import Image
    S = 152.4  # 6 inches in mm
    fig = plt.figure(figsize=(6, 6)); ax = fig.add_axes([0, 0, 1, 1]); mm_axes(ax, S, S)
    big, mid = 23, 21
    # row 1: inputs
    box(ax, 4, 104, 68, 44, fc="white"); T(ax, 38, 140, "Aging cohorts", big, fontweight="bold")
    for i, (ic, tint) in enumerate([("fibroblast", None), ("fibroblast", None), ("skin", None)]):
        place(ax, ic, 17 + i * 21, 121, 18 if ic == "fibroblast" else 15)
    for i, coh in enumerate(["Discovery", "GTEx fibroblasts", "GTEx skin"]):
        ax.add_patch(Circle((17 + i * 21, 108.5), 2.4, fc=COHORT_C[coh], ec="none", zorder=6))
    box(ax, 80.4, 104, 68, 44, fc="white"); T(ax, 114.4, 140, "CRISPRi screens", big, fontweight="bold")
    place(ax, "cas9", 98, 121, 28)
    for i, (c, ic) in enumerate([("K562", "myeloid_cell"), ("RPE1", "fibroblast_flat"), ("HepG2", "liver"), ("Jurkat", "lymphocyte")]):
        place(ax, ic, 121 + (i % 2) * 14, 126 - (i // 2) * 13, 11, tint=SCREEN_C[c])
    # row 2: score
    arrow(ax, 38, 101, 66, 86, lw=2.2, ms=18); arrow(ax, 114.4, 101, 86.4, 86, lw=2.2, ms=18)
    T(ax, 76.2, 79, "Reversal score", big, fontweight="bold")
    score_bar(ax, 36, 66, 80.4, 7, "", "", size=big)
    T(ax, 36, 61, "older", mid, ha="left", color=GRAY); T(ax, 116.4, 61, "younger", mid, ha="right", color=GRAY)
    for x in (30, 76.2, 122.4): arrow(ax, 76.2, 55, x, 47, lw=2.2, ms=18)
    # row 3: findings
    cards = [(4, "Cohort\ndecides\nranking"), (54.8, "Growth\ndrives\nscore"), (105.6, "Little left\nafter\nadjustment")]
    for x, lab in cards:
        box(ax, x, 3, 42.8, 42, fc="#f5f5f5"); T(ax, x + 21.4, 13.5, lab, mid, fontweight="bold", linespacing=1.0)
    # icon 1: two donor icons with crossed arrows
    place(ax, "old_man", 16, 35, 7); place(ax, "old_woman", 34, 35, 7.5)
    arrow(ax, 20, 38, 30, 29, col=INK, lw=1.6, ms=12); arrow(ax, 30, 38, 20, 29, col=INK, lw=1.6, ms=12)
    # icon 2: dividing cell -> slowed cell -> older end of the bar
    place(ax, "mitosis", 66, 34, 14); arrow(ax, 74.5, 34, 79, 34, col=INK, lw=1.6, ms=12)
    cc_wheel(ax, 84, 34, 4.2); arrow(ax, 89, 34, 92.5, 34, col=INK, lw=1.6, ms=12)
    score_bar(ax, 93.5, 32.5, 5, 3, "", "")
    # icon 3: near-flat bars
    for i, hgt in enumerate([1.6, 1.2, 1.9, 1.4, 1.7]):
        ax.add_patch(Rectangle((115 + i * 5, 29), 3.4, hgt, fc="#7f7f7f", ec="none", zorder=5))
    ax.plot([113, 141], [29, 29], color=INK, lw=1.2, zorder=6); ax.plot([113, 141], [38, 38], color=GRAY, lw=0.8, ls=":", zorder=6)
    fig.savefig(str(path) + ".pdf", facecolor="white"); fig.savefig(str(path) + ".png", dpi=300, facecolor="white")
    Image.open(str(path) + ".png").convert("RGB").save(str(path) + ".png", dpi=(300, 300))
    Image.open(str(path) + ".png").save(str(path) + ".tif", dpi=(300, 300), compression="tiff_lzw")
    Image.open(str(path) + ".png").resize((591, 591), Image.LANCZOS).save(str(path) + "_5cm_preview.png")
    plt.close(fig)
