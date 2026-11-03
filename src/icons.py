"""Licensed icons (figures/icons, see CREDITS.md) rasterized for matplotlib: icon(name, px) -> RGBA array; place(ax, name, x, y, w)
draws the icon centred at (x, y) with width w in data coordinates, keeping its aspect ratio."""
import io, os, re
from functools import lru_cache
import numpy as np
from PIL import Image
import resvg_py

ICON_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures", "icons")


def _svg_text(name):
    s = open(os.path.join(ICON_DIR, name if name.endswith(".svg") else name + ".svg"), encoding="utf-8", errors="ignore").read()
    head = re.search(r"<svg\b[^>]*>", s, flags=re.S).group(0)
    vb = re.search(r'viewBox\s*=\s*"([^"]+)"', head)
    if vb:  # pad the viewBox by 25% on every side so artwork that runs past it is never clipped; icon() then crops to the artwork
        x0, y0, w, h = [float(v) for v in re.split(r"[ ,]+", vb.group(1).strip())]
        px_, py_ = 0.25 * w, 0.25 * h; x0, y0, w, h = x0 - px_, y0 - py_, w + 2 * px_, h + 2 * py_
        new = re.sub(r'\s(width|height|viewBox)\s*=\s*"[^"]*"', "", head)
        new = new[:-1] + f' viewBox="{x0} {y0} {w} {h}" width="{w}" height="{h}">'
        s = s.replace(head, new, 1)
    return s


KNOCKOUT_WHITE = {"network"}   # icons drawn on a white rectangle: make that background transparent


@lru_cache(maxsize=64)
def icon(name, px=400, tint=None):
    png = bytes(resvg_py.svg_to_bytes(svg_string=_svg_text(name), width=int(px * 1.5)))
    im = np.asarray(Image.open(io.BytesIO(png)).convert("RGBA")).astype(float) / 255
    if name in KNOCKOUT_WHITE:
        im = im.copy(); im[..., 3] = np.where(im[..., :3].min(axis=2) > 0.97, 0, im[..., 3])
    ys, xs = np.nonzero(im[..., 3] > 0.02)   # crop to the artwork plus a 3% margin
    m = int(0.03 * max(np.ptp(ys), np.ptp(xs))) + 1
    im = im[max(ys.min() - m, 0):ys.max() + m + 1, max(xs.min() - m, 0):xs.max() + m + 1]
    if tint is not None:  # recolor all non-transparent pixels to one color, keeping shading through luminance
        from matplotlib.colors import to_rgb
        lum = im[..., :3].mean(axis=2, keepdims=True); c = np.array(to_rgb(tint))[None, None, :]
        im = np.concatenate([1 - (1 - c) * (1 - 0.55 * lum), im[..., 3:]], axis=2)
    return im


def place(ax, name, x, y, w, px=500, tint=None, zorder=5, transform=None):
    im = icon(name, px, tint); h = w * im.shape[0] / im.shape[1]
    if transform is None:
        # keep the icon's aspect in data units on an axes with equal aspect
        ax.imshow(im, extent=(x - w / 2, x + w / 2, y - h / 2, y + h / 2), zorder=zorder, interpolation="lanczos")
    return h
