"""Place a field label where it fits.

A label goes at the point deepest inside its field. If it does not fit there (it would cross a
drawn line, leave its field, touch another label or the legend, or leave the axes), smaller,
single-line and finally rotated versions are tried at the deepest few points. Fit is measured on
the rendered figure, so call this after the axes limits and layout are final.
"""
import numpy as np
from scipy.ndimage import distance_transform_edt
from matplotlib.collections import PolyCollection

# (single line?, font size, rotation) in order of preference
VARIANTS = [(False, 11.5, 0), (True, 11.5, 0), (False, 9.5, 0), (True, 9.5, 0),
            (True, 8.0, 0), (True, 9.5, 90), (True, 8.0, 90)]


def _occupied(ax, r):
    boxes = [t.get_window_extent(r) for t in ax.texts if t.get_text().strip()]
    if ax.get_legend() is not None:
        boxes.append(ax.get_legend().get_window_extent(r))
    return boxes


def _fits(ax, r, bb, field, X, Y, idx, taken):
    axbb = ax.get_window_extent(r)
    if not (axbb.contains(bb.x0, bb.y0) and axbb.contains(bb.x1, bb.y1)):
        return False
    if any(bb.overlaps(o) for o in taken):
        return False
    pts = ax.transData.inverted().transform(
        [(x, y) for x in np.linspace(bb.x0, bb.x1, 9) for y in np.linspace(bb.y0, bb.y1, 7)])
    ix = np.clip(np.rint((pts[:, 0] - X[0]) / (X[1] - X[0])).astype(int), 0, len(X) - 1)
    iy = np.clip(np.rint((pts[:, 1] - Y[0]) / (Y[1] - Y[0])).astype(int), 0, len(Y) - 1)
    if not np.all(field[iy, ix] == idx):
        return False
    for ln in ax.get_lines():
        if ln.get_path().transformed(ln.get_transform()).intersects_bbox(bb, filled=False):
            return False
    for col in ax.collections:
        filled = isinstance(col, PolyCollection)
        if any(p.transformed(col.get_transform()).intersects_bbox(bb, filled=filled) for p in col.get_paths()):
            return False
    return True


def place(ax, text, field, X, Y, idx, n_points=12, **style):
    """Label field index `idx` of `field` (shape len(Y) x len(X), regular grid). Returns the Text."""
    fig = ax.figure
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    p0, px = ax.transData.transform([(X[0], Y[0]), (X[1], Y[1])])
    # pad with a border of "outside" so the axes edge counts as a boundary; otherwise the deepest
    # point of a field that touches the edge sits on the edge itself
    mask = np.pad(field == idx, 1, constant_values=False)
    depth = distance_transform_edt(mask, sampling=(abs(px[1] - p0[1]), abs(px[0] - p0[0])))[1:-1, 1:-1]
    order = np.argsort(depth, axis=None)[::-1]
    picks = []                                   # deepest points, kept apart from each other
    for k in order[:20000]:
        j, i = np.unravel_index(k, depth.shape)
        if depth[j, i] <= 0:
            break
        if all((i - a) ** 2 + (j - b) ** 2 > 400 for a, b in picks):
            picks.append((i, j))
        if len(picks) == n_points:
            break
    taken = _occupied(ax, r)
    for single, size, rot in VARIANTS:
        label = text.replace("\n", " ") if single else text
        for i, j in picks:
            t = ax.text(X[i], Y[j], label, ha="center", va="center", fontsize=size, rotation=rot, **style)
            bb = t.get_window_extent(r).expanded(1.04, 1.10)
            if _fits(ax, r, bb, field, X, Y, idx, taken):
                return t
            t.remove()
    i, j = picks[0]
    print(f"label_fit: no clean fit for {text!r}; placed at the deepest point")
    return ax.text(X[i], Y[j], text.replace("\n", " "), ha="center", va="center", fontsize=8, **style)
