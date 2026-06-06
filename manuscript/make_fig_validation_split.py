#!/usr/bin/env python3
"""Figs. 3 & 4 (B&W, single-panel) — Engine validation, split into two figures.

Fig. 3: greigite SUPPRESSED (classical Fe-S, Dilner-Lee).
Fig. 4: greigite PRESENT (measured dHf) -> greigite is the stable intermediate.

Boundary curves from the cached engine output (match Table 2). Labels placed at
field medoids (no label crosses a boundary). No footnote on the figure; the
explanatory text lives in the manuscript caption. Black & white.
"""

import json
import sys
from pathlib import Path
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent  # = manuscript/
ROOT = HERE.parent  # = repo root
FIG = ROOT / "artifacts" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
ENGINE = ROOT / "engine"
sys.path.insert(0, str(HERE))
import bw_style as bw

bw.apply()


def series(dct, idx):
    Ts = np.array(sorted(float(t) for t in dct))
    return Ts, np.array([dct[str(int(t))][idx] for t in Ts])


sup = json.load(open(ENGINE / "fes_greigite_suppressed_fields.json"))
cen = json.load(open(ROOT / "artifacts" / "fes_engine_boundaries.json"))["central"][
    "fields"
]
Ts, fe_po = series(sup["PYRRHOTITE"], 0)
_, po_py = series(sup["PYRRHOTITE"], 1)
_, fe_gr = series(cen["GREIGITE"], 0)
_, gr_py = series(cen["GREIGITE"], 1)

YLO, YHI = -56.0, -2.0
Tg = np.linspace(300, 600, 400)
Yg = np.linspace(YLO, YHI, 400)
TT, YY = np.meshgrid(Tg, Yg)


def curve(pts):
    return np.polyval(np.polyfit(Ts, pts, 3), Tg)


def medoid(mask):
    cx, cy = TT[mask].mean(), YY[mask].mean()
    d = ((TT[mask] - cx) / 300.0) ** 2 + ((YY[mask] - cy) / (YHI - YLO)) ** 2
    j = np.argmin(d)
    return TT[mask][j], YY[mask][j]


def make(lower, upper, mid_label, title, out):
    lo, hi = curve(lower), curve(upper)
    fig, ax = plt.subplots(figsize=(7.0, 5.6))
    ax.plot(Tg, lo, color="black", lw=1.8)
    ax.plot(Tg, hi, color="black", lw=1.8)
    ax.set_xlim(300, 600)
    ax.set_ylim(YLO, YHI)
    ax.set_xlabel("Temperature, K")
    ax.set_ylabel("log $f$(S$_2$), bar")
    region = np.where(YY < lo[None, :], 0, np.where(YY < hi[None, :], 1, 2))
    for idx, lab in [(2, "FeS$_2$\npyrite"), (1, mid_label), (0, "Fe")]:
        m = region == idx
        if m.sum() > 50:
            x, y = medoid(m)
            ax.text(x, y, lab, ha="center", va="center", fontsize=13, fontweight="bold")
    ax.set_title(title, fontsize=12.5, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print("wrote", out)


make(
    fe_po,
    po_py,
    "Fe$_{1-x}$S\npyrrhotite",
    "Fig. 3 — greigite suppressed (classical Fe–S, Dilner–Lee)",
    str(FIG / "fig03_validation_suppressed.png"),
)
make(
    fe_gr,
    gr_py,
    "Fe$_3$S$_4$\ngreigite",
    "Fig. 4 — greigite present: greigite is the stable Fe–S intermediate",
    str(FIG / "fig04_validation_present.png"),
)
