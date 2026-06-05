#!/usr/bin/env python3
"""Fig. 5 (B&W) — Greigite stability field and its ±1σ error margin (engine).

Central greigite field (solid) bracketed by ±1σ enthalpy bounds (dashed); a
pyrrhotite sliver (hatched) appears only at the +1σ (least-stable) bound.
All boundaries are cached engine values, so the figure matches Table 2 / the
text margin (+19.5 ± 7.9 kJ, +2.5σ). White background, black lines, no colour.
"""

import json
import sys
from pathlib import Path
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent  # = manuscript/
ROOT = HERE.parent  # = repo root
FIG = ROOT / "artifacts" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
ENGINE = ROOT / "engine"
sys.path.insert(0, str(HERE))
import bw_style as bw

bw.apply()
OUT = str(FIG / "fig06_error_field.png")


def series(dct, idx):
    Ts = np.array(sorted(float(t) for t in dct))
    return Ts, np.array([dct[str(int(t))][idx] for t in Ts])


def smooth(T, v, deg=3, n=400):
    c = np.polyfit(T, v, deg)
    Tf = np.linspace(T.min(), T.max(), n)
    return Tf, np.polyval(c, Tf)


cen = json.load(open(ROOT / "artifacts" / "fes_engine_boundaries.json"))["central"][
    "fields"
]
f2b = json.load(open(HERE / "fig2b_fields.json"))

Ts, gc_lo = series(cen["GREIGITE"], 0)  # central Fe/Gr
_, gc_hi = series(cen["GREIGITE"], 1)  # central Gr/Py
_, gl_lo = series(f2b["lower"]["GREIGITE"], 0)  # -1sigma (wider) Fe/Gr
_, gl_hi = series(f2b["lower"]["GREIGITE"], 1)  # -1sigma Gr/Py
_, gu_lo = series(f2b["upper"]["GREIGITE"], 0)  # +1sigma (narrower) Gr lower (=Po/Gr)
_, gu_hi = series(f2b["upper"]["GREIGITE"], 1)  # +1sigma Gr/Py
_, po_lo = series(f2b["upper"]["PYRRHOTITE"], 0)  # +1sigma Fe/Po
_, po_hi = series(f2b["upper"]["PYRRHOTITE"], 1)  # +1sigma Po/Gr (= gu_lo)

fig, ax = plt.subplots(figsize=(9.5, 6.4))

# central greigite field (solid)
for v in (gc_lo, gc_hi):
    Tf, s = smooth(Ts, v)
    ax.plot(Tf, s, color="black", lw=2.0)
# +-1 sigma envelope (dashed)
for v in (gl_lo, gl_hi, gu_lo, gu_hi):
    Tf, s = smooth(Ts, v)
    ax.plot(Tf, s, color="black", lw=1.2, ls=(0, (5, 4)))

# pyrrhotite sliver that appears only at +1 sigma (hatched band Fe/Po -> Po/Gr)
Tf, plo = smooth(Ts, po_lo)
_, phi = smooth(Ts, po_hi)
ax.fill_between(Tf, plo, phi, **bw.hatch_kw("////"))

ax.set_xlim(300, 600)
ax.set_ylim(-56, -2)
ax.set_xlabel("Temperature, K")
ax.set_ylabel("log $f$(S$_2$), bar")

# --- physically-correct labels: medoid of each field on a grid ---
YLO, YHI = -56.0, -2.0
Tg = np.linspace(300, 600, 400)
Yg = np.linspace(YLO, YHI, 400)
TT, YY = np.meshgrid(Tg, Yg)


def _curve(pts):
    return np.polyval(np.polyfit(Ts, pts, 3), Tg)


def _medoid(mask):
    cx, cy = TT[mask].mean(), YY[mask].mean()
    dd = ((TT[mask] - cx) / 300.0) ** 2 + ((YY[mask] - cy) / (YHI - YLO)) ** 2
    j = np.argmin(dd)
    return TT[mask][j], YY[mask][j]


lo_c, hi_c = _curve(gc_lo), _curve(gc_hi)  # central greigite field
region = np.where(YY < lo_c[None, :], 0, np.where(YY < hi_c[None, :], 1, 2))
for idx, lab, fs in [
    (2, "FeS$_2$ pyrite", 12),
    (1, "Fe$_3$S$_4$\ngreigite", 13),
    (0, "Fe", 12),
]:
    m = region == idx
    if m.sum() > 50:
        x, y = _medoid(m)
        ax.text(x, y, lab, ha="center", va="center", fontsize=fs, fontweight="bold")
# Fe1-xS sliver (appears only at +1 sigma): label placed in clear space below
# the band, with an arrow pointing directly at the centre of the hatched sliver.
plo_c, phi_c = _curve(po_lo), _curve(po_hi)
ti = 345.0
k = int(np.argmin(np.abs(Tg - ti)))
ytip = 0.5 * (plo_c[k] + phi_c[k])
ax.annotate(
    "Fe$_{1-x}$S field\n(appears only at +1σ)",
    xy=(ti, ytip),
    xytext=(311, -54.3),
    ha="left",
    va="center",
    fontsize=9,
    style="italic",
    arrowprops=dict(arrowstyle="->", color="black", lw=1.0, shrinkA=2, shrinkB=2),
    bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.5", lw=0.6),
)

ax.set_title(
    "Fig. 6 — Greigite stability field and its ±1σ error margin (engine)\n"
    "pyrrhotite (Fe$_{1-x}$S) enters only if greigite sits at its +1σ (least-stable) bound",
    fontsize=11.5,
    fontweight="bold",
)

leg = [
    Line2D([0], [0], color="black", lw=2.0, label="greigite field (central ΔH$_f$)"),
    Line2D(
        [0], [0], color="black", lw=1.2, ls=(0, (5, 4)), label="greigite ±1σ margin"
    ),
    Patch(
        facecolor="none",
        edgecolor="black",
        hatch="////",
        label="Fe$_{1-x}$S appears only at greigite +1σ",
    ),
]
ax.legend(handles=leg, loc="lower right", frameon=True, edgecolor="black", fontsize=9.5)

fig.tight_layout(rect=[0, 0.03, 1, 1])
fig.savefig(OUT, dpi=200, bbox_inches="tight")
print("wrote", OUT)
