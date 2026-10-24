#!/usr/bin/env python3
"""Fig. 2 (B&W) — Greigite stability field and its error margin (engine).

Central greigite field (solid) bracketed by the two uncertainty bounds (dashed), which shift the
greigite enthalpy by its reported uncertainty, 21.9 kJ per Fe3S4 (Subramani et al. 2020, 2 SD of
the mean); pyrite and pyrrhotite stay at the assessed values. A pyrrhotite sliver (hatched)
appears only at the least-stable bound. Native-S saturation (dotted) is computed from the
database, and every boundary stops where it meets it. The central field is read from
artifacts/fes_engine_boundaries.json and the bounds from fig2b_fields.json
(explore_boundary_cases_2b.py probe lower/upper on the boundary databases). White background,
black lines, no colour.
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
from pycalphad import Database, calculate, variables as pv

HERE = Path(__file__).resolve().parent  # = manuscript/
ROOT = HERE.parent  # = repo root
FIG = ROOT / "artifacts" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE))
import bw_style as bw
import label_fit

bw.apply()
OUT = str(FIG / "Figure_2.png")


def series(dct, idx):
    Ts = np.array(sorted(float(t) for t in dct))
    return Ts, np.array([dct[str(int(t))][idx] for t in Ts])


cen = json.load(open(ROOT / "artifacts" / "fes_engine_boundaries.json"))["central"][
    "fields"
]
f2b = json.load(open(HERE / "fig2b_fields.json"))
DB = Database(str(ROOT / "artifacts" / "tdb" / "fes_greigite_v1.tdb"))


def _sym(name, t):
    return float(DB.symbols[name].subs({pv.T: float(t)}))


Ts, gc_lo = series(cen["GREIGITE"], 0)  # central Fe/Gr
_, gl_lo = series(f2b["lower"]["GREIGITE"], 0)  # most-stable bound (wider) Fe/Gr
_, gu_lo = series(f2b["upper"]["GREIGITE"], 0)  # least-stable bound (narrower) Gr lower (=Po/Gr)
_, po_lo = series(f2b["upper"]["PYRRHOTITE"], 0)  # least-stable bound Fe/Po
_, po_hi = series(f2b["upper"]["PYRRHOTITE"], 1)  # least-stable bound Po/Gr (= gu_lo)

YLO, YHI = -56.0, -2.0
Tg = np.linspace(300, 600, 400)
Yg = np.linspace(YLO, YHI, 400)
TT, YY = np.meshgrid(Tg, Yg)

# Native-S saturation, log f(S2) = (2 G(S, orthorhombic) - G(S2, gas)) / (RT ln 10), depends on
# sulfur alone, so it comes straight from the database. Above it native sulfur is stable.
sat_c = np.array([(2.0 * _sym("GHSERSS", t) - _sym("F15281T", t)) / (8.31451 * t * np.log(10.0))
                  for t in Tg])


def _gr_py(tdb):
    """Greigite/pyrite boundary from the Gibbs energies, cut at native-S saturation.

    Both are line compounds, so per mole Fe G(Fe3S4)/3 + mu(S2)/3 = G(FeS2) is exact. A fit
    through the probe points cannot place the end at saturation (it sits between grid points)."""
    db = Database(str(ROOT / "artifacts" / "tdb" / tdb))
    g = {ph: np.ravel(calculate(db, comps, ph, T=Tg, P=101325, N=1, output="GM").GM.values)
         for ph, comps in (("PYRITE", ["FE", "S"]), ("GREIGITE", ["FE", "S"]))}
    mu_s2 = 3.0 * (3.0 * g["PYRITE"] - 7.0 * g["GREIGITE"] / 3.0)
    s2_gas = np.array([_sym("F15281T", t) for t in Tg])
    c = (mu_s2 - s2_gas) / (8.31451 * Tg * np.log(10.0))
    above = np.nonzero(c >= sat_c)[0]
    if above.size:
        c[above[0]:] = np.nan
    return c


def _curve(pts):
    """Smooth boundary on Tg through the probe points below saturation, cut where it meets it.

    Used for the boundaries that involve Fe or non-stoichiometric pyrrhotite, which never reach
    saturation here."""
    keep = pts < np.interp(Ts, Tg, sat_c) - 0.02
    c = np.polyval(np.polyfit(Ts[keep], pts[keep], min(3, keep.sum() - 1)), Tg)
    above = np.nonzero(c >= sat_c)[0]
    if above.size:
        c[above[0]:] = np.nan
    return c


fig, ax = plt.subplots(figsize=(9.5, 6.4))

gp_c = _gr_py("fes_greigite_v1.tdb")              # central Gr/Py
gp_l = _gr_py("fes_greigite_boundary_lower.tdb")  # most-stable bound Gr/Py
gp_u = _gr_py("fes_greigite_boundary_upper.tdb")  # least-stable bound Gr/Py

# central greigite field (solid)
for c in (_curve(gc_lo), gp_c):
    ax.plot(Tg, c, color="black", lw=2.0)
# reported-uncertainty envelope (dashed)
for c in (_curve(gl_lo), gp_l, _curve(gu_lo), gp_u):
    ax.plot(Tg, c, color="black", lw=1.2, ls=(0, (5, 4)))
ax.plot(Tg, sat_c, color="black", lw=1.0, ls=(0, (1, 2)))

# pyrrhotite sliver that appears only at the least-stable bound (hatched band Fe/Po -> Po/Gr)
plo_c, phi_c = _curve(po_lo), _curve(po_hi)
ax.fill_between(Tg, plo_c, phi_c, **bw.hatch_kw("////"))

ax.set_xlim(300, 600)
ax.set_ylim(YLO, YHI)
ax.set_xlabel("Temperature, K")
ax.set_ylabel("log $f$(S$_2$), bar")

# fields of the central case: 0 = Fe, 1 = greigite, 2 = pyrite, 3 = native S
lo_c = _curve(gc_lo)
hi_c = np.fmin(np.nan_to_num(gp_c, nan=np.inf), sat_c)
region = np.where(YY < lo_c[None, :], 0, np.where(YY < hi_c[None, :], 1, 2))
region = np.where(YY > sat_c[None, :], 3, region)

m = region == 1
cx, cy = TT[m].mean(), YY[m].mean()
j = np.argmin(((TT[m] - cx) / 300.0) ** 2 + ((YY[m] - cy) / (YHI - YLO)) ** 2)
ax.text(TT[m][j], YY[m][j], "Fe$_3$S$_4$\ngreigite", ha="center", va="center", fontsize=13,
        fontweight="bold")

# The pyrite band is crossed by the most-stable dashed curve, so its label is a callout from the
# native-S wedge to the band centre at 315 K.
k = int(np.argmin(np.abs(Tg - 315.0)))
ax.annotate(
    "FeS$_2$ pyrite",
    xy=(Tg[k], 0.5 * (hi_c[k] + gp_l[k])),
    xytext=(330, -5.0),
    ha="left",
    va="center",
    fontsize=12,
    fontweight="bold",
    arrowprops=dict(arrowstyle="->", color="black", lw=1.0, shrinkA=2, shrinkB=2),
    bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.5", lw=0.6),
)
# Fe1-xS sliver (least-stable bound only): label placed in clear space below
# the band, with an arrow pointing directly at the centre of the hatched sliver.
k = int(np.argmin(np.abs(Tg - 345.0)))
ax.annotate(
    "Fe$_{1-x}$S field\n(least-stable bound only)",
    # position checked clear of every curve, the band, the legend and the axes edge
    xy=(345.0, 0.5 * (plo_c[k] + phi_c[k])),
    xytext=(345, -52.0),
    ha="left",
    va="center",
    fontsize=12,
    style="italic",
    fontweight="bold",
    arrowprops=dict(arrowstyle="->", color="black", lw=1.0, shrinkA=2, shrinkB=2),
    bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.5", lw=0.6),
)

leg = [
    Line2D([0], [0], color="black", lw=2.0, label="greigite field (central ΔH$_f$)"),
    Line2D(
        [0], [0], color="black", lw=1.2, ls=(0, (5, 4)), label="greigite ± reported uncertainty (2 SD)"
    ),
    Patch(
        facecolor="none",
        edgecolor="black",
        hatch="////",
        label="Fe$_{1-x}$S only at the least-stable bound",
    ),
    Line2D([0], [0], color="black", lw=1.0, ls=(0, (1, 2)), label="native-S saturation"),
]
ax.legend(handles=leg, loc="lower right", frameon=True, edgecolor="black", fontsize=13)

fig.tight_layout(rect=[0, 0.03, 1, 1])
# Fe and native-S labels last, once the legend and layout exist, each checked clear
label_fit.place(ax, "Fe", region, Tg, Yg, 0, fontweight="bold")
label_fit.place(ax, "native S", region, Tg, Yg, 3, fontweight="bold")
fig.savefig(OUT, dpi=200, bbox_inches="tight")
print("wrote", OUT)
