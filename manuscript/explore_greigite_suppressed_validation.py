#!/usr/bin/env python3
"""Fig. 3 (B&W) — Engine validation: greigite suppressed (classical Fe-S) vs
greigite present, log f(S2) vs T over 300-600 K.

Boundary curves are read from the cached ENGINE output so the figure matches
Table 2 exactly (no analytic drift, no pycalphad re-run needed):
  - suppressed panel : engine/fes_greigite_suppressed_fields.json (committed)
  - present panel    : artifacts/fes_engine_boundaries.json["central"] (built by validate_fes_engine.py)
The 7 engine temperatures are smoothed with a cubic polynomial for display.
White background, black boundary curves, in-field labels.
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
OUT = str(FIG / "explore_greigite_suppressed_validation.png")


def series(dct, idx):
    """Return (T array, value array) from a {T:[lo,hi]} dict, picking column idx."""
    Ts = np.array(sorted(float(t) for t in dct))
    vals = np.array([dct[str(int(t))][idx] for t in Ts])
    return Ts, vals


def smooth(T, v, deg=3, n=400):
    c = np.polyfit(T, v, deg)
    Tf = np.linspace(T.min(), T.max(), n)
    return Tf, np.polyval(c, Tf)


sup = json.load(open(ENGINE / "fes_greigite_suppressed_fields.json"))
cen = json.load(open(ROOT / "artifacts" / "fes_engine_boundaries.json"))["central"][
    "fields"
]

# suppressed: Fe/Po (PYRRHOTITE lower) and Po/Py (PYRRHOTITE upper)
Ts, fe_po = series(sup["PYRRHOTITE"], 0)
_, po_py = series(sup["PYRRHOTITE"], 1)
# present: Fe/Gr (GREIGITE lower) and Gr/Py (GREIGITE upper)
_, fe_gr = series(cen["GREIGITE"], 0)
_, gr_py = series(cen["GREIGITE"], 1)

YLO, YHI = -56.0, -2.0
Tg = np.linspace(300, 600, 400)
Yg = np.linspace(YLO, YHI, 400)
TT, YY = np.meshgrid(Tg, Yg)
fig, axes = plt.subplots(1, 2, figsize=(13, 6.0), sharey=True)


def curve(pts):
    return np.polyval(np.polyfit(Ts, pts, 3), Tg)


def medoid(mask):
    """In-field grid point nearest the field centroid (axis-normalised)."""
    cx, cy = TT[mask].mean(), YY[mask].mean()
    d = ((TT[mask] - cx) / 300.0) ** 2 + ((YY[mask] - cy) / (YHI - YLO)) ** 2
    j = np.argmin(d)
    return TT[mask][j], YY[mask][j]


def panel(ax, lower, upper, mid_label, title):
    lo, hi = curve(lower), curve(upper)
    ax.plot(Tg, lo, color="black", lw=1.8)
    ax.plot(Tg, hi, color="black", lw=1.8)
    ax.set_xlim(300, 600)
    ax.set_ylim(YLO, YHI)
    ax.set_xlabel("Temperature, K")
    # classify stable phase on the grid: 0 = Fe, 1 = intermediate, 2 = pyrite
    region = np.where(YY < lo[None, :], 0, np.where(YY < hi[None, :], 1, 2))
    for idx, lab in [(2, "FeS$_2$\npyrite"), (1, mid_label), (0, "Fe")]:
        m = region == idx
        if m.sum() > 50:
            x, y = medoid(m)
            ax.text(x, y, lab, ha="center", va="center", fontsize=12, fontweight="bold")
    ax.set_title(title, fontsize=11.5, fontweight="bold")


panel(
    axes[0],
    fe_po,
    po_py,
    "Fe$_{1-x}$S\npyrrhotite",
    "greigite SUPPRESSED\n= classical Fe–S assessment (Dilner–Lee)",
)
panel(
    axes[1],
    fe_gr,
    gr_py,
    "Fe$_3$S$_4$\ngreigite",
    "greigite PRESENT (measured ΔH$_f$)\ngreigite is the stable Fe–S intermediate",
)
axes[0].set_ylabel("log $f$(S$_2$), bar")

fig.suptitle(
    "Fig. 3 — Engine validation: suppressing greigite recovers the classical pyrrhotite–pyrite Fe–S diagram;\n"
    "with greigite, greigite becomes the stable intermediate (sulfide/pyrite boundary shifts −31.8 → −20.3 at 300 K)",
    fontsize=12,
    fontweight="bold",
)
fig.text(
    0.5,
    0.005,
    "Boundaries are the cached engine values (match Table 2). Left reproduces the Dilner-Lee assessment; "
    "on heating, DSC/XRD show greigite, pyrrhotite and pyrite coexist over the decomposition window.",
    ha="center",
    va="bottom",
    fontsize=8,
    color=bw.FOOT_GREY,
)
fig.tight_layout(rect=[0, 0.03, 1, 0.92])
fig.savefig(OUT, dpi=200, bbox_inches="tight")
print("wrote", OUT)
