#!/usr/bin/env python3
"""Fig. S6 (2026 re-measurement, B&W) — Simultaneous TG-DSC of bulk Fe3S4.

The new data taken six years after the 2020 DSC: blank-subtracted simultaneous
TG-DSC (Setaram Labsys EVO, Ar), read from the v2 datasets in data_dsc/. Left
axis: DSC heat flow (mW, exo up), heating solid / cooling dashed. Right axis: TG
mass (%, as exported) — the new information the 2026 run adds. Kelvin x-axis to
match the paper; no in-image title (the caption lives in the manuscript).

Companion to make_Figure_S4_dsc_2020.py (v1, 2020 DSC); see data_dsc/SOURCE.md.
"""

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
DATA = HERE / "data_dsc"
sys.path.insert(0, str(HERE))
import bw_style as bw

bw.apply()
OUT = str(FIG / "Figure_S6.png")


def load(fn, ncol):
    """Parse a v2 SETARAM .txt (UTF-16, space-padded). Returns an (N, ncol) array
    of the numeric columns after the 'Sample Temperature' header row."""
    lines = open(fn, "rb").read().decode("utf-16").splitlines()
    h = next(i for i, l in enumerate(lines) if l.startswith("Sample Temperature"))
    rows = []
    for l in lines[h + 1 :]:
        p = l.split()
        if len(p) < ncol:
            continue
        try:
            rows.append([float(x) for x in p[:ncol]])
        except ValueError:
            continue
    return np.array(rows)


# heating: T, TG%, HeatFlow(mW) ; cooling: T, HeatFlow(mW) (DSC only, no TG)
H = load(str(DATA / "Fe3S4_bulk_heating_v2.txt"), 3)
C = load(str(DATA / "Fe3S4_bulk_cooling_v2.txt"), 2)
Th, TGh, HFh = H[:, 0] + 273.15, H[:, 1], H[:, 2]
Tc, HFc = C[:, 0] + 273.15, C[:, 1]

fig, axL = plt.subplots(figsize=(8.2, 5.6))
axR = axL.twinx()

# DSC (left axis)
(lh,) = axL.plot(Th, HFh, color="black", lw=1.6, label="DSC heating")
(lc,) = axL.plot(Tc, HFc, color="black", lw=1.2, ls=(0, (5, 4)), label="DSC cooling")
axL.set_xlabel("Temperature, K")
axL.set_ylabel("DSC heat flow, mW  (exo up)")
axL.set_xlim(min(Th.min(), Tc.min()), max(Th.max(), Tc.max()))

# TG (right axis)
(lt,) = axR.plot(Th, TGh, color="0.45", lw=2.0, ls=(0, (1, 1)), label="TG (mass)")
axR.set_ylabel("TG mass, %", color="0.3")
axR.tick_params(axis="y", colors="0.3")

axL.text(
    0.03,
    0.06,
    "Fe$_3$S$_4$ — bulk (2026, blank-subtracted)",
    transform=axL.transAxes,
    fontsize=10.5,
    fontweight="bold",
    va="bottom",
    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=1.0),
)
axL.legend(
    handles=[lh, lc, lt],
    loc="upper right",
    frameon=True,
    edgecolor="black",
    fontsize=9.5,
)
fig.tight_layout()
fig.savefig(OUT, dpi=200, bbox_inches="tight")
print("wrote", OUT)
print(
    "heating %.0f-%.0f K | TG %.1f -> %.1f %% | endo trough %.0f K"
    % (Th.min(), Th.max(), TGh[0], TGh[-1], Th[int(np.argmin(HFh))])
)
print(
    "cooling %.0f-%.0f K | DSC exo peak %.0f K"
    % (Tc.min(), Tc.max(), Tc[int(np.argmax(HFc))])
)
