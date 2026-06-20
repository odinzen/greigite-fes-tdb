#!/usr/bin/env python3
"""Fig. 2 (B&W) — Greigite heat capacity, Debye-Einstein fit (Shumway 2022).

Reconstructs the manuscript Cp figure in black and white from the published
three-term Debye-Einstein fit in greigite_cp.py (per FeS1.33). Single black
curve, black calorimetric marker at 298.15 K, no colour.
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
sys.path.insert(0, str(HERE))
import bw_style as bw
import greigite_cp as g

bw.apply()
OUT = str(FIG / "fig02_cp.png")

T = np.linspace(5.0, 600.0, 600)
cp = g.cp_debye_einstein(T)

fig, ax = plt.subplots(figsize=(7.2, 5.2))
ax.plot(T, cp, color="black", lw=1.8, label="Debye–Einstein fit (Shumway 2022)")
ax.plot(
    298.15,
    g.CP_298,
    "o",
    color="black",
    ms=7,
    label="calorimetric C$_p$(298.15) = %.2f J/mol·K" % g.CP_298,
)

ax.set_xlabel("Temperature, K")
ax.set_ylabel("C$_p$ per FeS$_{1.33}$, J/mol·K")
ax.set_xlim(0, 600)
ax.set_ylim(0, 80)
ax.legend(loc="lower right", frameon=True, edgecolor="black", fontsize=10)

fig.tight_layout(rect=[0, 0.04, 1, 1])
fig.savefig(OUT, dpi=200, bbox_inches="tight")
print("wrote", OUT)
