"""Reproduce Waldner & Pelton (2005) Fig 9: log f(S2) vs X_S for single-phase
pyrrhotite at 700 C and 1100 C.  Validation gate for the WP model."""

import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TDB = ROOT / "artifacts" / "tdb"
FIG = ROOT / "artifacts" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE))

import pyrrhotite_wp as wp

# Sweep y_Fe so X_S = 1/(1+y_Fe) runs 0.50 -> 0.55.
# X_S=0.50 -> y_Fe=1 (singular, ln y_Va), so start just below.
XS_min, XS_max = 0.5001, 0.550
XS = np.linspace(XS_min, XS_max, 400)
yFe = wp.yFe_of_X_S(XS)
yFe = np.clip(yFe, 1e-6, 1.0 - 1e-9)

temps = {"700 C (973.15 K)": 973.15, "1100 C (1373.15 K)": 1373.15}
colors = {"700 C (973.15 K)": "C0", "1100 C (1373.15 K)": "C3"}

fig, ax = plt.subplots(figsize=(7, 5.5))
table_rows = []
for label, T in temps.items():
    lf = np.array([wp.log_fS2(y, T) for y in yFe])
    ax.plot(XS, lf, color=colors[label], lw=2, label=label)
    # sample points for the printed table
    for xs_t in (0.500, 0.505, 0.510, 0.520, 0.530, 0.540, 0.545):
        yt = (1.0 - xs_t) / xs_t
        yt = min(max(yt, 1e-6), 1.0 - 1e-9)
        table_rows.append((label, xs_t, wp.log_fS2(yt, T)))

ax.set_xlabel("mole fraction S,  X$_S$")
ax.set_ylabel("log$_{10}$ f(S$_2$)   [bar]")
ax.set_title(
    "Waldner & Pelton (2005) Fig 9 reproduction\nsingle-phase pyrrhotite Fe$_{1-x}$S"
)
ax.set_xlim(0.50, 0.55)
ax.set_ylim(-20, 6)
ax.axhline(0, color="gray", lw=0.5)
ax.grid(alpha=0.3)
ax.legend(loc="lower right")
fig.tight_layout()
fig.savefig(str(FIG / "pyrrhotite_wp_validation_fig9.png"), dpi=130)
print("saved pyrrhotite_wp_validation_fig9.png")

print("\n  T-label              X_S     log f(S2)")
print("  " + "-" * 44)
for lab, xs, lf in table_rows:
    print(f"  {lab:18s}  {xs:.3f}   {lf:8.3f}")
