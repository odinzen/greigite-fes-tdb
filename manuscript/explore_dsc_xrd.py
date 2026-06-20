#!/usr/bin/env python3
"""Fig. 8 (B&W) — Experimental decomposition of bulk greigite (DSC + powder XRD).

(a) DSC of bulk Fe3S4 re-plotted from the raw instrument export (heating to
    600 C, solid; cooling, dashed) in data_dsc/.
(b) Powder XRD of bulk Fe3S4 after DSC to 600 C — product is pyrrhotite-3T
    (Fe7S8); the measured pattern (data_dsc/) is shown verbatim.

All bulk, all from primary data (K. Lilova, co-author). Black-and-white.
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
OUT = str(FIG / "explore_dsc_xrd.png")


def parse_dsc(path):
    txt = open(path, "rb").read().decode("utf-16").splitlines()
    h = next(i for i, l in enumerate(txt) if l.startswith("Index"))
    T, HF = [], []
    for l in txt[h + 1 :]:
        p = l.split()
        if len(p) < 5:
            continue
        try:
            T.append(float(p[2]))
            HF.append(float(p[4]))
        except ValueError:
            continue
    return np.array(T), np.array(HF)


def smooth(y, w=31):
    k = np.ones(w) / w
    return np.convolve(y, k, mode="same")


Th, Hh = parse_dsc(str(DATA / "Fe3S4_bulk_heating.txt"))
Tc, Hc = parse_dsc(str(DATA / "Fe3S4_bulk_cooling.txt"))


# trim convolution edge artefacts
def trim(T, H, w=31):
    return T[w:-w], smooth(H, w)[w:-w]


Th, Hh = trim(Th, Hh)
Tc, Hc = trim(Tc, Hc)

fig, (axA, axB) = plt.subplots(
    1, 2, figsize=(13, 5.4), gridspec_kw={"width_ratios": [1.0, 1.05]}
)

# ---- panel (a): DSC ----
axA.plot(Th, Hh, color="black", lw=1.5, label="heating")
axA.plot(Tc, Hc, color="black", lw=1.2, ls=(0, (5, 4)), label="cooling")
axA.axvline(317, color="0.55", lw=0.9, ls=":")
axA.annotate(
    "pyritization\nonset ~317 °C",
    xy=(317, 1.2),
    xytext=(150, 2.0),
    fontsize=8.5,
    ha="left",
    arrowprops=dict(arrowstyle="->", color="0.4", lw=0.8),
)
axA.annotate(
    "decomposition\nendotherm ~540 °C",
    xy=(540, -4.6),
    xytext=(360, -3.6),
    fontsize=8.5,
    ha="left",
    arrowprops=dict(arrowstyle="->", color="0.4", lw=0.8),
)
axA.text(
    0.03,
    0.95,
    "Fe$_3$S$_4$ — bulk",
    transform=axA.transAxes,
    fontsize=11,
    fontweight="bold",
    va="top",
    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=1.0),
)
axA.set_xlabel("Temperature, °C")
axA.set_ylabel("Heat flow, µV  (exo up)")
axA.set_xlim(20, 600)
axA.legend(loc="lower left", frameon=True, edgecolor="black", fontsize=9.5)
axA.set_title(
    "(a) DSC of bulk Fe$_3$S$_4$ (heating + cooling)", fontsize=11, fontweight="bold"
)

# ---- panel (b): XRD image ----
img = plt.imread(str(DATA / "Fe3S4_bulk_XRD_postDSC600.png"))
axB.imshow(img)
axB.axis("off")
axB.set_title(
    "(b) PXRD of bulk Fe$_3$S$_4$ after DSC to 600 °C → pyrrhotite-3T (Fe$_7$S$_8$)",
    fontsize=11,
    fontweight="bold",
)

fig.suptitle(
    "Fig. 8 — Experimental decomposition of bulk greigite (DSC + powder XRD)",
    fontsize=12.5,
    fontweight="bold",
)
fig.text(
    0.5,
    0.005,
    "Greigite is stable to ~300 °C; surface-oxidation/pyritization exotherms 200–450 °C; "
    "decomposition endotherm ~540 °C; the 600 °C product is pyrrhotite-3T (Fe$_7$S$_8$). "
    "Data: K. Lilova (co-author).",
    ha="center",
    va="bottom",
    fontsize=8,
    color=bw.FOOT_GREY,
)
fig.tight_layout(rect=[0, 0.04, 1, 0.94])
fig.savefig(OUT, dpi=200, bbox_inches="tight")
print("wrote", OUT)
print(
    "heating T range %.0f-%.0f C ; cooling %.0f-%.0f C"
    % (Th.min(), Th.max(), Tc.min(), Tc.max())
)
