#!/usr/bin/env python3
"""Fig. 11 (B&W) — Powder XRD of bulk Fe3S4 after DSC to 600 C.

The measured pattern (data_dsc/, K. Lilova) shown verbatim; the product is
pyrrhotite-3T (Fe7S8). Single panel, no footnote (text in caption).
"""

import sys
from pathlib import Path
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
OUT = str(FIG / "fig11_xrd.png")

img = plt.imread(str(DATA / "Fe3S4_bulk_XRD_postDSC600.png"))
fig, ax = plt.subplots(figsize=(8.0, 6.0))
ax.imshow(img)
ax.axis("off")
ax.set_title(
    "Fig. 11 — PXRD of bulk Fe$_3$S$_4$ after DSC to 600 °C → pyrrhotite-3T (Fe$_7$S$_8$)",
    fontsize=12,
    fontweight="bold",
)
fig.tight_layout()
fig.savefig(OUT, dpi=200, bbox_inches="tight")
print("wrote", OUT)
