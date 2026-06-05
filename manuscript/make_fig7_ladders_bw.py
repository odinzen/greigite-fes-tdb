#!/usr/bin/env python3
"""Fig. 7 (B&W) — Greigite is the sulfur analogue of magnetite: parallel Fe-O
and Fe-S redox ladders at 298 K (schematic).

Pure schematic (no engine call): the Fe(II,III) spinel Fe3X4 is the stable
intermediate in BOTH systems — magnetite (Fe3O4) in Fe-O, greigite (Fe3S4) in
Fe-S — and both are ferrimagnetic. Boundary fugacities are the engine values
from Table 3. White background, black boxes and text, no colour.
"""

import sys
from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

HERE = Path(__file__).resolve().parent  # = manuscript/
ROOT = HERE.parent  # = repo root
FIG = ROOT / "artifacts" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE))
import bw_style as bw

bw.apply()
OUT = str(FIG / "fig09_ladders.png")

fig, ax = plt.subplots(figsize=(11.5, 7.4))
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis("off")

# column centres
xL, xR, wbox = 26, 74, 26
rows_y = [70, 44, 18]  # top, middle, bottom box centres
hbox = 14


def box(xc, yc, lines, shade=None, big=False):
    ax.add_patch(
        FancyBboxPatch(
            (xc - wbox / 2, yc - hbox / 2),
            wbox,
            hbox,
            boxstyle="round,pad=0.3,rounding_size=1.4",
            fc=(shade or "white"),
            ec="black",
            lw=1.6,
        )
    )
    ax.text(
        xc,
        yc,
        lines,
        ha="center",
        va="center",
        fontsize=12.5 if big else 11.5,
        fontweight="bold",
    )


# left ladder (Fe-O), right ladder (Fe-S)
box(xL, rows_y[0], "Fe$_2$O$_3$\nhematite")
box(xR, rows_y[0], "FeS$_2$\npyrite")
box(xL, rows_y[1], "Fe$_3$O$_4$\nmagnetite (spinel)", shade=bw.SHADE_LIGHT, big=True)
box(xR, rows_y[1], "Fe$_3$S$_4$\ngreigite (thiospinel)", shade=bw.SHADE_LIGHT, big=True)
box(xL, rows_y[2], "Fe metal")
box(xR, rows_y[2], "Fe metal")


# vertical arrows + boundary fugacities between rows
def ladder(xc, vtop, vbot):
    for ya, yb in ((rows_y[0], rows_y[1]), (rows_y[1], rows_y[2])):
        ax.annotate(
            "",
            xy=(xc, yb + hbox / 2 + 0.5),
            xytext=(xc, ya - hbox / 2 - 0.5),
            arrowprops=dict(arrowstyle="<|-|>", color="black", lw=1.4),
        )
    ax.text(
        xc + 5.5, (rows_y[0] + rows_y[1]) / 2, vtop, ha="left", va="center", fontsize=10
    )
    ax.text(
        xc + 5.5, (rows_y[1] + rows_y[2]) / 2, vbot, ha="left", va="center", fontsize=10
    )


ladder(xL, "log $f$(O$_2$)\n= −68.7", "log $f$(O$_2$)\n= −87.9")
ladder(xR, "log $f$(S$_2$)\n= −20.5", "log $f$(S$_2$)\n= −51.0")

# central dashed connector across the spinel row
ax.annotate(
    "",
    xy=(xR - wbox / 2, rows_y[1]),
    xytext=(xL + wbox / 2, rows_y[1]),
    arrowprops=dict(arrowstyle="<|-|>", color="black", lw=1.3, ls="--"),
)
ax.text(
    50,
    rows_y[1] + 8.6,
    "Fe(II,III) SPINEL  Fe$_3$X$_4$\nmagnetite ≙ greigite — both ferrimagnetic",
    ha="center",
    va="center",
    fontsize=10.5,
    fontweight="bold",
    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=1.0),
)

# column headers + endmember annotations
ax.text(
    xL,
    90,
    "Fe–O   (log $f$ O$_2$)",
    ha="center",
    va="center",
    fontsize=12,
    fontweight="bold",
)
ax.text(
    xR,
    90,
    "Fe–S   (log $f$ S$_2$)",
    ha="center",
    va="center",
    fontsize=12,
    fontweight="bold",
)
ax.text(
    50,
    rows_y[0],
    "fully oxidised /\nsulfidised:\nhematite ‖ pyrite",
    ha="center",
    va="center",
    fontsize=8.5,
    style="italic",
    color=bw.FOOT_GREY,
)
ax.text(
    50,
    rows_y[2],
    "monochalcogenide FeX\n(wüstite ‖ pyrrhotite,\nFe-rich / metastable)",
    ha="center",
    va="center",
    fontsize=8.5,
    style="italic",
    color=bw.FOOT_GREY,
)

# left vertical axis arrow
ax.annotate(
    "",
    xy=(4, 80),
    xytext=(4, 12),
    arrowprops=dict(arrowstyle="-|>", color="black", lw=1.6),
)
ax.text(
    1.5,
    46,
    "increasing chalcogen fugacity →",
    rotation=90,
    ha="center",
    va="center",
    fontsize=10,
)

ax.set_title(
    "Fig. 9 — Greigite is the sulfur analogue of magnetite: parallel Fe–O and Fe–S redox ladders (298 K)\n"
    "the Fe(II,III) spinel Fe$_3$X$_4$ is the stable intermediate in both systems",
    fontsize=11.5,
    fontweight="bold",
)
fig.tight_layout(rect=[0, 0.03, 1, 1])
fig.savefig(OUT, dpi=200, bbox_inches="tight")
print("wrote", OUT)
