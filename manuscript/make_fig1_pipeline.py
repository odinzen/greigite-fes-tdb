#!/usr/bin/env python
"""Fig. 1 - the four-step reproducible pipeline that produces this paper's
databases and figures.

No platform/product name: the system is described purely by the four steps it
runs. The Gibbs energy function (TDB format) is the common data currency passed
between steps; provenance is recorded at each step. This is the subset of a
larger research codebase that this study actually exercises, packaged on its own
under a permissive licence.
"""

from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

HERE = Path(__file__).resolve().parent  # = manuscript/
ROOT = HERE.parent  # = repo root
FIG = ROOT / "artifacts" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

OUT = str(FIG / "fig01_pipeline.png")

EN = "–"  # en dash
# (num, title, facecolor, edgecolor, what-it-does, tool, this-work output)
STEPS = [
    (
        "1",
        "Database lookup",
        "#cfe3f7",
        "#1f4e79",
        "Retrieve the assessed Fe"
        + EN
        + "S and\nCa"
        + EN
        + "Fe"
        + EN
        + "O"
        + EN
        + "S thermodynamic databases\nfrom the literature index.",
        "tool: TDB-DB lookup",
        "Dilner 2015 Fe"
        + EN
        + "S +\nDilner 2017 Ca"
        + EN
        + "Fe"
        + EN
        + "O"
        + EN
        + "S (.tdb)",
    ),
    (
        "2",
        "OCR extraction",
        "#cdeccd",
        "#1e5e1e",
        "Read the greigite enthalpy, entropy\nand heat capacity from the primary\ncalorimetry papers.",
        "tool: OCR of source papers",
        "Subramani 2020; Shumway 2022\n" + EN + "> machine-readable values",
    ),
    (
        "3",
        "TDB stitching",
        "#fbe3c4",
        "#8a4b00",
        "Graft greigite (and pyrite) into one\ndatabase, dedupe functions, build\ndeterministically.",
        "tool: Python stitch script",
        "single Fe" + EN + "S(" + EN + "O) .tdb,\nbyte-reproducible",
    ),
    (
        "4",
        "Predominance plotting",
        "#e3d7f0",
        "#4a2a70",
        "Minimise the grand potential over\nthe redox plane to find the stable\nphase at each point.",
        "tool: pycalphad 0.11.1 optimiser",
        "Fe"
        + EN
        + "S and Fe"
        + EN
        + "S"
        + EN
        + "O predominance\ndiagrams (the paper's figures)",
    ),
]

fig, ax = plt.subplots(figsize=(12, 8.4))
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis("off")

ax.text(
    48,
    97.5,
    "Fig. 1 " + EN + " The four-step reproducible pipeline",
    ha="center",
    va="top",
    fontsize=16,
    fontweight="bold",
)
ax.text(
    48,
    93.4,
    "The Gibbs energy function (TDB format) is the common data currency passed between steps; "
    "provenance is recorded at each step.",
    ha="center",
    va="top",
    fontsize=10.3,
)

box_x, box_w = 4, 80
top, bot = 88, 8
gap = 2.4
n = len(STEPS)
bh = (top - bot - gap * (n - 1)) / n
mid_x = box_x + 50  # divider between left body and right output lane
for i, (num, title, fc, ec, body, tool, out) in enumerate(STEPS):
    y = top - i * (bh + gap) - bh
    ax.add_patch(
        FancyBboxPatch(
            (box_x, y),
            box_w,
            bh,
            boxstyle="round,pad=0.4,rounding_size=1.8",
            fc="white",
            ec="black",
            lw=1.6,
        )
    )
    cy = y + bh
    # numbered badge + title (top-left)
    ax.add_patch(plt.Circle((box_x + 4.0, cy - 4.3), 2.4, color="black", zorder=5))
    ax.text(
        box_x + 4.0,
        cy - 4.3,
        num,
        ha="center",
        va="center",
        color="white",
        fontsize=13,
        fontweight="bold",
        zorder=6,
    )
    ax.text(
        box_x + 8.5,
        cy - 4.3,
        title,
        ha="left",
        va="center",
        fontsize=13,
        fontweight="bold",
        color="black",
    )
    # left lane: what it does
    ax.text(
        box_x + 8.5, cy - 9.4, body, ha="left", va="top", fontsize=9.2, color="black"
    )
    # right lane: tool (bold) + output (italic)
    ax.text(
        box_x + box_w - 1.5,
        cy - 4.3,
        tool,
        ha="right",
        va="center",
        fontsize=9.6,
        fontweight="bold",
        color="black",
    )
    ax.text(
        box_x + box_w - 1.5,
        cy - 9.4,
        out,
        ha="right",
        va="top",
        fontsize=8.8,
        style="italic",
        color="#444",
    )
    if i < n - 1:
        ax.annotate(
            "",
            xy=(box_x + box_w / 2, y - 0.3),
            xytext=(box_x + box_w / 2, y - gap + 0.3),
            arrowprops=dict(arrowstyle="-|>", color="black", lw=2.0),
        )

# common-currency vertical bar on the right
bar_x = box_x + box_w + 2.5
ax.add_patch(
    FancyBboxPatch(
        (bar_x, bot),
        9.5,
        top - bot,
        boxstyle="round,pad=0.3,rounding_size=1.8",
        fc="#ececec",
        ec="black",
        lw=1.3,
    )
)
ax.text(
    bar_x + 4.75,
    (top + bot) / 2,
    "Gibbs energy function  =  common data currency",
    ha="center",
    va="center",
    rotation=90,
    fontsize=11.5,
    fontweight="bold",
    color="black",
)

ax.text(
    48,
    4.0,
    "Open subset of a larger research codebase, packaged for reproduction - "
    "MIT-licensed; provenance-tracked.",
    ha="center",
    va="top",
    fontsize=9,
    color="#444",
)

fig.tight_layout()
fig.savefig(OUT, dpi=170, bbox_inches="tight")
print("wrote", OUT)
