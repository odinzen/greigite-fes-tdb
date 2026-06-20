#!/usr/bin/env python
"""Fig. 1 - the four-step reproducible pipeline that produces this paper's
databases and figures.

No platform/product name: the system is described purely by the four steps it
runs. The Gibbs energy function (TDB format) is the common data currency passed
between steps; provenance is recorded at each step. This is the subset of a
larger research codebase that this study actually exercises, packaged on its own
under a permissive licence.

Step 1 loads the typed knowledge-graph dump (engine/provenance_manifest.json):
the TDB fetch URLs, the measured greigite values, and the build recipe. The
later steps fetch the literature TDBs from those URLs and build from the
measured values. The upstream OCR + agentic extraction that populated the
knowledge graph lives in the internal codebase and is not re-run here; this open
subset builds from the dumped manifest, which keeps every value auditable back
to its primary source.
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
S = EN + "S"  # "–S"
CFOS = "Ca" + EN + "Fe" + EN + "O" + EN + "S"
FESO = "Fe" + EN + "S(" + EN + "O)"
# (num, title, what-it-does, tool, this-work output)
STEPS = [
    (
        "1",
        "Provenance manifest",
        "Load the typed knowledge-graph dump - TDB fetch\n"
        "URLs, the build recipe, and greigite values\n"
        "(Subramani 2020; Shumway 2022) OCR-extracted\n"
        "from the primary calorimetry papers.",
        "tool: provenance_manifest.json (typed KG export)",
        "TDB URLs + measured values\n+ build recipe (DAG)",
    ),
    (
        "2",
        "Database lookup",
        "Fetch the assessed Fe" + S + " (Dilner 2015) and\n"
        + CFOS + " (Dilner 2017) databases from\n"
        "the URLs recorded in the manifest.",
        "tool: TDB-DB fetch",
        "Dilner 2015 Fe" + S + " +\nDilner 2017 " + CFOS + " (.tdb)",
    ),
    (
        "3",
        "TDB stitching",
        "Build the greigite line compound from the\n"
        "measured values; dedupe and graft greigite\n"
        "(+ pyrite) onto a common SGTE91 reference.",
        "tool: Python build / graft",
        "single " + FESO + " .tdb\n(byte-reproducible)",
    ),
    (
        "4",
        "Predominance plotting",
        "Grand-potential minimisation (after Holland\n"
        "1959) with all-compound ±σ envelope and\n"
        "native-sulfur saturation capping.",
        "tool: pycalphad 0.11.1",
        "Fe" + S + " and Fe" + S + EN + "O predominance\n"
        "diagrams (the paper's figures)",
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
    "A single typed provenance manifest ties the steps together: every value is read from it and "
    "stays traceable to its primary source.",
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
for i, (num, title, body, tool, out) in enumerate(STEPS):
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

# provenance-manifest spine on the right: the typed record every step reads
# from and writes back to, keeping each value traceable to its primary source
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
    "provenance manifest  =  every value typed and traceable to its source",
    ha="center",
    va="center",
    rotation=90,
    fontsize=11,
    fontweight="bold",
    color="black",
)

ax.text(
    48,
    4.0,
    "Each step records typed provenance back to its primary source; packaged for reproduction - MIT-licensed.",
    ha="center",
    va="top",
    fontsize=9,
    color="#444",
)

fig.tight_layout()
fig.savefig(OUT, dpi=170, bbox_inches="tight")
print("wrote", OUT)
