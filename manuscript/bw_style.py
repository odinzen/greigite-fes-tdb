"""Shared black-and-white figure style for the greigite manuscript.

Matches the reviewer-preferred look: white background, black boundary lines,
in-field text labels, dashed/dotted lines for uncertainty bounds, hatching for
special regions. No colour. Import and call ``apply()`` at the top of each
figure script; use the GREY constants and ``hatch_kw`` helper for consistency.
"""

import matplotlib

PHASE_GREY = "#000000"  # boundary lines, labels
FIELD_FILL = "#ffffff"  # all fields white
SHADE_LIGHT = "#e6e6e6"  # optional light-grey fill for one emphasised field
SHADE_MID = "#bdbdbd"
FOOT_GREY = "#444444"  # footnote text


def apply():
    matplotlib.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
            "font.size": 12,
            "axes.edgecolor": "black",
            "axes.linewidth": 1.0,
            "axes.facecolor": "white",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.grid": False,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.color": "black",
            "ytick.color": "black",
            "text.color": "black",
            "axes.labelcolor": "black",
        }
    )


def hatch_kw(hatch="////"):
    """kwargs for a hatched, unfilled region (black hatch on white)."""
    return dict(facecolor="none", edgecolor="black", hatch=hatch, linewidth=0.0)
