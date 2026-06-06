#!/usr/bin/env python
"""Render Fig. 2b — the TWO labelled boundary cases, from the ENGINE.

Reads the all-compound boundary TDBs built by
engine/build_boundary_tdbs.py and draws a two-panel Fe-S
log f(S2)-vs-T predominance diagram:

  LEFT  = LOWER  (greigite -1s, pyrite +1s, pyrrhotite +1s) -> pyrrhotite eliminated
  RIGHT = UPPER  (greigite +1s, pyrite -1s, pyrrhotite -1s) -> max pyrrhotite field

The two panels ARE the on-figure uncertainty representation (Kristina's
"errors"): together they bound the +/-1 sigma all-compound envelope
(greigite +/-7.3, pyrite +/-6.0, pyrrhotite +/-3.5 kJ per f.u.).

Regression fixes implemented here:
  * boundary cases restored (both panels), built through pycalphad not the prototype;
  * error envelope shown (the two panels);
  * NO native-S cap on this Fe-S diagram (S2 is the axis variable; native-S
    saturation belongs on Fe-S-O / Fig. 3). Axis kept honest: a note states the
    S2 axis runs to native-S saturation, shown on Fig. 3, so pyrite's roof is
    not drawn as unbounded.
  * pyrrhotite composition: BOTH shown/labelled — equilibrium free composition
    (engine: troilite-end FeS1.00-1.14) AND the observed DSC/XRD product Fe7S8.
  * DSC decomposition markers: greigite stable to 300C/573K; pyritization onset
    317C/590K (Kristina DSC/XRD heating.docx).

Usage (probe each case, then plot):
    python render_fig2b.py probe lower     # probes LOWER, caches fields
    python render_fig2b.py probe upper     # probes UPPER, caches fields
    python render_fig2b.py plot            # draws fig2b_boundary_cases.png
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from pycalphad import Database, equilibrium, variables as v

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TDB = ROOT / "artifacts" / "tdb"
FIG = ROOT / "artifacts" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
CACHE = HERE / "fig2b_fields.json"

R = 8.31451
LN10 = math.log(10.0)
COMPS = ["FE", "S", "VA"]
CONDENSED = ["BCC_A2", "FCC_A1", "PYRRHOTITE", "PYRITE", "GREIGITE", "ORTHORHOMBIC_S"]
S2_REF = "F15281T"
T_GRID = [300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0]
# denser near the pyrite stoichiometry (X_S=2/3) so the thin pyrite roof is caught
XS_GRID = np.concatenate([np.linspace(0.34, 0.64, 42), np.linspace(0.642, 0.6664, 26)])

CASES = {
    "lower": dict(
        tdb="fes_greigite_boundary_lower.tdb",
        title="LOWER bound  (greigite −1σ, pyrrhotite +1σ; pyrite central)\npyrrhotite eliminated — greigite preempts",
    ),
    "upper": dict(
        tdb="fes_greigite_boundary_upper.tdb",
        title="UPPER bound  (greigite +1σ, pyrrhotite −1σ; pyrite central)\nmaximum pyrrhotite field",
    ),
}
COLORS = {
    "BCC_A2": "#d9d9d9",
    "PYRRHOTITE": "#e8d6a0",
    "GREIGITE": "#5fae77",
    "PYRITE": "#bcd4ec",
}
LABELS = {
    "BCC_A2": "Fe (bcc)",
    "PYRRHOTITE": "pyrrhotite",
    "GREIGITE": "Fe$_3$S$_4$ greigite",
    "PYRITE": "FeS$_2$ pyrite",
}


def _g_s2_ref(db, T):
    from pycalphad.variables import T as T_VAR

    return float(db.symbols[S2_REF].subs({T_VAR: float(T)}))


def probe(case_key):
    c = CASES[case_key]
    db = Database(str(TDB / c["tdb"]))
    fields = {}
    for T in T_GRID:
        gref = _g_s2_ref(db, T)
        for x in XS_GRID:
            eq = equilibrium(
                db, COMPS, CONDENSED, {v.X("S"): float(x), v.T: T, v.P: 101325, v.N: 1}
            )
            phases = sorted({p for p in np.ravel(eq.Phase.values) if p})
            try:
                muS = float(np.ravel(eq.MU.sel(component="S").values)[0])
            except Exception:
                muS = float(np.ravel(eq.MU.values)[-1])
            logfs2 = (2.0 * muS - gref) / (R * T * LN10)
            for ph in phases:
                fields.setdefault(ph, {}).setdefault(str(int(T)), []).append(logfs2)
    out = {
        ph: {T: [min(vals), max(vals)] for T, vals in d.items()}
        for ph, d in fields.items()
    }
    cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}
    cache[case_key] = out
    CACHE.write_text(json.dumps(cache, indent=2))
    print(f"probed {case_key}: phases = {sorted(out)}")


def plot():
    cache = json.loads(CACHE.read_text())
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 6.4), sharey=True)
    ymin, ymax = -56, 0
    for ax, key in zip(axes, ("lower", "upper")):
        fields = cache[key]
        for ph in ("BCC_A2", "PYRRHOTITE", "GREIGITE", "PYRITE"):
            if ph not in fields:
                continue
            Ts = sorted(int(t) for t in fields[ph])
            lo = np.array([fields[ph][str(t)][0] for t in Ts])
            hi = np.array([fields[ph][str(t)][1] for t in Ts])
            hi = np.where(
                ph == "PYRITE", np.minimum(hi, ymax), hi
            )  # cap roof at axis top
            ax.fill_between(
                Ts, lo, hi, color=COLORS[ph], alpha=0.9, lw=0, label=LABELS[ph]
            )
        # DSC/XRD decomposition markers (Kristina, heating.docx)
        ax.axvline(573, color="#b03030", ls=":", lw=1.3)
        ax.axvline(590, color="#b03030", ls="--", lw=1.3)
        ax.text(
            574,
            ymin + 2,
            "greigite stable\n≤300°C/573K",
            color="#b03030",
            fontsize=7.5,
            rotation=90,
            va="bottom",
            ha="left",
        )
        ax.text(
            591,
            ymin + 2,
            "pyritization onset\n317°C/590K",
            color="#b03030",
            fontsize=7.5,
            rotation=90,
            va="bottom",
            ha="left",
        )
        ax.set_title(CASES[key]["title"], fontsize=10)
        ax.set_xlabel("Temperature, K")
        ax.set_xlim(300, 600)
        ax.set_ylim(ymin, ymax)
    axes[0].set_ylabel(r"log $f$(S$_2$), bar")
    # legend incl. both pyrrhotite compositions
    handles = [
        Patch(facecolor=COLORS[p], label=LABELS[p])
        for p in ("BCC_A2", "PYRRHOTITE", "GREIGITE", "PYRITE")
    ]
    axes[1].legend(handles=handles, loc="lower right", fontsize=8.5, framealpha=0.95)
    fig.suptitle(
        "Fig. 2b — Fe–S equilibrium: two boundary cases bounding the ±1σ all-compound envelope\n"
        "(greigite ±7.3, pyrite ±6.0, pyrrhotite ±3.5 kJ/f.u.; pycalphad on Dilner-2015 + measured greigite)",
        fontsize=11,
    )
    note = (
        "Pyrrhotite composition: equilibrium (free) = troilite-end FeS$_{1.00}$–FeS$_{1.14}$ (engine, widens with T); "
        "observed DSC/XRD decomposition product = pyrrhotite-3T Fe$_7$S$_8$ (FeS$_{1.14}$).\n"
        "Sulfur is S$_2$(g) at all T — the f(S$_2$) axis runs up to native-S saturation, which is drawn on the Fe–S–O diagram (Fig. 3), not here. "
        "Pyrite roof clipped at axis top, not unbounded."
    )
    fig.text(0.5, 0.005, note, ha="center", va="bottom", fontsize=7.4, color="#333333")
    fig.tight_layout(rect=[0, 0.055, 1, 0.93])
    out = FIG / "fig2b_boundary_cases.png"
    fig.savefig(out, dpi=200)
    print(f"wrote {out}")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "plot"
    if cmd == "probe":
        probe(sys.argv[2])
    else:
        plot()
