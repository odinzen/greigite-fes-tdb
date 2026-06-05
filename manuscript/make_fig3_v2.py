#!/usr/bin/env python
"""Fig. 3 (Fe-S-O) — TWO boundary cases, native-S restored, large fonts.

HYBRID construction (honest label): Fe-bearing sulfide/oxide/sulfate fields are
a grand-potential argmin over log f(O2) x log f(S2) at 298.15 K using
formation data from thermo_data.py (NIST-JANAF oxides; GS1992 pyrrhotite
Fe0.877S; Subramani 2020 greigite; CRC FeSO4). This is NOT yet a single
Fe-S-O CALPHAD minimization — that needs the Dilner 2017 Ca-Fe-O-S TDB
(blocked on fetching FeCaOS.TDB.txt). Pyrrhotite is the provisional Fe0.877S
of the prototype; it will be unified to the engine free composition once the
single TDB lands.

Two panels mirror Fig. 2b:
  LEFT  LOWER  greigite -1 sigma (-21.9 kJ/Fe3S4) -> pyrrhotite eliminated
  RIGHT UPPER  greigite +1 sigma (+21.9 kJ/Fe3S4) -> pyrrhotite retained

Native-S saturation line at log f(S2) = -13.96 (298 K) is the ENGINE value
(pycalphad GHSERSS vs S2 reference F15281T), restored here on Fe-S-O ONLY.
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

HERE = Path(__file__).resolve().parent  # = manuscript/
ROOT = HERE.parent  # = repo root
FIG = ROOT / "artifacts" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE))
import reactions as rx

RR = rx.R
LN10 = rx.LN10
T = 298.15
LOGFS2_SAT = -13.96  # engine: GHSERSS vs F15281T at 298.15 K (traceable)

# FeSO4(cr): CRC Handbook (Wagman/NBS-1982): dHf=-928.4 kJ/mol, S=107.5 J/mol/K
rx.PHASE["FeSO4"] = (
    lambda Tx: np.full_like(np.atleast_1d(Tx), 100.6, float),
    -928400.0,
    107.5,
)

PHASES = ["Fe", "Po", "Gr", "FeS2", "Fe3O4", "Fe2O3", "FeSO4"]
LABELS = {
    "Fe": "Fe",
    "Po": "Fe$_{0.877}$S\npyrrhotite",
    "Gr": "Fe$_3$S$_4$\ngreigite",
    "FeS2": "FeS$_2$\npyrite",
    "Fe3O4": "Fe$_3$O$_4$\nmagnetite",
    "Fe2O3": "Fe$_2$O$_3$\nhematite",
    "FeSO4": "FeSO$_4$",
}
NU_S2 = {
    "Fe": 0.0,
    "Po": 0.5 / 0.877,
    "Gr": 2 / 3.0,
    "FeS2": 1.0,
    "Fe3O4": 0.0,
    "Fe2O3": 0.0,
    "FeSO4": 0.5,
}
NU_O2 = {
    "Fe": 0.0,
    "Po": 0.0,
    "Gr": 0.0,
    "FeS2": 0.0,
    "Fe3O4": 2 / 3.0,
    "Fe2O3": 0.75,
    "FeSO4": 2.0,
}
STOICH = {
    "Po": {"Fe": -1, "S2": -0.5 / 0.877, "Po": 1.0 / 0.877},
    "Gr": {"Fe": -1, "S2": -2 / 3.0, "Gr": 1 / 3.0},
    "FeS2": {"Fe": -1, "S2": -1.0, "FeS2": 1.0},
    "Fe3O4": {"Fe": -1, "O2": -2 / 3.0, "Fe3O4": 1 / 3.0},
    "Fe2O3": {"Fe": -1, "O2": -0.75, "Fe2O3": 0.5},
    "FeSO4": {"Fe": -1, "S2": -0.5, "O2": -2.0, "FeSO4": 1.0},
}
COLORS = {
    "Fe": "#cfd8dc",
    "Po": "#ffcc80",
    "Gr": "#a5d6a7",
    "FeS2": "#90caf9",
    "Fe3O4": "#b39ddb",
    "Fe2O3": "#ef9a9a",
    "FeSO4": "#ffe082",
}

O = np.linspace(-120.0, 0.0, 700)
S = np.linspace(-60.0, 0.0, 700)
OO, SS = np.meshgrid(O, S)

CASES = [
    (
        -21900.0,
        "LOWER bound\ngreigite −1σ (−21.9 kJ/Fe$_3$S$_4$)\npyrrhotite eliminated",
    ),
    (
        +21900.0,
        "UPPER bound\ngreigite +1σ (+21.9 kJ/Fe$_3$S$_4$)\nmaximum pyrrhotite field",
    ),
]

plt.rcParams.update({"font.size": 14})
fig, axes = plt.subplots(1, 2, figsize=(16.5, 8.6), sharey=True)

for ax, (off, title) in zip(axes, CASES):
    rx.set_greigite_dHf_offset(off)
    Ta = np.array([T])
    dGf = {"Fe": 0.0}
    for p in ["Po", "Gr", "FeS2", "Fe3O4", "Fe2O3", "FeSO4"]:
        dGf[p] = rx.reaction_dG(STOICH[p], Ta)[0]
    Phi = np.zeros((len(PHASES), SS.shape[0], SS.shape[1]))
    for i, p in enumerate(PHASES):
        Phi[i] = dGf[p] - NU_S2[p] * RR * T * LN10 * SS - NU_O2[p] * RR * T * LN10 * OO
    field = np.argmin(Phi, axis=0)
    cmap = ListedColormap([COLORS[p] for p in PHASES])
    norm = BoundaryNorm(np.arange(-0.5, len(PHASES) + 0.5, 1), cmap.N)
    ax.pcolormesh(OO, SS, field, cmap=cmap, norm=norm, shading="auto")
    ax.contour(
        OO,
        SS,
        field,
        levels=np.arange(0.5, len(PHASES) - 0.5 + 1, 1),
        colors="k",
        linewidths=1.1,
    )
    for i, p in enumerate(PHASES):
        mask = field == i
        if mask.sum() > 200:
            ax.text(
                OO[mask].mean(),
                SS[mask].mean(),
                LABELS[p],
                ha="center",
                va="center",
                fontsize=15,
                fontweight="bold",
            )
    # native-S saturation (engine value) restored on Fe-S-O
    ax.axhline(LOGFS2_SAT, color="#7a4a00", ls="--", lw=2.2)
    ax.fill_between(
        [-120, 0], LOGFS2_SAT, 0, color="#f4e3b0", alpha=0.45, hatch="xx", lw=0
    )
    ax.text(
        -60,
        (LOGFS2_SAT) / 2 + 0.5,
        "native S condenses  (S$_2$ saturation, log $f$(S$_2$)=−13.96, 298 K)",
        ha="center",
        va="center",
        fontsize=12.5,
        color="#5a3600",
        fontweight="bold",
    )
    ax.set_title(title, fontsize=15)
    ax.set_xlabel("log $f$(O$_2$), bar", fontsize=16)
    ax.tick_params(labelsize=14)
    ax.set_xlim(-120, 0)
    ax.set_ylim(-60, 0)
    rx.set_greigite_dHf_offset(0.0)

axes[0].set_ylabel("log $f$(S$_2$), bar", fontsize=16)
fig.suptitle(
    "Fig. 3 — Fe–S–O predominance at 298.15 K: two boundary cases (hybrid: measured/JANAF data)\n"
    "native-S saturation restored (engine value); greigite borders magnetite, pyrite borders hematite/FeSO$_4$",
    fontsize=15,
    fontweight="bold",
)
note = (
    "HYBRID, not yet a single-database minimization: sulfides (Subramani/Shumway greigite, GS1992 pyrrhotite Fe$_{0.877}$S — provisional), "
    "oxides NIST-JANAF, FeSO$_4$ CRC.  Pyrrhotite composition to be unified to the engine free comp. (FeS$_{1.00}$–FeS$_{1.14}$; observed product Fe$_7$S$_8$) "
    "once the Dilner-2017 Ca–Fe–O–S TDB is fetched.  Sulfur is S$_2$(g); native-S saturation = −13.96 (engine, 298 K)."
)
fig.text(0.5, 0.01, note, ha="center", va="bottom", fontsize=11, color="#222")
fig.tight_layout(rect=[0, 0.065, 1, 0.92])
fig.savefig(str(FIG / "fig3_FeSO_boundary_cases.png"), dpi=170)
print("wrote fig3_FeSO_boundary_cases.png")
