import sys
from pathlib import Path
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from itertools import combinations

HERE = Path(__file__).resolve().parent  # = manuscript/
ROOT = HERE.parent  # = repo root
FIG = ROOT / "artifacts" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE))
import reactions as rx

RR = rx.R
LN10 = rx.LN10
T = 298.15

# --- FeSO4(cr) added to the engine -----------------------------------------
# Anhydrous ferrous sulfate standard-state values, CRC Handbook of Chemistry
# and Physics (Wagman/NBS-1982 tabulation): dHf298 = -928.4 kJ/mol,
# S298 = 107.5 J/mol/K. Cp is immaterial at exactly T_REF (integral = 0).
FESO4_DHF = -928400.0
FESO4_S = 107.5
rx.PHASE["FeSO4"] = (
    lambda Tx: np.full_like(np.atleast_1d(Tx), 100.6, float),
    FESO4_DHF,
    FESO4_S,
)

rx.set_greigite_dHf_offset(+21900.0)  # upper bound -> pyrrhotite retained

PHASES = ["Fe", "Po", "Gr", "FeS2", "Fe3O4", "Fe2O3", "FeSO4"]
LABELS = {
    "Fe": "Fe",
    "Po": "Fe$_{(1-x)}$S\npyrrhotite",
    "Gr": "Fe$_3$S$_4$\ngreigite",
    "FeS2": "FeS$_2$\npyrite",
    "Fe3O4": "Fe$_3$O$_4$\nmagnetite",
    "Fe2O3": "Fe$_2$O$_3$\nhematite",
    "FeSO4": "FeSO$_4$\nferrous sulfate",
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

Ta = np.array([T])
dGf = {"Fe": 0.0}
for p in ["Po", "Gr", "FeS2", "Fe3O4", "Fe2O3", "FeSO4"]:
    dGf[p] = rx.reaction_dG(STOICH[p], Ta)[0]

# y to -60 so Fe appears at low f(S2); x to 0 so FeSO4 caps the
# high-f(O2)/high-f(S2) corner (the "orange above hematite" the reviewer asked
# for). Sulfur is S2 GAS throughout -- no native-S condensation line.
O = np.linspace(-120.0, 0.0, 901)
S = np.linspace(-60.0, 0.0, 901)
OO, SS = np.meshgrid(O, S)
Phi = np.zeros((len(PHASES), SS.shape[0], SS.shape[1]))
for i, p in enumerate(PHASES):
    Phi[i] = dGf[p] - NU_S2[p] * RR * T * LN10 * SS - NU_O2[p] * RR * T * LN10 * OO
field = np.argmin(Phi, axis=0)

cmap = ListedColormap([COLORS[p] for p in PHASES])
norm = BoundaryNorm(np.arange(-0.5, len(PHASES) + 0.5, 1), cmap.N)
fig, ax = plt.subplots(figsize=(9.2, 6.6))
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
    if mask.sum() > 60:
        ax.text(
            OO[mask].mean(),
            SS[mask].mean(),
            LABELS[p],
            ha="center",
            va="center",
            fontsize=10,
            fontweight="bold",
        )

# --- Locate all 3-phase junctions; check for mag/py/hem triple point --------
mag, py, hem = (PHASES.index("Fe3O4"), PHASES.index("FeS2"), PHASES.index("Fe2O3"))
junc = {}
for r in range(1, field.shape[0] - 1):
    for c in range(1, field.shape[1] - 1):
        w = frozenset(field[r - 1 : r + 2, c - 1 : c + 2].ravel())
        if len(w) >= 3:
            for tri in combinations(sorted(w), 3):
                junc.setdefault(tri, (OO[r, c], SS[r, c]))
mph = tuple(sorted((mag, py, hem)))
if mph in junc:
    ax.plot(*junc[mph], "k*", ms=15, mfc="yellow", mec="k", zorder=6)
    ax.annotate(
        "Fe$_3$O$_4$/FeS$_2$/Fe$_2$O$_3$ triple point",
        xy=junc[mph],
        xytext=(junc[mph][0] - 18, junc[mph][1] + 8),
        fontsize=8,
        arrowprops=dict(arrowstyle="->"),
    )
    tp_status = "present at %s" % (junc[mph],)
else:
    # Mark the oxide/sulfide junctions that DO exist (greigite intervenes).
    for tri in [
        tuple(sorted((PHASES.index("Gr"), py, hem))),
        tuple(sorted((PHASES.index("Gr"), mag, hem))),
    ]:
        if tri in junc:
            ax.plot(*junc[tri], "kP", ms=9, mfc="white", mec="k", zorder=6)
    tp_status = (
        "ABSENT - greigite separates magnetite from pyrite; the "
        "actual junctions are Gr/FeS2/Fe2O3 and Gr/Fe3O4/Fe2O3 "
        "(white markers)."
    )

ax.set_xlabel("log $f$(O$_2$)")
ax.set_ylabel("log $f$(S$_2$)")
ax.set_xlim(-120, 0)
ax.set_ylim(-60, 0)
ax.set_title(
    "Fe-S-O predominance at T = 298.15 K (greigite upper d$H_f$ "
    "bound, +21.9 kJ/Fe$_3$S$_4$)\nFeSO$_4$ included; sulfur treated "
    "as S$_2$ gas (no condensed-S phase)",
    fontsize=10,
)

note = (
    "No magnetite/pyrite/hematite triple point: greigite intervenes\n"
    "between pyrite and the oxides at 298 K (fixed-Fe$_{0.877}$S engine).\n"
    "Actual junctions: Gr/FeS$_2$/Fe$_2$O$_3$ and Gr/Fe$_3$O$_4$/Fe$_2$O$_3$ (+)."
)
ax.text(
    0.985,
    0.04,
    note,
    transform=ax.transAxes,
    ha="right",
    va="bottom",
    fontsize=7.0,
    style="italic",
    bbox=dict(boxstyle="round", fc="white", ec="0.6", alpha=0.9),
)

fig.text(
    0.5,
    0.005,
    "Sulfur treated as S$_2$ gas (no condensed-S phase); "
    "pyrrhotite/greigite boundaries provisional (fixed Fe$_{0.877}$S). "
    "FeSO$_4$(cr): $\\Delta_f H^\\circ=-928.4$ kJ/mol, $S^\\circ=107.5$ "
    "J/mol/K (CRC Handbook).\nOther phases: NIST-JANAF (Chase 1998); "
    "GS1992 pyrrhotite; Shumway 2022 / Subramani 2020 greigite.",
    ha="center",
    fontsize=6.6,
    style="italic",
)
fig.tight_layout(rect=[0, 0.04, 1, 1])
fig.savefig(str(FIG / "explore_feso_predominance_pyrrhotite.png"), dpi=150)
plt.close(fig)
rx.set_greigite_dHf_offset(0.0)

present = [PHASES[k] for k in sorted(set(field.ravel()))]
print("fields present:", present)
print("mag/py/hem triple point:", tp_status)
print(
    "area frac:",
    {
        p: round((field == i).mean(), 3)
        for i, p in enumerate(PHASES)
        if (field == i).any()
    },
)
