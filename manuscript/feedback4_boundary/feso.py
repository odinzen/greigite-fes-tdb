import sys, numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
MAN = HERE.parent
ROOT = MAN.parent
TDB = ROOT / "artifacts" / "tdb"
FIG = ROOT / "artifacts" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(MAN))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import reactions as rx, thermo_data as td

RR = rx.R
LN10 = rx.LN10
T = 298.15
# add FeSO4 (CRC)
rx.PHASE["FeSO4"] = (
    lambda Tx: np.full_like(np.atleast_1d(Tx), 100.6, float),
    -928400.0,
    107.5,
)

PHASES = ["Fe", "Po", "Gr", "FeS2", "Fe3O4", "Fe2O3", "FeSO4"]
LAB = {
    "Fe": "Fe",
    "Po": "Fe$_{1-x}$S\npyrrhotite",
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
STO = {
    "Po": {"Fe": -1, "S2": -0.5 / 0.877, "Po": 1.0 / 0.877},
    "Gr": {"Fe": -1, "S2": -2 / 3.0, "Gr": 1 / 3.0},
    "FeS2": {"Fe": -1, "S2": -1.0, "FeS2": 1.0},
    "Fe3O4": {"Fe": -1, "O2": -2 / 3.0, "Fe3O4": 1 / 3.0},
    "Fe2O3": {"Fe": -1, "O2": -0.75, "Fe2O3": 0.5},
    "FeSO4": {"Fe": -1, "S2": -0.5, "O2": -2.0, "FeSO4": 1.0},
}
COL = {
    "Fe": "#cfd8dc",
    "Po": "#ffcc80",
    "Gr": "#a5d6a7",
    "FeS2": "#90caf9",
    "Fe3O4": "#b39ddb",
    "Fe2O3": "#ef9a9a",
    "FeSO4": "#ffe082",
}

# native-S saturation: S2(g) <-> 2 S(orth);  S(orth) ref (dHf0, S 32.054); S2(g) dHf128600 S228.165
dG_S2cond = (
    -128600.0 + (2 * 32.054 - 228.165) * (-1) * T
)  # = -128600 +164.057*T  (J), ΔG° for S2->2S
logfS2_sat = dG_S2cond / (RR * T * LN10)
print("native-S saturation log f(S2) =", round(logfS2_sat, 2))

O = np.linspace(-120.0, 0.0, 701)
S = np.linspace(-60.0, 0.0, 701)
OO, SS = np.meshgrid(O, S)


def field_for(offset):
    rx.set_greigite_dHf_offset(offset)
    dGf = {"Fe": 0.0}
    for p in ["Po", "Gr", "FeS2", "Fe3O4", "Fe2O3", "FeSO4"]:
        dGf[p] = rx.reaction_dG(STO[p], np.array([T]))[0]
    rx.set_greigite_dHf_offset(0.0)
    Phi = np.zeros((len(PHASES),) + SS.shape)
    for i, p in enumerate(PHASES):
        Phi[i] = dGf[p] - NU_S2[p] * RR * T * LN10 * SS - NU_O2[p] * RR * T * LN10 * OO
    return np.argmin(Phi, axis=0)


cmap = ListedColormap([COL[p] for p in PHASES])
norm = BoundaryNorm(np.arange(-0.5, len(PHASES) + 0.5, 1), cmap.N)
fig, axes = plt.subplots(1, 2, figsize=(13, 5.7), sharey=True)
titles = {
    -21900.0: "Lower bound — greigite +1σ stable (−21.9 kJ/Fe₃S₄)\n→ no pyrrhotite field",
    21900.0: "Upper bound — greigite −1σ (+21.9 kJ/Fe₃S₄)\n→ pyrrhotite field retained",
}
for ax, off in zip(axes, [-21900.0, 21900.0]):
    field = field_for(off)
    ax.pcolormesh(OO, SS, field, cmap=cmap, norm=norm, shading="auto")
    ax.contour(
        OO,
        SS,
        field,
        levels=np.arange(0.5, len(PHASES) - 0.5 + 1, 1),
        colors="k",
        linewidths=1.0,
    )
    for i, p in enumerate(PHASES):
        m = field == i
        if m.sum() > 200:
            ax.text(
                OO[m].mean(),
                SS[m].mean(),
                LAB[p],
                ha="center",
                va="center",
                fontsize=9,
                fontweight="bold",
            )
    # native S saturation cap
    ax.axhspan(logfS2_sat, 0, color="#fff176", alpha=0.45, lw=0)
    ax.axhline(logfS2_sat, color="#b8860b", lw=1.8, ls="--")
    ax.text(
        -119,
        logfS2_sat + 0.6,
        "native S(s) saturation — condensed S₈ present above (log f(S₂)=%.1f, 298 K)"
        % logfS2_sat,
        fontsize=7.4,
        color="#7a5c00",
        va="bottom",
    )
    ax.set_xlabel("log $f$(O$_2$)")
    ax.set_title(titles[off], fontsize=9.5)
    ax.set_xlim(-120, 0)
    ax.set_ylim(-60, 0)
axes[0].set_ylabel("log $f$(S$_2$)")
fig.suptitle(
    "Fe–S–O predominance at 298 K — two error-bounding cases (hybrid: engine/prototype sulfides + JANAF/CRC oxides)\n"
    "Sulfur axis = log f(S₂); native-sulfur saturation reinstated (yellow cap)",
    fontsize=10.5,
)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(str(FIG / "fig_feso_bounds.png"), dpi=150)
print("SAVED fig_feso_bounds.png")
