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
import greigite_cp as G

RR = rx.R
LN10 = rx.LN10
DECOMP = G.GREIGITE_DECOMP_T  # 530 K (SHU2022)

# --- Decomposition cap: min(decomp zero-crossing, 530 K) -------------------
# dGdecomp(T) for Fe3S4 -> pyrrhotite + pyrite over 298-560 K.
Tdec = np.linspace(298.0, 560.0, 1500)
dGdec = rx.reaction_dG({"Gr": -1 / 3.0, "Po": 0.8886, "FeS2": 0.221}, Tdec)
sign = np.sign(dGdec)
cross = np.where(np.diff(sign) != 0)[0]
zero_cross = (
    np.interp(
        0, [dGdec[cross[0]], dGdec[cross[0] + 1]], [Tdec[cross[0]], Tdec[cross[0] + 1]]
    )
    if len(cross)
    else None
)
# Engine: dGdecomp stays POSITIVE (~+14.2 kJ at 298, rising). Greigite is
# never thermodynamically unstable vs Po+Py here; its field terminates at the
# SHU2022 kinetic decomposition (530 K) -> the triangle apex.
GR_CAP = min(zero_cross, DECOMP) if zero_cross else DECOMP

PHASES = ["Fe", "Po", "Gr", "FeS2"]
LABELS = {
    "Fe": "Fe",
    "Po": "Fe$_{(1-x)}$S\npyrrhotite",
    "Gr": "Fe$_3$S$_4$\ngreigite",
    "FeS2": "FeS$_2$\npyrite",
}
NU_S2 = {"Fe": 0.0, "Po": 0.5 / 0.877, "Gr": 2 / 3.0, "FeS2": 1.0}
STOICH = {
    "Po": {"Fe": -1, "S2": -0.5 / 0.877, "Po": 1.0 / 0.877},
    "Gr": {"Fe": -1, "S2": -2 / 3.0, "Gr": 1 / 3.0},
    "FeS2": {"Fe": -1, "S2": -1.0, "FeS2": 1.0},
}
COLORS = {"Fe": "#cfd8dc", "Po": "#ffcc80", "Gr": "#a5d6a7", "FeS2": "#90caf9"}

NOTE = (
    "Greigite field terminates at its ~500 K decomposition (kinetic).\n"
    "The sharp triangular pinch sketched by the reviewer emerges fully\n"
    "only with the continuous Fe$_{(1-x)}$S solution model (pending);\n"
    "the fixed-Fe$_{0.877}$S engine gives a narrowing wedge."
)


def make(offset, fname, title, ylo=-45.0):
    rx.set_greigite_dHf_offset(offset)
    Tg = np.linspace(300.0, 600.0, 601)
    Lg = np.linspace(ylo, 0.0, 481)
    dGf = {"Fe": np.zeros_like(Tg)}
    for p in ["Po", "Gr", "FeS2"]:
        dGf[p] = rx.reaction_dG(STOICH[p], Tg)
    TT, LL = np.meshgrid(Tg, Lg)
    Phi = np.zeros((len(PHASES), LL.shape[0], LL.shape[1]))
    for i, p in enumerate(PHASES):
        Phi[i] = dGf[p][None, :] - NU_S2[p] * RR * TT * LN10 * LL
    gi = PHASES.index("Gr")
    Phi[gi] = np.where(TT > GR_CAP, 1e30, Phi[gi])  # cap -> triangle apex
    field = np.argmin(Phi, axis=0)

    cmap = ListedColormap([COLORS[p] for p in PHASES])
    norm = BoundaryNorm(np.arange(-0.5, len(PHASES) + 0.5, 1), cmap.N)
    fig, ax = plt.subplots(figsize=(8.5, 6.2))
    ax.pcolormesh(TT, LL, field, cmap=cmap, norm=norm, shading="auto")
    ax.contour(
        TT,
        LL,
        field,
        levels=np.arange(0.5, len(PHASES) - 0.5 + 1, 1),
        colors="k",
        linewidths=1.3,
    )
    ax.axvline(GR_CAP, color="k", ls="--", lw=0.9, alpha=0.7)
    ax.text(
        GR_CAP + 3,
        ylo + 1.5,
        "Fe$_3$S$_4$ decomp.\n~%.0f K" % GR_CAP,
        fontsize=7.5,
        va="bottom",
    )
    for i, p in enumerate(PHASES):
        mask = field == i
        if mask.sum() > 0:
            mx, my = TT[mask].mean(), LL[mask].mean()
            if p == "Gr":
                ax.annotate(
                    LABELS[p],
                    xy=(mx, my),
                    xytext=(355, -3.5),
                    ha="center",
                    va="center",
                    fontsize=10,
                    fontweight="bold",
                    arrowprops=dict(arrowstyle="->", lw=1.0, color="k"),
                )
            elif p == "Fe":
                # Fe is a thin bottom sliver after the y-trim; label at bottom-mid.
                fx = 430.0
                fy = (
                    np.interp(
                        fx,
                        Tg,
                        np.minimum(
                            rx.reaction_dG(STOICH["Po"], Tg)
                            / (NU_S2["Po"] * RR * Tg * LN10),
                            0.0,
                        ),
                    )
                    - 1.2
                )
                ax.text(
                    fx,
                    max(fy, ylo + 0.8),
                    "Fe",
                    ha="center",
                    va="top",
                    fontsize=11,
                    fontweight="bold",
                )
            else:
                ax.text(
                    mx,
                    my,
                    LABELS[p],
                    ha="center",
                    va="center",
                    fontsize=11,
                    fontweight="bold",
                )
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel("log $f$(S$_2$)")
    ax.set_xlim(300, 600)
    ax.set_ylim(ylo, 0)
    ax.set_title(title, fontsize=11)
    ax.text(
        0.97,
        0.04,
        NOTE,
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
        "Sources: NIST-JANAF (Chase 1998); GS1992 pyrrhotite; "
        "Shumway 2022 greigite Cp/S; Subramani 2020 greigite $H_f$.",
        ha="center",
        fontsize=7,
        style="italic",
    )
    fig.tight_layout(rect=[0, 0.03, 1, 1])
    fig.savefig(fname, dpi=150)
    plt.close(fig)
    rx.set_greigite_dHf_offset(0.0)
    present = [PHASES[k] for k in sorted(set(field.ravel()))]
    fr = {p: round((field == i).mean(), 3) for i, p in enumerate(PHASES)}
    print("wrote", fname, "fields:", present, "area frac", fr)
    return present


print(
    "GR_CAP =",
    GR_CAP,
    "K   zero_cross =",
    zero_cross,
    "  dGdecomp(298)=%.0f J (positive => metastable, kinetic cap)" % dGdec[0],
)
make(
    +5951.0,
    str(FIG / "fig2_lower_admissible.png"),
    "Fe-S predominance, greigite d$H_f$ = $-$142.1 kJ/mol per FeS$_{1.33}$\n"
    "(lower admissible bound; offset +5951 J/Fe$_3$S$_4$)",
)
make(
    +21900.0,
    str(FIG / "fig2_upper_limit.png"),
    "Fe-S predominance, greigite d$H_f$ = $-$136.8 kJ/mol per FeS$_{1.33}$\n"
    "(upper +7.3 bound; offset +21900 J/Fe$_3$S$_4$)",
)
