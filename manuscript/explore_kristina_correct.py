import sys
from pathlib import Path
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TDB = ROOT / "artifacts" / "tdb"
FIG = ROOT / "artifacts" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE))
import pyrrhotite_wp as wp, greigite_cp as g

R = 8.314462618
LN10 = np.log(10)


def cp_py(TT):
    return 72.387 + 0.0088501 * TT + 7.3e-10 * TT**2 - 1.14279e6 / TT**2


def Gpy(T):
    gr = np.linspace(298.15, T, 800)
    return (
        -171048.03
        + np.trapezoid(cp_py(gr), gr)
        - T * (52.93 + np.trapezoid(cp_py(gr) / gr, gr))
    )


YFE = np.linspace(0.80, 0.99999, 160)
yv = 1.0 - YFE
XS = np.linspace(-65, 2, 600)


def fields(T, offs):
    gfes = wp.GFeS(T)
    gvas = wp.GVaS(T)
    Lv = wp.L_param(T)
    gpy = Gpy(T)
    G_S2 = wp.G_S2_gas(T)
    Gpo = (
        YFE * gfes
        + yv * gvas
        + R * T * (YFE * np.log(YFE) + yv * np.log(yv))
        + YFE * yv * Lv
    )
    muS = 0.5 * (G_S2 + R * T * LN10 * XS)
    om_po = np.min((Gpo[:, None] - muS[None, :]) / YFE[:, None], axis=0)
    Hg = g.enthalpy_increment(T)
    Sg = g.entropy(T)
    res = {}
    for off in offs:
        ggr = (-144100.0 + off) + Hg - T * Sg
        OM = np.vstack([np.zeros_like(XS), om_po, ggr - 1.33 * muS, gpy - 2.0 * muS])
        seq = np.argmin(OM, axis=0)

        def rng(p):
            idx = np.where(seq == p)[0]
            return (XS[idx[0]], XS[idx[-1]]) if len(idx) else (np.nan, np.nan)

        res[off] = (rng(2), rng(1))
    return res


T = np.linspace(300, 580, 26)
G = [fields(t, (0.0, 7300.0, -7300.0)) for t in T]
grc = np.array([x[0.0][0] for x in G])
gru = np.array([x[7300.0][0] for x in G])
grl = np.array([x[-7300.0][0] for x in G])
pou = np.array([x[7300.0][1] for x in G])
fig, ax = plt.subplots(figsize=(8.2, 6.4))
ax.fill_between(T, grl[:, 0], gru[:, 0], color="#a9d3b5", alpha=0.55, lw=0)
ax.fill_between(T, gru[:, 1], grl[:, 1], color="#a9d3b5", alpha=0.55, lw=0)
ax.fill_between(T, grc[:, 0], grc[:, 1], color="#5fae77", alpha=0.85, lw=0)
ax.plot(T, grc[:, 0], "k-", lw=1.4)
ax.plot(T, grc[:, 1], "k-", lw=1.4)
m = np.isfinite(pou[:, 0])
if m.any():
    ax.fill_between(
        T[m],
        pou[m, 0],
        pou[m, 1],
        facecolor="none",
        hatch="xxx",
        edgecolor="#b06a2e",
        lw=0,
    )
ax.text(330, -12, "FeS$_2$ pyrite", fontsize=12)
ax.text(360, -37, "Fe$_3$S$_4$ greigite", fontsize=13, weight="bold", color="#16431f")
ax.text(450, -56, "Fe", fontsize=12)
if m.any():
    ax.annotate(
        "pyrrhotite Fe$_{1-x}$S reappears only at the\n+7.3 kJ greigite bound, T > ~500 K",
        xy=(T[m][-1], pou[m, 0][-1]),
        xytext=(352, -52),
        fontsize=8,
        color="#8a5a23",
        arrowprops=dict(arrowstyle="->", color="#b06a2e", lw=0.8),
    )
ax.set_xlabel("Temperature, K", fontsize=12)
ax.set_ylabel(r"log $f$(S$_2$), bar", fontsize=12)
ax.set_xlim(300, 580)
ax.set_ylim(-60, 0)
ax.set_title(
    "Fe–S: continuous Waldner–Pelton pyrrhotite + greigite (equilibrium hull)\n"
    "at the measured greigite ΔH$_f$, greigite PREEMPTS pyrrhotite; Fe$_{1-x}$S is within the ±7.3 kJ band",
    fontsize=10,
)
ax.legend(
    handles=[
        Patch(facecolor="#5fae77", label="Fe$_3$S$_4$ field (central ΔH$_f$)"),
        Patch(facecolor="#a9d3b5", label="greigite boundary ±7.3 kJ/mol"),
        Patch(
            facecolor="white",
            hatch="xxx",
            edgecolor="#b06a2e",
            label="Fe$_{1-x}$S, only within +7.3 kJ",
        ),
    ],
    loc="lower right",
    fontsize=8.5,
    framealpha=0.95,
)
ax.text(
    0.5,
    -0.13,
    "Grand-potential minimisation over Fe + Waldner–Pelton 2005 Fe$_{1-x}$S (validated vs their Fig 9) + greigite + pyrite. "
    "Greigite is metastable; its place on the equilibrium hull is a question for review.",
    transform=ax.transAxes,
    ha="center",
    fontsize=7.2,
    color="0.4",
)
fig.tight_layout()
fig.savefig(str(FIG / "explore_kristina_correct.png"), dpi=200)
print("OK central300:", tuple(np.round(grc[0], 1)), "po@+7.3 appears:", bool(m.any()))
