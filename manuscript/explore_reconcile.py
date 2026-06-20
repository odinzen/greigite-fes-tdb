import numpy as np, sys
from pathlib import Path
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
    gr = np.linspace(298.15, T, 400)
    return (
        -171048.03
        + np.trapezoid(cp_py(gr), gr)
        - T * (52.93 + np.trapezoid(cp_py(gr) / gr, gr))
    )


YFE = np.linspace(0.80, 0.99999, 160)
yv = 1.0 - YFE
XS = np.linspace(-70, 2, 800)


def rng(seq, p):
    idx = np.where(seq == p)[0]
    return (XS[idx[0]], XS[idx[-1]]) if len(idx) else (np.nan, np.nan)


T = np.linspace(300, 560, 24)
OFF_TRI = 20000.0
A_gr = []
B_fe = []
B_po = []
B_py = []
B_gr = []
for t in T:
    gfes = wp.GFeS(t)
    gvas = wp.GVaS(t)
    Lv = wp.L_param(t)
    gpy = Gpy(t)
    G_S2 = wp.G_S2_gas(t)
    Gpo = (
        YFE * gfes
        + yv * gvas
        + R * t * (YFE * np.log(YFE) + yv * np.log(yv))
        + YFE * yv * Lv
    )
    muS = 0.5 * (G_S2 + R * t * LN10 * XS)
    om_po = np.min((Gpo[:, None] - muS[None, :]) / YFE[:, None], axis=0)
    om_fe = np.zeros_like(XS)
    om_py = gpy - 2.0 * muS
    Hg = g.enthalpy_increment(t)
    Sg = g.entropy(t)
    # A: off 0
    ggrA = (-144100.0) + Hg - t * Sg
    om_grA = ggrA - 1.33 * muS
    seqA = np.argmin(np.vstack([om_fe, om_po, om_py, om_grA]), axis=0)
    A_gr.append(rng(seqA, 3))
    # B: off +20kJ
    ggrB = (-144100.0 + OFF_TRI) + Hg - t * Sg
    om_grB = ggrB - 1.33 * muS
    seqB = np.argmin(np.vstack([om_fe, om_po, om_py, om_grB]), axis=0)
    B_fe.append(rng(seqB, 0))
    B_po.append(rng(seqB, 1))
    B_py.append(rng(seqB, 2))
    B_gr.append(rng(seqB, 3))
A_gr = np.array(A_gr)
B_fe = np.array(B_fe)
B_po = np.array(B_po)
B_py = np.array(B_py)
B_gr = np.array(B_gr)

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.8, 6.4), sharey=True)
axA.fill_between(T, A_gr[:, 0], A_gr[:, 1], color="#5fae77", alpha=0.85, lw=0)
axA.plot(T, A_gr[:, 0], "k-", lw=1.3)
axA.plot(T, A_gr[:, 1], "k-", lw=1.3)
axA.text(330, -12, "FeS$_2$ pyrite", fontsize=11)
axA.text(352, -37, "Fe$_3$S$_4$ greigite", fontsize=12, weight="bold", color="#16431f")
axA.text(460, -56, "Fe", fontsize=11)
axA.set_title(
    "A.  MEASURED greigite ΔH$_f$ (Subramani 2020, equilibrium hull)\n"
    "greigite preempts pyrrhotite — a WIDE field, no Fe$_{1-x}$S",
    fontsize=9.5,
)
axA.set_xlabel("Temperature, K", fontsize=11)
axA.set_ylabel(r"log $f$(S$_2$), bar", fontsize=11)

axB.fill_between(T, B_py[:, 0], B_py[:, 1], color="#bcd4ec", alpha=0.6, lw=0)
axB.fill_between(T, B_po[:, 0], B_po[:, 1], color="#e8d6a0", alpha=0.78, lw=0)
axB.fill_between(T, B_gr[:, 0], B_gr[:, 1], color="#5fae77", alpha=0.85, lw=0)
axB.fill_between(T, B_fe[:, 0], B_fe[:, 1], color="#d9d9d9", alpha=0.7, lw=0)
for b in (B_gr[:, 0], B_gr[:, 1], B_po[:, 0]):
    axB.plot(T, b, "k-", lw=1.0)
axB.text(330, -9, "FeS$_2$ pyrite", fontsize=11)
axB.annotate(
    "Fe$_3$S$_4$ greigite wedge",
    xy=(360, -30),
    xytext=(388, -19),
    fontsize=8.5,
    color="#16431f",
    arrowprops=dict(arrowstyle="->", color="#16431f", lw=0.8),
)
axB.text(425, -44, "Fe$_{1-x}$S pyrrhotite", fontsize=11)
axB.text(470, -58, "Fe", fontsize=11)
axB.set_title(
    "B.  greigite ΔH$_f$ set ~+20 kJ/mol-FeS$_{1.33}$ LESS stable (≈3σ)\n"
    "thin triangle with pyrrhotite below — your expected shape",
    fontsize=9.5,
)
axB.set_xlabel("Temperature, K", fontsize=11)
for ax in (axA, axB):
    ax.set_xlim(300, 560)
    ax.set_ylim(-62, 2)
fig.suptitle(
    "Fe–S: the expected triangle-with-pyrrhotite is NOT reproducible at the measured greigite enthalpy "
    "(continuous Waldner–Pelton pyrrhotite, validated vs their Fig 9)",
    fontsize=10.3,
    y=1.0,
)
fig.text(
    0.5,
    0.006,
    "Same engine, same data — only the greigite ΔH$_f$ differs. At the measured value (A) greigite is too stable to leave a pyrrhotite field. "
    "The triangle (B) needs greigite ≈3σ less stable than Subramani 2020 — OUTSIDE the published ±7.3 kJ/mol-FeS$_{1.33}$ band. "
    "The honest equilibrium answer is (A); (B) only locates the assumption your sketch implies.",
    ha="center",
    fontsize=7.2,
    color="0.35",
)
fig.tight_layout(rect=[0, 0.02, 1, 0.97])
fig.savefig(str(FIG / "explore_reconcile.png"), dpi=200)
print("saved explore_reconcile.png")
print("A_gr@300:", tuple(np.round(A_gr[0], 1)))
print(
    "B@300 py:",
    tuple(np.round(B_py[0], 1)),
    "gr:",
    tuple(np.round(B_gr[0], 1)),
    "po:",
    tuple(np.round(B_po[0], 1)),
    "fe:",
    tuple(np.round(B_fe[0], 1)),
)
