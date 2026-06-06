import json, numpy as np, sys
from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent  # = manuscript/
ROOT = HERE.parent  # = repo root
FIG = ROOT / "artifacts" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
ENGINE = ROOT / "engine"
sys.path.insert(0, str(HERE))
import pyrrhotite_wp as wp, greigite_cp as g

OUT = str(FIG)
R = 8.314462618
LN10 = np.log(10)

# ---- ENGINE band (Dilner basis) from the JSON pycalphad wrote ----
J = json.load(open(ROOT / "artifacts" / "fes_engine_boundaries.json"))
gf = J["central"]["fields"]["GREIGITE"]
Te = np.array(sorted(int(t) for t in gf))
e_lo = np.array([gf[str(t)][0] for t in Te])  # Fe/greigite boundary
e_hi = np.array([gf[str(t)][1] for t in Te])  # greigite/pyrite boundary


# ---- PROTOTYPE band (W-P/JANAF) recomputed by grand-potential, offset 0 ----
def cp_py(TT):
    return 72.387 + 0.0088501 * TT + 7.3e-10 * TT**2 - 1.14279e6 / TT**2


def Gpy(T):
    gr = np.linspace(298.15, T, 500)
    return (
        -171048.03
        + np.trapezoid(cp_py(gr), gr)
        - T * (52.93 + np.trapezoid(cp_py(gr) / gr, gr))
    )


YFE = np.linspace(0.80, 0.99999, 160)
yv = 1 - YFE
XS = np.linspace(-70, 2, 900)
Tp = np.linspace(300, 600, 16)
p_lo = []
p_hi = []
for T in Tp:
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
    ggr = (-144100.0) + Hg - T * Sg
    seq = np.argmin(
        np.vstack([np.zeros_like(XS), om_po, gpy - 2 * muS, ggr - 1.33 * muS]), axis=0
    )
    idx = np.where(seq == 3)[0]
    p_lo.append(XS[idx[0]] if len(idx) else np.nan)
    p_hi.append(XS[idx[-1]] if len(idx) else np.nan)
p_lo = np.array(p_lo)
p_hi = np.array(p_hi)

fig, (axE, axP) = plt.subplots(1, 2, figsize=(13.4, 6.2), sharey=True)
GREEN = "#5fae77"
BLUE = "#bcd4ec"
GREY = "#d9d9d9"


def panel(ax, T, lo, hi, title, ylo=-60, yhi=2):
    ax.fill_between(T, np.full_like(T, ylo), lo, color=GREY, alpha=0.7, lw=0)  # Fe
    ax.fill_between(T, lo, hi, color=GREEN, alpha=0.88, lw=0)  # greigite
    ax.fill_between(T, hi, np.full_like(T, yhi), color=BLUE, alpha=0.6, lw=0)  # pyrite
    ax.plot(T, lo, "k-", lw=1.2)
    ax.plot(T, hi, "k-", lw=1.2)
    ax.text(
        T[len(T) // 3], (ylo + lo[len(T) // 3]) / 2 - 3, "Fe", fontsize=12, ha="center"
    )
    ax.text(
        T[len(T) // 2],
        (lo[len(T) // 2] + hi[len(T) // 2]) / 2,
        "Fe$_3$S$_4$ greigite",
        fontsize=13,
        weight="bold",
        color="#16431f",
        ha="center",
        va="center",
    )
    ax.text(
        T[len(T) // 2],
        (hi[len(T) // 2] + yhi) / 2 + 2,
        "FeS$_2$ pyrite",
        fontsize=12,
        ha="center",
    )
    ax.set_title(title, fontsize=10.5)
    ax.set_xlabel("Temperature, K", fontsize=11)
    ax.set_xlim(300, 600)
    ax.set_ylim(ylo, yhi)


panel(
    axE,
    Te,
    e_lo,
    e_hi,
    "ENGINE — pycalphad on Dilner 2015 + measured greigite\n(real CALPHAD solver; published Fe–S basis)",
)
panel(
    axP,
    Tp,
    p_lo,
    p_hi,
    "PROTOTYPE — numpy grand-potential\n(Waldner–Pelton pyrrhotite + JANAF/greigite)",
)
axE.set_ylabel(r"log $f$(S$_2$), bar", fontsize=11)
axE.legend(
    handles=[
        Patch(facecolor=GREY, label="Fe (bcc)"),
        Patch(facecolor=GREEN, label="Fe$_3$S$_4$ greigite"),
        Patch(facecolor=BLUE, label="FeS$_2$ pyrite"),
    ],
    loc="lower right",
    fontsize=8.5,
    framealpha=0.95,
)
fig.suptitle(
    "Fe–S: independent confirmation — both methods give greigite PREEMPTING pyrrhotite (no Fe$_{1-x}$S field)",
    fontsize=11.5,
    y=0.99,
)
fig.text(
    0.5,
    0.005,
    "Same primary data for greigite (Subramani 2020 ΔH$_f$ + Shumway 2022 S°/Cp). Sulfur treated as S$_2$ gas at all T. Boundaries agree within ~1–2 log units; "
    "the small offset is the different pyrrhotite/pyrite assessment (Dilner/Lee vs Waldner–Pelton/JANAF), not a logic difference.",
    ha="center",
    fontsize=7.3,
    color="0.35",
)
fig.tight_layout(rect=[0, 0.02, 1, 0.97])
fig.savefig(f"{OUT}/fig_engine_vs_prototype.png", dpi=200)
print("saved fig_engine_vs_prototype.png")
print(
    "engine greigite @300:",
    round(e_lo[0], 1),
    round(e_hi[0], 1),
    " @600:",
    round(e_lo[-1], 1),
    round(e_hi[-1], 1),
)
print(
    "proto  greigite @300:",
    round(p_lo[0], 1),
    round(p_hi[0], 1),
    " @600:",
    round(p_lo[-1], 1),
    round(p_hi[-1], 1),
)
