#!/usr/bin/env python
"""Fig. 5 (A/B, B&W, single-panel) — Fe-S-O predominance at 298.15 K, split.

Fig. 5A: central dHf (greigite stable, borders the oxides).
Fig. 5B: greigite +1 sigma (pyrrhotite appears; pyrite borders hematite).

Single-database engine (one TDB). White fields, black boundary lines, in-field
medoid labels, native-S cap. No footnote on the figure (text goes in caption).
"""

import sys
from pathlib import Path
import numpy as np, warnings

warnings.filterwarnings("ignore")
import symengine as se
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pycalphad import Database, calculate
from pycalphad.variables import T as TV, P as PV

HERE = Path(__file__).resolve().parent  # = manuscript/
ROOT = HERE.parent  # = repo root
TDB = ROOT / "artifacts" / "tdb"
FIG = ROOT / "artifacts" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE))
import bw_style as bw

bw.apply()
DB = str(TDB / "fes_o_greigite_v1.tdb")
db = Database(DB)
T = 298.15
RTLN10 = 8.31451 * T * np.log(10.0)


def evf(name):
    e = db.symbols[name]
    for _ in range(60):
        ne = e.xreplace({se.Symbol(k): val for k, val in db.symbols.items()})
        if ne == e:
            break
        e = ne
    return float(e.xreplace({TV: T, PV: 101325.0}))


G_S2 = evf("F16023T")
G_O2 = evf("F14375T")
GHSERSS = evf("GHSERSS")
LFS2_SAT = (2 * GHSERSS - G_S2) / RTLN10


def minG_perFe(phase, target, atoms_per_Fe, tol=0.02):
    r = calculate(db, ["FE", "O", "S", "VA"], phase, T=T, P=101325, pdens=2000)
    X = {c: np.ravel(r.X.sel(component=c).values) for c in ("FE", "O", "S")}
    gm = np.ravel(r.GM.values)
    m = np.ones(len(gm), bool)
    for c, xt in target.items():
        m &= np.abs(X[c] - xt) < tol
    return gm[m].min() * atoms_per_Fe


SPEC = {
    "Fe": ("BCC_A2", {"FE": 1.0}, 1.0, 0.0, 0.0),
    "Po": ("PYRRHOTITE", {"FE": 0.5, "S": 0.5}, 2.0, 0.5, 0.0),
    "Gr": ("GREIGITE", {"FE": 3 / 7, "S": 4 / 7}, 7 / 3, 2 / 3, 0.0),
    "FeS2": ("PYRITE", {"FE": 1 / 3, "S": 2 / 3}, 3.0, 1.0, 0.0),
    "Fe3O4": ("SPINEL", {"FE": 3 / 7, "O": 4 / 7}, 7 / 3, 0.0, 2 / 3),
    "Fe2O3": ("CORUNDUM", {"FE": 2 / 5, "O": 3 / 5}, 5 / 2, 0.0, 3 / 4),
    "FeSO4": ("ANHYDRITE", {"FE": 1 / 6, "S": 1 / 6, "O": 4 / 6}, 6.0, 0.5, 2.0),
    "Fe2(SO4)3": (
        "FE2O12S3",
        {"FE": 2 / 17, "S": 3 / 17, "O": 12 / 17},
        17 / 2,
        0.75,
        3.0,
    ),
}
Gpf, nS2, nO2 = {}, {}, {}
for k, (ph, tgt, apf, ns, no) in SPEC.items():
    Gpf[k] = minG_perFe(ph, tgt, apf)
    nS2[k] = ns
    nO2[k] = no

PHASES = ["Fe", "Po", "Gr", "FeS2", "Fe3O4", "Fe2O3", "FeSO4", "Fe2(SO4)3"]
LAB = {
    "Fe": "Fe",
    "Po": "Fe$_{1-x}$S\npyrrhotite",
    "Gr": "Fe$_3$S$_4$\ngreigite",
    "FeS2": "FeS$_2$\npyrite",
    "Fe3O4": "Fe$_3$O$_4$\nmagnetite",
    "Fe2O3": "Fe$_2$O$_3$\nhematite",
    "FeSO4": "FeSO$_4$",
    "Fe2(SO4)3": "Fe$_2$(SO$_4$)$_3$",
}

O = np.linspace(-120.0, 0.0, 800)
S = np.linspace(-60.0, 0.0, 800)
OO, SS = np.meshgrid(O, S)
muS2 = G_S2 + RTLN10 * SS
muO2 = G_O2 + RTLN10 * OO


def make(dGr, dPo, title, out):
    Phi = np.zeros((len(PHASES),) + OO.shape)
    for i, p in enumerate(PHASES):
        g = Gpf[p] + (dGr if p == "Gr" else dPo if p == "Po" else 0.0)
        Phi[i] = g - nS2[p] * muS2 - nO2[p] * muO2
    field = np.argmin(Phi, axis=0)
    field = np.where(SS > LFS2_SAT, len(PHASES), field)
    fig, ax = plt.subplots(figsize=(8.2, 7.0))
    ax.contour(
        OO,
        SS,
        field,
        levels=np.arange(0.5, len(PHASES) + 0.5, 1),
        colors="black",
        linewidths=1.1,
    )
    for i, p in enumerate(PHASES):
        mk = field == i
        if mk.sum() > 300:
            cx, cy = OO[mk].mean(), SS[mk].mean()
            j = np.argmin((OO[mk] - cx) ** 2 + (SS[mk] - cy) ** 2)
            ax.text(
                OO[mk][j],
                SS[mk][j],
                LAB[p],
                ha="center",
                va="center",
                fontsize=11.5,
                fontweight="bold",
            )
    ax.axhline(LFS2_SAT, color="black", ls="--", lw=1.4)
    ax.text(
        -60,
        LFS2_SAT / 2.0,
        "native S  (S$_2$ saturation, log $f$(S$_2$) = %.2f)" % LFS2_SAT,
        ha="center",
        va="center",
        fontsize=10,
        fontweight="bold",
    )
    ax.set_xlabel("log $f$(O$_2$), bar")
    ax.set_ylabel("log $f$(S$_2$), bar")
    ax.set_xlim(-120, 0)
    ax.set_ylim(-60, 0)
    fig.tight_layout()
    fig.savefig(out, dpi=190, bbox_inches="tight")
    print("wrote", out)


make(
    0.0,
    0.0,
    "Fig. 7 — Fe–S–O predominance, 298.15 K, central ΔH$_f$ (greigite stable)",
    str(FIG / "Figure_5A.png"),
)
make(
    7300.0,
    -3500.0,
    "Fig. 8 — Fe–S–O predominance, 298.15 K, greigite +1σ (pyrrhotite appears)",
    str(FIG / "Figure_5B.png"),
)
