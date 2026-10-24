#!/usr/bin/env python
"""Fig. S6 (B&W, single-panel pair) — Fe-S-O predominance at 600 K.

Companion to the 300 K split (make_Figure_S5_predominance_300K.py, Fig. S5A/S5B). Same single TDB, same grand-potential construction, T = 600 K.
Two outputs (the manuscript stacks them as Fig. S6A / S6B):
  Figure_S6A.png     central dHf (greigite stable, narrower field)
  Figure_S6B.png  least-stable bound, greigite raised by its reported uncertainty
The per-Fe enthalpy offset is temperature-independent, so the greigite shift (+7300 J/mol-Fe)
carries over from the 300 K case; pyrite and pyrrhotite stay at the assessed values.
No in-image title (the caption lives in the manuscript).
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
import label_fit

bw.apply()
DB = str(TDB / "fes_o_greigite_v1.tdb")
db = Database(DB)
T = 600.0
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
    "FeS2": "FeS$_2$ pyrite",  # one line: at 600 K pyrite is only a thin band
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


def make(dGr, dPo, dPy, out):
    Phi = np.zeros((len(PHASES),) + OO.shape)
    for i, p in enumerate(PHASES):
        g = Gpf[p] + {"Gr": dGr, "Po": dPo, "FeS2": dPy}.get(p, 0.0)
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
    present = []
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
    # field labels last, once the layout is final, so each one is checked clear of every line
    for i, p in enumerate(PHASES):
        if (field == i).sum() > 300:
            label_fit.place(ax, LAB[p], field, O, S, i, fontweight="bold")
            present.append(p)
    fig.savefig(out, dpi=190, bbox_inches="tight")
    print("wrote", out, "| fields:", ",".join(present))


print("T=600 K | log f(S2) native-S sat = %.2f" % LFS2_SAT)
make(0.0, 0.0, 0.0, str(FIG / "Figure_S6A.png"))
make(7300.0, 0.0, 0.0, str(FIG / "Figure_S6B.png"))


def boundary(a, b, dGr=0.0):
    ga = Gpf[a] + (dGr if a == "Gr" else 0.0)
    gb = Gpf[b] + (dGr if b == "Gr" else 0.0)
    return ((ga - gb) / (nS2[a] - nS2[b]) - G_S2) / RTLN10


# Table 2 sulfide rows at 600 K, central and least-stable greigite bound
print("central Fe/Gr          log f(S2) = %.1f" % boundary("Fe", "Gr"))
print("central Gr/FeS2        log f(S2) = %.1f" % boundary("Gr", "FeS2"))
print("least stable Fe/Po     log f(S2) = %.1f" % boundary("Fe", "Po", 7300.0))
print("least stable Po/Gr     log f(S2) = %.1f" % boundary("Po", "Gr", 7300.0))
print("least stable Gr/FeS2   log f(S2) = %.1f" % boundary("Gr", "FeS2", 7300.0))
