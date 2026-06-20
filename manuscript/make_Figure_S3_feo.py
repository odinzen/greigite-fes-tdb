#!/usr/bin/env python
"""Fig. S3 (B&W) — COMPUTED Fe-O predominance vs temperature, 300-600 K.

Single-database engine. The Fe-O stability diagram (log f(O2) vs T): metallic
iron, the Fe(II,III) spinel magnetite (Fe3O4), and hematite (Fe2O3). This is the
oxide companion to the Fe-S diagram (manuscript Fig. 3) — both computed from the
same assembled database, so greigite (Fe3S4) is the sulfide analogue of magnetite.
White fields, black boundaries, medoid labels, no in-image title.
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
OUT = str(FIG / "Figure_S3.png")
db = Database(DB)


def evf(name, T):
    e = db.symbols[name]
    for _ in range(60):
        ne = e.xreplace({se.Symbol(k): val for k, val in db.symbols.items()})
        if ne == e:
            break
        e = ne
    return float(e.xreplace({TV: T, PV: 101325.0}))


def minG_perFe(phase, target, apf, T, tol=0.02):
    r = calculate(db, ["FE", "O", "S", "VA"], phase, T=T, P=101325, pdens=2000)
    X = {c: np.ravel(r.X.sel(component=c).values) for c in ("FE", "O", "S")}
    gm = np.ravel(r.GM.values)
    m = np.ones(len(gm), bool)
    for c, xt in target.items():
        m &= np.abs(X[c] - xt) < tol
    return gm[m].min() * apf


# Fe-O phases only: metal, magnetite (spinel), hematite (corundum).
SPEC = {
    "Fe": ("BCC_A2", {"FE": 1.0}, 1.0),
    "Fe3O4": ("SPINEL", {"FE": 3 / 7, "O": 4 / 7}, 7 / 3),
    "Fe2O3": ("CORUNDUM", {"FE": 2 / 5, "O": 3 / 5}, 5 / 2),
}
nO2 = {"Fe": 0.0, "Fe3O4": 2 / 3, "Fe2O3": 3 / 4}

Tg = np.array([300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0])
feo_lo, feo_hi = [], []
for T in Tg:
    RTLN10 = 8.31451 * T * np.log(10.0)
    G_O2 = evf("F14375T", T)
    Gpf = {k: minG_perFe(s[0], s[1], s[2], T) for k, s in SPEC.items()}
    feo_lo.append(
        ((Gpf["Fe"] - Gpf["Fe3O4"]) / (nO2["Fe"] - nO2["Fe3O4"]) - G_O2) / RTLN10
    )
    feo_hi.append(
        ((Gpf["Fe3O4"] - Gpf["Fe2O3"]) / (nO2["Fe3O4"] - nO2["Fe2O3"]) - G_O2) / RTLN10
    )
feo_lo, feo_hi = np.array(feo_lo), np.array(feo_hi)

Tf = np.linspace(300, 600, 400)
YLO, YHI = -95, -18
labels = {0: "Fe", 1: "Fe$_3$O$_4$\nmagnetite", 2: "Fe$_2$O$_3$\nhematite"}

fig, ax = plt.subplots(figsize=(8.2, 6.0))
L = np.polyval(np.polyfit(Tg, feo_lo, 3), Tf)
H = np.polyval(np.polyfit(Tg, feo_hi, 3), Tf)
ax.plot(Tf, L, "k", lw=1.8)
ax.plot(Tf, H, "k", lw=1.8)
ax.set_xlim(300, 600)
ax.set_ylim(YLO, YHI)
ax.set_xlabel("Temperature, K")
ax.set_ylabel("log $f$(O$_2$), bar")

# medoid label per region (0 below L, 1 between L and H, 2 above H)
Y = np.linspace(YLO, YHI, 400)[:, None]
reg = np.where(Y < L[None, :], 0, np.where(Y < H[None, :], 1, 2))
TT = np.broadcast_to(Tf[None, :], reg.shape)
YT = np.broadcast_to(Y, reg.shape)
for idx, lab in labels.items():
    m = reg == idx
    if m.sum() > 50:
        cx, cy = TT[m].mean(), YT[m].mean()
        d = ((TT[m] - cx) / 300.0) ** 2 + ((YT[m] - cy) / (YHI - YLO)) ** 2
        j = np.argmin(d)
        ax.text(
            TT[m][j],
            YT[m][j],
            lab,
            ha="center",
            va="center",
            fontsize=12,
            fontweight="bold",
        )

fig.tight_layout()
fig.savefig(OUT, dpi=190, bbox_inches="tight")
print("wrote", OUT)
print(
    "Fe/Fe3O4 %.1f..%.1f ; M-H %.1f..%.1f"
    % (feo_lo[0], feo_lo[-1], feo_hi[0], feo_hi[-1])
)
