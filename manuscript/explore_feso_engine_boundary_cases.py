#!/usr/bin/env python
"""Fig. 3 (Fe-S-O) — TRUE single-database engine predominance, two boundary cases.

All phase Gibbs energies come from ONE database (fes_o_greigite_v1.tdb =
Dilner & Selleby 2017 Ca-Fe-O-S, deduped, + grafted Dilner-2015 pyrite +
measured greigite). NO JANAF/CRC hybrid. Per-mol-Fe phase energies are taken
as the minimum GM (pycalphad `calculate`, robust vs the line-compound nan in
`equilibrium`) at each phase's stoichiometric composition; gas O2/S2 references
are the TDB GAS functions. Validated: oxide formation energies match literature
to ~1% (FeO -251.9, Fe3O4 -334.7/Fe, Fe2O3 -367.3 kJ).

Predominance = grand-potential argmin per mol-Fe over the (log fO2, log fS2)
plane at 298.15 K:
    Omega(phase) = G_perFe - nu_S2*mu_S2 - nu_O2*mu_O2
    mu_S2 = G_S2ref + RT ln10 logfS2 ;  mu_O2 = G_O2ref + RT ln10 logfO2

Two boundary cases mirror Fig. 2b (greigite +-1sigma, pyrrhotite -+1sigma per
mol-Fe; pyrite central): LOWER pyrrhotite eliminated, UPPER max pyrrhotite.
Native-S saturation (engine: 2 GHSERSS vs S2 ref = log fS2 -13.96 at 298 K)
restored on Fe-S-O ONLY, as a cap line.
"""

from pathlib import Path
import numpy as np, warnings

warnings.filterwarnings("ignore")
import symengine as se
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from pycalphad import Database, calculate
from pycalphad.variables import T as TV, P as PV

HERE = Path(__file__).resolve().parent  # = manuscript/
ROOT = HERE.parent  # = repo root
TDB = ROOT / "artifacts" / "tdb"
FIG = ROOT / "artifacts" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
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
LFS2_SAT = (2 * GHSERSS - G_S2) / RTLN10  # native-S saturation (engine)


def minG_perFe(phase, target, atoms_per_Fe, tol=0.02):
    r = calculate(db, ["FE", "O", "S", "VA"], phase, T=T, P=101325, pdens=2000)
    X = {c: np.ravel(r.X.sel(component=c).values) for c in ("FE", "O", "S")}
    gm = np.ravel(r.GM.values)
    m = np.ones(len(gm), bool)
    for c, xt in target.items():
        m &= np.abs(X[c] - xt) < tol
    return gm[m].min() * atoms_per_Fe


# phase: (calc-phase, target comp, atoms/Fe, nu_S2, nu_O2)
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
    print(f"  {k:10s} G/molFe={Gpf[k]:11.1f}  nS2={ns:.3f} nO2={no:.3f}")
print(f"  gas refs S2={G_S2:.1f} O2={G_O2:.1f} ; native-S sat log fS2={LFS2_SAT:.2f}")

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
COL = {
    "Fe": "#cfd8dc",
    "Po": "#ffcc80",
    "Gr": "#a5d6a7",
    "FeS2": "#90caf9",
    "Fe3O4": "#b39ddb",
    "Fe2O3": "#ef9a9a",
    "FeSO4": "#ffe082",
    "Fe2(SO4)3": "#fff0b3",
}

O = np.linspace(-120.0, 0.0, 800)
S = np.linspace(-60.0, 0.0, 800)
OO, SS = np.meshgrid(O, S)
muS2 = G_S2 + RTLN10 * SS
muO2 = G_O2 + RTLN10 * OO

CASES = [
    (
        0.0,
        0.0,
        "CENTRAL  ΔH$_f$ = −144.1 kJ/mol-FeS$_{1.33}$\ngreigite stable — greigite borders the oxides",
    ),
    (
        +7300.0,
        -3500.0,
        "UPPER bound  (greigite +1σ)\npyrrhotite appears; pyrite borders hematite",
    ),
]

plt.rcParams.update({"font.size": 14})
fig, axes = plt.subplots(1, 2, figsize=(17, 9), sharey=True)
for ax, (dGr, dPo, title) in zip(axes, CASES):
    Phi = np.zeros((len(PHASES),) + OO.shape)
    for i, p in enumerate(PHASES):
        g = Gpf[p] + (dGr if p == "Gr" else dPo if p == "Po" else 0.0)
        Phi[i] = g - nS2[p] * muS2 - nO2[p] * muO2
    field = np.argmin(Phi, axis=0)
    field = np.where(
        SS > LFS2_SAT, len(PHASES), field
    )  # native-S field above saturation
    PLOT_COL = [COL[p] for p in PHASES] + ["#efe3b0"]
    cmap = ListedColormap(PLOT_COL)
    norm = BoundaryNorm(np.arange(-0.5, len(PLOT_COL) + 0.5, 1), cmap.N)
    ax.pcolormesh(OO, SS, field, cmap=cmap, norm=norm, shading="auto")
    ax.contour(
        OO,
        SS,
        field,
        levels=np.arange(0.5, len(PHASES) + 0.5, 1),
        colors="k",
        linewidths=1.0,
    )
    for i, p in enumerate(PHASES):
        mk = field == i
        if mk.sum() > 300:
            cx, cy = OO[mk].mean(), SS[mk].mean()
            j = np.argmin(
                (OO[mk] - cx) ** 2 + (SS[mk] - cy) ** 2
            )  # in-field point nearest centroid
            ax.text(
                OO[mk][j],
                SS[mk][j],
                LAB[p],
                ha="center",
                va="center",
                fontsize=13,
                fontweight="bold",
            )
    adj = set()
    for a, b in [(field[:, :-1], field[:, 1:]), (field[:-1, :], field[1:, :])]:
        m = a != b
        for u, v in zip(a[m], b[m]):
            adj.add(tuple(sorted((int(u), int(v)))))
    NM = PHASES + ["nativeS"]
    print(
        title.split(chr(10))[0], "adjacencies:", sorted((NM[u], NM[v]) for u, v in adj)
    )
    ax.axhline(LFS2_SAT, color="#7a4a00", ls="--", lw=2.0)
    ax.text(
        -60,
        LFS2_SAT / 2.0,
        "native S\n(S$_2$ saturation, log $f$(S$_2$) = %.2f)" % LFS2_SAT,
        ha="center",
        va="center",
        fontsize=12,
        color="#5a3600",
        fontweight="bold",
    )
    ax.set_title(title, fontsize=14)
    ax.set_xlabel("log $f$(O$_2$), bar", fontsize=15)
    ax.tick_params(labelsize=13)
    ax.set_xlim(-120, 0)
    ax.set_ylim(-60, 0)
axes[0].set_ylabel("log $f$(S$_2$), bar", fontsize=16)
fig.suptitle(
    "Fig. 6 — Fe–S–O predominance at 298.15 K (single-database engine): central ΔH$_f$ (left) vs greigite +1σ (right)\n"
    "Dilner-2017 Ca–Fe–O–S (deduped) + grafted pyrite + measured greigite",
    fontsize=14,
    fontweight="bold",
)
note = (
    "Single database (no JANAF/CRC hybrid); oxides validated to ~1% vs literature; capped at native-S saturation "
    f"(log f(S$_2$)={LFS2_SAT:.2f}).  At central ΔH$_f$ greigite is stable and borders the oxides, so pyrite does NOT border "
    "hematite (greigite intervenes); the pyrite–hematite contact returns only at greigite +1σ (right) and in the greigite-free control (Fig. 4)."
)
fig.text(0.5, 0.02, note, ha="center", va="bottom", fontsize=9, color="#222")
fig.tight_layout(rect=[0, 0.055, 1, 0.92])
fig.savefig(str(FIG / "explore_feso_engine_boundary_cases.png"), dpi=170)
print("wrote explore_feso_engine_boundary_cases.png")
