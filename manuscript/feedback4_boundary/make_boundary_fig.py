import numpy as np, warnings, math, json
from pathlib import Path

warnings.filterwarnings("ignore")
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from pycalphad import Database, equilibrium, variables as v
from pycalphad.variables import T as T_VAR

HERE = Path(__file__).resolve().parent
MAN = HERE.parent
ROOT = MAN.parent
TDB = ROOT / "artifacts" / "tdb"
FIG = ROOT / "artifacts" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

R = 8.31451
LN10 = math.log(10.0)
COMPS = ["FE", "S", "VA"]
CONDENSED = ["BCC_A2", "FCC_A1", "PYRRHOTITE", "PYRITE", "GREIGITE", "ORTHORHOMBIC_S"]
T_GRID = np.array([300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0])
XS = np.linspace(0.40, 0.66, 46)
COL = {
    "BCC_A2": "#9aa0a6",
    "PYRRHOTITE": "#d4774e",
    "GREIGITE": "#4e79a7",
    "PYRITE": "#e6c34a",
    "ORTHORHOMBIC_S": "#b07aa1",
}
LAB = {
    "BCC_A2": "Fe (bcc)",
    "PYRRHOTITE": "pyrrhotite Fe$_{1-x}$S",
    "GREIGITE": "greigite Fe$_3$S$_4$",
    "PYRITE": "pyrite FeS$_2$",
    "ORTHORHOMBIC_S": "S",
}


def gref(db, T):
    return float(db.symbols["F15281T"].subs({T_VAR: float(T)}))


def fields(tdb):
    db = Database(tdb)
    out = {}
    pocomp = {}
    for T in T_GRID:
        g = gref(db, T)
        for xs in XS:
            eq = equilibrium(
                db,
                COMPS,
                CONDENSED,
                {v.X("S"): float(xs), v.T: float(T), v.P: 101325, v.N: 1},
            )
            ph = np.ravel(eq.Phase.values)
            Xs = eq.X.sel(component="S").values.ravel()
            muS = float(np.ravel(eq.MU.sel(component="S").values)[0])
            lf = (2 * muS - g) / (R * T * LN10)
            for p, x in zip(ph, Xs):
                if not p:
                    continue
                out.setdefault(p, {}).setdefault(T, []).append(lf)
                if p == "PYRRHOTITE" and np.isfinite(x):
                    pocomp.setdefault(float(T), []).append(float(x))
    bands = {p: {T: (min(L), max(L)) for T, L in d.items()} for p, d in out.items()}
    return bands, pocomp


fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.4), sharey=True)
titles = {
    "po_out": "Lower bound — greigite +1σ stable, pyrrhotite −1σ, pyrite +1σ\n→ pyrrhotite ELIMINATED",
    "po_in": "Upper bound — greigite −1σ, pyrrhotite +1σ, pyrite −1σ\n→ LARGEST pyrrhotite field",
}
order = ["BCC_A2", "PYRRHOTITE", "GREIGITE", "PYRITE", "ORTHORHOMBIC_S"]
summary = {}
for ax, (name, tdb) in zip(
    axes,
    [
        ("po_out", str(TDB / "fes_greigite_v1_dHf_lo.tdb")),
        ("po_in", str(TDB / "fes_greigite_v1_dHf_hi.tdb")),
    ],
):
    bands, poc = fields(tdb)
    summary[name] = {"pyrrhotite_present": "PYRRHOTITE" in bands}
    if poc:
        allc = [x for L in poc.values() for x in L]
        summary[name]["po_xS_range"] = [min(allc), max(allc)]
    for p in order:
        if p not in bands:
            continue
        Ts = sorted(bands[p])
        lo = [bands[p][t][0] for t in Ts]
        hi = [bands[p][t][1] for t in Ts]
        ax.fill_betweenx(Ts, lo, hi, color=COL[p], alpha=0.55, lw=0)
        ax.plot(lo, Ts, color=COL[p], lw=0.8)
        ax.plot(hi, Ts, color=COL[p], lw=0.8)
    # DSC/XRD pyritization band 590-723 K (cutoff justification)
    ax.axhspan(590, 600, color="grey", alpha=0.25, hatch="////", lw=0)
    ax.axhline(590, color="k", ls="--", lw=1)
    ax.text(
        ax.get_xlim()[0] + 0.3,
        592,
        "pyritization onset 317°C (590 K) — DSC/XRD",
        fontsize=7,
        va="bottom",
    )
    ax.set_title(titles[name], fontsize=9.5)
    ax.set_xlabel("log$_{10}$ f(S$_2$) / bar")
    ax.set_ylim(300, 600)
    ax.grid(alpha=0.25)
axes[0].set_ylabel("Temperature / K")
handles = [Patch(facecolor=COL[p], alpha=0.55, label=LAB[p]) for p in order]
axes[1].legend(handles=handles, loc="lower right", fontsize=8, framealpha=0.9)
fig.suptitle(
    "Fe–S equilibrium — two error-bounding cases (pycalphad engine, Dilner+greigite TDB)\nAll three sulfide ΔH$_f$ moved to their pyrrhotite-suppressing / pyrrhotite-maximizing 1σ limits",
    fontsize=10.5,
)
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(str(FIG / "fes_boundary_cases.png"), dpi=150)
print(json.dumps(summary, indent=2))
print("SAVED fes_boundary_cases.png")
