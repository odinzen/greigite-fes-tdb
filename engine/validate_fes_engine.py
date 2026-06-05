#!/usr/bin/env python
"""Engine validation: build the Fe-S log f(S2)-vs-T predominance diagram from
the REAL pycalphad solver on the greigite-grafted Dilner TDB, and compare the
greigite/pyrrhotite/pyrite topology against the numpy prototype.

Run (from repo root):
    python engine/validate_fes_engine.py

Reads the greigite TDB variants from artifacts/tdb/ (build them first with
build_greigite_tdb.py). Outputs:
    artifacts/figures/fes_engine_diagram.png  -- log f(S2) vs T fields
    artifacts/fes_engine_boundaries.json      -- machine-readable boundary table
Everything printed to stdout is also what to paste back for review.

Method
------
At each T we sweep bulk X(S) over condensed-phase-only equilibria and read the
S chemical potential mu_S = v.MU('S') of the stable assemblage, then convert to
the sulfur fugacity with the standard Kellogg-diagram convention:

    mu_S2(g) = G_S2_ref(T) + R*T*ln(p_S2/bar)   and   mu_S2(g) = 2*mu_S
    => log10 f(S2) = (2*mu_S - G_S2_ref(T)) / (R*T*ln10)

G_S2_ref(T) is the SER-referenced standard Gibbs energy of S2(g) from the TDB
GAS phase (FUNCTION F15281T; the +R*T*ln(1e-5*P) term is 0 at P=1 bar).
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pycalphad import Database, equilibrium, variables as v

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TDB_DIR = ROOT / "artifacts" / "tdb"
FIG_DIR = ROOT / "artifacts" / "figures"
R = 8.31451  # SGTE91 R, matches the TDB 'R' FUNCTION
LN10 = math.log(10.0)
COMPS = ["FE", "S", "VA"]
# condensed phases only -> mu_S reflects the condensed assemblage
CONDENSED = ["BCC_A2", "FCC_A1", "PYRRHOTITE", "PYRITE", "GREIGITE", "ORTHORHOMBIC_S"]
S2_REF_FUNCTION = "F15281T"  # G_S2_ref(T), SER-referenced, J/mol-S2
T_GRID = [300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0]
XS_GRID = np.linspace(0.34, 0.665, 70)  # mole-fraction S across Po..pyrite

VARIANTS = {
    "central": "fes_greigite_v1.tdb",
    "lo_dHf": "fes_greigite_v1_dHf_lo.tdb",
    "hi_dHf": "fes_greigite_v1_dHf_hi.tdb",
}


def _g_s2_ref(db, T: float) -> float:
    """Evaluate the S2(g) standard Gibbs energy FUNCTION at T (J/mol-S2)."""
    from pycalphad.variables import T as T_VAR

    return float(db.symbols[S2_REF_FUNCTION].subs({T_VAR: float(T)}))


def _phases_and_muS(db, T: float, x_s: float):
    """Return (list_of_phases, mu_S) at bulk X(S)=x_s, T (condensed only).

    The assemblage may be one or two phases (a tie-line). mu_S = v.MU('S') is
    the equilibrium sulfur chemical potential of that assemblage.
    """
    eq = equilibrium(
        db, COMPS, CONDENSED, {v.X("S"): float(x_s), v.T: T, v.P: 101325, v.N: 1}
    )
    phases = sorted({p for p in np.ravel(eq.Phase.values) if p})
    try:
        muS = float(np.ravel(eq.MU.sel(component="S").values)[0])
    except Exception:
        muS = float(np.ravel(eq.MU.values)[-1])
    return phases, muS


def build_fields(db, label: str):
    """Sweep (T, X_S); return dict: phase -> {T: (logfs2_min, logfs2_max)}.

    A line compound only ever appears in 2-phase tie-lines; the SET of mu_S
    (=> log f(S2)) values over which a phase is present in the assemblage spans
    exactly its predominance band (from its lower tie-line to its upper one).
    So we record log f(S2) for EVERY phase present in EACH assemblage and take
    the [min, max] span per phase per temperature.
    """
    rows = []
    fields: dict[str, dict[float, list]] = {}
    for T in T_GRID:
        gref = _g_s2_ref(db, T)
        for x in XS_GRID:
            phases, muS = _phases_and_muS(db, T, x)
            logfs2 = (2.0 * muS - gref) / (R * T * LN10)
            rows.append((T, float(x), "+".join(phases), logfs2))
            for ph in phases:
                fields.setdefault(ph, {}).setdefault(T, []).append(logfs2)
    out = {
        p: {T: (min(val), max(val)) for T, val in d.items()} for p, d in fields.items()
    }
    return out, rows


def consistency_sigma(db, T=298.15):
    """Greigite decomposition margin on this TDB basis (kJ, sigma)."""

    # Hf of pyrite and pyrrhotite(FeS1.092) and greigite at their compositions
    def hm(phase, x_s):
        eq = equilibrium(
            db,
            COMPS,
            [phase],
            {v.X("S"): x_s, v.T: T, v.P: 101325, v.N: 1},
            output="HM",
        )
        return float(np.ravel(eq.HM.values)[0])  # J/mol-atom

    try:
        H_py = hm("PYRITE", 2 / 3) * 3.0 / 1000.0  # per mol FeS2 (3 atoms)
        H_po = hm("PYRRHOTITE", 1.092 / 2.092) * 2.092 / 1000.0
        H_gr = hm("GREIGITE", 4 / 7) * 7.0 / 1000.0  # per mol Fe3S4
        # FeS1.33 -> 0.726 FeS1.092 + 0.274 FeS2 (per provenance)
        H_gr_FeS133 = H_gr / 3.0
        H_prod = 0.726 * H_po + 0.274 * H_py
        margin = H_prod - H_gr_FeS133  # >0 => greigite stable
        return {
            "H_py_kJ": H_py,
            "H_po_kJ": H_po,
            "H_gr_Fe3S4_kJ": H_gr,
            "decomp_margin_kJ_per_FeS133": margin,
            "sigma": margin / 7.3,
            "greigite_stable": margin > 0,
        }
    except Exception as e:  # HM path can vary by pycalphad build
        return {"error": f"{type(e).__name__}: {e}"}


def main():
    print("pycalphad engine validation — Fe-S greigite\n" + "=" * 56)
    results = {}
    central_fields = None
    for label, fname in VARIANTS.items():
        tdb = TDB_DIR / fname
        if not tdb.exists():
            print(f"  MISSING {fname} — skip")
            continue
        db = Database(str(tdb))
        print(f"\n[{label}]  {fname}")
        cons = consistency_sigma(db)
        print("  consistency @298.15K:", json.dumps(cons))
        # Full predominance sweep only for the central TDB (keeps runtime low);
        # lo/hi report the consistency margin only (shows the band stays stable).
        if label != "central":
            results[label] = {"consistency": cons}
            continue
        fields, rows = build_fields(db, label)
        # report greigite + neighbours boundaries
        for ph in ("BCC_A2", "PYRRHOTITE", "GREIGITE", "PYRITE"):
            if ph in fields:
                at = {
                    int(T): (round(r[0], 2), round(r[1], 2))
                    for T, r in sorted(fields[ph].items())
                }
                print(f"    {ph:11s} log f(S2) range by T: {at}")
            else:
                print(f"    {ph:11s} ABSENT (not stable anywhere on grid)")
        results[label] = {
            "consistency": cons,
            "fields": {
                p: {str(int(T)): list(r) for T, r in d.items()}
                for p, d in fields.items()
            },
        }
        central_fields = fields

    # ---- plot central diagram ----
    if central_fields:
        fig, ax = plt.subplots(figsize=(8.4, 6.2))
        colors = {
            "BCC_A2": "#d9d9d9",
            "PYRRHOTITE": "#e8d6a0",
            "GREIGITE": "#5fae77",
            "PYRITE": "#bcd4ec",
            "ORTHORHOMBIC_S": "#f0c0c0",
            "FCC_A1": "#cccccc",
        }
        T_arr = np.array(T_GRID)
        for ph, d in central_fields.items():
            Ts = sorted(d)
            lo = np.array([d[T][0] for T in Ts])
            hi = np.array([d[T][1] for T in Ts])
            ax.fill_between(
                Ts, lo, hi, color=colors.get(ph, "#eeeeee"), alpha=0.85, lw=0, label=ph
            )
        ax.set_xlabel("Temperature, K")
        ax.set_ylabel(r"log $f$(S$_2$), bar")
        ax.set_title(
            "Fe–S predominance from the pycalphad engine\n"
            "(Dilner 2015 base + measured greigite; central ΔH$_f$)"
        )
        ax.legend(loc="lower right", fontsize=8, framealpha=0.95)
        ax.set_xlim(min(T_GRID), max(T_GRID))
        fig.tight_layout()
        FIG_DIR.mkdir(parents=True, exist_ok=True)
        out_png = FIG_DIR / "fes_engine_diagram.png"
        fig.savefig(out_png, dpi=200)
        print(f"\nwrote {out_png}")

    out_json = ROOT / "artifacts" / "fes_engine_boundaries.json"
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(results, indent=2))
    print(f"wrote {out_json}")
    print(
        "\nDONE — paste the printed block above; outputs are in artifacts/ + engine/."
    )


if __name__ == "__main__":
    main()
