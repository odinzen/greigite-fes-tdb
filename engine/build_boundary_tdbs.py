#!/usr/bin/env python
"""Build the TWO all-compound boundary-case TDBs for the Fig. 2 envelope.

Earlier work shifted ONLY greigite (+/-21.9 kJ/Fe3S4) in _dHf_lo/_hi.tdb. The
boundary cases require shifting the Dilner PYRITE and PYRRHOTITE G-functions by
their own measured 1-sigma too, so the error envelope is propagated for ALL
compounds (greigite +/-7.3, pyrite +/-6.0, pyrrhotite +/-3.5 kJ per f.u.).

Two labelled cases (feedback 4):
  LOWER  (pyrrhotite eliminated / smallest pyrrhotite field):
      greigite at -1 sigma (MOST stable)  -> use _dHf_lo greigite function
      pyrite    at +1 sigma (LEAST stable) -> +6000 J on the -177763 constant
      pyrrhotite at +1 sigma (LEAST stable)-> +3500 J on GFES (-107518)
  UPPER  (largest pyrrhotite field):
      greigite at +1 sigma (LEAST stable) -> use _dHf_hi greigite function
      pyrite    at -1 sigma (MOST stable) -> -6000 J on the -177763 constant
      pyrrhotite at -1 sigma (MOST stable)-> -3500 J on GFES (-107518)

Method: pure, logged string surgery on the already-validated greigite +/-sigma
TDBs. The greigite functions themselves are NOT rebuilt (they reproduce
byte-for-byte from build_greigite_tdb.py); we only add the pyrite/pyrrhotite
1-sigma shifts to the Dilner REF2 (Lee et al.) functions.

Sigma sources (per f.u., from Subramani 2020 Table 1 / Xu & Navrotsky 2010):
  greigite -144.1 +/- 7.3 kJ/mol-FeS1.33   (x3 -> +/-21.9 kJ/Fe3S4)
  pyrite   -175.5 +/- 6.0 kJ/mol-FeS2
  pyrrhotite -106.2 +/- 3.5 kJ/mol-FeS1.092 (shift applied to the FeS unit GFES)

NOTE (traceability): the +/-3.5 is per FeS1.092; we apply it to GFES, the FeS
endmember energy of the (FE,VA):(S) sublattice model. This shifts the whole
pyrrhotite manifold by 3.5 kJ/mol-FeS, a defensible 1-sigma envelope on the
non-stoichiometric phase (documented, not a silent approximation).

I/O: reads the greigite +/-sigma TDBs from artifacts/tdb/ (build them first with
build_greigite_tdb.py), writes the two boundary TDBs back to artifacts/tdb/, and
writes the machine-readable boundary_cases_report.json to artifacts/ (a generated
report, not committed).
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
from pycalphad import Database, equilibrium, variables as v

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TDB_DIR = ROOT / "artifacts" / "tdb"

PYRITE_TOKEN = "-177763"  # Dilner/Lee pyrite FE:S constant (REF2)
GFES_TOKEN = "-107518"  # Dilner/Lee GFES (pyrrhotite FeS unit, REF2)

CASES = {
    "fes_greigite_boundary_lower.tdb": {
        "src": "fes_greigite_v1_dHf_lo.tdb",  # greigite -1 sigma (most stable)
        "d_pyrite": +6000,  # pyrite least stable
        "d_pyrrho": +3500,  # pyrrhotite least stable
        "tag": "LOWER: greigite -1s, pyrite +1s, pyrrhotite +1s (pyrrhotite minimised/eliminated)",
    },
    "fes_greigite_boundary_upper.tdb": {
        "src": "fes_greigite_v1_dHf_hi.tdb",  # greigite +1 sigma (least stable)
        "d_pyrite": -6000,  # pyrite most stable
        "d_pyrrho": -3500,  # pyrrhotite most stable
        "tag": "UPPER: greigite +1s, pyrite -1s, pyrrhotite -1s (pyrrhotite field maximised)",
    },
}

# ----- engine topology probe (mirrors validate_fes_engine.py conventions) ----
R = 8.31451
LN10 = math.log(10.0)
COMPS = ["FE", "S", "VA"]
CONDENSED = ["BCC_A2", "FCC_A1", "PYRRHOTITE", "PYRITE", "GREIGITE", "ORTHORHOMBIC_S"]
S2_REF_FUNCTION = "F15281T"
T_GRID = [300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0]
XS_GRID = np.linspace(0.34, 0.665, 70)


def _shift_token(text: str, token: str, delta: int) -> str:
    n = text.count(token)
    if n != 1:
        raise RuntimeError(f"token {token!r} appears {n}x (expected 1) — refusing")
    return text.replace(token, f"{token}{delta:+d}", 1)


def _g_s2_ref(db, T):
    from pycalphad.variables import T as T_VAR

    return float(db.symbols[S2_REF_FUNCTION].subs({T_VAR: float(T)}))


def _phases_and_muS(db, T, x_s):
    eq = equilibrium(
        db, COMPS, CONDENSED, {v.X("S"): float(x_s), v.T: T, v.P: 101325, v.N: 1}
    )
    phases = sorted({p for p in np.ravel(eq.Phase.values) if p})
    try:
        muS = float(np.ravel(eq.MU.sel(component="S").values)[0])
    except Exception:
        muS = float(np.ravel(eq.MU.values)[-1])
    return phases, muS


def probe(db):
    fields = {}
    for T in T_GRID:
        gref = _g_s2_ref(db, T)
        for x in XS_GRID:
            phases, muS = _phases_and_muS(db, T, x)
            logfs2 = (2.0 * muS - gref) / (R * T * LN10)
            for ph in phases:
                fields.setdefault(ph, {}).setdefault(T, []).append(logfs2)
    return {
        p: {T: (round(min(val), 2), round(max(val), 2)) for T, val in d.items()}
        for p, d in fields.items()
    }


def main():
    import sys, json as _json

    only = sys.argv[1] if len(sys.argv) > 1 else None  # optional case-name filter
    TDB_DIR.mkdir(parents=True, exist_ok=True)
    rep_path = ROOT / "artifacts" / "boundary_cases_report.json"
    report = _json.loads(rep_path.read_text()) if rep_path.exists() else {}
    print("Boundary-case TDB build (all-compound 1-sigma envelope)\n" + "=" * 58)
    for out_name, c in CASES.items():
        if only and only not in out_name:
            continue
        src = (TDB_DIR / c["src"]).read_text()
        txt = _shift_token(src, PYRITE_TOKEN, c["d_pyrite"])
        txt = _shift_token(
            txt, GFES_TOKEN, c["d_pyrro"] if "d_pyrro" in c else c["d_pyrrho"]
        )
        # banner so the file self-documents the edit
        banner = (
            f"$ BOUNDARY CASE — {c['tag']}\n"
            f"$ edits vs {c['src']}: PYRITE {c['d_pyrite']:+d} J on {PYRITE_TOKEN}; "
            f"GFES {c['d_pyrrho']:+d} J on {GFES_TOKEN}. Built by build_boundary_tdbs.py.\n"
        )
        txt = banner + txt
        out = TDB_DIR / out_name
        out.write_text(txt)
        sha = hashlib.sha256(txt.encode()).hexdigest()
        print(f"\n[{out_name}]  {c['tag']}")
        print(
            f"  pyrite {c['d_pyrite']:+d} J, pyrrhotite(GFES) {c['d_pyrrho']:+d} J  sha={sha[:12]}"
        )
        db = Database(str(out))  # parse check
        fields = probe(db)
        po = fields.get("PYRRHOTITE")
        print(f"  PYRRHOTITE: {'PRESENT' if po else 'ABSENT (eliminated)'}")
        for ph in ("BCC_A2", "PYRRHOTITE", "GREIGITE", "PYRITE"):
            if ph in fields:
                at = {int(T): list(r) for T, r in sorted(fields[ph].items())}
                print(f"    {ph:11s} log f(S2) by T: {at}")
        report[out_name] = {
            "tag": c["tag"],
            "sha256": sha,
            "pyrrhotite_present": bool(po),
            "fields": {
                p: {str(int(T)): list(r) for T, r in d.items()}
                for p, d in fields.items()
            },
        }
        rep_path.write_text(json.dumps(report, indent=2))  # persist per case
    print(
        "\nwrote artifacts/boundary_cases_report.json + the boundary TDB(s) in artifacts/tdb/"
    )


if __name__ == "__main__":
    main()
