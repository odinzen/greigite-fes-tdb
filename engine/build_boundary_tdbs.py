#!/usr/bin/env python
"""Build the two uncertainty-bound TDBs for the Fig. 2 envelope.

Only greigite carries a bound: its enthalpy of formation is shifted by the reported
uncertainty, 7.3 kJ per FeS1.33 (21.9 kJ per Fe3S4, 2 SD of the mean, Subramani et al.
2020 Table 1). Pyrite and pyrrhotite stay at the assessed Dilner (Lee) functions. The
6.0 and 3.5 kJ quoted beside them in Subramani's Table 1 are Xu & Navrotsky (2010)
drop-solution errors on Xu's own values (-175.5, -106.2 kJ), not on the database values
used here, so applying them would mix database values with experimental errors.

  LOWER  greigite at -u (most stable)   <- _dHf_lo
  UPPER  greigite at +u (least stable)  <- _dHf_hi

The files are the greigite +/-u TDBs from build_greigite_tdb.py with a banner, kept
under the boundary names so the probe and figure scripts read one pair of files.

I/O: reads artifacts/tdb/fes_greigite_v1_dHf_{lo,hi}.tdb, writes
artifacts/tdb/fes_greigite_boundary_{lower,upper}.tdb and the generated
artifacts/boundary_cases_report.json (not committed).
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

CASES = {
    "fes_greigite_boundary_lower.tdb": {
        "src": "fes_greigite_v1_dHf_lo.tdb",
        "tag": "LOWER: greigite -u (most stable), pyrite and pyrrhotite assessed",
    },
    "fes_greigite_boundary_upper.tdb": {
        "src": "fes_greigite_v1_dHf_hi.tdb",
        "tag": "UPPER: greigite +u (least stable), pyrite and pyrrhotite assessed",
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
    report = _json.loads(rep_path.read_text(encoding="utf-8")) if rep_path.exists() else {}
    print("Boundary-case TDB build (greigite-only uncertainty)\n" + "=" * 52)
    for out_name, c in CASES.items():
        if only and only not in out_name:
            continue
        src = (TDB_DIR / c["src"]).read_text(encoding="utf-8")
        banner = (
            f"$ BOUNDARY CASE - {c['tag']}\n"
            f"$ copy of {c['src']}. Built by build_boundary_tdbs.py.\n"
        )
        txt = banner + src
        out = TDB_DIR / out_name
        out.write_text(txt, encoding="utf-8")
        sha = hashlib.sha256(txt.encode()).hexdigest()
        print(f"\n[{out_name}]  {c['tag']}")
        print(f"  sha={sha[:12]}")
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
        rep_path.write_text(json.dumps(report, indent=2), encoding="utf-8")  # persist per case
    print(
        "\nwrote artifacts/boundary_cases_report.json + the boundary TDB(s) in artifacts/tdb/"
    )


if __name__ == "__main__":
    main()
