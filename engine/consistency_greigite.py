#!/usr/bin/env python
"""PHASE B.5 — greigite stability consistency check on the DILNER basis.

Greigite is stable iff its decomposition to the neighbouring phases is
endothermic (enthalpy-dominated, per Subramani 2020 / Shumway 2022):
    FeS1.33  ->  0.726 FeS1.092 (pyrrhotite) + 0.274 FeS2 (pyrite)        (per mole Fe)
The boundary greigite enthalpy (ΔH_decomp = 0) is
    ΔHf_bound = 0.726·ΔHf(pyrrhotite) + 0.274·ΔHf(pyrite).
Greigite must lie BELOW that bound; the margin is ΔH_decomp = ΔHf_bound − ΔHf(greigite),
and the σ-distance = ΔH_decomp / u(ΔHf,greigite)  with u = 7.3 kJ/mol (Subramani, absolute).

KEY CHANGE vs any earlier check: ΔHf(pyrrhotite) and ΔHf(pyrite) are taken from the
DILNER 2015 TDB (Lee/Selleby basis) — computed with pycalphad as the formation
enthalpy of those phases at 298.15 K — NOT from JANAF/Grønvold. The greigite ΔHf is
the measured Subramani value (provenance_manifest.json). The old 0.27σ figure (JANAF
basis) is deliberately not reused.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from pycalphad import Database, equilibrium, variables as v

from manifest import Manifest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TDB_DIR = ROOT / "artifacts" / "tdb"
MANIFEST = Manifest.load()
BASE = MANIFEST.tdb("tdb:dilner2015_femns").dest_path

# decomposition stoichiometry (per mole Fe), Shumway Eq. 6
A_PO, A_PY = 0.726, 0.274
X_PO = 1.092 / 2.092  # mole fraction S in FeS1.092
ATOMS_PO = 2.092
X_PY = 2.0 / 3.0  # mole fraction S in FeS2
ATOMS_PY = 3.0


def dHf_phase(db, phase, x_S, atoms_per_formula):
    """Formation enthalpy (J per formula) at 298.15 K from the Dilner TDB.

    At 298.15 K the SER element enthalpies are zero, so the phase HM (per mole
    atom) times atoms-per-formula is the standard enthalpy of formation.
    """
    eq = equilibrium(
        db,
        ["FE", "S", "VA"],
        [phase],
        {v.X("S"): x_S, v.T: 298.15, v.P: 101325, v.N: 1},
        output="HM",
    )
    hm_atom = float(np.ravel(eq.HM.values)[0])
    return hm_atom * atoms_per_formula


def main() -> None:
    MANIFEST.resolve_tdb("tdb:dilner2015_femns")  # fetch base if missing
    db = Database(str(BASE))
    dHf = MANIFEST.measurement("meas:subramani2020:dHf_greigite")
    dHf_grei = dHf.value  # kJ/mol FeS1.33
    u = dHf.uncertainty_value  # 7.3 kJ/mol

    dHf_py = dHf_phase(db, "PYRITE", X_PY, ATOMS_PY) / 1000.0  # kJ/mol FeS2
    dHf_po = dHf_phase(db, "PYRRHOTITE", X_PO, ATOMS_PO) / 1000.0  # kJ/mol FeS1.092

    bound = A_PO * dHf_po + A_PY * dHf_py  # boundary greigite ΔHf
    dH_decomp = bound - dHf_grei  # >0 => greigite stable
    sigma = dH_decomp / u

    print("=== greigite stability consistency check — DILNER 2015 basis ===")
    print(f"  ΔHf(pyrite  FeS2,    Dilner) = {dHf_py:8.2f} kJ/mol")
    print(
        "     (pycalphad SER formation enthalpy; NOT the bare -177.763 a-coefficient,"
    )
    print("      which is not ΔfH because GHSERFE/GHSERSS carry non-zero H at 298 K)")
    print(f"  ΔHf(pyrrhot FeS1.092,Dilner) = {dHf_po:8.2f} kJ/mol")
    print(f"  decomposition: FeS1.33 -> {A_PO} FeS1.092 + {A_PY} FeS2")
    print(f"  ΔHf boundary (ΔH_decomp=0)   = {bound:8.2f} kJ/mol")
    print(f"  ΔHf(greigite, Subramani)     = {dHf_grei:8.2f} +/- {u} kJ/mol")
    print(f"  ΔH_decomp (margin, >0=stable)= {dH_decomp:8.2f} kJ/mol")
    print(
        f"  σ-distance (margin / 7.3)    = {sigma:8.2f} σ  -> greigite "
        f"{'STABLE' if dH_decomp > 0 else 'UNSTABLE'} on the Dilner basis"
    )
    print()
    print("  greigite ΔHf limit variants (±3×7.3 not used here; per-FeS1.33 ±7.3):")
    for label, d in (
        ("lower (−7.3)", dHf_grei - u),
        ("nominal", dHf_grei),
        ("upper (+7.3)", dHf_grei + u),
    ):
        m = bound - d
        print(
            f"    {label:14s} ΔHf={d:7.1f}  margin={m:6.2f} kJ/mol  "
            f"({m / u:+.2f}σ, {'stable' if m > 0 else 'unstable'})"
        )

    out = {
        "basis": "Dilner 2015 (Lee/Selleby) pyrrhotite + pyrite, via pycalphad HM @298.15K",
        "dHf_pyrite_kJ": dHf_py,
        "dHf_pyrrhotite_FeS1p092_kJ": dHf_po,
        "dHf_boundary_kJ": bound,
        "dHf_greigite_kJ": dHf_grei,
        "u_kJ": u,
        "dH_decomp_kJ": dH_decomp,
        "sigma_distance": sigma,
        "greigite_stable": bool(dH_decomp > 0),
        "note": "supersedes the old JANAF-basis 0.27σ figure (not reused)",
    }
    out_path = ROOT / "artifacts" / "consistency_report.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2))
    print(f"  -> {out_path}")


if __name__ == "__main__":
    main()
