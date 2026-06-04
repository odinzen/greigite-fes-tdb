#!/usr/bin/env python
"""Sanity-check the greigite-grafted TDB.

1. Global equilibrium at the Fe3S4 composition (X(S)=4/7) over 298-600 K,
   solid phases only -> what is actually stable (greigite is expected to be
   METASTABLE w.r.t. pyrite + pyrrhotite, so it should NOT appear).
2. Metastability margin: G(greigite) minus the stable assemblage's Gibbs
   energy at the same composition (how far above the hull greigite sits).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from pycalphad import Database, equilibrium, variables as v

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TDB = ROOT / "artifacts" / "tdb" / "fes_greigite_v1.tdb"

COMPS = ["FE", "S", "VA"]
SOLIDS = [
    "BCC_A2",
    "FCC_A1",
    "PYRRHOTITE",
    "PYRITE",
    "ALABANDITE",
    "ORTHORHOMBIC_S",
    "GREIGITE",
]
X_FE3S4 = 4.0 / 7.0  # mole fraction S in Fe3S4


def main() -> None:
    db = Database(str(TDB))

    print("=== Global equilibrium at X(S)=4/7 (Fe3S4 composition), solids only ===")
    for T in (298.15, 400.0, 500.0, 600.0):
        eq = equilibrium(
            db, COMPS, SOLIDS, {v.X("S"): X_FE3S4, v.T: T, v.P: 101325, v.N: 1}
        )
        phases = sorted({p for p in np.ravel(eq.Phase.values) if p})
        gm = float(np.ravel(eq.GM.values)[0])
        print(f"  T={T:6.2f} K  stable: {phases}  GM={gm:.1f} J/mol-atom")

    print("\n=== Greigite metastability margin (per mole Fe3S4 = 7 atoms) ===")
    # G of pure greigite line compound vs stable hull at same composition.
    for T in (298.15, 600.0):
        eq_all = equilibrium(
            db, COMPS, SOLIDS, {v.X("S"): X_FE3S4, v.T: T, v.P: 101325, v.N: 1}
        )
        g_hull = float(np.ravel(eq_all.GM.values)[0])  # J/mol-atom of stable assemblage
        eq_g = equilibrium(
            db, COMPS, ["GREIGITE"], {v.X("S"): X_FE3S4, v.T: T, v.P: 101325, v.N: 1}
        )
        g_grei = float(np.ravel(eq_g.GM.values)[0])  # J/mol-atom of greigite
        dG_atom = g_grei - g_hull
        print(
            f"  T={T:6.2f} K  G(greigite)={g_grei:.1f}  G(hull)={g_hull:.1f}  "
            f"margin=+{dG_atom:.1f} J/mol-atom (+{dG_atom * 7 / 1000:.2f} kJ/mol-Fe3S4)"
        )


if __name__ == "__main__":
    main()
