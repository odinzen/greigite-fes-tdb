#!/usr/bin/env python
"""Parse-check + phase/provenance inspection of the fetched Fe-S TDBs.

For each *.tdb in artifacts/tdb/:
  1. pycalphad Database(path) parse check -> clean / error text.
  2. Elements, species, phases with sublattice constituents.
  3. Flag any non-stoichiometric Fe-S solution phase (pyrrhotite / FE1-xS).
  4. Dump REFERENCE / provenance lines from the raw text.
  5. Grep for FeS2 (pyrite) and pyrrhotite parameters + any G(FES2)/enthalpy.

Pure read-only; writes nothing. All output to stdout for the report.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TDB_DIR = HERE.parent / "artifacts" / "tdb"

from pycalphad import Database  # noqa: E402

SOLUTION_HINTS = ("PYRR", "FE1", "FEX", "MONOSULF", "TROIL", "FE_S")


def banner(s: str) -> None:
    print("\n" + "=" * 78)
    print(s)
    print("=" * 78)


def inspect(path: Path) -> None:
    banner(f"FILE: {path.name}  ({path.stat().st_size} bytes)")

    raw = path.read_text(errors="replace")

    # --- parse check -------------------------------------------------------
    try:
        db = Database(str(path))
        print("PARSE: OK (pycalphad Database loaded)")
    except Exception as exc:  # noqa: BLE001
        print(f"PARSE: FAILED -> {type(exc).__name__}: {exc}")
        db = None

    # --- phases ------------------------------------------------------------
    if db is not None:
        print(f"\nElements ({len(db.elements)}): {sorted(db.elements)}")
        print(f"Species  ({len(db.species)}): {sorted(s.name for s in db.species)}")
        print(f"\nPhases ({len(db.phases)}):")
        for name in sorted(db.phases):
            ph = db.phases[name]
            consts = [sorted(str(sp.name) for sp in subl) for subl in ph.constituents]
            model = " : ".join("(" + ", ".join(c) + ")" for c in consts)
            sites = list(ph.sublattices)
            print(f"  {name}")
            print(f"      sites: {sites}")
            print(f"      model: {model}")
        # solution-phase flag
        sol = [n for n in db.phases if any(h in n.upper() for h in SOLUTION_HINTS)]
        print(
            f"\nPossible non-stoichiometric Fe-S solution phases (by name): {sol or 'NONE'}"
        )

    # --- raw provenance ----------------------------------------------------
    banner(f"PROVENANCE / REFERENCE lines in {path.name}")
    prov_keys = (
        "waldner",
        "pelton",
        "dilner",
        "selleby",
        "mao",
        "quasichem",
        "mqm",
        "reference",
        "assess",
        "20",
    )
    shown = 0
    for ln in raw.splitlines():
        low = ln.lower()
        if any(
            k in low
            for k in (
                "waldner",
                "pelton",
                "dilner",
                "selleby",
                "mao",
                "quasichem",
                "mqmqa",
                "modified quasi",
            )
        ):
            print("  " + ln.strip()[:160])
            shown += 1
    if shown == 0:
        print("  (no author/assessment keyword hits in raw text)")

    # REFERENCE_FILE / REF blocks
    banner(f"REF blocks in {path.name}")
    for m in re.finditer(r"^\s*REF\w*.*$", raw, flags=re.MULTILINE | re.IGNORECASE):
        print("  " + m.group(0).strip()[:160])

    # --- FeS2 / pyrite + pyrrhotite params --------------------------------
    banner(f"FeS2 / pyrrhotite / FES parameter & function lines in {path.name}")
    for ln in raw.splitlines():
        up = ln.upper()
        if (
            "FES2" in up
            or "PYRIT" in up
            or "PYRRH" in up
            or "TROILIT" in up
            or re.search(r"\bFES\b", up)
        ):
            print("  " + ln.strip()[:200])


def main() -> None:
    files = sorted(TDB_DIR.glob("*.tdb"))
    if not files:
        print("No .tdb files in", TDB_DIR, file=sys.stderr)
        sys.exit(1)
    for f in files:
        inspect(f)


if __name__ == "__main__":
    main()
