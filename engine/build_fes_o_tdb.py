#!/usr/bin/env python
"""Build the Fe-S-O greigite TDB (artifacts/tdb/fes_o_greigite_v1.tdb).

This is the missing Fe-S-O builder. The four Fe-S-O figures
(make_fig4_feso_control_bw, make_fig6_feso_bw, make_fig_predominance_split,
make_fig3_engine) consume ``fes_o_greigite_v1.tdb``. This script fetches the
Dilner-2017 base at build time and applies two mechanical transforms (no new
CALPHAD modelling):

  STEP 0  fetch the Dilner & Selleby 2017 Ca-Fe-O-S database (itemid
          ``calphadj_1-s2.0-S0364591616301584-mmc1_1``) into artifacts/tdb/.
          The Elsevier zip member is ``FeCaOS.TDB.txt`` (matched by the
          ``.tdb.txt`` rule in tdbtools/tdbdb.py).

  STEP 1  dedupe -> artifacts/tdb/fes_o_dilner2017_clean.tdb
          Normalise CRLF->LF, then drop the SECOND (duplicate) occurrence of the
          three FUNCTION blocks pycalphad rejects (GCAOSOL, GWUSTITE, GFES). No
          thermodynamic values change (882 -> 876 content lines).

  STEP 2  graft -> artifacts/tdb/fes_o_greigite_v1.tdb
          Append the measured GREIGITE block (GREIGITE_GF function + GREIGITE
          phase), reused verbatim from artifacts/tdb/fes_greigite_v1.tdb (built
          by build_greigite_tdb.py), plus the Fe-only PYRITE block from the
          Dilner-2015 / Lee et al. assessment.

The repo therefore redistributes only our greigite + pyrite additions and this
transform code, never Dilner's database file.

Validation: the result parses cleanly under pycalphad. The boundary fugacities
(native-S sat log f(S2) = -13.96; Fe/Fe3S4 = -51.0) are checked downstream by
make_fig_predominance_split.py. Pass ``--verify <reference.tdb>`` to assert the
build reproduces a known-good file byte-for-byte (used during development; the
reference is not distributed).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from manifest import Manifest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TDB_DIR = ROOT / "artifacts" / "tdb"

MANIFEST = Manifest.load()
DILNER2017_TDB_ID = "tdb:dilner2017_cafeos"  # manifest artifact id
DILNER2017_ITEMID = MANIFEST.tdb(DILNER2017_TDB_ID).itemid
SOURCE_TDB = MANIFEST.tdb(DILNER2017_TDB_ID).dest_path
CLEAN_TDB = TDB_DIR / "fes_o_dilner2017_clean.tdb"
OUT_TDB = TDB_DIR / "fes_o_greigite_v1.tdb"
GREIGITE_TDB = ROOT / MANIFEST.derived["tdb:fes_greigite_v1"].output  # built upstream

# FUNCTIONs the source defines twice; pycalphad rejects duplicates. Keep the
# first occurrence of each, drop the rest.
DUP_FUNCS = ("GCAOSOL", "GWUSTITE", "GFES")

DEDUPE_BANNER = (
    f"$ DERIVED from {DILNER2017_ITEMID}.tdb (Dilner & Selleby 2017, Ca-Fe-O-S).\n"
    "$ Removed 3 BYTE-IDENTICAL duplicate FUNCTION blocks (GCAOSOL, GFES, "
    "GWUSTITE) that pycalphad rejects; no thermodynamic values changed. "
    "Dedupe only.\n"
)

# Fe-only PYRITE endmember from the Dilner-2015 / Lee et al. (REF2) assessment.
PYRITE_GRAFT = (
    " PHASE PYRITE  %  2 1 2 !\n"
    " CONSTITUENT PYRITE  :FE : S :  !\n"
    " PARAMETER G(PYRITE,FE:S;0)  2.98150E+02  +GHSERFE#+2*GHSERSS#-177763+48.567*T;"
    "  6.00000E+03  N REF0 !\n"
)


# --------------------------------------------------------------------------- #
# STEP 0 — fetch Dilner-2017
# --------------------------------------------------------------------------- #
def fetch_dilner2017() -> Path:
    """Fetch the Dilner-2017 Ca-Fe-O-S TDB into artifacts/tdb/ (per the manifest)."""
    return MANIFEST.resolve_tdb(DILNER2017_TDB_ID)


# --------------------------------------------------------------------------- #
# STEP 1 — dedupe
# --------------------------------------------------------------------------- #
def dedupe_functions(text_lf: str) -> str:
    """Drop the 2nd+ occurrence of each DUP_FUNCS FUNCTION statement.

    Operates on whole lines so all other bytes are preserved verbatim. A TDB
    FUNCTION statement runs from the ``FUNCTION`` keyword to its ``!`` terminator
    (possibly multi-line).
    """
    lines = text_lf.split("\n")
    seen: set[str] = set()
    drop: set[int] = set()
    i, n = 0, len(lines)
    while i < n:
        toks = lines[i].split()
        if len(toks) >= 2 and toks[0] == "FUNCTION":
            name = toks[1]
            j = i
            while j < n and "!" not in lines[j]:
                j += 1
            if name in DUP_FUNCS:
                if name in seen:
                    drop.update(range(i, j + 1))
                else:
                    seen.add(name)
            i = j + 1
        else:
            i += 1
    return "\n".join(ln for idx, ln in enumerate(lines) if idx not in drop)


def build_clean(source_path: Path) -> str:
    raw = source_path.read_bytes().decode("latin-1")
    text_lf = raw.replace("\r\n", "\n")
    return DEDUPE_BANNER + dedupe_functions(text_lf)


# --------------------------------------------------------------------------- #
# STEP 2 — graft GREIGITE (from fes_greigite_v1.tdb) + PYRITE
# --------------------------------------------------------------------------- #
def ensure_greigite_tdb() -> Path:
    """Build the Fe-S greigite TDB if it isn't in artifacts/tdb/ yet."""
    if GREIGITE_TDB.exists():
        return GREIGITE_TDB
    print(f"{GREIGITE_TDB.name} not found; running build_greigite_tdb.py ...")
    import build_greigite_tdb

    build_greigite_tdb.main()
    if not GREIGITE_TDB.exists():
        raise RuntimeError(f"build_greigite_tdb.py did not produce {GREIGITE_TDB}")
    return GREIGITE_TDB


def extract_greigite_lines(greigite_text: str) -> tuple[str, str, str, str]:
    """Pull the 4 grafted GREIGITE lines from a built fes_greigite_v1.tdb."""
    func = phase = const = param = None
    for ln in greigite_text.split("\n"):
        s = ln.strip()
        if s.startswith("FUNCTION GREIGITE_GF "):
            func = ln
        elif s.startswith("PHASE GREIGITE "):
            phase = ln
        elif s.startswith("CONSTITUENT GREIGITE "):
            const = ln
        elif s.startswith("PARAMETER G(GREIGITE,FE:S;0)"):
            param = ln
    missing = [
        n
        for n, val in (
            ("FUNCTION GREIGITE_GF", func),
            ("PHASE GREIGITE", phase),
            ("CONSTITUENT GREIGITE", const),
            ("PARAMETER G(GREIGITE)", param),
        )
        if val is None
    ]
    if missing:
        raise RuntimeError(f"GREIGITE lines not found in greigite TDB: {missing}")
    if not func.rstrip().endswith("!"):
        raise RuntimeError(
            "GREIGITE_GF FUNCTION is not single-line; cannot graft verbatim"
        )
    return func, phase, const, param


def build_greigite_o(clean_text: str, greigite_text: str) -> str:
    func, phase, const, param = extract_greigite_lines(greigite_text)
    graft = (
        "\n".join(
            [
                "",
                "$ ===== GREIGITE graft (Subramani 2020 dHf + Shumway 2022 S/Cp; "
                "same as fes_greigite_v1) =====",
                func,
                phase,
                const,
                param,
                "",
                "$ ===== PYRITE graft (Dilner 2015 / Lee et al. REF2; Fe-only, Mn dropped) =====",
            ]
        )
        + "\n"
        + PYRITE_GRAFT
    )
    return clean_text + graft


# --------------------------------------------------------------------------- #
def main() -> int:
    ap = argparse.ArgumentParser(description="Build the Fe-S-O greigite TDB.")
    ap.add_argument(
        "--verify-clean",
        type=Path,
        default=None,
        help="assert the dedupe output matches this reference byte-for-byte",
    )
    ap.add_argument(
        "--verify",
        type=Path,
        default=None,
        help="assert the final TDB matches this reference byte-for-byte",
    )
    args = ap.parse_args()

    TDB_DIR.mkdir(parents=True, exist_ok=True)

    # STEP 0
    src = fetch_dilner2017()
    print(f"[0] source: {src}  ({src.stat().st_size} B)")

    # STEP 1
    clean_text = build_clean(src)
    CLEAN_TDB.write_text(clean_text)
    print(
        f"[1] dedupe: {CLEAN_TDB}  ({len(clean_text.encode('latin-1'))} B, "
        f"{clean_text.count(chr(10))} lines)"
    )
    if args.verify_clean:
        ref = args.verify_clean.read_bytes()
        ok = clean_text.encode("latin-1") == ref
        print(
            f"    verify-clean vs {args.verify_clean.name}: "
            f"{'BYTE-IDENTICAL' if ok else 'MISMATCH'}"
        )
        if not ok:
            return 2

    # STEP 2
    greigite_text = ensure_greigite_tdb().read_text()
    out_text = build_greigite_o(clean_text, greigite_text)
    OUT_TDB.write_text(out_text)
    print(
        f"[2] graft:  {OUT_TDB}  ({len(out_text.encode('latin-1'))} B, "
        f"{out_text.count(chr(10))} lines)"
    )
    if args.verify:
        ref = args.verify.read_bytes()
        ok = out_text.encode("latin-1") == ref
        print(
            f"    verify vs {args.verify.name}: "
            f"{'BYTE-IDENTICAL' if ok else 'MISMATCH'}"
        )
        if not ok:
            return 2

    # parse check
    try:
        from pycalphad import Database

        db = Database(str(OUT_TDB))
        print(
            f"[OK] pycalphad parse: {len(db.phases)} phases, "
            f"{len(db.elements)} elements; GREIGITE present: "
            f"{'GREIGITE' in db.phases}, PYRITE present: {'PYRITE' in db.phases}"
        )
    except ImportError:
        print("[--] pycalphad not installed; skipped parse check")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
