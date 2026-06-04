#!/usr/bin/env python
"""Pull full (multi-line) FUNCTION + PARAMETER blocks of interest from the
fetched Fe-S TDBs, and peek inside the 2017 Ca-Fe-O-S zip that had no .tdb.

TDB statements run from a keyword to the terminating '!'. We reassemble those
multi-line statements so truncated constants (e.g. the pyrite enthalpy) are
shown in full.
"""

from __future__ import annotations

import io
import re
import urllib.request
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
TDB_DIR = HERE.parent / "artifacts" / "tdb"

WANTED_FUNCS = ("GFES", "GFE1S1", "GFE1S2", "GFES2", "GMNS", "GHSERSS", "GPYRRH")
WANTED_PARAM_SUBSTR = ("PYRITE", "PYRRHOTITE")


def statements(raw: str):
    """Yield TDB statements (keyword ... '!') as single whitespace-joined strings."""
    # Strip $-comments first (TDB comment char), keep statement structure.
    lines = []
    for ln in raw.splitlines():
        # comments start at '$'
        idx = ln.find("$")
        lines.append(ln if idx < 0 else ln[:idx])
    text = "\n".join(lines)
    for stmt in text.split("!"):
        s = " ".join(stmt.split())
        if s:
            yield s


def dump(path: Path) -> None:
    print("\n" + "#" * 78)
    print("#", path.name)
    print("#" * 78)
    raw = path.read_text(errors="replace")
    stmts = list(statements(raw))

    print("\n--- FUNCTIONs of interest ---")
    for s in stmts:
        if s.upper().startswith("FUNCTION"):
            name = s.split()[1].upper() if len(s.split()) > 1 else ""
            if name in WANTED_FUNCS:
                print("  " + s[:400])

    print("\n--- PYRITE / PYRRHOTITE PARAMETERs (full) ---")
    for s in stmts:
        up = s.upper()
        if up.startswith("PARAMETER") and any(w in up for w in WANTED_PARAM_SUBSTR):
            print("  " + s[:400])

    print("\n--- TYPE_DEFINITION / magnetic / SPECIES FE1S2 etc ---")
    for s in stmts:
        up = s.upper()
        if up.startswith("SPECIES") and ("S2" in up or "S1" in up):
            print("  " + s[:200])


def peek_2017_zip() -> None:
    url = "http://ars.els-cdn.com/content/image/1-s2.0-S0364591616301584-mmc1.zip"
    print("\n" + "#" * 78)
    print("# 2017 Ca-Fe-O-S supplementary zip contents (Dilner & Selleby 2017)")
    print("#  url:", url)
    print("#" * 78)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        data = urllib.request.urlopen(req, timeout=30).read()  # nosec B310
        zf = zipfile.ZipFile(io.BytesIO(data))
        for n in zf.namelist():
            print(f"  member: {n}  ({zf.getinfo(n).file_size} bytes)")
        # If there is a .dat / FactSage file, sniff its head for the liquid model.
        for n in zf.namelist():
            if n.lower().endswith((".dat", ".tdb", ".txt")):
                head = zf.read(n)[:1200].decode("latin-1", "replace")
                print(f"\n  --- head of {n} ---")
                for ln in head.splitlines()[:25]:
                    print("    " + ln)
    except Exception as exc:  # noqa: BLE001
        print("  PEEK FAILED:", type(exc).__name__, exc)


def main() -> None:
    for f in sorted(TDB_DIR.glob("*.tdb")):
        dump(f)
    peek_2017_zip()


if __name__ == "__main__":
    main()
