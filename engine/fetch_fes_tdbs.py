#!/usr/bin/env python
"""Fetch the campaign's literature TDBs from the public NIMS TDBDB into artifacts/tdb/.

This is the *fetch* step of the build: the base databases are downloaded into
``artifacts/tdb/`` on first run (``artifacts/`` is gitignored).

The set of TDBs to fetch — their itemids, source URLs, and destinations — is
declared in the campaign manifest (``provenance_manifest.json``); this script
just resolves each ``tdb_from_tdbdb`` artifact via ``manifest.resolve_tdb`` (a
live download + zip-extract through the vendored, self-contained ``tdbtools``
package). The builders also auto-fetch their own base on demand, so running this
explicitly is optional.

Everything generated is written under the repo-root ``artifacts/`` dir:
  artifacts/tdb/<itemid>.tdb     - the fetched databases
  artifacts/fetch_report.json    - machine-readable fetch outcome
Provenance is recorded in engine/provenance.md and the manifest.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

from manifest import Manifest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
ARTIFACTS = ROOT / "artifacts"

# Live download required to actually pull bytes from TDBDB / Elsevier.
os.environ.setdefault("TDBDB_LIVE", "1")


def main() -> None:
    m = Manifest.load()
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    report = {"campaign": m.campaign.get("id"), "num_tdbs": len(m.tdbs), "records": []}

    for tid, art in m.tdbs.items():
        entry = {
            "id": tid,
            "itemid": art.itemid,
            "role": art.role,
            "elements": art.elements,
            "url": art.fetch.get("url"),
            "fetch": None,
        }
        try:
            path = m.resolve_tdb(tid)
            entry["fetch"] = {
                "ok": True,
                "path": str(path),
                "bytes": path.stat().st_size,
            }
        except Exception as exc:  # noqa: BLE001 - report, don't crash the sweep
            entry["fetch"] = {"ok": False, "error": f"{type(exc).__name__}: {exc}"}
        report["records"].append(entry)

    (ARTIFACTS / "fetch_report.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
