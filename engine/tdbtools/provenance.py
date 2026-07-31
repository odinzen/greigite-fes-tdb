"""Structured provenance for every datum that flows through the build.

Every value (a literature TDB from TDBDB, a measured quantity from a paper)
carries a `Provenance` block recording:

  - WHERE the value came from (source + source-specific ID)
  - HOW to cite it (human-readable citation string + canonical URL)
  - HOW it was computed (functional, model version, table edition, etc.)
  - WHEN we accessed it

The blocks aggregate into the `provenance_manifest.json` index, so a downstream
user (or a paper reviewer) can audit every value in a fitted database back to
its original source.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from typing import Any


def _now_iso() -> str:
    """ISO 8601 UTC timestamp, second precision."""
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


@dataclass
class Provenance:
    """Where a single datum came from and how to cite it.

    Machine-readable (via `to_dict`) and human-readable (for citation in a
    paper). Attached to every literature TDB and measured value the build
    consumes.
    """

    source: str  # e.g. "literature_tdb"
    source_id: str  # source-specific unique identifier
    citation: str  # human-readable citation, BibTeX-friendly when possible
    url: str | None = None  # canonical URL where the data can be re-fetched
    method: dict[str, Any] = field(default_factory=dict)  # source-specific method info
    accessed_at: str = field(default_factory=_now_iso)

    def to_dict(self) -> dict:
        return asdict(self)


CITATION_TDBDB = (
    "TDBDB: a free, open-access database of thermodynamic data files. "
    "Brown University, https://avdwgroup.engin.brown.edu/"
)


def from_literature_tdb(
    tdbdb_itemid: str,
    elements: list[str],
    authoryear: str,
    paper_doi: str | None = None,
    tdb_url: str | None = None,
) -> Provenance:
    """Provenance for a literature CALPHAD database fetched via TDBDB."""
    citation = f"{authoryear}"
    if paper_doi:
        citation += f", doi:{paper_doi}"
    citation += f" | indexed via {CITATION_TDBDB}"
    return Provenance(
        source="literature_tdb",
        source_id=tdbdb_itemid,
        citation=citation,
        url=tdb_url,
        method={
            "elements": list(elements),
            "tdbdb_itemid": tdbdb_itemid,
            "paper_doi": paper_doi,
        },
    )
