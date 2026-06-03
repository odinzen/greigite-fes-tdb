#!/usr/bin/env python
"""Typed read API over the campaign provenance manifest.

`provenance_manifest.json` is the single, machine-readable, provenance-tracked
description of the greigite TDB campaign: the literature TDBs to fetch, the
measured values taken from papers, the collaborator experimental datasets, and
the recipe DAG that turns those into the derived TDBs and figures. The build
scripts read from here instead of from ad-hoc files, so every number is
auditable back to a typed source.

This module is the accessor the scripts (and agents) use. Stdlib only.

    from manifest import Manifest
    m = Manifest.load()
    dHf = m.measurement("meas:subramani2020:dHf_greigite")      # -> Measurement
    base = m.resolve_tdb("tdb:dilner2015_femns")                # -> Path (fetched if missing)
    dsc = m.experimental("exp:dsc_heating").resolved_path       # -> Path in the repo
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent  # engine/
ROOT = HERE.parent  # repo root
MANIFEST_PATH = HERE / "provenance_manifest.json"


# --------------------------------------------------------------------------- #
# Typed nodes
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Provenance:
    """Where a datum came from + how to cite it (mirrors the KG Provenance block)."""

    source: str
    source_id: str
    citation: str = ""
    url: str | None = None
    method: dict[str, Any] = field(default_factory=dict)
    accessed_at: str | None = None

    @classmethod
    def from_dict(cls, d: dict | None) -> "Provenance | None":
        if not d:
            return None
        return cls(
            source=d.get("source", ""),
            source_id=d.get("source_id", ""),
            citation=d.get("citation", ""),
            url=d.get("url"),
            method=d.get("method", {}) or {},
            accessed_at=d.get("accessed_at"),
        )


@dataclass(frozen=True)
class Paper:
    id: str
    doi: str | None
    title: str
    authors: list[str]
    year: int | None
    venue: str | None
    url: str | None
    kg: dict[str, Any]
    raw: dict[str, Any] = field(default_factory=dict, repr=False)


@dataclass(frozen=True)
class Measurement:
    """A measured physical quantity, associated with a Paper.

    The new abstraction the KG needs: a typed value (value/unit/basis/
    uncertainty/conditions) with the measurement method, the in-paper locator,
    and a provenance block. `model` carries a fitted form (e.g. Debye-Einstein
    Cp coefficients) when the "measurement" is a model rather than a scalar.
    """

    id: str
    paper: str
    quantity: str
    phase: str | None = None
    value: float | None = None
    unit: str | None = None
    basis: str | None = None
    uncertainty: dict[str, Any] | None = None
    conditions: dict[str, Any] = field(default_factory=dict)
    measurement_method: str | None = None
    locator: dict[str, Any] = field(default_factory=dict)
    derived_values: list[dict[str, Any]] = field(default_factory=list)
    model: dict[str, Any] | None = None
    role: str | None = None
    provenance: Provenance | None = None
    raw: dict[str, Any] = field(default_factory=dict, repr=False)

    @property
    def uncertainty_value(self) -> float | None:
        return None if not self.uncertainty else self.uncertainty.get("value")

    def coeff(self, name: str) -> float:
        """A fitted-model coefficient (e.g. Debye-Einstein 'ThetaD_K')."""
        if not self.model or "coefficients" not in self.model:
            raise KeyError(f"measurement {self.id} has no model coefficients")
        return float(self.model["coefficients"][name])

    @property
    def locator_str(self) -> str:
        t, p = self.locator.get("table"), self.locator.get("page")
        bits = [b for b in (t, f"p.{p}" if p is not None else None) if b]
        return ", ".join(bits)


@dataclass(frozen=True)
class TdbArtifact:
    id: str
    role: str | None
    elements: list[str]
    fetch: dict[str, Any]
    provenance: Provenance | None = None
    raw: dict[str, Any] = field(default_factory=dict, repr=False)

    @property
    def itemid(self) -> str:
        return self.fetch["itemid"]

    @property
    def dest_path(self) -> Path:
        return ROOT / self.fetch["dest"]


@dataclass(frozen=True)
class ExperimentalDataset:
    id: str
    technique: str | None
    run: str | None
    path: str
    provenance: Provenance | None = None
    raw: dict[str, Any] = field(default_factory=dict, repr=False)

    @property
    def resolved_path(self) -> Path:
        return ROOT / self.path


@dataclass(frozen=True)
class Derived:
    id: str
    type: str
    produced_by: str | None
    derived_from: list[str]
    output: Any = None
    method: str | None = None
    raw: dict[str, Any] = field(default_factory=dict, repr=False)


# --------------------------------------------------------------------------- #
# Manifest container
# --------------------------------------------------------------------------- #
class Manifest:
    def __init__(self, data: dict[str, Any]):
        self.data = data
        self.campaign: dict[str, Any] = data.get("campaign", {})
        self.papers: dict[str, Paper] = {}
        self.measurements: dict[str, Measurement] = {}
        self.tdbs: dict[str, TdbArtifact] = {}
        self.experimental: dict[str, ExperimentalDataset] = {}
        self.derived: dict[str, Derived] = {}

        for p in data.get("papers", []):
            self.papers[p["id"]] = Paper(
                id=p["id"],
                doi=p.get("doi"),
                title=p.get("title", ""),
                authors=p.get("authors", []),
                year=p.get("year"),
                venue=p.get("venue"),
                url=p.get("url"),
                kg=p.get("kg", {}),
                raw=p,
            )

        for a in data.get("artifacts", []):
            t = a.get("type")
            if t == "measurement":
                self.measurements[a["id"]] = Measurement(
                    id=a["id"],
                    paper=a["paper"],
                    quantity=a["quantity"],
                    phase=a.get("phase"),
                    value=a.get("value"),
                    unit=a.get("unit"),
                    basis=a.get("basis"),
                    uncertainty=a.get("uncertainty"),
                    conditions=a.get("conditions", {}),
                    measurement_method=a.get("measurement_method"),
                    locator=a.get("locator", {}),
                    derived_values=a.get("derived_values", []),
                    model=a.get("model"),
                    role=a.get("role"),
                    provenance=Provenance.from_dict(a.get("provenance")),
                    raw=a,
                )
            elif t == "tdb_from_tdbdb":
                self.tdbs[a["id"]] = TdbArtifact(
                    id=a["id"],
                    role=a.get("role"),
                    elements=a.get("elements", []),
                    fetch=a.get("fetch", {}),
                    provenance=Provenance.from_dict(a.get("provenance")),
                    raw=a,
                )
            elif t == "experimental_dataset":
                self.experimental[a["id"]] = ExperimentalDataset(
                    id=a["id"],
                    technique=a.get("technique"),
                    run=a.get("run"),
                    path=a["path"],
                    provenance=Provenance.from_dict(a.get("provenance")),
                    raw=a,
                )

        for d in data.get("derived", []):
            self.derived[d["id"]] = Derived(
                id=d["id"],
                type=d.get("type", ""),
                produced_by=d.get("produced_by"),
                derived_from=d.get("derived_from", []),
                output=d.get("output"),
                method=d.get("method"),
                raw=d,
            )

    # ---- loading ---------------------------------------------------------- #
    @classmethod
    def load(cls, path: Path | str = MANIFEST_PATH) -> "Manifest":
        return cls(json.loads(Path(path).read_text()))

    # ---- accessors -------------------------------------------------------- #
    def paper(self, pid: str) -> Paper:
        return self.papers[pid]

    def measurement(self, mid: str) -> Measurement:
        return self.measurements[mid]

    def tdb(self, tid: str) -> TdbArtifact:
        return self.tdbs[tid]

    def experimental_dataset(self, eid: str) -> ExperimentalDataset:
        return self.experimental[eid]

    def measurements_for(self, paper_id: str) -> list[Measurement]:
        return [m for m in self.measurements.values() if m.paper == paper_id]

    # ---- TDB fetch -------------------------------------------------------- #
    def resolve_tdb(self, tid: str, *, live: bool = True) -> Path:
        """Return the local path to a literature TDB, fetching it if missing.

        Uses the vendored tdbtools to live-download from the itemid/url recorded
        in the manifest into the artifact's `dest` (under the gitignored
        artifacts/ tree).
        """
        art = self.tdb(tid)
        dest = art.dest_path
        if dest.exists():
            return dest
        if not live:
            raise FileNotFoundError(f"{dest} missing and live fetch disabled")

        os.environ.setdefault("TDBDB_LIVE", "1")
        from tdbtools import tdbdb  # vendored, self-contained

        dest.parent.mkdir(parents=True, exist_ok=True)
        record = tdbdb.TdbdbRecord(
            itemid=art.itemid,
            elements=list(art.elements),
            authoryear=(art.provenance.citation if art.provenance else ""),
            tdb_url=art.fetch.get("url"),
        )
        out = Path(tdbdb.record_tdb_fixture(record, dest.parent))
        if out != dest and out.exists() and not dest.exists():
            out.replace(dest)
        if not dest.exists():
            raise RuntimeError(
                f"fetch of {tid} ({art.itemid}) did not produce {dest}; "
                "needs network access to NIMS TDBDB / Elsevier."
            )
        return dest


if __name__ == "__main__":
    m = Manifest.load()
    print(f"campaign: {m.campaign.get('name')}")
    print(f"papers: {list(m.papers)}")
    print(f"measurements: {list(m.measurements)}")
    print(f"tdbs: {list(m.tdbs)}  experimental: {list(m.experimental)}")
    print(f"derived: {list(m.derived)}")
