"""Structured provenance for every datum that flows through the pipeline.

Every result from a data source (Materials Project, NIST-JANAF, pySIPFENN,
literature TDBs from TDBDB) carries a `Provenance` block recording:

  - WHERE the value came from (source + source-specific ID)
  - HOW to cite it (human-readable citation string + canonical URL)
  - HOW it was computed (functional, model version, table edition, etc.)
  - WHEN we accessed it

The synthesis tool aggregates these into a `provenance.json` manifest emitted
alongside every dataset directory and TDB it produces, so a downstream user
(or a paper reviewer) can audit every formation energy in a fitted database
back to its original source.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def _now_iso() -> str:
    """ISO 8601 UTC timestamp, second precision."""
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


@dataclass
class Provenance:
    """Where a single datum came from and how to cite it.

    Stored alongside every synthesised phase and emitted as part of the
    manifest from `build_dataset_directory`. Designed to be machine-readable
    (for the admin console / agent) AND human-readable (for citation in a
    materials science paper).
    """

    source: str  # "materials_project" | "nist_janaf" | "pysipfenn" | "literature_tdb"
    source_id: str  # source-specific unique identifier
    citation: str  # human-readable citation, BibTeX-friendly when possible
    url: str | None = None  # canonical URL where the data can be re-fetched
    method: dict[str, Any] = field(default_factory=dict)  # source-specific method info
    accessed_at: str = field(default_factory=_now_iso)

    def to_dict(self) -> dict:
        return asdict(self)


# ---------------------------------------------------------------------------
# Citation strings for the data sources used here.
#
# These are the canonical references that should appear in any paper or report
# built from the pipeline's output. Per-datum source IDs (e.g. mp-942733-GGA) get
# substituted into the source_id field at construction time.
# ---------------------------------------------------------------------------


CITATION_MATERIALS_PROJECT = (
    "Jain, A. et al. The Materials Project: A materials genome approach to "
    "accelerating materials innovation. APL Materials 1, 011002 (2013). "
    "doi:10.1063/1.4812323"
)

CITATION_PYMATGEN = (
    "Ong, S. P. et al. Python Materials Genomics (pymatgen): A robust, open-source "
    "python library for materials analysis. Comput. Mater. Sci. 68, 314–319 (2013). "
    "doi:10.1016/j.commatsci.2012.10.028"
)

CITATION_NIST_JANAF = (
    "Chase, M. W., Jr. NIST-JANAF Thermochemical Tables, Fourth Edition. "
    "J. Phys. Chem. Ref. Data Monograph 9 (1998). https://janaf.nist.gov/"
)

CITATION_PYSIPFENN = (
    "Krajewski, A. M., Siegel, J. W., Xu, J. & Liu, Z.-K. Extensible Structure-Informed "
    "Prediction of Formation Energy with improved accuracy and usability employing neural "
    "networks. Comput. Mater. Sci. 208, 111254 (2022). doi:10.1016/j.commatsci.2022.111254"
)

CITATION_PYSIPFENN_MODELS = (
    "pySIPFENN pretrained model weights, Zenodo. doi:10.5281/zenodo.7373089"
)

CITATION_TDBDB = (
    "TDBDB: a free, open-access database of thermodynamic data files. "
    "Brown University, https://avdwgroup.engin.brown.edu/"
)

CITATION_HALLSTEDT_1992 = (
    "Hallstedt, B. Thermodynamic assessment of the system MgO-Al2O3. "
    "Calphad 16, 53-61 (1992). doi:10.1016/0364-5916(92)90015-Q | "
    "GHSER(1/2 H2,g) coefficients per Dinsdale, A. T. SGTE data for pure "
    "elements. Calphad 15, 317-425 (1991). doi:10.1016/0364-5916(91)90030-N"
)


# ---------------------------------------------------------------------------
# Convenience builders for each data source.
# ---------------------------------------------------------------------------


def from_materials_project(
    material_id: str, functional: str | None = None
) -> Provenance:
    """Provenance for a single Materials Project entry.

    Args:
        material_id: MP material ID (e.g. "mp-942733-GGA"). The trailing
            "-GGA" / "-GGA+U" / "-SCAN" suffix encodes the DFT functional.
        functional: Override the functional string (otherwise inferred from
            the material_id suffix).
    """
    if functional is None and "-" in material_id:
        # Extract everything after the last "-" as the functional code
        functional = material_id.rsplit("-", 1)[-1]

    return Provenance(
        source="materials_project",
        source_id=material_id,
        citation=CITATION_MATERIALS_PROJECT + " | retrieved via " + CITATION_PYMATGEN,
        url=f"https://materialsproject.org/materials/{material_id.split('-GGA')[0]}",
        method={"functional": functional} if functional else {},
    )


def from_nist_janaf(formula: str, phase_code: str, T: float) -> Provenance:
    """Provenance for a NIST-JANAF tabulated value at a specific temperature."""
    return Provenance(
        source="nist_janaf",
        source_id=f"janaf:{formula}({phase_code})@{T}K",
        citation=CITATION_NIST_JANAF,
        url=f"https://janaf.nist.gov/tables/{formula}-{phase_code}.html",
        method={"phase_code": phase_code, "temperature_K": float(T)},
    )


def from_pysipfenn(model: str, descriptor: str, formula: str) -> Provenance:
    """Provenance for a pySIPFENN ML formation-energy prediction."""
    return Provenance(
        source="pysipfenn",
        source_id=f"sipfenn:{model}:{formula}",
        citation=CITATION_PYSIPFENN + " | model weights: " + CITATION_PYSIPFENN_MODELS,
        url="https://doi.org/10.5281/zenodo.7373089",
        method={"model": model, "descriptor": descriptor},
    )


def from_hallstedt_1992() -> Provenance:
    """Provenance for the SGTE/Dinsdale 1991 hydrogen unary lattice stability.

    This entry is injected into ESPEI's SGTE91 reference tables at fit time
    so that H-bearing phases (hydroxides, hydrides) can be parameter-selected.
    The numerical GHSERHH coefficients come from Dinsdale 1991; we attribute
    the H-O assessment context to Hallstedt 1992 per the project convention.
    """
    return Provenance(
        source="dinsdale_1991_sgte",
        source_id="sgte91:H:1/2_MOLE_H2(G)",
        citation=CITATION_HALLSTEDT_1992,
        url="https://doi.org/10.1016/0364-5916(91)90030-N",
        method={
            "reference_phase": "1/2_MOLE_H2(G)",
            "ref_state": "SGTE91",
            "valid_range_K": [298.15, 6000.0],
        },
    )


def from_paper(
    *,
    paper_id: int,
    sha256: str,
    title: str | None,
    authors: list[str],
    year: int | None,
    doi: str | None = None,
    journal: str | None = None,
    page_number: int | None = None,
    table_or_figure: str | None = None,
    confidence: str = "medium",
    extraction_notes: str = "",
) -> Provenance:
    """Provenance for a thermodynamic data point extracted from a research paper.

    The paper lives in the local `PaperLibrary` (slice 14). The agent reads the
    paper, identifies the data point, and produces an `ExtractedDataPoint` —
    this builder turns the paper's library metadata into a Provenance block
    that travels with the resulting ESPEI dataset entry.
    """
    first_author = authors[0] if authors else "Unknown"
    parts: list[str] = []
    if first_author and len(authors) > 1:
        parts.append(f"{first_author} et al.")
    else:
        parts.append(first_author)
    if year is not None:
        parts.append(f"({year})")
    citation = " ".join(parts)
    if title:
        citation += f". {title}"
    if journal:
        citation += f". {journal}"
    if doi:
        citation += f". doi:{doi}"

    return Provenance(
        source="paper",
        source_id=f"paper:{paper_id}:sha256:{sha256[:12]}",
        citation=citation,
        url=f"https://doi.org/{doi}" if doi else None,
        method={
            "paper_id": int(paper_id),
            "sha256": sha256,
            "page_number": page_number,
            "table_or_figure": table_or_figure,
            "confidence": confidence,
            "extraction_notes": extraction_notes,
        },
    )


def from_literature_calorimetry(
    *,
    reference: str,
    components: list[str],
    phase: str,
    paper_id: int | None = None,
    doi: str | None = None,
) -> Provenance:
    """Provenance for mixing data extracted from published calorimetry.

    Used by the B1 (literature mining) sub-path of solution-phase synthesis.
    The data originates from research papers reporting enthalpy-of-mixing
    measurements (drop calorimetry, solution calorimetry, etc.).

    Args:
        reference: Human-readable citation string (e.g. "Sommer 1982").
        components: Element list for the mixing system (e.g. ["CU", "MG"]).
        phase: CALPHAD phase name (e.g. "LIQUID").
        paper_id: Optional PaperLibrary paper ID if the data was extracted
            via the paper ingestion pipeline.
        doi: Optional DOI for the source publication.
    """
    source_id = f"calorimetry:{'-'.join(sorted(c.upper() for c in components))}:{phase}"
    if paper_id is not None:
        source_id += f":paper:{paper_id}"
    return Provenance(
        source="literature_calorimetry",
        source_id=source_id,
        citation=reference,
        url=f"https://doi.org/{doi}" if doi else None,
        method={
            "extraction_type": "literature_calorimetry",
            "components": [c.upper() for c in components],
            "phase": phase,
            "paper_id": paper_id,
        },
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


# ---------------------------------------------------------------------------
# Manifest writer
# ---------------------------------------------------------------------------


@dataclass
class ProvenanceManifest:
    """Aggregated provenance for a synthesised dataset directory or TDB.

    A manifest collects every Provenance entry that contributed to a single
    output (a dataset directory, a fitted TDB, or a phase diagram), plus
    summary statistics about the source mix.
    """

    target: str  # path or label of what this manifest describes
    chemsys: list[str]  # element list the synthesis was run for
    entries: list[dict]  # list of {"phase": ..., "provenance": ...}
    generated_at: str = field(default_factory=_now_iso)

    @classmethod
    def from_phases(cls, target: str, chemsys: list[str], phases) -> ProvenanceManifest:
        entries = []
        for p in phases:
            prov = getattr(p, "provenance", None)
            entry = {
                "phase": getattr(p, "phase_name", None) or getattr(p, "name", "?"),
                "formula": getattr(p, "formula", None),
                "source": getattr(prov, "source", None) if prov else None,
                "provenance": prov.to_dict() if prov else None,
            }
            # Per-phase Cp source tag (set by espei_datasets when JANAF
            # provides a temperature-dependent Gibbs polynomial; otherwise
            # "neumann_kopp" to flag a constant-ΔH_f-only line compound).
            cp_source = getattr(p, "cp_source", None)
            if cp_source is not None:
                entry["cp_source"] = cp_source
            cp_fit = getattr(p, "cp_fit", None)
            if cp_fit is not None and hasattr(cp_fit, "to_dict"):
                entry["cp_fit"] = cp_fit.to_dict()
            # Slice 07: ThermoModel type and fit diagnostics
            thermo_model = getattr(p, "thermo_model", None)
            if thermo_model is not None and hasattr(thermo_model, "thermo_type"):
                entry["thermo_type"] = thermo_model.thermo_type()
            thermo_fit = getattr(p, "thermo_model_fit_result", None)
            if thermo_fit is not None and hasattr(thermo_fit, "diagnostics"):
                entry["fit_diagnostics"] = thermo_fit.diagnostics.to_dict()
            entries.append(entry)
        return cls(target=str(target), chemsys=list(chemsys), entries=entries)

    def write(self, path: str | Path) -> Path:
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(asdict(self), indent=2))
        return p

    def summary_lines(self) -> list[str]:
        """Human-readable summary suitable for printing or putting in a tool response."""
        by_source: dict[str, int] = {}
        for e in self.entries:
            src = e.get("source") or "unknown"
            by_source[src] = by_source.get(src, 0) + 1
        lines = [
            f"Provenance manifest: {self.target}",
            f"  Chemsys:    {'-'.join(self.chemsys)}",
            f"  Total phases: {len(self.entries)}",
        ]
        for src, n in sorted(by_source.items()):
            lines.append(f"    {src}: {n}")
        return lines
