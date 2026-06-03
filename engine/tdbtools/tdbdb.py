"""TDBDB — literature CALPHAD database index at Brown University.

Queries the public TDBDB index (https://avdwgroup.engin.brown.edu/) for free
literature TDB files. Returns structured records with author, year, paper DOI,
TDB download URL, and the elements covered.

This is the data source used for Phase 1.5 grafting (Path A) — finding
existing solution-phase backbones to merge with our synthesised line compounds.
Every record carries enough metadata to populate a `Provenance` block, so
downstream uses (in the synthesis tool, in publications) can attribute the
literature TDB properly.

Like the Materials Project data source, this module supports an offline
fixture mode: when `TDBDB_FIXTURES` is set (or `use_fixtures()` called),
queries read from `tests/fixtures/tdbdb_responses/*.json` instead of hitting
the live API.
"""

from __future__ import annotations

import json
import os
import shutil
import tempfile
import urllib.parse
import urllib.request
import zipfile
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

from . import provenance as prov_mod

TDBDB_BASE_URL = "https://avdwgroup.engin.brown.edu/getdbid.php"
USER_AGENT = "tdbdb-fetch/1.0 (research)"

# Artifact type used when archiving a tdbdb-fetched .tdb. Distinct from
# "tdb" (fitted/grafted databases) so list_available_tdbs() can
# bucket literature TDBs separately and apply the tdbdb:<itemid>
# reference scheme to them.
TDBDB_ARTIFACT_TYPE = "tdb_literature"

# Source name for the RunStore query_cache entry pointing at each
# fetched .tdb's archive path.
TDBDB_CACHE_SOURCE = "tdbdb_fetch"


@dataclass
class TdbdbRecord:
    """One literature TDB record returned by TDBDB."""

    itemid: str
    elements: list[str]
    authoryear: str
    paper_doi: str | None = None
    tdb_url: str | None = None
    is_kinetic: bool = False
    has_no_valid_tdb: bool = False

    @classmethod
    def from_raw(cls, raw: dict) -> TdbdbRecord:
        elements_raw = raw.get("element", [])
        if isinstance(elements_raw, str):
            # Sometimes returned as a string instead of a list
            try:
                elements_raw = json.loads(elements_raw)
            except json.JSONDecodeError:
                elements_raw = [elements_raw]

        paper_doi = raw.get("paperdoi") or None
        if paper_doi is not None and str(raw.get("hasdoi", "0")) != "1":
            paper_doi = None

        tdb_url = raw.get("tdburl") or None
        if tdb_url:
            tdb_url = tdb_url.replace(
                "TDBDB:", "https://avdwgroup.engin.brown.edu/tdb/"
            )

        return cls(
            itemid=str(raw.get("itemid", "")),
            elements=list(elements_raw)
            if isinstance(elements_raw, (list, tuple))
            else [],
            authoryear=str(raw.get("authoryear", "")),
            paper_doi=paper_doi,
            tdb_url=tdb_url,
            is_kinetic=str(raw.get("kinetic", "0")) == "1",
            has_no_valid_tdb=str(raw.get("novalidtdb", "0")) == "1",
        )

    def to_provenance(self) -> prov_mod.Provenance:
        return prov_mod.from_literature_tdb(
            tdbdb_itemid=self.itemid,
            elements=self.elements,
            authoryear=self.authoryear,
            paper_doi=self.paper_doi,
            tdb_url=self.tdb_url,
        )

    def to_dict(self) -> dict:
        return asdict(self)


# ---------------------------------------------------------------------------
# Fixture mode (mirrors data_sources/materials_project.py)
# ---------------------------------------------------------------------------


_FIXTURE_DIR: Path | None = None


def use_fixtures(path: Path | str | None) -> None:
    """Make TDBDB queries read from disk fixtures instead of the live API."""
    global _FIXTURE_DIR
    _FIXTURE_DIR = Path(path) if path is not None else None


def _maybe_use_env_fixtures() -> None:
    global _FIXTURE_DIR
    if os.environ.get("TDBDB_LIVE"):
        return
    if _FIXTURE_DIR is None:
        env_path = os.environ.get("TDBDB_FIXTURES")
        if env_path:
            _FIXTURE_DIR = Path(env_path)


def _chemsys_key(elements: list[str]) -> str:
    return "-".join(sorted({e.strip().capitalize() for e in elements}))


def _fixture_path(elements: list[str]) -> Path:
    assert _FIXTURE_DIR is not None
    return _FIXTURE_DIR / f"tdbdb_{_chemsys_key(elements)}.json"


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def search(elements: list[str]) -> list[TdbdbRecord]:
    """Query TDBDB for literature CALPHAD databases containing ALL of these elements.

    Returns parsed `TdbdbRecord` objects, each carrying enough metadata to
    populate a `Provenance` block via `.to_provenance()`. De-duplicated by
    itemid: TDBDB indexes multi-binary databases under several chemsys
    aliases (e.g. `nims_gasbzn_der` shows up for Ga-Sb, Ga-Zn, and Sb-Zn),
    and upstream search() would return the same record three times for a
    Ga-Zn query. Dedup preserves first-seen order so the ranking logic
    downstream stays stable.
    """
    _maybe_use_env_fixtures()

    if _FIXTURE_DIR is not None:
        path = _fixture_path(elements)
        if not path.exists():
            return []
        raw_list = json.loads(path.read_text())
    else:
        raw_list = _live_query(elements)

    records: list[TdbdbRecord] = []
    seen: set[str] = set()
    for raw in raw_list:
        rec = TdbdbRecord.from_raw(raw)
        if rec.itemid in seen:
            continue
        seen.add(rec.itemid)
        records.append(rec)
    return records


def _live_query(elements: list[str]) -> list[dict]:
    """Hit the live TDBDB endpoint and return the raw JSON record list."""
    normalised = [e.strip().capitalize() for e in elements]
    query = ",".join(normalised)
    url = f"{TDBDB_BASE_URL}?element={urllib.parse.quote(query)}"
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    # TDBDB_BASE_URL is a literal HTTPS endpoint to tdbdb.nims.go.jp;
    # user input is URL-encoded into the query string only.
    with urllib.request.urlopen(req, timeout=15) as resp:  # nosec B310
        raw = resp.read().decode("utf-8").strip()

    # TDBDB sometimes wraps responses in parentheses
    if raw.startswith("(") and raw.endswith(")"):
        raw = raw[1:-1]
    return json.loads(raw)


def record_fixture(elements: list[str], output_dir: Path | str) -> Path:
    """Hit the live API and write the raw response to a fixture file.

    Used by `tests/fixtures/record_tdbdb_fixtures.py`. Always bypasses fixture
    mode (we want fresh data).
    """
    raw_list = _live_query(elements)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"tdbdb_{_chemsys_key(elements)}.json"
    path.write_text(json.dumps(raw_list, indent=2))
    return path


# ---------------------------------------------------------------------------
# TDB file fetch (Phase 1.5 / Slice 12)
#
# `search()` returns metadata; `fetch()` actually downloads the .tdb file the
# record points at, with the same fixture-mode + cache-on-disk pattern that
# JANAF uses for its pickle cache. Live downloads are gated on
# `TDBDB_LIVE` exactly like `search()`.
# ---------------------------------------------------------------------------


_FIXTURE_TDB_DIR: Path | None = None
_BROWSER_UA = "Mozilla/5.0 (compatible; tdbdb-fetch/1.0)"


def use_tdb_fixtures(path: Path | str | None) -> None:
    """Make `fetch()` read .tdb files from a local fixture directory.

    Mirrors `use_fixtures()` (which controls `search()`) but for the actual
    downloaded TDB files. Tests and the offline demo set this to
    `tests/fixtures/literature_tdbs/`.
    """
    global _FIXTURE_TDB_DIR
    _FIXTURE_TDB_DIR = Path(path) if path is not None else None


def _maybe_use_env_tdb_fixtures() -> None:
    global _FIXTURE_TDB_DIR
    if os.environ.get("TDBDB_LIVE"):
        return
    if _FIXTURE_TDB_DIR is None:
        env_path = os.environ.get("TDBDB_TDB_FIXTURES")
        if env_path:
            _FIXTURE_TDB_DIR = Path(env_path)


def _open_with_ua(url: str, timeout: float = 30.0):
    """urllib open() with a browser User-Agent. Narrowly scoped to this module.

    Mirrors janaf.py's `_NistUrllibShim` philosophy: never monkeypatch a
    process-global opener; just build a one-off Request object with the UA we
    need. Some upstream hosts (NIST, Elsevier) reject Python's default UA.
    """
    req = urllib.request.Request(url, headers={"User-Agent": _BROWSER_UA})
    # URL comes from the TDBDB record we just fetched; it's a tdbdb.nims.go.jp
    # redirect target (TDB download endpoint), not user-controlled input.
    return urllib.request.urlopen(req, timeout=timeout)  # nosec B310


def _fixture_tdb_path(itemid: str) -> Path | None:
    """Look up a fixture .tdb (or .zip) for an itemid; return None if missing."""
    assert _FIXTURE_TDB_DIR is not None
    direct = _FIXTURE_TDB_DIR / f"{itemid}.tdb"
    if direct.exists():
        return direct
    zipped = _FIXTURE_TDB_DIR / f"{itemid}.zip"
    if zipped.exists():
        return zipped
    return None


def _extract_tdb_from_zip(zip_path: Path, itemid: str, dest_dir: Path) -> Path:
    """Extract the first/best .tdb member from a zip into dest_dir.

    If multiple .tdb files exist, prefer one whose stem contains `itemid`.
    Some Elsevier supplementary zips name the database ``*.TDB.txt`` (e.g.
    Dilner & Selleby 2017, member ``FeCaOS.TDB.txt``), so match that too.
    """
    dest_dir.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path) as zf:
        candidates = [
            n for n in zf.namelist() if n.lower().endswith((".tdb", ".tdb.txt"))
        ]
        if not candidates:
            raise ValueError(f"No .tdb file inside {zip_path}")
        # Prefer a member whose stem contains the itemid
        candidates.sort(key=lambda n: (itemid.lower() not in Path(n).stem.lower(), n))
        member = candidates[0]
        out_path = dest_dir / f"{itemid}.tdb"
        with zf.open(member) as src, open(out_path, "wb") as dst:
            shutil.copyfileobj(src, dst)
        return out_path


# ---------------------------------------------------------------------------
# RunStore-backed cache (replaces the legacy ~/.cache/literature_tdbs layout)
#
# Layout per fetched record:
#   - The .tdb file lives at the local artifact cache
#     (content-addressed; deduped across runs).
#   - A row in `run_artifacts` (artifact_type=TDBDB_ARTIFACT_TYPE) ties the
#     archived file to a singleton "tdbdb_fetch" run for this process.
#   - A row in `query_cache` (source=TDBDB_CACHE_SOURCE, query_key=itemid)
#     stores the archive_path so callers can resolve `tdbdb:<itemid>` ->
#     archive path without rescanning run_artifacts.
# ---------------------------------------------------------------------------


_TDBDB_RUN_ID: str | None = None
_TDB_PATH_CACHE: dict[str, Path] = {}


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def _get_or_create_tdbdb_run() -> str:
    """Lazily start a singleton tdbdb_fetch run for this process and return its id."""
    global _TDBDB_RUN_ID
    if _TDBDB_RUN_ID is None:
        # The RunStore-backed cache is not part of this self-contained build.
        raise RuntimeError(
            "The RunStore-backed TDB cache is unavailable in this "
            "self-contained build. Use fetch(record, cache_dir=...) or "
            "record_tdb_fixture(record, output_dir) for the directory cache."
        )
    return _TDBDB_RUN_ID


def cached_path_for(itemid: str) -> Path | None:
    """Resolve a tdbdb itemid to its cached archive path, if any.

    Returns ``None`` if no cache entry exists or the archived file no longer
    lives at the recorded path. Used by ``_formatting._resolve_tdbdb`` to
    serve ``tdbdb:<itemid>`` references via the RunStore-backed cache.
    """
    if itemid in _TDB_PATH_CACHE:
        path = _TDB_PATH_CACHE[itemid]
        if path.exists():
            return path
        del _TDB_PATH_CACHE[itemid]

    # The persistent (RunStore-backed) cache tier is not part of this
    # self-contained build; only the in-memory tier above is available.
    return None


def _record_tdb_in_runstore(record: TdbdbRecord, source_path: Path) -> Path:
    """Not available in the self-contained build (see module docstring)."""
    raise RuntimeError(
        "The RunStore-backed TDB archive is unavailable in this self-contained "
        "build. Use fetch(record, cache_dir=...) or record_tdb_fixture(record, "
        "output_dir), which write the .tdb directly under a directory you pass."
    )


def fetch(record: TdbdbRecord, *, cache_dir: Path | None = None) -> Path:
    """Download (or read from fixtures) the .tdb file referenced by `record`.

    Returns a local path to a usable .tdb file. Behaviour:

    1. If TDB fixture mode is active, look up `<fixture_dir>/<itemid>.tdb` (or
       `<itemid>.zip`); raise `FileNotFoundError` if absent. This is what
       tests + the offline demo run on.
    2. If ``cache_dir`` is explicitly passed, route the fetch through the
       legacy directory-based layout (``cache_dir/<itemid>.tdb``). This is
       used by ``record_tdb_fixture`` to populate ``tests/fixtures/literature_tdbs/``;
       normal callers should leave ``cache_dir`` unset.
    3. Otherwise consult the RunStore-backed cache (in-memory tier ->
       ``query_cache``); return the archived file if present.
    4. Otherwise (live mode only — gated on ``TDBDB_LIVE``) GET
       ``record.tdb_url`` with a browser UA, unpack zips if needed, archive
       the .tdb via ``record_artifact``, and stamp a ``query_cache`` pointer
       row.
    """
    _maybe_use_env_tdb_fixtures()

    if not record.itemid:
        raise ValueError("TdbdbRecord has no itemid; cannot fetch")
    if record.has_no_valid_tdb:
        raise ValueError(
            f"TDBDB record {record.itemid!r} is flagged as having no valid TDB"
        )

    # Fixture mode — never hit the network
    if _FIXTURE_TDB_DIR is not None:
        fix = _fixture_tdb_path(record.itemid)
        if fix is None:
            raise FileNotFoundError(
                f"No fixture TDB for itemid={record.itemid!r} under {_FIXTURE_TDB_DIR}"
            )
        if fix.suffix.lower() == ".zip":
            return _extract_tdb_from_zip(fix, record.itemid, _FIXTURE_TDB_DIR)
        return fix

    # Legacy directory mode (only when callers explicitly opt in via cache_dir).
    # record_tdb_fixture uses this to populate the on-disk fixture corpus.
    if cache_dir is not None:
        return _fetch_to_dir(record, Path(cache_dir))

    # RunStore-backed cache hit (in-memory -> query_cache).
    cached = cached_path_for(record.itemid)
    if cached is not None:
        return cached

    if not os.environ.get("TDBDB_LIVE"):
        raise RuntimeError(
            f"TDB fetch for {record.itemid!r} requires TDBDB_LIVE=1 "
            f"(or set fixtures via use_tdb_fixtures())"
        )
    if not record.tdb_url:
        raise ValueError(f"TDBDB record {record.itemid!r} has no tdb_url")

    # Stream the URL to a temp file; we don't know yet whether it's .zip or .tdb
    with _open_with_ua(record.tdb_url) as resp:
        ctype = resp.headers.get("Content-Type", "").lower()
        looks_zipped = record.tdb_url.lower().endswith(".zip") or "zip" in ctype
        suffix = ".zip" if looks_zipped else ".tdb"
        with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as tmp:
            shutil.copyfileobj(resp, tmp)
        tmp_path = Path(tmp.name)

    try:
        if looks_zipped:
            extract_dir = Path(tempfile.mkdtemp())
            try:
                tdb_path = _extract_tdb_from_zip(tmp_path, record.itemid, extract_dir)
                return _record_tdb_in_runstore(record, tdb_path)
            finally:
                shutil.rmtree(extract_dir, ignore_errors=True)
        # Plain .tdb body — rename to <itemid>.tdb so record_artifact archives
        # it under that filename inside the local artifact cache.
        renamed = tmp_path.parent / f"{record.itemid}.tdb"
        shutil.move(str(tmp_path), str(renamed))
        tmp_path = renamed
        return _record_tdb_in_runstore(record, tmp_path)
    finally:
        tmp_path.unlink(missing_ok=True)


def _fetch_to_dir(record: TdbdbRecord, cache_dir: Path) -> Path:
    """Legacy directory-mode fetch used by ``record_tdb_fixture``.

    Unlike the RunStore-backed path, this writes directly under ``cache_dir``
    so the fixture-recording script can populate ``tests/fixtures/literature_tdbs/``
    with a flat ``<itemid>.tdb`` layout that the test suite expects.
    """
    cached = cache_dir / f"{record.itemid}.tdb"
    if cached.exists():
        return cached
    if not os.environ.get("TDBDB_LIVE"):
        raise RuntimeError(
            f"TDB fetch for {record.itemid!r} requires TDBDB_LIVE=1 "
            f"(or set fixtures via use_tdb_fixtures())"
        )
    if not record.tdb_url:
        raise ValueError(f"TDBDB record {record.itemid!r} has no tdb_url")

    try:
        cache_dir.mkdir(parents=True, exist_ok=True)
        cache_writable = True
    except OSError:
        cache_writable = False

    with _open_with_ua(record.tdb_url) as resp:
        ctype = resp.headers.get("Content-Type", "").lower()
        looks_zipped = record.tdb_url.lower().endswith(".zip") or "zip" in ctype
        suffix = ".zip" if looks_zipped else ".tdb"
        with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as tmp:
            shutil.copyfileobj(resp, tmp)
        tmp_path = Path(tmp.name)

    if looks_zipped:
        dest_dir = cache_dir if cache_writable else Path(tempfile.mkdtemp())
        try:
            out = _extract_tdb_from_zip(tmp_path, record.itemid, dest_dir)
        finally:
            tmp_path.unlink(missing_ok=True)
        return out

    if cache_writable:
        try:
            shutil.move(str(tmp_path), str(cached))
            return cached
        except OSError:
            pass
    return tmp_path


def record_tdb_fixture(record: TdbdbRecord, output_dir: Path | str) -> Path:
    """Live-download a TDB and copy it under `output_dir/<itemid>.tdb`.

    Used by the fixture-recording script when bootstrapping
    `tests/fixtures/literature_tdbs/`. Always bypasses fixture mode.
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    # Force a live fetch even if TDBDB_LIVE isn't set
    prev_live = os.environ.get("TDBDB_LIVE")
    prev_fix = _FIXTURE_TDB_DIR
    use_tdb_fixtures(None)
    os.environ["TDBDB_LIVE"] = "1"
    try:
        path = fetch(record, cache_dir=out_dir)
    finally:
        if prev_live is None:
            os.environ.pop("TDBDB_LIVE", None)
        else:
            os.environ["TDBDB_LIVE"] = prev_live
        use_tdb_fixtures(prev_fix)
    return path
