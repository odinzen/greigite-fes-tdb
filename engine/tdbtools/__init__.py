"""Self-contained, vendored helpers for fetching literature CALPHAD TDBs.

These helpers keep the repository self-contained: the only requirement for the
TDB-fetch path is the Python standard library. ``tdbdb`` queries the public
NIMS TDBDB service and downloads the referenced ``.tdb`` files; ``provenance``
provides the small dataclass used to stamp where each record came from.
"""
