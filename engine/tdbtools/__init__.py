"""Self-contained, vendored helpers for fetching literature CALPHAD TDBs.

This package was extracted from the internal assessment toolkit so this
repository is fully self-contained: the only third-party requirement for the
TDB-fetch path is the Python standard library. ``tdbdb`` queries the public
NIMS TDBDB service and downloads the referenced ``.tdb`` files; ``provenance``
provides the small dataclass used to stamp where each record came from.
"""
