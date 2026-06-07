# Fe–S(–O) TDB provenance — greigite

Where the CALPHAD inputs come from and how the greigite TDBs are built from them.

## Base CALPHAD databases
The Fe–S base is the **Dilner, Mao & Selleby (2015) Fe–Mn–S** assessment; the
Fe–S–O base is **Dilner & Selleby (2017) Ca–Fe–O–S**. Both are fetched from
**TDBDB** (Brown University's literature-TDB index, `avdwgroup.engin.brown.edu`).
Their Fe–S binary edge traces to **Lee et al. (ISIJ)**, with a partially-ionic
two-sublattice liquid (`IONIC_LIQ`). The Waldner & Pelton (2005) Fe–S assessment
is not in TDBDB and is not used here.

## How the files were obtained
The TDBs are downloaded straight from TDB-DB by the self-contained fetch helper
(`fetch_fes_tdbs.py`, using the vendored `tdbtools` package — stdlib only):
`tdbdb.search(["Fe","S"])` then `tdbdb.record_tdb_fixture(record, artifacts/tdb/)`
with `TDBDB_LIVE=1`. That path (`fetch(..., cache_dir=...)` → `_fetch_to_dir`) does
a live HTTP GET + zip extract into a flat directory — no database or service.
The base TDBs land in the gitignored `artifacts/tdb/` on first run (the builders
auto-fetch if absent). pycalphad version **0.11.1**; environment setup is in
`reproduce.md`.

## TDB-DB Fe–S records used

| itemid | reference (TDB-DB) | elements | source URL | outcome |
|---|---|---|---|---|
| `calphadj_1-s2.0-S0364591614000960-mmc1_1` | Dilner, Mao & Selleby (2015) | Fe-Mn-S | http://ars.els-cdn.com/content/image/1-s2.0-S0364591614000960-mmc1.zip | **fetched + parses** → `artifacts/tdb/calphadj_1-s2.0-S0364591614000960-mmc1_1.tdb` (15 502 B) — the **base** for the greigite build |
| `calphadj_1-s2.0-S0364591616300256-mmc1_1` | Dilner (2016) | Ca-Fe-Mg-Mn-S | http://ars.els-cdn.com/content/image/1-s2.0-S0364591616300256-mmc1.zip | **fetched + parses** → `artifacts/tdb/…256….tdb` (24 925 B); fetched during exploration, **not used** |
| `calphadj_1-s2.0-S0364591616301584-mmc1_1` | Dilner & Selleby (2017) | Ca-Fe-O-S | http://ars.els-cdn.com/content/image/1-s2.0-S0364591616301584-mmc1.zip | the zip member is `FeCaOS.TDB.txt` (44 674 B), not `*.tdb`; the `tdbtools` matcher accepts `.tdb.txt`, so `build_fes_o_tdb.py` fetches it and builds the **Fe–S–O** TDB on top |

## Base TDB source & authorship (the file our build actually uses)
Our greigite build reads exactly one base TDB:
**`artifacts/tdb/calphadj_1-s2.0-S0364591614000960-mmc1_1.tdb`** (set as `BASE`
in `build_greigite_tdb.py`, auto-fetched if missing). The companion
`…S0364591616300256…` (Ca‑Fe‑Mg‑Mn‑S) was fetched during exploration but is
**not** used. The Fe–S–O build (`build_fes_o_tdb.py`) instead uses the 2017
`…S0364591616301584…` Ca‑Fe‑O‑S file as its base.

What the `.zip` archives reveal (fetched + read from `ars.els-cdn.com`):
- Our base zip `…S0364591614000960-mmc1.zip` contains the single member
  **`mmc1/FeMnS.TDB`** (15 502 B) — the **Fe–Mn–S** assessment.
- The other zip `…S0364591616300256-mmc1.zip` contains member `mmc1.TDB` (24 925 B).
- Both TDB headers read `Database file written 2014-2-14` / `From database: SSUB4`
  — i.e. exported from Thermo‑Calc against SGTE's SSUB4 substances DB. The file's
  own `REF` block names only its data sources (Dinsdale SGTE91; Lee–Sundman–Kim–Chin
  ISIJ for Fe–S; Huang for Fe–Mn; SSUB5 gas) — it does **not** self-name its authors.

**Author of our base TDB: Daniella Dilner, Huahai Mao, and Malin Selleby**
(KTH Royal Institute of Technology, Stockholm), Fe–Mn–S assessment, *Calphad* (2015).
Attribution comes from two independent sources (NOT from reading the paywalled
article): (a) the TDB‑DB harvest metadata labels the PII "Dilner, Mao and Selleby
(2015)"; (b) the companion 2016 file's `REF5` explicitly cites
*"D. Dilner, H. Mao, M. Selleby Calphad … 95–105; Fe‑Mn‑S."* The other file is
**Dilner (2016)**, Ca‑Fe‑Mg‑Mn‑S.

URLs:
- **Open (fetched + read here):** the Elsevier supplementary archives —
  `http://ars.els-cdn.com/content/image/1-s2.0-S0364591614000960-mmc1.zip` (ours)
  and `…1-s2.0-S0364591616300256-mmc1.zip`.
- **Paywalled (ScienceDirect article pages, deterministic PII pattern — NOT
  opened/read):** `https://www.sciencedirect.com/science/article/pii/S0364591614000960`
  and `…/pii/S0364591616300256`. (No DOI quoted — not reliably derivable from the
  PII without a lookup, and the articles were not read.)

For the greigite **corpus papers** (Subramani 389, Shumway 384, Waldner–Pelton 388)
and their full metadata table, see `reproduce.md`.

## Key extracted values (from the two parsed Dilner files; identical Fe–S edge)

Pyrrhotite phase — **non-stoichiometric Fe(1-x)S sublattice CEF**:
`PHASE PYRRHOTITE % 2 1 1` → model `(FE,MN,VA)_1 : (S)_1`. Fe/Va mixing on the
metal sublattice gives the Fe-deficiency. Parameters (REF2 = Lee et al., ISIJ):
```
G(PYRRHOTITE,FE:S;0)   = +GFES# +GHSERFE# +GHSERSS#
G(PYRRHOTITE,VA:S;0)   = +GHSERSS# +258600
G(PYRRHOTITE,FE,VA:S;0)= -407000 +10*T
G(PYRRHOTITE,FE,VA:S;1)= +60000 +20*T
FUNCTION GFES = -107518 -18.19*T +1.78*T*LN(T)   (298.15–6000 K)   → ΔfH(FeS) ≈ -107.518 kJ/mol
```

Pyrite (FeS2), `PHASE PYRITE % 2 1 2` → `(FE,MN):(S)`:
```
G(PYRITE,FE:S;0) = +GHSERFE# +2*GHSERSS# -177763 +48.567*T   (REF2 = Lee et al., ISIJ)
   → enthalpy a-coefficient (≈ ΔfH298, SER-referenced) = -177.763 kJ/mol
```
Compared to the stand-in: pyrite **-177.763** vs stand-in -171.544 / -167.36
kJ/mol → Dilner/Lee pyrite is ~6.2 / ~10.4 kJ **more negative**.
Pyrrhotite basis here = Lee/Dilner `GFES` (NOT JANAF, NOT Grønvold–Stølen).

## Reference blocks (raw, from the TDB text)
2015 (Fe-Mn-S): REF1 Dinsdale SGTE91; **REF2 = B.-J. Lee, B. Sundman, S. Kim,
K.-W. Chin, ISIJ International** (the Fe–S source); REF3 Huang Fe-Mn; REF4 SSUB5 gas.
2016 (Ca-Fe-Mg-Mn-S): adds REF5 Dilner, Mao, Selleby Calphad 2015; REF7 Selleby &
Sundman Ca-Fe-O; REF10 Dilner, Kjellqvist & Selleby; REF11 Hallstedt; REF12 Tibballs.
**No "Waldner" or "Pelton" string occurs in either file.**

## Parse check (pycalphad 0.11.1)
Both base TDBs load cleanly under `pycalphad.Database(...)`; their `IONIC_LIQ`
liquid model is parsed natively.

## ---------------------------------------------------------------------------
## Greigite graft (bolt greigite onto a working Fe(1-x)S base)
## ---------------------------------------------------------------------------

Rather than transcribe Waldner–Pelton, greigite is added to a base that already
parses and already has a non-stoichiometric Fe(1-x)S phase: Dilner 2015
(`artifacts/tdb/calphadj_1-s2.0-S0364591614000960-mmc1_1.tdb`).

### Greigite thermodynamic source — MEASURED data (Subramani 2020 + Shumway 2022)
Greigite thermodynamics come from two papers, recorded as typed `measurement`
artifacts in the manifest (each citing its paper + table + page):
- **Subramani 2020** (PNAS, Table 1, p.2) —
  ΔfH° = **−144.1 ± 7.3 kJ/mol per FeS1.33** (⇒ **−432.3 kJ/mol Fe3S4**).
  ±7.3 is **absolute** (kJ/mol, drop-solution calorimetry), confirmed against the
  Shumway 2022 Table 6 caption (p.5).
- **Shumway 2022** (J. Chem. Thermodyn.) —
  **S°(298.15) = 71.334 J/mol·K** (Table 5, p.4) and the **Debye–Einstein** high-T
  Cp (Table 4, p.3): m=0.9740, Θ_D=237.85 K, n=1.1329, Θ_E=405.26 K, A=0.0398.
- **Waldner & Pelton 2005** is **NOT** used; the Dilner 2015 assessment stays the
  base.

The build's audit chain runs through the campaign manifest (every value is a
typed `measurement` carrying its paper + table + page + a provenance block):
```
provenance_manifest.json  (measurement artifacts)  →  artifacts/tdb/fes_greigite_v1.tdb
```
- The committed manifest keeps the build turnkey: it carries the cited numeric
  values and their source locators.

### Authored TDB (greigite G = measured T-dependent line compound)
`artifacts/tdb/fes_greigite_v1.tdb` = Dilner 2015 base + a stoichiometric **(FE)3(S)4** GREIGITE
line compound whose Gibbs energy is built by `build_greigite_tdb.py` from the
manifest `measurement` artifacts (deterministic, idempotent within an environment;
across environments the lstsq fit coefficients move in the ~8th sig fig only):
```
G_greigite(T) = ΔfH°(298) + ∫₂₉₈ᵀ Cp dT − T·[ S°(298) + ∫₂₉₈ᵀ (Cp/T) dT ]
PARAMETER G(GREIGITE,FE:S;0) = +3*GHSERFE# +4*GHSERSS# +GREIGITE_GF#   (298.15–600 K)
```
- Cp = Shumway Debye–Einstein (×3 → Fe3S4). The **DE form** is used, NOT the cubic /
  mid-T polynomial (which diverges >300 K). Cp is fit to an SGTE polynomial over
  **298–530 K**; the SER remainder fit has **max residual 0.03 J/mol**.
- G is written **+3 GHSERFE +4 GHSERSS + remainder**, reusing Dilner's SER element
  references so greigite sits on the same scale as the base.
- **Extrapolation flag:** Cp 298–530 K extrapolates a fit valid to ~303 K, through
  the ~530 K kinetic decomposition of greigite; the TDB ceiling is 600 K. This
  shifts fugacity boundaries <0.1 log unit.
- Two ΔfH limit variants for the upper/lower Fig. 2 figures:
  `fes_greigite_v1_dHf_lo.tdb` (−432.3 − 21.9 kJ) and `_dHf_hi.tdb` (−432.3 + 21.9 kJ),
  where 21.9 = 3 × 7.3 kJ/mol (per Fe3S4).

### Validation (pycalphad 0.11.1)
- **Internal Cp/S checks** (`build.log`): DE Cp(298.15) = 59.687 vs paper 59.686
  (Δ 0.001); S(300) from ∫DE/T·0→300 = 72.15 vs paper 71.704 (+0.44, 0.6% — the DE
  form slightly overestimates sub-56 K entropy where Shumway used dedicated low-T
  fits; this does NOT enter the TDB, which anchors on the measured S°298 and
  integrates Cp only from 298 K up).
- **Parse:** all three TDB variants load; `GREIGITE` present as (FE)3(S)4.
- **Idempotent:** rerunning the build reproduces identical sha256.
- **Physics:** with the measured data greigite is the **stable** phase at its
  composition across 298–600 K (Subramani's central result).
- **Consistency on the Dilner basis** (`consistency_greigite.py`): the
  decomposition FeS1.33 → 0.726 FeS1.092 + 0.274 FeS2 using Dilner's
  pycalphad-computed ΔHf(pyrite)=−168.61, ΔHf(pyrrhotite,FeS1.092)=−105.19 kJ/mol
  gives a stability margin **+21.5 kJ/mol = +2.95σ** (σ = 7.3 kJ); stable across
  the whole ±7.3 band (+1.95σ to +3.95σ). **Supersedes the old JANAF-basis 0.27σ
  figure** (not reused). Note: the bare −177.763 pyrite a-coefficient is NOT ΔfH —
  the SGTE GHSER functions carry non-zero enthalpy at 298 K, so the true SER
  formation enthalpy (pycalphad) is −168.61 kJ/mol.

## ---------------------------------------------------------------------------
## Fe–S–O builder + all-compound boundary engine
## ---------------------------------------------------------------------------

### Fe–S–O greigite TDB (`build_fes_o_tdb.py` → `artifacts/tdb/fes_o_greigite_v1.tdb`)
The Fe–S–O predominance figures need a greigite-bearing Ca‑Fe‑O‑S TDB. It is
built mechanically (no new CALPHAD modelling) from the fetched Dilner‑2017 base:
1. **fetch** Dilner & Selleby 2017 (`…S0364591616301584…`; zip member
   `FeCaOS.TDB.txt`) into `artifacts/tdb/`.
2. **dedupe** → `fes_o_dilner2017_clean.tdb`: normalise CRLF→LF and drop the
   second (duplicate) occurrence of three FUNCTIONs pycalphad rejects
   (`GCAOSOL`, `GWUSTITE`, `GFES`). No thermodynamic values change.
3. **graft** → `fes_o_greigite_v1.tdb`: append the measured **GREIGITE** block
   (the same `GREIGITE_GF` function + GREIGITE phase produced by
   `build_greigite_tdb.py`, reused verbatim) and the Fe-only **PYRITE** block
   from the Dilner-2015 / Lee et al. assessment
   (`+GHSERFE# +2*GHSERSS# −177763 +48.567*T`).
Both transforms reproduce the known-good file byte-for-byte (`--verify`); the
result parses with `GREIGITE` + `PYRITE` present and oxide formation energies
within ~1% of the literature (FeO −251.9, Fe₃O₄ −334.7/Fe, Fe₂O₃ −367.3 kJ).

### All-compound boundary TDBs (`build_boundary_tdbs.py` → `artifacts/tdb/`)
Builds the LOWER/UPPER boundary-case TDBs for the ±1σ envelope by shifting the
Dilner pyrite (±6 kJ on the −177763 constant) and pyrrhotite (±3.5 kJ on `GFES`)
functions alongside the greigite ±σ variants — pure logged string surgery on the
already-validated greigite ±σ TDBs. `validate_fes_engine.py` runs the pycalphad
predominance sweep on the central/lo/hi variants and writes the boundary table to
`artifacts/fes_engine_boundaries.json` + a diagram to `artifacts/figures/`.

## Provenance manifest (the campaign's machine-readable index)
`provenance_manifest.json` is the single typed, provenance-tracked description of
the campaign; `manifest.py` is the read API. Scripts (and agents) drive the whole
build from this one file:
- `campaign` — name, description, chemical systems.
- `papers[]` — reference entities (DOI, authors, year) + a `kg` block
  (`corpus_id`, `tags`) for the metadata that's useful outside the DB.
- `artifacts[]` — typed nodes, each with a `provenance` block
  (`source`/`source_id`/`citation`/`url`/`method`):
  - `tdb_from_tdbdb` — a literature TDB to fetch (`itemid`, `url`, `dest`).
  - `measurement` — a measured quantity tied to a paper (value/unit/basis/
    uncertainty/conditions/method/locator; `model` holds fitted forms like the
    Debye-Einstein Cp). These ARE the numbers the build consumes.
  - `experimental_dataset` — collaborator data (with `author` + CRediT
    `contributor_role`), an in-repo `path`
    instead of a URL.
- `derived[]` — the recipe DAG: each output lists `produced_by` (script) +
  `derived_from` (artifact ids).

## Files in this folder
- `provenance_manifest.json` — the campaign manifest (above); **committed**, the source of the measured values.
- `manifest.py`             — typed read API over the manifest (Campaign / Measurement / Tdb / ExperimentalDataset + Provenance).
- `fetch_fes_tdbs.py`       — fetch the manifest's literature TDBs into `artifacts/tdb/`; writes `artifacts/fetch_report.json`.
- `inspect_fes_tdbs.py`     — pycalphad parse-check + phase/sublattice/reference dump of `artifacts/tdb/`.
- `extract_functions.py`    — full multi-line FUNCTION/PARAMETER blocks + 2017-zip peek.
- `build_greigite_tdb.py`   — measured-data builder (reads the manifest) → the 3 greigite TDB variants in `artifacts/tdb/` + `artifacts/build_report.json`.
- `build_fes_o_tdb.py`      — Fe–S–O builder (fetch + dedupe + graft) → `artifacts/tdb/fes_o_greigite_v1.tdb`.
- `build_boundary_tdbs.py`  — all-compound ±1σ boundary TDBs → `artifacts/tdb/`; report → `artifacts/boundary_cases_report.json`.
- `validate_fes_engine.py`  — predominance sweep → `artifacts/fes_engine_boundaries.json` + `artifacts/figures/fes_engine_diagram.png`.
- `consistency_greigite.py` — Dilner-basis greigite stability σ-check → `artifacts/consistency_report.json`.
- `validate_greigite.py`    — equilibrium / metastability-margin sanity check.
- `fes_greigite_suppressed_fields.json` — **committed** reference fields (no generator) the validation figures read.
- **Not committed** (built/fetched into `artifacts/` on first run): all `*.tdb`
  (base + greigite + Fe–S–O + boundary), figures, and the generated result JSONs
  (`fes_engine_boundaries.json`, `boundary_cases_report.json`,
  build/consistency/fetch reports).

### Reproduce from scratch
```
python engine/build_greigite_tdb.py      # auto-fetches the Dilner-2015 base (per the manifest)
python engine/build_fes_o_tdb.py         # fetches Dilner-2017; dedupe + graft
python engine/build_boundary_tdbs.py     # ±1σ boundary TDBs
python engine/validate_fes_engine.py     # predominance sweep + diagram
python engine/consistency_greigite.py    # Dilner-basis σ-distance
```
The build + checks run from the committed `provenance_manifest.json`.
