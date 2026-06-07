# Reproducing the greigite Fe–S(–O) TDBs

How to rebuild the CALPHAD databases from the files in this repo. The base Dilner
TDBs are fetched, and every derived TDB is built, into the gitignored
**`artifacts/`** tree on first run. The whole campaign is described by a single
machine-readable **provenance manifest**
(`provenance_manifest.json`, read via `manifest.py`): typed, provenance-tracked
artifacts — the literature TDBs to fetch, the measured values taken from papers,
the collaborator experimental datasets — plus the recipe DAG that turns them
into the derived TDBs and figures. The builders read their cited measured values
from the manifest, so the headline greigite build is **turnkey without the
source PDFs**.

I/O layout:
- **`artifacts/`** — ALL generated output (gitignored): fetched + built TDBs in
  `artifacts/tdb/`, figures in `artifacts/figures/`, build/consistency/fetch
  reports in `artifacts/`.

## Prerequisites

- **Python 3.12**, via conda / miniforge (conda-forge channel). Any environment
  manager works; the commands below use conda.
- **Network access on the first build** — the builders fetch the base Dilner
  TDBs from NIMS TDBDB / Elsevier into `artifacts/tdb/`. Later runs are offline.
- **Python packages** (installed in step 1): `pycalphad=0.11.1`, `numpy`, `scipy`,
  `matplotlib` for the core build + figures; `reportlab` only for the optional
  `md_to_pdf.py` step. No other services or databases.

## 1. Environment (conda)

```bash
# miniforge / conda-forge community scientific tooling:
conda create -n greigite-fes-tdb -c conda-forge -y python=3.12 \
    "pycalphad=0.11.1" numpy scipy matplotlib reportlab
conda activate greigite-fes-tdb
```

| package    | used by                                          | notes |
|------------|--------------------------------------------------|-------|
| pycalphad  | build_greigite / build_fes_o / build_boundary / validate / consistency | `Database`, `equilibrium` (v0.11.1) |
| numpy      | builds + validation + figures                    | grids, least-squares fit |
| scipy      | build_greigite (`scipy.integrate.quad`)          | Debye/Einstein + entropy integrals |
| symengine  | build_greigite (`Symbol`)                        | ships with pycalphad; evaluates GHSER |
| matplotlib | validate_fes_engine + manuscript figures         | `Agg` backend |
| reportlab  | md_to_pdf.py                                      | optional .md → PDF |
| (stdlib)   | fetch (`tdbtools`)                                | no third-party deps |

The fetch step (`fetch_fes_tdbs.py` / the auto-fetch inside the builders) uses
the self-contained, stdlib-only `tdbtools` package vendored here plus network
access to NIMS TDBDB / Elsevier.

## 2. Core build path (fetch → build → validate → figures)

Run from the repo root (scripts resolve their own paths via `__file__`, so the
working directory does not matter). Each step fetches what it needs if missing.

```bash
conda activate greigite-fes-tdb

# (a) Fe–S greigite TDBs. Auto-fetches the Dilner-2015 base into artifacts/tdb/
#     if absent, reads the measured values from provenance_manifest.json, and
#     writes fes_greigite_v1.tdb + fes_greigite_v1_dHf_{lo,hi}.tdb to artifacts/tdb/
#     (+ artifacts/build.log, artifacts/build_report.json).
python engine/build_greigite_tdb.py

# (b) Fe–S–O greigite TDB. Fetches Dilner-2017 (Ca-Fe-O-S), dedupes duplicate
#     FUNCTIONs, grafts the measured GREIGITE block (reused from step a) + the
#     Fe-only PYRITE block -> artifacts/tdb/fes_o_greigite_v1.tdb.
python engine/build_fes_o_tdb.py

# (c) All-compound 1-sigma boundary TDBs (LOWER/UPPER) -> artifacts/tdb/
#     (+ writes artifacts/boundary_cases_report.json).
python engine/build_boundary_tdbs.py

# (d) checks (optional but recommended)
python engine/validate_fes_engine.py    # predominance sweep -> artifacts/fes_engine_boundaries.json
                                         #   + artifacts/figures/fes_engine_diagram.png
python engine/consistency_greigite.py    # Dilner-basis stability sigma-distance -> artifacts/
python engine/validate_greigite.py       # parse + equilibrium sanity

# (e) manuscript figures (read TDBs from artifacts/tdb/, engine JSONs from
#     engine/, write PNGs to artifacts/figures/)
python manuscript/make_fig4_feso_control_bw.py    # etc. — any make_fig*.py / render_fig2b.py
```

The greigite build is **deterministic and idempotent within an environment** —
rerunning reproduces the TDBs byte-for-byte (sha256 in `artifacts/build_report.json`).
Across environments the SGTE least-squares fit coefficients can differ in the
~8th significant figure (BLAS/scipy build differences); this is thermodynamically
negligible (fit residual stays 0.03 J/mol). `build_fes_o_tdb.py` is a mechanical
dedupe + graft: its STEP-1 dedupe is byte-identical anywhere, and STEP 2 grafts
the GREIGITE block from the freshly-built `fes_greigite_v1.tdb` (pass
`--verify <ref.tdb>` to byte-compare against a known-good file from the same env).
Expected checkpoints:

- DE `Cp(298.15) = 59.687` vs paper `59.686` (Δ 0.001)
- SER-remainder → SGTE fit `max|resid| = 0.03 J/mol` over 298–530 K
- all 3 greigite variants parse, `GREIGITE` present as `(FE)3(S)4`
- greigite is the **stable** phase at its composition across 298–600 K
- Dilner-basis stability margin **+21.5 kJ/mol = +2.95σ**

## 3. Where the measured values come from

The build reads its measured values from the `measurement` artifacts in
`provenance_manifest.json`. Each one carries a precise **source locator** (paper +
table + page) and a `provenance` block, so every number points back to an exact
spot in the primary source — a verification map, not a copy of the paper.

The base TDBs are fetched automatically by the builders; to fetch them explicitly:

```bash
python engine/fetch_fes_tdbs.py          # manifest TDBs -> artifacts/tdb/ + artifacts/fetch_report.json
```

## 4. Source papers used for the greigite values

| authors | year | title | journal | DOI | role in this build |
|---------|------|-------|---------|-----|--------------------|
| T. Subramani, K. Lilova, M. Abramchuk, K.D. Leinenweber, A. Navrotsky | 2020 | Greigite (Fe₃S₄) is thermodynamically stable: Implications for its terrestrial and planetary occurrence | Proc. Natl. Acad. Sci. (PNAS) | 10.1073/pnas.2017312117 | **ΔfH° = −144.1 ± 7.3 kJ/mol per FeS₁.₃₃** (Table 1) |
| S.G. Shumway, J. Wilson, K. Lilova, T. Subramani, A. Navrotsky, B.F. Woodfield | 2022 | The low-temperature heat capacity and thermodynamic properties of greigite (Fe₃S₄) | J. Chem. Thermodyn. | 10.1016/j.jct.2022.106836 | **S°₂₉₈ = 71.334 J/mol·K (Table 5) + Debye–Einstein Cp (Table 4)** |
| P. Waldner, A.D. Pelton | 2005 | Thermodynamic Modeling of the Fe-S System | J. Phase Equilib. Diffus. | 10.1361/15477030522455 | **NOT used** (noted only; Dilner 2015 stays the base) |

The build's audit chain runs through the manifest:
`provenance_manifest.json` (`measurement` artifacts, citing paper + table + page)
→ `artifacts/tdb/fes_greigite_v1.tdb`.

The **base** CALPHAD TDBs (Dilner, Mao & Selleby 2015, Fe–Mn–S; and Dilner &
Selleby 2017, Ca–Fe–O–S) are a separate provenance chain — see `provenance.md`.
