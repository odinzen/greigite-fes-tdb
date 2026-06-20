# greigite-fes-tdb

Reproducible CALPHAD thermodynamic databases for **greigite (Fe₃S₄)** in the
Fe–S and Fe–S–O systems — the scripts that build them and the figures for the
accompanying manuscript.

No database files are committed: a published CALPHAD base is fetched and the
greigite databases are built from a small set of cited measured values on first
run.

## Quick start

Needs Python 3.12; the first build needs network access to fetch the base
database. (Full environment notes: [`engine/reproduce.md`](engine/reproduce.md).)

```bash
conda create -n greigite-fes-tdb -c conda-forge -y python=3.12 \
    "pycalphad=0.11.1" numpy scipy matplotlib
conda activate greigite-fes-tdb

python engine/build_greigite_tdb.py
```

This fetches the base TDB and writes the greigite **Fe–S** database (plus its
±1σ enthalpy variants) to **`artifacts/tdb/fes_greigite_v1.tdb`**. The rest of
the campaign:

```bash
python engine/build_fes_o_tdb.py         # -> artifacts/tdb/fes_o_greigite_v1.tdb  (Fe–S–O)
python engine/build_boundary_tdbs.py     # -> artifacts/tdb/fes_greigite_boundary_{lower,upper}.tdb
python engine/validate_greigite.py       # parse + equilibrium sanity checks
```

## Figures

The manuscript figures read the built databases from `artifacts/tdb/` and write
PNGs to `artifacts/figures/`. Build the databases first (above), then:

```bash
python engine/validate_fes_engine.py     # predominance sweep -> artifacts/figures/fes_engine_diagram.png
python manuscript/make_fig3_engine.py    # Fe–S–O single-database predominance
# the remaining manuscript/make_fig*.py and render_fig2b.py regenerate the rest
```

Each script resolves its own paths, so the working directory doesn't matter.

## What you get

Everything lands in the gitignored `artifacts/` tree:

| Output | Path |
|--------|------|
| Fe–S greigite database (+ `_dHf_lo` / `_dHf_hi` variants) | `artifacts/tdb/fes_greigite_v1.tdb` |
| Fe–S–O greigite database | `artifacts/tdb/fes_o_greigite_v1.tdb` |
| All-compound ±1σ boundary databases | `artifacts/tdb/fes_greigite_boundary_{lower,upper}.tdb` |
| Manuscript figures | `artifacts/figures/*.png` |

## How it works

The whole campaign is described by one machine-readable file,
[`engine/provenance_manifest.json`](engine/provenance_manifest.json): typed,
provenance-tracked artifacts — the literature TDB to fetch, the measured values
taken from the source papers (each with its citation, table, and page), the
experimental datasets — plus the recipe that turns them into the derived
databases and figures. The builders read it through
[`engine/manifest.py`](engine/manifest.py), so every number in a database traces
back to a cited source. Details: [`engine/provenance.md`](engine/provenance.md).

## Repository layout

- **`engine/`** — the provenance manifest + its reader, the TDB build/validation
  scripts, and the vendored TDBDB fetch helper.
- **`manuscript/`** — figure-generation scripts, DSC/XRD source data
  (`data_dsc/`), and the journal `submission/` bundle.
- **`artifacts/`** — all generated output (gitignored).

## Data sources & copyright

This repository fetches the source CALPHAD databases rather than including them
(for copyright and fairness reasons). The Fe–S and Ca–Fe–O–S base databases are
journal supplementary material, indexed by
[TDBDB](https://avdwgroup.engin.brown.edu/); the builders download them from the
publisher at build time, and our greigite/pyrite additions are grafted on
locally.

## License

Dual-licensed by component:

- **Code** (the `engine/` and `manuscript/` scripts) — MIT, see [`LICENSE`](LICENSE).
- **Data, figures, and manuscript** (e.g. the experimental data under
  `manuscript/data_dsc/`, and generated figures) — Creative Commons Attribution 4.0
  International (CC-BY-4.0), see
  [`manuscript/data_dsc/LICENSE-CC-BY-4.0.txt`](manuscript/data_dsc/LICENSE-CC-BY-4.0.txt).

Copyright 2026 Odinzen LLC (Michael Bustamante, Gabriel Bustamante) and the authors.
