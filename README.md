# greigite-fes-tdb

Reproducible computational basis for the **greigite (Fe₃S₄)** stability work, in three parts:

- **Thermodynamics** — CALPHAD databases for the Fe–S and Fe–S–O systems and the anhydrous
  predominance figures (`engine/`, `manuscript/`). Greigite enters as a stable phase built from the
  measured calorimetric enthalpy and entropy.
- **Porewater** — the aqueous Eh–pH stability diagram computed with PHREEQC (`aqueous/`), which
  places greigite's equilibrium field in reducing, near-neutral to alkaline sediment porewater and
  carries the calorimetric uncertainty as a ±1σ band.
- **Kinetics** — the solid-state diffusion falsification (`kinetics/`): bulk sulfur diffusion cannot
  convert greigite at diagenetic temperature, so the observed conversion to pyrite is
  solution-mediated and gated by porewater oxidant supply.

No database or Word files are committed: a published CALPHAD base is fetched and the greigite
databases are built from a small set of cited measured values on first run; everything generated
lands in the gitignored `artifacts/` tree.

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
PNGs to `artifacts/figures/`. Build the databases first (above), then regenerate
every published figure with one command:

```bash
python manuscript/make_all_figures.py    # -> artifacts/figures/Figure_*.png
```

Each figure also has its own script (resolving its own paths, so the working
directory doesn't matter), named for the published figure it produces:

| Published figure | Script | Output PNG |
|---|---|---|
| Fig. 1 & Fig. 3 | `manuscript/make_Figure_1_and_3_validation.py` | `Figure_1.png`, `Figure_3.png` |
| Fig. 2 | `manuscript/make_Figure_2_feso_control.py` | `Figure_2.png` |
| Fig. 4 | `manuscript/make_Figure_4_errorfield.py` | `Figure_4.png` |
| Fig. 5A/5B | `manuscript/make_Figure_5_predominance.py` | `Figure_5A.png`, `Figure_5B.png` |
| Fig. 6A/6B | `manuscript/make_Figure_6_predominance_600K.py` | `Figure_6A.png`, `Figure_6B.png` |
| Fig. S1 | `manuscript/make_Figure_S1_pipeline.py` | `Figure_S1.png` |
| Fig. S2 | `manuscript/make_Figure_S2_cp.py` | `Figure_S2.png` |
| Fig. S3 | `manuscript/make_Figure_S3_feo.py` | `Figure_S3.png` |
| Fig. S4 | `manuscript/make_Figure_S4_dsc_2020.py` | `Figure_S4.png` |
| Fig. S6 | `manuscript/make_Figure_S6_dsc_2026.py` | `Figure_S6.png` |

Figures S5 and S7 (powder XRD of the post-DSC products) were produced in external
software from raw diffraction data not included in this repository. The
`manuscript/explore_*.py` scripts are exploratory/superseded drafts, not part of
the published figure set.

## Porewater and kinetics

The aqueous and kinetic results sit alongside the CALPHAD databases (details in each folder's
README):

```bash
pip install phreeqpython numpy scipy matplotlib

python aqueous/derive_greigite_logk.py   # greigite dissolution log_K (-68.95; ±1σ -65.12/-72.79)
python aqueous/build_ehph_diagram.py     # PHREEQC Eh-pH predominance -> artifacts/aqueous/
python aqueous/make_ehph_figure.py       # -> artifacts/figures/fig_ehph_greigite.png
python kinetics/make_kinetics_figure.py  # -> artifacts/figures/fig_kinetics_falsification.png
```

The greigite dissolution constant is derived from the same measured thermochemistry as the CALPHAD
databases, on a database-consistent aqueous reference, and cross-checked by reproducing the
reference pyrite constant to within 0.04 log units
([`aqueous/LOGK_DERIVATION.md`](aqueous/LOGK_DERIVATION.md)). The kinetic conversion times were
computed with the proprietary Odinzen equilibrium engine, which is not included; the tabulated
result is ([`kinetics/README.md`](kinetics/README.md)).

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
- **`aqueous/`** — PHREEQC porewater Eh–pH diagram and the greigite dissolution log_K derivation.
- **`kinetics/`** — tabulated solid-state diffusion conversion times and the falsification figure.
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

- **Code** (the `engine/`, `manuscript/`, `aqueous/` and `kinetics/` scripts) — MIT, see
  [`LICENSE`](LICENSE).
- **Data, figures, and manuscript** (e.g. the experimental data under
  `manuscript/data_dsc/`, the tabulated result `kinetics/kinetics_conversion_time.csv`, and
  generated figures) — Creative Commons Attribution 4.0 International (CC-BY-4.0), see
  [`manuscript/data_dsc/LICENSE-CC-BY-4.0.txt`](manuscript/data_dsc/LICENSE-CC-BY-4.0.txt).

Copyright 2026 Odinzen LLC (Michael Bustamante, Gabriel Bustamante) and the authors.
