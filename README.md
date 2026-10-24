# greigite-fes-tdb

Reproducible computational basis for the **greigite (Fe₃S₄)** stability work, in three parts:

- **Thermodynamics** — CALPHAD databases for the Fe–S and Fe–S–O systems and the anhydrous
  predominance figures (`engine/`, `manuscript/`). Greigite enters as a stable phase built from the
  measured calorimetric enthalpy and entropy.
- **Porewater** — the aqueous Eh–pH stability diagram computed with PHREEQC (`aqueous/`), which
  places greigite's equilibrium field in reducing porewater along the lower stability limit of water
  from pH 2 to 12 and carries the reported calorimetric uncertainty (2 SD of the mean) as a band.
- **Kinetics** — the solid-state diffusion falsification (`kinetics/`): bulk sulfur diffusion cannot
  convert greigite at diagenetic temperature, so the observed conversion to pyrite runs through
  faster pathways that need an electron acceptor supplied by the porewater.

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
enthalpy variants at ± the reported uncertainty, 2 SD of the mean) to **`artifacts/tdb/fes_greigite_v1.tdb`**. The rest of
the campaign:

```bash
python engine/build_fes_o_tdb.py         # -> artifacts/tdb/fes_o_greigite_v1.tdb  (Fe–S–O)
python engine/build_boundary_tdbs.py     # -> artifacts/tdb/fes_greigite_boundary_{lower,upper}.tdb
python engine/validate_fes_engine.py     # predominance sweep -> artifacts/fes_engine_boundaries.json
python engine/validate_greigite.py       # parse + equilibrium sanity checks
```

`validate_fes_engine.py` is required before the figures: Figs. 1, 2 and S1 read the
boundary curves it writes.

## Figures

The manuscript figures read the built databases from `artifacts/tdb/` and write
PNGs to `artifacts/figures/`. Build the databases first (above), then regenerate
every manuscript figure with one command (Fig. 3 also needs `phreeqpython`; see
below):

```bash
python manuscript/make_all_figures.py    # -> artifacts/figures/Figure_*.png
```

Each figure also has its own script (resolving its own paths, so the working
directory doesn't matter), named for the manuscript figure it produces:

| Manuscript figure | Script | Output PNG |
|---|---|---|
| Fig. 1 & Fig. S1 | `manuscript/make_Figure_1_and_S1_fes.py` | `Figure_1.png`, `Figure_S1.png` |
| Fig. 2 | `manuscript/make_Figure_2_errorfield.py` | `Figure_2.png` |
| Fig. 3 | `aqueous/make_ehph_figure.py` (after `aqueous/build_ehph_diagram.py`) | `Figure_3.png` |
| Fig. 4 | `kinetics/make_kinetics_figure.py` | `Figure_4.png` |
| Fig. S2 | `manuscript/make_Figure_S2_feso_control.py` | `Figure_S2.png` |
| Fig. S3 | `manuscript/make_Figure_S3_cp.py` | `Figure_S3.png` |
| Fig. S4 | `manuscript/make_Figure_S4_feo.py` | `Figure_S4.png` |
| Fig. S5A/S5B | `manuscript/make_Figure_S5_predominance_300K.py` | `Figure_S5A.png`, `Figure_S5B.png` |
| Fig. S6A/S6B | `manuscript/make_Figure_S6_predominance_600K.py` | `Figure_S6A.png`, `Figure_S6B.png` |
| Figs. S5, S6 as printed (A above B) | `manuscript/stack_Figure_S5_S6.py` | `Figure_S5.png`, `Figure_S6.png` |
| Fig. S7A (2020 DSC, raw export) | `manuscript/make_Figure_S7_dsc.py` | `Figure_S7.png` |
| Fig. S8A (2026 DSC, raw export) | `manuscript/make_dsc_2026_remeasurement.py` | `dsc_2026_remeasurement.png` |

The anhydrous Fe–S–O diagrams (Figs. S2, S5) are computed at 300 K, the temperature of
the manuscript's boundary-fugacity tables. The uncertainty bounds (Fig. 2 and panels S5B, S6B) shift only
the greigite enthalpy, by its reported uncertainty (7.3 kJ per mole Fe, 2 SD of the mean;
Subramani et al. 2020). Pyrite and pyrrhotite are assessed database values and stay fixed. Figs. S7 and S8 are
measured data, not calculations: panel A of each is the DSC run (2020, and the same sample again in
2026) and panel B the powder XRD of the post-DSC product. The two DSC scripts re-plot the raw
instrument exports in `manuscript/data_dsc/` (heat flow in µV); the manuscript panels show the same
runs as plotted by the instrument software (calibrated mW). The XRD panels were produced in
external software from raw diffraction data not included in this repository.

One further script is kept for reference but is not in the current manuscript:
`manuscript/make_pipeline_diagram.py` (the build-pipeline schematic). The
`manuscript/explore_*.py` scripts are exploratory/superseded drafts.

## Porewater and kinetics

The aqueous and kinetic results sit alongside the CALPHAD databases (details in each folder's
README):

```bash
pip install phreeqpython numpy scipy matplotlib

python aqueous/derive_greigite_logk.py   # greigite dissolution log_K (-68.95; ±2 SD -65.12/-72.79)
python aqueous/build_ehph_diagram.py     # PHREEQC Eh-pH predominance -> artifacts/aqueous/
python aqueous/make_ehph_figure.py       # -> artifacts/figures/Figure_3.png
python kinetics/make_kinetics_figure.py  # -> artifacts/figures/Figure_4.png
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
| All-compound uncertainty-bound databases | `artifacts/tdb/fes_greigite_boundary_{lower,upper}.tdb` |
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
  (`data_dsc/`), and the manuscript companion note (`submission/`).
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
