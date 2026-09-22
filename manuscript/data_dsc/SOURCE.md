README
# Bulk Fe3S4 (greigite) decomposition data — source provenance

Experimental data measured by Tamilarasan Subramani (co-author) — Investigation
(sample preparation and calorimetric/thermal measurements).

- `Fe3S4_bulk_heating.txt`  — DSC/TG, bulk greigite, heating to 600 C (SETARAM export, UTF-16).
- `Fe3S4_bulk_cooling.txt`  — DSC/TG, bulk greigite, cooling run.
  Columns: Index, Time (s), Sample Temperature (C), TG (mg, baseline-corrected), HeatFlow (uV).
- `Fe3S4_bulk_XRD_postDSC600.png` — powder XRD of bulk greigite after DSC to 600 C;
  product is pyrrhotite-3T (Fe7S8). Extracted from "DSC and XRD heating.docx".

## v2 (blank-subtracted, 2026-06-09) — 2026 re-measurement
Re-measured on a fresh bulk Fe3S4 pellet (Setaram Labsys EVO, Ar, ~10 C/min to ~600 C),
blank-subtracted against an empty alumina crucible. Converted from the instrument export
to this repo's UTF-16 SETARAM `.txt` layout.

- `Fe3S4_bulk_heating_v2.txt` — TG-DSC heating. Initial mass 23.9313 mg.
  Columns: Sample Temperature (C), TG %, HeatFlow (mW). TG 100% -> ~84.4%.
- `Fe3S4_bulk_cooling_v2.txt` — DSC cooling, ~20 C/min. DSC only (NO TG channel).
  Columns: Sample Temperature (C), HeatFlow (mW).

Licence: code is MIT (repo-root `LICENSE`); the data and figures in this directory are
licensed CC-BY 4.0 — full text in `LICENSE-CC-BY-4.0.txt`. The v1 (2020) run is the one in
manuscript Fig. S7 (`../make_Figure_S7_dsc.py`); v2 is plotted by
`../make_dsc_2026_remeasurement.py` and is not in the current manuscript.
