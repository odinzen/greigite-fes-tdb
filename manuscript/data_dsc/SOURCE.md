README
# Bulk Fe3S4 (greigite) decomposition data — source provenance

Experimental data measured by K. Lilova (co-author), Navrotsky group, ASU.

- `Fe3S4_bulk_heating.txt`  — DSC/TG, bulk greigite, heating to 600 C (SETARAM export, UTF-16).
- `Fe3S4_bulk_cooling.txt`  — DSC/TG, bulk greigite, cooling run.
  Columns: Index, Time (s), Sample Temperature (C), TG (mg, baseline-corrected), HeatFlow (uV).
- `Fe3S4_bulk_XRD_postDSC600.png` — powder XRD of bulk greigite after DSC to 600 C;
  product is pyrrhotite-3T (Fe7S8). Extracted from "DSC and XRD heating.docx".

Licence: CC-BY 4.0 (figure/data). Used to build Fig. 8 via make_fig8_dscxrd_bw.py.
