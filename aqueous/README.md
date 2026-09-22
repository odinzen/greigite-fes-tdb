# aqueous — porewater Eh-pH stability of greigite (PHREEQC)

The wet-sediment companion to the anhydrous CALPHAD diagrams: greigite's equilibrium field in real
porewater, on a redox (Eh) vs pH predominance diagram, computed with PHREEQC.

The greigite dissolution constant is built from the same measured thermochemistry as the CALPHAD
databases (Subramani 2020 enthalpy, Shumway 2022 entropy) on a database-consistent aqueous
reference, then PHREEQC returns the predominant Fe mineral over a pH x pe grid at diagenetic
dissolved-S, -Fe and -C levels. The grid is run for the nominal greigite log_K and its ±1σ bounds
so the field carries the calorimetric uncertainty.

## Run

```bash
pip install phreeqpython numpy scipy matplotlib

python aqueous/derive_greigite_logk.py     # log_K derivation + pyrite verifier (stdlib only)
python aqueous/build_ehph_diagram.py       # -> artifacts/aqueous/ehph_fields.npz
python aqueous/make_ehph_figure.py         # -> artifacts/figures/Figure_3.png
```

`build_ehph_diagram.py` uses phreeqpython's bundled `phreeqc.dat` and grafts on the Fe-S(-O) phases
it lacks (greigite, pyrrhotite, magnetite). Mackinawite is excluded as a metastable precursor.

## Files

- `derive_greigite_logk.py` — self-contained derivation of the greigite dissolution log_K
  (-68.95, ±1σ -65.12/-72.79) from cited constants, with the pyrite verifier. Details in
  [`LOGK_DERIVATION.md`](LOGK_DERIVATION.md).
- `build_ehph_diagram.py` — the PHREEQC predominance calculation over the pH x pe grid.
- `make_ehph_figure.py` — greyscale render (greigite field with its ±1σ range, water-stability
  window masked).

The greigite field sits in reducing, near-neutral to alkaline porewater and is robust across the
calorimetric uncertainty.
