# greigite-feso — an open CALPHAD workflow for Fe–S(–O) predominance

Reproducible pipeline that grafts measured greigite (Fe₃S₄) thermodynamics onto an
Fe–S(–O) CALPHAD database (pycalphad) and computes Fe–S log *f*(S₂)–T and
Fe–S–O log *f*(O₂)–log *f*(S₂) predominance diagrams by Gibbs-energy minimisation.

Companion code for the manuscript *"An open, reproducible CALPHAD workflow for Fe–S(–O)
predominance: greigite (Fe₃S₄) is a stable Fe(II,III) thiospinel and the sulfur analogue
of magnetite"* (Bustamante, Bustamante & Lilova), submitted for publication.

## What it does
- Builds a greigite line-compound Gibbs energy from measured ΔfH° (Subramani et al. 2020),
  S° and Debye–Einstein Cp (Shumway et al. 2022) — deterministic, byte-reproducible.
- Grafts greigite (and pyrite) onto the Dilner–Mao–Selleby (2015) Fe–S and Dilner–Selleby
  (2017) Ca–Fe–O–S assessments on a common SGTE91 reference.
- Computes predominance diagrams (grand-potential argmin, after Holland 1959) with
  all-compound ±1σ uncertainty envelopes and native-S saturation capping.
- Validation: suppressing greigite recovers the classical pyrrhotite–pyrite assessment.

## Requirements
Python 3.10+, `pycalphad==0.11.1`, NumPy, SciPy, Matplotlib. Install:
```
python3 -m venv .venv && .venv/bin/pip install pycalphad==0.11.1 matplotlib
```

## Reproduce the figures
```
python scripts/greigite/build_greigite_tdb.py          # build greigite TDB variants
python scripts/greigite/build_boundary_tdbs.py         # all-compound boundary-case TDBs
python scripts/greigite/validate_fes_engine.py         # central Fe–S engine diagram
python scripts/fe_s_greigite/make_fig3_engine.py       # Fe–S–O single-database predominance
# (Fig. 1, 2, 4, 5, 6 scripts likewise; see manuscript figure list)
```

## Key results
- Greigite is thermodynamically stable; decomposition margin **+19.5 ± 7.9 kJ·mol⁻¹ (+2.5σ)**.
- Native-S saturation at 298.15 K: log *f*(S₂) = −13.96 (engine).
- Oxide formation energies reproduce literature to ~1%.

## License
**MIT License** (see `LICENSE`). OSI-approved, permissive.

## Citation
Bustamante, M., Bustamante, G., Lilova, K. (2026). An open, reproducible CALPHAD workflow
for Fe–S(–O) predominance. Manuscript submitted for publication. Archived: Zenodo
(extends DOI 10.5281/zenodo.19835550).

## Acknowledgements
Built on the open-source pycalphad and ESPEI projects. Calorimetric data: Navrotsky group,
Arizona State University.
