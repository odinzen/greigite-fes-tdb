# Manuscript companion

Companion code for the manuscript *"Greigite is thermodynamically stable and its pyritization
is gated by porewater oxidant, not solid-state diffusion"* (Michael E. Bustamante, Gabriel
Bustamante, Tamilarasan Subramani, Kristina Lilova and Alexandra Navrotsky), submitted for
publication.

Everything needed to reproduce the manuscript figures and boundary-fugacity tables is described
in the repository [README](../../README.md): the CALPHAD databases and anhydrous diagrams
(`engine/`, `manuscript/`), the PHREEQC porewater diagram (`aqueous/`) and the solid-state
diffusion falsification (`kinetics/`).

## Key results

- Greigite is thermodynamically stable: its reaction to pyrite and pyrrhotite is endothermic by
  +19.5 ± 7.9 kJ/mol per FeS₁.₃₃ (+2.5σ).
- Suppressing greigite recovers the classical Fe–S and Fe–S–O diagrams; at 300 K the pyrite
  boundary sits at log f(S₂) = −31.8 and native-S saturation at −13.82. With greigite, the pyrite
  onset moves to log f(S₂) = −20.3.
- The greigite dissolution constant on the CODATA aqueous reference is log K = −68.95
  (±1σ: −65.12 / −72.79); the same reference set reproduces the database pyrite constant to
  within 0.04 log units.
- Converting a one-micrometre greigite crystal by bulk sulfur diffusion takes of order 10¹⁷ years
  at 25 °C, so the observed conversion to pyrite is solution-mediated.

## License

Code MIT, data and figures CC BY 4.0; see the repository [`LICENSE`](../../LICENSE).

## Acknowledgements

Built on the open-source pycalphad, ESPEI and PHREEQC projects. Calorimetric data: Navrotsky
group, Arizona State University.
