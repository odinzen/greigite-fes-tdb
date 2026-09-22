# Greigite aqueous dissolution log_K derivation

Produced by [`derive_greigite_logk.py`](derive_greigite_logk.py) (standard library only). Purpose:
give the porewater Eh-pH diagram a greigite dissolution constant built from the same measured
thermochemistry as the CALPHAD databases — Subramani et al. 2020 formation enthalpy and Shumway et
al. 2022 entropy — on a database-consistent aqueous reference, replacing the legacy value carried in
older aqueous databases.

## Result (298.15 K)

| Phase | Reaction | log_K | source |
|---|---|---|---|
| Greigite (nominal) | Fe3S4 + 4H+ = 2Fe+3 + Fe+2 + 4HS- | **-68.95** | this work |
| Greigite (+1σ, least stable) | " | -65.12 | this work |
| Greigite (-1σ, most stable) | " | -72.79 | this work |
| Greigite (legacy) | " | -45.035 | older aqueous database (replaced) |
| Pyrrhotite/FeS (added) | FeS + H+ = Fe+2 + HS- | -6.03 | this work |

greigite dGf(CODATA, 298 K) = -433.48 kJ/mol Fe3S4 (nominal); -411.58 (+1σ); -455.38 (-1σ).

## Why it is trustworthy

1. **Verifier pass.** The same aqueous set reproduces the reference pyrite dissolution constant
   (FeS2 + 2H+ + 2e- = Fe+2 + 2HS-) at log_K = -18.435 against the database value of -18.479, a
   residual of 0.04 log units. The aqueous species set (Robie-Hemingway 1995 / CODATA: Fe+2 =
   -78.90, Fe+3 = -4.60, HS- = +12.05 kJ/mol) is therefore database-consistent.
2. **Entropy basis confirmed.** S(298 K) = 215.11 J/mol/K = 3 × Shumway's per-FeS1.33 value
   (71.704), so the description is per-Fe3S4 (no factor-of-three error).
3. **Reference reconciliation.** A CALPHAD solid Gibbs energy sits on the SGTE element reference
   and is not directly mixable with CODATA aqueous species. The aqueous log_K uses
   dGf = dHf − T(S_greigite − 3 S_Fe − 4 S_S) with CODATA element entropies (S_Fe = 27.28,
   S_S = 32.05), giving -433.48 kJ/mol. Using a raw SGTE-referenced formation term would be
   several log units wrong.

## Robustness

Across the full ±1σ calorimetric enthalpy band the aqueous log_K stays between -65.1 and -72.8, at
least 20 log units more stable than the legacy value. The ±1σ uncertainty that flips the pyrrhotite
field in the anhydrous f(S2)-f(O2) diagram does not flip the porewater conclusion; the greigite
field is robust to the calorimetric uncertainty. The remaining exposure is a systematic offset in
the single calorimetric enthalpy.

## Inputs (all cited constants)

- dHf(greigite) = -431.97 ± 21.9 kJ/mol Fe3S4, CODATA reference — Subramani et al. 2020.
- S(greigite, 298 K) = 215.11 J/mol/K — Shumway et al. 2022 (3 × per-FeS1.33).
- CODATA element entropies S_Fe = 27.28, S_S = 32.05 J/mol/K.
- Aqueous species dGf (kJ/mol): Fe+2 = -78.90, Fe+3 = -4.60, HS- = +12.05 — Robie-Hemingway 1995.
- Pyrite anchor: dHf = -171.5 kJ/mol, S = 52.9 J/mol/K — Robie-Hemingway 1995.
