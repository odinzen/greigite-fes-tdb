# kinetics — solid-state diffusion falsification

Why a thermodynamically stable mineral is still observed to become pyrite: not by solid-state
diffusion. The characteristic time to convert a one-micrometre greigite crystal by bulk sulfur
diffusion is `L^2 / D*(T)`, with `D*(T) = D0 exp(-Ea/RT)` and the measured pyrite sulfur
self-diffusion coefficient (Ea = 132 ± 12.5 kJ/mol; Watson, Cherniak & Frank 2009). It is of order
10^17 years at 25 °C, many orders of magnitude slower than the observed transformation, which
completes in 60 to 110 minutes at 125-200 °C (Hunger & Benning 2007, Table 2; the 100 °C run
was stopped after 311 minutes, before completion). Bulk sulfur diffusion
therefore cannot control the conversion; the faster pathways need an electron acceptor supplied by
the porewater (zero-valent sulfur in polysulfides, or dissolved H2S).

## Provenance and scope

The conversion times were computed with the **Odinzen equilibrium engine**, which is proprietary
software and is **not** part of this repository. Only the tabulated result is distributed here
(`kinetics_conversion_time.csv`), which is sufficient to regenerate the figure. The engine is
available under license or collaboration by request to the corresponding author.

## Files

- `kinetics_conversion_time.csv` holds the conversion time (years) against temperature (°C),
  nominal and the band from the activation-energy uncertainty (132 ± 12.5 kJ/mol).
- `hunger_benning_2007_table2.csv` holds the observed total reaction times transcribed from
  Hunger & Benning (2007) Table 2 (p. 8), with completion status; the 100 °C time is a lower bound.
- `make_kinetics_figure.py` renders `artifacts/figures/Figure_4.png` from the two tables.

```bash
pip install numpy matplotlib
python kinetics/make_kinetics_figure.py
```
