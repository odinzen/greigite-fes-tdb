# kinetics — solid-state diffusion falsification

Why a thermodynamically stable mineral is still observed to become pyrite: not by solid-state
diffusion. The characteristic time to convert a one-micrometre greigite crystal by bulk sulfur
diffusion is `L^2 / D*(T)`, with `D*(T) = D0 exp(-Ea/RT)` and the measured pyrite sulfur
self-diffusion coefficient (Ea = 132 ± 12.5 kJ/mol; Watson, Cherniak & Frank 2009). It is of order
10^17 years at 25 °C — many orders of magnitude slower than the observed transformation, which
completes in minutes to hours at 100-200 °C (Hunger & Benning 2007). The operative pathway is
therefore solution-mediated and gated by porewater oxidant supply, not intracrystalline diffusion.

## Provenance and scope

The conversion times were computed with the **Odinzen equilibrium engine**, which is proprietary
software and is **not** part of this repository. Only the tabulated result is distributed here
(`kinetics_conversion_time.csv`), which is sufficient to regenerate the figure. The engine is
available under license or collaboration by request to the corresponding author.

## Files

- `kinetics_conversion_time.csv` — conversion time (years) vs temperature (°C), nominal and the
  ±1σ diffusion-coefficient band.
- `make_kinetics_figure.py` — renders `artifacts/figures/fig_kinetics_falsification.png` from the
  CSV against the observed transformation points.

```bash
pip install numpy matplotlib
python kinetics/make_kinetics_figure.py
```
