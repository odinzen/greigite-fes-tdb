"""Kinetic falsification figure: solid-state S diffusion vs the observed conversion.

Plots the solid-state sulfur-diffusion conversion time of a one-micrometre greigite crystal against
the observed mackinawite->greigite->pyrite transformation times (Hunger & Benning 2007, minutes to
hours at 100-200 C). The conversion time is the characteristic diffusion time L^2/D*(T) with the
measured pyrite sulfur self-diffusion coefficient (Ea = 132 kJ/mol); it was computed with the
Odinzen equilibrium engine (proprietary; only the tabulated result is distributed here, in
kinetics_conversion_time.csv). The many-order-of-magnitude gap, plus the ~2x activation-energy
mismatch against the measured transformation, rules out bulk solid-state diffusion and leaves the
solution-mediated, oxidant-limited pathway.

No title or finding text is baked onto the image (it belongs in the manuscript caption).
Writes artifacts/figures/fig_kinetics_falsification.png.
"""
import csv, pathlib, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = pathlib.Path(__file__).resolve().parents[1]
ART = REPO / "artifacts" / "figures"
ART.mkdir(parents=True, exist_ok=True)
HERE = pathlib.Path(__file__).resolve().parent

Tc, nom, hi, lo = [], [], [], []
with open(HERE / "kinetics_conversion_time.csv") as f:
    for r in csv.DictReader(f):
        Tc.append(float(r["Tc"])); nom.append(float(r["t_nom_yr"]))
        hi.append(float(r["t_hi_yr"])); lo.append(float(r["t_lo_yr"]))
Tc, nom, hi, lo = map(np.array, (Tc, nom, hi, lo))

plt.rcParams.update({"font.size": 12})
fig, ax = plt.subplots(figsize=(8.2, 6.6))

# solid-state diffusion band (Ea +/- 12.5 kJ/mol) + nominal
ax.fill_between(Tc, lo, hi, color="0.75", alpha=0.7, lw=0, zorder=2)
ax.plot(Tc, nom, "k-", lw=2.6, zorder=3, label="solid-state S diffusion (E$_a$=132 kJ mol$^{-1}$)")

# observed transformation (Hunger & Benning 2007): minutes to ~1.3 h over 100-200 C
obs_T = np.array([100, 125, 150, 175, 200])
obs_t = np.array([1.46e-4, 5.7e-5, 2.9e-5, 1.9e-5, 9.5e-6])   # yr (77 min -> minutes)
ax.fill_between([100, 200], [1.7e-4, 1.7e-4], [8e-6, 8e-6], color="0.35", alpha=0.30, lw=0, zorder=1)
ax.plot(obs_T, obs_t, "o", ms=9, mfc="0.15", mec="k", zorder=4,
        label="observed transformation (Hunger & Benning 2007, E$_a$≈66 kJ mol$^{-1}$)")

# reference timescale lines
for yv, txt in [(1 / 8760, "1 hour"), (1, "1 year"), (1e3, "1 kyr"), (1e6, "1 Myr"),
                (1e9, "1 Gyr"), (1.37e10, "age of universe")]:
    ax.axhline(yv, color="0.85", lw=0.8, zorder=0)
    ax.text(249, yv * 1.3, txt, fontsize=7.5, ha="right", color="0.4")

ax.axvspan(0, 25, color="0.92", zorder=0)
ax.text(2, 3e-3, "diagenetic\nsediment", fontsize=8.5, style="italic")

# separation between the two: arrow spanning the gap (magnitude only, no conclusion text)
ax.annotate("", xy=(150, 3.7e10), xytext=(150, 2.9e-5),
            arrowprops=dict(arrowstyle="<->", color="k", lw=1.4), zorder=5)
ax.text(153, 3e3, "≥ 10$^{13}$×", fontsize=10, va="center")

ax.set_yscale("log")
ax.set_xlabel("Temperature (°C)", fontsize=14)
ax.set_ylabel("Greigite → pyrite conversion time (years)", fontsize=14)
ax.set_xlim(0, 250); ax.set_ylim(1e-6, 1e22)
ax.legend(fontsize=9.5, loc="upper right", framealpha=0.95)
fig.tight_layout()
out = ART / "fig_kinetics_falsification.png"
fig.savefig(out, dpi=300)
print("saved", out)
