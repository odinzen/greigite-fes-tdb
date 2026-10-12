"""Fig. 4 — kinetic falsification: solid-state S diffusion vs the observed conversion.

Plots the solid-state sulfur-diffusion conversion time of a one-micrometre greigite crystal against
the observed total reaction times of the mackinawite -> greigite -> pyrite transformation, read
from Hunger & Benning (2007) Table 2 (hunger_benning_2007_table2.csv). The 100 C run was stopped
before completion, so it is drawn as a lower bound. The gap label is computed from the data at the
temperature where the gap is smallest.

The conversion time is the characteristic diffusion time L^2/D*(T) with the measured pyrite
sulfur self-diffusion coefficient (Ea = 132 kJ/mol); it was computed with the Odinzen equilibrium
engine (proprietary; only the tabulated result is distributed here, in
kinetics_conversion_time.csv).

No title or finding text is baked onto the image (it belongs in the manuscript caption).
Writes artifacts/figures/Figure_4.png.
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

# observed total reaction times, Hunger & Benning (2007) Table 2 (p. 8): mackinawite -> greigite -> pyrite
MIN_PER_YR = 365.25 * 24 * 60
obs = list(csv.DictReader(open(HERE / "hunger_benning_2007_table2.csv", encoding="utf-8")))
obs_T = np.array([float(r["T_C"]) for r in obs])
obs_t = np.array([float(r["reaction_time_min"]) for r in obs]) / MIN_PER_YR
done = np.array([r["complete"] == "yes" for r in obs])
ax.fill_between([obs_T.min(), obs_T.max()], [obs_t.max() * 1.25] * 2, [obs_t.min() * 0.8] * 2,
                color="0.35", alpha=0.30, lw=0, zorder=1)
ax.plot(obs_T[done], obs_t[done], "o", ms=9, mfc="0.15", mec="k", zorder=4,
        label="observed reaction time, mackinawite → greigite → pyrite\n"
              "(Hunger & Benning 2007, Table 2; E$_a$ = 67.5 kJ mol$^{-1}$ for mackinawite\n"
              "consumption, 36 to 37 kJ mol$^{-1}$ for pyrite growth)")
ax.plot(obs_T[~done], obs_t[~done], "o", ms=9, mfc="white", mec="k", mew=1.4, zorder=4,
        label="reaction incomplete when stopped (lower bound)")
for T, t in zip(obs_T[~done], obs_t[~done]):
    ax.annotate("", xy=(T, t * 60), xytext=(T, t * 1.6),
                arrowprops=dict(arrowstyle="->", color="k", lw=1.2), zorder=4)

# reference timescale lines
for yv, txt in [(1 / 8760, "1 hour"), (1, "1 year"), (1e3, "1 kyr"), (1e6, "1 Myr"),
                (1e9, "1 Gyr"), (1.37e10, "age of universe")]:
    ax.axhline(yv, color="0.85", lw=0.8, zorder=0)
    ax.text(249, yv * 1.3, txt, fontsize=7.5, ha="right", color="0.4")

ax.axvspan(0, 25, color="0.92", zorder=0)
ax.text(2, 3e-3, "diagenetic\nsediment", fontsize=8.5, style="italic")

# gap label: smallest central-line gap over the completed runs (slowest run at that T)
logt_nom = lambda T: np.interp(T, Tc, np.log10(nom))
gap = 10 ** logt_nom(obs_T[done]) / obs_t[done]
k = int(np.argmin(gap))
Tg, tg = obs_T[done][k], obs_t[done][k]
ax.annotate("", xy=(Tg, 10 ** logt_nom(Tg)), xytext=(Tg, tg),
            arrowprops=dict(arrowstyle="<->", color="k", lw=1.4, shrinkA=6), zorder=5)
mant, expo = f"{gap[k]:.0e}".split("e")
ax.text(Tg + 3, 1e-1, f"{mant} × 10$^{{{int(expo)}}}$×\nat {Tg:.0f} °C", fontsize=10, va="center")

ax.set_yscale("log")
ax.set_xlabel("Temperature (°C)", fontsize=14)
ax.set_ylabel("Greigite → pyrite conversion time (years)", fontsize=14)
ax.set_xlim(0, 250); ax.set_ylim(1e-6, 1e22)
ax.legend(fontsize=9.5, loc="upper right", framealpha=0.95)
fig.tight_layout()
out = ART / "Figure_4.png"
fig.savefig(out, dpi=300)
print("saved", out)
