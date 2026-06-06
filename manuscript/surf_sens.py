from pathlib import Path
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
FIG = ROOT / "artifacts" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

# Engine-computed decomposition margin (Dilner basis, pycalphad) -- the anchor.
MARGIN0 = 21.529  # kJ per mol FeS1.33 (engine consistency_report)
SIGMA = 7.3  # kJ per mol FeS1.33 (Subramani 2020 uncertainty)

# Greigite molar volume (Fe3S4) -> per FeS1.33.
M_Fe3S4 = 3 * 55.845 + 4 * 32.06  # g/mol = 295.775
RHO = 4.079  # g/cm3 (greigite)
Vm_Fe3S4 = M_Fe3S4 / RHO * 1e-6  # m3/mol  = 7.25e-5
Vm_FeS133 = Vm_Fe3S4 / 3.0  # m3/mol per FeS1.33

print(
    f"Vm(Fe3S4) = {Vm_Fe3S4 * 1e6:.2f} cm3/mol ; Vm(FeS1.33) = {Vm_FeS133 * 1e6:.3f} cm3/mol"
)


# Surface enthalpy contribution per mol FeS1.33 for a spherical particle of
# diameter d (m), surface energy gamma (J/m2):  dH_surf = gamma * 6 * Vm / d
def dH_surf_kJ(d_nm, gamma):
    d = d_nm * 1e-9
    return gamma * 6.0 * Vm_FeS133 / d / 1000.0  # kJ/mol


# critical diameter where the surface penalty equals the margin (margin->0)
def d_crit_nm(gamma):
    # gamma*6*Vm/d = MARGIN0*1000  -> d = gamma*6*Vm/(MARGIN0*1000)
    return gamma * 6.0 * Vm_FeS133 / (MARGIN0 * 1000.0) * 1e9


for g in (1.0, 1.5, 2.0):
    print(
        f"gamma={g} J/m2:  dH_surf(10nm)={dH_surf_kJ(10, g):.1f} kJ  d_crit(margin=0)={d_crit_nm(g):.1f} nm"
    )

# ---- figure: margin vs particle size, three surface energies ----
d = np.linspace(3, 80, 400)
fig, ax = plt.subplots(figsize=(8.6, 5.8))
colors = {1.0: "#7aa9d6", 1.5: "#3b7dbf", 2.0: "#1f4e85"}
for g in (1.0, 1.5, 2.0):
    margin = MARGIN0 - dH_surf_kJ(d, g)
    ax.plot(
        d,
        margin,
        lw=2,
        color=colors[g],
        label=f"γ = {g:.1f} J/m²  (d* = {d_crit_nm(g):.0f} nm)",
    )
ax.axhline(MARGIN0, ls=":", color="#16431f", lw=1.4)
ax.text(
    70,
    MARGIN0 + 0.4,
    "bulk greigite  (+2.95σ, engine)",
    color="#16431f",
    fontsize=9,
    ha="right",
)
ax.axhline(0, color="0.3", lw=1.2)
ax.text(
    70,
    0.5,
    "metastability threshold (pyrrhotite returns)",
    color="0.3",
    fontsize=9,
    ha="right",
)
ax.fill_between(d, -8, 0, color="#f3d6d6", alpha=0.5)
# secondary sigma axis
ax2 = ax.twinx()
ax2.set_ylim(np.array(ax.get_ylim()) / SIGMA)
ax2.set_ylabel("margin in σ (σ = 7.3 kJ)", fontsize=10)
ax.set_xlabel("greigite particle diameter, nm", fontsize=11)
ax.set_ylabel("greigite decomposition margin, kJ/mol-FeS$_{1.33}$", fontsize=11)
ax.set_title(
    "Sensitivity of greigite stability to nanoparticle surface enthalpy\n"
    "(anchored on the engine margin +21.5 kJ; surface term = γ·6V$_m$/d)",
    fontsize=10.5,
)
ax.set_xlim(3, 80)
ax.set_ylim(-6, 26)
ax.legend(loc="lower right", fontsize=9)
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(str(FIG / "fig_surface_sensitivity.png"), dpi=200)
print("saved fig_surface_sensitivity.png")
