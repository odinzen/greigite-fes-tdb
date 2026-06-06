import sys, numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
MAN = HERE.parent
ROOT = MAN.parent
TDB = ROOT / "artifacts" / "tdb"
FIG = ROOT / "artifacts" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(MAN))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import thermo_data as td, greigite_cp as g, reactions as rx

# Decomposition per FeS1.33: FeS1.33 -> 0.726 Po(Fe0.877S) + 0.274 FeS2  (Shumway Eq.6 basis)
A_PO, A_PY = 0.726, 0.274
# Gr in rx.PHASE is per Fe3S4; use -1/3 to get per-FeS1.33
STO = {"Gr": -1 / 3.0, "Po": A_PO, "FeS2": A_PY}
T = np.linspace(298.15, 600.0, 80)

dG = rx.reaction_dG(STO, T) / 1000.0  # kJ per mol FeS1.33
# entropy of reaction at 298 from S298 values
S_gr = g.S_298
S_po = td.FeS_pyrrhotite_0877["S298"]
S_py = td.FeS2_pyrite["S298"]
dS298 = A_PO * S_po + A_PY * S_py - S_gr
dH298 = (
    A_PO * td.FeS_pyrrhotite_0877["dHf298"]
    + A_PY * td.FeS2_pyrite["dHf298"]
    - 3 * (-144100) / 3.0
) / 1000.0  # kJ (greigite per FeS1.33 = -144.1)
print(
    f"decomposition margin dH298 = {dH298:+.2f} kJ/mol-FeS1.33  (Subramani measured +18.9)"
)
print(
    f"dS_rxn(298)               = {dS298:+.2f} J/mol/K   -> T*dS at 600K = {600 * dS298 / 1000:+.2f} kJ"
)
print(f"dG(298)                   = {dG[0]:+.2f} kJ ; dG(600) = {dG[-1]:+.2f} kJ")

# --- Uncertainty budget (per FeS1.33) ---
sig_gr = 7.3  # Subramani greigite dHf
sig_po = 3.5  # Xu&Navrotsky 2010 pyrrhotite (~Fe0.875S)
sig_py = 0.9  # Gronvold-Westrum pyrite (direct-reaction calorimetry)
# enthalpy contributions to the decomposition margin (coeff * sigma)
c_gr, c_po, c_py = 1.0 * sig_gr, A_PO * sig_po, A_PY * sig_py
sig_indep = (c_gr**2 + c_po**2 + c_py**2) ** 0.5  # independent quadrature
sig_worst = c_gr + c_po + c_py  # correlated worst-case (linear)
sig_shumway = 8.6  # Shumway 2022 published, fully propagated
# entropy uncertainty: ~1% of phase entropies (Shumway). Contribution to margin uncertainty:
sigS = 0.01 * np.array([S_gr, S_po, S_py])
sigTS_600 = (
    600
    * ((A_PO * sigS[1]) ** 2 + (A_PY * sigS[2]) ** 2 + (1 * sigS[0]) ** 2) ** 0.5
    / 1000.0
)
print(
    f"\nenthalpy sigma: greigite {c_gr:.2f} | pyrrhotite {c_po:.2f} | pyrite {c_py:.2f} kJ"
)
print(
    f"  independent (quadrature) = {sig_indep:.2f} kJ ; worst-case (linear) = {sig_worst:.2f} kJ ; Shumway published = {sig_shumway} kJ"
)
print(
    f"entropy-uncertainty contribution at 600 K (T*sigma_dS) = {sigTS_600:.2f} kJ  -> negligible vs enthalpy"
)
print(
    f"margin / sigma:  {dH298 / sig_shumway:+.2f} sigma (Shumway u) ; {dG[0] / sig_indep:+.2f} sigma (indep)"
)

# --- Cp-extrapolation sensitivity: perturb greigite linear Cp term A_LIN by +-30% above 530 K ---
base = g.A_LIN


def margin_with_Alin(factor):
    g.A_LIN = base * factor
    val = rx.reaction_dG(STO, np.array([600.0]))[0] / 1000.0
    g.A_LIN = base
    return val


dG600_base = dG[-1]
dG600_hi = margin_with_Alin(1.3)
dG600_lo = margin_with_Alin(0.7)
print(
    f"\nCp-extrapolation test (greigite A_lin +-30%): dG(600) = {dG600_base:.2f} [{dG600_lo:.2f}, {dG600_hi:.2f}] kJ  -> spread {abs(dG600_hi - dG600_lo):.2f} kJ"
)

# --- Figure: margin vs T with both envelopes ---
fig, ax = plt.subplots(figsize=(8.4, 5.2))
ax.axhline(0, color="k", lw=0.8)
ax.fill_between(
    T,
    dG - sig_worst,
    dG + sig_worst,
    color="#e0a0a0",
    alpha=0.35,
    label="worst-case envelope (±%.1f kJ, correlated)" % sig_worst,
)
ax.fill_between(
    T,
    dG - sig_shumway,
    dG + sig_shumway,
    color="#6a9bd1",
    alpha=0.45,
    label="independent / published envelope (±8.6 kJ, Shumway 2022)",
)
ax.plot(T, dG, color="#1f3864", lw=2.2, label="decomposition margin ΔG(T)")
ax.axvspan(590, 600, color="grey", alpha=0.25, hatch="////", lw=0)
ax.axvline(590, color="k", ls="--", lw=1)
ax.text(
    588,
    ax.get_ylim()[0],
    "pyritization onset 590 K",
    rotation=90,
    va="bottom",
    ha="right",
    fontsize=8,
)
ax.set_xlabel("Temperature / K")
ax.set_ylabel("Decomposition margin ΔG  /  kJ·mol⁻¹ (per FeS$_{1.33}$)")
ax.set_title(
    "Greigite decomposition margin with full error propagation\nFeS$_{1.33}$ → 0.726 Fe$_{1-x}$S + 0.274 FeS$_2$  (margin>0 ⇒ greigite stable)"
)
ax.set_xlim(298, 600)
ax.legend(loc="upper left", fontsize=8.5)
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(str(FIG / "fig_propagation_margin.png"), dpi=150)
print("\nSAVED fig_propagation_margin.png")
