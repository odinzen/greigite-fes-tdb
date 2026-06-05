"""Heat capacity, entropy and enthalpy increment of greigite (Fe3S4).

All quantities in this module are expressed PER MOLE OF FeS1.33, i.e. per
one-third of a formula unit of Fe3S4 (Fe3S4 == 3 x FeS1.33). Working in the
FeS1.33 sub-formula keeps the numbers on the same per-metal footing as the
mono-/di-sulfides in thermo_data.py and matches how Shumway 2022 report S298.

Cp(T) above the low-T calorimetry is the Debye-Einstein-plus-linear fit of
Shumway, Wilson, Lilova, Subramani, Navrotsky & Woodfield,
J. Chem. Thermodyn. 173 (2022) 106836 (their Table 4). S298 and Cp298 are the
calorimetric values from the same paper.

Greigite is metastable: it decomposes (to pyrrhotite + pyrite) on heating
near GREIGITE_DECOMP_T, so the fit is only physically meaningful below that.

Units: Cp in J/mol/K (per FeS1.33) ; S in J/mol/K ; H increment in J/mol ;
T in K.
"""

import numpy as np

# --- Constants -------------------------------------------------------------
R = 8.314462618  # J/mol/K, CODATA 2018
T_REF = 298.15  # K, standard reference temperature

# Greigite is metastable; it decomposes around this temperature on heating.
GREIGITE_DECOMP_T = 530.0  # K (SHU2022)

# Calorimetric standard-state values, per mole FeS1.33 (SHU2022).
S_298 = 71.334  # J/mol/K, third-law entropy at 298.15 K
CP_298 = 59.686  # J/mol/K, measured Cp at 298.15 K

# Shumway 2022 Table 4: Debye-Einstein + linear high-T fit coefficients.
# Cp = m*Debye(THETA_D) + n*Einstein(THETA_E) + A_LIN*T  (per FeS1.33).
m = 0.9740  # number of Debye oscillators
THETA_D = 237.85  # K, Debye temperature
n = 1.1329  # number of Einstein oscillators
THETA_E = 405.26  # K, Einstein temperature
A_LIN = 0.0398  # J/mol/K^2, linear (electronic/dilation) term


def _debye_cp(theta_over_T):
    """Per-mole-oscillator Debye heat capacity (J/mol/K) -> 3R at high T.

        C = 9R (T/theta)^3 * integral_0^{theta/T} u^4 e^u / (e^u - 1)^2 du

    Evaluated by numerical (trapezoid) integration on a fine grid. For very
    small theta/T (i.e. T >> theta) the classical Dulong-Petit limit 3R is
    returned directly to avoid 0/0 in the integrand.
    """
    x = theta_over_T
    if x < 1e-8:
        return 3.0 * R
    u = np.linspace(0.0, x, 4000)
    # integrand u^4 e^u / (e^u - 1)^2 ; the u=0 endpoint is a 0/0 limit -> 0.
    integrand = np.zeros_like(u)
    nz = u > 1e-12
    eu = np.exp(u[nz])
    integrand[nz] = u[nz] ** 4 * eu / (eu - 1.0) ** 2
    integral = np.trapezoid(integrand, u)
    return 9.0 * R * (1.0 / x) ** 3 * integral


def _einstein_cp(x):
    """Per-mole-oscillator Einstein heat capacity (J/mol/K).

    C = 3R x^2 e^x / (e^x - 1)^2 ,  with x = THETA_E / T.
    """
    ex = np.exp(x)
    return 3.0 * R * x**2 * ex / (ex - 1.0) ** 2


def cp_debye_einstein(T):
    """Greigite Cp(T) per FeS1.33 from the SHU2022 three-term fit (J/mol/K).

    Accepts a scalar or a numpy array of temperatures.
    """
    T = np.asarray(T, dtype=float)
    scalar = T.ndim == 0
    Tf = np.atleast_1d(T)
    out = np.empty_like(Tf)
    for i, Ti in enumerate(Tf):
        debye = m * _debye_cp(THETA_D / Ti)
        einst = n * _einstein_cp(THETA_E / Ti)
        out[i] = debye + einst + A_LIN * Ti
    return out[0] if scalar else out


def entropy(T):
    """Third-law entropy S(T) per FeS1.33 (J/mol/K).

        S(T) = S_298 + integral_{298.15}^{T} Cp/T dT

    Integrated by the trapezoid rule on a fine grid from 298.15 K to T.
    """
    grid = np.linspace(T_REF, float(T), 2000)
    cp = cp_debye_einstein(grid)
    return S_298 + np.trapezoid(cp / grid, grid)


def enthalpy_increment(T):
    """Enthalpy increment H(T) - H(298.15) per FeS1.33 (J/mol).

        H(T) - H(298) = integral_{298.15}^{T} Cp dT

    Integrated by the trapezoid rule on a fine grid from 298.15 K to T.
    """
    grid = np.linspace(T_REF, float(T), 2000)
    cp = cp_debye_einstein(grid)
    return np.trapezoid(cp, grid)


if __name__ == "__main__":
    # Validation against SHU2022 published values.
    cp298 = cp_debye_einstein(298.15)
    s300 = entropy(300.0)
    print(f"Cp(298.15 K) = {cp298:.3f} J/mol/K   (expect ~59.69)")
    print(f"S(300.0 K)   = {s300:.3f} J/mol/K   (Shumway Table 5: 71.704)")
    assert abs(cp298 - 59.69) < 0.1, cp298
    assert abs(s300 - 71.704) < 0.1, s300
    print("greigite_cp self-checks passed.")
