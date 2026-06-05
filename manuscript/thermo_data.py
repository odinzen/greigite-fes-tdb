"""Sourced standard thermodynamic data for the Fe-S(-O) system.

EVERY value below is transcribed from a named primary source. Nothing is
estimated or fabricated. Where the literature is internally inconsistent
(e.g. two NIST values for pyrite formation enthalpy) the alternatives are
documented in comments and the geochemically-standard choice is used.

Reference key
-------------
JANAF  : NIST-JANAF Thermochemical Tables, 4th ed. (Chase, M.W. Jr.,
         J. Phys. Chem. Ref. Data, Monograph 9, 1998). Table numbers cited.
GW1962 : Gronvold & Westrum, Inorg. Chem. 1 (1962) 36-48 (pyrite Cp, S).
GS1992 : Gronvold & Stolen, J. Chem. Thermodyn. (1992) (pyrrhotite Fe0.877S).
SUB2020: Subramani, Lilova, Abramchuk, Leinenweber, Navrotsky,
         PNAS 117 (2020) 28645  (greigite enthalpy of formation).
SHU2022: Shumway, Wilson, Lilova, Subramani, Navrotsky, Woodfield,
         J. Chem. Thermodyn. 173 (2022) 106836 (greigite low-T Cp, S298).

Units: dHf298 in J/mol ; S298 in J/mol/K ; Cp in J/mol/K ; T in K.
All sulfide/oxide values are PER FORMULA UNIT AS WRITTEN in the key below.
Greigite is handled in greigite_cp.py (measured + extrapolated Cp).
"""

import numpy as np

# --- Universal constants ---------------------------------------------------
R = 8.314462618  # J/mol/K, CODATA 2018
T_REF = 298.15  # K, standard reference temperature


def shomate_cp(T, A, B, C, D, E):
    """NIST Shomate equation for Cp (J/mol/K).

    Cp = A + B t + C t^2 + D t^3 + E / t^2,  with t = T/1000 (T in K).
    Coefficients A..E are the JANAF Shomate parameters for each phase.
    """
    t = T / 1000.0
    return A + B * t + C * t**2 + D * t**3 + E / t**2


# --- Gas-phase species -----------------------------------------------------
# Diatomic sulfur gas. dHf and S298 from JANAF S-012 (Chase 1998).
S2_g = {
    "dHf298": 128600.0,  # J/mol  (JANAF S-012 / Chase 1998)
    "S298": 228.165,  # J/mol/K
    "shomate": (33.51313, 5.065360, -1.059670, 0.089905, -0.211911),
    "source": "JANAF S-012 / Chase 1998",
}

# Dioxygen gas. JANAF O-029 (Chase 1998). Reference element -> dHf = 0.
O2_g = {
    "dHf298": 0.0,
    "S298": 205.147,
    "shomate": (31.32234, -20.23531, 57.86644, -36.50624, -0.007374),
    "source": "JANAF O-029 / Chase 1998",
}

# --- Condensed phases ------------------------------------------------------
# Body-centred-cubic alpha-iron. Reference element -> dHf = 0.
Fe_bcc = {
    "dHf298": 0.0,
    "S298": 27.28,  # JANAF Fe / Chase 1998
    "shomate": (18.42868, 24.64301, -8.913720, 9.664706, -0.012643),
    "source": "JANAF Fe / Chase 1998",
}

# Pyrite FeS2. dHf298 = -171544 J/mol (JANAF Fe-028). S298 from GW1962.
# (An alternative JANAF value of -178200 J/mol exists; the GW1962-consistent
#  -171544 is the geochemically-standard choice and used here.)
FeS2_pyrite = {
    "dHf298": -171544.0,
    "S298": 52.916,  # GW1962
    "shomate": (82.30388, -23.72240, 35.14598, -9.598515, -1.396524),
    "source": "JANAF Fe-028 ; S298 from GW1962",
}

# Pyrrhotite Fe0.877S (the monoclinic 4C superstructure composition).
# dHf298, S298 and the measured Cp table are from GS1992 / JANAF Fe-002.
FeS_pyrrhotite_0877 = {
    "formula_FexS": 0.877,
    "dHf298": -105437.0,
    "S298": 60.799,
    "cp_table_T": [298.15, 300.0, 400.0, 500.0, 598.0, 600.0, 700.0],
    "cp_table_Cp": [49.883, 50.000, 56.149, 62.383, 58.175, 58.116, 55.564],
    "source": "GS1992 / JANAF Fe-002 (Fe0.877S)",
}

# Stoichiometric troilite FeS. JANAF Fe-002.
FeS_troilite = {
    "dHf298": -101671.0,
    "S298": 60.321,
    "source": "JANAF Fe-002 (FeS troilite)",
}

# --- Iron oxides -----------------------------------------------------------
# Magnetite Fe3O4. JANAF Fe-023.
Fe3O4_magnetite = {
    "dHf298": -1120894.0,
    "S298": 145.2,
    "shomate": (104.2096, 178.5108, 10.61510, 1.132534, -0.994202),
    "source": "JANAF Fe-023",
}

# Hematite Fe2O3. JANAF Fe-022.
Fe2O3_hematite = {
    "dHf298": -825503.0,
    "S298": 87.28,
    "shomate": (93.43834, 108.3577, -50.86447, 25.58683, -1.611330),
    "source": "JANAF Fe-022",
}

# Wustite FeO. JANAF Fe-021.
FeO_wustite = {
    "dHf298": -272044.0,
    "S298": 60.75,
    "shomate": (45.75120, 18.78553, -5.952201, 0.852779, -0.081265),
    "source": "JANAF Fe-021",
}


# --- Accessors -------------------------------------------------------------
def cp_shomate(phase, T):
    """Cp(T) for any phase carrying a 'shomate' tuple (J/mol/K)."""
    return shomate_cp(T, *phase["shomate"])


def cp_pyrrhotite(T):
    """Cp(T) for Fe0.877S by linear interpolation of the measured table.

    GS1992 measured Cp on a coarse grid; geochemical use only needs values
    between 298 and 700 K, so linear interpolation is adequate. np.interp
    clamps at the table endpoints (no extrapolation).
    """
    return np.interp(
        T,
        FeS_pyrrhotite_0877["cp_table_T"],
        FeS_pyrrhotite_0877["cp_table_Cp"],
    )
