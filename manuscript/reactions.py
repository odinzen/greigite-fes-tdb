"""Reaction engine for Fe-S(-O) phase-boundary sulfur fugacities.

Computes dG(T) for solid<->solid(+S2) reactions and converts them into
log10 f(S2) phase boundaries, using:

  * standard-state data (dHf298, S298) and Cp(T) functions from thermo_data,
  * the greigite Debye-Einstein Cp fit from greigite_cp.

For every reaction:

    dG(T) = dH298 + INT_298^T dCp dT  -  T * ( dS298 + INT_298^T dCp/T dT )

where dH298, dS298 are the stoichiometric sums at 298.15 K and dCp(T) is the
stoichiometric sum of phase heat capacities. The Cp integrals are evaluated
by cumulative trapezoid integration on a shared temperature grid, so the
whole function is vectorised over an array of target temperatures T.

Sign convention: stoichiometric coefficients nu are NEGATIVE for reactants
and POSITIVE for products. log f(S2) boundaries assume S2(g) is the only
fugacity-bearing species, with stoichiometric coefficient nu_S2 (reactant,
so nu_S2 here is the POSITIVE magnitude consumed per reaction).

Greigite enthalpy carries a published +-7300 J/mol-FeS1.33 uncertainty
(Subramani 2020); set_greigite_dHf_offset / greigite_boundaries expose it so
admissible-band envelopes can be drawn.

Units: dG in J/mol-reaction ; T in K ; log f(S2) dimensionless (log10).
"""

import numpy as np

import thermo_data as td
import greigite_cp as g

# --- Constants -------------------------------------------------------------
R = 8.314462618  # J/mol/K, CODATA 2018
LN10 = np.log(10.0)  # natural log of 10, for dG -> log10 f(S2)
T_REF = 298.15  # K, standard reference temperature


# --- Per-phase Cp(T) functions (all J/mol/K, per formula unit) -------------
def cp_Fe(T):
    return td.cp_shomate(td.Fe_bcc, T)


def cp_FeS2(T):
    return td.cp_shomate(td.FeS2_pyrite, T)


def cp_Po(T):
    # Pyrrhotite Fe0.877S via interpolation of the measured GS1992 table.
    return td.cp_pyrrhotite(T)


def cp_S2(T):
    return td.cp_shomate(td.S2_g, T)


def cp_O2(T):
    return td.cp_shomate(td.O2_g, T)


def cp_Fe3O4(T):
    return td.cp_shomate(td.Fe3O4_magnetite, T)


def cp_Fe2O3(T):
    return td.cp_shomate(td.Fe2O3_hematite, T)


def cp_FeO(T):
    return td.cp_shomate(td.FeO_wustite, T)


def cp_Gr(T):
    # Greigite Cp PER Fe3S4 = 3 x Cp per FeS1.33 (greigite_cp works in FeS1.33).
    return 3.0 * g.cp_debye_einstein(T)


# --- Phase registry: name -> (Cp_func, dHf_per_formula, S_per_formula) ------
# Greigite Fe3S4: dHf298 = 3 * (-144100) J/mol-FeS1.33 (Subramani 2020),
# S298 = 3 * g.S_298 (Shumway 2022). Per the Fe3S4 formula unit.
PHASE = {
    "Fe": (cp_Fe, td.Fe_bcc["dHf298"], td.Fe_bcc["S298"]),
    "FeS2": (cp_FeS2, td.FeS2_pyrite["dHf298"], td.FeS2_pyrite["S298"]),
    "Po": (cp_Po, td.FeS_pyrrhotite_0877["dHf298"], td.FeS_pyrrhotite_0877["S298"]),
    "S2": (cp_S2, td.S2_g["dHf298"], td.S2_g["S298"]),
    "O2": (cp_O2, td.O2_g["dHf298"], td.O2_g["S298"]),
    "Fe3O4": (cp_Fe3O4, td.Fe3O4_magnetite["dHf298"], td.Fe3O4_magnetite["S298"]),
    "Fe2O3": (cp_Fe2O3, td.Fe2O3_hematite["dHf298"], td.Fe2O3_hematite["S298"]),
    "FeO": (cp_FeO, td.FeO_wustite["dHf298"], td.FeO_wustite["S298"]),
    "Gr": (cp_Gr, 3 * (-144100), 3 * g.S_298),
}


def reaction_dG(stoich, T):
    """dG(T) for a reaction (J/mol-reaction), vectorised over T.

    Parameters
    ----------
    stoich : dict {phase_name: nu}
        nu < 0 for reactants, nu > 0 for products.
    T : float or array of K
        Target temperature(s).

    Returns
    -------
    float or numpy array (matching T) of dG in J.
    """
    T = np.asarray(T, dtype=float)
    scalar = T.ndim == 0
    Tf = np.atleast_1d(T)

    # Standard-state reaction enthalpy and entropy at 298.15 K.
    dH298 = sum(nu * PHASE[p][1] for p, nu in stoich.items())
    dS298 = sum(nu * PHASE[p][2] for p, nu in stoich.items())

    # Shared integration grid from 298.15 K up to the hottest target T.
    Tmax = max(float(Tf.max()), T_REF)
    grid = np.linspace(T_REF, Tmax, 1500)

    # Stoichiometric dCp(T) along the grid.
    dCp = np.zeros_like(grid)
    for p, nu in stoich.items():
        dCp = dCp + nu * PHASE[p][0](grid)

    # Cumulative trapezoid integrals: INT dCp dT and INT dCp/T dT.
    def _cumtrap(y, x):
        dx = np.diff(x)
        seg = 0.5 * (y[1:] + y[:-1]) * dx
        out = np.empty_like(x)
        out[0] = 0.0
        out[1:] = np.cumsum(seg)
        return out

    H_cum = _cumtrap(dCp, grid)  # INT_298^grid dCp dT
    S_cum = _cumtrap(dCp / grid, grid)  # INT_298^grid dCp/T dT

    # Interpolate the increments at the requested temperatures.
    Hinc = np.interp(Tf, grid, H_cum)
    Sinc = np.interp(Tf, grid, S_cum)

    dG = (dH298 + Hinc) - Tf * (dS298 + Sinc)
    return dG[0] if scalar else dG


def log_fS2(stoich, nu_S2, T):
    """log10 f(S2) along a boundary (dimensionless), vectorised over T.

        log10 f(S2) = dG(T) / (nu_S2 * ln(10) * R * T)

    nu_S2 is the POSITIVE magnitude of S2 consumed per reaction as written.
    """
    return reaction_dG(stoich, T) / (nu_S2 * LN10 * R * np.asarray(T, dtype=float))


# --- Boundary reactions ----------------------------------------------------
# Pyrrhotite composition Fe(X)S with X = 0.877 (the GS1992 4C phase).
X = 0.877

# Fe + 0.5 S2 -> Fe0.877S  (iron / pyrrhotite boundary; nu_S2 = 0.5)
RX_FE_PO = {"Fe": -X, "S2": -0.5, "Po": 1}
NU_FE_PO = 0.5

# (1/X) Fe0.877S + nu S2 -> FeS2  (pyrrhotite / pyrite),
#   a = 1/X (mol Po giving 1 Fe), nu = (2 - a)/2  S2 to make FeS2.
_a_po_py = 1.0 / X
_nu_po_py = (2.0 - _a_po_py) / 2.0
RX_PO_PY = {"Po": -_a_po_py, "S2": -_nu_po_py, "FeS2": 1}
NU_PO_PY = _nu_po_py

# (3/X) Fe0.877S + nu S2 -> Fe3S4  (pyrrhotite / greigite),
#   a = 3/X (mol Po giving 3 Fe), nu = (4 - a)/2  S2 to make Fe3S4.
_a_po_gr = 3.0 / X
_nu_po_gr = (4.0 - _a_po_gr) / 2.0
RX_PO_GR = {"Po": -_a_po_gr, "S2": -_nu_po_gr, "Gr": 1}
NU_PO_GR = _nu_po_gr

# Fe3S4 + S2 -> 3 FeS2  (greigite / pyrite; nu_S2 = 1)
RX_GR_PY = {"Gr": -1, "S2": -1, "FeS2": 3}
NU_GR_PY = 1

# Published greigite enthalpy uncertainty, per Fe3S4 (Subramani 2020:
# +-7300 J/mol-FeS1.33 -> x3 for the Fe3S4 formula unit).
GR_DHF_UNCERTAINTY_PER_Fe3S4 = 3 * 7300


def set_greigite_dHf_offset(off):
    """Shift greigite dHf298 (per Fe3S4) by `off` J, keeping S298 unchanged.

    Used to draw the admissible band around the central greigite enthalpy.
    """
    PHASE["Gr"] = (cp_Gr, 3 * (-144100) + off, PHASE["Gr"][2])


def greigite_boundaries(T, offset_per_Fe3S4=0.0):
    """Return (log f(S2) Po/Gr, log f(S2) Gr/Py) at temperature(s) T.

    Applies the given greigite-enthalpy offset (per Fe3S4) for the duration
    of the calculation, then resets it so the module state is unchanged.
    """
    set_greigite_dHf_offset(offset_per_Fe3S4)
    try:
        po_gr = log_fS2(RX_PO_GR, NU_PO_GR, T)
        gr_py = log_fS2(RX_GR_PY, NU_GR_PY, T)
    finally:
        set_greigite_dHf_offset(0.0)
    return po_gr, gr_py


if __name__ == "__main__":
    Ts = np.array([300.0, 500.0])
    po_gr, gr_py = greigite_boundaries(Ts, 0.0)
    print("T            =", Ts)
    print("log fS2 Po/Gr =", po_gr, " (expect ~ -54.6 / -28.4)")
    print("log fS2 Gr/Py =", gr_py, " (expect ~ -21.9 /  -7.1)")
    decomp = (
        reaction_dG({"Gr": -1 / 3.0, "Po": 0.8886, "FeS2": 0.221}, np.array([298.15]))
        / 1000.0
    )
    print("greigite decomposition dG (kJ) =", decomp, " (expect ~ +14)")
