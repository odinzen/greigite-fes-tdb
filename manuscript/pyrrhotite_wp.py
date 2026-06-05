"""Waldner & Pelton (2005) continuous Fe(1-x)S pyrrhotite solution model.

Reference
---------
WP2005 : P. Waldner & A.D. Pelton, "Thermodynamic Modeling of the Fe-S
         System", J. Phase Equilib. Diffus. 26 (2005) 23-38.

The pyrrhotite phase is modeled as a two-sublattice (Fe,Va)1(S)1 compound
energy formalism phase (their Eq 9-10).  Per mole of formula unit = per
mole of S, the metal sublattice carries site fractions y_Fe and y_Va with
y_Fe + y_Va = 1.  The composition is Fe(1-x)S with x = y_Va, and the mole
fraction of sulfur is  X_S = 1 / (1 + y_Fe).

G_po(y_Fe,T) = y_Fe*GFeS + y_Va*GVaS
             + R*T*(y_Fe*ln y_Fe + y_Va*ln y_Va)
             + y_Fe*y_Va*L
with  L = -225830.67 + 26.359*T  (J/mol).

GFeS and GVaS are END-MEMBER apparent Gibbs energies, referenced to the
elements, computed the same way reactions.py builds apparent-G:

    G(T) = dHf298 + INT_298^T Cp dT - T*(S298 + INT_298^T Cp/T dT)

Sulfur fugacity (KEY relation, their model):
    mu_S    = G_po - y_Fe * dG_po/dy_Fe
    dG/dyFe = GFeS - GVaS + R*T*(ln y_Fe - ln y_Va) + (1 - 2*y_Fe)*L
    mu_S2(gas) = 2*mu_S = G_S2_gas(T) + R*T*ln f(S2)
    log10 f(S2) = (2*mu_S - G_S2_gas(T)) / (R*T*ln10)

Units: J/mol, K; log f(S2) is log10, dimensionless.
"""

import numpy as np

import thermo_data as td

R = 8.314462618
LN10 = np.log(10.0)
T_REF = 298.15

# Interaction parameter L (J/mol), WP2005 Eq 10.
L_A = -225830.67
L_B = 26.359


def L_param(T):
    """Regular-solution interaction parameter L(T) (J/mol)."""
    return L_A + L_B * T


# ---------------------------------------------------------------------------
# End-member heat capacities (J/mol/K)
# ---------------------------------------------------------------------------
def cp_FeS_endmember(T):
    """Cp of the (Fe)1(S)1 end-member (WP2005), piecewise in T (J/mol/K).

    Includes the 590 K transition only as an enthalpy/entropy increment of
    290 J/mol (handled in the integrators), NOT in Cp itself.
    """
    T = np.asarray(T, dtype=float)
    out = np.empty_like(T)
    # 298-420
    m = (T >= 298.15) & (T < 420.0)
    out[m] = 2437.135 - 9.901929 * T[m] + 0.01156762 * T[m] ** 2 - 41123870 * T[m] ** -2
    # 420-440
    m = (T >= 420.0) & (T < 440.0)
    out[m] = 83.0
    # 440-590
    m = (T >= 440.0) & (T < 590.0)
    out[m] = 344.055 - 1.1307420 * T[m] + 0.001173304 * T[m] ** 2 - 13870 * T[m] ** -2
    # 590-1200
    m = (T >= 590.0) & (T < 1200.0)
    out[m] = 36.401 + 0.0236417 * T[m] - 5.53585e-6 * T[m] ** 2 + 3740990 * T[m] ** -2
    # 1200-1400
    m = (T >= 1200.0) & (T <= 1400.0)
    out[m] = 47.203 + 0.0082851 * T[m] + 3.6613e-7 * T[m] ** 2 + 2482690 * T[m] ** -2
    # below 298.15 (shouldn't occur) clamp to the 298-420 expression
    m = T < 298.15
    out[m] = 2437.135 - 9.901929 * T[m] + 0.01156762 * T[m] ** 2 - 41123870 * T[m] ** -2
    return out


# 590 K transition enthalpy of the (Fe)1(S)1 end-member (J/mol).
FES_TRANS_T = 590.0
FES_TRANS_H = 290.0


def cp_VaS_endmember(T):
    """Cp of the (Va)1(S)1 end-member = Cp of orthorhombic sulfur.

    The SGTE/Dinsdale orthorhombic-S Cp is not reproduced from memory with
    confidence here, so we use the DOCUMENTED linear approximation
        Cp_S(orth) ~ 14.795 + 0.02436*T   (J/mol/K)
    (flagged in USED_S_ORTH_APPROX).  y_Va is small except for the most
    S-rich pyrrhotite, so this term is second-order.
    """
    T = np.asarray(T, dtype=float)
    return 14.795 + 0.02436 * T


USED_S_ORTH_APPROX = True  # we used the documented linear S(orth) Cp


# End-member standard-state data (WP2005).
FES_dHf298 = -96291.00
FES_S298 = 69.429
VAS_dHf298 = 140049.39
VAS_S298 = 32.054


# ---------------------------------------------------------------------------
# Apparent Gibbs-energy integrator (general)
# ---------------------------------------------------------------------------
def _apparent_G(cp_func, dHf298, S298, T, trans_T=None, trans_H=0.0):
    """Apparent G(T) = dHf298 + INT Cp dT - T*(S298 + INT Cp/T dT).

    Vectorised over T.  If a first-order transition (trans_T, trans_H) is
    given, its enthalpy is added to the H integral above trans_T, and
    trans_H/trans_T is added to the S integral above trans_T.
    """
    T = np.asarray(T, dtype=float)
    scalar = T.ndim == 0
    Tf = np.atleast_1d(T)
    Tmax = max(float(Tf.max()), T_REF)
    grid = np.linspace(T_REF, Tmax, 4000)
    cp = cp_func(grid)

    def _cumtrap(y, x):
        dx = np.diff(x)
        seg = 0.5 * (y[1:] + y[:-1]) * dx
        out = np.empty_like(x)
        out[0] = 0.0
        out[1:] = np.cumsum(seg)
        return out

    H_cum = _cumtrap(cp, grid)
    S_cum = _cumtrap(cp / grid, grid)

    Hinc = np.interp(Tf, grid, H_cum)
    Sinc = np.interp(Tf, grid, S_cum)

    if trans_T is not None and trans_H:
        above = Tf >= trans_T
        Hinc = Hinc + np.where(above, trans_H, 0.0)
        Sinc = Sinc + np.where(above, trans_H / trans_T, 0.0)

    G = (dHf298 + Hinc) - Tf * (S298 + Sinc)
    return G[0] if scalar else G


def GFeS(T):
    """Apparent G of the (Fe)1(S)1 end-member (J/mol)."""
    return _apparent_G(
        cp_FeS_endmember,
        FES_dHf298,
        FES_S298,
        T,
        trans_T=FES_TRANS_T,
        trans_H=FES_TRANS_H,
    )


def GVaS(T):
    """Apparent G of the (Va)1(S)1 end-member (J/mol)."""
    return _apparent_G(cp_VaS_endmember, VAS_dHf298, VAS_S298, T)


def G_S2_gas(T):
    """Apparent G of S2(g) (J/mol), from td.S2_g Shomate Cp."""
    return _apparent_G(
        lambda TT: td.cp_shomate(td.S2_g, TT), td.S2_g["dHf298"], td.S2_g["S298"], T
    )


# ---------------------------------------------------------------------------
# Pyrrhotite solution thermodynamics
# ---------------------------------------------------------------------------
def G_po(y_Fe, T):
    """Molar Gibbs energy of pyrrhotite per mole formula (=per mole S)."""
    y_Fe = np.asarray(y_Fe, dtype=float)
    y_Va = 1.0 - y_Fe
    gfes = GFeS(T)
    gvas = GVaS(T)
    Lv = L_param(T)
    mix = R * T * (_xlnx(y_Fe) + _xlnx(y_Va))
    return y_Fe * gfes + y_Va * gvas + mix + y_Fe * y_Va * Lv


def _xlnx(y):
    """y*ln(y) with the y->0 limit (0) handled."""
    y = np.asarray(y, dtype=float)
    out = np.zeros_like(y)
    nz = y > 0.0
    out[nz] = y[nz] * np.log(y[nz])
    return out


def dG_po_dyFe(y_Fe, T):
    """Derivative dG_po/dy_Fe at fixed T (J/mol)."""
    y_Fe = np.asarray(y_Fe, dtype=float)
    y_Va = 1.0 - y_Fe
    gfes = GFeS(T)
    gvas = GVaS(T)
    Lv = L_param(T)
    return gfes - gvas + R * T * (np.log(y_Fe) - np.log(y_Va)) + (1.0 - 2.0 * y_Fe) * Lv


def mu_S(y_Fe, T):
    """Partial molar Gibbs energy of S in pyrrhotite (J/mol)."""
    return G_po(y_Fe, T) - y_Fe * dG_po_dyFe(y_Fe, T)


def log_fS2(y_Fe, T):
    """log10 f(S2) in equilibrium with pyrrhotite of site fraction y_Fe."""
    T = float(T)
    return (2.0 * mu_S(y_Fe, T) - G_S2_gas(T)) / (R * T * LN10)


def X_S_of_yFe(y_Fe):
    """Mole fraction of S for a pyrrhotite of metal site fraction y_Fe."""
    return 1.0 / (1.0 + np.asarray(y_Fe, dtype=float))


def yFe_of_X_S(X_S):
    """Inverse: metal site fraction y_Fe for a target mole fraction of S."""
    X_S = np.asarray(X_S, dtype=float)
    return (1.0 - X_S) / X_S


if __name__ == "__main__":
    # quick smoke test
    for T in (973.15, 1373.15):
        for XS in (0.500, 0.510, 0.520, 0.530, 0.540):
            yfe = yFe_of_X_S(XS)
            if yfe >= 1.0:
                yfe = 1.0 - 1e-6
            print(f"T={T:.1f} XS={XS:.3f} yFe={yfe:.4f} log fS2={log_fS2(yfe, T):.3f}")


# ---------------------------------------------------------------------------
# Scalar-T caches: the end-member apparent-G integrals are expensive (each
# rebuilds a 4000-point grid).  In composition solves we hold T fixed and
# vary y_Fe, so cache G(T) per scalar temperature.
# ---------------------------------------------------------------------------
_GFES_CACHE = {}
_GVAS_CACHE = {}
_GS2_CACHE = {}


def GFeS_cached(T):
    key = round(float(T), 6)
    if key not in _GFES_CACHE:
        _GFES_CACHE[key] = float(GFeS(T))
    return _GFES_CACHE[key]


def GVaS_cached(T):
    key = round(float(T), 6)
    if key not in _GVAS_CACHE:
        _GVAS_CACHE[key] = float(GVaS(T))
    return _GVAS_CACHE[key]


def G_S2_gas_cached(T):
    key = round(float(T), 6)
    if key not in _GS2_CACHE:
        _GS2_CACHE[key] = float(G_S2_gas(T))
    return _GS2_CACHE[key]


def G_po_fast(y_Fe, T, gfes, gvas):
    y_Fe = float(y_Fe)
    y_Va = 1.0 - y_Fe
    Lv = L_param(T)
    mix = (
        R
        * T
        * (
            (y_Fe * np.log(y_Fe) if y_Fe > 0 else 0.0)
            + (y_Va * np.log(y_Va) if y_Va > 0 else 0.0)
        )
    )
    return y_Fe * gfes + y_Va * gvas + mix + y_Fe * y_Va * Lv


def dG_po_dyFe_fast(y_Fe, T, gfes, gvas):
    y_Fe = float(y_Fe)
    y_Va = 1.0 - y_Fe
    Lv = L_param(T)
    return gfes - gvas + R * T * (np.log(y_Fe) - np.log(y_Va)) + (1.0 - 2.0 * y_Fe) * Lv


def mu_S_fast(y_Fe, T, gfes, gvas):
    return G_po_fast(y_Fe, T, gfes, gvas) - y_Fe * dG_po_dyFe_fast(y_Fe, T, gfes, gvas)
