"""Derive the aqueous dissolution log_K for greigite (Fe3S4) on a database-consistent reference.

The porewater diagram needs a greigite dissolution constant that is built from the SAME measured
thermochemistry as the CALPHAD databases (Subramani 2020 formation enthalpy, Shumway 2022 entropy),
not from the legacy value in older aqueous databases. Two traps are checked before any number is
trusted:

  1. Reference reconciliation. A CALPHAD solid Gibbs energy sits on the SGTE element reference and
     is NOT directly mixable with CODATA/SUPCRT aqueous species. The aqueous log_K uses the CODATA
     reference throughout: dGf = dHf - T*(S_greigite - 3 S_Fe - 4 S_S) with CODATA element
     entropies. Using the raw SGTE formation term would be several log units wrong.
  2. Verifier pass. The same aqueous set reproduces the reference pyrite dissolution constant from
     Robie-Hemingway solid data to within ~0.04 log units, so the machinery is trusted for greigite.

Self-contained: every input is a cited constant (below), so this reproduces the values in
LOGK_DERIVATION.md with only the standard library.
"""
import math

R = 8.314462618            # J/mol/K
T = 298.15                 # K
LN10RT = math.log(10) * R * T   # J/mol per log10 unit = 5708.0

# --- Standard-state formation enthalpy and entropy of greigite (per Fe3S4) ---
# dHf: Subramani et al. 2020 (oxide-melt drop-solution calorimetry), CODATA scale.
# S(298): Shumway et al. 2022 low-temperature heat capacity, 3 x per-FeS1.33 value (71.704).
dHf_greigite = -431.97     # kJ/mol Fe3S4 (CODATA reference)
dHf_sigma    = 21.9        # kJ/mol, 1 sigma calorimetric uncertainty
S_greigite   = 215.11      # J/mol/K, = 3 x 71.704 (per-Fe3S4 basis confirmed)

# --- CODATA element entropies (J/mol/K) for the formation-entropy term ---
S_FE = 27.28               # bcc-Fe
S_S  = 32.05               # orthorhombic S

# --- Aqueous species formation Gibbs energies (kJ/mol), Robie-Hemingway 1995 / CODATA ---
dGf_aq = {"Fe+2": -78.90, "Fe+3": -4.60, "HS-": 12.05, "H+": 0.0, "e-": 0.0}

# --- Solid data (dHf kJ/mol, S J/mol/K) for the pyrite verifier, Robie-Hemingway 1995 ---
PYRITE = (-171.5, 52.9)    # FeS2


def dGf_from_dHf_S(dHf_kJ, S, n_Fe, n_S):
    """Formation Gibbs energy on the element reference: dGf = dHf - T*(S - sum S_elements)."""
    dS_form = S - (n_Fe * S_FE + n_S * S_S)      # J/mol/K
    return dHf_kJ - T * dS_form / 1000.0         # kJ/mol


def logK(dGr_kJ):
    """log10 K from the reaction Gibbs energy (kJ/mol)."""
    return -(dGr_kJ * 1000.0) / LN10RT


# --- Verifier: reproduce the reference pyrite dissolution constant --------------------------
# FeS2 + 2 H+ + 2 e- = Fe+2 + 2 HS-   (e- has dGf = 0 by convention)
dGf_pyrite = dGf_from_dHf_S(*PYRITE, n_Fe=1, n_S=2)
pyr_dGr = (dGf_aq["Fe+2"] + 2 * dGf_aq["HS-"]) - dGf_pyrite
print("verifier  pyrite  FeS2 + 2H+ + 2e- = Fe+2 + 2HS-")
print(f"          dGf(pyrite) = {dGf_pyrite:8.2f} kJ/mol   log_K = {logK(pyr_dGr):8.3f}"
      f"   (reference phreeqc.dat: -18.479; residual 0.04)")

# --- Greigite dissolution constant, nominal and +/- 1 sigma ---------------------------------
# Fe3S4 + 4 H+ = 2 Fe+3 + Fe+2 + 4 HS-   (congruent, non-redox)
prod = 2 * dGf_aq["Fe+3"] + dGf_aq["Fe+2"] + 4 * dGf_aq["HS-"]
print("\ngreigite  Fe3S4 + 4H+ = 2Fe+3 + Fe+2 + 4HS-")
for tag, dHf in (("nominal", dHf_greigite),
                 ("+1sigma (least stable)", dHf_greigite + dHf_sigma),
                 ("-1sigma (most stable)",  dHf_greigite - dHf_sigma)):
    dGf_g = dGf_from_dHf_S(dHf, S_greigite, n_Fe=3, n_S=4)
    print(f"  [{tag:22}] dHf={dHf:8.2f}  dGf={dGf_g:8.2f} kJ/mol  ->  log_K = {logK(prod - dGf_g):8.3f}")

print("\n  (legacy aqueous-database greigite log_K = -45.035; this derivation replaces it)")
