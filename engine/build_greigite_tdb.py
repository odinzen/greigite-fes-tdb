#!/usr/bin/env python
"""PHASE B.4 — build the greigite (Fe3S4) CALPHAD description from MEASURED data.

Deterministic + rerunnable. Reads the greigite measurements (dHf, S°, the
Debye-Einstein Cp coefficients) from the committed campaign manifest
(provenance_manifest.json, via manifest.py) — NOT hardcoded literals — and
writes a temperature-dependent GREIGITE line compound into the Dilner-2015 base.

The Dilner-2015 base TDB is read from artifacts/tdb/ and auto-fetched there (per
the manifest's TDBDB record) on first run if missing.

Thermodynamic construction (per the agreed formula):
    G_greigite(T) = ΔfH°(298) + ∫₂₉₈^T Cp dT − T·[ S°(298) + ∫₂₉₈^T (Cp/T) dT ]
  inputs, per FeS1.33 (×3 → Fe3S4):
    ΔfH°(298) = -144.1 kJ/mol           Subramani 2020 (id 389, Table 1)
    S°(298)   = 71.334 J/mol/K          Shumway 2022  (id 384, Table 5)
    Cp(T)     = m·D(Θ_D/T)+n·E(Θ_E/T)+A·T   Shumway 2022 (id 384, Table 4)
                m=0.9740 Θ_D=237.85 n=1.1329 Θ_E=405.26 A=0.0398
  The Debye-Einstein (DE) high-T Cp is used (NOT the cubic / mid-T polynomial,
  which diverges >300 K). Cp over 298–530 K is fit to an SGTE polynomial; the
  greigite G is written as +3*GHSERFE +4*GHSERSS +<remainder fit>, reusing
  Dilner's SER element references so the phase sits on the same scale as the base.

Outputs (artifacts/tdb/):
    fes_greigite_v1.tdb            nominal ΔfH
    fes_greigite_v1_dHf_lo.tdb     ΔfH − 3×7.3 kJ (lower limit, more negative)
    fes_greigite_v1_dHf_hi.tdb     ΔfH + 3×7.3 kJ (upper limit, less negative)
    artifacts/build.log            every consumed value + fit coeffs + residual
    artifacts/build_report.json    inputs, fit coeffs, output TDB sha256s

±7.3 kJ/mol is per FeS1.33 (absolute, drop-solution) → ±21.9 kJ/mol per Fe3S4.
Cp 298–530 K is an extrapolation of a fit valid to ~303 K, through the ~530 K
kinetic decomposition of greigite; this shifts fugacity boundaries <0.1 log unit.
"""

from __future__ import annotations

import hashlib
import json
import logging
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from pycalphad import Database
from symengine import Symbol

from manifest import Manifest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TDB_DIR = ROOT / "artifacts" / "tdb"
BASE_TDB_ID = "tdb:dilner2015_femns"  # manifest artifact id (Dilner 2015 Fe-S base)
LOG = ROOT / "artifacts" / "build.log"
REPORT = ROOT / "artifacts" / "build_report.json"

MANIFEST = Manifest.load()
BASE = MANIFEST.tdb(BASE_TDB_ID).dest_path  # artifacts/tdb/<itemid>.tdb


def ensure_base() -> Path:
    """Fetch the Dilner-2015 base TDB into artifacts/tdb/ (per the manifest) if missing."""
    return MANIFEST.resolve_tdb(BASE_TDB_ID)


R = 8.314462618
T0 = 298.15
TFIT_LO, TFIT_HI = 298.15, 530.0  # SGTE fit window (instruction)
TDB_HI = 600.0  # TDB validity ceiling (>530 flagged extrapolation)

LOG.parent.mkdir(parents=True, exist_ok=True)  # artifacts/ for the log + reports
logger = logging.getLogger("greigite_build")
logger.setLevel(logging.INFO)
logger.handlers.clear()
for _h in (logging.FileHandler(LOG, mode="w"), logging.StreamHandler()):
    _h.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(_h)


# --------------------------------------------------------------------------- #
# 1. Pull measured values (with provenance) from the campaign manifest
# --------------------------------------------------------------------------- #
def load_values() -> dict:
    """Read the greigite measurements from provenance_manifest.json.

    Returns rows shaped {parsed, unit, source, id} so the build + report logic
    is source-agnostic. dHf, S°, the internal-check Cp/S, and the Debye-Einstein
    Cp coefficients all trace back to a typed Measurement node.
    """
    dHf = MANIFEST.measurement("meas:subramani2020:dHf_greigite")
    s298 = MANIFEST.measurement("meas:shumway2022:S298_greigite")
    cp298 = MANIFEST.measurement("meas:shumway2022:Cp298_greigite")
    s300 = MANIFEST.measurement("meas:shumway2022:S300_check")
    de = MANIFEST.measurement("meas:shumway2022:cp_debye_einstein")

    def src(m):
        return f"{m.paper} {m.locator_str}".strip()

    def row(parsed, unit, m):
        return {"parsed": parsed, "unit": unit, "source": src(m), "id": m.id}

    v = {
        "dHf_FeS133": row(dHf.value, dHf.unit, dHf),
        "dHf_unc": row(
            dHf.uncertainty_value, dHf.uncertainty.get("unit", dHf.unit), dHf
        ),
        "S298": row(s298.value, s298.unit, s298),
        "Cp298": row(cp298.value, cp298.unit, cp298),
        "S300_check": row(s300.value, s300.unit, s300),
        "m": row(de.coeff("m"), "mol", de),
        "ThetaD": row(de.coeff("ThetaD_K"), "K", de),
        "n": row(de.coeff("n"), "mol", de),
        "ThetaE": row(de.coeff("ThetaE_K"), "K", de),
        "A": row(de.coeff("A_per_K"), "K^-1", de),
    }
    logger.info("=== consumed measured values (source: provenance_manifest.json) ===")
    for k, r in v.items():
        logger.info(
            "  %-12s = %-10s %-8s  <- %s  [%s]",
            k,
            r["parsed"],
            r["unit"],
            r["source"],
            r["id"],
        )
    logger.info("")
    return v


# --------------------------------------------------------------------------- #
# 2. Debye-Einstein Cp (per FeS1.33) and the internal S(300) recovery check
# --------------------------------------------------------------------------- #
def make_cp(m, ThetaD, n, ThetaE, A):
    def _dint(t):  # Debye integrand, overflow-safe via exp(-t)
        if t < 1e-8:
            return t * t  # t^4 e^t/(e^t-1)^2 -> t^2 as t->0
        em = math.exp(-t)
        return t**4 * em / (1.0 - em) ** 2

    def debye(theta, T):
        x = theta / T
        integ, _ = quad(_dint, 0, x, limit=200)
        return 9 * R * (T / theta) ** 3 * integ  # -> 3R as T->inf (per oscillator)

    def einstein(theta, T):
        x = theta / T
        em = math.exp(-x)  # overflow-safe: -> 0 for large x
        return 3 * R * x * x * em / (1.0 - em) ** 2

    def cp_FeS133(T):
        return m * debye(ThetaD, T) + n * einstein(ThetaE, T) + A * T

    return cp_FeS133


def internal_checks(cp_FeS133, cp298_expected, s300_expected):
    cp298 = cp_FeS133(T0)
    # S(300) from the DE fit integrated from ~0 K (DE form is valid to 0 K).
    s300, _ = quad(lambda T: cp_FeS133(T) / T, 0.5, 300.0, limit=200)
    logger.info("=== internal checks (per FeS1.33) ===")
    logger.info(
        "  Cp(298.15) DE   = %.3f  vs paper %.3f   (diff %.3f)",
        cp298,
        cp298_expected,
        cp298 - cp298_expected,
    )
    logger.info(
        "  S(300) ∫DE/T 0→300 = %.3f  vs paper %.3f (diff %.3f)",
        s300,
        s300_expected,
        s300 - s300_expected,
    )
    logger.info("")
    return cp298, s300


# --------------------------------------------------------------------------- #
# 3. SGTE element references (evaluate Dilner GHSERFE / GHSERSS numerically)
# --------------------------------------------------------------------------- #
def make_ghser(db):
    symmap = {Symbol(k): val for k, val in db.symbols.items()}

    def ev(name, T):
        e = db.symbols[name]
        for _ in range(25):
            ne = e.xreplace(symmap)
            if ne == e:
                break
            e = ne
        return float(e.xreplace({Symbol("T"): T, Symbol("P"): 101325.0}))

    return ev


# --------------------------------------------------------------------------- #
# 4. Build greigite G(T) per Fe3S4 and fit the SER remainder to an SGTE poly
# --------------------------------------------------------------------------- #
SGTE_BASIS = [
    ("1", lambda T: np.ones_like(T)),
    ("T", lambda T: T),
    ("T*LN(T)", lambda T: T * np.log(T)),
    ("T**2", lambda T: T**2),
    ("T**(-1)", lambda T: 1.0 / T),
    ("T**3", lambda T: T**3),
]


def build_one(dHf_Fe3S4, S298_Fe3S4, cp_FeS133, ev):
    """Return (coeffs dict, max_resid_J, grid diagnostics) for one ΔfH choice."""
    Tg = np.linspace(TFIT_LO, TFIT_HI, 120)
    cp_Fe3S4 = np.array([3.0 * cp_FeS133(T) for T in Tg])  # ×3: FeS1.33 -> Fe3S4

    # H(T)-H(298) and S(T)-S(298) by cumulative integration of Cp / (Cp/T)
    def cumint(y):
        out = np.zeros_like(Tg)
        out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(Tg))
        return out

    Hrel = cumint(cp_Fe3S4)
    Srel = cumint(cp_Fe3S4 / Tg)
    G_abs = dHf_Fe3S4 + Hrel - Tg * (S298_Fe3S4 + Srel)  # greigite absolute G
    remainder = G_abs - np.array(
        [3 * ev("GHSERFE", T) + 4 * ev("GHSERSS", T) for T in Tg]
    )

    M = np.column_stack([f(Tg) for _, f in SGTE_BASIS])
    coef, *_ = np.linalg.lstsq(M, remainder, rcond=None)
    fit = M @ coef
    max_resid = float(np.max(np.abs(fit - remainder)))
    coeffs = {name: float(c) for (name, _), c in zip(SGTE_BASIS, coef)}
    return coeffs, max_resid, (Tg, G_abs, cp_Fe3S4)


def fmt_func(coeffs) -> str:
    sym = {
        "1": "",
        "T": "*T",
        "T*LN(T)": "*T*LN(T)",
        "T**2": "*T**2",
        "T**(-1)": "*T**(-1)",
        "T**3": "*T**3",
    }
    terms = []
    for name, _ in SGTE_BASIS:
        c = coeffs[name]
        terms.append(f"{c:+.8E}{sym[name]}")
    return "".join(terms)


def write_tdb(
    out: Path, base_text: str, func_name: str, expr: str, dHf_Fe3S4: float, tag: str
) -> str:
    block = f"""
$ ============================================================================
$ GREIGITE (Fe3S4) — measured-data build ({tag})
$ Source: Subramani 2020 (corpus id 389, Table 1) ΔfH; Shumway 2022 (id 384,
$ Tables 4/5) S°298 + Debye-Einstein Cp. Built by build_greigite_tdb.py from
$ provenance_manifest.json. ΔfH(298, Fe3S4) used = {dHf_Fe3S4:.1f} J/mol.
$ G = +3 GHSERFE +4 GHSERSS + {func_name} (SER-referenced; Cp 298–530 K fit,
$ extrapolated to {TDB_HI:.0f} K — >530 K is through greigite decomposition).
$ ============================================================================
 FUNCTION {func_name}  2.98150E+02  {expr};  {TDB_HI:.5E}  N !
 PHASE GREIGITE  %  2 3 4 !
 CONSTITUENT GREIGITE  :FE : S :  !
 PARAMETER G(GREIGITE,FE:S;0)  2.98150E+02  +3*GHSERFE#+4*GHSERSS#+{func_name}#;  {TDB_HI:.5E}  N REF0 !
"""
    text = base_text.rstrip() + "\n" + block
    out.write_text(text)
    return hashlib.sha256(text.encode()).hexdigest()


def main() -> None:
    v = load_values()
    dHf_FeS133 = v["dHf_FeS133"]["parsed"] * 1000.0  # kJ->J
    dHf_unc_FeS133 = v["dHf_unc"]["parsed"] * 1000.0
    S298_FeS133 = v["S298"]["parsed"]
    cp_FeS133 = make_cp(
        v["m"]["parsed"],
        v["ThetaD"]["parsed"],
        v["n"]["parsed"],
        v["ThetaE"]["parsed"],
        v["A"]["parsed"],
    )

    cp298, s300 = internal_checks(
        cp_FeS133, v["Cp298"]["parsed"], v["S300_check"]["parsed"]
    )

    ensure_base()
    TDB_DIR.mkdir(parents=True, exist_ok=True)
    base_text = BASE.read_text()
    db = Database(str(BASE))
    ev = make_ghser(db)

    S298_Fe3S4 = 3.0 * S298_FeS133
    variants = {
        "fes_greigite_v1.tdb": (3.0 * dHf_FeS133, "GREIGITE_GF", "nominal ΔfH"),
        "fes_greigite_v1_dHf_lo.tdb": (
            3.0 * (dHf_FeS133 - dHf_unc_FeS133),
            "GREIGITE_GFLO",
            "lower limit ΔfH−3×7.3 kJ",
        ),
        "fes_greigite_v1_dHf_hi.tdb": (
            3.0 * (dHf_FeS133 + dHf_unc_FeS133),
            "GREIGITE_GFHI",
            "upper limit ΔfH+3×7.3 kJ",
        ),
    }

    report = {
        "inputs": {
            k: {
                "parsed": r["parsed"],
                "unit": r["unit"],
                "source": r["source"],
                "measurement_id": r["id"],
            }
            for k, r in v.items()
        },
        "internal_checks": {
            "Cp298_DE": cp298,
            "Cp298_paper": v["Cp298"]["parsed"],
            "S300_DE_from0": s300,
            "S300_paper": v["S300_check"]["parsed"],
        },
        "S298_Fe3S4": S298_Fe3S4,
        "fit_window_K": [TFIT_LO, TFIT_HI],
        "tdb_ceiling_K": TDB_HI,
        "variants": {},
    }

    logger.info("=== greigite G(T) build (per Fe3S4; SER = 3 GHSERFE + 4 GHSERSS) ===")
    for fname, (dHf, func, tag) in variants.items():
        coeffs, max_resid, _ = build_one(dHf, S298_Fe3S4, cp_FeS133, ev)
        expr = fmt_func(coeffs)
        sha = write_tdb(TDB_DIR / fname, base_text, func, expr, dHf, tag)
        logger.info(
            "  %-28s ΔfH(Fe3S4)=%.1f J  fit max|resid|=%.2f J  sha=%s",
            fname,
            dHf,
            max_resid,
            sha[:12],
        )
        logger.info("      %s = %s", func, expr)
        report["variants"][fname] = {
            "dHf_Fe3S4_J": dHf,
            "tag": tag,
            "function": func,
            "expr": expr,
            "fit_max_resid_J": max_resid,
            "coeffs": coeffs,
            "sha256": sha,
        }
    logger.info("")

    REPORT.write_text(json.dumps(report, indent=2))
    logger.info(
        "wrote %d TDB variants -> %s + artifacts/build_report.json",
        len(variants),
        TDB_DIR,
    )


if __name__ == "__main__":
    main()
