"""Fe-S log f(S2) vs T diagram with CONTINUOUS (variable-composition)
pyrrhotite, built on the Waldner-Pelton (2005) solution model.

At each T the pyrrhotite single-phase field spans the f(S2) range swept by
its composition, bounded by:
  * lower : Fe-saturated pyrrhotite (coexists with bcc-Fe)
  * upper : pyrite-saturated pyrrhotite (coexists with FeS2)

WP2005 consistent line compounds (Table 4):
  FeS2 pyrite : dHf = -171048.03 , S298 = 52.93 ,
                Cp(298-1350) = 72.387 + 0.0088501 T + 7.3e-10 T^2 - 1.14279e6 T^-2
Greigite Fe3S4 (from greigite_cp, 3xFeS1.33): dHf = -432300 , S298 = 3*71.334.
bcc-Fe : td.Fe_bcc (JANAF).
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TDB = ROOT / "artifacts" / "tdb"
FIG = ROOT / "artifacts" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(HERE))

import thermo_data as td
import pyrrhotite_wp as wp
import greigite_cp as g

R = wp.R
LN10 = wp.LN10


def cp_pyrite_wp(T):
    T = np.asarray(T, dtype=float)
    return 72.387 + 0.0088501 * T + 7.3e-10 * T**2 - 1.14279e6 * T**-2


_FEBCC = {}
_FES2 = {}
_FE3S4 = {}


def G_Fe_bcc(T):
    k = round(float(T), 6)
    if k not in _FEBCC:
        _FEBCC[k] = float(
            wp._apparent_G(
                lambda TT: td.cp_shomate(td.Fe_bcc, TT),
                td.Fe_bcc["dHf298"],
                td.Fe_bcc["S298"],
                T,
            )
        )
    return _FEBCC[k]


def G_FeS2_wp(T):
    k = round(float(T), 6)
    if k not in _FES2:
        _FES2[k] = float(wp._apparent_G(cp_pyrite_wp, -171048.03, 52.93, T))
    return _FES2[k]


def _apparent_G_coarse(cp_func, dHf, S298, T, npts=200):
    """Apparent G with a coarse integration grid (for smooth Cp like greigite)."""
    grid = np.linspace(wp.T_REF, max(float(T), wp.T_REF), npts)
    cp = cp_func(grid)
    dx = np.diff(grid)
    Hinc = np.sum(0.5 * (cp[1:] + cp[:-1]) * dx)
    Sinc = np.sum(0.5 * (cp[1:] / grid[1:] + cp[:-1] / grid[:-1]) * dx)
    return (dHf + Hinc) - float(T) * (S298 + Sinc)


def G_Fe3S4(T):
    k = round(float(T), 6)
    if k not in _FE3S4:
        _FE3S4[k] = float(
            _apparent_G_coarse(
                lambda TT: 3.0 * g.cp_debye_einstein(TT),
                -432300.0,
                3 * 71.334,
                T,
                npts=150,
            )
        )
    return _FE3S4[k]


def mu_Fe_fast(y, T, gfes, gvas):
    return wp.G_po_fast(y, T, gfes, gvas) + (1.0 - y) * wp.dG_po_dyFe_fast(
        y, T, gfes, gvas
    )


def mu_FeS_fast(y, T, gfes, gvas):
    return mu_Fe_fast(y, T, gfes, gvas) + wp.mu_S_fast(y, T, gfes, gvas)


def _bisect(f, a, b, tol=1e-10, nmax=200):
    fa, fb = f(a), f(b)
    if fa * fb > 0:
        return None
    for _ in range(nmax):
        m = 0.5 * (a + b)
        fm = f(m)
        if abs(b - a) < tol:
            return m
        if fa * fm <= 0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    return 0.5 * (a + b)


def log_fS2_from_y(y, T, gfes, gvas, gs2):
    return (2.0 * wp.mu_S_fast(y, T, gfes, gvas) - gs2) / (R * T * LN10)


if __name__ == "__main__":
    Ts = np.linspace(300.0, 900.0, 81)

    lower = []
    upper = []
    yfe_lo = []
    yfe_hi = []
    po_gr = []
    for T in Ts:
        gfes = wp.GFeS_cached(T)
        gvas = wp.GVaS_cached(T)
        gs2 = wp.G_S2_gas_cached(T)
        gfe_bcc = G_Fe_bcc(T)
        gpy = G_FeS2_wp(T)
        gr = G_Fe3S4(T)

        # Fe-saturated boundary: the most Fe-rich pyrrhotite is the one whose
        # solution-model log f(S2) matches the Fe(bcc)+1/2 S2 -> FeS reaction.
        # Equivalently mu_S(po) = mu_S set by bcc-Fe coexistence:
        #   2*mu_S = GFeS_endmember-referenced fS2 ... we solve log fS2(y)=target.
        lfs2_fesat = (gfes - gfe_bcc - 0.5 * gs2) / (0.5 * LN10 * R * T)
        ylo = _bisect(
            lambda y: log_fS2_from_y(y, T, gfes, gvas, gs2) - lfs2_fesat,
            1e-6,
            1.0 - 1e-9,
        )
        # pyrite-saturated: mu_Fe + 2 mu_S = G_FeS2  (FeS2 = "FeS"(po)+1/2 S2)
        yhi = _bisect(
            lambda y: (
                mu_Fe_fast(y, T, gfes, gvas)
                + 2.0 * wp.mu_S_fast(y, T, gfes, gvas)
                - gpy
            ),
            1e-6,
            1.0 - 1e-9,
        )
        lfs2_fesat_arr = lfs2_fesat
        yfe_lo.append(ylo)
        yfe_hi.append(yhi)
        lower.append(log_fS2_from_y(ylo, T, gfes, gvas, gs2) if ylo else np.nan)
        upper.append(log_fS2_from_y(yhi, T, gfes, gvas, gs2) if yhi else np.nan)

        # pyrrhotite/greigite: 3 FeS(po@py-sat) + 0.5 S2 -> Fe3S4
        ysat = yhi if yhi else 0.85
        dG = gr - 3.0 * mu_FeS_fast(ysat, T, gfes, gvas) - 0.5 * gs2
        po_gr.append(dG / (0.5 * LN10 * R * T))

    lower = np.array(lower)
    upper = np.array(upper)
    po_gr = np.array(po_gr)

    # greigite/pyrite line: Fe3S4 + S2 -> 3 FeS2
    gr_py = []
    for T in Ts:
        dG = 3.0 * G_FeS2_wp(T) - G_Fe3S4(T) - wp.G_S2_gas_cached(T)
        gr_py.append(dG / (1.0 * LN10 * R * T))
    gr_py = np.array(gr_py)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.fill_between(
        Ts,
        lower,
        upper,
        color="C0",
        alpha=0.22,
        label="pyrrhotite Fe$_{1-x}$S (continuous field)",
    )
    ax.plot(Ts, lower, "C0-", lw=1.5, label="Fe-saturated po (po / bcc-Fe)")
    ax.plot(Ts, upper, "C0--", lw=1.5, label="pyrite-saturated po (po / FeS$_2$)")
    ax.plot(Ts, gr_py, "C3-", lw=1.8, label="greigite / pyrite")
    ax.plot(
        Ts,
        po_gr,
        "C2-",
        lw=1.5,
        alpha=0.85,
        label="pyrrhotite / greigite (provisional)",
    )

    valid = (po_gr < gr_py) & (Ts <= g.GREIGITE_DECOMP_T)
    if valid.any():
        ax.fill_between(
            Ts,
            po_gr,
            gr_py,
            where=valid,
            color="C2",
            alpha=0.18,
            label="greigite Fe$_3$S$_4$ field",
        )

    ax.set_xlabel("T  [K]")
    ax.set_ylabel("log$_{10}$ f(S$_2$)  [bar]")
    ax.set_title(
        "Fe-S: continuous pyrrhotite (Waldner-Pelton 2005)\n"
        "with greigite / pyrite line compounds  (greigite provisional)"
    )
    ax.set_xlim(300, 900)
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right", fontsize=8)

    # Inset: zoom on the greigite triangular wedge (po/gr vs gr/py).
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

    axin = inset_axes(ax, width="38%", height="34%", loc="upper left", borderpad=2.2)
    axin.plot(Ts, gr_py, "C3-", lw=1.6)
    axin.plot(Ts, po_gr, "C2-", lw=1.4)
    axin.fill_between(Ts, po_gr, gr_py, where=valid, color="C2", alpha=0.30)
    # pyrite-saturated pyrrhotite upper edge for context
    axin.plot(Ts, upper, "C0--", lw=1.0)
    axin.set_xlim(300, 560)
    yvis = np.concatenate([gr_py[Ts <= 560], po_gr[Ts <= 560]])
    axin.set_ylim(np.nanmin(yvis) - 1, np.nanmax(yvis) + 1)
    axin.set_title("greigite wedge (pinches ~525 K)", fontsize=7)
    axin.tick_params(labelsize=6)
    axin.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(str(FIG / "explore_continuous_pyrrhotite.png"), dpi=130)
    print("saved explore_continuous_pyrrhotite.png")

    print("\n  T     yFe(Fe-sat) yFe(py-sat) logfS2_lo logfS2_hi  gr/py   po/gr")
    for i in range(0, len(Ts), 10):
        lo = yfe_lo[i] if yfe_lo[i] else float("nan")
        hi = yfe_hi[i] if yfe_hi[i] else float("nan")
        print(
            f"  {Ts[i]:5.0f}  {lo:.4f}     {hi:.4f}    {lower[i]:8.2f} "
            f"{upper[i]:8.2f} {gr_py[i]:7.2f} {po_gr[i]:7.2f}"
        )
    gw = np.where(valid)[0]
    if gw.size:
        print(
            f"\n  greigite wedge T={Ts[gw[0]]:.0f}..{Ts[gw[-1]]:.0f} K; "
            f"pinches near T={Ts[gw[-1]]:.0f} K (po/gr meets gr/py) -> TRIANGULAR"
        )
    else:
        print("\n  greigite NOT triangular (po/gr never below gr/py in range)")
