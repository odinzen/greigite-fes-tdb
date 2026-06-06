"""Publication figure for K. Lilova: Fe-S log f(S2)-T diagram built on the
VALIDATED continuous Waldner-Pelton (2005) (Fe,Va)1S1 pyrrhotite model.

Four fields by increasing f(S2): Fe (bottom) | Fe(1-x)S pyrrhotite (large
continuous field) | Fe3S4 greigite (triangular wedge, capped at the ~525 K
metastable decomposition) | FeS2 pyrite (top).

Greigite carries the published +-7.3 kJ/mol-FeS1.33 enthalpy uncertainty
(= +-21900 J / Fe3S4); we draw the central boundary plus a light envelope
between the -7.3 (greigite-stabilized, wide) and +7.3 (destabilized) bounds.

Reuses make_fig2_continuous machinery (pyrrhotite_wp, greigite_cp).
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

import pyrrhotite_wp as wp
import greigite_cp as g
import make_fig2_continuous as m

R = wp.R
LN10 = wp.LN10
DECOMP = g.GREIGITE_DECOMP_T  # ~530 K metastable decomposition cap
GR_PINCH = 525.0  # documented ~525 K greigite pinch/cap
GR_OFFSET = 3 * 7300  # +-7.3 kJ/mol-FeS1.33 -> +-21900 J/Fe3S4


def G_Fe3S4_off(T, off):
    """Apparent G of Fe3S4 with dHf offset `off` (J / Fe3S4)."""
    grid = np.linspace(wp.T_REF, max(float(T), wp.T_REF), 150)
    cp = 3.0 * g.cp_debye_einstein(grid)
    dx = np.diff(grid)
    Hinc = np.sum(0.5 * (cp[1:] + cp[:-1]) * dx)
    Sinc = np.sum(0.5 * (cp[1:] / grid[1:] + cp[:-1] / grid[:-1]) * dx)
    return (-432300.0 + off + Hinc) - float(T) * (3 * 71.334 + Sinc)


def fields_at_T(T, off):
    """Return (lower, upper, po_gr, gr_py) log f(S2) at T for greigite dHf
    offset `off`.  lower=Fe-saturated po; upper=pyrite-saturated po;
    po_gr=pyrrhotite/greigite; gr_py=greigite/pyrite."""
    gfes = wp.GFeS_cached(T)
    gvas = wp.GVaS_cached(T)
    gs2 = wp.G_S2_gas_cached(T)
    gfe_bcc = m.G_Fe_bcc(T)
    gpy = m.G_FeS2_wp(T)
    gr = G_Fe3S4_off(T, off)
    # Fe-saturated pyrrhotite lower edge
    lfs2_fesat = (gfes - gfe_bcc - 0.5 * gs2) / (0.5 * LN10 * R * T)
    ylo = m._bisect(
        lambda y: m.log_fS2_from_y(y, T, gfes, gvas, gs2) - lfs2_fesat, 1e-6, 1 - 1e-9
    )
    # pyrite-saturated pyrrhotite upper edge
    yhi = m._bisect(
        lambda y: (
            m.mu_Fe_fast(y, T, gfes, gvas) + 2.0 * wp.mu_S_fast(y, T, gfes, gvas) - gpy
        ),
        1e-6,
        1 - 1e-9,
    )
    lower = m.log_fS2_from_y(ylo, T, gfes, gvas, gs2) if ylo else np.nan
    upper = m.log_fS2_from_y(yhi, T, gfes, gvas, gs2) if yhi else np.nan
    ysat = yhi if yhi else 0.85
    po_gr = (gr - 3.0 * m.mu_FeS_fast(ysat, T, gfes, gvas) - 0.5 * gs2) / (
        0.5 * LN10 * R * T
    )
    gr_py = (3.0 * gpy - gr - gs2) / (1.0 * LN10 * R * T)
    return lower, upper, po_gr, gr_py


if __name__ == "__main__":
    Ts = np.linspace(300.0, 600.0, 121)
    lower = np.empty_like(Ts)
    upper = np.empty_like(Ts)
    pg_c = np.empty_like(Ts)
    gp_c = np.empty_like(Ts)  # central
    pg_lo = np.empty_like(Ts)
    gp_lo = np.empty_like(Ts)  # -7.3 (greigite wide)
    pg_hi = np.empty_like(Ts)
    gp_hi = np.empty_like(Ts)  # +7.3 (greigite narrow)
    for i, T in enumerate(Ts):
        lower[i], upper[i], pg_c[i], gp_c[i] = fields_at_T(T, 0.0)
        _, _, pg_lo[i], gp_lo[i] = fields_at_T(T, -GR_OFFSET)
        _, _, pg_hi[i], gp_hi[i] = fields_at_T(T, +GR_OFFSET)

    # Greigite exists only up to its metastable decomposition (~525-530 K).
    cap = Ts <= GR_PINCH

    fig, ax = plt.subplots(figsize=(9.2, 7.0))
    ax.set_facecolor("white")
    fig.patch.set_facecolor("white")

    # y-range chosen to frame all four fields + the greigite wedge.
    YLO, YHI = -55.0, 2.0

    # --- Fe (bottom): everything below the Fe-saturated pyrrhotite edge ---
    ax.fill_between(Ts, YLO, lower, color="#d9d9d9", alpha=0.9, zorder=0)
    # --- Pyrrhotite continuous field: Fe-sat (lower) to pyrite-sat (upper) ---
    ax.fill_between(Ts, lower, upper, color="#9ecae1", alpha=0.75, zorder=1)
    # --- Pyrite (top): above the greigite/pyrite line (capped) and above
    #     pyrite-saturated pyrrhotite where greigite is gone ---
    top_edge = np.where(cap, np.maximum(gp_c, upper), upper)
    ax.fill_between(Ts, top_edge, YHI, color="#fdd0a2", alpha=0.85, zorder=1)

    # --- Greigite wedge: ENVELOPE between -7.3 (wide) and +7.3 (narrow) ---
    # Light envelope: union of admissible greigite fields over the +-7.3 band.
    env_lo = np.minimum(pg_c, pg_lo)  # lowest po/gr across the band
    env_hi = np.maximum(gp_c, gp_lo)  # highest gr/py across the band
    ax.fill_between(
        Ts,
        env_lo,
        env_hi,
        where=cap & (env_hi > env_lo),
        color="#74c476",
        alpha=0.45,
        zorder=2,
        label=r"greigite Fe$_3$S$_4$ ($\pm$7.3 kJ envelope)",
    )
    # Central greigite field: darker core wedge
    ax.fill_between(
        Ts, pg_c, gp_c, where=cap & (gp_c > pg_c), color="#238b45", alpha=0.9, zorder=3
    )

    # Boundary lines
    ax.plot(Ts, lower, "-", color="#08519c", lw=1.6, zorder=4)
    ax.plot(Ts, upper, "-", color="#08519c", lw=1.6, zorder=4)
    ax.plot(Ts[cap], gp_c[cap], "-", color="#006d2c", lw=1.8, zorder=5)
    ax.plot(Ts[cap], pg_c[cap], "-", color="#006d2c", lw=1.8, zorder=5)
    # Dashed bounding lines for the +-7.3 band on greigite boundaries
    ax.plot(Ts[cap], gp_lo[cap], "--", color="#006d2c", lw=0.9, alpha=0.7, zorder=4)
    ax.plot(Ts[cap], pg_lo[cap], "--", color="#006d2c", lw=0.9, alpha=0.7, zorder=4)

    # --- ~525 K greigite pinch marker ---
    yp = float(np.interp(GR_PINCH, Ts, env_hi))
    ax.axvline(GR_PINCH, color="#444444", lw=1.0, ls=":", zorder=6)
    ax.annotate(
        "greigite pinch /\nmetastable decomp.\n~525 K",
        xy=(GR_PINCH, yp),
        xytext=(540, -30),
        fontsize=9.5,
        color="#222222",
        ha="left",
        va="center",
        arrowprops=dict(arrowstyle="->", color="#444444", lw=1.0),
        zorder=8,
    )

    # --- Field labels INSIDE the fields ---
    ax.text(
        450,
        -50,
        "Fe  (bcc)",
        fontsize=17,
        weight="bold",
        color="#525252",
        ha="center",
        va="center",
        zorder=7,
    )
    ax.text(
        440,
        -28,
        r"Fe$_{1-x}$S  pyrrhotite",
        fontsize=18,
        weight="bold",
        color="#08306b",
        ha="center",
        va="center",
        zorder=7,
    )
    ax.text(
        415,
        0.6,
        r"FeS$_2$  pyrite",
        fontsize=17,
        weight="bold",
        color="#a63603",
        ha="center",
        va="center",
        zorder=7,
    )
    # greigite label on the green wedge, with a short leader
    yg = float(np.interp(400, Ts, 0.5 * (env_lo + env_hi)))
    ax.annotate(
        r"Fe$_3$S$_4$" + "\ngreigite",
        xy=(400, yg),
        xytext=(330, -32),
        fontsize=12.5,
        weight="bold",
        color="#00441b",
        ha="center",
        va="center",
        arrowprops=dict(arrowstyle="->", color="#00441b", lw=1.1),
        zorder=8,
    )

    ax.set_xlim(300, 600)
    ax.set_ylim(YLO, YHI)
    ax.set_xlabel("Temperature  [K]", fontsize=13)
    ax.set_ylabel(r"log$_{10}$ $f$(S$_2$)   [bar]", fontsize=13)
    ax.set_title(
        "Fe-S with continuous pyrrhotite (Waldner-Pelton 2005) + greigite",
        fontsize=14,
        weight="bold",
        pad=12,
    )
    ax.grid(True, which="major", color="#e8e8e8", lw=0.6, zorder=0)
    ax.tick_params(labelsize=11)
    ax.legend(loc="lower right", fontsize=9, framealpha=0.9)

    foot = (
        "Pyrrhotite Fe$_{1-x}$S = validated Waldner-Pelton 2005 solution "
        "model (reproduces their Fig 9).  Greigite Fe$_3$S$_4$ boundary "
        "provisional\n(metastable; shaded = $\\pm$7.3 kJ/mol).  "
        "Sources: Waldner-Pelton 2005; Shumway 2022; Subramani 2020."
    )
    fig.text(0.012, 0.012, foot, fontsize=7.6, color="#555555", va="bottom")

    fig.subplots_adjust(left=0.085, right=0.975, top=0.93, bottom=0.13)
    fig.savefig(str(FIG / "fig2_kristina_continuous.png"), dpi=150, facecolor="white")
    print("saved fig2_kristina_continuous.png")

    # Report numbers: greigite triangle corners (central model)
    print(
        f"base T=300 K: po/gr={pg_c[0]:.2f}  gr/py={gp_c[0]:.2f}  "
        f"width={gp_c[0] - pg_c[0]:.2f}"
    )
    ip = np.argmin(np.abs(Ts - GR_PINCH))
    print(f"pinch T={Ts[ip]:.0f} K: po/gr={pg_c[ip]:.2f}  gr/py={gp_c[ip]:.2f}")
    print(f"-7.3 base T=300: po/gr={pg_lo[0]:.2f} gr/py={gp_lo[0]:.2f}")
