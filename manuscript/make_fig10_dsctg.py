#!/usr/bin/env python3
"""Fig. 10 (B&W) — Simultaneous TG-DSC of bulk greigite on heating to 600 C.

Re-plotted from the raw instrument export in data_dsc/. Left axis: DSC heat
flow (exo up). Right axis: TG mass change (%, relative to the 15.57 mg initial
mass; the startup transient < ~60 s is trimmed). Heating solid, cooling dashed.
No footnote text on the figure itself (it lives in the caption). Black & white.
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent  # = manuscript/
ROOT = HERE.parent  # = repo root
FIG = ROOT / "artifacts" / "figures"
FIG.mkdir(parents=True, exist_ok=True)
DATA = HERE / "data_dsc"
sys.path.insert(0, str(HERE))
import bw_style as bw

bw.apply()
OUT = str(FIG / "fig10_dsc_tg.png")
M0 = 15.57  # mg, initial sample mass (from the instrument header)


def parse(fn):
    txt = open(fn, "rb").read().decode("utf-16").splitlines()
    h = next(i for i, l in enumerate(txt) if l.startswith("Index"))
    rows = []
    for l in txt[h + 1 :]:
        p = l.split()
        if len(p) < 5:
            continue
        try:
            rows.append(
                [float(p[1]), float(p[2]), float(p[3]), float(p[4])]
            )  # time,T,TG,HF
        except ValueError:
            continue
    a = np.array(rows)
    a = a[a[:, 0] > 60.0]  # trim startup transient (first ~60 s)
    return a[:, 1], a[:, 2], a[:, 3]  # T, TG(mg), HF(uV)


def smooth(y, w=31):
    k = np.ones(w) / w
    s = np.convolve(y, k, mode="same")
    return s[w:-w]


def prep(fn):
    T, TG, HF = parse(fn)
    return T[31:-31], smooth(TG), smooth(HF)


Th, TGh, HFh = prep(str(DATA / "Fe3S4_bulk_heating.txt"))
Tc, TGc, HFc = prep(str(DATA / "Fe3S4_bulk_cooling.txt"))
# TG: drop the low-T startup transient (T is monotonic on heating) and convert to
# % mass change vs the post-transient baseline (mean over 90-150 C).
mTG = Th > 70.0
Th_tg, TGh = Th[mTG], TGh[mTG]
base = TGh[(Th_tg > 90) & (Th_tg < 150)].mean()
TGh_pct = (TGh - base) / M0 * 100.0

fig, axL = plt.subplots(figsize=(8.2, 5.6))
axR = axL.twinx()

# DSC (left axis)
(lh,) = axL.plot(Th, HFh, color="black", lw=1.6, label="DSC heating")
(lc,) = axL.plot(Tc, HFc, color="black", lw=1.2, ls=(0, (5, 4)), label="DSC cooling")
axL.set_xlabel("Temperature, °C")
axL.set_ylabel("DSC heat flow, µV  (exo up)")
axL.set_xlim(40, 600)

# TG (right axis)
(lt,) = axR.plot(
    Th_tg, TGh_pct, color="0.45", lw=2.0, ls=(0, (1, 1)), label="TG heating (mass)"
)
axR.set_ylabel("TG mass change, %", color="0.3")
axR.tick_params(axis="y", colors="0.3")

axL.axvline(540, color="0.6", lw=0.8, ls=":")
axL.text(
    0.03,
    0.05,
    "Fe$_3$S$_4$ — bulk",
    transform=axL.transAxes,
    fontsize=11,
    fontweight="bold",
    va="bottom",
    bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=1.0),
)

axL.legend(
    handles=[lh, lc, lt],
    loc="upper left",
    frameon=True,
    edgecolor="black",
    fontsize=9.5,
)
axL.set_title(
    "Fig. 10 — Simultaneous TG–DSC of bulk Fe$_3$S$_4$ (heating to 600 °C)",
    fontsize=12,
    fontweight="bold",
)
fig.tight_layout()
fig.savefig(OUT, dpi=200, bbox_inches="tight")
print("wrote", OUT)
print(
    "TG net at 600C: %.2f%% ; range %.2f..%.2f"
    % (TGh_pct[-1], TGh_pct.min(), TGh_pct.max())
)
