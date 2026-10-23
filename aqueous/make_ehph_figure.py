"""Fig. 3, the porewater Eh-pH stability diagram, from artifacts/aqueous/ehph_fields.npz.

Fields are drawn from the per-phase scores (SI per Fe atom) that build_ehph_diagram.py stores, not
from the cell labels: each scores array is interpolated onto a fine grid and a phase boundary is
where two scores cross, so the lines follow the chemistry instead of the calculation grid. The
predominant phase on the fine grid is checked against the coarse labels before anything is drawn.
Fields are coloured and labelled in place; greigite keeps its hatch and
bold outline, with its range under the reported enthalpy uncertainty (2 SD of the mean) shaded.
No title or caption is baked onto the image. Writes artifacts/figures/Figure_3.png.
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.interpolate import RegularGridInterpolator
import pathlib

REPO = pathlib.Path(__file__).resolve().parents[1]
ART = REPO / "artifacts"
(ART / "figures").mkdir(parents=True, exist_ok=True)
d = np.load(ART / "aqueous" / "ehph_fields.npz", allow_pickle=True)
pH, pe, labels = d["pH"], d["pe"], list(d["labels"])
cands = list(d["candidates"])
VPE = 0.059159                      # V per pe unit at 25 C
Eh = pe * VPE

COLOR = {"Hematite": "#f4b394", "Pyrite": "#efc85a", "Greigite": "#7fb0ea",
         "Magnetite": "#7fd3b3", "Pyrrhotite": "#c9a0dc", "Goethite": "#e3a36f",
         "Siderite": "#b8c96a", "Fe(OH)3(a)": "#f7d2c0"}
BAND = "#cfe0f7"

# fine grid: 8x in each direction
UP = 8
fpH = np.linspace(pH.min(), pH.max(), (len(pH) - 1) * UP + 1)
fpe = np.linspace(pe.min(), pe.max(), (len(pe) - 1) * UP + 1)
PE, PH = np.meshgrid(fpe, fpH, indexing="ij")
fEh = fpe * VPE
pts = np.column_stack([PE.ravel(), PH.ravel()])

def fine_scores(variant):
    """Interpolated SI/n_Fe per candidate on the fine grid; undefined SI never wins."""
    s = np.where(np.isnan(d["score_" + variant]), -1e3, d["score_" + variant])
    return {m: RegularGridInterpolator((pe, pH), s[k])(pts).reshape(PE.shape)
            for k, m in enumerate(cands)}

def field_fn(sc, phase):
    """>0 inside the phase field: its score minus the best competitor (aqueous = 0)."""
    others = np.maximum.reduce([sc[m] for m in cands if m != phase] + [np.zeros(PE.shape)])
    return sc[phase] - others

def coarse_check(sc, variant):
    """The fine-grid winner at the coarse nodes must equal the stored labels."""
    stack = np.stack([sc[m] for m in cands])
    best = np.argmax(stack, axis=0)
    win = np.where(stack.max(axis=0) >= 0, best + 1, 0)       # +1: labels index 0 is aqueous
    node = win[::UP, ::UP]
    want = np.vectorize(lambda v: {labels.index(m): k + 1 for k, m in enumerate(cands)}.get(v, 0))(d[variant])
    bad = int((node != want).sum())
    assert bad == 0, f"{variant}: {bad} coarse nodes disagree with the stored labels"

S = {v: fine_scores(v) for v in ("nominal", "hi", "lo")}
for v in S:
    coarse_check(S[v], v)

upper = 1.229 - VPE * PH
lower = -VPE * PH
EH = PE * VPE
inside = (EH <= upper) & (EH >= lower)

plt.rcParams.update({"font.size": 14, "axes.linewidth": 1.1, "hatch.linewidth": 0.8})
fig, ax = plt.subplots(figsize=(7.8, 6.4))
ax.set_facecolor("white")

present = []
for m in cands:
    f = np.where(inside, field_fn(S["nominal"], m), -1)
    if (f > 0).sum() > 50 * UP * UP:          # ignore slivers of a few coarse cells
        present.append(m)
        ax.contourf(fpH, fEh, f, levels=[0, 1e9], colors=[COLOR[m]],
                    zorder=1.6 if m == "Greigite" else 1)   # greigite sits above the band

# greigite uncertainty band: inside the most-stable field but outside the least-stable one
g_lo = np.where(inside, field_fn(S["lo"], "Greigite"), -1)
g_hi = np.where(inside, field_fn(S["hi"], "Greigite"), -1)
g_nom = np.where(inside, field_fn(S["nominal"], "Greigite"), -1)
band = np.where((g_lo > 0) & (g_hi <= 0), 1.0, 0.0)
ax.contourf(fpH, fEh, band, levels=[0.5, 1.5], colors=[BAND], zorder=1.5)
ax.contourf(fpH, fEh, g_nom, levels=[0, 1e9], colors="none", hatches=["//"], zorder=2)

# phase boundaries between the drawn fields (thin), then the greigite curves
for m in present:
    if m != "Greigite":
        ax.contour(fpH, fEh, np.where(inside, field_fn(S["nominal"], m), -1), levels=[0],
                   colors="0.2", linewidths=0.9, zorder=2.5)
ax.contour(fpH, fEh, g_nom, levels=[0], colors="k", linewidths=2.6, zorder=3)
for g in (g_hi, g_lo):
    ax.contour(fpH, fEh, g, levels=[0], colors="#1c3f73", linewidths=1.4,
               linestyles=[(0, (5, 2))], zorder=3)

# outside the stability limits of water: stippled, no mineral fields
out = (~inside).astype(float)
ax.contourf(fpH, fEh, out, levels=[0.5, 1.5], colors=["#efefef"], zorder=0.5)
ax.contourf(fpH, fEh, out, levels=[0.5, 1.5], colors="none", hatches=["...."], zorder=0.6)
ax.plot(pH, 1.229 - VPE * pH, "k-", lw=1.3, zorder=4)
ax.plot(pH, -VPE * pH, "k-", lw=1.3, zorder=4)
LBOX = dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85)
ax.text(6.4, 0.94, "O$_2$ degassing", fontsize=12, style="italic", ha="center", bbox=LBOX, zorder=5)
ax.text(4.0, -0.58, "H$_2$ degassing", fontsize=12, style="italic", ha="center", bbox=LBOX, zorder=5)

# direct field labels at the deepest point of each field (distance to its edge, on the fine grid)
from scipy.ndimage import distance_transform_edt
LABEL_TXT = {"Greigite": "Greigite", "Pyrite": "Pyrite", "Hematite": "Hematite",
             "Magnetite": "Magnetite", "Pyrrhotite": "Pyrrhotite"}
# narrow fields get a hand-placed label (checked on the rendered figure); the rest go deepest-in-field
FIXED = {"Pyrite": (2.6, 0.085)}
for m in present:
    if m in FIXED:
        ax.text(*FIXED[m], LABEL_TXT.get(m, m), fontsize=13, ha="center", va="center", weight="bold", zorder=5)
        continue
    mask = np.where(inside, field_fn(S["nominal"], m), -1) > 0
    if m == "Greigite":
        mask &= g_hi > 0        # keep the label off the uncertainty curves
    mask[[0, -1], :] = False    # the plot frame counts as an edge
    mask[:, [0, -1]] = False
    dist = distance_transform_edt(mask, sampling=(6.4 / (len(fEh) - 1), 7.8 / (len(fpH) - 1)))
    j, i = np.unravel_index(np.argmax(dist), dist.shape)
    ax.text(fpH[i], fEh[j], LABEL_TXT.get(m, m), fontsize=13, ha="center", va="center",
            weight="bold", zorder=5,
            bbox=LBOX if m == "Greigite" else None)
ax.text(2.8, 0.29, "aqueous\nFe", fontsize=13, ha="center", va="center", weight="bold", zorder=5)

ax.set_xlabel("pH", fontsize=16)
ax.set_ylabel("Eh (V)", fontsize=16)
ax.tick_params(labelsize=14)
handles = [plt.Line2D([0], [0], color="k", lw=2.6, label="greigite, mean"),
           Patch(facecolor=BAND, edgecolor="#1c3f73", ls="--", label="greigite, ±2 SD")]
ax.legend(handles=handles, fontsize=13, loc="upper right", framealpha=0.95)
ax.set_xlim(pH.min(), pH.max())
ax.set_ylim(Eh.min(), Eh.max())
fig.tight_layout()
out_png = ART / "figures" / "Figure_3.png"
fig.savefig(out_png, dpi=320)
print("saved", out_png, "| fields drawn:", present)
