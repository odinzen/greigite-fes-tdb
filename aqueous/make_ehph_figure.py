"""Fig. 3 — render the porewater Eh-pH stability diagram from artifacts/aqueous/ehph_fields.npz.

Greyscale house style; greigite is the hero phase (mid grey + hatch + bold outline) with its +/-1
sigma range shaded. No title or caption is baked onto the image (it belongs in the manuscript
caption); axis labels only. Writes artifacts/figures/Figure_3.png.
"""
import numpy as np, matplotlib, pathlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.ndimage import median_filter, binary_opening, binary_closing

REPO = pathlib.Path(__file__).resolve().parents[1]
ART = REPO / "artifacts"
(ART / "figures").mkdir(parents=True, exist_ok=True)


def denoise_labels(a, size=3):
    return median_filter(a.astype(int), size=size)


def denoise_mask(m, size=2):
    return binary_closing(binary_opening(m, iterations=1), iterations=1)


d = np.load(ART / "aqueous" / "ehph_fields.npz", allow_pickle=True)
pH, pe, labels = d["pH"], d["pe"], list(d["labels"])
Eh = pe * 0.059159  # V at 25 C

# grey levels (0=black..1=white); greigite gets a mid grey + hatch as the hero phase
GREY = {"aqueous": 1.0, "Hematite": 0.86, "Goethite": 0.80, "Magnetite": 0.45,
        "Greigite": 0.62, "Pyrite": 0.30, "Pyrrhotite": 0.15, "Siderite": 0.72,
        "Fe(OH)3(a)": 0.92, "error": 1.0}
order = ["aqueous", "Hematite", "Goethite", "Fe(OH)3(a)", "Siderite", "Magnetite",
         "Greigite", "Pyrite", "Pyrrhotite", "error"]
cmap = matplotlib.colors.ListedColormap([str(GREY[n]) for n in order])
idx_remap = {labels.index(n): i for i, n in enumerate(order) if n in labels}

plt.rcParams.update({"font.size": 12, "axes.linewidth": 1.1})
fig, ax = plt.subplots(figsize=(7.8, 6.4))
grid = d["nominal"].astype(float)

# water-stability window: outside it there is no liquid water (H2 or O2 degassing) -> mask
Eh2d, pH2d = np.meshgrid(Eh, pH, indexing="ij")
upper = 1.229 - 0.059159 * pH2d
lower = -0.059159 * pH2d
outside = (Eh2d > upper) | (Eh2d < lower)

grid = denoise_labels(grid)                          # despeckle the field map
remapped = np.vectorize(lambda v: idx_remap.get(int(v), 0))(grid).astype(float)
remapped_in = np.where(outside, np.nan, remapped)   # mineral fields only inside the window
ax.pcolormesh(pH, Eh, remapped_in, cmap=cmap, vmin=0, vmax=len(order) - 1, shading="nearest")

# grey-hatch the degassing regions (no stable water)
ax.contourf(pH, Eh, outside.astype(float), levels=[0.5, 1.5], colors="none", hatches=["...."])
ax.pcolormesh(pH, Eh, np.where(outside, 1.0, np.nan),
              cmap=matplotlib.colors.ListedColormap(["0.93"]), shading="nearest", alpha=0.55)

# greigite: nominal field (bold outline + hatch) and a single +/-1sigma uncertainty band
gi = labels.index("Greigite")
g_nom = denoise_mask((d["nominal"] == gi) & ~outside)
g_hi = denoise_mask((d["hi"] == gi) & ~outside)     # +1s least stable = smallest field
g_lo = denoise_mask((d["lo"] == gi) & ~outside)     # -1s most stable = largest field
band = (g_lo & ~g_hi)                                # where greigite presence is uncertain
ax.pcolormesh(pH, Eh, np.where(band, 1.0, np.nan),
              cmap=matplotlib.colors.ListedColormap(["0.72"]), shading="nearest", alpha=0.7, zorder=1.5)
ax.contourf(pH, Eh, g_nom.astype(float), levels=[0.5, 1.5], colors="none", hatches=["///"], zorder=2)
ax.contour(pH, Eh, g_nom.astype(float), levels=[0.5], colors="k", linewidths=2.6, zorder=3)     # nominal
ax.contour(pH, Eh, g_hi.astype(float), levels=[0.5], colors="0.25", linewidths=1.4, linestyles=[(0, (5, 2))], zorder=3)
ax.contour(pH, Eh, g_lo.astype(float), levels=[0.5], colors="0.25", linewidths=1.4, linestyles=[(0, (5, 2))], zorder=3)

# water stability limits (bold) + labels
ax.plot(pH, 1.229 - 0.059159 * pH, "k-", lw=1.3)
ax.plot(pH, -0.059159 * pH, "k-", lw=1.3)
# O2 label sits left of the legend (upper right) so the two never overlap
ax.text(5.0, 1.229 - 0.059159 * 5.0 + 0.02, "O$_2$ degassing", fontsize=7, style="italic", ha="left")
ax.text(11.3, -0.059159 * 11.3 - 0.06, "H$_2$ degassing", fontsize=7, style="italic", ha="right")

ax.set_xlabel("pH", fontsize=14)
ax.set_ylabel("Eh (V)", fontsize=14)
ax.tick_params(labelsize=12)
present = [n for n in order if n in labels and n != "error"
          and ((grid == labels.index(n)) & ~outside).any()]
handles = [Patch(facecolor=str(GREY[n]), edgecolor="k",
                 hatch=("////" if n == "Greigite" else None), label=n) for n in present]
handles += [plt.Line2D([0], [0], color="k", lw=2.6, label="greigite (nominal)"),
            Patch(facecolor="0.72", alpha=0.7, edgecolor="0.25", ls="--",
                  label="greigite ±1σ range")]
ax.legend(handles=handles, fontsize=9.5, loc="upper right", framealpha=0.95)
ax.set_xlim(pH.min(), pH.max())
ax.set_ylim(Eh.min(), Eh.max())
fig.tight_layout()
out = ART / "figures" / "Figure_3.png"
fig.savefig(out, dpi=320)
print("saved", out)
