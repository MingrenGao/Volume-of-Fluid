# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# =============================================================================
# GLOBAL FONT SETTINGS (Times-like, works on Linux)
# =============================================================================
plt.rcParams.update({
    "font.family": "STIXGeneral",
    "mathtext.fontset": "stix",
    "font.size": 10,
})

# =============================================================================
# CONFIG (put all user controls here)
# =============================================================================
CONFIG = {
    # input files
    "KAPPA_FILE": "kappa_field.dat",
    "SEG_FILE": "interface_segments.dat",

    # output
    "OUT_PNG": "kappa_gridcolored_with_redline.png",
    "DPI": 400,

    # kappa colormap on grid
    "CMAP": "RdBu_r",
    "VMIN": -1,
    "VMAX":  1,

    # grid line appearance
    "DRAW_GRID_LINES": True,
    "GRID_LW": 0.6,
    "GRID_ALPHA": 1.0,

    # analytic/interface line (red)
    "LINE_COLOR": "r",
    "LINE_LW": 1.5,

    # border
    "SPINE_LW": 1.5,

    # axes
    "HIDE_TICKS": True,

    # layout
    "FIGSIZE": (5, 3.8),
    "TIGHT_BBOX": True,
    "PAD_INCHES": 0.02,

    # ---------------- colorbar style ----------------
    "CBAR_FRACTION": 0.045,
    "CBAR_PAD": 0.04,
    "CBAR_OUTLINE_LW": 2.0,
    "CBAR_SHOW_TICKS": False,
    "CBAR_ADD_KMIN_KMAX": True,
    "CBAR_FONT_SIZE": 14,
    "CBAR_WEIGHT": "bold",
}

# =============================================================================
# load kappa field
# =============================================================================
def load_kappa_field(fname):
    data = np.loadtxt(Path(fname), comments="#")
    i = data[:, 0].astype(int)
    j = data[:, 1].astype(int)
    x = data[:, 2]
    y = data[:, 3]
    k = data[:, 4]

    nx = i.max() + 1
    ny = j.max() + 1

    xg = np.full((ny, nx), np.nan)
    yg = np.full((ny, nx), np.nan)
    kg = np.full((ny, nx), np.nan)

    xg[j, i] = x
    yg[j, i] = y
    kg[j, i] = k
    return xg, yg, kg


def cell_edges_from_centers(xg, yg):
    dx = xg[0, 1] - xg[0, 0]
    dy = yg[1, 0] - yg[0, 0]
    x0 = xg[0, 0] - 0.5 * dx
    y0 = yg[0, 0] - 0.5 * dy
    nx = xg.shape[1]
    ny = xg.shape[0]
    xe = x0 + dx * np.arange(nx + 1)
    ye = y0 + dy * np.arange(ny + 1)
    return xe, ye


# =============================================================================
# load interface segments
# =============================================================================
def load_interface_segments_xy(fname):
    data = np.loadtxt(Path(fname), comments="#")
    x1 = data[:, 2]; y1 = data[:, 3]
    x2 = data[:, 4]; y2 = data[:, 5]
    return x1, y1, x2, y2


# =============================================================================
# draw structured grid lines
# =============================================================================
def draw_structured_mesh_lines(ax, xe, ye, lw=0.6, alpha=1.0):
    for x in xe:
        ax.plot([x, x], [ye[0], ye[-1]],
                color="k", lw=lw, alpha=alpha, zorder=3)
    for y in ye:
        ax.plot([xe[0], xe[-1]], [y, y],
                color="k", lw=lw, alpha=alpha, zorder=3)


# =============================================================================
# main
# =============================================================================
xg, yg, kg = load_kappa_field(CONFIG["KAPPA_FILE"])
xe, ye = cell_edges_from_centers(xg, yg)
x1, y1, x2, y2 = load_interface_segments_xy(CONFIG["SEG_FILE"])

fig, ax = plt.subplots(figsize=CONFIG["FIGSIZE"])

ax.set_xlim(xe[0], xe[-1])
ax.set_ylim(ye[0], ye[-1])
ax.set_aspect("equal")

if CONFIG["HIDE_TICKS"]:
    ax.set_xticks([])
    ax.set_yticks([])

for sp in ax.spines.values():
    sp.set_linewidth(CONFIG["SPINE_LW"])

im = ax.pcolormesh(
    xe, ye, kg,
    shading="flat",
    cmap=CONFIG["CMAP"],
    vmin=CONFIG["VMIN"],
    vmax=CONFIG["VMAX"],
    zorder=1,
)

if CONFIG["DRAW_GRID_LINES"]:
    draw_structured_mesh_lines(
        ax, xe, ye,
        lw=CONFIG["GRID_LW"],
        alpha=CONFIG["GRID_ALPHA"],
    )

for a, b, c, d in zip(x1, y1, x2, y2):
    ax.plot([a, c], [b, d],
            color=CONFIG["LINE_COLOR"],
            lw=CONFIG["LINE_LW"],
            zorder=4)

# =============================================================================
# COLORBAR (κmax / κmin on the right)
# =============================================================================
cbar = fig.colorbar(
    im, ax=ax,
    fraction=CONFIG["CBAR_FRACTION"],
    pad=CONFIG["CBAR_PAD"],
)

FS = CONFIG["CBAR_FONT_SIZE"]
WT = CONFIG["CBAR_WEIGHT"]

cbar.set_label(
    r"Curvature $\kappa$",
    fontsize=FS,
    weight=WT,
    labelpad=10,
)

if not CONFIG["CBAR_SHOW_TICKS"]:
    cbar.set_ticks([])
    cbar.ax.tick_params(length=0)

cbar.outline.set_linewidth(CONFIG["CBAR_OUTLINE_LW"])
for sp in cbar.ax.spines.values():
    sp.set_linewidth(CONFIG["CBAR_OUTLINE_LW"])

if CONFIG["CBAR_ADD_KMIN_KMAX"]:
    cbar.ax.text(
        1.5, 1.00, r"$\kappa_{\max}$",
        transform=cbar.ax.transAxes,
        ha="left", va="center",
        fontsize=FS,
        weight=WT,
    )
    cbar.ax.text(
        1.5, 0.00, r"$\kappa_{\min}$",
        transform=cbar.ax.transAxes,
        ha="left", va="center",
        fontsize=FS,
        weight=WT,
    )

# =============================================================================
# save
# =============================================================================
fig.tight_layout()

save_kwargs = dict(dpi=CONFIG["DPI"])
if CONFIG["TIGHT_BBOX"]:
    save_kwargs.update(
        dict(bbox_inches="tight",
             pad_inches=CONFIG["PAD_INCHES"])
    )

fig.savefig(CONFIG["OUT_PNG"], **save_kwargs)
plt.close(fig)

print(f"Saved: {Path(CONFIG['OUT_PNG']).resolve()}")
