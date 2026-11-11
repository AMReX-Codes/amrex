"""Draw extend_domain_face.png for the EB chapter of the AMReX User's Guide.

Usage: python extend_domain_face.py [extend_domain_face.png]
Requires numpy and matplotlib. Uses the Noto Sans Mono and Lato fonts when
they are installed.

A smooth body crosses the lower domain face. The embedded boundary is drawn
as AMReX represents it: one straight segment in each cut cell, through the
points where the level set changes sign on the cell edges. With
extend_domain_face = true, nodes outside the domain take the level set value
on the domain face, so the body outside the domain is an extrusion of the
face. With extend_domain_face = false, the actual surface is used.
"""
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle

MONO = ["Noto Sans Mono", "DejaVu Sans Mono"]
SANS = ["Lato", "DejaVu Sans"]

NX, NY, NOUT = 10, 6, 3         # domain cells in x and y, cells shown below the face

BODY = "#a6a6a6"                # body inside the domain, as in the other EB figures
OUTBODY = "#c8c8c8"             # body outside the domain
EDGE = "#1a1a1a"                # embedded boundary
OUTSIDE = "#f2f2f2"             # cells outside the domain
GRID = "#bdbdbd"
DOMAIN = "#303030"
TEXT = "#303030"


def radius(th):
    return 2.75 * (1 + 0.07 * np.cos(3 * th) + 0.05 * np.sin(2 * th))


def levelset(x, y, extend):
    """Negative inside the body."""
    if extend:
        y = np.maximum(y, 0.0)
    dx, dy = x - 5.0, (y - 1.3) / 1.1
    return np.hypot(dx, dy) - radius(np.arctan2(dy, dx))


def cut_cell(x0, y0, phi):
    """Body part of the cell with lower left corner (x0, y0) and the
    embedded boundary segment in it (marching squares on the node values)."""
    corners = [(x0, y0), (x0 + 1, y0), (x0 + 1, y0 + 1), (x0, y0 + 1)]
    vals = [phi(*c) for c in corners]
    poly, cuts = [], []
    for m in range(4):
        a, b = corners[m], corners[(m + 1) % 4]
        va, vb = vals[m], vals[(m + 1) % 4]
        if va < 0:
            poly.append(a)
        if (va < 0) != (vb < 0):
            t = va / (va - vb)
            p = (a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]))
            poly.append(p)
            cuts.append(p)
    return poly, cuts


def draw(ax, extend):
    ax.add_patch(Rectangle((0, -NOUT), NX, NOUT, facecolor=OUTSIDE, edgecolor="none",
                           zorder=0))

    phi = lambda x, y: levelset(x, y, extend)
    for i in range(NX):
        for j in range(-NOUT, NY):
            poly, cuts = cut_cell(i, j, phi)
            if poly:
                color = BODY if j >= 0 else OUTBODY
                # an edge in the fill color hides the seams between cells
                ax.add_patch(Polygon(poly, closed=True, facecolor=color, edgecolor=color,
                                     lw=0.6, zorder=1))
            if len(cuts) == 2:
                ax.plot(*np.array(cuts).T, color=EDGE, lw=1.8, zorder=3,
                        solid_capstyle="round")

    for x in range(NX + 1):
        ax.plot([x, x], [-NOUT, NY], color=GRID, lw=0.6, zorder=2)
    for y in range(-NOUT, NY + 1):
        ax.plot([0, NX], [y, y], color=GRID, lw=0.6, zorder=2)

    if extend:
        # the actual surface outside the domain, which is not used
        th = np.linspace(0, 2 * np.pi, 721)
        x = 5.0 + radius(th) * np.cos(th)
        y = 1.3 + 1.1 * radius(th) * np.sin(th)
        y = np.where(y < 0, y, np.nan)
        ax.plot(x, y, color=EDGE, lw=1.0, ls=(0, (3, 2.5)), alpha=0.6, zorder=3)

    ax.add_patch(Rectangle((0, 0), NX, NY, facecolor="none", edgecolor=DOMAIN, lw=1.8,
                           zorder=5))
    ax.text(NX - 1.2, -NOUT / 2, "ghost cells", ha="center", va="center", fontsize=12,
            family=SANS, color=TEXT, zorder=6,
            bbox=dict(boxstyle="round,pad=0.25", facecolor=OUTSIDE, edgecolor="none"))

    ax.set_title(f"extend_domain_face = {'true' if extend else 'false'}",
                 fontsize=13, family=MONO, color=TEXT, pad=8)
    ax.set_xlim(-0.05, NX + 0.05)
    ax.set_ylim(-NOUT - 0.05, NY + 0.05)
    ax.set_aspect("equal")
    ax.axis("off")


fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.8))
for ax, extend in zip(axes, (True, False)):
    draw(ax, extend)
fig.tight_layout(w_pad=3.0)

out = sys.argv[1] if len(sys.argv) > 1 else "extend_domain_face.png"
fig.savefig(out, dpi=150, facecolor="white", bbox_inches="tight", pad_inches=0.08)
