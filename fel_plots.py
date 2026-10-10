"""
fel_plots.py -- styled plotting for the PCA / FEL pipeline + style previews.
============================================================================

All five pipeline figures are drawn here, driven by a style dict from
fel_styles.py. The style previews are made by THESE SAME functions on a
synthetic landscape, so what you see in the picker is exactly what the
real run produces (only the data differ).

Pure numpy + matplotlib (scipy optional, for smoothing) -- no Schrodinger
imports, so it also runs in plain python3.

Generate the preview images used by the PyMOL style picker:
    python3 fel_plots.py --previews fel_style_previews
    $SCHRODINGER/run python3 fel_plots.py --previews fel_style_previews

Make a single contact sheet of all 20 styles:
    python3 fel_plots.py --gallery fel_style_gallery.png

License: MIT (see LICENSE file).
"""

import os
import sys
import json
import logging
import argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")  # non-interactive backend, safe for batch/cluster runs
import matplotlib.pyplot as plt
import matplotlib.patheffects as patheffects
from matplotlib.colors import LinearSegmentedColormap, to_rgba

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fel_styles import get_style, list_styles, STYLES  # noqa: E402

# Missing optional fonts (e.g. Times New Roman on Linux) just fall back to
# the next one in the list -- don't spam the log about it.
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

KB = 0.0019872041          # Boltzmann constant, kcal/(mol*K)


# ----------------------------------------------------------------------
# Style plumbing
# ----------------------------------------------------------------------

def get_cmap(spec, name="custom"):
    """A matplotlib colormap name, or a list of colours (low -> high)."""
    if isinstance(spec, (list, tuple)):
        return LinearSegmentedColormap.from_list(name, list(spec), N=512)
    return plt.get_cmap(spec)


def style_rc(style):
    """rcParams for a style -- applied with plt.rc_context() around each figure."""
    fs = style["font_scale"]
    fg = style["fg"]
    rc = {
        "figure.dpi": 150,
        "font.size": 13 * fs,
        "font.family": style["font_family"],
        "font.serif": ["Times New Roman", "Times", "Nimbus Roman", "Liberation Serif",
                       "DejaVu Serif"],
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
        "font.monospace": ["DejaVu Sans Mono", "Liberation Mono", "Courier New"],
        "mathtext.fontset": "dejavuserif" if style["font_family"] == "serif" else "dejavusans",
        "axes.titlesize": 16 * fs,
        "axes.titleweight": style["title_weight"],
        "axes.labelsize": 14 * fs,
        "axes.labelweight": style["label_weight"],
        "axes.linewidth": 1.2,
        "axes.facecolor": style["bg"],
        "axes.edgecolor": fg,
        "axes.labelcolor": fg,
        "axes.titlecolor": fg,
        "text.color": fg,
        "xtick.color": fg,
        "ytick.color": fg,
        "xtick.labelsize": 11 * fs,
        "ytick.labelsize": 11 * fs,
        "grid.color": fg,
        "figure.facecolor": style["bg"],
        "savefig.facecolor": style["bg"],
        "savefig.edgecolor": style["bg"],
        "savefig.bbox": "tight",
        "legend.facecolor": style["bg"],
        "legend.edgecolor": fg,
        "legend.labelcolor": fg,
    }
    if style["open_spines"]:
        rc["axes.spines.top"] = False
        rc["axes.spines.right"] = False
    return rc


def _style_axes(ax, style):
    ax.set_facecolor(style["bg"])
    if style["grid"]:
        ax.grid(True, alpha=style["grid_alpha"], linewidth=0.6)
        ax.set_axisbelow(True)


def _title(ax, style, text, pad=12):
    if style["show_title"]:
        ax.set_title(text, pad=pad)


def _style_colorbar(cbar, style, label, vrange=None):
    if vrange is not None:
        # contourf colorbars default to ticks on contour levels (1.766, 1.589...);
        # use round numbers instead.
        from matplotlib.ticker import MaxNLocator
        lo, hi = vrange
        ticks = [t for t in MaxNLocator(6).tick_values(lo, hi) if lo - 1e-9 <= t <= hi + 1e-9]
        cbar.set_ticks(ticks)
    cbar.set_label(label, fontsize=13 * style["font_scale"], fontweight=style["label_weight"],
                   color=style["fg"])
    cbar.ax.tick_params(labelsize=10 * style["font_scale"], colors=style["fg"])
    cbar.outline.set_edgecolor(style["fg"])
    cbar.outline.set_linewidth(0.8)


def _save(fig, out_png, dpi, what):
    fig.savefig(out_png, dpi=dpi)
    plt.close(fig)
    print(f"Saved {what} to {out_png} ({dpi} dpi)")


def pc_axis_labels(variance_explained):
    """'PC1 (25.6%)' style labels when variance is known, else 'PC1'/'PC2'."""
    if variance_explained is None:
        return "PC1", "PC2"
    pc1_pct, pc2_pct = variance_explained
    return f"PC1 ({pc1_pct:.1f}%)", f"PC2 ({pc2_pct:.1f}%)"


# ----------------------------------------------------------------------
# FEL computation (numpy only)
# ----------------------------------------------------------------------

def compute_fel(pc1, pc2, n_bins=60, temperature=300.0):
    """
    2D free energy landscape from PC1/PC2 projections:
        FEL(x, y) = -kB*T * ln( P(x, y) / P_max )
    Returns bin centres (x, y) and the FEL grid in kcal/mol (unsampled = NaN).
    """
    hist, xedges, yedges = np.histogram2d(pc1, pc2, bins=n_bins)
    prob = hist / hist.sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        fel = -KB * temperature * np.log(prob / prob.max())
    fel[hist == 0] = np.nan
    x_centers = 0.5 * (xedges[:-1] + xedges[1:])
    y_centers = 0.5 * (yedges[:-1] + yedges[1:])
    return x_centers, y_centers, fel.T   # transpose for correct x/y orientation


def find_local_minima(fel, min_distance=8, max_minima=5, energy_cutoff=1.5):
    """
    Distinct local minima (basins) of a 2D FEL grid, lowest first, as
    (energy, row_idx, col_idx). Only minima within `energy_cutoff` kcal/mol
    of the global minimum are kept, and minima closer than `min_distance`
    grid points to an accepted one are discarded.
    """
    ny, nx = fel.shape
    finite = np.isfinite(fel)
    global_min = np.nanmin(fel)

    candidates = []
    for i in range(ny):
        for j in range(nx):
            if not finite[i, j]:
                continue
            val = fel[i, j]
            if val > global_min + energy_cutoff:
                continue
            neighborhood = fel[max(0, i - 1):min(ny, i + 2), max(0, j - 1):min(nx, j + 2)]
            if val <= np.nanmin(neighborhood):
                candidates.append((val, i, j))

    candidates.sort(key=lambda t: t[0])
    selected = []
    for val, i, j in candidates:
        if any(abs(si - i) < min_distance and abs(sj - j) < min_distance
               for _, si, sj in selected):
            continue
        selected.append((val, i, j))
        if len(selected) >= max_minima:
            break
    return selected


def _simple_box_blur(a, iterations=2):
    """Minimal separable blur without scipy, used only if scipy is unavailable."""
    kernel = np.array([0.25, 0.5, 0.25])
    out = a.copy()
    for _ in range(iterations):
        out = np.apply_along_axis(lambda m: np.convolve(m, kernel, mode="same"), axis=0, arr=out)
        out = np.apply_along_axis(lambda m: np.convolve(m, kernel, mode="same"), axis=1, arr=out)
    return out


def smooth_and_upsample_fel(x_centers, y_centers, fel, upsample_factor=6, sigma=1.2):
    """
    Smooth the FEL grid and upsample it onto a finer mesh (blocky histogram ->
    continuous surface). Unsampled (NaN) regions stay blank. sigma=0 and
    upsample_factor=1 return the raw grid untouched.
    """
    if (not sigma or sigma <= 0) and (not upsample_factor or upsample_factor <= 1):
        return x_centers, y_centers, fel

    nan_mask = ~np.isfinite(fel)
    fel_filled = np.where(nan_mask, np.nanmax(fel), fel)
    upsample_factor = max(1, int(upsample_factor or 1))

    try:
        from scipy.ndimage import gaussian_filter, zoom
        if sigma and sigma > 0:
            fel_smooth = gaussian_filter(fel_filled, sigma=sigma)
            mask_smooth = gaussian_filter(nan_mask.astype(float), sigma=sigma)
        else:
            fel_smooth, mask_smooth = fel_filled, nan_mask.astype(float)
        if upsample_factor > 1:
            fel_fine = zoom(fel_smooth, upsample_factor, order=3)
            mask_fine = zoom(mask_smooth, upsample_factor, order=1)
        else:
            fel_fine, mask_fine = fel_smooth, mask_smooth
        x_fine = np.linspace(x_centers[0], x_centers[-1], fel_fine.shape[1])
        y_fine = np.linspace(y_centers[0], y_centers[-1], fel_fine.shape[0])
    except ImportError:
        fel_fine = _simple_box_blur(fel_filled, iterations=3)
        mask_fine = _simple_box_blur(nan_mask.astype(float), iterations=3)
        x_fine, y_fine = x_centers, y_centers

    fel_fine = np.where(mask_fine > 0.5, np.nan, fel_fine)
    return x_fine, y_fine, fel_fine


def _prepared_fel(x_centers, y_centers, fel, style):
    return smooth_and_upsample_fel(x_centers, y_centers, fel,
                                   upsample_factor=style["upsample"],
                                   sigma=style["smooth_sigma"])


def _levels(fel):
    finite_vals = fel[np.isfinite(fel)]
    if finite_vals.size == 0:
        raise ValueError(
            "FEL grid is entirely empty -- check your PC1/PC2 projections "
            "and n_bins (too many bins for too few frames can cause this)."
        )
    vmin, vmax = finite_vals.min(), finite_vals.max()
    if vmin == vmax:
        vmax = vmin + 1e-6   # guard against a perfectly flat landscape
    return vmin, vmax


# ----------------------------------------------------------------------
# Drawing primitives (draw_* draw onto a given axes; plot_* make + save a figure)
# ----------------------------------------------------------------------

def draw_fel_2d(fig, ax, x_centers, y_centers, fel, style, minima=None,
                variance_explained=None, title="Free Energy Landscape (PC1 vs PC2)",
                colorbar=True, labels=True):
    """Filled 2D FEL (+ optional iso-lines, minima markers, colorbar)."""
    x, y, f = _prepared_fel(x_centers, y_centers, fel, style)
    vmin, vmax = _levels(f)
    cmap = get_cmap(style["fel_cmap"], "fel")
    cmap.set_bad(style["bg"])

    if style["fel_render"] == "pcolormesh":
        mappable = ax.pcolormesh(x, y, np.ma.masked_invalid(f), cmap=cmap,
                                 vmin=vmin, vmax=vmax, shading="nearest")
    else:
        levels = np.linspace(vmin, vmax, max(3, int(style["n_levels"])))
        mappable = ax.contourf(x, y, f, levels=levels, cmap=cmap, extend="neither")
        if style["contour_lines"]:
            every = max(1, int(style["contour_every"]))
            ax.contour(x, y, f, levels=levels[::every], colors=style["contour_color"],
                       linewidths=style["contour_lw"], alpha=style["contour_alpha"])

    if colorbar:
        cbar = fig.colorbar(mappable, ax=ax, pad=0.02)
        _style_colorbar(cbar, style, "Free energy (kcal/mol)", vrange=(vmin, vmax))

    if minima:
        for m in minima:
            mx, my = m["pc1"], m["pc2"]
            ax.plot(mx, my, marker=style["minima_marker"], markersize=11,
                    markeredgecolor=style["minima_edge"], markerfacecolor=style["minima_face"],
                    markeredgewidth=2.0, zorder=5, linestyle="none")
            txt = ax.annotate(str(m["rank"]), xy=(mx, my), xytext=(6, 6),
                              textcoords="offset points", fontsize=13 * style["font_scale"],
                              fontweight="bold", color=style["minima_text"], zorder=6)
            txt.set_path_effects([patheffects.withStroke(linewidth=3,
                                                         foreground=style["minima_halo"])])

    if labels:
        xl, yl = pc_axis_labels(variance_explained)
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        _title(ax, style, title)
    _style_axes(ax, style)
    if style["fel_render"] == "pcolormesh":
        ax.grid(False)   # grid lines over raw pixels look like bin edges
    return mappable


def draw_fel_3d(fig, ax, x_centers, y_centers, fel, style, minima=None,
                variance_explained=None, colorbar=True, labels=True):
    """3D FEL surface with optional projected contour floor and minima markers."""
    from mpl_toolkits.mplot3d import proj3d

    # 3D always uses a smoothed surface -- raw bins look broken in 3D.
    sigma = style["smooth_sigma"] if style["smooth_sigma"] > 0 else 1.0
    up = style["upsample"] if style["upsample"] > 1 else 4
    x, y, f = smooth_and_upsample_fel(x_centers, y_centers, fel, upsample_factor=up, sigma=sigma)
    X, Y = np.meshgrid(x, y)
    vmin, vmax = _levels(f)
    f_filled = np.where(np.isfinite(f), f, vmax)
    cmap = get_cmap(style["fel_cmap"], "fel3d")

    edge = style["surface_edge"]
    surf = ax.plot_surface(
        X, Y, f_filled, cmap=cmap, rstride=3 if edge else 2, cstride=3 if edge else 2,
        linewidth=0.12 if edge else 0,
        edgecolor=to_rgba(edge, 0.35) if edge else "none",
        antialiased=True, alpha=style["surface_alpha"], vmin=vmin, vmax=vmax,
    )

    floor_z = vmin - 0.15 * (vmax - vmin)
    if style["floor"]:
        ax.contourf(X, Y, f_filled, zdir="z", offset=floor_z,
                    levels=np.linspace(vmin, vmax, max(3, int(style["n_levels"]))), cmap=cmap)

    ax.set_zlim(floor_z, vmax)
    ax.set_facecolor(style["bg"])
    if labels:
        xl, yl = pc_axis_labels(variance_explained)
        ax.set_xlabel(xl, labelpad=10)
        ax.set_ylabel(yl, labelpad=10)
        ax.set_zlabel("Free energy (kcal/mol)", labelpad=10)
        if style["show_title"]:
            ax.set_title("3D Free Energy Landscape", pad=20)
    ax.view_init(elev=style["elev"], azim=style["azim"])
    pane_rgba = to_rgba(style["fg"], style["pane_alpha"])
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_pane_color(pane_rgba)
        axis.pane.set_edgecolor(to_rgba(style["fg"], 0.35))
        axis._axinfo["grid"]["color"] = to_rgba(style["fg"], 0.15)
        axis.label.set_color(style["fg"])
    ax.tick_params(colors=style["fg"])

    if minima:
        for m in minima:
            # Single annotate() per minimum (circle bbox + number) -- see the
            # original pipeline notes: a separate ax.plot() marker on Axes3D
            # would be re-projected and land somewhere else.
            x2, y2, _ = proj3d.proj_transform(m["pc1"], m["pc2"], floor_z, ax.get_proj())
            ax.annotate(
                str(m["rank"]), xy=(x2, y2), ha="center", va="center",
                fontsize=9 * style["font_scale"], fontweight="bold",
                color=style["minima_text"], zorder=100,
                bbox=dict(boxstyle="circle,pad=0.15", facecolor=style["minima_face"],
                          edgecolor=style["minima_edge"], linewidth=1.4),
            )

    if colorbar:
        cbar = fig.colorbar(surf, ax=ax, shrink=0.6, pad=0.14)
        _style_colorbar(cbar, style, "Free energy (kcal/mol)", vrange=(vmin, vmax))
    return surf


def _point_density(pc1, pc2, bins=60):
    """Fast per-point density via a 2D histogram lookup (no KDE / scipy needed)."""
    hist, xe, ye = np.histogram2d(pc1, pc2, bins=bins)
    ix = np.clip(np.searchsorted(xe, pc1, side="right") - 1, 0, bins - 1)
    iy = np.clip(np.searchsorted(ye, pc2, side="right") - 1, 0, bins - 1)
    return hist[ix, iy]


def draw_scatter(fig, ax, pc1, pc2, style, variance_explained=None, colorbar=True, labels=True,
                 time_label="Frame index (simulation time)"):
    """PC1 vs PC2 scatter, coloured by time or by local density."""
    cmap = get_cmap(style["scatter_cmap"], "scatter")
    if style["scatter_color"] == "density":
        c = _point_density(pc1, pc2)
        order = np.argsort(c)                 # densest points drawn last (on top)
        pc1, pc2, c = pc1[order], pc2[order], c[order]
        cbar_label = "Point density (frames per bin)"
    else:
        c = np.arange(len(pc1))
        cbar_label = time_label

    sc = ax.scatter(pc1, pc2, c=c, cmap=cmap, s=style["scatter_size"],
                    alpha=style["scatter_alpha"], edgecolors="none")
    if colorbar:
        cbar = fig.colorbar(sc, ax=ax, pad=0.02)
        _style_colorbar(cbar, style, cbar_label)
    if labels:
        xl, yl = pc_axis_labels(variance_explained)
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        _title(ax, style, "PCA Projection (PC1 vs PC2)")
    ax.set_facecolor(style["bg"])
    ax.grid(True, alpha=max(style["grid_alpha"], 0.1) if style["grid"] else 0.15, linewidth=0.6)
    ax.set_axisbelow(True)
    return sc


def draw_scree(fig, ax1, eigvals, style, n_show=10, legend=True, labels=True):
    """Variance explained per PC (bars) + cumulative (line, right axis)."""
    eigvals = np.asarray(eigvals, dtype=float)
    n_show = min(n_show, len(eigvals))
    variance_pct = eigvals / eigvals.sum() * 100
    cumulative_pct = np.cumsum(variance_pct)
    x = np.arange(1, n_show + 1)

    ax1.bar(x, variance_pct[:n_show], color=style["bar_color"], edgecolor=style["fg"],
            linewidth=0.8, label="Individual")
    ax1.set_xticks(x)
    ax1.set_facecolor(style["bg"])

    ax2 = ax1.twinx()
    ax2.plot(x, cumulative_pct[:n_show], color=style["line_color"], marker="o",
             markersize=5, linewidth=2, label="Cumulative")
    ax2.set_ylim(0, 105)
    if style["open_spines"]:
        ax2.spines["top"].set_visible(False)
        ax1.spines["right"].set_visible(True)
    ax2.tick_params(colors=style["fg"])

    if labels:
        ax1.set_xlabel("Principal component")
        ax1.set_ylabel("Variance explained (%)")
        ax2.set_ylabel("Cumulative variance explained (%)")
        _title(ax1, style, "PCA Scree Plot")
    if legend:
        b, bl = ax1.get_legend_handles_labels()
        l, ll = ax2.get_legend_handles_labels()
        ax1.legend(b + l, bl + ll, loc="center right")
    if style["grid"]:
        ax1.grid(True, axis="y", alpha=style["grid_alpha"])
        ax1.set_axisbelow(True)
    return ax2


# ----------------------------------------------------------------------
# Full-figure wrappers used by pca_fel_analysis.py
# ----------------------------------------------------------------------

def plot_fel(x_centers, y_centers, fel, style, variance_explained=None,
             out_png="pca_fel.png", dpi=1200):
    with plt.rc_context(style_rc(style)):
        fig, ax = plt.subplots(figsize=(8, 6.5))
        draw_fel_2d(fig, ax, x_centers, y_centers, fel, style,
                    variance_explained=variance_explained)
        fig.tight_layout()
        _save(fig, out_png, dpi, "2D FEL plot")


def plot_fel_with_minima_labels(x_centers, y_centers, fel, minima, style,
                                variance_explained=None, out_png="pca_fel_labeled.png", dpi=1200):
    with plt.rc_context(style_rc(style)):
        fig, ax = plt.subplots(figsize=(8, 6.5))
        draw_fel_2d(fig, ax, x_centers, y_centers, fel, style, minima=minima,
                    variance_explained=variance_explained,
                    title="Free Energy Landscape with Minima Labeled")
        fig.tight_layout()
        _save(fig, out_png, dpi, "labeled FEL plot")


def plot_fel_3d(x_centers, y_centers, fel, style, minima=None, variance_explained=None,
                out_png="pca_fel_3d.png", dpi=1200):
    with plt.rc_context(style_rc(style)):
        fig = plt.figure(figsize=(9, 8))
        ax = fig.add_subplot(111, projection="3d")
        draw_fel_3d(fig, ax, x_centers, y_centers, fel, style, minima=minima,
                    variance_explained=variance_explained)
        fig.tight_layout()
        _save(fig, out_png, dpi, "3D FEL plot")


def plot_pca_scatter(pc1, pc2, style, variance_explained=None, out_png="pca_scatter.png",
                     dpi=1200):
    with plt.rc_context(style_rc(style)):
        fig, ax = plt.subplots(figsize=(8, 6.5))
        draw_scatter(fig, ax, pc1, pc2, style, variance_explained=variance_explained)
        fig.tight_layout()
        _save(fig, out_png, dpi, "PCA scatter plot")


def plot_scree(eigvals, style, n_show=10, out_png="pca_scree.png", dpi=1200):
    with plt.rc_context(style_rc(style)):
        fig, ax = plt.subplots(figsize=(8, 6))
        draw_scree(fig, ax, eigvals, style, n_show=n_show)
        fig.tight_layout()
        _save(fig, out_png, dpi, "scree plot")


# ----------------------------------------------------------------------
# Synthetic demo landscape (for previews only -- never used in a real run)
# ----------------------------------------------------------------------

def demo_data(seed=7, n_frames=6000):
    """
    A believable 3-basin trajectory in PC space: the 'protein' starts in
    basin A, hops to B, visits C, and returns -- so time-colouring and
    density-colouring both look like a real run.
    """
    rng = np.random.default_rng(seed)
    basins = [((-12.0, 3.0), (3.2, 2.4)), ((6.0, -4.0), (3.6, 2.8)),
              ((10.0, 9.0), (2.4, 2.2)), ((-12.0, 3.0), (3.0, 2.6))]
    weights = [0.38, 0.30, 0.14, 0.18]
    pcs = []
    for (center, spread), w in zip(basins, weights):
        n = int(n_frames * w)
        walk = np.cumsum(rng.normal(0, 0.35, size=(n, 2)), axis=0)
        walk -= np.linspace(0, 1, n)[:, None] * walk[-1]        # pull walk back to centre
        pts = np.array(center) + rng.normal(0, 1, size=(n, 2)) * spread * 0.85 + walk * 0.25
        pcs.append(pts)
    # short transition segments between consecutive basins
    full = []
    for k, seg in enumerate(pcs):
        full.append(seg)
        if k + 1 < len(pcs):
            a, b = seg[-1], pcs[k + 1][0]
            t = np.linspace(0, 1, 60)[:, None]
            full.append(a + (b - a) * t + rng.normal(0, 1.2, size=(60, 2)))
    proj = np.vstack(full)
    eigvals = np.array([31.0, 17.5, 9.8, 6.9, 5.1, 3.9, 3.1, 2.5, 2.0, 1.7, 1.4, 1.2] + [0.6] * 20)
    var = eigvals / eigvals.sum() * 100
    return proj[:, 0], proj[:, 1], eigvals, (var[0], var[1])


def _demo_minima(pc1, pc2, x_c, y_c, fel):
    out = []
    for rank, (e, r, c) in enumerate(find_local_minima(fel), start=1):
        d2 = (pc1 - x_c[c]) ** 2 + (pc2 - y_c[r]) ** 2
        i = int(np.argmin(d2))
        out.append({"rank": rank, "energy_kcal_mol": float(e), "pc1": float(pc1[i]),
                    "pc2": float(pc2[i]), "frame_idx": i})
    return out


# ----------------------------------------------------------------------
# Preview generation
# ----------------------------------------------------------------------

def render_style_preview(style_name, out_png, data=None, dpi=110):
    """2x2 composite (labeled 2D FEL, 3D FEL, PCA scatter, scree) for one style."""
    style = get_style(style_name)
    pc1, pc2, eigvals, var = data or demo_data()
    x_c, y_c, fel = compute_fel(pc1, pc2)
    minima = _demo_minima(pc1, pc2, x_c, y_c, fel)

    preview = dict(style)
    preview["font_scale"] = style["font_scale"] * 0.8
    with plt.rc_context(style_rc(preview)):
        fig = plt.figure(figsize=(13, 10.5))
        ax1 = fig.add_subplot(2, 2, 1)
        draw_fel_2d(fig, ax1, x_c, y_c, fel, preview, minima=minima, variance_explained=var,
                    title="Free Energy Landscape")
        ax2 = fig.add_subplot(2, 2, 2, projection="3d")
        draw_fel_3d(fig, ax2, x_c, y_c, fel, preview, minima=minima, variance_explained=var)
        ax3 = fig.add_subplot(2, 2, 3)
        draw_scatter(fig, ax3, pc1, pc2, preview, variance_explained=var)
        ax4 = fig.add_subplot(2, 2, 4)
        draw_scree(fig, ax4, eigvals, preview)
        fig.suptitle(f"{style['number']:02d} · {style['label']}", fontsize=18,
                     fontweight="bold", color=style["fg"], y=0.995)
        fig.tight_layout(rect=(0, 0, 1, 0.97))
        fig.savefig(out_png, dpi=dpi)
        plt.close(fig)


def render_style_thumbnail(style_name, out_png, data=None, dpi=100):
    """Small, label-free 2D FEL tile for the picker grid."""
    style = get_style(style_name)
    pc1, pc2, eigvals, var = data or demo_data()
    x_c, y_c, fel = compute_fel(pc1, pc2)
    minima = _demo_minima(pc1, pc2, x_c, y_c, fel)
    thumb = dict(style)
    thumb["font_scale"] = 0.55
    with plt.rc_context(style_rc(thumb)):
        fig, ax = plt.subplots(figsize=(2.4, 2.0))
        draw_fel_2d(fig, ax, x_c, y_c, fel, thumb, minima=minima, colorbar=False, labels=False)
        for line in ax.lines:
            line.set_markersize(5)
            line.set_markeredgewidth(1.0)
        for t in ax.texts:
            t.set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        fig.subplots_adjust(0.02, 0.02, 0.98, 0.98)
        fig.savefig(out_png, dpi=dpi, bbox_inches=None)
        plt.close(fig)


def render_gallery(out_png, data=None, dpi=130):
    """One contact sheet with all styles (4 columns x 5 rows) -- quick overview."""
    pc1, pc2, eigvals, var = data or demo_data()
    x_c, y_c, fel = compute_fel(pc1, pc2)
    minima = _demo_minima(pc1, pc2, x_c, y_c, fel)
    keys = list(STYLES.keys())
    ncol = 4
    nrow = int(np.ceil(len(keys) / ncol))
    fig = plt.figure(figsize=(ncol * 3.6, nrow * 3.4), facecolor="#f3f3f3")
    for i, key in enumerate(keys):
        style = get_style(key)
        tile = dict(style)
        tile["font_scale"] = 0.6
        with plt.rc_context(style_rc(tile)):
            ax = fig.add_subplot(nrow, ncol, i + 1)
            draw_fel_2d(fig, ax, x_c, y_c, fel, tile, minima=minima, colorbar=False,
                        labels=False)
            ax.set_xticks([])
            ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(True)
                sp.set_color("#888888")
            ax.set_title(f"{i + 1}. {style['label']}\n({key})", fontsize=9.5,
                         color="#222222", fontweight="bold", pad=5,
                         fontfamily="sans-serif")
    fig.suptitle("PCA / FEL pipeline — 20 figure styles  (use --style <name or number>)",
                 fontsize=15, fontweight="bold", color="#222222")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(out_png, dpi=dpi, facecolor="#f3f3f3")
    plt.close(fig)
    print(f"Saved style gallery to {out_png}")


def generate_previews(out_dir):
    """Write <key>_preview.png + <key>_thumb.png for every style, a gallery, and styles.json."""
    os.makedirs(out_dir, exist_ok=True)
    data = demo_data()
    for i, (key, label, desc, cat) in enumerate(list_styles(), start=1):
        render_style_thumbnail(key, os.path.join(out_dir, f"{key}_thumb.png"), data=data)
        render_style_preview(key, os.path.join(out_dir, f"{key}_preview.png"), data=data)
        print(f"  [{i:2d}/{len(STYLES)}] {key}")
    render_gallery(os.path.join(out_dir, "fel_style_gallery.png"), data=data)
    with open(os.path.join(out_dir, "styles.json"), "w") as f:
        json.dump([{"key": k, "label": l, "description": d, "category": c}
                   for k, l, d, c in list_styles()], f, indent=2)
    print(f"Previews written to {out_dir}")


def main(argv=None):
    ap = argparse.ArgumentParser(description="Generate FEL style previews.")
    ap.add_argument("--previews", metavar="DIR", help="write thumbnails + previews for all styles")
    ap.add_argument("--gallery", metavar="PNG", help="write one contact sheet of all styles")
    ap.add_argument("--one", nargs=2, metavar=("STYLE", "PNG"),
                    help="write the 2x2 preview of a single style")
    args = ap.parse_args(argv)
    if not (args.previews or args.gallery or args.one):
        ap.print_help()
        return
    if args.previews:
        generate_previews(args.previews)
    if args.gallery:
        render_gallery(args.gallery)
    if args.one:
        render_style_preview(args.one[0], args.one[1])
        print(f"Saved {args.one[1]}")


if __name__ == "__main__":
    main()
