"""
PCA-based Free Energy Landscape (FEL) analysis for a Desmond MD trajectory.

Computes PCA on protein backbone atoms, projects the trajectory onto
PC1/PC2, and constructs a 2D free energy landscape from the resulting
probability distribution.

Usage:
    $SCHRODINGER/run python3 pca_fel_analysis.py system-out.cms system_trj

Notes:
    - Requires a Schrodinger Python environment (run via $SCHRODINGER/run).
    - Method names on the PCA object (getEigvecs, getEigvals, calcModes)
      follow ProDy's NMA API, which schrodinger.trajectory.prody wraps.
      If your Schrodinger version differs slightly, run `dir(p)` on the
      PCA instance to confirm the exact accessor names.
    - Adjust TEMPERATURE to match the temperature used in your simulation.

License: MIT (see LICENSE file). Note that this code calls the Schrodinger
Suite Python API and PyMOL, which are separate third-party software not
covered by this license -- see LICENSE for details.
"""

import sys
import os
import json
import shutil
import subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")  # non-interactive backend, safe for batch/cluster runs
import matplotlib.pyplot as plt
import matplotlib.patheffects as patheffects

from schrodinger.application.desmond.packages import topo, traj
from schrodinger.trajectory.prody import pca

# ----------------------------------------------------------------------
# Plot styling
# ----------------------------------------------------------------------
DPI = 1200
FEL_CMAP = "RdYlBu_r"       # blue (low energy) -> red (high energy), smooth gradient -- used
                             # consistently across pca_fel.png, pca_fel_3d.png, and
                             # pca_fel_labeled.png so all three FEL views match visually
SCATTER_CMAP = "plasma"     # vivid time-progression colormap (frame index, not energy --
                             # deliberately different from FEL_CMAP since it encodes a
                             # different quantity; let me know if you'd rather it match too)
N_LEVELS = 60                # more contour levels = smoother color gradient, less banding

plt.rcParams.update({
    "figure.dpi": 150,          # on-screen/preview only; savefig always uses DPI below
    "savefig.dpi": DPI,
    "font.size": 13,
    "font.family": "sans-serif",
    "axes.titlesize": 16,
    "axes.titleweight": "bold",
    "axes.labelsize": 14,
    "axes.labelweight": "bold",
    "axes.linewidth": 1.2,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.bbox": "tight",
})

# ----------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------
KB = 0.0019872041          # Boltzmann constant, kcal/(mol*K)
TEMPERATURE = 300.0        # K -- set to match your simulation
N_BINS = 60                 # bins per PC axis for the FEL histogram
MAX_MINIMA = 5               # how many distinct free-energy basins to extract
MINIMA_ENERGY_CUTOFF = 1.5   # kcal/mol above the global min -- basins above this are ignored
MIN_BASIN_SEPARATION = 8     # grid points apart two minima must be to count as distinct

CMS_FILE = sys.argv[1] if len(sys.argv) > 1 else "system-out.cms"
TRJ_DIR = sys.argv[2] if len(sys.argv) > 2 else "system_trj"


def load_trajectory(cms_file, trj_dir):
    """Load the .cms model and trajectory frames."""
    msys_model, cms_model = topo.read_cms(cms_file)
    tr = traj.read_traj(trj_dir)
    return msys_model, cms_model, tr


def get_backbone_coordsets(msys_model, cms_model, tr):
    """
    Extract backbone atom coordinates for every frame.
    Returns an (n_frames, n_atoms, 3) numpy array plus the atom ids used.
    """
    aslexpr = "backbone and protein"
    backbone_aids = cms_model.select_atom(aslexpr)

    if not backbone_aids:
        raise ValueError(
            "No backbone atoms found -- check your ASL selection or "
            "confirm the system contains a protein chain."
        )

    gids = topo.aids2gids(cms_model, backbone_aids, include_pseudoatoms=False)

    n_frames = len(tr)
    n_atoms = len(gids)
    coordsets = np.zeros((n_frames, n_atoms, 3), dtype=np.float64)

    for i, frame in enumerate(tr):
        coordsets[i] = frame.pos(gids)

    return coordsets, backbone_aids


def find_local_minima(fel, min_distance=MIN_BASIN_SEPARATION,
                       max_minima=MAX_MINIMA, energy_cutoff=MINIMA_ENERGY_CUTOFF):
    """
    Find distinct local minima (basins) in a 2D free energy grid using a
    simple neighborhood comparison (no scipy dependency).

    Only minima within `energy_cutoff` kcal/mol of the global minimum are
    considered (filters out noise far up the walls of the landscape), and
    minima closer than `min_distance` grid points to an already-accepted
    one are discarded so nearby noisy pixels don't count as separate basins.

    Returns a list of (energy, row_idx, col_idx) tuples, lowest energy first.
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
            i0, i1 = max(0, i - 1), min(ny, i + 2)
            j0, j1 = max(0, j - 1), min(nx, j + 2)
            neighborhood = fel[i0:i1, j0:j1]
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


def nearest_frame(pc1, pc2, target_x, target_y):
    """Return the index of the trajectory frame closest to (target_x, target_y) in PC space."""
    d2 = (pc1 - target_x) ** 2 + (pc2 - target_y) ** 2
    return int(np.argmin(d2))


def extract_structure(cms_model, tr, frame_idx, out_path):
    """
    Write out the full-system structure (protein, ligand, waters, etc.) at
    a given trajectory frame index. Format is inferred from out_path's
    extension (e.g. '.mae' or '.pdb').
    """
    fsys_ct = cms_model.fsys_ct.copy()
    topo.update_ct(fsys_ct, cms_model, tr[frame_idx])
    fsys_ct.write(out_path)


def extract_minima_structures(cms_model, tr, x_centers, y_centers, fel, pc1, pc2):
    """
    Locate the distinct free-energy basins in the FEL, find the trajectory
    frame closest to each basin's center, and save that frame's full
    structure to disk. Returns a list of dicts describing each saved minimum.
    """
    minima = find_local_minima(fel)
    results = []

    for rank, (energy, row, col) in enumerate(minima, start=1):
        target_x = x_centers[col]
        target_y = y_centers[row]
        frame_idx = nearest_frame(pc1, pc2, target_x, target_y)

        out_mae = f"minimum_{rank}_frame{frame_idx}.mae"
        out_pdb = f"minimum_{rank}_frame{frame_idx}.pdb"
        extract_structure(cms_model, tr, frame_idx, out_mae)
        extract_structure(cms_model, tr, frame_idx, out_pdb)

        results.append({
            "rank": rank,
            "energy_kcal_mol": float(energy),
            "frame_idx": frame_idx,
            "pc1": float(pc1[frame_idx]),
            "pc2": float(pc2[frame_idx]),
            "mae_file": out_mae,
            "pdb_file": out_pdb,
        })
        print(f"  Minimum #{rank}: {energy:.2f} kcal/mol -> frame {frame_idx} "
              f"(PC1={pc1[frame_idx]:.1f}, PC2={pc2[frame_idx]:.1f}) "
              f"saved to {out_mae} / {out_pdb}")

    with open("minima_summary.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"  Saved minima summary to minima_summary.json "
          f"(used by the PyMOL Minima Viewer plugin)")

    return results


def _clean_subprocess_env():
    """
    Build a subprocess environment stripped of Schrodinger-specific
    variables. $SCHRODINGER/run rewrites PYTHONHOME/PYTHONPATH and
    prepends its own bin directories to PATH for its own Python process --
    if that same environment is inherited by an external tool's launcher
    script (like PyMOL's, which itself invokes `python3`), the external
    tool ends up accidentally running under Schrodinger's bundled Python
    instead of its own, and fails with a missing-module error even though
    it's installed correctly. Stripping these lets PyMOL resolve its own
    interpreter normally.

    It also sets LD_LIBRARY_PATH to point at its own bundled shared
    libraries (including an older libstdc++.so.6) -- if PyMOL's compiled
    extensions inherit that, the dynamic linker finds Schrodinger's older
    libstdc++ before the system's newer one, and PyMOL fails with a
    GLIBCXX version-not-found ImportError even though the system library
    it actually needs is present. Stripping LD_LIBRARY_PATH (and
    LD_PRELOAD, for the same reason) lets the system's default linker
    search path resolve the correct libraries.
    """
    env = os.environ.copy()
    for var in ("PYTHONHOME", "PYTHONPATH", "LD_LIBRARY_PATH", "LD_PRELOAD"):
        env.pop(var, None)

    system_dirs = ["/usr/local/bin", "/usr/bin", "/bin"]
    current_path = env.get("PATH", "")
    env["PATH"] = os.pathsep.join(system_dirs + [current_path])
    return env


def render_structure_pymol(pdb_path, out_png, width=1000, height=1000, dpi=600,
                            color_scheme="rainbow"):
    """
    Render a cartoon image of a structure using PyMOL in headless batch mode
    (`pymol -cq`), run as a subprocess so it never has to be imported into
    the Schrodinger Python environment (which avoids version/DLL conflicts
    between the two separate Python installations).

    :param pdb_path: input structure file (the .pdb written by extract_structure).
    :param out_png: output image path.
    :param color_scheme: 'rainbow' (N-to-C spectrum, matches the reference
        figure's coloring) or 'chain' (color by chain instead).
    """
    # $SCHRODINGER/run launches Python in its own isolated environment, which
    # does NOT inherit your normal shell's PATH -- so even if `pymol` works
    # fine in a regular terminal, shutil.which() may not find it here. Set
    # PYMOL_EXECUTABLE to override with an absolute path if that happens:
    #     PYMOL_EXECUTABLE=/full/path/to/pymol $SCHRODINGER/run python3 ...
    pymol_exe = os.environ.get("PYMOL_EXECUTABLE") or shutil.which("pymol")
    if pymol_exe is None:
        raise RuntimeError(
            "Could not find the 'pymol' executable. This usually means "
            "$SCHRODINGER/run's isolated environment doesn't see the same "
            "PATH as your normal shell. Run `which pymol` in a regular "
            "terminal (outside $SCHRODINGER/run) to get its full path, then "
            "re-run this script as: "
            "PYMOL_EXECUTABLE=/full/path/to/pymol $SCHRODINGER/run python3 ..."
        )

    color_cmd = (
        "spectrum count, rainbow, mol and polymer and name CA"
        if color_scheme == "rainbow" else
        "util.cbc mol"
    )

    script = f"""
load {pdb_path}, mol
remove solvent
hide everything
bg_color white
set ray_opaque_background, 0
set antialias, 2
set ray_shadows, 0
set ambient, 0.35
set specular, 0.25
set cartoon_fancy_helices, 1
show cartoon, polymer
{color_cmd}
show sticks, organic
set stick_radius, 0.22, organic
color magenta, organic and elem C
util.cnc organic
orient mol
zoom mol, 3
ray {width}, {height}
png {out_png}, dpi={dpi}
quit
"""
    script_path = os.path.splitext(out_png)[0] + "_render.pml"
    with open(script_path, "w") as f:
        f.write(script)

    log_path = os.path.splitext(out_png)[0] + "_pymol.log"
    result = subprocess.run(
        [pymol_exe, "-cq", script_path],
        env=_clean_subprocess_env(),
        capture_output=True,
        text=True,
    )
    with open(log_path, "w") as f:
        f.write("--- stdout ---\n" + result.stdout + "\n--- stderr ---\n" + result.stderr)

    # PyMOL can exit with code 0 even when a command inside the script
    # failed (e.g. a headless-display issue during `ray`), so checking the
    # return code alone isn't reliable -- verify the PNG actually exists.
    if result.returncode != 0 or not os.path.exists(out_png):
        raise RuntimeError(
            f"PyMOL did not produce {out_png} (exit code {result.returncode}). "
            f"Full PyMOL output saved to {log_path} -- check it for the real "
            f"error. A common cause on headless servers is PyMOL needing a "
            f"virtual display for 'ray'; try prefixing the whole command with "
            f"`xvfb-run -a` (e.g. install with `apt install xvfb` / "
            f"`conda install -c conda-forge xorg-x11-server-xvfb`), or run "
            f"`{pymol_exe} -cq {script_path}` by hand to see the error directly."
        )


def render_all_minima_structures(minima, color_scheme="rainbow"):
    """
    Render every extracted minimum's structure with PyMOL and return a
    {rank: png_path} dict ready to hand to plot_fel_with_structures().
    Minima whose PyMOL render fails are skipped (with a printed warning)
    rather than aborting the whole run.
    """
    structure_image_paths = {}
    for m in minima:
        out_png = f"structure_{m['rank']}_frame{m['frame_idx']}.png"
        try:
            render_structure_pymol(m["pdb_file"], out_png, color_scheme=color_scheme)
            structure_image_paths[m["rank"]] = out_png
            print(f"  Rendered minimum #{m['rank']} -> {out_png}")
        except (RuntimeError, subprocess.CalledProcessError) as e:
            print(f"  Warning: PyMOL render failed for minimum #{m['rank']}: {e}")

    return structure_image_paths


def kabsch_align(mobile, target):
    """
    Superpose `mobile` (N,3) onto `target` (N,3) using the Kabsch algorithm
    (optimal least-squares rotation + translation, no scaling).
    Returns the aligned mobile coordinates.
    """
    mobile_center = mobile.mean(axis=0)
    target_center = target.mean(axis=0)
    mobile_c = mobile - mobile_center
    target_c = target - target_center

    H = mobile_c.T @ target_c
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T

    return (R @ mobile_c.T).T + target_center


def align_trajectory(coordsets, n_iter=2):
    """
    Remove overall translation/rotation from every frame so PCA reflects
    only internal conformational motion, not rigid-body drift/tumbling.

    Aligns all frames to the first frame, then re-aligns to the resulting
    mean structure (repeated n_iter times) for a stable reference -- this
    is the standard iterative-mean-structure approach used by MDAnalysis,
    GROMACS, etc.
    """
    aligned = coordsets.copy()
    reference = aligned[0]

    for _ in range(n_iter):
        for i in range(aligned.shape[0]):
            aligned[i] = kabsch_align(aligned[i], reference)
        reference = aligned.mean(axis=0)   # update reference to converged mean

    return aligned


def run_pca(coordsets):
    """Build the covariance matrix and diagonalize to get PCA modes."""
    p = pca.PCA()
    p.buildCovariance(coordsets)   # covariance about the mean structure
    p.calcModes()                  # diagonalize -> eigenvectors/eigenvalues
    return p


def project_onto_pcs(p, coordsets, n_modes=2):
    """
    Project each frame's mean-deviation onto the first n_modes eigenvectors.
    Returns an (n_frames, n_modes) array of projections.
    """
    n_frames = coordsets.shape[0]
    mean = coordsets.mean(axis=0)
    flat_dev = (coordsets - mean).reshape(n_frames, -1)   # (n_frames, 3N)

    eigvecs = p.getEigvecs()[:, :n_modes]                 # (3N, n_modes)
    projections = flat_dev.dot(eigvecs)                   # (n_frames, n_modes)
    return projections


def compute_fel(pc1, pc2, n_bins=N_BINS, temperature=TEMPERATURE):
    """
    Compute a 2D free energy landscape from PC1/PC2 projections.

        FEL(x, y) = -kB*T * ln( P(x, y) / P_max )

    Returns bin centers (x, y) and the FEL grid in kcal/mol, with
    unsampled bins set to NaN so they're left blank when plotted.
    """
    hist, xedges, yedges = np.histogram2d(pc1, pc2, bins=n_bins)

    prob = hist / hist.sum()

    # Suppress the expected log(0)/divide-by-zero warnings for empty bins;
    # those positions are overwritten with NaN on the next line anyway.
    with np.errstate(divide="ignore", invalid="ignore"):
        fel = -KB * temperature * np.log(prob / prob.max())
    fel[hist == 0] = np.nan   # unsampled bins -> blank, not -inf/inf

    x_centers = 0.5 * (xedges[:-1] + xedges[1:])
    y_centers = 0.5 * (yedges[:-1] + yedges[1:])

    return x_centers, y_centers, fel.T   # transpose for correct x/y orientation


def _pc_axis_labels(variance_explained):
    """
    Build x/y axis label strings, including the % variance explained by
    each PC when available -- e.g. "PC1 (25.6%)" instead of just "PC1".
    :param variance_explained: (pc1_pct, pc2_pct) tuple, or None.
    """
    if variance_explained is None:
        return "PC1", "PC2"
    pc1_pct, pc2_pct = variance_explained
    return f"PC1 ({pc1_pct:.1f}%)", f"PC2 ({pc2_pct:.1f}%)"


def plot_scree(eigvals, n_show=10, out_png="pca_scree.png"):
    """
    Plot a scree plot: variance explained by each principal component
    (bars) plus cumulative variance explained (line) -- shows at a glance
    how many components actually matter, not just PC1/PC2 in isolation.

    :param eigvals: full eigenvalue array from p.getEigvals().
    :param n_show: how many leading components to plot (default 10).
    """
    eigvals = np.asarray(eigvals)
    n_show = min(n_show, len(eigvals))
    variance_pct = eigvals / eigvals.sum() * 100
    cumulative_pct = np.cumsum(variance_pct)

    fig, ax1 = plt.subplots(figsize=(8, 6))

    x = np.arange(1, n_show + 1)
    ax1.bar(x, variance_pct[:n_show], color="#4472C4", edgecolor="black",
            linewidth=0.8, label="Individual")
    ax1.set_xlabel("Principal component")
    ax1.set_ylabel("Variance explained (%)", fontweight="bold")
    ax1.set_xticks(x)
    ax1.set_facecolor("white")

    ax2 = ax1.twinx()
    ax2.plot(x, cumulative_pct[:n_show], color="#C00000", marker="o",
              markersize=5, linewidth=2, label="Cumulative")
    ax2.set_ylabel("Cumulative variance explained (%)", fontweight="bold")
    ax2.set_ylim(0, 105)

    # Combined legend from both axes
    bars, bar_labels = ax1.get_legend_handles_labels()
    lines, line_labels = ax2.get_legend_handles_labels()
    ax1.legend(bars + lines, bar_labels + line_labels, loc="center right")

    ax1.set_title("PCA Scree Plot", pad=12)
    fig.tight_layout()
    fig.savefig(out_png, dpi=DPI)
    plt.close(fig)
    print(f"Saved scree plot to {out_png} ({DPI} dpi)")


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
    Smooth the FEL grid and upsample it onto a much finer mesh, turning the
    blocky raw histogram into a crisp, continuous-looking surface -- the
    difference between a pixelated plot and a polished, publication-style one.
    Unsampled (NaN) regions are preserved as blank after upsampling.
    """
    nan_mask = ~np.isfinite(fel)
    fel_filled = np.where(nan_mask, np.nanmax(fel), fel)

    try:
        from scipy.ndimage import gaussian_filter, zoom

        fel_smooth = gaussian_filter(fel_filled, sigma=sigma)
        mask_smooth = gaussian_filter(nan_mask.astype(float), sigma=sigma)

        fel_fine = zoom(fel_smooth, upsample_factor, order=3)
        mask_fine = zoom(mask_smooth, upsample_factor, order=1)
        x_fine = np.linspace(x_centers[0], x_centers[-1], fel_fine.shape[1])
        y_fine = np.linspace(y_centers[0], y_centers[-1], fel_fine.shape[0])
    except ImportError:
        # Fallback: blur only (no upsampling) if scipy isn't available.
        fel_fine = _simple_box_blur(fel_filled, iterations=3)
        mask_fine = _simple_box_blur(nan_mask.astype(float), iterations=3)
        x_fine, y_fine = x_centers, y_centers

    fel_fine = np.where(mask_fine > 0.5, np.nan, fel_fine)
    return x_fine, y_fine, fel_fine


def plot_fel(x_centers, y_centers, fel, variance_explained=None, out_png="pca_fel.png"):
    """
    Plot the FEL as a filled contour map (2D, top-down view).

    :param variance_explained: optional (pc1_pct, pc2_pct) tuple -- if
        given, axis labels show "PC1 (25.6%)" instead of just "PC1".
    """
    x_centers, y_centers, fel = smooth_and_upsample_fel(x_centers, y_centers, fel)

    fig, ax = plt.subplots(figsize=(8, 6.5))

    finite_vals = fel[np.isfinite(fel)]
    if finite_vals.size == 0:
        raise ValueError(
            "FEL grid is entirely empty -- check your PC1/PC2 projections "
            "and n_bins (too many bins for too few frames can cause this)."
        )
    vmin, vmax = finite_vals.min(), finite_vals.max()
    if vmin == vmax:
        vmax = vmin + 1e-6   # guard against a perfectly flat landscape
    levels = np.linspace(vmin, vmax, N_LEVELS)

    cf = ax.contourf(x_centers, y_centers, fel, levels=levels, cmap=FEL_CMAP, extend="neither")
    ax.contour(x_centers, y_centers, fel, levels=levels[::6], colors="black",
               linewidths=0.25, alpha=0.35)

    cbar = fig.colorbar(cf, ax=ax, pad=0.02)
    cbar.set_label("Free energy (kcal/mol)", fontsize=13, fontweight="bold")
    cbar.ax.tick_params(labelsize=10)

    xlabel, ylabel = _pc_axis_labels(variance_explained)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("Free Energy Landscape (PC1 vs PC2)", pad=12)
    ax.set_facecolor("white")

    fig.tight_layout()
    fig.savefig(out_png, dpi=DPI)
    plt.close(fig)
    print(f"Saved 2D FEL plot to {out_png} ({DPI} dpi)")


def plot_fel_with_minima_labels(x_centers, y_centers, fel, minima, variance_explained=None,
                                 out_png="pca_fel_labeled.png"):
    """
    Plot the 2D FEL with each minimum's exact (PC1, PC2) site marked by a
    highlighted point and labeled with its rank number (1, 2, 3, ...) --
    no arrows, no structure thumbnails, just the sites called out directly
    on the landscape.

    :param minima: list of dicts as returned by extract_minima_structures(),
        each containing 'rank', 'pc1', 'pc2', 'energy_kcal_mol'.
    :param variance_explained: optional (pc1_pct, pc2_pct) tuple -- if
        given, axis labels show "PC1 (25.6%)" instead of just "PC1".
    """
    x_fine, y_fine, fel_fine = smooth_and_upsample_fel(x_centers, y_centers, fel)

    fig, ax = plt.subplots(figsize=(8, 6.5))

    finite_vals = fel_fine[np.isfinite(fel_fine)]
    vmin, vmax = finite_vals.min(), finite_vals.max()
    if vmin == vmax:
        vmax = vmin + 1e-6
    levels = np.linspace(vmin, vmax, N_LEVELS)

    cf = ax.contourf(x_fine, y_fine, fel_fine, levels=levels, cmap=FEL_CMAP, extend="neither")
    cbar = fig.colorbar(cf, ax=ax, pad=0.02)
    cbar.set_label("Free energy (kcal/mol)", fontsize=13, fontweight="bold")
    cbar.ax.tick_params(labelsize=10)

    for m in minima:
        x, y = m["pc1"], m["pc2"]

        # Highlighted marker at the exact site
        ax.plot(x, y, marker="o", markersize=11, markeredgecolor="black",
                 markerfacecolor="white", markeredgewidth=2.0, zorder=5)

        # Rank number placed just beside the marker, with a small white
        # halo (path effect) so it stays legible against any color in the
        # landscape underneath it.
        txt = ax.annotate(
            str(m["rank"]), xy=(x, y), xytext=(6, 6),
            textcoords="offset points",
            fontsize=13, fontweight="bold", color="black", zorder=6,
        )
        txt.set_path_effects([
            patheffects.withStroke(linewidth=3, foreground="white")
        ])

    xlabel, ylabel = _pc_axis_labels(variance_explained)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("Free Energy Landscape with Minima Labeled", pad=12)
    ax.set_facecolor("white")

    fig.tight_layout()
    fig.savefig(out_png, dpi=DPI)
    plt.close(fig)
    print(f"Saved labeled FEL plot to {out_png} ({DPI} dpi)")


def plot_fel_3d(x_centers, y_centers, fel, minima=None, variance_explained=None,
                 out_png="pca_fel_3d.png"):
    """
    Plot the FEL as a 3D surface with a projected 2D contour floor,
    in the style commonly used for MD free-energy-landscape figures
    (mountain-like surface above, colored 'island' contour below).

    :param minima: optional list of dicts as returned by
        extract_minima_structures() (each with 'rank', 'pc1', 'pc2',
        'energy_kcal_mol'). If given, each site is marked on the floor
        projection with a circled rank number.
    :param variance_explained: optional (pc1_pct, pc2_pct) tuple -- if
        given, axis labels show "PC1 (25.6%)" instead of just "PC1".
    """
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 (registers 3D projection)
    from mpl_toolkits.mplot3d import proj3d

    x_centers, y_centers, fel = smooth_and_upsample_fel(x_centers, y_centers, fel)

    X, Y = np.meshgrid(x_centers, y_centers)

    finite_vals = fel[np.isfinite(fel)]
    if finite_vals.size == 0:
        raise ValueError("FEL grid is entirely empty -- cannot make 3D plot.")
    vmin, vmax = finite_vals.min(), finite_vals.max()

    # Replace NaN (unsampled) bins with the max energy so the surface
    # doesn't have holes in it -- purely cosmetic, doesn't affect the data.
    fel_filled = np.where(np.isfinite(fel), fel, vmax)

    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(111, projection="3d")
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    surf = ax.plot_surface(
        X, Y, fel_filled,
        cmap=FEL_CMAP, rstride=2, cstride=2,
        linewidth=0, antialiased=True, alpha=0.97,
        vmin=vmin, vmax=vmax,
    )

    # Projected 2D contour "floor" beneath the surface, like the blue/green
    # island plots in published FEL figures.
    floor_z = vmin - 0.15 * (vmax - vmin)
    ax.contourf(X, Y, fel_filled, zdir="z", offset=floor_z,
                levels=np.linspace(vmin, vmax, N_LEVELS), cmap=FEL_CMAP)

    # Axis/view setup MUST happen before we compute any 2D screen
    # projections below -- proj3d.proj_transform() uses the axes' current
    # projection matrix (ax.get_proj()), which depends on the view angle
    # and axis limits already being finalized.
    ax.set_zlim(floor_z, vmax)
    xlabel, ylabel = _pc_axis_labels(variance_explained)
    ax.set_xlabel(xlabel, labelpad=10)
    ax.set_ylabel(ylabel, labelpad=10)
    ax.set_zlabel("Free energy (kcal/mol)", labelpad=10)
    ax.set_title("3D Free Energy Landscape", pad=20)
    ax.view_init(elev=35, azim=-60)   # adjust viewing angle to taste
    ax.xaxis.pane.set_edgecolor("gray")
    ax.yaxis.pane.set_edgecolor("gray")
    ax.zaxis.pane.set_edgecolor("gray")
    ax.xaxis.pane.set_alpha(0.05)
    ax.yaxis.pane.set_alpha(0.05)
    ax.zaxis.pane.set_alpha(0.05)

    if minima:
        for m in minima:
            x, y = m["pc1"], m["pc2"]
            label = str(m["rank"])

            # The circle and number are drawn together in a SINGLE
            # ax.annotate() call, using a circular bbox around the text,
            # rather than a separate marker + label. Earlier versions used
            # ax.plot(..., marker="o") for the circle and ax.annotate() for
            # the number, both fed the same proj3d-projected (x2, y2) --
            # but ax.plot() on an Axes3D does NOT treat those as
            # already-projected screen coordinates the way ax.annotate()
            # does; it re-runs them through the full 3D pipeline as if they
            # were raw (PC1, PC2, z=0) data, landing the circle somewhere
            # different from its own number. A single annotate() call has
            # no such mismatch, since there's only one coordinate transform
            # involved for both the shape and the text.
            x2, y2, _ = proj3d.proj_transform(x, y, floor_z, ax.get_proj())

            ax.annotate(
                label, xy=(x2, y2), ha="center", va="center",
                fontsize=9, fontweight="bold", color="black", zorder=100,
                bbox=dict(boxstyle="circle,pad=0.15", facecolor="white",
                          edgecolor="black", linewidth=1.4),
            )

    cbar = fig.colorbar(surf, ax=ax, shrink=0.6, pad=0.1)
    cbar.set_label("Free energy (kcal/mol)", fontsize=13, fontweight="bold")
    cbar.ax.tick_params(labelsize=10)

    fig.tight_layout()
    fig.savefig(out_png, dpi=DPI)
    plt.close(fig)
    print(f"Saved 3D FEL plot to {out_png} ({DPI} dpi)")


def plot_pca_scatter(pc1, pc2, variance_explained=None, out_png="pca_scatter.png"):
    """
    Scatter plot of every frame's PC1/PC2 projection, colored by
    simulation time (frame index), the standard companion plot to a FEL.

    :param variance_explained: optional (pc1_pct, pc2_pct) tuple -- if
        given, axis labels show "PC1 (25.6%)" instead of just "PC1".
    """
    fig, ax = plt.subplots(figsize=(8, 6.5))

    frame_idx = np.arange(len(pc1))
    sc = ax.scatter(pc1, pc2, c=frame_idx, cmap=SCATTER_CMAP, s=10,
                     alpha=0.65, edgecolors="none")

    cbar = fig.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label("Frame index (simulation time)", fontsize=13, fontweight="bold")
    cbar.ax.tick_params(labelsize=10)

    xlabel, ylabel = _pc_axis_labels(variance_explained)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title("PCA Projection (PC1 vs PC2)", pad=12)
    ax.grid(alpha=0.15)
    ax.set_facecolor("white")

    fig.tight_layout()
    fig.savefig(out_png, dpi=DPI)
    plt.close(fig)
    print(f"Saved PCA scatter plot to {out_png} ({DPI} dpi)")


def main():
    print(f"Loading {CMS_FILE} / {TRJ_DIR} ...")
    msys_model, cms_model, tr = load_trajectory(CMS_FILE, TRJ_DIR)

    print("Extracting backbone coordinates from trajectory ...")
    coordsets, backbone_aids = get_backbone_coordsets(msys_model, cms_model, tr)
    print(f"  {coordsets.shape[0]} frames, {coordsets.shape[1]} backbone atoms")

    print("Aligning frames (removing translation/rotation) ...")
    coordsets = align_trajectory(coordsets)

    print("Running PCA ...")
    p = run_pca(coordsets)

    eigvals = p.getEigvals()
    variance_pct = eigvals / eigvals.sum() * 100
    pc_variance = (variance_pct[0], variance_pct[1])   # (PC1 %, PC2 %) for axis labels
    print(f"  PC1 explains {variance_pct[0]:.1f}% of variance")
    print(f"  PC2 explains {variance_pct[1]:.1f}% of variance")

    print("Plotting scree plot (variance explained per component) ...")
    plot_scree(eigvals)

    print("Projecting trajectory onto PC1/PC2 ...")
    projections = project_onto_pcs(p, coordsets, n_modes=2)
    pc1, pc2 = projections[:, 0], projections[:, 1]

    np.savetxt(
        "pc1_pc2_projections.csv",
        projections,
        delimiter=",",
        header="PC1,PC2",
        comments="",
    )
    print("  Saved PC1/PC2 projections to pc1_pc2_projections.csv")

    print("Computing free energy landscape ...")
    x_centers, y_centers, fel = compute_fel(pc1, pc2)

    plot_fel(x_centers, y_centers, fel, variance_explained=pc_variance)
    plot_pca_scatter(pc1, pc2, variance_explained=pc_variance)

    print("Finding free-energy minima and extracting representative structures ...")
    minima = extract_minima_structures(cms_model, tr, x_centers, y_centers, fel, pc1, pc2)
    if not minima:
        print("  No distinct minima found -- try loosening MINIMA_ENERGY_CUTOFF "
              "or MIN_BASIN_SEPARATION at the top of the script.")
        plot_fel_3d(x_centers, y_centers, fel, variance_explained=pc_variance)
        return

    print("Labeling minima on the FEL (2D and 3D) ...")
    plot_fel_with_minima_labels(x_centers, y_centers, fel, minima, variance_explained=pc_variance)
    plot_fel_3d(x_centers, y_centers, fel, minima=minima, variance_explained=pc_variance)

    print("Rendering minima structures with PyMOL ...")
    structure_image_paths = render_all_minima_structures(minima)
    if not structure_image_paths:
        print("  No structures were rendered -- check that PyMOL is on PATH "
              "('pymol -cq' should run without error).")


if __name__ == "__main__":
    main()
