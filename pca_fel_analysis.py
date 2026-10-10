"""
PCA-based Free Energy Landscape (FEL) analysis for a Desmond MD trajectory.

Computes PCA on protein backbone atoms, projects the trajectory onto
PC1/PC2, and constructs a 2D free energy landscape from the resulting
probability distribution.

Usage:
    $SCHRODINGER/run python3 pca_fel_analysis.py system-out.cms system_trj
    $SCHRODINGER/run python3 pca_fel_analysis.py system-out.cms system_trj --style neon_cyber
    $SCHRODINGER/run python3 pca_fel_analysis.py system-out.cms system_trj --style 18 --dpi 600
    $SCHRODINGER/run python3 pca_fel_analysis.py --list-styles

    Run from a terminal WITHOUT --style and you get a numbered menu of the
    20 figure styles to choose from first (see fel_style_previews/
    fel_style_gallery.png for what each looks like). From the PyMOL plugin,
    the style is picked in the GUI with a live preview.

Notes:
    - Requires a Schrodinger Python environment (run via $SCHRODINGER/run).
    - Method names on the PCA object (getEigvecs, getEigvals, calcModes)
      follow ProDy's NMA API, which schrodinger.trajectory.prody wraps.
      If your Schrodinger version differs slightly, run `dir(p)` on the
      PCA instance to confirm the exact accessor names.
    - Set --temperature to match the temperature used in your simulation.
    - Needs fel_styles.py and fel_plots.py in the same folder as this script.

License: MIT (see LICENSE file). Note that this code calls the Schrodinger
Suite Python API and PyMOL, which are separate third-party software not
covered by this license -- see LICENSE for details.
"""

import sys
import os
import json
import shutil
import argparse
import subprocess
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fel_styles import get_style, format_style_table, DEFAULT_STYLE, STYLES  # noqa: E402
import fel_plots  # noqa: E402  (sets the matplotlib Agg backend)
from fel_plots import compute_fel as _compute_fel, find_local_minima as _find_local_minima  # noqa: E402

from schrodinger.application.desmond.packages import topo, traj
from schrodinger.trajectory.prody import pca

# ----------------------------------------------------------------------
# Defaults (all overridable from the command line)
# ----------------------------------------------------------------------
DPI = 1200
TEMPERATURE = 300.0        # K -- set to match your simulation (--temperature)
N_BINS = 60                 # bins per PC axis for the FEL histogram (--bins)
MAX_MINIMA = 5               # how many distinct free-energy basins to extract
MINIMA_ENERGY_CUTOFF = 1.5   # kcal/mol above the global min -- basins above this are ignored
MIN_BASIN_SEPARATION = 8     # grid points apart two minima must be to count as distinct

STYLE = get_style(DEFAULT_STYLE)   # replaced in main() by --style / menu choice


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
    """Distinct FEL basins, lowest first, as (energy, row, col). See fel_plots."""
    return _find_local_minima(fel, min_distance=min_distance, max_minima=max_minima,
                              energy_cutoff=energy_cutoff)


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
                            color_scheme="rainbow", bg_color="white"):
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
bg_color {bg_color}
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
            render_structure_pymol(m["pdb_file"], out_png, color_scheme=color_scheme,
                                   bg_color=STYLE["pymol_bg"])
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


def compute_fel(pc1, pc2, n_bins=None, temperature=None):
    """FEL(x, y) = -kB*T ln(P/Pmax) on PC1/PC2 (kcal/mol, unsampled = NaN)."""
    return _compute_fel(pc1, pc2, n_bins=n_bins or N_BINS,
                        temperature=temperature or TEMPERATURE)


# ----------------------------------------------------------------------
# Plot wrappers -- all drawing lives in fel_plots.py and follows STYLE.
# ----------------------------------------------------------------------

def plot_scree(eigvals, n_show=10, out_png="pca_scree.png"):
    fel_plots.plot_scree(eigvals, STYLE, n_show=n_show, out_png=out_png, dpi=DPI)


def plot_fel(x_centers, y_centers, fel, variance_explained=None, out_png="pca_fel.png"):
    fel_plots.plot_fel(x_centers, y_centers, fel, STYLE,
                       variance_explained=variance_explained, out_png=out_png, dpi=DPI)


def plot_fel_with_minima_labels(x_centers, y_centers, fel, minima, variance_explained=None,
                                out_png="pca_fel_labeled.png"):
    fel_plots.plot_fel_with_minima_labels(x_centers, y_centers, fel, minima, STYLE,
                                          variance_explained=variance_explained,
                                          out_png=out_png, dpi=DPI)


def plot_fel_3d(x_centers, y_centers, fel, minima=None, variance_explained=None,
                out_png="pca_fel_3d.png"):
    fel_plots.plot_fel_3d(x_centers, y_centers, fel, STYLE, minima=minima,
                          variance_explained=variance_explained, out_png=out_png, dpi=DPI)


def plot_pca_scatter(pc1, pc2, variance_explained=None, out_png="pca_scatter.png"):
    fel_plots.plot_pca_scatter(pc1, pc2, STYLE, variance_explained=variance_explained,
                               out_png=out_png, dpi=DPI)


# ----------------------------------------------------------------------
# Command line / style selection
# ----------------------------------------------------------------------

def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        description="PCA + free energy landscape analysis of a Desmond trajectory.")
    ap.add_argument("cms_file", nargs="?", default="system-out.cms", help="*-out.cms file")
    ap.add_argument("trj_dir", nargs="?", default="system_trj", help="matching *_trj folder")
    ap.add_argument("--style", "-s", default=None,
                    help="figure style: name (e.g. neon_cyber) or number 1-%d" % len(STYLES))
    ap.add_argument("--list-styles", action="store_true", help="list the styles and exit")
    ap.add_argument("--dpi", type=int, default=DPI, help="output DPI (default %(default)s)")
    ap.add_argument("--temperature", "-T", type=float, default=TEMPERATURE,
                    help="simulation temperature in K (default %(default)s)")
    ap.add_argument("--bins", type=int, default=N_BINS,
                    help="FEL histogram bins per axis (default %(default)s)")
    return ap.parse_args(argv)


def choose_style_interactively():
    """Numbered style menu for terminal runs started without --style."""
    here = os.path.dirname(os.path.abspath(__file__))
    gallery = os.path.join(here, "fel_style_previews", "fel_style_gallery.png")
    print("\nChoose a figure style for all plots:")
    print(format_style_table())
    if os.path.exists(gallery):
        print(f"\n  (Visual overview of all styles: {gallery})")
    while True:
        try:
            ans = input(f"\nStyle number or name [Enter = 1, {DEFAULT_STYLE}]: ").strip()
        except EOFError:
            ans = ""
        try:
            return get_style(ans or DEFAULT_STYLE)
        except KeyError as e:
            print(f"  {e}")


def main(argv=None):
    global STYLE, DPI, TEMPERATURE, N_BINS
    args = parse_args(argv)

    if args.list_styles:
        print(format_style_table())
        return

    if args.style is not None:
        try:
            STYLE = get_style(args.style)
        except KeyError as e:
            sys.exit(f"ERROR: {e}")
    elif sys.stdin is not None and sys.stdin.isatty():
        STYLE = choose_style_interactively()
    else:
        STYLE = get_style(DEFAULT_STYLE)

    DPI, TEMPERATURE, N_BINS = args.dpi, args.temperature, args.bins
    print(f"Figure style: {STYLE['number']}. {STYLE['label']} ({STYLE['key']}) | "
          f"{DPI} dpi | T = {TEMPERATURE} K | {N_BINS} bins")

    with open("run_settings.json", "w") as f:
        json.dump({"style": STYLE["key"], "dpi": DPI, "temperature_K": TEMPERATURE,
                   "bins": N_BINS, "cms_file": os.path.abspath(args.cms_file),
                   "trj_dir": os.path.abspath(args.trj_dir)}, f, indent=2)

    print(f"Loading {args.cms_file} / {args.trj_dir} ...")
    msys_model, cms_model, tr = load_trajectory(args.cms_file, args.trj_dir)

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
