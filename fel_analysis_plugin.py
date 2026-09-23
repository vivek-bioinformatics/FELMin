"""
PyMOL Plugin: FEL Analysis (full pipeline orchestrator)
=========================================================

IMPORTANT ARCHITECTURAL NOTE -- read this first:

PyMOL's own Python interpreter cannot run the PCA/FEL computation itself.
The trajectory loading (topo.read_cms, traj.read_traj) and PCA
(schrodinger.trajectory.prody.pca) used by pca_fel_analysis.py depend on
Schrodinger's proprietary Python packages, which only exist inside
Schrodinger's own bundled Python environment ($SCHRODINGER/run python3) --
they are not, and cannot be, importable from PyMOL's separate Python.
Desmond's -out.cms / _trj trajectory format also isn't something PyMOL
can read natively.

So this plugin works as an ORCHESTRATOR:
  1. It launches pca_fel_analysis.py as a subprocess through
     $SCHRODINGER/run, exactly as you would from the command line.
  2. That script does the PCA, builds the FEL, finds minima, extracts
     structures, and (via its own PyMOL subprocess call) renders them.
  3. Once it finishes, THIS plugin loads the resulting structures and
     images back into your live PyMOL session automatically, so you get
     a one-click "run everything and show me the results" experience
     inside PyMOL, without ever leaving it or touching a terminal.

Install like any other PyMOL plugin: Plugin > Plugin Manager > Install
New Plugin > Choose file..., select this .py file, restart PyMOL, then
use Plugin > FEL Analysis. It also works headlessly (no GUI) via the
commands registered near the bottom of this file.

Place pca_fel_analysis.py in the SAME folder as this plugin file (the
default ANALYSIS_SCRIPT path below assumes that); or point at it
explicitly with the analysis_script argument / GUI field.

License: MIT (see LICENSE file). Note that this plugin calls the
Schrodinger Suite Python API and PyMOL, which are separate third-party
software not covered by this license -- see LICENSE for details.
"""

import os
import sys
import glob
import subprocess

from pymol import cmd

# Default location of the analysis script: same folder as this plugin file.
ANALYSIS_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "pca_fel_analysis.py")

RESULT_IMAGES = [
    "pca_fel.png", "pca_fel_3d.png", "pca_fel_labeled.png",
    "pca_scatter.png", "pca_scree.png",
]


# ----------------------------------------------------------------------
# Environment handling (same PATH/LD_LIBRARY_PATH cleanup pca_fel_analysis.py
# needed for its own PyMOL subprocess call -- kept here too since we're now
# going the OTHER direction: PyMOL calling OUT to Schrodinger's environment,
# which has the equivalent opposite problem if not handled carefully).
# ----------------------------------------------------------------------

def _schrodinger_run_executable(schrodinger_home=None):
    """Resolve the path to $SCHRODINGER/run, from an explicit arg or the
    SCHRODINGER environment variable."""
    home = schrodinger_home or os.environ.get("SCHRODINGER")
    if not home:
        raise RuntimeError(
            "Schrodinger installation path not found. Either set the "
            "SCHRODINGER environment variable before launching PyMOL, or "
            "pass schrodinger_home explicitly (e.g. "
            "run_fel_pipeline(..., schrodinger_home='/opt/schrodingerYYYY-N')."
        )
    run_exe = os.path.join(home, "run")
    if not os.path.exists(run_exe):
        raise RuntimeError(f"Could not find '{run_exe}' -- check the Schrodinger path.")
    return run_exe


# ----------------------------------------------------------------------
# Pipeline orchestration
# ----------------------------------------------------------------------

def run_fel_pipeline(cms_file, trj_dir, output_dir=None, schrodinger_home=None,
                      analysis_script=None, load_results=True):
    """
    Run the full PCA/FEL pipeline (pca_fel_analysis.py) via Schrodinger's
    Python, then load the results back into the current PyMOL session.

    USAGE (from the PyMOL command line or another script):
        run_fel_pipeline system-out.cms, system_trj
        run_fel_pipeline system-out.cms, system_trj, output_dir=/path/to/results

    :param cms_file: path to the *-out.cms file.
    :param trj_dir: path to the matching *_trj trajectory directory.
    :param output_dir: directory to run in / write results to (default:
        the cms file's own directory).
    :param schrodinger_home: path to the Schrodinger installation
        (default: the SCHRODINGER environment variable).
    :param analysis_script: path to pca_fel_analysis.py (default: assumes
        it sits next to this plugin file).
    :param load_results: if True (default), automatically load the
        extracted minima structures and open the result images once the
        pipeline finishes.
    """
    cms_file = os.path.abspath(cms_file)
    trj_dir = os.path.abspath(trj_dir)
    output_dir = os.path.abspath(output_dir) if output_dir else os.path.dirname(cms_file)
    analysis_script = analysis_script or ANALYSIS_SCRIPT

    if not os.path.exists(analysis_script):
        raise RuntimeError(
            f"Could not find pca_fel_analysis.py at '{analysis_script}'. Place "
            f"it next to this plugin file, or pass analysis_script explicitly."
        )

    run_exe = _schrodinger_run_executable(schrodinger_home)

    print(f"[FEL Analysis] Running pipeline: {analysis_script}")
    print(f"[FEL Analysis]   cms_file  = {cms_file}")
    print(f"[FEL Analysis]   trj_dir   = {trj_dir}")
    print(f"[FEL Analysis]   output_dir = {output_dir}")

    cmd_args = [run_exe, "python3", analysis_script, cms_file, trj_dir]
    result = subprocess.run(
        cmd_args,
        cwd=output_dir,
        capture_output=True,
        text=True,
    )

    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr)
        raise RuntimeError(
            f"pca_fel_analysis.py exited with code {result.returncode}. "
            f"See the output above for details."
        )

    print(f"[FEL Analysis] Pipeline finished successfully. Results in {output_dir}")

    if load_results:
        load_fel_results(output_dir)

    return output_dir


def load_fel_results(output_dir):
    """
    Load a completed pipeline run's results into the current PyMOL session:
    every minimum_*_frame*.pdb is loaded as its own object with the same
    cartoon+ligand styling used elsewhere in this plugin, and the FEL/PCA
    PNG images are opened with the system's default image viewer.

    USAGE:
        load_fel_results /path/to/results_dir
    """
    output_dir = os.path.abspath(output_dir)

    pdb_files = sorted(glob.glob(os.path.join(output_dir, "minimum_*_frame*.pdb")))
    for pdb_path in pdb_files:
        obj_name = os.path.splitext(os.path.basename(pdb_path))[0]
        obj_name = obj_name.replace("-", "_").replace(".", "_")
        cmd.load(pdb_path, obj_name)
        cmd.remove(f"{obj_name} and solvent")
        cmd.hide("everything", obj_name)
        cmd.show("cartoon", f"{obj_name} and polymer")
        cmd.spectrum("count", "rainbow", f"{obj_name} and polymer and name CA")
        cmd.show("sticks", f"{obj_name} and organic")
        cmd.color("magenta", f"{obj_name} and organic and elem C")
        cmd.util.cnc(f"{obj_name} and organic")
        cmd.disable(obj_name)   # loaded but hidden by default -- enable ones you want to view

    if pdb_files:
        cmd.orient()
        print(f"[FEL Analysis] Loaded {len(pdb_files)} minima structure(s) into PyMOL "
              f"(disabled by default -- enable individual objects to view them).")
    else:
        print(f"[FEL Analysis] No minimum_*_frame*.pdb files found in {output_dir}")

    for img_name in RESULT_IMAGES:
        img_path = os.path.join(output_dir, img_name)
        if os.path.exists(img_path):
            _open_with_system_viewer(img_path)


def _open_with_system_viewer(path):
    """Open an image with the OS's default viewer (best-effort, non-blocking)."""
    try:
        if sys.platform.startswith("darwin"):
            subprocess.Popen(["open", path])
        elif os.name == "nt":
            os.startfile(path)  # noqa: this only runs on Windows
        else:
            subprocess.Popen(["xdg-open", path])
        print(f"[FEL Analysis] Opened {path}")
    except Exception as e:
        print(f"[FEL Analysis] Could not auto-open {path} ({e}) -- open it manually.")


# ----------------------------------------------------------------------
# Standalone rendering commands (for re-rendering structures without
# re-running the whole pipeline -- e.g. after tweaking colors/orientation)
# ----------------------------------------------------------------------

def render_fel_minimum(pdb_path, out_png=None, width=1000, height=1000,
                        dpi=600, color_scheme="rainbow"):
    """
    Render a single minimum-energy structure as a cartoon (+ ligand) image.

    USAGE:
        render_fel_minimum minimum_1_frame5294.pdb
    """
    width, height, dpi = int(width), int(height), int(dpi)
    if out_png is None:
        out_png = os.path.splitext(pdb_path)[0] + ".png"

    obj_name = "render_" + os.path.splitext(os.path.basename(pdb_path))[0]
    obj_name = obj_name.replace("-", "_").replace(".", "_")

    cmd.load(pdb_path, obj_name)
    cmd.remove(f"{obj_name} and solvent")
    cmd.hide("everything", obj_name)
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 0)
    cmd.set("antialias", 2)
    cmd.set("ray_shadows", 0)
    cmd.set("ambient", 0.35)
    cmd.set("specular", 0.25)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.show("cartoon", f"{obj_name} and polymer")

    if color_scheme == "rainbow":
        cmd.spectrum("count", "rainbow", f"{obj_name} and polymer and name CA")
    else:
        cmd.util.cbc(obj_name)

    cmd.show("sticks", f"{obj_name} and organic")
    cmd.set("stick_radius", 0.22, f"{obj_name} and organic")
    cmd.color("magenta", f"{obj_name} and organic and elem C")
    cmd.util.cnc(f"{obj_name} and organic")

    cmd.orient(obj_name)
    cmd.zoom(obj_name, 3)
    cmd.ray(width, height)
    cmd.png(out_png, dpi=dpi)
    cmd.delete(obj_name)

    print(f"[FEL Analysis] Saved {out_png}")
    return out_png


def render_fel_minima_batch(input_dir, pattern="minimum_*_frame*.pdb",
                             output_dir=None, color_scheme="rainbow"):
    """
    Re-render every minimum-energy structure found in input_dir.

    USAGE:
        render_fel_minima_batch /path/to/minima_dir
    """
    output_dir = output_dir or input_dir
    pdb_files = sorted(glob.glob(os.path.join(input_dir, pattern)))
    if not pdb_files:
        print(f"[FEL Analysis] No files matching '{pattern}' found in {input_dir}")
        return []

    results = []
    for pdb_path in pdb_files:
        base = os.path.splitext(os.path.basename(pdb_path))[0]
        out_png = os.path.join(output_dir, base + ".png")
        results.append(render_fel_minimum(pdb_path, out_png, color_scheme=color_scheme))

    print(f"[FEL Analysis] Rendered {len(results)} structure(s) to {output_dir}")
    return results


cmd.extend("run_fel_pipeline", run_fel_pipeline)
cmd.extend("load_fel_results", load_fel_results)
cmd.extend("render_fel_minimum", render_fel_minimum)
cmd.extend("render_fel_minima_batch", render_fel_minima_batch)


# ----------------------------------------------------------------------
# GUI (only used when PyMOL loads this as a plugin through the Plugin
# Manager -- not needed for the headless commands above)
# ----------------------------------------------------------------------

def __init_plugin__(app=None):
    from pymol.plugins import addmenuitemqt
    addmenuitemqt("FEL Analysis", run_gui)


def run_gui():
    from pymol.Qt import QtWidgets

    dialog = QtWidgets.QDialog()
    dialog.setWindowTitle("FEL Analysis")
    layout = QtWidgets.QFormLayout(dialog)

    def _file_row(label, is_dir=False):
        line = QtWidgets.QLineEdit()
        btn = QtWidgets.QPushButton("Browse...")
        row = QtWidgets.QHBoxLayout()
        row.addWidget(line)
        row.addWidget(btn)
        container = QtWidgets.QWidget()
        container.setLayout(row)
        layout.addRow(label, container)

        def _browse():
            if is_dir:
                path = QtWidgets.QFileDialog.getExistingDirectory(dialog, label)
            else:
                path, _ = QtWidgets.QFileDialog.getOpenFileName(dialog, label)
            if path:
                line.setText(path)

        btn.clicked.connect(_browse)
        return line

    cms_line = _file_row("-out.cms file")
    trj_line = _file_row("_trj folder", is_dir=True)
    out_line = _file_row("Output folder (optional)", is_dir=True)

    schrodinger_line = QtWidgets.QLineEdit(os.environ.get("SCHRODINGER", ""))
    layout.addRow("Schrodinger home", schrodinger_line)

    log_box = QtWidgets.QPlainTextEdit()
    log_box.setReadOnly(True)
    log_box.setMinimumHeight(200)
    layout.addRow(log_box)

    run_btn = QtWidgets.QPushButton("Run Full FEL Analysis")
    layout.addRow(run_btn)

    def _run():
        cms_file = cms_line.text().strip()
        trj_dir = trj_line.text().strip()
        output_dir = out_line.text().strip() or None
        schrodinger_home = schrodinger_line.text().strip() or None

        if not cms_file or not trj_dir:
            QtWidgets.QMessageBox.warning(dialog, "FEL Analysis",
                                           "Please select both the -out.cms file and the _trj folder.")
            return

        log_box.appendPlainText("Running pipeline -- this can take a while for large trajectories...")
        QtWidgets.QApplication.processEvents()

        try:
            result_dir = run_fel_pipeline(
                cms_file, trj_dir,
                output_dir=output_dir,
                schrodinger_home=schrodinger_home,
            )
            log_box.appendPlainText(f"Done. Results in {result_dir}")
            QtWidgets.QMessageBox.information(dialog, "FEL Analysis",
                                               f"Pipeline finished.\nResults loaded from:\n{result_dir}")
        except Exception as e:
            log_box.appendPlainText(f"ERROR: {e}")
            QtWidgets.QMessageBox.critical(dialog, "FEL Analysis", str(e))

    run_btn.clicked.connect(_run)

    dialog.resize(600, 420)
    dialog.exec_()
