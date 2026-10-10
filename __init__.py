"""
PyMOL Plugin: FEL Analysis (full pipeline orchestrator + 20 figure styles)
===========================================================================

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
  1. You pick one of 20 figure styles in the GUI, with a live preview of
     all four figure types (2D FEL, 3D FEL, PCA scatter, scree plot).
  2. It launches pca_fel_analysis.py --style <your choice> as a subprocess
     through $SCHRODINGER/run, streaming its log into the dialog.
  3. That script does the PCA, builds the FEL, finds minima, extracts
     structures, draws every figure in the chosen style, and renders the
     minima with PyMOL.
  4. When it finishes, THIS plugin loads the minima structures and opens
     the result images automatically.

INSTALL (as a package, so the helper files come along):
    Plugin > Plugin Manager > Install New Plugin > Choose file...
    -> select fel_analysis_plugin.zip, restart PyMOL, then Plugin > FEL Analysis.
The zip contains this file (as fel_analysis/__init__.py) together with
pca_fel_analysis.py, fel_styles.py, fel_plots.py and fel_style_previews/.

Headless use from the PyMOL command line:
    list_fel_styles
    fel_style_gallery
    run_fel_pipeline system-out.cms, system_trj, style=neon_cyber
    run_fel_pipeline system-out.cms, system_trj, style=14, dpi=600, temperature=310

License: MIT (see LICENSE file). Note that this plugin calls the
Schrodinger Suite Python API and PyMOL, which are separate third-party
software not covered by this license -- see LICENSE for details.
"""

import os
import sys
import glob
import subprocess

from pymol import cmd

PLUGIN_DIR = os.environ.get("FEL_PIPELINE_DIR") or os.path.dirname(os.path.abspath(__file__))
if PLUGIN_DIR not in sys.path:
    sys.path.insert(0, PLUGIN_DIR)

from fel_styles import get_style, list_styles, format_style_table, DEFAULT_STYLE  # noqa: E402

ANALYSIS_SCRIPT = os.path.join(PLUGIN_DIR, "pca_fel_analysis.py")
PLOTS_SCRIPT = os.path.join(PLUGIN_DIR, "fel_plots.py")
PREVIEW_DIR = os.path.join(PLUGIN_DIR, "fel_style_previews")

RESULT_IMAGES = [
    "pca_fel.png", "pca_fel_3d.png", "pca_fel_labeled.png",
    "pca_scatter.png", "pca_scree.png",
]


# ----------------------------------------------------------------------
# Helpers
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
    if os.name == "nt" and not os.path.exists(run_exe):
        run_exe = os.path.join(home, "run.exe")
    if not os.path.exists(run_exe):
        raise RuntimeError(f"Could not find '{run_exe}' -- check the Schrodinger path.")
    return run_exe


def _pipeline_args(run_exe, cms_file, trj_dir, style, dpi, temperature, analysis_script):
    return [run_exe, "python3", analysis_script, cms_file, trj_dir,
            "--style", str(style), "--dpi", str(int(dpi)),
            "--temperature", str(float(temperature))]


def _open_with_system_viewer(path):
    """Open a file with the OS's default viewer (best-effort, non-blocking)."""
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
# Headless commands
# ----------------------------------------------------------------------

def list_fel_styles():
    """
    Print the 20 available figure styles.

    USAGE:
        list_fel_styles
    """
    print("[FEL Analysis] Available figure styles (use the name or the number):")
    print(format_style_table())


def fel_style_gallery():
    """
    Open the contact sheet showing all 20 styles side by side.

    USAGE:
        fel_style_gallery
    """
    path = os.path.join(PREVIEW_DIR, "fel_style_gallery.png")
    if os.path.exists(path):
        _open_with_system_viewer(path)
    else:
        print(f"[FEL Analysis] Gallery not found at {path}. Generate it with:\n"
              f"  $SCHRODINGER/run python3 {PLOTS_SCRIPT} --previews {PREVIEW_DIR}")


def run_fel_pipeline(cms_file, trj_dir, output_dir=None, style=DEFAULT_STYLE, dpi=1200,
                     temperature=300.0, schrodinger_home=None, analysis_script=None,
                     load_results=True):
    """
    Run the full PCA/FEL pipeline (pca_fel_analysis.py) via Schrodinger's
    Python in the chosen figure style, then load the results into PyMOL.
    (Blocking -- the GUI uses a non-blocking version of the same call.)

    USAGE:
        run_fel_pipeline system-out.cms, system_trj
        run_fel_pipeline system-out.cms, system_trj, style=neon_cyber
        run_fel_pipeline system-out.cms, system_trj, style=14, dpi=600, temperature=310

    :param style: style name or number 1-20 (see list_fel_styles).
    :param dpi: output resolution of the figures.
    :param temperature: simulation temperature in K (used for the FEL).
    """
    style_info = get_style(style)            # validates early, before the long run
    load_results = str(load_results).lower() not in ("0", "false", "no")
    cms_file = os.path.abspath(cms_file)
    trj_dir = os.path.abspath(trj_dir)
    output_dir = os.path.abspath(output_dir) if output_dir else os.path.dirname(cms_file)
    analysis_script = analysis_script or ANALYSIS_SCRIPT
    os.makedirs(output_dir, exist_ok=True)

    if not os.path.exists(analysis_script):
        raise RuntimeError(
            f"Could not find pca_fel_analysis.py at '{analysis_script}'. Install the plugin "
            f"from the zip (so the helper files are copied too), or pass analysis_script."
        )

    run_exe = _schrodinger_run_executable(schrodinger_home)
    args = _pipeline_args(run_exe, cms_file, trj_dir, style_info["key"], dpi, temperature,
                          analysis_script)

    print(f"[FEL Analysis] Style: {style_info['number']}. {style_info['label']}")
    print(f"[FEL Analysis] Running: {' '.join(args)}")
    result = subprocess.run(args, cwd=output_dir, capture_output=True, text=True,
                            stdin=subprocess.DEVNULL)
    print(result.stdout)
    if result.returncode != 0:
        print(result.stderr)
        raise RuntimeError(f"pca_fel_analysis.py exited with code {result.returncode}. "
                           f"See the output above for details.")

    print(f"[FEL Analysis] Pipeline finished successfully. Results in {output_dir}")
    if load_results:
        load_fel_results(output_dir)
    return output_dir


def load_fel_results(output_dir, open_images=True):
    """
    Load a completed run's results into PyMOL: every minimum_*_frame*.pdb as
    its own object (cartoon + ligand sticks), and open the FEL/PCA images.

    USAGE:
        load_fel_results /path/to/results_dir
    """
    output_dir = os.path.abspath(output_dir)
    open_images = str(open_images).lower() not in ("0", "false", "no")

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
        cmd.disable(obj_name)   # loaded but hidden -- enable the ones you want to view

    if pdb_files:
        cmd.orient()
        print(f"[FEL Analysis] Loaded {len(pdb_files)} minima structure(s) into PyMOL "
              f"(disabled by default -- enable individual objects to view them).")
    else:
        print(f"[FEL Analysis] No minimum_*_frame*.pdb files found in {output_dir}")

    if open_images:
        for img_name in RESULT_IMAGES:
            img_path = os.path.join(output_dir, img_name)
            if os.path.exists(img_path):
                _open_with_system_viewer(img_path)


# ----------------------------------------------------------------------
# Standalone rendering commands (re-render structures without re-running)
# ----------------------------------------------------------------------

def render_fel_minimum(pdb_path, out_png=None, width=1000, height=1000,
                       dpi=600, color_scheme="rainbow", bg_color="white"):
    """
    Render a single minimum-energy structure as a cartoon (+ ligand) image.

    USAGE:
        render_fel_minimum minimum_1_frame5294.pdb
        render_fel_minimum minimum_1_frame5294.pdb, bg_color=black
    """
    width, height, dpi = int(width), int(height), int(dpi)
    if out_png is None:
        out_png = os.path.splitext(pdb_path)[0] + ".png"

    obj_name = "render_" + os.path.splitext(os.path.basename(pdb_path))[0]
    obj_name = obj_name.replace("-", "_").replace(".", "_")

    cmd.load(pdb_path, obj_name)
    cmd.remove(f"{obj_name} and solvent")
    cmd.hide("everything", obj_name)
    cmd.bg_color(bg_color)
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
                            output_dir=None, color_scheme="rainbow", bg_color="white"):
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
        results.append(render_fel_minimum(pdb_path, out_png, color_scheme=color_scheme,
                                          bg_color=bg_color))

    print(f"[FEL Analysis] Rendered {len(results)} structure(s) to {output_dir}")
    return results


cmd.extend("run_fel_pipeline", run_fel_pipeline)
cmd.extend("load_fel_results", load_fel_results)
cmd.extend("render_fel_minimum", render_fel_minimum)
cmd.extend("render_fel_minima_batch", render_fel_minima_batch)
cmd.extend("list_fel_styles", list_fel_styles)
cmd.extend("fel_style_gallery", fel_style_gallery)


# ----------------------------------------------------------------------
# GUI: Step 1 choose style (with preview) -> Step 2 inputs -> Run
# ----------------------------------------------------------------------

def __init_plugin__(app=None):
    from pymol.plugins import addmenuitemqt
    addmenuitemqt("FEL Analysis", run_gui)


_dialog_ref = None   # keep the (non-modal) dialog alive


def run_gui():
    global _dialog_ref
    from pymol.Qt import QtWidgets, QtCore, QtGui

    class ScaledImage(QtWidgets.QLabel):
        """QLabel that keeps its pixmap scaled to fit, preserving aspect ratio."""

        def __init__(self, *a):
            super().__init__(*a)
            self._pix = None
            self.setAlignment(QtCore.Qt.AlignCenter)
            self.setMinimumSize(420, 340)
            self.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                               QtWidgets.QSizePolicy.Expanding)

        def set_image(self, path):
            self._pix = QtGui.QPixmap(path) if path and os.path.exists(path) else None
            if self._pix is None or self._pix.isNull():
                self._pix = None
                self.setText("No preview image found.\nClick 'Generate previews' below.")
            self._rescale()

        def _rescale(self):
            if self._pix is not None:
                self.setPixmap(self._pix.scaled(self.size(), QtCore.Qt.KeepAspectRatio,
                                                QtCore.Qt.SmoothTransformation))

        def resizeEvent(self, event):
            self._rescale()
            super().resizeEvent(event)

    settings = QtCore.QSettings("QuantaCalculus", "FEL_Analysis")
    styles = list_styles()

    dialog = QtWidgets.QDialog()
    dialog.setWindowTitle("FEL Analysis — choose a figure style, then run")
    dialog.setAttribute(QtCore.Qt.WA_DeleteOnClose, False)
    root = QtWidgets.QVBoxLayout(dialog)

    # ---------------- Step 1: style picker ----------------
    step1 = QtWidgets.QGroupBox("Step 1 — Choose a figure style (applies to all 5 figures)")
    s1 = QtWidgets.QHBoxLayout(step1)

    style_list = QtWidgets.QListWidget()
    style_list.setViewMode(QtWidgets.QListView.IconMode)
    style_list.setIconSize(QtCore.QSize(150, 125))
    style_list.setGridSize(QtCore.QSize(172, 168))
    style_list.setResizeMode(QtWidgets.QListView.Adjust)
    style_list.setMovement(QtWidgets.QListView.Static)
    style_list.setWordWrap(True)
    style_list.setUniformItemSizes(True)
    style_list.setMinimumWidth(372)
    style_list.setMaximumWidth(372)

    def _fill_list():
        style_list.clear()
        for i, (key, label, desc, cat) in enumerate(styles, start=1):
            thumb = os.path.join(PREVIEW_DIR, f"{key}_thumb.png")
            item = QtWidgets.QListWidgetItem(f"{i}. {label}")
            if os.path.exists(thumb):
                item.setIcon(QtGui.QIcon(thumb))
            item.setData(QtCore.Qt.UserRole, key)
            item.setToolTip(f"{key}  [{cat}]\n{desc}")
            style_list.addItem(item)

    _fill_list()
    s1.addWidget(style_list)

    right = QtWidgets.QVBoxLayout()
    style_title = QtWidgets.QLabel()
    f = style_title.font()
    f.setPointSize(f.pointSize() + 5)
    f.setBold(True)
    style_title.setFont(f)
    style_desc = QtWidgets.QLabel()
    style_desc.setWordWrap(True)
    preview = ScaledImage()
    preview.setStyleSheet("QLabel { background: #2b2b2b; border-radius: 6px; }")
    right.addWidget(style_title)
    right.addWidget(style_desc)
    right.addWidget(preview, 1)

    pv_buttons = QtWidgets.QHBoxLayout()
    btn_gallery = QtWidgets.QPushButton("Open all-styles gallery")
    btn_full = QtWidgets.QPushButton("Open this preview full-size")
    btn_regen = QtWidgets.QPushButton("Generate previews")
    btn_regen.setToolTip("Re-create the preview images with Schrodinger's Python "
                         "(needed only if fel_style_previews/ is missing or you edited styles).")
    for b in (btn_gallery, btn_full, btn_regen):
        pv_buttons.addWidget(b)
    pv_buttons.addStretch(1)
    right.addLayout(pv_buttons)
    s1.addLayout(right, 1)
    root.addWidget(step1, 3)

    # ---------------- Step 2: inputs ----------------
    step2 = QtWidgets.QGroupBox("Step 2 — Input trajectory and settings")
    form = QtWidgets.QFormLayout(step2)

    def _file_row(label, is_dir=False, key=None):
        line = QtWidgets.QLineEdit(settings.value(key, "") if key else "")
        btn = QtWidgets.QPushButton("Browse...")
        row = QtWidgets.QHBoxLayout()
        row.addWidget(line)
        row.addWidget(btn)
        container = QtWidgets.QWidget()
        container.setLayout(row)
        row.setContentsMargins(0, 0, 0, 0)
        form.addRow(label, container)

        def _browse():
            start = line.text().strip() or os.path.expanduser("~")
            if is_dir:
                path = QtWidgets.QFileDialog.getExistingDirectory(dialog, label, start)
            else:
                path, _ = QtWidgets.QFileDialog.getOpenFileName(
                    dialog, label, start, "Desmond CMS (*-out.cms *.cms);;All files (*)")
            if path:
                line.setText(path)

        btn.clicked.connect(_browse)
        return line

    cms_line = _file_row("-out.cms file", key="cms")
    trj_line = _file_row("_trj folder", is_dir=True, key="trj")
    out_line = _file_row("Output folder (optional)", is_dir=True, key="out")

    schrodinger_line = QtWidgets.QLineEdit(
        settings.value("schrodinger", "") or os.environ.get("SCHRODINGER", ""))
    form.addRow("Schrodinger home", schrodinger_line)

    opts = QtWidgets.QHBoxLayout()
    temp_box = QtWidgets.QDoubleSpinBox()
    temp_box.setRange(1.0, 2000.0)
    temp_box.setDecimals(1)
    temp_box.setSuffix(" K")
    temp_box.setValue(float(settings.value("temperature", 300.0)))
    dpi_box = QtWidgets.QComboBox()
    for d in ("300", "600", "1200"):
        dpi_box.addItem(d)
    dpi_box.setCurrentText(str(settings.value("dpi", "1200")))
    opts.addWidget(QtWidgets.QLabel("Temperature:"))
    opts.addWidget(temp_box)
    opts.addSpacing(20)
    opts.addWidget(QtWidgets.QLabel("Figure DPI:"))
    opts.addWidget(dpi_box)
    opts.addStretch(1)
    opts_w = QtWidgets.QWidget()
    opts.setContentsMargins(0, 0, 0, 0)
    opts_w.setLayout(opts)
    form.addRow("Options", opts_w)
    root.addWidget(step2)

    # ---------------- Log + run ----------------
    log_box = QtWidgets.QPlainTextEdit()
    log_box.setReadOnly(True)
    log_box.setMinimumHeight(110)
    log_box.setStyleSheet("QPlainTextEdit { font-family: monospace; }")
    root.addWidget(log_box, 1)

    run_row = QtWidgets.QHBoxLayout()
    run_btn = QtWidgets.QPushButton()
    run_btn.setMinimumHeight(36)
    stop_btn = QtWidgets.QPushButton("Stop")
    stop_btn.setEnabled(False)
    run_row.addWidget(run_btn, 1)
    run_row.addWidget(stop_btn)
    root.addLayout(run_row)

    state = {"key": None, "proc": None, "out_dir": None}

    def _select_style(key):
        s = get_style(key)
        state["key"] = key
        style_title.setText(f"{s['number']}. {s['label']}   ({key})")
        style_desc.setText(f"[{s['category']}]  {s['description']}")
        preview.set_image(os.path.join(PREVIEW_DIR, f"{key}_preview.png"))
        run_btn.setText(f"▶  Run full FEL analysis with style: {s['label']}")

    def _on_item_changed(current, _prev):
        if current is not None:
            _select_style(current.data(QtCore.Qt.UserRole))

    style_list.currentItemChanged.connect(_on_item_changed)

    start_key = settings.value("style", DEFAULT_STYLE)
    keys = [k for k, *_ in styles]
    style_list.setCurrentRow(keys.index(start_key) if start_key in keys else 0)

    def _log(text):
        log_box.appendPlainText(text.rstrip("\n"))
        log_box.verticalScrollBar().setValue(log_box.verticalScrollBar().maximum())

    def _set_running(running):
        run_btn.setEnabled(not running)
        stop_btn.setEnabled(running)
        btn_regen.setEnabled(not running)
        style_list.setEnabled(not running)

    def _start_process(program_args, cwd, on_done):
        proc = QtCore.QProcess(dialog)
        proc.setWorkingDirectory(cwd)
        proc.setProcessChannelMode(QtCore.QProcess.MergedChannels)
        proc.readyReadStandardOutput.connect(
            lambda: _log(bytes(proc.readAllStandardOutput()).decode(errors="replace")))

        def _finished(code, _status):
            state["proc"] = None
            _set_running(False)
            on_done(code)

        proc.finished.connect(_finished)
        state["proc"] = proc
        _set_running(True)
        proc.start(program_args[0], program_args[1:])
        if not proc.waitForStarted(10000):
            _log(f"ERROR: could not start {program_args[0]}")
            state["proc"] = None
            _set_running(False)

    def _schrodinger_or_warn():
        try:
            return _schrodinger_run_executable(schrodinger_line.text().strip() or None)
        except RuntimeError as e:
            QtWidgets.QMessageBox.critical(dialog, "FEL Analysis", str(e))
            return None

    def _run():
        cms_file = cms_line.text().strip()
        trj_dir = trj_line.text().strip()
        if not cms_file or not trj_dir:
            QtWidgets.QMessageBox.warning(dialog, "FEL Analysis",
                                          "Please select both the -out.cms file and the _trj folder.")
            return
        if not os.path.exists(ANALYSIS_SCRIPT):
            QtWidgets.QMessageBox.critical(dialog, "FEL Analysis",
                                           f"pca_fel_analysis.py not found in {PLUGIN_DIR}.\n"
                                           f"Install the plugin from the zip file.")
            return
        run_exe = _schrodinger_or_warn()
        if run_exe is None:
            return

        cms_file, trj_dir = os.path.abspath(cms_file), os.path.abspath(trj_dir)
        out_dir = out_line.text().strip()
        out_dir = os.path.abspath(out_dir) if out_dir else os.path.dirname(cms_file)
        os.makedirs(out_dir, exist_ok=True)
        state["out_dir"] = out_dir

        for k, v in (("cms", cms_line.text()), ("trj", trj_line.text()),
                     ("out", out_line.text()), ("schrodinger", schrodinger_line.text()),
                     ("temperature", temp_box.value()), ("dpi", dpi_box.currentText()),
                     ("style", state["key"])):
            settings.setValue(k, v)

        args = _pipeline_args(run_exe, cms_file, trj_dir, state["key"],
                              dpi_box.currentText(), temp_box.value(), ANALYSIS_SCRIPT)
        log_box.clear()
        _log(f"Style: {get_style(state['key'])['label']}  |  output: {out_dir}")
        _log("Running pipeline -- large trajectories can take a while; PyMOL stays usable.\n")

        def _done(code):
            if code == 0:
                _log(f"\nDone. Results in {out_dir}")
                try:
                    load_fel_results(out_dir)
                except Exception as e:  # noqa: BLE001
                    _log(f"Could not load results into PyMOL: {e}")
                QtWidgets.QMessageBox.information(dialog, "FEL Analysis",
                                                  f"Pipeline finished.\nResults in:\n{out_dir}")
            else:
                _log(f"\nERROR: pca_fel_analysis.py exited with code {code}.")
                QtWidgets.QMessageBox.critical(dialog, "FEL Analysis",
                                               f"Pipeline failed (exit code {code}). "
                                               f"See the log in this window.")

        _start_process(args, out_dir, _done)

    def _stop():
        proc = state["proc"]
        if proc is not None:
            _log("\nStopping...")
            proc.kill()

    def _regen():
        run_exe = _schrodinger_or_warn()
        if run_exe is None:
            return
        log_box.clear()
        _log("Generating style previews (about a minute)...")

        def _done(code):
            if code == 0:
                _fill_list()
                keys_now = [k for k, *_ in styles]
                style_list.setCurrentRow(keys_now.index(state["key"]) if state["key"] in keys_now else 0)
                _select_style(state["key"] or DEFAULT_STYLE)
                _log("Previews ready.")
            else:
                _log(f"Preview generation failed (exit code {code}).")

        _start_process([run_exe, "python3", PLOTS_SCRIPT, "--previews", PREVIEW_DIR],
                       PLUGIN_DIR, _done)

    run_btn.clicked.connect(_run)
    stop_btn.clicked.connect(_stop)
    btn_regen.clicked.connect(_regen)
    btn_gallery.clicked.connect(fel_style_gallery)
    btn_full.clicked.connect(lambda: _open_with_system_viewer(
        os.path.join(PREVIEW_DIR, f"{state['key']}_preview.png")))

    dialog.resize(1180, 900)
    _dialog_ref = dialog
    dialog.show()
    dialog.raise_()
