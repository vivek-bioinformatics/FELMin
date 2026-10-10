# FELMin

**FELMin: A PyMOL Plugin for PCA and Free Energy Landscape Analysis of Desmond MD Trajectories with Automated Minima Detection**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22909334.svg)](https://doi.org/10.5281/zenodo.22909334)
![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)

FELMin performs principal component analysis (PCA) on the protein backbone of a
Desmond MD trajectory, builds a 2D free energy landscape (FEL) on PC1/PC2,
detects the distinct free-energy minima automatically, extracts a
representative structure for each minimum, and renders publication-ready
figures — all from a single dialog inside PyMOL or one command in a terminal.

## Authors
- **Dr. Vivek Dhar Dwivedi** — Raja Shankar Shah University, Chhindwara, Madhya Pradesh, India
- **Dr. Amaresh Kumar Sahoo** — Department of Applied Sciences, Indian Institute of Information Technology Allahabad, Prayagraj 211015, India

## Features
- Backbone PCA with iterative Kabsch alignment (rigid-body motion removed)
- FEL: G(PC1, PC2) = −k<sub>B</sub>T ln(P / P<sub>max</sub>), in kcal/mol
- Automated detection of distinct minima, with the closest trajectory frame saved as `.mae` and `.pdb`
- Five figures: 2D FEL, labelled 2D FEL, 3D FEL surface, PCA scatter, scree plot
- **20 ready-made figure styles**, chosen at start-up with a live preview
- PyMOL rendering of every minimum structure; results load into PyMOL automatically
- Non-blocking run with a live log — PyMOL stays usable during long analyses
- `run_settings.json` records style, DPI, temperature and bins for reproducibility

## Requirements
- Schrödinger Suite with Desmond (the analysis runs in Schrödinger's Python via `$SCHRODINGER/run`)
- PyMOL 2.x or later (open-source or Incentive) for the plugin and structure rendering
- Input: a Desmond `*-out.cms` file and its matching `*_trj` folder

## Files
| File | Role |
|---|---|
| `__init__.py` | PyMOL plugin: style picker GUI and pipeline orchestrator |
| `pca_fel_analysis.py` | PCA, FEL, minima detection and structure extraction (Schrödinger Python) |
| `fel_styles.py` | The 20 style definitions — edit or add styles here |
| `fel_plots.py` | All plotting, shared by the real run and the style previews |
| `fel_style_previews/` | Thumbnails, 2×2 previews and `fel_style_gallery.png` |

## Installation
In PyMOL: **Plugin → Plugin Manager → Install New Plugin → Choose file…** →
select `fel_analysis_plugin.zip` → restart PyMOL → **Plugin → FEL Analysis**.

Install the **zip**, not a single `.py` file, so that all helper files are copied together.

## Usage

### PyMOL GUI
1. **Step 1 — Choose a style.** Click a thumbnail; the large 2×2 preview (2D FEL, 3D FEL, PCA scatter, scree plot) shows exactly how your figures will look.
2. **Step 2 — Inputs.** Select the `-out.cms` file, the `_trj` folder, the Schrödinger installation, the simulation temperature and the figure DPI.
3. Click **Run**. The log streams live; when the run finishes, the minima structures load into PyMOL and the figures open. Your last style and paths are remembered.

### Terminal
```bash
$SCHRODINGER/run python3 pca_fel_analysis.py --list-styles
$SCHRODINGER/run python3 pca_fel_analysis.py sys-out.cms sys_trj --style neon_cyber
$SCHRODINGER/run python3 pca_fel_analysis.py sys-out.cms sys_trj --style 14 --dpi 600 -T 310
$SCHRODINGER/run python3 pca_fel_analysis.py sys-out.cms sys_trj    # numbered style menu
```

| Option | Default | Meaning |
|---|---|---|
| `--style`, `-s` | `classic_rdylbu` | style name or number 1–20 |
| `--dpi` | 1200 | figure resolution |
| `--temperature`, `-T` | 300 | simulation temperature (K) |
| `--bins` | 60 | FEL histogram bins per axis |

### PyMOL command line
```
list_fel_styles
fel_style_gallery
run_fel_pipeline sys-out.cms, sys_trj, style=neon_cyber, dpi=600
load_fel_results /path/to/results
render_fel_minima_batch /path/to/results, bg_color=black
```

## Output
| File | Content |
|---|---|
| `pca_fel.png`, `pca_fel_labeled.png`, `pca_fel_3d.png` | Free energy landscape (2D, with minima, 3D) |
| `pca_scatter.png`, `pca_scree.png` | PC1/PC2 projection and variance explained |
| `minimum_N_frameX.mae / .pdb` | Representative structure for each minimum |
| `structure_N_frameX.png` | PyMOL rendering of each minimum |
| `minima_summary.json` | Energy, frame, PC1/PC2 of each minimum |
| `pc1_pc2_projections.csv` | Per-frame PC1/PC2 values |
| `run_settings.json` | Style, DPI, temperature and bins used |

## The 20 figure styles
| # | Style | # | Style |
|---|---|---|---|
| 1 | `classic_rdylbu` | 11 | `coolwarm_minimal` |
| 2 | `nature_viridis` | 12 | `cividis_accessible` (colour-blind safe) |
| 3 | `jet_classic` | 13 | `grayscale_print` |
| 4 | `gromacs_raw` (unsmoothed bins) | 14 | `serif_journal` |
| 5 | `turbo_vivid` | 15 | `parula_matlab` |
| 6 | `magma_dark` | 16 | `sunset_glow` |
| 7 | `inferno_density` | 17 | `emerald_forest` |
| 8 | `ocean_deep` | 18 | `neon_cyber` |
| 9 | `terrain_topo` | 19 | `pastel_soft` |
| 10 | `spectral_rainbow` | 20 | `discrete_bands` |

See `fel_style_previews/fel_style_gallery.png` for all styles side by side.
Each style sets the colormap, smooth or raw-bin rendering, contour bands and
iso-lines, light or dark background, fonts, axis style, minima markers,
scatter colouring (time or density), 3D view, scree colours and the PyMOL
background for structure renders.

### Adding your own style
Copy any entry in `STYLES` in `fel_styles.py`, give it a new key and change the
fields you want (unset fields fall back to `DEFAULTS`). Then regenerate the previews:
```bash
$SCHRODINGER/run python3 fel_plots.py --previews fel_style_previews
```
or click **Generate previews** in the plugin.

## Troubleshooting
- **`pymol` not found during structure rendering:** run with
  `PYMOL_EXECUTABLE=/full/path/to/pymol $SCHRODINGER/run python3 ...`
- **Rendering fails on a headless server:** prefix the command with `xvfb-run -a`.
- **No minima found:** loosen `MINIMA_ENERGY_CUTOFF` or `MIN_BASIN_SEPARATION` at the top of `pca_fel_analysis.py`.

## Citation
If you use FELMin in your work, please cite:

> Dwivedi, V. D., & Sahoo, A. K. (2026). *FELMin: A PyMOL Plugin for PCA and Free Energy Landscape Analysis of Desmond MD Trajectories with Automated Minima Detection*. Zenodo. https://doi.org/10.5281/zenodo.22909334

```bibtex
@software{dwivedi_sahoo_2026_felmin,
  author    = {Dwivedi, Vivek Dhar and Sahoo, Amaresh Kumar},
  title     = {{FELMin: A PyMOL Plugin for PCA and Free Energy Landscape
                Analysis of Desmond MD Trajectories with Automated Minima
                Detection}},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.22909334},
  url       = {https://doi.org/10.5281/zenodo.22909334}
}
```

## License
MIT License — see `LICENSE`. FELMin calls the Schrödinger Suite Python API and
PyMOL, which are separate third-party software not covered by this license and
are subject to their own licence terms.
