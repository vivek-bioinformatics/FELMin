# PCA / Free Energy Landscape Analysis for Desmond Trajectories

Scripts for running principal component analysis on a Desmond MD trajectory, building a free energy landscape from the PC1/PC2 projection, locating the distinct energy minima, and pulling out the structures at those minima for rendering in PyMOL.

Two files here:

- `pca_fel_analysis.py` — does the actual work: reads the trajectory, aligns frames, runs PCA, builds the FEL, finds minima, extracts structures, renders them with PyMOL, and produces the plots.
- `fel_analysis_plugin.py` — a PyMOL plugin that wraps the above into a one-click workflow from inside PyMOL's GUI.

## Why two separate environments

This is worth explaining up front, because it trips people up. The trajectory reading and PCA depend on Schrodinger's own Python packages (`schrodinger.application.desmond.packages`, `schrodinger.trajectory.prody`), which only exist inside Schrodinger's bundled Python — not inside PyMOL's Python, and not in a normal system Python either. So `pca_fel_analysis.py` has to be run through `$SCHRODINGER/run`.

PyMOL, meanwhile, is used purely for rendering the extracted structures as cartoon images. The main script calls it as a separate subprocess rather than trying to import it, since mixing the two Python environments in one process causes library conflicts (mismatched `libstdc++`, wrong Python picked up by PyMOL's own launcher, that sort of thing — I ran into all of these while building this and the fixes are baked into the script already).

The plugin (`fel_analysis_plugin.py`) exists for people who'd rather not touch a terminal at all. It launches the main script through `$SCHRODINGER/run` from inside PyMOL, waits for it to finish, then loads the results straight into your PyMOL session.

## Requirements

- A Schrodinger Suite license, with `$SCHRODINGER` pointing at your install directory.
- PyMOL, reachable either on your normal shell `PATH` or via an explicit path (see below — `$SCHRODINGER/run`'s environment doesn't inherit your regular shell PATH, which is a separate gotcha from the Python environment issue above).
- Python packages that ship with Schrodinger's own Python already cover most of what's needed (numpy, matplotlib). If `scipy` isn't present, the smoothing step falls back to a cruder blur automatically rather than failing outright.

Neither Schrodinger nor PyMOL is included in this repository — you need your own valid installations of both.

## Running it from the command line

```bash
$SCHRODINGER/run python3 pca_fel_analysis.py system-out.cms system_trj
```

If PyMOL isn't found automatically (common when `$SCHRODINGER/run`'s environment doesn't see your normal PATH), point at it directly:

```bash
PYMOL_EXECUTABLE=/path/to/pymol $SCHRODINGER/run python3 pca_fel_analysis.py system-out.cms system_trj
```

## Running it from inside PyMOL

Install `fel_analysis_plugin.py` through PyMOL's Plugin Manager (Plugin → Plugin Manager → Install New Plugin), keeping it in the same folder as `pca_fel_analysis.py`. Restart PyMOL, then use Plugin → FEL Analysis to point at your `.cms` file and trajectory folder and run the whole thing from a small dialog.

The same commands are also available headlessly if you'd rather script it:

```
run fel_analysis_plugin.py
run_fel_pipeline system-out.cms, system_trj
```

## What it does, roughly in order

1. Loads the trajectory and pulls out backbone atom coordinates for every frame.
2. Aligns every frame to a converged mean structure (Kabsch superposition), removing rigid-body translation and rotation. Skipping this step is the single most common reason people get a useless, noise-dominated FEL — the rigid-body motion completely swamps the actual internal conformational fluctuations if it's left in.
3. Runs PCA on the aligned coordinates and projects the trajectory onto PC1/PC2.
4. Bins the projection into a 2D histogram and converts it to a free energy surface via ΔG = −k_BT ln(P/P_max).
5. Smooths and upsamples the resulting grid so the plots look like a continuous surface rather than blocky histogram bins.
6. Finds the distinct local minima in the landscape (not just the single global minimum — the actual separate basins), matches each one to its nearest real trajectory frame, and extracts that frame's full structure.
7. Hands each extracted structure to PyMOL for rendering: cartoon representation, rainbow N-to-C coloring, ligand shown as sticks if present, solvent stripped out.
8. Produces the final set of plots — 2D and 3D free energy landscapes (with minima labeled), a PCA scree plot, and a raw PC1/PC2 scatter colored by simulation time.

## Output files

Running the pipeline on a trajectory produces, in the working directory:

- `pc1_pc2_projections.csv` — raw PC1/PC2 value for every frame
- `pca_fel.png`, `pca_fel_3d.png`, `pca_fel_labeled.png` — the free energy landscape, in three views
- `pca_scatter.png` — PC1/PC2 scatter colored by frame index
- `pca_scree.png` — variance explained per component
- `minimum_N_frameXXXX.mae` / `.pdb` — the extracted structure at each minimum
- `structure_N_frameXXXX.png` — the PyMOL-rendered image of each

## A few settings worth knowing about

Near the top of `pca_fel_analysis.py`:

- `TEMPERATURE` — set this to whatever temperature your simulation actually ran at. It's used directly in the free energy calculation.
- `N_BINS` — how finely the PC1/PC2 space gets binned before smoothing. Too high relative to how many frames you have and the histogram gets sparse and noisy.
- `MAX_MINIMA`, `MINIMA_ENERGY_CUTOFF`, `MIN_BASIN_SEPARATION` — control how many distinct minima get reported and how far apart they need to be to count as separate basins rather than the same one picked up twice.
- `DPI` — resolution for all the output images.

## License

MIT — see `LICENSE`. That covers this code only. Schrodinger's software and PyMOL are separate third-party tools with their own licensing, not included here.
