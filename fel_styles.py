"""
fel_styles.py -- 20 ready-made figure styles for the PCA / FEL pipeline.
=========================================================================

This module is PURE PYTHON DATA (no matplotlib / numpy import), so it can be
imported from anywhere -- PyMOL's own Python (for the style-picker GUI),
Schrodinger's Python (for the analysis run), or plain python3.

Every style controls ALL five output figures consistently:
    pca_fel.png, pca_fel_labeled.png, pca_fel_3d.png,
    pca_scatter.png, pca_scree.png
plus the background colour used when PyMOL renders the minima structures.

How to use
----------
    $SCHRODINGER/run python3 pca_fel_analysis.py sys-out.cms sys_trj --style neon_cyber
    $SCHRODINGER/run python3 pca_fel_analysis.py --list-styles

Adding your own style: copy any entry in STYLES, give it a new key, change
what you want. Any field you leave out falls back to DEFAULTS below.
A colormap can be a matplotlib name ("viridis") or a list of hex colours
(low free energy first) -- the list is turned into a smooth gradient.

License: MIT (see LICENSE file).
"""

DEFAULT_STYLE = "classic_rdylbu"

# ----------------------------------------------------------------------
# Every field a style can set, with its default value.
# ----------------------------------------------------------------------
DEFAULTS = {
    "label": "",
    "description": "",
    "category": "Light",

    # --- FEL colouring / rendering ---
    "fel_cmap": "RdYlBu_r",       # low free energy = first colour
    "fel_render": "contourf",     # "contourf" (smooth) or "pcolormesh" (raw bins, GROMACS-like)
    "n_levels": 60,               # contour levels (low = discrete bands, high = smooth)
    "smooth_sigma": 1.2,          # Gaussian smoothing of the FEL grid (0 = none)
    "upsample": 6,                # interpolation factor (1 = none)
    "contour_lines": True,        # thin iso-energy lines over the filled map
    "contour_color": "black",
    "contour_alpha": 0.35,
    "contour_lw": 0.25,
    "contour_every": 6,           # draw a line every N filled levels

    # --- figure / axes look ---
    "bg": "white",                # figure + axes background (unsampled FEL regions show this)
    "fg": "black",                # text, ticks, spines
    "grid": False,
    "grid_alpha": 0.15,
    "open_spines": False,         # True = remove top/right axis lines (journal look)
    "font_family": "sans-serif",  # "sans-serif" | "serif" | "monospace"
    "title_weight": "bold",
    "label_weight": "bold",
    "show_title": True,
    "font_scale": 1.0,

    # --- minima markers ---
    "minima_face": "white",
    "minima_edge": "black",
    "minima_text": "black",
    "minima_halo": "white",
    "minima_marker": "o",

    # --- PCA scatter ---
    "scatter_cmap": "plasma",
    "scatter_color": "time",      # "time" (frame index) or "density" (local point density)
    "scatter_size": 10,
    "scatter_alpha": 0.65,

    # --- 3D surface ---
    "surface_alpha": 0.97,
    "surface_edge": None,         # colour for a wireframe mesh over the surface, or None
    "floor": True,                # projected 2D contour "floor" under the surface
    "elev": 35,
    "azim": -60,
    "pane_alpha": 0.05,

    # --- scree plot ---
    "bar_color": "#4472C4",
    "line_color": "#C00000",

    # --- PyMOL rendering of minima structures ---
    "pymol_bg": "white",
}


# ----------------------------------------------------------------------
# The 20 styles. Order here = order shown in the picker.
# ----------------------------------------------------------------------
STYLES = {
    # 1
    "classic_rdylbu": {
        "label": "Classic RdYlBu",
        "description": "The original pipeline look: blue basins rising through yellow to red "
                       "barriers, soft iso-energy lines. Safe default for any paper.",
    },
    # 2
    "nature_viridis": {
        "label": "Nature Viridis",
        "description": "Perceptually uniform viridis, open axes, no title -- the clean "
                       "single-panel look Nature/Cell-family journals favour.",
        "fel_cmap": "viridis", "scatter_cmap": "viridis",
        "open_spines": True, "show_title": False,
        "contour_color": "white", "contour_alpha": 0.25,
        "bar_color": "#3b528b", "line_color": "#5ec962",
    },
    # 3
    "jet_classic": {
        "label": "Jet Classic",
        "description": "The familiar blue-to-red 'jet' landscape seen in countless MD papers, "
                       "smoothed and contoured.",
        "fel_cmap": "jet", "scatter_cmap": "jet",
        "contour_alpha": 0.45,
        "bar_color": "#1f3fbf", "line_color": "#d62828",
    },
    # 4
    "gromacs_raw": {
        "label": "GROMACS Raw Bins",
        "description": "Unsmoothed histogram bins drawn as pixels (like gmx sham + xpm2ps). "
                       "Honest view of exactly what was sampled -- good for SI figures.",
        "fel_render": "pcolormesh", "fel_cmap": "jet", "scatter_cmap": "jet",
        "smooth_sigma": 0, "upsample": 1, "contour_lines": False,
        "font_family": "monospace", "title_weight": "normal",
        "bar_color": "#1f3fbf", "line_color": "#d62828",
    },
    # 5
    "turbo_vivid": {
        "label": "Turbo Vivid",
        "description": "Google's Turbo rainbow: as vivid as jet but smoother and more "
                       "even in perceived brightness. Great for slides.",
        "fel_cmap": "turbo", "scatter_cmap": "turbo",
        "contour_alpha": 0.3,
        "bar_color": "#30123b", "line_color": "#f05b12",
    },
    # 6
    "magma_dark": {
        "label": "Magma Dark",
        "description": "Black background with glowing magma basins -- striking for "
                       "presentations and graphical abstracts.",
        "category": "Dark",
        "fel_cmap": "magma_r", "scatter_cmap": "magma",
        "bg": "#0d0d12", "fg": "#e8e8ee",
        "contour_color": "black", "contour_alpha": 0.4,
        "minima_face": "#0d0d12", "minima_edge": "#ffd166",
        "minima_text": "#ffd166", "minima_halo": "#0d0d12",
        "bar_color": "#b73779", "line_color": "#fcfdbf",
        "pymol_bg": "black",
    },
    # 7
    "inferno_density": {
        "label": "Inferno Density",
        "description": "Dark inferno landscape; the PCA scatter is coloured by local point "
                       "density instead of time, so hot-spots match the basins.",
        "category": "Dark",
        "fel_cmap": "inferno_r", "scatter_cmap": "inferno", "scatter_color": "density",
        "scatter_alpha": 0.85, "scatter_size": 8,
        "bg": "#000000", "fg": "#f2f2f2",
        "contour_color": "black", "contour_alpha": 0.5,
        "minima_face": "black", "minima_edge": "#fca50a",
        "minima_text": "#fca50a", "minima_halo": "black",
        "bar_color": "#bc3754", "line_color": "#fcffa4",
        "pymol_bg": "black",
    },
    # 8
    "ocean_deep": {
        "label": "Ocean Deep",
        "description": "Basins as deep navy trenches shoaling to sand-coloured barriers -- "
                       "a calm, elegant gradient.",
        "fel_cmap": ["#03045e", "#023e8a", "#0077b6", "#00b4d8", "#90e0ef", "#caf0f8", "#fefae0"],
        "scatter_cmap": ["#03045e", "#0077b6", "#00b4d8", "#90e0ef"],
        "contour_color": "#03045e", "contour_alpha": 0.3,
        "bar_color": "#0077b6", "line_color": "#e76f51",
    },
    # 9
    "terrain_topo": {
        "label": "Terrain Topographic",
        "description": "Reads like a topographic map: basins are lakes, barriers are "
                       "mountains, with dense contour lines.",
        "fel_cmap": "terrain", "scatter_cmap": "gist_earth",
        "contour_every": 3, "contour_alpha": 0.55, "contour_lw": 0.35,
        "bar_color": "#3d8b37", "line_color": "#8c5a2b",
    },
    # 10
    "spectral_rainbow": {
        "label": "Spectral Rainbow",
        "description": "ColorBrewer Spectral: a softer rainbow that avoids jet's harsh "
                       "bands. Colourful but still print-friendly.",
        "fel_cmap": "Spectral_r", "scatter_cmap": "Spectral",
        "contour_alpha": 0.3,
        "bar_color": "#3288bd", "line_color": "#d53e4f",
    },
    # 11
    "coolwarm_minimal": {
        "label": "Coolwarm Minimal",
        "description": "Diverging cool-warm map, no contour lines, open axes, light grid, "
                       "no title. Modern and uncluttered.",
        "fel_cmap": "coolwarm", "scatter_cmap": "coolwarm",
        "contour_lines": False, "open_spines": True, "grid": True, "grid_alpha": 0.12,
        "show_title": False, "title_weight": "normal", "label_weight": "normal",
        "bar_color": "#3b4cc0", "line_color": "#b40426",
    },
    # 12
    "cividis_accessible": {
        "label": "Cividis Accessible",
        "description": "Optimised for colour-vision deficiency (reads correctly for "
                       "deuteranopia/protanopia) and converts cleanly to greyscale.",
        "fel_cmap": "cividis", "scatter_cmap": "cividis",
        "open_spines": True, "contour_color": "white", "contour_alpha": 0.3,
        "bar_color": "#00224e", "line_color": "#d9b33a",
    },
    # 13
    "grayscale_print": {
        "label": "Grayscale Print",
        "description": "Pure greyscale with strong black contours -- for journals that "
                       "charge for colour or for black-and-white printing.",
        "fel_cmap": "gray", "scatter_cmap": "gray_r",
        "n_levels": 30, "contour_every": 2, "contour_alpha": 0.7, "contour_lw": 0.4,
        "font_family": "serif",
        "minima_face": "white", "minima_edge": "black",
        "bar_color": "#7f7f7f", "line_color": "#000000",
    },
    # 14
    "serif_journal": {
        "label": "Serif Journal",
        "description": "Times-style serif fonts, red-blue map, no title, regular-weight "
                       "labels -- matches classic ACS/Elsevier typesetting.",
        "fel_cmap": "RdBu_r", "scatter_cmap": "RdBu_r",
        "font_family": "serif", "show_title": False,
        "title_weight": "normal", "label_weight": "normal",
        "open_spines": True, "contour_alpha": 0.25,
        "bar_color": "#2166ac", "line_color": "#b2182b",
    },
    # 15
    "parula_matlab": {
        "label": "Parula",
        "description": "MATLAB's default blue-to-yellow gradient, for figures that need "
                       "to match plots made in MATLAB.",
        "fel_cmap": ["#352a87", "#0f5cdd", "#1481d6", "#06a4ca", "#2eb7a4",
                     "#87bf77", "#d1bb59", "#fec832", "#f9fb0e"],
        "scatter_cmap": ["#352a87", "#0f5cdd", "#06a4ca", "#87bf77", "#f9fb0e"],
        "contour_alpha": 0.3,
        "bar_color": "#0f5cdd", "line_color": "#d1bb59",
    },
    # 16
    "sunset_glow": {
        "label": "Sunset Glow",
        "description": "Deep violet basins warming through magenta and orange to pale "
                       "gold barriers. Eye-catching cover-art style.",
        "fel_cmap": ["#2d0a4e", "#6a1b9a", "#c2185b", "#f4511e", "#ffb300", "#fff3b0"],
        "scatter_cmap": ["#2d0a4e", "#c2185b", "#ffb300"],
        "contour_color": "#2d0a4e", "contour_alpha": 0.3,
        "minima_face": "#fff3b0", "minima_edge": "#2d0a4e",
        "bar_color": "#6a1b9a", "line_color": "#f4511e",
    },
    # 17
    "emerald_forest": {
        "label": "Emerald Forest",
        "description": "Dark emerald basins fading to pale lime -- a natural, "
                       "single-hue-family look.",
        "fel_cmap": ["#00332a", "#00684a", "#2e9e5b", "#8cc84b", "#e2ef9a", "#fbfbe9"],
        "scatter_cmap": ["#00332a", "#2e9e5b", "#e2ef9a"],
        "contour_color": "#00332a", "contour_alpha": 0.3,
        "bar_color": "#2e9e5b", "line_color": "#00332a",
    },
    # 18
    "neon_cyber": {
        "label": "Neon Cyber",
        "description": "Dark navy canvas, electric-cyan basins through magenta into the "
                       "dark, cyan wireframe on the 3D surface. Built for talks.",
        "category": "Dark",
        "fel_cmap": ["#00f5ff", "#00b3ff", "#7b2cff", "#ff00c8", "#5a0050", "#12062b"],
        "scatter_cmap": ["#00f5ff", "#7b2cff", "#ff00c8"],
        "bg": "#0a0a1f", "fg": "#c8f7ff",
        "contour_color": "#00f5ff", "contour_alpha": 0.25,
        "surface_edge": "#00f5ff",
        "grid": True, "grid_alpha": 0.08,
        "font_family": "monospace",
        "minima_face": "#0a0a1f", "minima_edge": "#00f5ff",
        "minima_text": "#00f5ff", "minima_halo": "#0a0a1f",
        "bar_color": "#7b2cff", "line_color": "#00f5ff",
        "pymol_bg": "black",
    },
    # 19
    "pastel_soft": {
        "label": "Pastel Soft",
        "description": "Muted pastel palette on a warm off-white background -- gentle, "
                       "friendly, good for posters and teaching material.",
        "fel_cmap": ["#5e81ac", "#88c0d0", "#a3be8c", "#ebcb8b", "#d08770", "#e8b4b8"],
        "scatter_cmap": ["#5e81ac", "#a3be8c", "#d08770"],
        "bg": "#fbfaf7", "fg": "#3b4252",
        "contour_color": "#3b4252", "contour_alpha": 0.2,
        "open_spines": True, "title_weight": "normal",
        "bar_color": "#88c0d0", "line_color": "#d08770",
    },
    # 20
    "discrete_bands": {
        "label": "Discrete Bands",
        "description": "Twelve sharp energy bands with a line on every boundary -- each "
                       "colour step is a readable kcal/mol interval.",
        "fel_cmap": "viridis", "scatter_cmap": "viridis",
        "n_levels": 13, "contour_every": 1, "contour_alpha": 0.6, "contour_lw": 0.5,
        "bar_color": "#31688e", "line_color": "#fde725",
    },
}


def list_styles():
    """Return [(key, label, description, category), ...] in picker order."""
    out = []
    for key, s in STYLES.items():
        merged = get_style(key)
        out.append((key, merged["label"] or key, merged["description"], merged["category"]))
    return out


def get_style(name=None):
    """
    Return a complete style dict (DEFAULTS overlaid with the chosen style).
    Accepts a style key ("neon_cyber"), its 1-based number ("18"), or None
    (-> DEFAULT_STYLE). Raises KeyError with the valid names if unknown.
    """
    if name is None or str(name).strip() == "":
        name = DEFAULT_STYLE
    name = str(name).strip()
    keys = list(STYLES.keys())
    if name.isdigit():
        idx = int(name) - 1
        if not 0 <= idx < len(keys):
            raise KeyError(f"Style number must be 1-{len(keys)}, got {name}")
        name = keys[idx]
    if name not in STYLES:
        raise KeyError(f"Unknown style '{name}'. Valid styles: {', '.join(keys)}")
    style = dict(DEFAULTS)
    style.update(STYLES[name])
    style["key"] = name
    style["number"] = keys.index(name) + 1
    if not style["label"]:
        style["label"] = name
    return style


def format_style_table():
    """Human-readable numbered list, used by --list-styles and the terminal prompt."""
    lines = []
    for i, (key, label, desc, cat) in enumerate(list_styles(), start=1):
        lines.append(f"  {i:2d}. {key:<20s} {label:<22s} [{cat}]")
    return "\n".join(lines)


if __name__ == "__main__":
    print(format_style_table())
