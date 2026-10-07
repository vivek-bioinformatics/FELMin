#!/usr/bin/env python3
"""
fel_minima_figure.py
--------------------
Multi-panel publication figure: several FEL plots (e.g. 4 systems), each with
its minima structures placed around it and connected by arrows from the
numbered minimum markers (1, 2, 3 ...) on the FEL.

What the script does automatically
  * finds the circled minimum markers on each FEL and reads their number
    (template matching; you can give coordinates manually instead)
  * moves the colorbar to the outer edge so arrows never run over it
  * routes every arrow along a path that avoids axis labels / tick text
  * puts each structure in a labelled box ("Minimum k", "Frame N")
  * builds a grid of panels (A), (B), (C), (D) at an exact physical width
    and saves PNG, TIFF (LZW) and PDF at the requested DPI (default 600)

Requirements:  pip install numpy scipy pillow matplotlib

HOW TO USE (no editing needed)
  FEL_Minima/                      <- main folder
  |-- fel_minima_figure.py         <- this script
  |-- A/   FEL image + its minima structure images
  |-- B/   ...
  |-- C/   ...
  |-- D/   ...
  Run:  python fel_minima_figure.py      (or double-click it on Windows)
  Output goes into the main folder: FEL_minima_figure.png / .tif / .pdf

  * Sub-folders are used in name order -> panels (A), (B), (C), (D).
    Name them A, B, C, D (or 1_Apo, 2_Cpd1 ...) so the order is right.
  * FEL image  = the file whose name contains "fel" (else the largest image).
  * Structures = all other images in the folder. The minimum number is read
    from the name: structure_3..., min3, minimum_3, 3.png, 3_xxx.png ...
    ("frame644" in the name is shown as "Frame 644"). If no number is found,
    files are numbered 1, 2, 3 ... in name order.
  * Optional markers.txt in a folder ("1 3979 5656" per line = minimum, x, y
    pixel in the FEL image) if a circled number is not detected automatically.
  * Folders starting with "_" or "." are ignored.
"""

import glob
import math
import os
import re
import sys
import traceback

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy import ndimage
from matplotlib import font_manager

# =============================== USER SETTINGS ===============================
# Main folder = folder of this script (or give a path:  python fel_minima_figure.py D:\FEL_Minima)
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
IMG_EXT  = (".png", ".jpg", ".jpeg", ".tif", ".tiff", ".bmp")

TITLE_FROM_FOLDER = False   # True -> folder name ("Apo_HisA") printed next to (A)
NCOLS          = 2          # panels per row (2 -> 2x2 grid for 4 systems; 1 system -> 1)
FIG_WIDTH_MM   = 180        # final printed width (Elsevier/ACS double column ~174-190 mm)
DPI            = 600
OUT_BASENAME   = "FEL_minima_figure"   # writes .png, .tif, .pdf

FONT_FAMILY    = "DejaVu Sans"   # "Arial" if installed; DejaVu ships with matplotlib
PANEL_LETTERS  = True            # (A), (B), (C) ...
SHOW_TITLES    = True            # used only when TITLE_FROM_FOLDER = True
LETTER_PT      = 10
TITLE_PT       = 8
MIN_LABEL_PT   = 5.5             # "Minimum k"
FRAME_PT       = 5               # "Frame N"   (set None to hide)
ONE_LINE_LABEL = True            # "Minimum k  |  Frame N" on one line
ARROW_PT       = 0.7             # arrow line width in points
ARROW_COLOR    = (20, 20, 20)
BOX_COLOR      = (70, 70, 70)
MOVE_COLORBAR  = True            # move colorbar to the outer side of each panel
SAVE_DEBUG     = True            # saves _marker_checks/<folder>.png in the main folder
PAUSE_AT_END   = True            # keep the window open when double-clicked (Windows)
# ============================================================================

Image.MAX_IMAGE_PIXELS = None


# ----------------------------------------------------------------- helpers --
def font_file(bold=True, family=None):
    fp = font_manager.FontProperties(family=family or FONT_FAMILY,
                                     weight="bold" if bold else "normal")
    return font_manager.findfont(fp, fallback_to_default=True)


def get_font(px, bold=True):
    return ImageFont.truetype(font_file(bold), max(int(round(px)), 6))


def flatten_white(im):
    im = im.convert("RGBA")
    bg = Image.new("RGBA", im.size, (255, 255, 255, 255))
    bg.alpha_composite(im)
    return bg.convert("RGB")


def content_box(im, thr=245, pad=0):
    a = np.asarray(im)
    ys, xs = np.where(a.min(axis=2) < thr)
    if len(xs) == 0:
        return (0, 0, im.width, im.height)
    return (max(xs.min() - pad, 0), max(ys.min() - pad, 0),
            min(xs.max() + pad + 1, im.width), min(ys.max() + pad + 1, im.height))


def parse_structure_name(path):
    name = os.path.splitext(os.path.basename(path))[0]
    m = (re.search(r"(?:structure|struct|minimum|minima|min|state|conf)[ _\-]*(\d+)", name, re.I)
         or re.match(r"(\d+)(?:$|[ _\-.])", name))
    f = re.search(r"frame[ _\-]*(\d+)", name, re.I)
    return (int(m.group(1)) if m else None), (f.group(1) if f else None)


def natural_key(s):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", s)]


def discover_systems(base):
    """Every sub-folder with images = one panel. Returns list of cfg dicts."""
    systems = []
    subs = sorted((d for d in os.listdir(base)
                   if os.path.isdir(os.path.join(base, d)) and not d.startswith(("_", "."))),
                  key=natural_key)
    for d in subs:
        folder = os.path.join(base, d)
        imgs = sorted((os.path.join(folder, f) for f in os.listdir(folder)
                       if f.lower().endswith(IMG_EXT) and "markers_check" not in f.lower()),
                      key=lambda p: natural_key(os.path.basename(p)))
        if len(imgs) < 2:
            print(f"  skip folder '{d}': needs a FEL image + at least one structure")
            continue
        named = [p for p in imgs if "fel" in os.path.basename(p).lower()]
        if len(named) == 1:
            fel = named[0]
        else:   # largest image (FEL plots are much bigger than structure renders)
            cand = named or imgs
            fel = max(cand, key=lambda p: Image.open(p).size[0] * Image.open(p).size[1])
        structs = [p for p in imgs if p != fel]

        nums = [parse_structure_name(p)[0] for p in structs]
        if None in nums or len(set(nums)) != len(nums):
            print(f"  [i] '{d}': minimum number not found in all file names -> "
                  f"numbered 1..{len(structs)} in name order")
            smap = {i + 1: (p, parse_structure_name(p)[1]) for i, p in enumerate(structs)}
        else:
            smap = {n: (p, parse_structure_name(p)[1]) for n, p in zip(nums, structs)}

        markers = None
        mt = os.path.join(folder, "markers.txt")
        if os.path.isfile(mt):
            markers = {}
            for line in open(mt):
                v = line.replace(",", " ").split()
                if len(v) >= 3 and v[0].isdigit():
                    markers[int(v[0])] = (float(v[1]), float(v[2]))
            print(f"  '{d}': using markers.txt for minima {sorted(markers)}")

        systems.append(dict(name=d, title=d if TITLE_FROM_FOLDER else "",
                            fel=fel, structures=smap, markers=markers))
    return systems


# ------------------------------------------------------- colorbar handling --
def split_colorbar(img):
    """Return (plot_without_colorbar, colorbar_image or None)."""
    a = np.asarray(img).astype(int)
    H, W, _ = a.shape
    dark = a.max(axis=2) < 70
    x0 = int(W * 0.6)
    cols = np.where(dark[:, x0:].sum(axis=0) > 0.3 * H)[0] + x0
    if len(cols) < 2:
        return img, None
    # group adjacent columns into lines
    groups = np.split(cols, np.where(np.diff(cols) > 3)[0] + 1)
    centers = [g.mean() for g in groups]
    for i in range(len(centers) - 1):
        L, R = groups[i][0], groups[i + 1][-1]
        if not (0.008 * W < R - L < 0.09 * W):
            continue
        inner = a[int(H * 0.3):int(H * 0.7), groups[i][-1] + 3:groups[i + 1][0] - 2]
        if inner.size == 0:
            continue
        sat = (inner.max(axis=2) - inner.min(axis=2)).mean()
        if sat < 40:            # not a coloured bar
            continue
        # cut at the white gap just left of the bar
        nonwhite = (a.min(axis=2) < 240).any(axis=0)
        x = L - 1
        while x > 0 and nonwhite[x]:
            x -= 1
        gap_r = x
        while x > 0 and not nonwhite[x]:
            x -= 1
        if gap_r - x < 0.004 * W:      # no clean white gap -> leave it
            return img, None
        cut = (x + gap_r) // 2
        plot = img.crop((0, 0, cut, H))
        cb = img.crop((cut, 0, W, H))
        cb = cb.crop(content_box(cb, pad=10))
        return plot, cb
    return img, None


# --------------------------------------------------------- marker finding --
_TEMPL = {}


def _digit_templates():
    if _TEMPL:
        return _TEMPL
    for bold in (True, False):
        f = ImageFont.truetype(font_file(bold, "DejaVu Sans"), 120)
        for d in "123456789":
            im = Image.new("L", (200, 200), 255)
            ImageDraw.Draw(im).text((40, 20), d, font=f, fill=0)
            m = np.asarray(im) < 128
            ys, xs = np.where(m)
            m = m[ys.min():ys.max() + 1, xs.min():xs.max() + 1]
            _TEMPL.setdefault(int(d), []).append(_norm(m))
    return _TEMPL


def _norm(mask, n=32):
    im = Image.fromarray((mask * 255).astype(np.uint8)).resize((n, n), Image.BILINEAR)
    return np.asarray(im) > 127


def _read_digit(mask):
    ys, xs = np.where(mask)
    if len(xs) < 10:
        return None, 0
    m = _norm(mask[ys.min():ys.max() + 1, xs.min():xs.max() + 1])
    best, score = None, -1
    for d, temps in _digit_templates().items():
        for t in temps:
            s = (m & t).sum() / max((m | t).sum(), 1)
            if s > score:
                best, score = d, s
    return best, score


def detect_markers(img):
    """Find white circles with a dark ring and a digit inside -> {num: (x, y)}."""
    a = np.asarray(img).astype(int)
    H, W, _ = a.shape
    white = a.min(axis=2) > 230
    dark = a.max(axis=2) < 110
    lab, _ = ndimage.label(white)
    lo, hi = 0.006 * W, 0.035 * W
    hits = []
    for i, sl in enumerate(ndimage.find_objects(lab)):
        if sl is None:
            continue
        h, w = sl[0].stop - sl[0].start, sl[1].stop - sl[1].start
        if not (lo < h < hi and lo < w < hi and 0.85 < h / w < 1.18):
            continue
        comp = lab[sl] == i + 1
        if comp.sum() < 0.45 * h * w:
            continue
        cy, cx = (sl[0].start + sl[0].stop) / 2, (sl[1].start + sl[1].stop) / 2
        r = (h + w) / 4
        ang = np.linspace(0, 2 * np.pi, 48, endpoint=False)
        ring = 0
        for f in (1.05, 1.12, 1.2):
            yy = np.clip((cy + f * r * np.sin(ang)).astype(int), 0, H - 1)
            xx = np.clip((cx + f * r * np.cos(ang)).astype(int), 0, W - 1)
            ring = max(ring, dark[yy, xx].mean())
        if ring < 0.6:
            continue
        Y, X = np.mgrid[sl[0], sl[1]]
        inside = (Y - cy) ** 2 + (X - cx) ** 2 < (0.85 * r) ** 2
        digit = inside & ~comp & dark[sl]
        num, score = _read_digit(digit)
        if num is not None and score > 0.35:
            hits.append((num, score, (cx, cy), r))
    out, radius = {}, {}
    for num, score, xy, r in sorted(hits, key=lambda t: -t[1]):
        if num not in out:
            out[num], radius[num] = xy, r
    return out, (np.median(list(radius.values())) if radius else 0.012 * W)


# ------------------------------------------------------------ arrow paths --
def route(dark_small, f, markers, edge_x, H, r, all_markers=None):
    """Choose a via point on the plot edge for each marker so that the first
    arrow segment crosses as little text/axes as possible; keeps order so
    arrows on one side never cross each other."""
    Hs, Ws = dark_small.shape
    vias, last = {}, -1e9
    minsep = 0.035 * H
    for k, (mx, my) in sorted(markers.items(), key=lambda t: t[1][1]):
        best = None
        for ty in np.arange(max(my - 0.15 * H, 0.03 * H),
                            min(my + 0.15 * H, 0.97 * H), H / 300):
            if ty < last + minsep:
                continue
            L = math.hypot(edge_x - mx, ty - my)
            t0 = (1.4 * r) / max(L, 1)
            ts = np.linspace(t0, 1, max(int(L / f), 2))
            xs = ((mx + (edge_x - mx) * ts) / f).astype(int).clip(0, Ws - 1)
            ys = ((my + (ty - my) * ts) / f).astype(int).clip(0, Hs - 1)
            cost = dark_small[ys, xs].sum() + 0.08 * abs(ty - my) / f
            for kk, (ox, oy) in (all_markers or {}).items():   # keep clear of other markers
                if kk == k:
                    continue
                dmin = np.min(np.hypot(xs * f - ox, ys * f - oy))
                if dmin < 2.2 * r:
                    cost += 2000 * (2.2 * r - dmin) / r
            if best is None or cost < best[0]:
                best = (cost, ty)
        ty = best[1] if best else max(my, last + minsep)
        vias[k] = ty
        last = ty
    return vias


def draw_arrow(d, pts, color, w, head):
    (x0, y0), (x1, y1) = pts[-2], pts[-1]
    ang = math.atan2(y1 - y0, x1 - x0)
    bx, by = x1 - head * math.cos(ang), y1 - head * math.sin(ang)
    path = list(pts[:-1]) + [(bx, by)]
    halo = max(int(w * 1.8), w + 4)
    d.line(path, fill="white", width=halo, joint="curve")
    d.line(path, fill=color, width=w, joint="curve")

    def tri(hl, ext):
        hx, hy = x1 + ext * math.cos(ang), y1 + ext * math.sin(ang)
        cx, cy = x1 - hl * math.cos(ang), y1 - hl * math.sin(ang)
        pw = hl * 0.5
        return [(hx, hy), (cx + pw * math.sin(ang), cy - pw * math.cos(ang)),
                (cx - pw * math.sin(ang), cy + pw * math.cos(ang))]
    d.polygon(tri(head + halo * 0.6, halo * 0.3), fill="white")
    d.polygon(tri(head, 0), fill=color)


# ------------------------------------------------------------ one panel ----
def build_block(cfg, block_w_final):
    raw = flatten_white(Image.open(cfg["fel"]))
    bx0, by0, bx1, by1 = content_box(raw, pad=20)
    fel = raw.crop((bx0, by0, bx1, by1))

    plot, cbar = split_colorbar(fel) if MOVE_COLORBAR else (fel, None)
    H = plot.height

    if cfg.get("markers"):
        markers = {int(k): (x - bx0, y - by0) for k, (x, y) in cfg["markers"].items()}
        r = 0.012 * plot.width
    else:
        markers, r = detect_markers(plot)

    # structures
    structs = cfg["structures"]

    missing = sorted(set(structs) - set(markers))
    if missing:
        print(f"  [!] {cfg['name']}: no circled number found on the FEL for minima {missing} "
              f"-> skipped (add them to {cfg['name']}/markers.txt)")
    ks = sorted(set(structs) & set(markers))
    extra = sorted(set(markers) - set(structs))
    if extra:
        print(f"  [i] '{cfg['name']}': markers {extra} on the FEL have no structure file - ignored")
    print(f"  FEL = {os.path.basename(cfg['fel'])}; minima {ks} at "
          + ", ".join(f"{k}:({markers[k][0]+bx0:.0f},{markers[k][1]+by0:.0f})" for k in ks))

    if SAVE_DEBUG:
        dbg = plot.copy()
        dd = ImageDraw.Draw(dbg)
        fdb = get_font(r * 2.2)
        for k in ks:
            x, y = markers[k]
            dd.ellipse([x - 2 * r, y - 2 * r, x + 2 * r, y + 2 * r], outline="red", width=max(int(r / 5), 2))
            dd.text((x + 2.2 * r, y - 2.6 * r), str(k), font=fdb, fill="red")
        s = 1600 / dbg.width
        ddir = os.path.join(BASE_DIR, "_marker_checks")
        os.makedirs(ddir, exist_ok=True)
        dbg.resize((1600, int(dbg.height * s))).save(os.path.join(ddir, cfg["name"] + ".png"))

    # split minima into left / right columns, balance counts
    cx_plot = plot.width / 2
    left = [k for k in ks if markers[k][0] < 0.45 * plot.width]
    right = [k for k in ks if k not in left]
    while len(right) > math.ceil(len(ks) / 2):
        m = min(right, key=lambda k: markers[k][0]); right.remove(m); left.append(m)
    while len(left) > math.ceil(len(ks) / 2):
        m = max(left, key=lambda k: markers[k][0]); left.remove(m); right.append(m)

    # geometry (native pixels)
    SW = int(0.40 * H)                 # side column width
    G = int(0.035 * H)                 # gutter between plot and boxes
    nL, nR = len(left), len(right)
    cbw = cbar.width if cbar is not None else 0
    W = (SW + G if nL else 0) + plot.width + (G + SW if nR else 0) + (int(0.03 * H) + cbw if cbar is not None else 0)
    s = block_w_final / W               # final scale factor
    px = lambda pt: pt * DPI / 72 / s   # points -> native pixels

    canvas = Image.new("RGB", (W, H), "white")
    FX = SW + G if nL else 0
    canvas.paste(plot, (FX, 0))
    if cbar is not None:
        cy = max((H - cbar.height) // 2, 0)
        canvas.paste(cbar, (W - cbw, cy))
    d = ImageDraw.Draw(canvas)

    f_min = get_font(px(MIN_LABEL_PT), True)
    f_frm = get_font(px(FRAME_PT), False) if FRAME_PT else None
    lw = max(int(px(ARROW_PT)), 2)
    head = lw * 5.5
    box_lw = max(int(px(0.5)), 2)

    def column(keys, x_left, edge_x, side):
        if not keys:
            return
        keys = sorted(keys, key=lambda k: markers[k][1])
        gap = int(0.03 * H)
        ph = min(int(SW * 0.82), (H - (len(keys) - 1) * gap) // len(keys))
        total = len(keys) * ph + (len(keys) - 1) * gap
        y = (H - total) // 2
        dark = (np.asarray(plot).astype(int).max(axis=2) < 100)
        f = 4
        dark_small = ndimage.binary_dilation(dark[::f, ::f], iterations=4)
        mk = {k: markers[k] for k in keys}
        vias = route(dark_small, f, mk, edge_x - FX, H, r, all_markers={kk: markers[kk] for kk in ks})
        for k in keys:
            path, frame = structs[k]
            x0, y0 = x_left, y
            d.rounded_rectangle([x0, y0, x0 + SW, y0 + ph], radius=int(0.04 * ph),
                                outline=BOX_COLOR, width=box_lw, fill="white")
            title = f"Minimum {k}"
            ty = y0 + int(0.04 * ph)
            sub = f"  |  Frame {frame}"
            one_line = (ONE_LINE_LABEL and f_frm and frame and
                        d.textlength(title, font=f_min) + d.textlength(sub, font=f_frm) < 0.92 * SW)
            if one_line:
                w1 = d.textlength(title, font=f_min); w2 = d.textlength(sub, font=f_frm)
                xs = x0 + (SW - w1 - w2) / 2
                d.text((xs, ty), title, font=f_min, fill="black")
                d.text((xs + w1, ty + (f_min.size - f_frm.size) * 0.75), sub, font=f_frm, fill=(80, 80, 80))
                ty += int(f_min.size * 1.3)
            else:
                tw = d.textlength(title, font=f_min)
                d.text((x0 + (SW - tw) / 2, ty), title, font=f_min, fill="black")
                ty += int(f_min.size * 1.2)
            if f_frm and frame and not one_line:
                sub = f"Frame {frame}"
                sw_ = d.textlength(sub, font=f_frm)
                d.text((x0 + (SW - sw_) / 2, ty), sub, font=f_frm, fill=(80, 80, 80))
                ty += int(f_frm.size * 1.25)
            st = flatten_white(Image.open(path))
            st = st.crop(content_box(st, pad=10))
            bw, bh = SW * 0.92, (y0 + ph - ty) * 0.95
            sc = min(bw / st.width, bh / st.height)
            st = st.resize((max(int(st.width * sc), 1), max(int(st.height * sc), 1)), Image.LANCZOS)
            canvas.paste(st, (int(x0 + (SW - st.width) / 2), int(ty + (y0 + ph - ty - st.height) / 2)))

            mx, my = markers[k]; mx += FX
            vx, vy = edge_x, vias[k]
            tx = x0 - lw if side == "R" else x0 + SW + lw
            tyy = y0 + ph / 2
            a = math.atan2(vy - my, vx - mx)
            start = (mx + 1.2 * r * math.cos(a), my + 1.2 * r * math.sin(a))
            draw_arrow(d, [start, (vx, vy), (tx, tyy)], ARROW_COLOR, lw, head)
            y += ph + gap

    column(right, FX + plot.width + G, FX + plot.width - int(0.005 * plot.width), "R")
    column(left, 0, FX + int(0.005 * plot.width), "L")
    return canvas, s


# --------------------------------------------------------------- assemble --
def main():
    global BASE_DIR, NCOLS
    if len(sys.argv) > 1:
        BASE_DIR = os.path.abspath(sys.argv[1])
    print("Main folder:", BASE_DIR)
    systems = discover_systems(BASE_DIR)
    if not systems:
        raise SystemExit("No sub-folders with images found next to the script.")
    print("Panels:", ", ".join(f"({chr(65 + i)})={c['name']}" for i, c in enumerate(systems)))
    NCOLS = min(NCOLS, len(systems))

    fig_w = int(round(FIG_WIDTH_MM / 25.4 * DPI))
    margin = int(0.01 * fig_w)
    hgap = int(0.025 * fig_w)
    vgap = int(0.02 * fig_w)
    bw = (fig_w - 2 * margin - (NCOLS - 1) * hgap) // NCOLS

    blocks = []
    for cfg in systems:
        print(f"\nBuilding panel from folder '{cfg['name']}'")
        blk, s = build_block(cfg, bw)
        blk = blk.resize((bw, int(round(blk.height * bw / blk.width))), Image.LANCZOS)
        blocks.append((cfg, blk))

    pt = lambda p: p * DPI / 72
    head_h = int(pt(max(LETTER_PT, TITLE_PT)) * 1.5) if (PANEL_LETTERS or SHOW_TITLES) else 0
    rows = [blocks[i:i + NCOLS] for i in range(0, len(blocks), NCOLS)]
    row_h = [head_h + max(b.height for _, b in r) for r in rows]
    fig_h = 2 * margin + sum(row_h) + (len(rows) - 1) * vgap

    fig = Image.new("RGB", (fig_w, fig_h), "white")
    d = ImageDraw.Draw(fig)
    f_let = get_font(pt(LETTER_PT), True)
    f_tit = get_font(pt(TITLE_PT), True)
    y = margin
    n = 0
    for r, rh in zip(rows, row_h):
        x = margin
        for cfg, blk in r:
            tx = x
            if PANEL_LETTERS:
                lab = f"({chr(65 + n)})"
                d.text((tx, y), lab, font=f_let, fill="black")
                tx += d.textlength(lab, font=f_let) + pt(3)
            if SHOW_TITLES and cfg.get("title"):
                d.text((tx, y + (f_let.size - f_tit.size) * 0.8), cfg["title"],
                       font=f_tit, fill="black")
            fig.paste(blk, (x, y + head_h + (rh - head_h - blk.height) // 2))
            x += bw + hgap
            n += 1
        y += rh + vgap

    out = os.path.join(BASE_DIR, OUT_BASENAME)
    fig.save(out + ".png", dpi=(DPI, DPI))
    fig.save(out + ".tif", dpi=(DPI, DPI), compression="tiff_lzw")
    fig.save(out + ".pdf", resolution=DPI)
    print(f"\nSaved {out}.png/.tif/.pdf  -  {fig_w} x {fig_h} px "
          f"= {FIG_WIDTH_MM:.0f} x {fig_h / DPI * 25.4:.0f} mm at {DPI} dpi")


if __name__ == "__main__":
    try:
        main()
    except SystemExit as e:
        print(e)
    except Exception:
        traceback.print_exc()
    if PAUSE_AT_END and os.name == "nt":
        input("\nPress Enter to close...")
