#!/usr/bin/env python
"""
On-demand RGB composite of a lens's wavebands with every band's contaminant mask outlined,
for judging masks by COLOUR without opening the GUI.

When vetting whether a masked (or unmasked) feature is a contaminant or a lensed image, the
discriminant is its colour relative to the deflector (AGENTS.md, *Masks*: "compare colours as
a ratio against the deflector's own"; the multiple images of one source share a colour).
This writes one PNG per lens into diagnostics/masks/ -- the folder for "is the mask boundary
right?" figures -- so that question can be answered from a picture.

THE COMPOSITE IS "DEFLECTOR-GREY". Each band has its sky (median) subtracted and is divided
by its mean inside a small central aperture, so the deflector has unit flux in every band;
the bands then share ONE arcsinh stretch. The deflector renders neutral grey and anything
bluer or redder than it tints blue or red. Units (WFPC2 DN/s vs ACS e/s), depths and zero
points cancel. It is a finding aid, not a calibrated colour image. The helper is a copy of
`autogalaxy.gui.display_util.rgb_from_bands` (PyAutoGalaxy branch feature/scribbler-rgb);
the project venv pins autogalaxy 2026.7.29.2, which predates it.

BANDS ARE CHOSEN BY WAVELENGTH: the longest available is red, the shortest blue, the one
between them green. With two bands the green channel is the mean of red and blue; with four
or more the middle band is used (override with --bands). One-band lenses are skipped.

ALL BANDS ARE REPROJECTED ONTO THE REFERENCE BAND'S WCS GRID (`detect_arcs.regrid`,
bilinear) before combining -- the optical bands are not index-aligned (sub-pixel offsets) and
f160W is a different grid altogether (0.06"/px, 200 px). Each band's `{prefix}_mask.fits`
is reprojected the same way (`make_masks.reproject_mask_bool`, nearest) and drawn as an
outline in that band's colour; the arc mask, where present, as a dashed white outline.

THE FIGURE IS IN THE GUI FRAME: pixel axes, origin 'upper' (row 0 at the top), so it matches
what make_masks.py shows -- East LEFT, North DOWN -- with a 1" scale bar. Give positions
read off it as (column, row).

Usage:
    uv run python scripts/make_rgb.py --sample gallery --lens J0029+2544
    uv run python scripts/make_rgb.py --sample slacs_gold --lens J0008-0004   # two bands
    uv run python scripts/make_rgb.py --sample sl2s_dinos_a                   # whole sample
    uv run python scripts/make_rgb.py --lens J0959+0206 --sample sl2s_dinos_a --bands f814W,f606W,f475X
"""

import argparse
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

ws_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mast_target_names
import stamp_paths
import make_masks
import detect_arcs

# Central wavelength (nm) of every filter in the stamp tree; unknown filters are dropped
# with a note. f606W_v2 is the shorter second visit of the split-visit lens J0822+2652 (AGENTS.md).
WAVELENGTH_NM = {
    'f275W': 275, 'f438W': 438, 'f475X': 475, 'f555W': 555, 'f606W': 606,
    'f606W_v2': 606, 'f600LP': 750, 'f814W': 814, 'f160W': 1600,
}

# Outline colour per band, fixed by wavelength rank so the same band is always the same
# colour across lenses: red-ish for the reddest band through to blue-ish for the bluest.
OUTLINE_COLOURS = ['#ff5555', '#ffcc00', '#55ff55', '#55ddff', '#dd88ff']


# ----------------------------------------------------------------------------- composite

def _arcsinh_stretch_shared(channels, clip_percentiles=(0.5, 99.5)):
    """arcsinh stretch of several images with ONE set of parameters (sky, noise sigma and
    clip limits measured on the pooled pixels), so equal inputs map to equal outputs.
    Copy of autogalaxy.gui.display_util._arcsinh_stretch_shared."""
    arrays = [np.asarray(c, dtype=float) for c in channels]
    pooled = np.concatenate([a[np.isfinite(a)] for a in arrays])
    if not pooled.size:
        return [np.zeros_like(a) for a in arrays]
    median = np.median(pooled)
    sigma = 1.4826 * np.median(np.abs(pooled - median))
    if not sigma > 0:
        sigma = np.std(pooled)
    if not sigma > 0:
        return [np.zeros_like(a) for a in arrays]
    stretched = [np.arcsinh((np.where(np.isfinite(a), a, median) - median) / sigma)
                 for a in arrays]
    pooled = np.arcsinh((pooled - median) / sigma)
    lo, hi = np.percentile(pooled, clip_percentiles)
    if not hi > lo:
        lo, hi = pooled.min(), pooled.max()
    if not hi > lo:
        return [np.zeros_like(a) for a in arrays]
    out = []
    for a, st in zip(arrays, stretched):
        o = np.clip((st - lo) / (hi - lo), 0.0, 1.0)
        o[~np.isfinite(a)] = 0.0
        out.append(o)
    return out


def rgb_from_bands(red, green, blue, aperture_radius_px=None, clip_percentiles=(0.5, 99.5)):
    """Deflector-grey RGB of two or three bands already on one grid: (ny, nx, 3) in [0, 1].
    `green=None` synthesises green as the mean of red and blue. Copy of
    autogalaxy.gui.display_util.rgb_from_bands (deflector_grey=True)."""
    bands = {'red': red, 'green': green, 'blue': blue}
    arrays = {k: np.asarray(b, dtype=float) for k, b in bands.items() if b is not None}
    shapes = {k: a.shape for k, a in arrays.items()}
    if any(len(s) != 2 for s in shapes.values()) or len(set(shapes.values())) != 1:
        raise ValueError(f'every band must be 2D and on the same grid, got {shapes}')
    ny, nx = next(iter(shapes.values()))
    if aperture_radius_px is None:
        aperture_radius_px = max(3, int(round(0.05 * min(ny, nx))))
    yy, xx = np.mgrid[0:ny, 0:nx]
    aperture = np.hypot(yy - (ny - 1) / 2.0, xx - (nx - 1) / 2.0) <= aperture_radius_px
    normalised = []
    for a in arrays.values():
        finite = np.isfinite(a)
        sky = np.median(a[finite]) if finite.any() else 0.0
        a = a - sky
        inside = aperture & finite
        mean = float(a[inside].mean()) if inside.any() else 0.0
        if mean > 0:
            a = a / mean
        normalised.append(a)
    stretched = dict(zip(arrays.keys(), _arcsinh_stretch_shared(normalised, clip_percentiles)))
    if 'green' not in stretched:
        stretched['green'] = 0.5 * (stretched['red'] + stretched['blue'])
    return np.clip(np.dstack([stretched['red'], stretched['green'], stretched['blue']]), 0.0, 1.0)


# ----------------------------------------------------------------------------- bands

def pick_rgb_bands(filts, requested=None):
    """(red, green, blue) band names by wavelength, green None for two bands; None if fewer
    than two bands have a known wavelength. `requested` is an explicit 'R,G,B' (G may be
    'none')."""
    if requested:
        parts = [p.strip() for p in requested.split(',')]
        if len(parts) != 3:
            raise SystemExit(f'--bands wants R,G,B (G may be none), got {requested!r}')
        red, green, blue = parts
        green = None if green.lower() == 'none' else green
        missing = [b for b in (red, green, blue) if b is not None and b not in filts]
        if missing:
            raise SystemExit(f'--bands names band(s) this lens does not have: {missing}')
        return red, green, blue
    known = sorted((f for f in filts if f in WAVELENGTH_NM), key=lambda f: WAVELENGTH_NM[f])
    unknown = [f for f in filts if f not in WAVELENGTH_NM]
    if unknown:
        print(f'  ignoring band(s) with no wavelength in WAVELENGTH_NM: {unknown}')
    if len(known) < 2:
        return None
    red, blue = known[-1], known[0]
    green = known[len(known) // 2] if len(known) >= 3 else None
    return red, green, blue


STAMP_SIZE = None   # main() sets it from --size; None: each lens's modelling size


def load_band(stamp_dir, drizzle_pass):
    """sci / header / contaminant mask / arc region for one band, via detect_arcs.load_band
    when the band has a PSF stamp (it centres the kernel) and directly otherwise."""
    prefix = stamp_paths.find_prefix(stamp_dir, drizzle_pass, STAMP_SIZE)
    if prefix is None:
        return None
    if os.path.exists(os.path.join(stamp_dir, f'psf.fits')):
        band = detect_arcs.load_band(stamp_dir, prefix)
        sci, hdr, cmask = band['sci'], band['hdr'], band['cmask']
    else:
        sci_path = stamp_paths.find_stamp(stamp_dir, prefix, 'sci')
        sci = np.nan_to_num(fits.getdata(sci_path).astype(float))
        hdr = fits.getheader(sci_path)
        mask_path = os.path.join(stamp_dir, f'{prefix}_mask.fits')
        cmask = (np.asarray(fits.getdata(mask_path), dtype=bool) if os.path.exists(mask_path)
                 else np.zeros(sci.shape, dtype=bool))
    arc_path = os.path.join(stamp_dir, f'{prefix}_mask_arcs.fits')
    arcs = (~np.asarray(fits.getdata(arc_path), dtype=bool) if os.path.exists(arc_path)
            else None)   # the file stores the inverse: True = NOT arc
    return dict(prefix=prefix, sci=sci, hdr=hdr, cmask=cmask, arcs=arcs)


# ----------------------------------------------------------------------------- figure

def make_rgb_png(lens, bands, rgb_names, ref_filt, out_path, with_arcs=True, dpi=110):
    """Reproject every band onto `ref_filt`'s grid, build the composite, outline each band's
    mask and write the PNG. Returns the path."""
    ref_hdr = bands[ref_filt]['hdr']
    ref_wcs = WCS(ref_hdr).celestial
    ref_shape = bands[ref_filt]['sci'].shape
    ps_ref = make_masks.pixel_scale_from_header(ref_hdr)

    on_ref = {f: detect_arcs.regrid(b['sci'], b['hdr'], ref_hdr, order=1) for f, b in bands.items()}
    red, green, blue = rgb_names
    rgb = rgb_from_bands(on_ref[red], None if green is None else on_ref[green], on_ref[blue])

    fig, ax = plt.subplots(figsize=(8.5, 8.5))
    ax.imshow(rgb, origin='upper', interpolation='none')

    by_wavelength = sorted(bands, key=lambda f: -WAVELENGTH_NM.get(f, 0))
    handles = []
    for i, filt in enumerate(by_wavelength):
        b = bands[filt]
        if not b['cmask'].any():
            handles.append(Line2D([], [], color='none', label=f'{filt}: no mask on disk'))
            continue
        mask = make_masks.reproject_mask_bool(b['cmask'], WCS(b['hdr']).celestial, ref_wcs, ref_shape)
        colour = OUTLINE_COLOURS[i % len(OUTLINE_COLOURS)]
        # No extent / origin on contour: Z[0, 0] lands at data (0, 0), where imshow draws
        # pixel [0, 0] -- the same alignment reasoning as the Scribbler's proposal outline.
        ax.contour(mask.astype(float), levels=[0.5], colors=[colour], linewidths=1.4)
        role = {red: ' (R)', green: ' (G)', blue: ' (B)'}.get(filt, '')
        handles.append(Line2D([], [], color=colour, label=f'{filt} mask{role}  {int(mask.sum())} px'))
        if with_arcs and b['arcs'] is not None:
            arcs = make_masks.reproject_mask_bool(b['arcs'], WCS(b['hdr']).celestial, ref_wcs, ref_shape)
            ax.contour(arcs.astype(float), levels=[0.5], colors='w', linewidths=1.0, linestyles='dashed')
            handles.append(Line2D([], [], color='w', linestyle='dashed', label=f'{filt} arc mask'))

    ax.add_artist(AnchoredSizeBar(ax.transData, 1.0 / ps_ref, '1"', 'lower left', pad=0.4,
                                  color='w', frameon=False, size_vertical=1.0))
    ax.legend(handles=handles, loc='upper right', fontsize=8, framealpha=0.6)
    ax.set_xlabel('column (px)')
    ax.set_ylabel('row (px)   [row 0 at top: East LEFT, North DOWN]')
    green_label = green if green is not None else 'mean(R,B)'
    ax.set_title(f'{lens}   R={red}  G={green_label}  B={blue}   ref {ref_filt} @ {ps_ref:.3f}"/px   '
                 f'deflector-grey arcsinh', fontsize=10)
    ax.grid(color='w', alpha=0.15, linewidth=0.5)
    fig.tight_layout()
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)
    return out_path


# ----------------------------------------------------------------------------- driver

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--sample', default=mast_target_names.DEFAULT_SAMPLE)
    p.add_argument('--lens', default=None, help='one lens (default: every lens in the sample)')
    p.add_argument('--pass', dest='drizzle_pass', choices=['auto', 'cr', 'nocrrej'], default='auto')
    p.add_argument('--size', type=float, default=None,
                   help='stamp size in arcsec (default: each lens\'s modelling size, its stamp_size in '
                        'info/samples/lens_samples.json, else 12"). Only modelling-size stamps are '
                        'tracked, so another size is skipped unless --include-ignored')
    p.add_argument('--ref-band', default=None,
                   help='band whose grid the composite is on (default: make_masks priority, f475X then f814W)')
    p.add_argument('--bands', default=None, help='explicit R,G,B band names; G may be "none"')
    p.add_argument('--out-dir', default=os.path.join(ws_path, 'diagnostics', 'masks'))
    p.add_argument('--no-arcs', action='store_true', help='do not outline arc masks')
    p.add_argument('--include-ignored', action='store_true')
    args = p.parse_args()
    global STAMP_SIZE
    STAMP_SIZE = args.size

    root = stamp_paths.stamps_root(ws_path)
    per_lens = {}
    for lens, filt, stamp_dir in make_masks.discover_targets(root, args.sample, args.lens,
                                                                include_ignored=args.include_ignored):
        per_lens.setdefault(lens, {})[filt] = stamp_dir
    if not per_lens:
        raise SystemExit(f'no stamps under {root}/{args.sample}' + (f'/{args.lens}' if args.lens else ''))

    for lens, filt_dirs in per_lens.items():
        if len(filt_dirs) < 2:
            print(f'{lens}: only {", ".join(filt_dirs)}, no colour -- skipped')
            continue
        bands = {}
        for filt, stamp_dir in filt_dirs.items():
            band = load_band(stamp_dir, args.drizzle_pass)
            if band is not None:
                bands[filt] = band
        rgb_names = pick_rgb_bands(list(bands), args.bands)
        if rgb_names is None:
            print(f'{lens}: fewer than two bands with a known wavelength ({list(bands)}) -- skipped')
            continue
        ref_filt = make_masks.pick_display_filt(list(bands), args.ref_band)
        if ref_filt is None:
            print(f'{lens}: --ref-band {args.ref_band} not among {list(bands)} -- skipped')
            continue
        out_path = os.path.join(args.out_dir, f'{lens}_rgb_masks.png')
        make_rgb_png(lens, bands, rgb_names, ref_filt, out_path, with_arcs=not args.no_arcs)
        red, green, blue = rgb_names
        print(f'{lens}: R={red} G={green or "mean(R,B)"} B={blue} on {ref_filt} -> '
              f'{os.path.relpath(out_path, ws_path)}')


if __name__ == '__main__':
    main()
