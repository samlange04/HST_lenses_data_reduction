#!/usr/bin/env python
"""
Interactive GUI mask-making, ONE draw per band, seeded by a REVIEWED mask from another band.

Merged 2026-11-28 across the HST and JWST repos: the file is identical in both. The repo
profile (lens_samples / mast_target_names) is the import seam below; the stamp layout is each
repo's own `stamp_paths`. On screen North is DOWN and East is LEFT in every panel, the one
convention of both repos (display_look.ORIGIN; the JWST side drew North-up until the merge).

For each lens in a sample, launches PyAutoLens's `Scribbler` GUI over the lens's
best-available stamp -- the same tool as
autolens_workspace:scripts/imaging/data_preparation/gui/mask.py -- so you can scribble the
areas to REMOVE from the fit (painted = excluded, everything unpainted = kept; this is the
opposite polarity to the autolens_workspace GUI, which scribbles the region to keep), then
writes the result as stamp<S>_mask.fits alongside that stamp's sci/noise/psf products.

ONE DRAW PER BAND. Each run draws on ONE band -- the lens's PRIMARY band (its
`best_band_visual` / `best_band_arc_snr`, see choose_primary_band.py; the fixed order below is
the fallback) unless --filt forces it -- and writes the mask for THAT BAND
ONLY. Run again with --filt <band> to do the next band; the skip check is per-band, so a
sample sweep resumes cleanly. Per-band drawing is the rule because what a mask should exclude
is NOT in practice band-independent: a contaminant can be bright in one filter and absent in
another, and each band's own depth and PSF change where the sensible boundary falls.

THE NEXT BAND STARTS FROM A REVIEWED PROPOSAL, NOT A BLANK CANVAS (--propose-from, default
'auto'). A mask already drawn for another band of this lens is reprojected onto the band
about to be drawn and OUTLINED for approval over the usual two panels:

    LEFT   radial-subtracted -- where the arcs are visible at all, so you can see what the
           inherited mask may be clipping;
    RIGHT  as-observed -- the contaminant's real extent and the galaxy envelope in THIS band.

The outline is all the review needs; a filled-in "mask applied" panel was tried and dropped,
since it showed nothing the boundary does not already carry while hiding the pixels being
judged. A multi-band lens also gets a COLOUR-INDEX panel in the same window (log10 of the
bluest/reddest band ratio against the deflector's own, PSF-matched; blue = bluer than the
deflector): the lensed-source vs contaminant discriminant made visible. You then edit the
proposal with two brushes --
'1' WHITE adds to the mask, '2' BLACK erases from it, scribbling on whichever panel you like --
and on closing the GUI choose to [a]pply proposal+edits, keep only what you [d]rew (rejecting
the proposal outright), or [s]kip the band. Source: the band's own existing mask when
--force-redrawing it (so refining a mask does not mean redrawing it); for a split-visit band
the primary visit's mask in the SAME filter; else the lens's primary band's mask, then the
highest-priority other band that has one. --propose-from <band> forces the source,
--propose-from none draws on a blank canvas.

--broadcast (mutually exclusive with the above) is the old unreviewed route: write the one
drawn mask to every band of the lens, sight-unseen. Both share the same transfer -- a rigorous
per-pixel WCS reprojection (nearest-neighbour, source sci-header WCS -> target sci-header WCS),
NOT an array-index copy: bands that share a plate scale share a grid to <1px so it is nearly
identity, but a genuinely different grid (HST's f160W at 0.06"/px vs the 0.05" optical bands;
JWST's SW/LW 2:1 nesting) needs the reprojection to place the mask on the same sky.

Defaults to each lens's MODELLING-SIZE stamp (stamp_paths.modelling_size) -- the size the
hand-drawn masks are meant for and the size tracked in git (see .gitignore) precisely because
it carries them, which no script can regenerate. (The stamp pixels themselves ship through
Zenodo in the HST repo; a fresh clone fetches them with scripts/fetch_data.py, which the tools
print when they find no stamps.) Pass --size 20 to mask another size instead.

One stamp dir per band, so one mask per band, written beside the sci it was drawn on.
On the HST repo's `main` that sci is the standard drizzle; the bad-column-filled variant lives
only on its `bcfill` branch, under a tagged stamp name in the same dir (stamp_paths.py), and
shares the grid to <1e-3 px, so a mask drawn on either is valid on both. Either way there is
nothing to choose and no second copy of the mask to keep in step. A (lens, filt) with a mask
already on disk is SKIPPED; --force redraws it.

Every mask written queues that band's `{prefix}_dataset.png` QC subplot and, at the END of the
run, the dataset subplots and the changed lenses' band overviews are re-rendered in ONE
DETACHED background process (rerender_overviews; log in data/run_logs/). Inline rendering put
an isophote fit between one GUI and the next, and -- worse -- a pyplot figure after a Scribbler
has run flashes up as a blank window (the Scribbler switches pyplot to interactive mode).
--no-dataset-subplot skips the queue. The renders are best-effort QC: the mask is the product
and is never lost to a plotting error.

This is a manual, one-lens-at-a-time tool (not a batch driver): each lens blocks on its own
Tk window until you press Esc. A lens that already has a mask is skipped so a run resumes
across a sample; --force redraws.

Usage:
    uv run python scripts/make_masks.py --sample slacs                # best band per lens
    uv run python scripts/make_masks.py --lens J0008-0004                  # one lens, best band
    uv run python scripts/make_masks.py --lens J0008-0004 --filt f606W     # force the draw band
    uv run python scripts/make_masks.py --lens J0008-0004 --filt f160W --force  # redo one band
    uv run python scripts/make_masks.py --sample slacs --filt f555W    # review f814W's
    uv run python scripts/make_masks.py --lens J0008-0004 --propose-from none  # blank canvas
    uv run python scripts/make_masks.py --lens J0008-0004 --broadcast --propose-from none
    uv run python scripts/make_masks.py --sample slacs --force        # redraw everything
"""

import argparse
import glob
import os
import sys
import warnings

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.visualization import (AsinhStretch, LogStretch, SqrtStretch,
                                   LinearStretch, ManualInterval)

# autonerves prints a workspace-version-mismatch UserWarning on import in this repo (it's
# not an autolens_workspace checkout) -- harmless, silence it so it doesn't repeat per lens.
warnings.filterwarnings('ignore', category=UserWarning, module='autonerves')

import autolens as al
import autolens.plot as aplt
# The display helpers this repo used to carry its own copies of. Re-exported from here on
# purpose: make_positions.py and make_arc_masks.py call `make_masks.radial_median_subtract`
# / `make_masks.mask_boundary`, and both are now the upstream functions (same algorithm, same
# `exclude=` keyword, same inner-boundary / border_value=1 semantics as the copies they replace).
from autogalaxy.gui.display_util import radial_median_subtract, mask_boundary  # noqa: F401

ws_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import display_look
import info_json
import stamp_paths
try:
    import lens_samples as mast_target_names   # the JWST repo's profile module
except ImportError:
    import mast_target_names                   # the HST repo's

MASKS_JSON = os.path.join(ws_path, 'info', 'provenance', 'lens_masks.json')


#: blank columns between al.Scribbler's two side-by-side display panels (its `panel_gap`).
_PANEL_GAP = 6


# The GUI itself is al.Scribbler, unmodified. Everything this file once had to build around
# it -- a subclass guarding out-of-axes mouse motion and forcing brush redraws, x1.4 brush
# resizing with a 1 px floor, a title smuggled in through a class attribute, a hand-hstacked
# side-by-side display with a NaN gutter, a 1-px proposal outline burned into the pixels, and
# the fold-back of strokes from the composite onto the band grid -- is native in autolens
# >= 2026.10.7 (`proposal=`, `subtract_radial=`, `side_by_side=`, `panel_gap=`, `positions=`,
# `mask_overlay=`, `title=`, `mask_from()`), so draw_mask_gui just passes the arguments through.

_STRETCHES = {
    'asinh': AsinhStretch,   # bright cores + faint structure at once; --asinh-a tunes it
    'log': LogStretch,
    'sqrt': SqrtStretch,
    'linear': LinearStretch,
}


def stretched_display(data_native, pixel_scales, stretch='asinh',
                      vmin_percent=5.0, vmax_percent=99.5, asinh_a=0.1,
                      exclude=None, vmax_value=None):
    """Return an al.Array2D holding a contrast-stretched COPY of `data_native`, for
    display in make_positions.py's Clicker and in the QC PNGs only.

    NOT used by the Scribbler GUIs any more: al.Scribbler takes display_look.SCRIBBLER as its
    `cmap` (each panel normalised with its own PercentileInterval(99) + AsinhStretch(0.1) on
    inferno), so draw_mask_gui hands it the raw base array. This percentile + astropy
    AsinhStretch(a) version is kept for the Clicker, which imshows its own axes, and for the
    QC PNGs, so those products render exactly as before.

    Positions are built from click coordinates, never from pixel values, so remapping the
    displayed values cannot change the saved product -- this is purely so bright cores and
    faint/problematic structure are both visible while clicking. Values are clipped to the
    [vmin_percent, vmax_percent] percentiles then mapped through `stretch` onto [0, 1] (a
    plain linear imshow of that IS the stretch).

    `exclude` (bool array, True = hide) blanks pixels to the display floor AND drops them
    from the percentile statistics -- the second half is the point: a deflector core left
    in the stats pins vmax far above the arcs, so blanking it also re-scales everything
    else. `vmax_value` sets the upper clip in the base array's own units (S/N, or flux)
    instead of by percentile, so anything brighter saturates at the top of the colourmap.
    """
    vals = np.asarray(data_native, dtype=float)
    keep = np.isfinite(vals)
    if exclude is not None:
        keep &= ~np.asarray(exclude, dtype=bool)
    finite = vals[keep]
    if finite.size == 0:
        return al.Array2D.no_mask(values=np.nan_to_num(vals), pixel_scales=pixel_scales).native
    vmin, vmax = np.percentile(finite, [vmin_percent, vmax_percent])
    if vmax_value is not None:
        vmax = float(vmax_value)
    if not vmax > vmin:                      # flat/degenerate frame -- avoid /0 in interval
        vmax = vmin + (abs(vmin) or 1.0)
    st = AsinhStretch(a=asinh_a) if stretch == 'asinh' else _STRETCHES[stretch]()
    disp = np.nan_to_num(st(ManualInterval(vmin=vmin, vmax=vmax)(vals)))  # clip -> stretch
    if exclude is not None:
        disp[np.asarray(exclude, dtype=bool)] = 0.0    # hidden region -> colourmap floor
    return al.Array2D.no_mask(values=disp, pixel_scales=pixel_scales).native


# `radial_median_subtract(data_native, exclude=None)` and `mask_boundary(m)` used to be
# defined here; they are now imported from `autogalaxy.gui.display_util` at the top of this
# file (identical algorithms -- the upstream versions were lifted from these) and are still
# reachable as `make_masks.radial_median_subtract` / `make_masks.mask_boundary` for
# make_positions.py's Clicker panels. The Scribbler does its own radial subtraction
# (`subtract_radial=True`) and proposal outlining (`proposal=`) internally.


def central_disc(shape, pixel_scales, radius_arcsec):
    """Bool array, True inside a disc of `radius_arcsec` at the stamp centre.

    Stamps are recentred on the deflector (`make_stamps.py --center-band`), so the stamp
    centre IS the lens galaxy -- this is the region to hide when the central light drowns
    the arcs. Returns None for a non-positive radius (feature off).
    """
    if not radius_arcsec or radius_arcsec <= 0:
        return None
    n_y, n_x = shape
    yy, xx = np.mgrid[0:n_y, 0:n_x]
    cy, cx = (n_y - 1) / 2.0, (n_x - 1) / 2.0
    r_px = float(radius_arcsec) / float(pixel_scales)
    return ((yy - cy) ** 2 + (xx - cx) ** 2) <= r_px ** 2


def pixel_scale_from_header(hdr):
    """arcsec/pixel from a stamp FITS header's WCS -- Array2D.from_fits wants a plain
    float, not the WCS itself."""
    return float(proj_plane_pixel_scales(WCS(hdr).celestial)[0] * 3600.0)


# The stamp size the GUI tools work on. Every size shares one band directory, so the size
# picks the prefix (stamp12_*, stamp20_*) rather than a tree; each tool's main() sets this
# from its --size before discovering targets.
STAMP_SIZE = None   # None: each lens's modelling size


def find_prefix(stamp_dir, drizzle_pass='auto'):
    """Delegates to stamp_paths.find_prefix at STAMP_SIZE: 'stamp12' for this band (or
    'stamp12_nocrrej' for a lone comparison stamp), whichever reduction (standard / bcfill
    / ...) its stamp is from -- the stamp names are variant-tagged, so this is the one
    place that knows how to find them."""
    return stamp_paths.find_prefix(stamp_dir, drizzle_pass, STAMP_SIZE)


def stamp_pass_of(stamp_dir, prefix):
    """'cr' or 'nocrrej' for the sci stamp under `prefix` (stamp_paths.stamp_pass)."""
    return stamp_paths.stamp_pass(stamp_paths.find_stamp(stamp_dir, prefix, 'sci'))


def gitignored(paths):
    """Delegates to stamp_paths.gitignored -- see there for why it lives in that module."""
    return stamp_paths.gitignored(paths, ws_path)


def discover_targets(root, sample, lens=None, filt=None, include_ignored=False):
    """Yield (lens, filt, stamp_dir) per band with a stamp under `root` -- the dir both
    scribbled on and written to -- sorted for a reproducible run order.

    GITIGNORED PRODUCT DIRS ARE SKIPPED (2026-09-22, and this is the guard the repo was
    missing; test updated 2026-11-28). A band whose hand-made product would be gitignored
    is, by construction, not a science product -- the three 420 s SLACS SNAP f814W
    diagnostics (J1134+6027, J1403+0006, J1538+5817) are the reason the rule exists.
    Before this, the sweep tools happily offered them, and a hand-drawn mask or a set of
    hand-marked positions would be written into a directory git refuses to store:
    non-regenerable work that no clone would ever see. It already happened -- J1538+5817's
    image positions were MARKED on its SNAP f814W band and broadcast from there to f606W
    across a 688 mas (13.8 px) band misregistration, and both had to be deleted on
    2026-09-22.

    The test is the prospective `{prefix}_mask.fits` path, not the sci stamp (2026-11-28):
    the sci FITS is regenerable -- in the HST repo it ships through Zenodo and is ALWAYS
    gitignored since the 2026-11-27 untracking, so testing it would skip everything -- while
    a mask written into a gitignored directory (a non-modelling stamp size) is exactly the
    product that must be tracked.

    `include_ignored=True` opts back in, for deliberately working on a diagnostic product.
    The skip is announced, never silent, so a short target list is explicable.
    """
    found, sci_of, prefixes = {}, {}, {}
    # Lenses come from the roster, so every stamp layout is found; the band dirs are each
    # repo's own (stamp_paths.lens_dir).
    names = [lens] if lens else mast_target_names.lenses(sample)
    for this_lens in sorted(names):
        lens_dir = stamp_paths.lens_dir(ws_path, 'stamps', sample, this_lens)
        for stamp_dir in sorted(glob.glob(os.path.join(lens_dir, filt or '*'))):
            if not os.path.isdir(stamp_dir):
                continue
            prefix = find_prefix(stamp_dir, 'auto')
            if prefix is None:
                continue
            this_filt = os.path.basename(stamp_dir)
            found[(this_lens, this_filt)] = stamp_dir
            prefixes[(this_lens, this_filt)] = prefix
            sci_of[(this_lens, this_filt)] = stamp_paths.find_stamp(stamp_dir, prefix, 'sci')
    if not include_ignored:
        # Test the prospective hand-made product (the contaminant-mask path), NOT the sci
        # stamp: the sci FITS is regenerable -- in the HST repo it ships through Zenodo and is
        # ALWAYS gitignored since 2026-11-27, so testing it skips everything -- while a mask
        # or positions file written into a gitignored dir is non-regenerable work no clone
        # would see (in the JWST repo a non-modelling stamp size is exactly that). The guard
        # the 2026-09-22 rule meant, kept honest about what is actually tracked.
        prospect = {k: os.path.join(v, f'{prefixes[k]}_mask.fits') for k, v in found.items()}
        ignored = gitignored(list(prospect.values()))
        skipped = sorted(k for k, v in prospect.items() if os.path.abspath(v) in ignored)
        for key in skipped:
            del found[key]
        if skipped:
            print(f"  skipping {len(skipped)} gitignored band(s), not science products: "
                  + ', '.join(f'{l}/{f}' for l, f in skipped)
                  + "  (--include-ignored to use them anyway)")
    if not found:
        # A fresh clone of the HST repo carries no data products (they live on Zenodo, not in
        # git). Say so with the next step, instead of an empty run that reads like a bug.
        what = lens or f"sample {sample}"
        print(f"  no stamps found for {what} under {root}")
        hint = stamp_paths.fetch_hint(sample, filt)
        if hint:
            print(f"  if {what} has data products, fetch them from Zenodo with:")
            print(f"      {hint}")
    for key in sorted(found):
        yield key[0], key[1], found[key]


def load_display_base(display_dir, prefix, display, pixel_scales):
    """Return (base_array, source_label): the array to stretch and scribble on.

    display='snr' returns signal/noise (the recommended base -- the mask should track
    statistical significance, and this pipeline's noise map is spatially non-uniform, so
    brightness != significance); 'sci' returns the raw signal. S/N is set to 0 where the
    noise map is <=0 or non-finite (zero-coverage edges), which correctly suppresses
    those regions rather than blowing them up. Falls back to signal, with a note, if the
    noise map is missing.
    """
    sci = al.Array2D.from_fits(
        file_path=stamp_paths.find_stamp(display_dir, prefix, 'sci'),
        pixel_scales=pixel_scales).native
    if display == 'sci':
        return sci, 'signal'
    noise_path = stamp_paths.find_stamp(display_dir, prefix, 'noise')
    if noise_path is None:
        print(f"  NOTE: no {prefix}[_<variant>]_noise.fits -- displaying signal instead of S/N")
        return sci, 'signal'
    noise = np.asarray(al.Array2D.from_fits(
        file_path=noise_path, pixel_scales=pixel_scales).native, dtype=float)
    snr = np.zeros_like(noise)
    good = np.isfinite(noise) & (noise > 0)
    snr[good] = np.asarray(sci, dtype=float)[good] / noise[good]
    return al.Array2D.no_mask(values=snr, pixel_scales=pixel_scales).native, 'S/N'


# ── Lupton RGB colour composite (positions/masks GUIs, 2026-11-17) ──────────────────
# Colour is the discriminant between lensed images (blue) and the deflector/contaminants
# (red), so both GUIs show a true-colour panel. The bands are resampled onto the finest
# band's grid through their WCS (map_coordinates, bilinear), sky-neutralised per channel,
# and combined with astropy's Lupton asinh. A 2-band lens uses G = sqrt(R*B) (Lupton's
# own 2-band recipe); >3 bands average the middle ones into G. Display only, everywhere.

def band_pivot_wavelength(filt):
    """Pivot wavelength (um) of a band key, split-visit suffixes stripped (f606W_visit2 ->
    f606W), for wavelength ordering. The table is make_band_overview's FILTER_UM (one table for
    both repos' filters; locally imported: keep the module graph one-directional)."""
    import re as _re
    from make_band_overview import FILTER_UM
    key = _re.sub(r'_visit\d+$', '', filt)
    if key not in FILTER_UM:
        raise KeyError(f'unknown filter {filt!r}; add its pivot wavelength to '
                       f'make_band_overview.FILTER_UM')
    return FILTER_UM[key]


def lupton_rgb(bands, stretch_scale=0.15, Q=8.0, balance=False):
    """(rgb float array (ny,nx,3), base_index) from a list of band dicts
    {'filt', 'data', 'hdr', 'ps'} (2+ entries; returns (None, None) for fewer).

    `bands` are resampled onto the finest pixel scale's grid through the WCS. `data` should
    be FLUX (sci), or whatever per-band quantity the caller wants coloured (the positions
    GUI colours both the raw sci and the isophote-subtracted panels). Per-channel the border
    median is subtracted so the sky stays neutral; the asinh stretch is
    `stretch_scale` x the 99.9th percentile of the red channel. NaNs (zero-coverage edges,
    the odd NaN pixel in a stamp) are zeroed -- one NaN pixel otherwise blacks out the whole
    Lupton image (found on J0841+3824 f160W, 2026-11-17). `balance=True` first divides each
    channel by its own 99.9th percentile -- for the residual colour panel, where a band whose
    isophote fit failed carries its full unsubtracted flux and would otherwise own the red
    channel; for flux (the default) the true ratios are the point and are kept.
    """
    from astropy.visualization import make_lupton_rgb
    from scipy.ndimage import map_coordinates

    usable = [b for b in bands if b is not None]
    if len(usable) < 2:
        return None, None
    order = sorted(range(len(usable)),
                   key=lambda i: band_pivot_wavelength(usable[i]['filt']))
    usable = [usable[i] for i in order]              # blue -> red
    base = min(range(len(usable)), key=lambda i: usable[i]['ps'])   # finest grid
    tgt = usable[base]
    ny, nx = np.asarray(tgt['data']).shape
    yy, xx = np.mgrid[:ny, :nx]
    wcs_t = WCS(tgt['hdr'])
    sky = wcs_t.pixel_to_world(xx, yy)

    channels = []
    for b in usable:
        d = np.nan_to_num(np.asarray(b['data'], float))
        if b is tgt:
            channels.append(d)
            continue
        wcs_b = WCS(b['hdr'])
        bx, by = wcs_b.world_to_pixel(sky)
        channels.append(np.nan_to_num(map_coordinates(d, [by, bx], order=1, mode='constant',
                                                      cval=0.0)))
    # sky neutralisation: subtract each channel's border median
    border = 10
    for i, c in enumerate(channels):
        rim = np.concatenate([c[:border].ravel(), c[-border:].ravel(),
                              c[:, :border].ravel(), c[:, -border:].ravel()])
        channels[i] = c - np.nanmedian(rim)
    channels = [np.clip(c, 0.0, None) for c in channels]
    if balance:
        channels = [c / max(np.percentile(c, 99.9), 1e-30) for c in channels]

    R = channels[-1]
    B = channels[0]
    if len(channels) == 2:
        G = np.sqrt(R * B)
    else:
        G = np.mean(channels[1:-1], axis=0) if len(channels) > 3 else channels[1]
    hi = np.percentile(R, 99.9)
    rgb = make_lupton_rgb(R, G, B, stretch=max(stretch_scale * hi, 1e-30), Q=Q,
                          minimum=0.0)
    return rgb, base


def centre_kernel(psf):
    """Return `psf` shifted so its own centre sits on the array centre, renormalised.

    REQUIRED before convolving one band with another's kernel, because these kernels are NOT
    centred (measured over the HST stamp PSFs: median offset 0.2-0.5 px by band, the worst
    ~1.7 px). Convolving with an off-centre kernel TRANSLATES the image by that offset, and a
    kernel swap then prints a bright-plus-hole dipole at the core (the fictitious "2.2 px
    astrometric offset" of detect_arcs.py was exactly this). Anything that convolves these
    kernels should centre them first or know why it does not. (Same algorithm as
    detect_arcs.centre_kernel; carried here so the GUIs do not depend on the detector, which
    is a placeholder in the JWST repo.)
    """
    from photutils.centroids import centroid_com
    from scipy.ndimage import shift as ndshift
    k = np.nan_to_num(np.asarray(psf, dtype=float))
    ny, nx = k.shape
    py, px = np.unravel_index(np.argmax(k), k.shape)
    h = 3
    y0, x0 = max(py - h, 0), max(px - h, 0)
    box = np.clip(k[y0:py + h + 1, x0:px + h + 1], 0, None)
    if box.sum() <= 0:
        return k / k.sum()
    bx, by = centroid_com(box)
    dy, dx = (y0 + by) - (ny - 1) / 2.0, (x0 + bx) - (nx - 1) / 2.0
    if np.hypot(dy, dx) > 0.02 and np.isfinite(dy) and np.isfinite(dx):
        k = ndshift(k, (-dy, -dx), order=3, mode='constant', cval=0.0)
    total = k.sum()
    return k / total if total else k


def colour_index_map(bands, deflector_radius=0.3, snr_floor=3.0, smooth_px=1.0):
    """Log colour-ratio map of the BLUEST over the REDDEST band, relative to the deflector.

    Returns (map, extent_base_band_dict): log10[(B/R) / (B/R)_deflector], so the deflector's
    own colour reads 0 and a pixel reads blue (negative) / red (positive) against it -- the
    "two images of one source share a colour" test made visual. (B/R)_deflector is the
    median ratio within `deflector_radius` arcsec of the stamp centre. Pixels with red-band
    S/N < `snr_floor` are NaN (a ratio of noise is meaningless).

    PSF-MATCHED (2026-11-17, the detect_arcs.py rule): each band is convolved with the OTHER
    band's kernel, so both carry PSF_B * PSF_R and the ratio is not a sharp-minus-blurry
    artefact at every steep gradient (an un-matched f450W/f160W map prints a red ring around
    the deflector core). The kernels are centred first (centre_kernel; an off-centre kernel
    translates whatever it convolves) and resampled onto the finest grid when the pixel
    scales differ. The map is lightly smoothed (1 px). Needs >=2 bands with 'psf' arrays
    (returns (None, None) otherwise)."""
    from scipy.ndimage import map_coordinates, gaussian_filter
    from scipy.signal import fftconvolve

    usable = [b for b in bands if b is not None and b.get('psf') is not None]
    if len(usable) < 2:
        return None, None
    usable = sorted(usable, key=lambda b: band_pivot_wavelength(b['filt']))
    blue, red = usable[0], usable[-1]
    base = min(usable, key=lambda b: b['ps'])
    ps_t = base['ps']
    ny, nx = np.asarray(base['data']).shape
    yy, xx = np.mgrid[:ny, :nx]
    sky = WCS(base['hdr']).pixel_to_world(xx, yy)

    def onto(band, key='data'):
        d = np.nan_to_num(np.asarray(band[key], float))
        if band is base:
            return d
        bx, by = WCS(band['hdr']).world_to_pixel(sky)
        return np.nan_to_num(map_coordinates(d, [by, bx], order=1, mode='constant',
                                             cval=0.0))

    def kernel_onto(band):
        k = np.nan_to_num(np.asarray(band['psf'], float))
        zoom = band['ps'] / ps_t            # kernel px per target px
        if abs(zoom - 1.0) > 1e-3:
            h = np.array(k.shape)
            c = (h - 1) / 2.0
            n_new = (np.maximum(h / zoom, 5) // 2 * 2 + 1).astype(int)
            c_new = (n_new - 1) / 2.0
            ky, kx = np.mgrid[:n_new[0], :n_new[1]]
            k = map_coordinates(k, [c[0] + (ky - c_new[0]) * zoom,
                                    c[1] + (kx - c_new[1]) * zoom],
                                order=3, mode='constant', cval=0.0)
        k = centre_kernel(k)
        return k / k.sum()

    kb, kr = kernel_onto(blue), kernel_onto(red)
    B = fftconvolve(onto(blue), kr, mode='same')          # blue carries PSF_R * PSF_B
    R = fftconvolve(onto(red), kb, mode='same')           # red carries PSF_B * PSF_R
    N_R = onto(red, 'noise') if red.get('noise') is not None else None
    if N_R is not None:
        good = (N_R > 0) & np.isfinite(N_R)
        # Fill invalid pixels with the median valid sigma^2 and mask afterwards by the
        # convolved invalid weight: a huge fill value would ring through the FFT and
        # destroy the whole frame (the 1e8-sentinel trap; measured 2026-11-17 on J0008:
        # 1e30 fill gave var_R = -1e11 at the CENTRE).
        fill = np.median(N_R[good] ** 2) if good.any() else 1.0
        var_R = fftconvolve(np.where(good, N_R ** 2, fill), kb ** 2, mode='same')
        bad_w = fftconvolve((~good).astype(float), kb ** 2, mode='same')
        var_R = np.where(bad_w > 0.5, np.nan, var_R)
        N_R = np.where(np.isfinite(var_R) & (var_R > 0), var_R ** 0.5, np.inf)
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = B / np.where(R > 0, R, np.nan)
        snr_R = R / N_R if N_R is not None else np.full_like(R, np.inf)
    cy, cx = (ny - 1) / 2.0, (nx - 1) / 2.0
    rr = np.hypot(yy - cy, xx - cx) * ps_t        # arcsec
    centre_ok = (rr < deflector_radius) & np.isfinite(ratio) & (snr_R > 2)
    ref = np.nanmedian(ratio[centre_ok])
    if not np.isfinite(ref) or ref <= 0:
        return None, None
    m = np.log10(ratio / ref)
    m = np.where((snr_R >= snr_floor) & np.isfinite(m), m, np.nan)
    # smooth on the finite values only (masked pixels must not bleed in)
    filled = np.where(np.isfinite(m), m, 0.0)
    w = np.isfinite(m).astype(float)
    sm = gaussian_filter(filled, smooth_px) / np.clip(gaussian_filter(w, smooth_px), 1e-9, None)
    sm[w < 0.5] = np.nan
    return sm, base


# FALLBACK order for "best band to draw the mask on" when --filt is not given, used when the
# lens has no primary band set (mast_target_names.primary_band: no best_band_visual and no
# measured best_band_arc_snr). The union of both repos' orders; a repo only ever holds its own
# filters, so the cross-ordering never bites. HST: f814W is the highest-S/N SLACS band, the
# stamp --center-band and BELLS's only band; f606W is the gallery primary; f475X leads
# (2026-10-05): only SL2S has it, and there it is the deep GO-17130 band. JWST: reddest first
# (the arcs are red sources; a fixed order alone picked the wrong band for 13 of 20 lenses).
_BAND_PRIORITY = ['f475X', 'f814W', 'f606W', 'f606W_visit2', 'f555W', 'f160W', 'f438W', 'f450W',
                  'f275W', 'f225W', 'f600LP',
                  'F444W', 'F356W', 'F410M', 'F277W', 'F200W', 'F150W', 'F115W', 'F090W',
                  'F070W']


def _band_rank(filt):
    try:
        return _BAND_PRIORITY.index(filt)
    except ValueError:
        return len(_BAND_PRIORITY)


def band_order(filts, lens=None, sample=None):
    """`filts` in preference order: the lens's primary band first (if it has one and it is
    among `filts`), then _BAND_PRIORITY. make_positions and make_arc_masks share it, so the
    GUIs agree on which band is primary."""
    order = sorted(filts, key=lambda f: (_band_rank(f), f))
    if lens is not None:
        primary, _ = mast_target_names.primary_band(lens, sample)
        if primary in order:
            order.remove(primary)
            order.insert(0, primary)
    return order


def pick_display_filt(filts, requested, lens=None, sample=None):
    """Choose which band to draw on. `requested` (--filt) wins if present among `filts`;
    else the lens's PRIMARY band (mast_target_names.primary_band: the by-eye
    `best_band_visual`, else choose_primary_band.py's measured `best_band_arc_snr`) when it is
    among `filts`; else the highest-priority available band (band_order). Returns None if
    the requested band is absent. Shared with make_positions and make_arc_masks so all three
    GUIs pick the same primary band per lens.
    """
    if requested is not None:
        return requested if requested in filts else None
    return band_order(filts, lens, sample)[0]


def reproject_mask_bool(src_bool, src_wcs, dst_wcs, dst_shape):
    """Rigorously reproject a boolean mask from the source stamp's WCS grid onto a target
    stamp's WCS grid (nearest-neighbour), returning a bool array of shape `dst_shape`.

    Kept over autolens's `autogalaxy.gui.display_util.mask_2d_regridded_from`, which regrids
    between two centred uniform grids by pixel scale alone: this one goes through each band's
    own sci-header WCS, so it also absorbs any band-to-band offset of the stamp centres.

    Each target pixel's sky position (via `dst_wcs`) is mapped back to a source pixel (via
    `src_wcs`) and the nearest source mask value taken -- so the mask lands on the SAME sky
    region regardless of pixel scale/size. For the 0.05" optical bands this is near-identity
    (they share a grid to <1px); for f160W (0.06", 200px) it is the only correct transfer (a
    common sky point is ~20px away by array index). Nearest-neighbour (order 0) is right for a
    keep/exclude boolean -- no interpolation across the True/False boundary. Target pixels that
    fall outside the source footprint (only sub-pixel stamp-edge slivers, since every band cuts
    the same 12" angular size) default to False = kept.
    """
    from scipy.ndimage import map_coordinates
    ny, nx = dst_shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    sky = dst_wcs.pixel_to_world(xx.ravel(), yy.ravel())      # 0-based (x=col, y=row)
    sx, sy = src_wcs.world_to_pixel(sky)                      # float source (col, row)
    vals = map_coordinates(np.asarray(src_bool, dtype=float),
                           np.vstack([sy, sx]),               # map_coordinates wants [row, col]
                           order=0, mode='constant', cval=0.0)
    return vals.reshape(ny, nx) > 0.5


def find_proposal_mask(filt_dirs, display_filt, display_hdr, drizzle_pass, propose_from,
                       force=False, lens=None, sample=None):
    """Find a mask already drawn for this lens to PROPOSE as the starting point for the band
    about to be drawn, reprojected onto that band's grid. Returns (source_filt, bool array)
    or (None, None).

    This is the reviewed alternative to `--broadcast`: instead of writing another band's mask
    out blind, it is shown for approval and can be edited or rejected (see draw_mask_gui and
    process_lens_mask). Source preference, `propose_from='auto'`:
      1. the draw band's OWN existing mask, when there is one -- i.e. a `--force` redraw
         starts from the mask being replaced, so refining it does not mean redrawing it;
      2. for a split-visit band (`f450W_visit2`), the PRIMARY visit's mask in the SAME
         filter (`f450W`): same field, same PSF, same depth -- only the grid moved, and the
         WCS reprojection below lands it exactly. Beats the cross-filter primary band
         (2026-11-25);
      3. otherwise the lens's primary band, then the highest-priority OTHER band that has a
         mask (band_order).
    An explicit `propose_from='<band>'` restricts the search to that band; 'none' is handled
    by the caller (no search at all).

    The transfer is the same rigorous per-pixel WCS reprojection `--broadcast` uses -- near
    identity across the bands that share a plate scale, a real regrid onto a different one.
    """
    if propose_from in (None, 'none'):
        return None, None
    if propose_from == 'auto':
        others = [f for f in band_order(filt_dirs, lens, sample) if f != display_filt]
        base_filt = display_filt.split('_visit')[0]
        if base_filt != display_filt and base_filt in others:
            others.remove(base_filt)
            others.insert(0, base_filt)
        candidates = ([display_filt] if force else []) + others
    else:
        candidates = [propose_from]

    dst_wcs = WCS(display_hdr).celestial
    dst_shape = (display_hdr['NAXIS2'], display_hdr['NAXIS1'])
    for cand in candidates:
        stamp_dir = filt_dirs.get(cand)
        if stamp_dir is not None:
            prefix = find_prefix(stamp_dir, drizzle_pass)
            if prefix is None:
                continue
            mask_path = os.path.join(stamp_dir, f'{prefix}_mask.fits')
            if not os.path.exists(mask_path):
                continue
            src_bool = np.asarray(fits.getdata(mask_path), dtype=bool)
            if cand == display_filt:
                return cand, src_bool                 # same grid -- no reprojection needed
            # The mask FITS is a bare array (aplt.fits_array), so its band's sci header
            # carries the WCS -- exactly as the --broadcast path resolves the target grid.
            src_hdr = fits.getheader(stamp_paths.find_stamp(stamp_dir, prefix, 'sci'))
            return cand, reproject_mask_bool(src_bool, WCS(src_hdr).celestial,
                                             dst_wcs, dst_shape)
    return None, None


def confirm_proposal(lens, filt, source_filt, n_prop, n_add, n_erase, n_final, n_drawn):
    """Ask what to do with a reviewed proposal, AFTER the GUI has been closed -- so the
    decision is made having actually seen the mask over this band's own image.

    Three outcomes, and an empty answer means the common one (apply):
      apply  -- write proposal + white additions - black erasures (the reviewed mask);
      drawn  -- reject the proposal outright and keep only what was painted white;
      skip   -- write NOTHING, leaving the band pending for a later run, which is also
                what a non-interactive stdin or Ctrl-C gets (never silently write a mask
                nobody approved).
    """
    print(f"\n  {lens} {filt}: proposal from {source_filt} = {n_prop}px; "
          f"you added {n_add}px, erased {n_erase}px")
    print(f"    [a] apply proposal + your edits  -> {n_final}px excluded   (default)")
    print(f"    [d] drawn only, reject the proposal -> {n_drawn}px excluded")
    print(f"    [s] skip this band, write nothing")
    try:
        ans = input("  choose [A]/d/s: ").strip().lower()
    except (EOFError, KeyboardInterrupt):
        print()
        return 'skip'
    if ans in ('', 'a', 'apply', 'y'):
        return 'apply'
    if ans in ('d', 'drawn', 'r', 'reject'):
        return 'drawn'
    return 'skip'


def load_colour_bands(filt_dirs, drizzle_pass='auto'):
    """Band dicts ({filt, data, hdr, ps, noise, psf}, raw sci flux) for lupton_rgb /
    colour_index_map, for every band of a lens with a stamp for the pass. Shared by
    make_masks.py / make_arc_masks.py (the Scribbler's colour panel) and make_positions.py
    (the colour column). psf/noise are loaded when present (colour_index_map needs them for
    the PSF-matching and the S/N floor; lupton_rgb ignores them)."""
    bands = []
    for filt in sorted(filt_dirs):
        prefix = find_prefix(filt_dirs[filt], drizzle_pass)
        if prefix is None:
            continue
        sci_path = stamp_paths.find_stamp(filt_dirs[filt], prefix, 'sci')
        hdr = fits.getheader(sci_path)
        noise_path = stamp_paths.find_stamp(filt_dirs[filt], prefix, 'noise')
        psf_path = os.path.join(filt_dirs[filt], 'psf.fits')
        bands.append({'filt': filt, 'data': fits.getdata(sci_path), 'hdr': hdr,
                      'ps': pixel_scale_from_header(hdr),
                      'noise': (np.asarray(fits.getdata(noise_path), float)
                                if noise_path else None),
                      'psf': (np.nan_to_num(np.asarray(fits.getdata(psf_path), float))
                              if os.path.exists(psf_path) else None)})
    return bands


def draw_mask_gui(display_dir, display_prefix, lens, filt, brush_radius, brush_width,
                  display, subtract_radial=False, side_by_side=True, prompt=None,
                  proposal=None, proposal_label='proposed mask', positions=None,
                  position_marker_size=None, mask_overlay=None, title=None,
                  colour_bands=None):
    """Run the Scribbler GUI once on the chosen band, blocking until it is closed. Returns
    (scribbler, pixel_scales, sci_header, brush_width, start_radius, source_label); read the
    closed `scribbler` back with `get_scribble_masks()` (each brush on this band's grid) and
    `mask_from()` (the reviewed mask, see below).

    The window is al.Scribbler's own (autolens >= 2026.10.7) -- this function only loads the
    base array, derives the brush size and passes everything through. What it shows:

    TWO BRUSHES. '1' WHITE (`al.Scribbler.ADD_SEGMENT`) adds to the mask, '2' BLACK
    (`ERASE_SEGMENT`) erases from it; press the key to switch, the brush ring follows the
    cursor in the active colour. `get_scribble_masks()[ADD_SEGMENT]` is the white scribble
    and `[ERASE_SEGMENT]` the black one, each already on the band grid (a stroke on either
    side-by-side panel is clipped to the panel it was drawn in and folded onto the one
    image). `mask_from()` is `(proposal | painted) & ~erased`, with `proposal` the one the
    window was opened with. Erasing is what makes an inherited mask editable rather than
    all-or-nothing.

    `proposal` (bool array on this band's grid, optional) is a mask inherited from another
    band (or from the arc detector), passed to `al.Scribbler(proposal=)`: its boundary is
    drawn as a contour (black on a white halo) on every panel -- where that edge falls
    against this band's own structure being the whole question. Only an outline: a filled-in
    "mask applied" panel was tried and dropped, since blanking the interior showed nothing the
    boundary does not already carry while hiding the pixels being judged.

    `subtract_radial` shows the deflector's radial-median-subtracted image (the Scribbler's
    `subtract_radial=True`, `display_util.radial_median_subtract`) -- DISPLAY ONLY, the
    scribble is read back as brush positions so it cannot change what a stroke masks. Off by
    default here (a contaminant mask is drawn against the real sky), on by default in
    make_arc_masks.py, where the arcs are the whole subject.

    `side_by_side` (default True, only meaningful with subtract_radial) shows BOTH views at
    once -- subtracted LEFT, as-observed RIGHT, separated by a blank gutter -- because each
    answers a different question: the subtracted panel is where the arcs are visible at all,
    the as-observed panel is where a contaminant's real extent and the galaxy envelope are.
    Scribble on EITHER; strokes land on the one mask at the same sky position.

    STRETCH. The raw base array (S/N or signal, --display) is handed over untouched and the
    Scribbler maps each panel through display_look.SCRIBBLER -- the dataset-subplot look:
    inferno with the panel's own PercentileInterval(99) + AsinhStretch(0.1) (2026-11-28; before
    the merge the HST side passed bare cmap='inferno' with the Scribbler's own median/MAD
    arcsinh). The repo's `stretched_display` and its
    --stretch/--asinh-a/--vmin-percent/--vmax-percent levers no longer apply to this GUI; they
    remain for the Clicker and the QC PNGs.

    `positions` (list of (y, x) arcsec, e.g. make_positions.py's) are marked with a small
    '+' cross on every panel (`position_marker_size` px arms; None = the Scribbler's own,
    1.5% of the stamp); `mask_overlay` (an al.Mask2D) has its edge scattered as black 'x's.
    Both display only. `title` is the first line of the window title (default names the
    lens, band and polarity); the Scribbler appends the panel layout, a legend for every
    marker it is showing, and the keys. `prompt` replaces the printed instruction lines (so
    the arc tool can state its own, opposite, polarity).

    POLARITY. Here painted = REMOVE from the fit (painted = excluded, unpainted = kept --
    al.Mask2D's own convention, so the scribble IS the mask, no inversion; opposite to
    autolens_workspace's mask.py, which scribbles the region to keep).
    """
    sci_path = stamp_paths.find_stamp(display_dir, display_prefix, 'sci')
    with fits.open(sci_path) as hdul:
        sci_hdr = hdul[0].header
    pixel_scales = pixel_scale_from_header(sci_hdr)

    print(f"\n{lens}  drawing on {filt}  [{display_prefix} pass]  "
          f"pixel_scale={pixel_scales:.4f}\"/pix")
    print(f"  {sci_path}")
    for line in (prompt or [
            "  Scribbler GUI: scribble the areas to REMOVE from the fit (contaminants, field",
            "    sources, bad regions). Everything you do NOT paint is kept in the fit."]):
        print(line)

    base, source_label = load_display_base(display_dir, display_prefix, display, pixel_scales)
    # The panels the Scribbler will build (its own names), for the provenance label.
    if subtract_radial:
        panel_names = ['radial-subtracted'] + (['as-observed'] if side_by_side else [])
    else:
        panel_names = ['as-observed']
    source_label += ', ' + ' + '.join(panel_names) + ' (Scribbler arcsinh stretch)'
    if proposal is not None:
        source_label += f' ({proposal_label} outlined)'
    # al.Scribbler sizes the brush as max(int(image_height * brush_width), 1), a fraction. To
    # get a stamp-independent default (--brush-radius px) we derive the fraction from this
    # stamp's own height; an explicit --brush-width fraction (if given) overrides it.
    n_pix = base.shape_native[0]
    if brush_width is None:
        brush_width = brush_radius / n_pix
    start_radius = max(int(n_pix * brush_width), 1)
    print(f"  display: {source_label}")
    if len(panel_names) > 1:
        print(f"  2 PANELS side by side -- LEFT {panel_names[0]}, RIGHT {panel_names[1]}. "
              f"Scribble on either; both land on the one mask (same sky position).")
    if proposal is not None:
        print(f"  the proposed mask ({proposal_label}) is OUTLINED on every panel.")
    print("    keys '1' = WHITE brush, ADD to the mask   |   '2' = BLACK brush, ERASE from it")
    print(f"  brush start radius = {start_radius} px")

    if title is None:
        title = f"{lens} {filt}: scribble what to REMOVE from the fit"
        if proposal is not None:
            title += f"  --  reviewing the {proposal_label}"

    # Colour panel (2026-11-17, --colour default on; 2026-11-25 the COLOUR-INDEX map, not the
    # raw-flux RGB): log10(bluest/reddest vs the deflector's own ratio, PSF-matched), on RdBu
    # -- blue = bluer than the deflector, red = redder, deflector white. The raw Lupton RGB
    # was uninformative on a disk-dominated lens (the red deflector owns every pixel: J0841,
    # J0822) and, being RGB, can take no colourbar; the index map is scalar, so it has one.
    # No ticks or axis labels on the panel (Sam, 2026-11-25).
    # Shown IN the Scribbler's window as a left-hand panel -- one window total (Sam's call
    # 2026-11-24; the separate companion window was refused). The Scribbler only paints on
    # self.ax (its mouse callbacks guard on event.inaxes), so the added axes is display-only
    # and safe. Its own rgb_image= hook is NOT usable: it draws the RGB on subplot(121) and
    # the image on subplot(111) OVER it, hiding the colour panel entirely (verified
    # 2026-11-24). So: build with block=False, shrink self.ax to the right, add the colour
    # axes on the left, then start().
    cmap_data = None
    if colour_bands:
        cmap_data, cbase = colour_index_map(colour_bands)

    # Blocks at start() until Esc / q closes the window.
    scribbler = al.Scribbler(
        image=base, brush_width=brush_width, proposal=proposal,
        subtract_radial=subtract_radial, side_by_side=side_by_side, panel_gap=_PANEL_GAP,
        positions=positions or None, position_marker_size=position_marker_size,
        mask_overlay=mask_overlay, title=title, block=False,
        figsize=((20, 10) if cmap_data is not None else (15, 15)),
        cmap=display_look.SCRIBBLER)   # the dataset-subplot look (inferno, per-panel 99% asinh)
    # The Scribbler scatters the mask-overlay edge in black, invisible on inferno's dark sky --
    # recolour it cyan.
    from matplotlib.collections import PathCollection
    for coll in scribbler.ax.collections:
        if type(coll) is PathCollection:
            coll.set_color('cyan')
    scribbler.ax.set_title(scribbler.title_text + '\nNorth down, East left', fontsize=10)
    if cmap_data is not None:
        names = sorted((b['filt'] for b in colour_bands), key=band_pivot_wavelength)
        n_y, n_x = cmap_data.shape
        hw = int(n_x / 2) * cbase['ps']
        scribbler.ax.set_position([0.335, 0.03, 0.65, 0.90])
        axc = scribbler.figure.add_axes([0.03, 0.10, 0.27, 0.72])
        imc = axc.imshow(cmap_data, extent=[-hw, hw, -hw, hw], cmap='RdBu',
                         vmin=-0.3, vmax=0.3)
        axc.set_xticks([])
        axc.set_yticks([])
        axc.set_title(f'colour index ({names[0]}/{names[-1]})\nblue = bluer than the '
                      f'deflector', fontsize=10)
        cb = scribbler.figure.colorbar(imc, ax=axc, fraction=0.046)
        cb.ax.tick_params(labelsize=8)
        cb.set_label(f'log10({names[0]}/{names[-1]} vs deflector)', fontsize=8)
        print(f"  colour index is the LEFT panel of the mask window ({names[0]}/{names[-1]} "
              f"vs the deflector; blue = bluer) -- reference only, paint on the data panels")
    scribbler.start()
    _drain_tk()
    return scribbler, pixel_scales, sci_hdr, brush_width, start_radius, source_label


def _drain_tk():
    """Tidy up after the Scribbler window closes: switch pyplot's interactive mode back OFF, and
    flush Tk's pending events now with Tcl background errors silenced.

    Interactive mode: the Scribbler's `start()` calls `plt.ion()` and (in autogalaxy
    <= 2026.10.12.1) never switches it off. Every figure drawn later in this process -- the
    dataset QC subplot rebuilt right after "a" -- then pops up as a window for a moment before
    it is saved and closed (2026-10-14: one FigureManagerTk.show() per rebuild with ion on,
    none with it off).

    Pending events: see below.

    On Esc the Scribbler calls `canvas.draw_idle()` and then `plt.close()`, which leaves an
    `idle_draw` callback queued for a canvas that no longer exists. Python's Tk input hook
    services Tk events while `input()` waits, so that stale callback fired at the prompt and
    printed `invalid command name "...idle_draw"` over it (J0217-0513 f160W, 2026-10-14). It is
    harmless (the canvas is gone), so it is run here, quietly, while nothing is waiting on stdin.
    The real fix is upstream: Scribbler.on_keypress should not `draw_idle()` right before
    `plt.close()`."""
    import matplotlib.pyplot as plt
    plt.ioff()
    try:
        import tkinter
        root = tkinter._default_root
        if root is None:
            return
        root.tk.eval('proc bgerror {args} {}')
        root.update()
    except Exception:
        pass


#: (lens, filt, size) of every band whose mask this run wrote: their dataset subplots are drawn
#: by rerender_overviews' background process at the end of the run.
PENDING_SUBPLOTS = []


def rerender_overviews(sample, lenses, subplots=()):
    """Re-render, at the END of the run and in ONE DETACHED background process, the dataset
    subplot of every band in `subplots` ((lens, filt, size); make_dataset_subplots.py --force)
    and then the band overview (make_band_overview.py) of every lens in `lenses`.

    Neither may run in the GUI process. Done after each lens, the overview put an isophote fit
    (~10 s) between one GUI and the next; done in this process at the end, it left the GUI app
    alive but busy and not answering macOS, which showed as a blank window and a spinning
    cursor, like a hang. The dataset subplot drawn inline also flashed up as a blank window:
    al.Scribbler switches pyplot to interactive mode (plt.ion), so any later pyplot figure,
    autolens's subplot included, is shown as it is made. Both are best-effort QC renders (the
    mask is the product); the log says what failed."""
    lenses = list(lenses)
    subplots = list(subplots)
    if not lenses and not subplots:
        return
    import shlex
    import subprocess
    import time
    try:
        import matplotlib.pyplot as plt
        plt.close('all')
    except Exception:
        pass
    here = os.path.dirname(os.path.abspath(__file__))
    py = shlex.quote(sys.executable)
    steps = [f'{py} {shlex.quote(os.path.join(here, "make_dataset_subplots.py"))} --sample {sample} '
             f'--lens {shlex.quote(l)} --filt {f} --force' + (f' --size {s:g}' if s is not None else '')
             for l, f, s in subplots]
    if lenses:
        steps.append(f'{py} {shlex.quote(os.path.join(here, "make_band_overview.py"))} --sample {sample} --lens '
                     + ' '.join(shlex.quote(l) for l in lenses))
    log = os.path.join(ws_path, 'data', 'run_logs', f'band_overview_{sample}_{time.strftime("%Y%m%d_%H%M%S")}.log')
    os.makedirs(os.path.dirname(log), exist_ok=True)
    with open(log, 'w') as fh:
        subprocess.Popen(['/bin/sh', '-c', '; '.join(steps) + '; echo "=== renders done"'],
                         stdout=fh, stderr=subprocess.STDOUT, start_new_session=True,
                         env=dict(os.environ, MPLBACKEND='Agg'))
    what = ' + '.join(x for x in (f'{len(subplots)} dataset subplot(s)' if subplots else '',
                                  f'{len(lenses)} band overview(s)' if lenses else '') if x)
    print(f"\n{what} rendering IN THE BACKGROUND; log: {os.path.relpath(log, ws_path)}")


def scribbles_from(scribbler, proposal=None):
    """Read a closed Scribbler back as (painted, erased, drawn_only, reviewed), all bool
    arrays on the band grid:

      painted / erased  -- the WHITE and BLACK brushes (`get_scribble_masks()`), for the
                           provenance counts;
      drawn_only        -- `painted & ~erased`, what you drew with the proposal ignored;
      reviewed          -- `scribbler.mask_from()` = `(proposal | painted) & ~erased`, the
                           proposal as edited (equal to drawn_only when there is none).
    """
    segments = scribbler.get_scribble_masks()
    painted = np.asarray(segments[al.Scribbler.ADD_SEGMENT], dtype=bool)
    erased = np.asarray(segments[al.Scribbler.ERASE_SEGMENT], dtype=bool)
    reviewed = np.asarray(scribbler.mask_from(proposal), dtype=bool)
    return painted, erased, painted & ~erased, reviewed


def process_lens_mask(lens, filt_dirs, sample, requested_filt, drizzle_pass, force,
                      broadcast=False, brush_radius=6, brush_width=None, display='snr',
                      subtract_radial=None, side_by_side=True,
                      size=None, dataset_subplot=True,
                      propose_from='auto', colour=True):
    """Draw a mask once for one lens, on its best/forced band. By default the mask is written
    for THAT BAND ONLY; `broadcast=True` also writes it to the lens's other bands by WCS
    reprojection. `filt_dirs` maps filt -> stamp_dir from discover_targets; the mask FITS
    goes into that dir, beside the sci it describes. Returns True if a mask was written.

    `propose_from` ('auto' by default, 'none' to disable, or an explicit band) turns this into
    a REVIEW of a mask already drawn for another band of this lens instead of a blank canvas:
    the inherited mask is reprojected onto this band and outlined on every panel, and you then
    add (white) / erase (black) parts of it before choosing to apply, reject, or skip. It is
    the deliberate alternative to `--broadcast`, which writes the
    same mask to every band unseen -- what a mask should exclude is not band-independent in
    practice, so the transfer is worth looking at in the band it is landing on.
    """
    filts = sorted(filt_dirs)
    display_filt = pick_display_filt(filts, requested_filt, lens=lens, sample=sample)
    if display_filt is None:
        print(f"{lens}: requested --filt {requested_filt} not present "
              f"(have {', '.join(filts)}), skipping")
        return False

    display_dir = filt_dirs[display_filt]                # stamp dir of the draw band
    display_prefix = find_prefix(display_dir, drizzle_pass)
    if display_prefix is None:
        print(f"{lens} {display_filt}: no stamp sci for --pass {drizzle_pass}, skipping")
        return False

    # Skip if the draw band already has a mask, unless --force.
    if (os.path.exists(os.path.join(display_dir, f'{display_prefix}_mask.fits'))
            and not force):
        print(f"{lens}: mask already exists ({display_filt}), skipping (--force to redraw)")
        return False

    display_hdr = fits.getheader(stamp_paths.find_stamp(display_dir, display_prefix, 'sci'))
    proposal_from, proposal = find_proposal_mask(
        filt_dirs, display_filt, display_hdr, drizzle_pass, propose_from, force=force,
        lens=lens, sample=sample)
    if proposal is not None:
        print(f"{lens} {display_filt}: reviewing the mask from {proposal_from} "
              f"({int(proposal.sum())}px) -- edit it in the GUI, then apply/reject")
    # 'auto' radial subtraction: on when there is a proposal to review (you need to see the
    # arcs the inherited mask may be clipping), off otherwise -- a contaminant mask drawn from
    # scratch is normally judged against the real sky. An explicit --(no-)subtract-radial wins.
    if subtract_radial is None:
        subtract_radial = proposal is not None

    (scribbler, src_ps, src_hdr, brush_width, start_radius,
     source_label) = draw_mask_gui(
        display_dir, display_prefix, lens, display_filt, brush_radius, brush_width, display,
        subtract_radial=subtract_radial, side_by_side=side_by_side,
        proposal=proposal, proposal_label=f'{proposal_from} mask' if proposal_from else None,
        colour_bands=(load_colour_bands(filt_dirs, drizzle_pass)
                      if colour and len(filt_dirs) > 1 else None))

    painted, erased, drawn_only, reviewed = scribbles_from(scribbler)
    n_add, n_erase = int(painted.sum()), int(erased.sum())

    # Nothing painted = you skipped this lens: write NO file. An all-False mask would be a
    # legitimate "exclude nothing" mask, which is exactly the ambiguity to avoid -- an absent
    # file means "not done yet" and the next run re-offers the lens, where an empty one would
    # silently mark it finished. The same rule holds after a review: an approved mask that
    # ends up empty (everything erased) is not written either.
    if proposal is None:
        scribbled, source_tag = drawn_only, 'drawn'
    else:
        # reviewed = scribbler.mask_from() = (proposal | painted) & ~erased
        choice = confirm_proposal(lens, display_filt, proposal_from, int(proposal.sum()),
                                  n_add, n_erase, int(reviewed.sum()), int(drawn_only.sum()))
        if choice == 'skip':
            print(f"  skipped -- no mask written for {lens} {display_filt} (still pending)")
            return False
        if choice == 'drawn':
            scribbled, source_tag = drawn_only, 'drawn'
            print(f"  proposal from {proposal_from} REJECTED -- keeping only what you painted")
        else:
            scribbled = reviewed
            source_tag = (f'edited_from_{proposal_from}' if (n_add or n_erase)
                          else f'accepted_from_{proposal_from}')
    if not np.asarray(scribbled).any():
        print(f"  nothing drawn -- no mask written for {lens} (still pending)")
        return False
    src_wcs = WCS(src_hdr).celestial

    # Default: the draw band only, using the scribbled array verbatim. --broadcast additionally
    # writes every other band, WCS-reprojected onto its own grid (near-identity for the 0.05"
    # optical bands, a real ~20px regrid for f160W).
    targets = (filts if broadcast else [display_filt])
    for filt in targets:
        stamp_dir = filt_dirs[filt]
        prefix = find_prefix(stamp_dir, drizzle_pass)
        if prefix is None:
            continue
        band_hdr = fits.getheader(stamp_paths.find_stamp(stamp_dir, prefix, 'sci'))
        band_ps = pixel_scale_from_header(band_hdr)
        is_draw = (stamp_dir == display_dir)
        if is_draw:
            out_bool = scribbled
        else:
            out_bool = reproject_mask_bool(scribbled, src_wcs, WCS(band_hdr).celestial,
                                           (band_hdr['NAXIS2'], band_hdr['NAXIS1']))
        mask_path = os.path.join(stamp_dir, f'{prefix}_mask.fits')
        aplt.fits_array(array=al.Mask2D(mask=out_bool, pixel_scales=band_ps),
                        file_path=mask_path, overwrite=True)
        n_excl = int(np.asarray(out_bool).sum())
        print(f"  wrote {mask_path}  ({source_tag if is_draw else f'reprojected<-{display_filt}'}"
              f", {n_excl}px excluded)")
        entry = {
            'prefix': prefix,
            'drizzle_pass': stamp_pass_of(stamp_dir, prefix) or 'single',
            'pixel_scale_arcsec': round(band_ps, 6),
            'display': source_label,
            # Which reduction the stamp was drawn on, from its own header (stamp_paths.py).
            'bcfill': bool(band_hdr.get('BCFILL', False)),
            'source': source_tag if is_draw else f'reprojected_from_{display_filt}',
            'n_excluded_px': n_excl,
        }
        if band_hdr.get('REDUCVAR'):            # the JWST repo's variant card
            entry['variant'] = str(band_hdr['REDUCVAR'])
        if is_draw:
            entry['brush_width'] = round(brush_width, 6)
            entry['brush_radius_px'] = start_radius
            if proposal is not None:
                # What the reviewed mask inherited vs what this session changed by hand --
                # 'accepted_from_x' with 0 edits is a real, deliberate outcome and needs to
                # be distinguishable from a hand-drawn mask that merely resembles one.
                entry['proposal_from'] = proposal_from
                entry['proposal_px'] = int(proposal.sum())
                entry['n_added_px'] = n_add
                entry['n_erased_px'] = n_erase
        info_json.update(MASKS_JSON, sample, lens, filt, entry)
        if dataset_subplot:
            PENDING_SUBPLOTS.append((lens, filt, size))   # drawn in the background at the end
    return True


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--sample', type=mast_target_names.sample_arg, default=None,
                   help='sample subdirectory of data/stamps/. Default: with --lens, the '
                        'lens\'s own sample (each lens is in exactly one); without --lens, '
                        f'{mast_target_names.DEFAULT_SAMPLE}')
    p.add_argument('--lens', default=None,
                   help='restrict to one lens (no --sample needed -- it is resolved from '
                        'the roster); default every lens with stamps in --sample')
    p.add_argument('--filt', default=None,
                   help='force which band you DRAW on (e.g. f606W); default the best available '
                        'band per lens (f475X>f814W>f606W>f555W>f160W>...). Only that band gets a '
                        'mask unless --broadcast is passed')
    p.add_argument('--broadcast', action=argparse.BooleanOptionalAction, default=False,
                   help='ALSO write the drawn mask to the lens\'s other bands, WCS-reprojected '
                        'onto each grid. Default OFF: one draw covers one band, and you run '
                        'again with --filt for the next, because what a mask should exclude is '
                        'not in practice band-independent (a contaminant can be bright in one '
                        'filter and absent in another, and each band\'s depth and PSF move the '
                        'sensible boundary)')
    p.add_argument('--propose-from', default='auto',
                   help="review a mask already drawn for this lens instead of starting from a "
                        "blank canvas: 'auto' (default) proposes the band's own existing mask "
                        'when --force-redrawing it, else the highest-priority other band that '
                        'has one; a band name (e.g. f814W) forces the source; \'none\' disables '
                        'it. The proposal is reprojected onto the drawn band and OUTLINED on '
                        'every panel, and you add (white brush, key 1) / erase (black brush, '
                        'key 2) parts of it before choosing to apply, reject, or skip -- the '
                        'reviewed alternative to --broadcast')
    p.add_argument('--colour', action=argparse.BooleanOptionalAction, default=True,
                   help='show a COLOUR-INDEX panel in the mask window (default on, '
                        '2026-11-25; inactive for single-band lenses): log10(bluest/reddest '
                        'band vs the deflector\'s own ratio, PSF-matched) -- the discriminant '
                        'between blue lensed images and red contaminants. Reference only -- '
                        'paint on the data panels. --no-colour disables')
    p.add_argument('--pass', dest='drizzle_pass', choices=['auto', 'cr', 'nocrrej'],
                   default='auto',
                   help="which stamp pass to mask, matching make_stamps.py's --pass: "
                        "'auto' (default) uses the band's science stamp (stamp<S>_*: the CR pass, "
                        "or the no-CR pass where none exists, F160W); 'cr'/'nocrrej' ask for one pass")
    p.add_argument('--size', type=float, default=None,
                   help='stamp size in arcsec (default: each lens\'s modelling size, its stamp_size in '
                        'info/samples/lens_samples.json, else 12"). Only modelling-size stamps are '
                        'tracked, so another size is skipped unless --include-ignored')
    p.add_argument('--force', action='store_true', default=False,
                   help='redraw a mask that already exists (default: skip it)')
    p.add_argument('--brush-radius', type=int, default=6,
                   help='starting brush radius in PIXELS (default 6). Stamp-independent: the '
                        'fraction al.Scribbler wants is derived from each stamp height. '
                        "Resize live in the GUI with '='/'-'. Overridden by --brush-width")
    p.add_argument('--brush-width', type=float, default=None,
                   help='explicit brush width as a FRACTION of image height (brush_radius = '
                        'int(N_pix * brush_width)); overrides --brush-radius when given. '
                        'Default: unset, so --brush-radius (px) sets the size')
    p.add_argument('--display', choices=['snr', 'sci'], default='snr',
                   help="image to scribble on: 'snr' (default) = signal/noise, which tracks "
                        "statistical significance -- the right criterion for what to keep, "
                        "given this pipeline's non-uniform noise map; 'sci' = raw signal, a "
                        "morphology cross-check. Choice does not affect the saved mask. The "
                        "display stretch is al.Scribbler's own arcsinh about the sky level "
                        "(the former --stretch/--asinh-a/--vmin-percent/--vmax-percent "
                        "levers are gone from this tool)")
    p.add_argument('--side-by-side', action=argparse.BooleanOptionalAction, default=True,
                   help='with --subtract-radial, show BOTH views at once -- subtracted left, '
                        'as-observed right -- and accept scribbles on either panel, combining '
                        'them onto the one mask (default on; ignored without --subtract-radial, '
                        'where there is nothing to compare)')
    p.add_argument('--subtract-radial', action=argparse.BooleanOptionalAction, default=None,
                   help="subtract the deflector's azimuthally-averaged radial profile from "
                        'the DISPLAYED image, so the arcs stand clear of the galaxy envelope '
                        '(same lever as make_positions.py, where it is ON by default). Default '
                        'is AUTO: on when a mask is being reviewed (--propose-from), where you '
                        'need to see the arcs the inherited mask may be clipping, and off when '
                        'drawing from scratch, since a contaminant mask is normally judged '
                        'against the real sky. Display only -- the mask is built from brush '
                        'positions, so it cannot change what a stroke masks')
    p.add_argument('--dataset-subplot', action=argparse.BooleanOptionalAction, default=True,
                   help='after writing each mask, regenerate that band\'s '
                        '{prefix}_dataset.png QC subplot (scripts/make_dataset_subplots.py) '
                        'so it shows the mask just drawn rather than the previous one '
                        '(default on; best-effort -- a failure never costs the mask)')
    p.add_argument('--include-ignored', action='store_true', default=False,
                   help='also offer bands whose stamp directory is GITIGNORED. Off by '
                        'default: such a band is not a science product (the three 420 s '
                        'SLACS SNAP f814W diagnostics, since deleted, are why the rule exists), and hand-drawn '
                        'work written there is non-regenerable and invisible to every clone. '
                        'Use only to work on a diagnostic deliberately.')
    a = p.parse_args()
    mast_target_names.resolve_sample(p, a, default=mast_target_names.DEFAULT_SAMPLE)

    # --broadcast writes one draw to every band unseen; --propose-from reviews another band's
    # mask in the band it is landing on. They are two answers to the same question, and
    # combining them would write a reviewed-for-f555W mask back over f814W's own.
    if a.broadcast and a.propose_from != 'none':
        p.error('--broadcast and --propose-from are alternatives: pass --propose-from none '
                'to broadcast one draw to every band unseen, or drop --broadcast to review '
                'the inherited mask per band')

    global STAMP_SIZE
    STAMP_SIZE = a.size
    root = stamp_paths.stamps_root(ws_path)

    # discover_targets yields per (lens, filt); regroup to filt -> stamp_dir per lens so we
    # know every band of a lens (filt=None: consider all bands; pick_display_filt picks one).
    lens_filts = {}
    for lens, filt, stamp_dir in discover_targets(root, a.sample, a.lens, None,
                                                   include_ignored=a.include_ignored):
        lens_filts.setdefault(lens, {})[filt] = stamp_dir

    if not lens_filts:
        raise SystemExit(f"no stamps found under {root} for sample {a.sample} matching "
                         f"lens={a.lens!r}")

    print(f"{len(lens_filts)} lens(es) to process")
    made = skipped = 0
    changed = []
    for lens in sorted(lens_filts):
        if process_lens_mask(lens, lens_filts[lens], a.sample, a.filt, a.drizzle_pass,
                             a.force, a.broadcast, brush_radius=a.brush_radius,
                             brush_width=a.brush_width, display=a.display,
                             subtract_radial=a.subtract_radial,
                             side_by_side=a.side_by_side, size=a.size,
                             dataset_subplot=a.dataset_subplot,
                             propose_from=a.propose_from, colour=a.colour):
            made += 1
            changed.append(lens)
        else:
            skipped += 1

    print(f"\nDone: {made} lens(es) masked, {skipped} skipped")
    rerender_overviews(a.sample, changed, PENDING_SUBPLOTS)


if __name__ == '__main__':
    main()
