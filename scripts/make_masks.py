#!/usr/bin/env python
"""
Interactive GUI mask-making, ONE draw per band, seeded by a REVIEWED mask from another band.

For each lens in a sample, launches PyAutoLens's `Scribbler` GUI over the lens's
best-available cutout -- the same tool as
autolens_workspace:scripts/imaging/data_preparation/gui/mask.py -- so you can scribble the
areas to REMOVE from the fit (painted = excluded, everything unpainted = kept; this is the
opposite polarity to the autolens_workspace GUI, which scribbles the region to keep), then
writes the result as cutout_[cr_]mask.fits alongside that cutout's sci/noise/psf products.

ONE DRAW PER BAND. Each run draws on ONE band -- the highest-S/N one available (f814W >
f606W > f555W > f160W > ...) unless --filt forces it -- and writes the mask for THAT BAND
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
judged. You then edit the proposal with two brushes --
'1' WHITE adds to the mask, '2' BLACK erases from it, scribbling on whichever panel you like --
and on closing the GUI choose to [a]pply proposal+edits, keep only what you [d]rew (rejecting
the proposal outright), or [s]kip the band. Source: the band's own existing mask when
--force-redrawing it (so refining a mask does not mean redrawing it), else the highest-priority
other band that has one. --propose-from <band> forces the source, --propose-from none draws on
a blank canvas.

--broadcast (mutually exclusive with the above) is the old unreviewed route: write the one
drawn mask to every band of the lens, sight-unseen. Both share the same transfer -- a rigorous
per-pixel WCS reprojection (nearest-neighbour, source sci-header WCS -> target sci-header WCS),
NOT an array-index copy: the 0.05" optical bands (f814W/f606W/f555W) share a grid to <1px so it
is nearly identity, but f160W is a genuinely different grid (0.06"/px, 200px vs 240px for a 12"
stamp -- a common sky point sits ~20px away by index), and only reprojection places the mask
correctly there.

Defaults to the pipeline's standard 12" cutout tree (data/cutouts/, cutout_paths.DEFAULT_SIZE)
-- the tree these hand-drawn masks are meant for and the only size-variant tree tracked in
git (see .gitignore) precisely because it carries them, which no script can regenerate.
Pass --size 20 to mask the (untracked, regenerable) 20" tree instead.

One cutout dir per band, so one mask per band, written beside the sci it was drawn on.
The bad-column-filled re-drizzle (--bcfill) supersedes the standard one in place rather
than starting a parallel tree (cutout_paths.py), so for ACS (f814W/f555W) + WFPC2 (f606W)
you are scribbling on the bcfill image -- the same geometry with its dead-column noise
stripes filled, i.e. a cleaner image to draw on -- and for f160W and the gallery, which
have no bcfill, on the standard one. Either way there is nothing to choose and no second
copy of the mask to keep in step. A (lens, filt) with a mask already on disk is SKIPPED;
--force redraws it.

Every mask written also REGENERATES that band's `{prefix}_dataset.png` QC subplot
(scripts/make_dataset_subplots.py), so the one PNG that shows the mask laid over
data/noise/S-N always describes the mask currently on disk instead of a previous draw. It is
best-effort: a failure (e.g. no PSF kernel for that band yet) is reported and the mask is
kept. --no-dataset-subplot skips it.

This is a manual, one-lens-at-a-time tool (not a batch driver): each lens blocks on its own
Tk window until you press Esc. A lens that already has a mask is skipped so a run resumes
across a sample; --force redraws.

Usage:
    uv run python scripts/make_masks.py --sample slacs_gold                # best band per lens
    uv run python scripts/make_masks.py --lens J0008-0004                  # one lens, best band
    uv run python scripts/make_masks.py --lens J0008-0004 --filt f606W     # force the draw band
    uv run python scripts/make_masks.py --lens J0008-0004 --filt f160W --force  # redo one band
    uv run python scripts/make_masks.py --sample slacs_gold --filt f555W    # review f814W's
    uv run python scripts/make_masks.py --lens J0008-0004 --propose-from none  # blank canvas
    uv run python scripts/make_masks.py --lens J0008-0004 --broadcast --propose-from none
    uv run python scripts/make_masks.py --sample slacs_gold --force        # redraw everything
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
import mast_target_names
import info_json
import cutout_paths

MASKS_JSON = os.path.join(ws_path, 'info', 'lens_masks.json')


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

    NOT used by the Scribbler GUIs any more: al.Scribbler applies its own arcsinh stretch
    (median / MAD-scaled, clipped at the 0.5 / 99.5 percentiles of the stretched values --
    `autogalaxy.gui.display_util.arcsinh_stretch`) to each panel it shows, so draw_mask_gui
    hands it the raw base array. This percentile + astropy AsinhStretch(a) version is kept
    for the Clicker, which imshows its own axes, and for the arc-mask QC PNG, so those
    products render exactly as before.

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

    Cutouts are recentred on the deflector (`make_cutouts.py --center-band`), so the stamp
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
    """arcsec/pixel from a cutout FITS header's WCS -- Array2D.from_fits wants a plain
    float, not the WCS itself."""
    return float(proj_plane_pixel_scales(WCS(hdr).celestial)[0] * 3600.0)


def find_prefix(cutout_dir, drizzle_pass='auto'):
    """Delegates to cutout_paths.find_prefix: 'cutout_cr' or 'cutout' for this band, whichever
    reduction (standard / bcfill / ...) its stamp is from -- the stamp names are variant-
    tagged, so this is the one place that knows how to find them."""
    return cutout_paths.find_prefix(cutout_dir, drizzle_pass)


def gitignored(paths):
    """Delegates to cutout_paths.gitignored -- see there for why it lives in that module."""
    return cutout_paths.gitignored(paths, ws_path)


def discover_targets(root, sample, lens=None, filt=None, include_ignored=False):
    """Yield (lens, filt, cutout_dir) per band with a cutout under `root` -- the dir both
    scribbled on and written to -- sorted for a reproducible run order.

    GITIGNORED CUTOUT DIRS ARE SKIPPED (2026-09-22, and this is the guard the repo was
    missing). A band whose directory is gitignored is, by construction, not a science
    product -- the three 420 s SLACS SNAP f814W diagnostics (J1134+6027, J1403+0006,
    J1538+5817) are the reason the rule exists. Before this, the sweep tools happily offered
    them, and a hand-drawn mask or a set of hand-marked positions would be written into a
    directory git refuses to store: non-regenerable work that no clone would ever see. It
    already happened -- J1538+5817's image positions were MARKED on its SNAP f814W band and
    broadcast from there to f606W across a 688 mas (13.8 px) band misregistration, and both
    had to be deleted on 2026-09-22.

    `include_ignored=True` opts back in, for deliberately working on a diagnostic product.
    The skip is announced, never silent, so a short target list is explicable.
    """
    found = {}
    pattern = os.path.join(root, sample, lens or '*', filt or '*')
    for cutout_dir in sorted(glob.glob(pattern)):
        if find_prefix(cutout_dir, 'auto') is None:
            continue
        this_filt = os.path.basename(cutout_dir)
        this_lens = os.path.basename(os.path.dirname(cutout_dir))
        found[(this_lens, this_filt)] = cutout_dir
    if not include_ignored:
        ignored = gitignored(list(found.values()))
        skipped = sorted(k for k, v in found.items() if os.path.abspath(v) in ignored)
        for key in skipped:
            del found[key]
        if skipped:
            print(f"  skipping {len(skipped)} gitignored band(s), not science products: "
                  + ', '.join(f'{l}/{f}' for l, f in skipped)
                  + "  (--include-ignored to use them anyway)")
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
        file_path=cutout_paths.find_stamp(display_dir, prefix, 'sci'),
        pixel_scales=pixel_scales).native
    if display == 'sci':
        return sci, 'signal'
    noise_path = cutout_paths.find_stamp(display_dir, prefix, 'noise')
    if noise_path is None:
        print(f"  NOTE: no {prefix}[_<variant>]_noise.fits -- displaying signal instead of S/N")
        return sci, 'signal'
    noise = np.asarray(al.Array2D.from_fits(
        file_path=noise_path, pixel_scales=pixel_scales).native, dtype=float)
    snr = np.zeros_like(noise)
    good = np.isfinite(noise) & (noise > 0)
    snr[good] = np.asarray(sci, dtype=float)[good] / noise[good]
    return al.Array2D.no_mask(values=snr, pixel_scales=pixel_scales).native, 'S/N'


# Preference order for "best band to draw the mask on" when --filt is not given. Same list
# used by make_positions.pick_display_filt (which imports this) so the two GUIs agree on
# which band is primary. f814W is the highest-S/N SLACS band and the cutout --center-band;
# f606W is the gallery/BELLS primary; the rest follow by typical arc S/N. f475X leads (Sam,
# 2026-10-05): only SL2S has it, and there it is the deep GO-17130 band (2366 s, UVIS
# sampling, bright blue arcs), so every Dinos lens -- including J0959+0206, which also has
# F814W -- opens on F475X; snapshot-only lenses open on F475X over F600LP. No SLACS, GALLERY
# or BELLS lens has either band, so their primary is unchanged.
_BAND_PRIORITY = ['f475X', 'f814W', 'f606W', 'f606W_v2', 'f555W', 'f160W', 'f438W', 'f275W',
                  'f225W', 'f600LP']


def _band_rank(filt):
    try:
        return _BAND_PRIORITY.index(filt)
    except ValueError:
        return len(_BAND_PRIORITY)


def pick_display_filt(filts, requested):
    """Choose which band to draw on. `requested` (--filt) wins if present among `filts`;
    else the highest-priority available band. Returns None if the requested band is absent.
    Shared with make_positions so both tools pick the same primary band per lens.
    """
    if requested is not None:
        return requested if requested in filts else None
    return sorted(filts, key=lambda f: (_band_rank(f), f))[0]


def reproject_mask_bool(src_bool, src_wcs, dst_wcs, dst_shape):
    """Rigorously reproject a boolean mask from the source cutout's WCS grid onto a target
    cutout's WCS grid (nearest-neighbour), returning a bool array of shape `dst_shape`.

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
                       force=False):
    """Find a mask already drawn for this lens to PROPOSE as the starting point for the band
    about to be drawn, reprojected onto that band's grid. Returns (source_filt, bool array)
    or (None, None).

    This is the reviewed alternative to `--broadcast`: instead of writing another band's mask
    out blind, it is shown for approval and can be edited or rejected (see draw_mask_gui and
    process_lens_mask). Source preference, `propose_from='auto'`:
      1. the draw band's OWN existing mask, when there is one -- i.e. a `--force` redraw
         starts from the mask being replaced, so refining it does not mean redrawing it;
      2. otherwise the highest-priority OTHER band of this lens that has a mask
         (`_BAND_PRIORITY`, so f814W's mask is what usually seeds f555W/f160W).
    An explicit `propose_from='<band>'` restricts the search to that band; 'none' is handled
    by the caller (no search at all).

    The transfer is the same rigorous per-pixel WCS reprojection `--broadcast` uses -- near
    identity across the 0.05" optical bands, a real ~20px regrid onto f160W's 0.06" grid.
    """
    if propose_from in (None, 'none'):
        return None, None
    if propose_from == 'auto':
        others = sorted((f for f in filt_dirs if f != display_filt),
                        key=lambda f: (_band_rank(f), f))
        candidates = ([display_filt] if force else []) + others
    else:
        candidates = [propose_from]

    dst_wcs = WCS(display_hdr).celestial
    dst_shape = (display_hdr['NAXIS2'], display_hdr['NAXIS1'])
    for cand in candidates:
        cutout_dir = filt_dirs.get(cand)
        if cutout_dir is not None:
            prefix = find_prefix(cutout_dir, drizzle_pass)
            if prefix is None:
                continue
            mask_path = os.path.join(cutout_dir, f'{prefix}_mask.fits')
            if not os.path.exists(mask_path):
                continue
            src_bool = np.asarray(fits.getdata(mask_path), dtype=bool)
            if cand == display_filt:
                return cand, src_bool                 # same grid -- no reprojection needed
            # The mask FITS is a bare array (aplt.fits_array), so its band's sci header
            # carries the WCS -- exactly as the --broadcast path resolves the target grid.
            src_hdr = fits.getheader(cutout_paths.find_stamp(cutout_dir, prefix, 'sci'))
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


def draw_mask_gui(display_dir, display_prefix, lens, filt, brush_radius, brush_width,
                  display, subtract_radial=False, side_by_side=True, prompt=None,
                  proposal=None, proposal_label='proposed mask', positions=None,
                  position_marker_size=None, mask_overlay=None, title=None):
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
    Scribbler maps each panel through its own arcsinh stretch about the panel's sky level
    (median / MAD-scaled, clipped at the 0.5 / 99.5 percentiles of the stretched values).
    The repo's `stretched_display` and its --stretch/--asinh-a/--vmin-percent/--vmax-percent
    levers no longer apply to this GUI; they remain for the Clicker and the QC PNGs.

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
    sci_path = cutout_paths.find_stamp(display_dir, display_prefix, 'sci')
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

    # Blocks here until Esc / q closes the window.
    scribbler = al.Scribbler(
        image=base, brush_width=brush_width, proposal=proposal,
        subtract_radial=subtract_radial, side_by_side=side_by_side,
        positions=positions or None, position_marker_size=position_marker_size,
        mask_overlay=mask_overlay, title=title)
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


def regenerate_dataset_subplot(lens, filt, sample, size, drizzle_pass):
    """Rebuild this band's `{prefix}_dataset.png` (scripts/make_dataset_subplots.py) so the
    QC subplot -- the only pipeline PNG that shows the MASK laid over the data/noise/S-N a
    fit actually consumes -- always matches the mask just written, rather than silently
    continuing to describe the previous draw.

    Best-effort and never fatal. The hand-drawn mask is the non-regenerable product here and
    the subplot is not, so any failure (missing PSF kernel for that band, an autolens
    plotting error) is reported with the command to retry and then swallowed -- it must never
    cost a mask that has just been drawn by hand.

    """
    try:
        import make_dataset_subplots as mds
        root = cutout_paths.cutouts_root(ws_path, size)
        for t_lens, t_filt, cutout_dir in mds.discover_targets(root, sample, lens, filt):
            mds.process(t_lens, t_filt, sample, cutout_dir, drizzle_pass,
                        force=True, check_noise_map=False)
    except Exception as exc:                    # deliberately broad -- see docstring
        print(f"  NOTE: dataset subplot not regenerated for {lens} {filt} ({exc!r}). "
              f"The mask is written; rebuild the subplot with\n"
              f"    uv run python scripts/make_dataset_subplots.py "
              f"--sample {sample} --lens {lens} --filt {filt} --force")


def process_lens_mask(lens, filt_dirs, sample, requested_filt, drizzle_pass, force,
                      broadcast=False, brush_radius=6, brush_width=None, display='snr',
                      subtract_radial=None, side_by_side=True,
                      size=cutout_paths.DEFAULT_SIZE, dataset_subplot=True,
                      propose_from='auto'):
    """Draw a mask once for one lens, on its best/forced band. By default the mask is written
    for THAT BAND ONLY; `broadcast=True` also writes it to the lens's other bands by WCS
    reprojection. `filt_dirs` maps filt -> cutout_dir from discover_targets; the mask FITS
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
    display_filt = pick_display_filt(filts, requested_filt)
    if display_filt is None:
        print(f"{lens}: requested --filt {requested_filt} not present "
              f"(have {', '.join(filts)}), skipping")
        return False

    display_dir = filt_dirs[display_filt]                # cutout dir of the draw band
    display_prefix = find_prefix(display_dir, drizzle_pass)
    if display_prefix is None:
        print(f"{lens} {display_filt}: no cutout sci for --pass {drizzle_pass}, skipping")
        return False

    # Skip if the draw band already has a mask, unless --force.
    if (os.path.exists(os.path.join(display_dir, f'{display_prefix}_mask.fits'))
            and not force):
        print(f"{lens}: mask already exists ({display_filt}), skipping (--force to redraw)")
        return False

    display_hdr = fits.getheader(cutout_paths.find_stamp(display_dir, display_prefix, 'sci'))
    proposal_from, proposal = find_proposal_mask(
        filt_dirs, display_filt, display_hdr, drizzle_pass, propose_from, force=force)
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
        proposal=proposal, proposal_label=f'{proposal_from} mask' if proposal_from else None)

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
        cutout_dir = filt_dirs[filt]
        prefix = find_prefix(cutout_dir, drizzle_pass)
        if prefix is None:
            continue
        band_hdr = fits.getheader(cutout_paths.find_stamp(cutout_dir, prefix, 'sci'))
        band_ps = pixel_scale_from_header(band_hdr)
        is_draw = (cutout_dir == display_dir)
        if is_draw:
            out_bool = scribbled
        else:
            out_bool = reproject_mask_bool(scribbled, src_wcs, WCS(band_hdr).celestial,
                                           (band_hdr['NAXIS2'], band_hdr['NAXIS1']))
        mask_path = os.path.join(cutout_dir, f'{prefix}_mask.fits')
        aplt.fits_array(array=al.Mask2D(mask=out_bool, pixel_scales=band_ps),
                        file_path=mask_path, overwrite=True)
        n_excl = int(np.asarray(out_bool).sum())
        print(f"  wrote {mask_path}  ({source_tag if is_draw else f'reprojected<-{display_filt}'}"
              f", {n_excl}px excluded)")
        entry = {
            'prefix': prefix,
            'drizzle_pass': 'cr' if prefix == 'cutout_cr' else 'nocrrej',
            'pixel_scale_arcsec': round(band_ps, 6),
            'display': source_label,
            # Which reduction the stamp was drawn on, from its own header (cutout_paths.py).
            'bcfill': bool(band_hdr.get('BCFILL', False)),
            'source': source_tag if is_draw else f'reprojected_from_{display_filt}',
            'n_excluded_px': n_excl,
        }
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
            regenerate_dataset_subplot(lens, filt, sample, size, drizzle_pass)
    return True


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--sample', default=mast_target_names.DEFAULT_SAMPLE,
                   help=f'sample subdirectory of data/cutouts/ (default '
                        f'{mast_target_names.DEFAULT_SAMPLE})')
    p.add_argument('--lens', default=None,
                   help='restrict to one lens; default every lens with cutouts in --sample')
    p.add_argument('--filt', default=None,
                   help='force which band you DRAW on (e.g. f606W); default the best available '
                        'band per lens (f814W>f606W>f555W>f160W>...). Only that band gets a '
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
    p.add_argument('--pass', dest='drizzle_pass', choices=['auto', 'cr', 'nocrrej'],
                   default='auto',
                   help="which cutout pass to mask, matching make_cutouts.py's --pass: "
                        "'auto' (default) prefers cutout_cr_*, falling back to cutout_* "
                        "where no CR pass exists (F160W)")
    p.add_argument('--size', type=float, default=cutout_paths.DEFAULT_SIZE,
                   help=f'cutout tree to draw masks for (default {cutout_paths.DEFAULT_SIZE:g}", '
                        'i.e. data/cutouts/ -- the tree this tool is meant for and '
                        'the only size-variant tree tracked in git, see .gitignore; pass '
                        '--size 20 for the untracked, regenerable 20" tree)')
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
                   help='also offer bands whose cutout directory is GITIGNORED. Off by '
                        'default: such a band is not a science product (the three 420 s '
                        'SLACS SNAP f814W diagnostics are why the rule exists), and hand-drawn '
                        'work written there is non-regenerable and invisible to every clone. '
                        'Use only to work on a diagnostic deliberately.')
    a = p.parse_args()

    # --broadcast writes one draw to every band unseen; --propose-from reviews another band's
    # mask in the band it is landing on. They are two answers to the same question, and
    # combining them would write a reviewed-for-f555W mask back over f814W's own.
    if a.broadcast and a.propose_from != 'none':
        p.error('--broadcast and --propose-from are alternatives: pass --propose-from none '
                'to broadcast one draw to every band unseen, or drop --broadcast to review '
                'the inherited mask per band')

    root = cutout_paths.cutouts_root(ws_path, a.size)

    # discover_targets yields per (lens, filt); regroup to filt -> cutout_dir per lens so we
    # know every band of a lens (filt=None: consider all bands; pick_display_filt picks one).
    lens_filts = {}
    for lens, filt, cutout_dir in discover_targets(root, a.sample, a.lens, None,
                                                   include_ignored=a.include_ignored):
        lens_filts.setdefault(lens, {})[filt] = cutout_dir

    if not lens_filts:
        raise SystemExit(f"no cutouts found under {root} for sample {a.sample} matching "
                         f"lens={a.lens!r}")

    print(f"{len(lens_filts)} lens(es) to process")
    made = skipped = 0
    for lens in sorted(lens_filts):
        if process_lens_mask(lens, lens_filts[lens], a.sample, a.filt, a.drizzle_pass,
                             a.force, a.broadcast, brush_radius=a.brush_radius,
                             brush_width=a.brush_width, display=a.display,
                             subtract_radial=a.subtract_radial,
                             side_by_side=a.side_by_side, size=a.size,
                             dataset_subplot=a.dataset_subplot,
                             propose_from=a.propose_from):
            made += 1
        else:
            skipped += 1

    print(f"\nDone: {made} lens(es) masked, {skipped} skipped")


if __name__ == '__main__':
    main()
