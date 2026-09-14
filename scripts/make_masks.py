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
'1' GREEN adds to the mask, '2' RED erases from it, scribbling on whichever panel you like --
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

ws_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mast_target_names
import info_json
import cutout_paths

MASKS_JSON = os.path.join(ws_path, 'info', 'lens_masks.json')


class _GuardedScribbler(al.Scribbler):
    """al.Scribbler with two behavioural fixes for this pipeline's masking.

    1. Motion handler ignores events outside the image axes, and forces a redraw so
       the brush circle always renders. Upstream `on_mouse_motion` (a) unconditionally
       sets the brush-circle centre to (event.xdata, event.ydata), which are None when
       the cursor is over non-image space (padding, colour-bar gutter, off-figure) --
       the first such event raises `TypeError: unsupported operand type(s) for +:
       'float' and 'NoneType'` inside matplotlib's `add_patch` (harmless -- the mask is
       built from mouse-down scribbles, not motion events -- but it spams a traceback);
       and (b) on the reposition path (`self.brush.center = center`) never asks the
       canvas to redraw, so the blue brush circle only paints when some *other* event
       forces a draw. On a fresh figure an idle-draw happens to service it, but from the
       SECOND lens onward -- after show_mask()'s `plt.ioff()` and a figure close -- no
       such draw fires, so the brush patch exists but is never painted and only the OS
       cursor shows (scribbling still works, since mouse-down draws itself). Skipping
       out-of-axes motion suppresses (a); a `draw_idle()` on every in-axes motion fixes
       (b) so the circle tracks the cursor on every lens.

    2. MID-RUN brush resizing is multiplicative, not a fixed +/-2px step. Press '='
       to grow and '-' to shrink the brush (upstream keybindings, unchanged); each
       press scales the radius by ~1.4x (at least +/-1px), so you can drop to a 1px
       fine brush for detail and ramp up to a thick brush for the outskirts in a few
       presses instead of dozens. The radius floor drops from 4px to 1px, and the
       current radius is echoed to the console on every change.
    """

    _RESIZE_FACTOR = 1.4
    _MIN_RADIUS = 1

    #: set before construction to title the figure (used to label side-by-side panels).
    #: al.Scribbler builds its figure and then BLOCKS inside __init__, so there is no
    #: post-construction hook -- the title has to go on from inside, hence this class attr.
    title = None

    def __init__(self, *args, **kwargs):
        title, type(self).title = type(self).title, None
        self._pending_title = title
        super().__init__(*args, **kwargs)

    def on_mouse_motion(self, event):
        if event.inaxes != self.ax or event.xdata is None or event.ydata is None:
            return
        pending = getattr(self, '_pending_title', None)
        if pending is not None:                  # first in-axes motion: label the panels
            self._pending_title = None
            self.ax.set_title(pending, fontsize=11)
        super().on_mouse_motion(event)
        # Upstream repositions the brush patch but never redraws on the motion path, so
        # the circle can stop painting on later lenses (see class docstring). Force a
        # coalesced redraw so it always tracks the cursor. draw_idle (not draw) keeps
        # this cheap -- one paint per event-loop idle, not one per motion event.
        self.figure.canvas.draw_idle()

    def _apply_radius(self, radius):
        self.brush_radius = max(self._MIN_RADIUS, int(round(radius)))
        if self.brush is not None:
            self.brush.radius = self.brush_radius
            self.figure.canvas.draw()
        print(f"  brush radius = {self.brush_radius} px", flush=True)

    def enlarge_brush(self):
        # at least +1px, and faster (x1.4) the bigger it already is
        self._apply_radius(max(self.brush_radius + 1,
                               self.brush_radius * self._RESIZE_FACTOR))

    def shrink_brush(self):
        # at least -1px, so it always moves even at small radii
        self._apply_radius(min(self.brush_radius - 1,
                               self.brush_radius / self._RESIZE_FACTOR))


#: blank columns between the two side-by-side display panels (see draw_mask_gui).
_PANEL_GAP = 6

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
    display in the Scribbler/Clicker only.

    The mask is built from the pixel positions of scribble circles (and positions from
    click coordinates), never from pixel values, so remapping the displayed values cannot
    change the saved product -- this is purely so bright cores and faint/problematic
    structure are both visible while scribbling. Values are clipped to the
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


def radial_median_subtract(data_native):
    """Subtract the azimuthally-averaged (median) radial profile about the stamp centre.

    A DISPLAY transform for finding lensed arcs: the deflector is a smooth, near-circular
    elliptical, so its light is almost entirely a function of radius, while the arcs are
    not -- subtracting the per-radius median removes the galaxy and leaves the source
    images standing out, including the ones buried inside the galaxy envelope that a
    central blank cannot reach without also hiding them. It is not a galaxy FIT: real
    ellipticity leaves a quadrupole residual, and the arcs themselves bias the median at
    their own radius (self-subtraction), so the result is a finding aid only -- never a
    photometric product. Clicks snap on the untouched flux array regardless.
    """
    a = np.asarray(data_native, dtype=float)
    n_y, n_x = a.shape
    yy, xx = np.mgrid[0:n_y, 0:n_x]
    r_bin = np.hypot(yy - (n_y - 1) / 2.0, xx - (n_x - 1) / 2.0).astype(int)
    prof = np.zeros(r_bin.max() + 1)
    for i in range(prof.size):
        ring = a[r_bin == i]
        ring = ring[np.isfinite(ring)]
        if ring.size:
            prof[i] = np.median(ring)
    return a - prof[r_bin]


def mask_boundary(m):
    """1-pixel inner boundary of a boolean mask, for outlining a proposed mask in the display.

    Inner (`m & ~erosion(m)`), so the outline lies ON masked pixels -- it therefore survives
    the proposal panel, where the mask interior is blanked to the display floor. `border_value
    =1` treats outside-the-array as masked, so a mask running off the stamp edge is not
    outlined along the frame itself.
    """
    from scipy.ndimage import binary_erosion
    m = np.asarray(m, dtype=bool)
    if not m.any():
        return m
    return m & ~binary_erosion(m, np.ones((3, 3), dtype=bool), border_value=1)


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
    """Pick 'cutout_cr' or 'cutout', mirroring make_cutouts.py's --pass auto logic: prefer
    the CR-rejected pass, fall back to no-CR (e.g. F160W, which has no CR pass). Returns
    None if the requested pass has no sci file in this cutout_dir.
    """
    has_cr = os.path.exists(os.path.join(cutout_dir, 'cutout_cr_sci.fits'))
    has_nocr = os.path.exists(os.path.join(cutout_dir, 'cutout_sci.fits'))
    if drizzle_pass == 'cr':
        return 'cutout_cr' if has_cr else None
    if drizzle_pass == 'nocrrej':
        return 'cutout' if has_nocr else None
    return 'cutout_cr' if has_cr else ('cutout' if has_nocr else None)


def discover_targets(root, sample, lens=None, filt=None):
    """Yield (lens, filt, cutout_dir) per band with a cutout under `root` -- the dir both
    scribbled on and written to -- sorted for a reproducible run order.
    """
    found = {}
    pattern = os.path.join(root, sample, lens or '*', filt or '*')
    for cutout_dir in sorted(glob.glob(pattern)):
        if find_prefix(cutout_dir, 'auto') is None:
            continue
        this_filt = os.path.basename(cutout_dir)
        this_lens = os.path.basename(os.path.dirname(cutout_dir))
        found[(this_lens, this_filt)] = cutout_dir
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
        file_path=os.path.join(display_dir, f'{prefix}_sci.fits'),
        pixel_scales=pixel_scales).native
    if display == 'sci':
        return sci, 'signal'
    noise_path = os.path.join(display_dir, f'{prefix}_noise.fits')
    if not os.path.exists(noise_path):
        print(f"  NOTE: no {prefix}_noise.fits -- displaying signal instead of S/N")
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
# f606W is the gallery/BELLS primary; the rest follow by typical arc S/N.
_BAND_PRIORITY = ['f814W', 'f606W', 'f606W_v2', 'f555W', 'f160W', 'f438W', 'f275W', 'f225W']


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
            src_hdr = fits.getheader(os.path.join(cutout_dir, f'{prefix}_sci.fits'))
            return cand, reproject_mask_bool(src_bool, WCS(src_hdr).celestial,
                                             dst_wcs, dst_shape)
    return None, None


def confirm_proposal(lens, filt, source_filt, n_prop, n_add, n_erase, n_final, n_drawn):
    """Ask what to do with a reviewed proposal, AFTER the GUI has been closed -- so the
    decision is made having actually seen the mask over this band's own image.

    Three outcomes, and an empty answer means the common one (apply):
      apply  -- write proposal + green additions - red erasures (the reviewed mask);
      drawn  -- reject the proposal outright and keep only what was painted green;
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


#: positional words used to label the composited panels in the figure title.
_PANEL_POSITIONS = {1: [''], 2: ['LEFT', 'RIGHT'], 3: ['LEFT', 'MIDDLE', 'RIGHT'],
                    4: ['1st', '2nd', '3rd', '4th']}


def fold_panels(scribbled, n_panels, panel_n_x):
    """Fold a scribble drawn on the horizontally-composited display back onto one band grid.

    Every panel shows the SAME sky at the same pixel grid, so a stroke on any of them means
    the same thing: the per-panel slices are OR-ed together. Panels are laid out at a stride
    of `panel_n_x + _PANEL_GAP` (see draw_mask_gui). Returns (folded, px_per_panel).
    """
    scribbled = np.asarray(scribbled, dtype=bool)
    if n_panels <= 1:
        return scribbled, [int(scribbled.sum())]
    step = panel_n_x + _PANEL_GAP
    out = np.zeros((scribbled.shape[0], panel_n_x), dtype=bool)
    per_panel = []
    for i in range(n_panels):
        s = scribbled[:, i * step: i * step + panel_n_x]
        per_panel.append(int(s.sum()))
        out |= s
    return out, per_panel


def draw_mask_gui(display_dir, display_prefix, lens, filt, brush_radius, brush_width,
                  display, stretch, vmin_percent, vmax_percent, asinh_a,
                  subtract_radial=False, side_by_side=True, prompt=None, overlay=None,
                  proposal=None, proposal_label='proposed mask'):
    """Run the Scribbler GUI once on the chosen band. Returns
    (painted, erased, pixel_scales, sci_header, brush_width, start_radius, source_label).

    TWO BRUSHES. al.Scribbler carries two independent scribble segments; this pipeline uses
    them as ADD (segment '1', GREEN, the default) and ERASE (segment '2', RED) -- press '1'
    and '2' in the GUI to switch. `painted` is the green scribble and `erased` the red one;
    the caller decides what they mean (make_masks: mask = (proposal | painted) & ~erased).
    Erasing is what makes an inherited mask editable rather than all-or-nothing.

    `proposal` (bool array on this band's grid, optional) is a mask inherited from another
    band, shown for review as a 1-px OUTLINE burned into every panel -- where its edge falls
    against this band's own structure being the whole question. Only an outline: a filled-in
    "mask applied" panel was tried and dropped, since blanking the interior showed nothing the
    boundary does not already carry while hiding the pixels being judged.

    `subtract_radial` shows the deflector's radial-median-subtracted image (see
    radial_median_subtract) -- DISPLAY ONLY, the scribble is read back as pixel positions so
    it cannot change what a given brush stroke masks. Off by default here (a contaminant mask
    is drawn against the real sky), on by default in make_arc_masks.py, where the arcs are
    the whole subject.

    `side_by_side` (default True, and only meaningful with subtract_radial) shows BOTH views
    at once -- subtracted on the left, as-observed on the right, separated by a blank gutter
    -- because each answers a different question: the subtracted panel is where the arcs are
    visible at all, the as-observed panel is where a contaminant's real extent and the
    galaxy envelope are. You may scribble on EITHER panel: the two halves are read back and
    UNIONed onto the single-band mask, so a stroke on the right lands at the same sky
    position as the same stroke on the left. Each panel is stretched independently (their
    dynamic ranges differ by orders of magnitude, so a shared scale would flatten one).

    `prompt` replaces the printed instruction lines (so the arc tool can state its own,
    opposite, polarity) and `overlay` is a callable applied to each stretched panel just
    before they are composited.
    The scribbled array marks True = REMOVE from the fit (painted = excluded, unpainted =
    kept -- al.Mask2D's own convention, so the scribbled array IS the mask, no inversion;
    opposite polarity to autolens_workspace's mask.py, which scribbles the region to keep).
    """
    sci_path = os.path.join(display_dir, f'{display_prefix}_sci.fits')
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
    print("    keys:  '=' bigger brush | '-' smaller brush | 'z' undo last | Esc/q done")

    base, source_label = load_display_base(display_dir, display_prefix, display, pixel_scales)
    # Which view(s) to show. Extra panels only when there is genuinely something to compare.
    panels = []
    if subtract_radial:
        panels.append(('radial-subtracted', radial_median_subtract(base)))
    if not subtract_radial or side_by_side:
        panels.append(('as-observed', base))

    def _panel(values):
        d = stretched_display(values, pixel_scales, stretch=stretch, vmin_percent=vmin_percent,
                              vmax_percent=vmax_percent, asinh_a=asinh_a)
        return np.asarray((overlay(d, pixel_scales) if overlay is not None else d).native)

    rendered = [_panel(v) for _, v in panels]
    panel_names = [name for name, _ in panels]
    title_extra = ''
    if proposal is not None:
        # OUTLINE ONLY, on every panel -- deliberately no extra "mask applied" panel. Blanking
        # the interior showed nothing the boundary does not already carry, and it hid the very
        # pixels you are judging the mask against while costing every other panel a third of
        # the window. Where the edge falls against this band's structure is the question.
        edge = mask_boundary(np.asarray(proposal, dtype=bool))
        for r in rendered:
            r[edge] = 1.0
        title_extra = f'   --  {proposal_label} OUTLINED'

    panel_n_x = rendered[0].shape[1]
    if len(rendered) == 1:
        composite, title = rendered[0], None
    else:
        gutter = np.zeros((rendered[0].shape[0], _PANEL_GAP))
        parts = []
        for i, r in enumerate(rendered):
            if i:
                parts.append(gutter)
            parts.append(r)
        composite = np.hstack(parts)
        pos = _PANEL_POSITIONS.get(len(rendered),
                                   [str(i + 1) for i in range(len(rendered))])
        title = '   |   '.join(f'{p}: {n}' for p, n in zip(pos, panel_names))
        title += title_extra + '   --  scribble on any panel'
        if proposal is not None:
            title += "   ('1' green = ADD, '2' red = ERASE)"
    disp_array = al.Array2D.no_mask(values=composite, pixel_scales=pixel_scales).native
    source_label += ', ' + ' + '.join(panel_names)
    if proposal is not None:
        source_label += f' ({proposal_label} outlined)'
    # al.Scribbler sizes the brush as int(image_height * brush_width), a fraction. To get a
    # stamp-independent default (--brush-radius px) we derive the fraction from this stamp's
    # own height; an explicit --brush-width fraction (if given) overrides it.
    n_pix = disp_array.shape_native[0]
    if brush_width is None:
        brush_width = brush_radius / n_pix
    start_radius = int(n_pix * brush_width)
    print(f"  display: {source_label}, {stretch} stretch"
          + (f" (a={asinh_a})" if stretch == 'asinh' else '')
          + f", clip [{vmin_percent:g}, {vmax_percent:g}] percentile")
    if len(rendered) > 1:
        print(f"  {len(rendered)} PANELS side by side -- "
              + ", ".join(f"{p} {n}" for p, n in
                          zip(_PANEL_POSITIONS.get(len(rendered), []), panel_names)) + ".")
        print(f"    Scribble on any of them; the panels are combined onto the one mask "
              f"(same sky position).")
    if proposal is not None:
        print(f"  the proposed mask ({proposal_label}) is OUTLINED on every panel.")
        print(f"    keys '1' = GREEN brush, ADD to the mask   |   "
              f"'2' = RED brush, ERASE from it")
    print(f"  brush start radius = {start_radius} px")
    _GuardedScribbler.title = title
    scribbler = _GuardedScribbler(image=disp_array, brush_width=brush_width)
    scribbler.show_mask()             # returns the ADD segment; also does the plt.ioff()
    segments = scribbler.get_scribble_masks()     # both segments: '1' = add, '2' = erase

    n_panels = len(rendered)
    painted, add_per_panel = fold_panels(segments.get('1'), n_panels, panel_n_x)
    erased, _ = fold_panels(segments.get('2'), n_panels, panel_n_x)
    if n_panels > 1 and sum(1 for n in add_per_panel if n) > 1:
        print("  combined panels (" + " + ".join(f"{n}px {p}" for n, p in
              zip(add_per_panel, _PANEL_POSITIONS.get(n_panels, []))) + ")")
    return painted, erased, pixel_scales, sci_hdr, brush_width, start_radius, source_label


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
                      stretch='asinh', vmin_percent=5.0, vmax_percent=99.5, asinh_a=0.1,
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
    add (green) / erase (red) parts of it before choosing to apply, reject, or skip. It is the deliberate alternative to `--broadcast`, which writes the
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

    display_hdr = fits.getheader(os.path.join(display_dir, f'{display_prefix}_sci.fits'))
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

    (painted, erased, src_ps, src_hdr, brush_width, start_radius,
     source_label) = draw_mask_gui(
        display_dir, display_prefix, lens, display_filt, brush_radius, brush_width,
        display, stretch, vmin_percent, vmax_percent, asinh_a,
        subtract_radial=subtract_radial, side_by_side=side_by_side,
        proposal=proposal, proposal_label=f'{proposal_from} mask' if proposal_from else None)

    painted = np.asarray(painted, dtype=bool)
    erased = np.asarray(erased, dtype=bool)
    drawn_only = painted & ~erased
    n_add, n_erase = int(painted.sum()), int(erased.sum())

    # Nothing painted = you skipped this lens: write NO file. An all-False mask would be a
    # legitimate "exclude nothing" mask, which is exactly the ambiguity to avoid -- an absent
    # file means "not done yet" and the next run re-offers the lens, where an empty one would
    # silently mark it finished. The same rule holds after a review: an approved mask that
    # ends up empty (everything erased) is not written either.
    if proposal is None:
        scribbled, source_tag = drawn_only, 'drawn'
    else:
        reviewed = (proposal | painted) & ~erased
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
        band_hdr = fits.getheader(os.path.join(cutout_dir, f'{prefix}_sci.fits'))
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
                        'every panel, and you add (green) / erase '
                        '(red) parts of it before choosing to apply, reject, or skip -- the '
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
                        "morphology cross-check. Choice does not affect the saved mask")
    p.add_argument('--stretch', choices=list(_STRETCHES), default='asinh',
                   help='intensity stretch for the DISPLAYED image only (does not affect the '
                        'saved mask): asinh (default) shows bright cores and faint structure '
                        'together; log/sqrt/linear also available')
    p.add_argument('--asinh-a', type=float, default=0.1,
                   help='asinh softening parameter (smaller = stronger stretch, more faint '
                        'detail; only used with --stretch asinh; default 0.1)')
    p.add_argument('--vmin-percent', type=float, default=5.0,
                   help='lower clip percentile for the display stretch (default 5.0)')
    p.add_argument('--vmax-percent', type=float, default=99.5,
                   help='upper clip percentile for the display stretch; lower it to brighten '
                        'and reveal problematic bright areas (default 99.5)')
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
    for lens, filt, cutout_dir in discover_targets(root, a.sample, a.lens, None):
        lens_filts.setdefault(lens, {})[filt] = cutout_dir

    if not lens_filts:
        raise SystemExit(f"no cutouts found under {root} for sample {a.sample} matching "
                         f"lens={a.lens!r}")

    print(f"{len(lens_filts)} lens(es) to process")
    made = skipped = 0
    for lens in sorted(lens_filts):
        if process_lens_mask(lens, lens_filts[lens], a.sample, a.filt, a.drizzle_pass,
                             a.force, a.broadcast, brush_radius=a.brush_radius,
                             brush_width=a.brush_width, display=a.display, stretch=a.stretch,
                             vmin_percent=a.vmin_percent, vmax_percent=a.vmax_percent,
                             asinh_a=a.asinh_a, subtract_radial=a.subtract_radial,
                             side_by_side=a.side_by_side, size=a.size,
                             dataset_subplot=a.dataset_subplot,
                             propose_from=a.propose_from):
            made += 1
        else:
            skipped += 1

    print(f"\nDone: {made} lens(es) masked, {skipped} skipped")


if __name__ == '__main__':
    main()
