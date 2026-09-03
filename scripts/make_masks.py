#!/usr/bin/env python
"""
Interactive GUI mask-making, ONE draw per lens, broadcast (WCS-reprojected) to every band.

For each lens in a sample, launches PyAutoLens's `Scribbler` GUI over the lens's
best-available cutout -- the same tool as
autolens_workspace:scripts/imaging/data_preparation/gui/mask.py -- so you can scribble the
areas to REMOVE from the fit (painted = excluded, everything unpainted = kept; this is the
opposite polarity to the autolens_workspace GUI, which scribbles the region to keep), then
writes the result as cutout_[cr_]mask.fits alongside that cutout's sci/noise/psf products.

ONE DRAW PER LENS, THEN BROADCAST. A mask marks a sky region (the deflector+arcs to keep,
contaminants to exclude), which is the same physical region in every band. So you draw ONCE,
on the highest-S/N band available (f814W > f606W > f555W > f160W > ...), and the mask is
broadcast to the lens's other bands. The broadcast is a rigorous per-pixel WCS reprojection
(nearest-neighbour, source sci-header WCS -> target sci-header WCS), NOT an array-index copy:
the 0.05" optical bands (f814W/f606W/f555W) share a grid to <1px so it is nearly identity,
but f160W is a genuinely different grid (0.06"/px, 200px vs 240px for a 12" stamp -- a common
sky point sits ~20px away by index), and only reprojection places the mask correctly there.
`--filt` forces which band you draw on; `--no-broadcast` writes only that one band (the escape
hatch for a mask that must genuinely differ per band, e.g. a contaminant bright in only one
filter).

Defaults to the pipeline's standard 12" cutout tree (data/cutouts/, cutout_paths.DEFAULT_SIZE)
-- the tree these hand-drawn masks are meant for and the only size-variant tree tracked in
git (see .gitignore) precisely because it carries them, which no script can regenerate.
Pass --size 20 to mask the (untracked, regenerable) 20" tree instead.

bcfill variant (--variant, default 'auto'). The bad-column-filled re-drizzle
(data/cutouts_bcfill/, also tracked) shares crop geometry with the standard tree EXACTLY
-- identical NAXIS/CRPIX/CRVAL -- so a mask drawn on one is pixel-valid on the other, and
bcfill differs only by having its dead-column noise stripes filled (a cleaner image to
scribble on). It covers ACS (f814W/f555W) + WFPC2 (f606W) only -- there is no bcfill for
f160W (WFC3/IR has no bad columns) or the gallery sample. So the default is:
    --variant auto  (default): per (lens, filt), a single mask is drawn on and written to
        the bcfill cutout if it exists, else the standard cutout. Because the geometry is
        shared, that one mask serves whichever reduction is modelled -- it is NOT duplicated
        across trees. A (lens, filt) is SKIPPED if a mask already exists in EITHER tree
        (bcfill or standard); with --force it is redrawn, again into the priority tree.
    --variant bcfill: bcfill tree only (skip lenses/filters with no bcfill cutout).
    --variant standard: standard tree only (the historic behaviour, data/cutouts/).

This is a manual, one-lens-at-a-time tool (not a batch driver): each lens blocks on its own
Tk window until you press Esc. A lens that already has a mask is skipped so a run resumes
across a sample; --force redraws.

Usage:
    uv run python scripts/make_masks.py --sample slacs_gold                # best band per lens
    uv run python scripts/make_masks.py --lens J0008-0004                  # one lens, best band
    uv run python scripts/make_masks.py --lens J0008-0004 --filt f606W     # force the draw band
    uv run python scripts/make_masks.py --lens J0008-0004 --filt f160W --no-broadcast --force
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

ws_path = '/Users/samlange/Code/HST_lenses_data_reduction'
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

    def on_mouse_motion(self, event):
        if event.inaxes != self.ax or event.xdata is None or event.ydata is None:
            return
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


def discover_targets(trees, sample, lens=None, filt=None):
    """Yield (lens, filt, display_dir, write_dirs) per lens/filt that has a cutout in any
    of `trees`, sorted for a reproducible run order.

    `trees` is an ordered list of (variant, cutouts_root_dir) in display-preference order
    (bcfill before standard for --variant auto). For each (lens, filt) present in any tree:
      - write_dirs   = every tree's cutout_dir that has a matching sci file (both trees
                       share crop geometry exactly, so one mask is valid for all of them),
                       tagged with its variant;
      - display_dir  = the first such tree in preference order (the image scribbled on).
    """
    seen = {}
    for variant, root in trees:
        pattern = os.path.join(root, sample, lens or '*', filt or '*')
        for cutout_dir in sorted(glob.glob(pattern)):
            if find_prefix(cutout_dir, 'auto') is None:
                continue
            this_filt = os.path.basename(cutout_dir)
            this_lens = os.path.basename(os.path.dirname(cutout_dir))
            seen.setdefault((this_lens, this_filt), []).append((variant, cutout_dir))
    for (this_lens, this_filt) in sorted(seen):
        write_dirs = seen[(this_lens, this_filt)]        # already in preference order
        display_variant, display_dir = write_dirs[0]
        yield this_lens, this_filt, display_dir, write_dirs


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


def draw_mask_gui(display_dir, display_prefix, lens, filt, brush_radius, brush_width,
                  display, stretch, vmin_percent, vmax_percent, asinh_a):
    """Run the Scribbler GUI once on the chosen band. Returns
    (scribbled_bool, pixel_scales, sci_header, brush_width, start_radius, source_label).
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
    print("  Scribbler GUI: scribble the areas to REMOVE from the fit (contaminants, field")
    print("    sources, bad regions). Everything you do NOT paint is kept in the fit.")
    print("    keys:  '=' bigger brush | '-' smaller brush | 'z' undo last | Esc/q done")

    base, source_label = load_display_base(display_dir, display_prefix, display, pixel_scales)
    disp_array = stretched_display(base, pixel_scales, stretch=stretch,
                                   vmin_percent=vmin_percent, vmax_percent=vmax_percent,
                                   asinh_a=asinh_a)
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
    print(f"  brush start radius = {start_radius} px")
    scribbler = _GuardedScribbler(image=disp_array, brush_width=brush_width)
    scribbled = scribbler.show_mask()
    return scribbled, pixel_scales, sci_hdr, brush_width, start_radius, source_label


def process_lens_mask(lens, filt_dirs, sample, requested_filt, drizzle_pass, force,
                      no_broadcast, brush_radius=6, brush_width=None, display='snr',
                      stretch='asinh', vmin_percent=5.0, vmax_percent=99.5, asinh_a=0.1):
    """Draw a mask once for one lens (on its best/forced band) and broadcast it to every band
    by WCS reprojection. `filt_dirs` maps filt -> write_dirs (ordered (variant, cutout_dir),
    bcfill before standard) from discover_targets. One mask FITS is written per band into that
    band's priority tree (the two variants share a grid, so one serves both -- the skip check
    covers either). Returns True if a mask was written.
    """
    filts = sorted(filt_dirs)
    display_filt = pick_display_filt(filts, requested_filt)
    if display_filt is None:
        print(f"{lens}: requested --filt {requested_filt} not present "
              f"(have {', '.join(filts)}), skipping")
        return False

    display_dir = filt_dirs[display_filt][0][1]          # priority-tree dir of the draw band
    display_prefix = find_prefix(display_dir, drizzle_pass)
    if display_prefix is None:
        print(f"{lens} {display_filt}: no cutout sci for --pass {drizzle_pass}, skipping")
        return False

    # Skip if the draw band already has a mask in EITHER variant, unless --force (one mask per
    # band serves both variants; matches the historic skip-if-either behaviour).
    existing = []
    for variant, cutout_dir in filt_dirs[display_filt]:
        prefix = find_prefix(cutout_dir, drizzle_pass)
        if prefix and os.path.exists(os.path.join(cutout_dir, f'{prefix}_mask.fits')):
            existing.append(variant or 'standard')
    if existing and not force:
        print(f"{lens}: mask already exists ({display_filt} in [{', '.join(existing)}]), "
              f"skipping (--force to redraw)")
        return False

    scribbled, src_ps, src_hdr, brush_width, start_radius, source_label = draw_mask_gui(
        display_dir, display_prefix, lens, display_filt, brush_radius, brush_width,
        display, stretch, vmin_percent, vmax_percent, asinh_a)
    src_wcs = WCS(src_hdr).celestial

    # Broadcast: write one mask per band into that band's priority tree. The draw band uses the
    # scribbled array verbatim (identity); every other band is WCS-reprojected onto its own grid
    # (near-identity for the 0.05" optical bands, a real ~20px regrid for f160W). --no-broadcast
    # restricts to the draw band only (the escape hatch for a genuinely band-specific mask).
    targets = ([display_filt] if no_broadcast else filts)
    for filt in targets:
        variant, cutout_dir = filt_dirs[filt][0]         # priority tree for this band
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
        print(f"  wrote {mask_path}  ({'drawn' if is_draw else f'reprojected<-{display_filt}'}, "
              f"{n_excl}px excluded)")
        entry = {
            'prefix': prefix,
            'drizzle_pass': 'cr' if prefix == 'cutout_cr' else 'nocrrej',
            'pixel_scale_arcsec': round(band_ps, 6),
            'display': source_label,
            'variant': variant or 'standard',
            'source': 'drawn' if is_draw else f'reprojected_from_{display_filt}',
            'n_excluded_px': n_excl,
        }
        if is_draw:
            entry['brush_width'] = round(brush_width, 6)
            entry['brush_radius_px'] = start_radius
        info_json.update(MASKS_JSON, sample, lens, filt, entry)
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
                        'band per lens (f814W>f606W>f555W>f160W>...). The mask is broadcast to '
                        'every band regardless (unless --no-broadcast)')
    p.add_argument('--no-broadcast', action='store_true', default=False,
                   help='write the mask only to the drawn band, not reprojected to the others '
                        '-- the escape hatch for a mask that must differ per band (e.g. a '
                        'contaminant bright in only one filter). Pair with --filt and --force '
                        'to refine one band without touching the rest')
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
    p.add_argument('--variant', choices=['auto', 'bcfill', 'standard'], default='auto',
                   help="which reduction's cutouts to draw on: 'auto' (default) draws on the "
                        "bcfill sci where it exists (cleaner image, shared geometry) else the "
                        "standard sci, skipping if a mask exists in EITHER tree; 'bcfill' uses "
                        "only data/cutouts_bcfill/; 'standard' uses only data/cutouts/")
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
    a = p.parse_args()

    # Ordered by display preference: bcfill first (cleaner image), standard second.
    if a.variant == 'bcfill':
        variants = ['bcfill']
    elif a.variant == 'standard':
        variants = ['']
    else:  # auto
        variants = ['bcfill', '']
    trees = [(v, cutout_paths.cutouts_root(ws_path, a.size, variant=v)) for v in variants]

    # discover_targets yields per (lens, filt); regroup to filt -> write_dirs per lens so we
    # draw once and broadcast across the lens's bands (filt=None: consider all bands, pick best).
    lens_filts = {}
    for lens, filt, _display_dir, write_dirs in discover_targets(trees, a.sample, a.lens, None):
        lens_filts.setdefault(lens, {})[filt] = write_dirs

    if not lens_filts:
        roots = ', '.join(r for _, r in trees)
        raise SystemExit(f"no cutouts found under [{roots}] for sample {a.sample} matching "
                         f"lens={a.lens!r}")

    print(f"{len(lens_filts)} lens(es) to process")
    made = skipped = 0
    for lens in sorted(lens_filts):
        if process_lens_mask(lens, lens_filts[lens], a.sample, a.filt, a.drizzle_pass,
                             a.force, a.no_broadcast, brush_radius=a.brush_radius,
                             brush_width=a.brush_width, display=a.display, stretch=a.stretch,
                             vmin_percent=a.vmin_percent, vmax_percent=a.vmax_percent,
                             asinh_a=a.asinh_a):
            made += 1
        else:
            skipped += 1

    print(f"\nDone: {made} lens(es) masked, {skipped} skipped")


if __name__ == '__main__':
    main()
