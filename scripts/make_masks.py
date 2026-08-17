#!/usr/bin/env python
"""
Interactive GUI mask-making, cycling through a sample's cutouts.

For each (lens, filt) cutout on disk, launches PyAutoLens's `Scribbler` GUI over the
cutout's science image -- the same tool as
autolens_workspace:scripts/imaging/data_preparation/gui/mask.py -- so you can scribble the
areas to REMOVE from the fit (painted = excluded, everything unpainted = kept; this is the
opposite polarity to the autolens_workspace GUI, which scribbles the region to keep), then
writes the result as cutout_[cr_]mask.fits alongside that cutout's sci/noise/psf products.

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

This is a manual, one-image-at-a-time tool (not a batch driver): each cutout blocks on its
own Tk window until you press Esc. Already-masked cutouts are skipped so a run can be
resumed across lenses; --force redraws.

Usage:
    uv run python scripts/make_masks.py --sample slacs_gold
    uv run python scripts/make_masks.py --lens J0008-0004 --filt f814W
    uv run python scripts/make_masks.py --sample slacs_gold --filt f814W   # one band, every lens
    uv run python scripts/make_masks.py --sample slacs_gold --force
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
                      vmin_percent=5.0, vmax_percent=99.5, asinh_a=0.1):
    """Return an al.Array2D holding a contrast-stretched COPY of `data_native`, for
    display in the Scribbler only.

    The mask is built from the pixel positions of scribble circles, never from pixel
    values, so remapping the displayed values cannot change the saved mask -- this is
    purely so bright cores and faint/problematic structure are both visible while
    scribbling. Values are clipped to the [vmin_percent, vmax_percent] percentiles then
    mapped through `stretch` onto [0, 1] (a plain linear imshow of that IS the stretch).
    """
    vals = np.asarray(data_native, dtype=float)
    finite = vals[np.isfinite(vals)]
    if finite.size == 0:
        return al.Array2D.no_mask(values=vals, pixel_scales=pixel_scales).native
    vmin, vmax = np.percentile(finite, [vmin_percent, vmax_percent])
    if not vmax > vmin:                      # flat/degenerate frame -- avoid /0 in interval
        vmax = vmin + (abs(vmin) or 1.0)
    st = AsinhStretch(a=asinh_a) if stretch == 'asinh' else _STRETCHES[stretch]()
    disp = st(ManualInterval(vmin=vmin, vmax=vmax)(vals))  # clip -> stretch -> [0, 1]
    return al.Array2D.no_mask(values=np.nan_to_num(disp), pixel_scales=pixel_scales).native


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


def make_mask_for(display_dir, write_dirs, lens, filt, sample, drizzle_pass, force,
                  brush_radius=6, brush_width=None, display='snr', stretch='asinh',
                  vmin_percent=5.0, vmax_percent=99.5, asinh_a=0.1):
    """Run the Scribbler GUI once and write ONE mask FITS, into the priority tree only
    (`write_dirs[0]` -- bcfill where its cutout exists, else standard).

    The trees share crop geometry exactly, so a single mask serves whichever reduction is
    modelled; it is not duplicated across trees. Skips (draws nothing) if a mask already
    exists in ANY of `write_dirs`, unless --force. Returns True if a mask was written.
    """
    display_prefix = find_prefix(display_dir, drizzle_pass)
    if display_prefix is None:
        print(f"{lens} {filt}: no cutout sci file for --pass {drizzle_pass}, skipping")
        return False

    # Existing masks anywhere across the trees -> skip (one mask suffices for both).
    existing = []
    for variant, cutout_dir in write_dirs:
        prefix = find_prefix(cutout_dir, drizzle_pass)
        if prefix is None:
            continue
        if os.path.exists(os.path.join(cutout_dir, f'{prefix}_mask.fits')):
            existing.append(variant or 'standard')
    if existing and not force:
        print(f"{lens} {filt}: mask already exists in [{', '.join(existing)}], "
              f"skipping (--force to redraw)")
        return False

    # Draw on, and write to, the priority tree only.
    write_variant = write_dirs[0][0] or 'standard'
    mask_path = os.path.join(display_dir, f'{display_prefix}_mask.fits')
    sci_path = os.path.join(display_dir, f'{display_prefix}_sci.fits')
    with fits.open(sci_path) as hdul:
        sci_hdr = hdul[0].header
    pixel_scales = pixel_scale_from_header(sci_hdr)

    print(f"\n{lens} {filt}  [{display_prefix} pass, {write_variant}]  "
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
    # Scribbling here means REMOVE: the painted region is excluded from the fit, everything
    # unpainted is kept. Scribbler marks the painted region True, and al.Mask2D's convention
    # already has True = excluded from the fit, so the scribbled array IS the mask -- no
    # inversion. (This is the opposite polarity from
    # autolens_workspace:scripts/imaging/data_preparation/gui/mask.py, which scribbles the
    # region to KEEP and inverts.)
    mask = al.Mask2D(mask=scribbled, pixel_scales=pixel_scales)

    aplt.fits_array(array=mask, file_path=mask_path, overwrite=True)
    print(f"  wrote {mask_path}")

    info_json.update(MASKS_JSON, sample, lens, filt, {
        'prefix': display_prefix,
        'drizzle_pass': 'cr' if display_prefix == 'cutout_cr' else 'nocrrej',
        'pixel_scale_arcsec': round(pixel_scales, 6),
        'brush_width': round(brush_width, 6),
        'brush_radius_px': start_radius,
        'display': source_label,
        'variant': write_variant,
    })
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
                   help='restrict to one filter; default every filter cutout the lens has')
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
                   help="which reduction's cutouts to mask: 'auto' (default) draws one mask "
                        "on the bcfill sci where it exists (cleaner image, shared geometry) "
                        "else the standard sci, skipping if a mask exists in EITHER tree; "
                        "'bcfill' uses only data/cutouts_bcfill/; 'standard' uses only "
                        "data/cutouts/")
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

    targets = list(discover_targets(trees, a.sample, a.lens, a.filt))
    if not targets:
        roots = ', '.join(r for _, r in trees)
        raise SystemExit(f"no cutouts found under [{roots}] for sample {a.sample} matching "
                         f"lens={a.lens!r} filt={a.filt!r}")

    print(f"{len(targets)} lens/filter cutout(s) to process")
    made = skipped = 0
    for lens, filt, display_dir, write_dirs in targets:
        if make_mask_for(display_dir, write_dirs, lens, filt, a.sample, a.drizzle_pass,
                         a.force, brush_radius=a.brush_radius, brush_width=a.brush_width,
                         display=a.display, stretch=a.stretch,
                         vmin_percent=a.vmin_percent, vmax_percent=a.vmax_percent,
                         asinh_a=a.asinh_a):
            made += 1
        else:
            skipped += 1

    print(f"\nDone: {made} mask(s) drawn, {skipped} skipped")


if __name__ == '__main__':
    main()
