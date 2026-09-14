#!/usr/bin/env python
"""
Interactive GUI ARC mask-making: paint the arcs / multiple images, keep only those.

The second of the two hand-drawn mask tools, and the OPPOSITE polarity to `make_masks.py`:

    make_masks.py       paint what to REMOVE  -> cutout_[cr_]mask.fits       (contaminants)
    make_arc_masks.py   paint the ARCS to KEEP -> cutout_[cr_]mask_arcs.fits (this script)

Here you scribble ONLY the lensed source -- the arcs and the multiple images -- and
everything you do not paint is masked out, deflector light included. The product is a mask
that isolates the source: the region a source-only fit, an arc S/N measurement, or a
source-plane analysis should look at.

WHAT IS ON DISK, AND ITS POLARITY. The saved array is in `al.Mask2D`'s own convention --
True = EXCLUDED from the fit -- exactly like `cutout_[cr_]mask.fits`, so both files load
through the same reader and can never be mixed up by polarity:

    mask = al.Mask2D.from_fits(file_path=".../cutout_cr_mask_arcs.fits", pixel_scales=0.05)
    # mask == True  everywhere EXCEPT the arcs   -> fitting under it fits only the arcs
    # ~mask         is the arc region you painted (what MASKARC/NARCPX count)

So what you PAINT is inverted before writing (`arc_region = painted`, `saved = ~painted`).
The header records this explicitly: MASKTYPE='ARCS', MASKPOL, NARCPX. `make_masks.py`'s
scribble, by contrast, is saved verbatim -- there the painted region IS the excluded region.

THE DISPLAY SUBTRACTS THE DEFLECTOR BY DEFAULT (`--subtract-radial`, on here, off in
make_masks.py). The arcs are the subject of this mask and on most lenses they are invisible
under the galaxy envelope until its azimuthally-averaged radial profile is removed. It is a
DISPLAY transform only -- the mask is read back as brush pixel positions, so it cannot
change what a stroke masks -- and a finding aid, not photometry (an elliptical leaves a
quadrupole residual; an arc partly self-subtracts at its own radius). Pass
--no-subtract-radial to draw against the raw image.

ALREADY-MARKED POSITIONS ARE SHOWN. Where `make_positions.py` has recorded
cutout_[cr_]positions.json for the band, each marked image is ringed in the display so the
arc mask can be drawn around the same multiple images. --no-show-positions turns it off.

Everything else -- ONE DRAW PER BAND by default (--broadcast opts into sharing the draw
across the lens's bands by WCS reprojection), the band-priority order, --size tree routing,
the skip/--force behaviour -- is `make_masks.py`'s, imported from it rather than
re-implemented, so the two tools stay in lock-step. Where --broadcast IS used, note the
reprojection is applied to the ARC REGION (True = arc), not to the saved array:
outside-footprint pixels then default to 'not arc', which is the safe side -- reprojecting
the inverted array would instead leave stamp-edge slivers unmasked.

Usage:
    uv run python scripts/make_arc_masks.py --sample slacs_gold             # best band per lens
    uv run python scripts/make_arc_masks.py --lens J0330-0020               # one lens
    uv run python scripts/make_arc_masks.py --lens J0330-0020 --filt f606W  # force the band
    uv run python scripts/make_arc_masks.py --lens J0330-0020 --no-subtract-radial --force
    uv run python scripts/make_arc_masks.py --lens J0330-0020 --broadcast  # share one draw
"""

import argparse
import os
import sys
import warnings

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

warnings.filterwarnings('ignore', category=UserWarning, module='autonerves')

import autolens as al
import autolens.plot as aplt
from matplotlib import pyplot as plt

ws_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mast_target_names
import info_json
import cutout_paths
# Every display/tree/discovery/broadcast helper is make_masks'. This script contributes the
# opposite polarity, the arc-specific display defaults, and its own product name.
import make_masks

ARC_MASKS_JSON = os.path.join(ws_path, 'info', 'lens_arc_masks.json')

ARC_SUFFIX = 'mask_arcs'          # cutout_[cr_]mask_arcs.fits, beside cutout_[cr_]mask.fits


def load_positions(cutout_dir, prefix):
    """Return the band's marked image positions as a list of (y, x) arcsec, or [] if
    make_positions.py has not run for it. Read-only -- this tool never writes positions."""
    path = os.path.join(cutout_dir, f'{prefix}_positions.json')
    if not os.path.exists(path):
        return []
    try:
        return [(float(p[0]), float(p[1])) for p in al.from_json(file_path=path)]
    except Exception as exc:                       # a hand-edited/partial file must not block
        print(f"  NOTE: could not read {os.path.basename(path)} ({exc}); no position rings")
        return []


def ring_overlay(positions, ring_radius_px=6):
    """Return an `overlay(disp_array, pixel_scales)` callable for make_masks.draw_mask_gui
    that rings each marked image position in the DISPLAYED array.

    Rings are burned in at the display MINIMUM (dark), which reads clearly against both the
    mid-tone background and the bright arcs of a radial-subtracted image -- a maximum-valued
    ring would be indistinguishable from arc flux, which is the one thing it must not hide.
    Display only: the Scribbler returns brush positions, so this cannot enter the mask.
    """
    def _overlay(disp_array, pixel_scales):
        a = np.array(disp_array.native, dtype=float)
        n_y, n_x = a.shape
        cy, cx = (n_y - 1) / 2.0, (n_x - 1) / 2.0
        yy, xx = np.mgrid[0:n_y, 0:n_x]
        for (y_as, x_as) in positions:
            # arcsec -> pixel. +x arcsec is +column, but +y arcsec is -ROW: autoarray's
            # native grid puts row 0 at the TOP (+y), verified against
            # al.Grid2D.uniform(...).native[0, 0] == (+5.975, -5.975) for a 240px/0.05"
            # stamp. Getting this sign wrong silently rings the mirror image of each
            # position -- it lands plausibly near the lens and on nothing.
            py, px = cy - y_as / pixel_scales, cx + x_as / pixel_scales
            r = np.hypot(yy - py, xx - px)
            a[(r >= ring_radius_px - 0.8) & (r <= ring_radius_px + 0.8)] = float(a.min())
        return al.Array2D.no_mask(values=a, pixel_scales=pixel_scales).native
    return _overlay


def save_arc_qc_png(cutout_dir, prefix, arc_bool, pixel_scales, out_path, lens, filt,
                    stretch='asinh', vmin_percent=5.0, vmax_percent=99.5, asinh_a=0.1):
    """QC overlay: the arc region outlined over the band's own (radial-subtracted) image.

    Written per band, so a broadcast can be checked where it actually matters -- the f160W
    regrid (0.06"/px), where an array-index copy would land the region ~20px off.
    """
    base, _ = make_masks.load_display_base(cutout_dir, prefix, 'snr', pixel_scales)
    disp = make_masks.stretched_display(make_masks.radial_median_subtract(base), pixel_scales,
                                        stretch=stretch, vmin_percent=vmin_percent,
                                        vmax_percent=vmax_percent, asinh_a=asinh_a)
    a = np.asarray(disp.native)
    hw = int(a.shape[1] / 2) * pixel_scales
    ext = [-hw, hw, -hw, hw]
    fig = plt.figure(figsize=(7, 7))
    plt.imshow(a, cmap='jet', extent=ext)
    # origin='upper' is REQUIRED: with an extent, imshow draws row 0 at the top but
    # contour defaults to row 0 at the BOTTOM, so without it the outline is drawn
    # vertically mirrored -- a QC image that looks fine and checks nothing (verified).
    plt.contour(np.asarray(arc_bool, dtype=float), levels=[0.5], colors='white',
                linewidths=1.2, extent=ext, origin='upper')
    plt.title(f'{lens}  {filt}  --  arc mask region (white outline)')
    plt.xlabel('arcsec'); plt.ylabel('arcsec')
    plt.tight_layout()
    fig.savefig(out_path, dpi=110)
    plt.close(fig)


def write_arc_mask(arc_bool, pixel_scales, path):
    """Write the arc mask in al.Mask2D convention (True = EXCLUDED, so only the painted arc
    region survives) and stamp the polarity into the header, since this file's meaning is
    the inverse of what was painted."""
    excluded = ~np.asarray(arc_bool, dtype=bool)
    aplt.fits_array(array=al.Mask2D(mask=excluded, pixel_scales=pixel_scales),
                    file_path=path, overwrite=True)
    with fits.open(path, mode='update') as hdul:
        hdr = hdul[0].header
        hdr['MASKTYPE'] = ('ARCS', 'hand-drawn lensed-source (arc/multiple-image) mask')
        hdr['MASKPOL'] = ('True=excluded', 'al.Mask2D convention; ~mask is the arc region')
        hdr['NARCPX'] = (int(np.asarray(arc_bool).sum()), 'pixels inside the arc region')


def process_lens_arc_mask(lens, filt_dirs, sample, requested_filt, drizzle_pass, force,
                          broadcast, brush_radius, brush_width, display, stretch,
                          vmin_percent, vmax_percent, asinh_a, subtract_radial,
                          side_by_side, show_positions, ring_radius_px):
    """Draw the arc mask once for one lens, on its best/forced band only unless
    `broadcast` is set (then also WCS-reprojected to the lens's other bands). Mirrors
    make_masks.process_lens_mask; differs only in polarity, product name, and the arc-specific
    display (radial subtraction + position rings). Returns True if a mask was written."""
    filts = sorted(filt_dirs)
    display_filt = make_masks.pick_display_filt(filts, requested_filt)
    if display_filt is None:
        print(f"{lens}: requested --filt {requested_filt} not present "
              f"(have {', '.join(filts)}), skipping")
        return False

    display_dir = filt_dirs[display_filt]
    display_prefix = make_masks.find_prefix(display_dir, drizzle_pass)
    if display_prefix is None:
        print(f"{lens} {display_filt}: no cutout sci for --pass {drizzle_pass}, skipping")
        return False

    # Skip if the draw band already has an arc mask, unless --force.
    if (os.path.exists(os.path.join(display_dir, f'{display_prefix}_{ARC_SUFFIX}.fits'))
            and not force):
        print(f"{lens}: arc mask already exists ({display_filt}), skipping (--force to redraw)")
        return False

    positions = load_positions(display_dir, display_prefix) if show_positions else []
    overlay = ring_overlay(positions, ring_radius_px) if positions else None
    if positions:
        print(f"  {len(positions)} marked position(s) ringed in the display "
              f"(from {display_prefix}_positions.json)")

    prompt = [
        "  Scribbler GUI: paint ONLY the arcs / multiple images of the source.",
        "    Everything you do NOT paint is masked OUT (the deflector included) -- this is",
        "    the OPPOSITE polarity to make_masks.py, which paints what to remove.",
        "    keys '1' = GREEN brush, paint arc   |   '2' = RED brush, un-paint (erase)",
    ]
    (painted, erased, src_ps, src_hdr, brush_width, start_radius,
     source_label) = make_masks.draw_mask_gui(
        display_dir, display_prefix, lens, display_filt, brush_radius, brush_width,
        display, stretch, vmin_percent, vmax_percent, asinh_a,
        subtract_radial=subtract_radial, side_by_side=side_by_side,
        prompt=prompt, overlay=overlay)

    # painted = the arcs to KEEP; the red brush trims an over-painted stroke back off it.
    arc_region = np.asarray(painted, dtype=bool) & ~np.asarray(erased, dtype=bool)
    if not arc_region.any():
        print(f"  nothing painted -- no arc mask written for {lens} (an empty arc region "
              f"would mask out the entire stamp)")
        return False
    src_wcs = WCS(src_hdr).celestial

    targets = (filts if broadcast else [display_filt])
    for filt in targets:
        cutout_dir = filt_dirs[filt]
        prefix = make_masks.find_prefix(cutout_dir, drizzle_pass)
        if prefix is None:
            continue
        band_hdr = fits.getheader(os.path.join(cutout_dir, f'{prefix}_sci.fits'))
        band_ps = make_masks.pixel_scale_from_header(band_hdr)
        is_draw = (cutout_dir == display_dir)
        # Reproject the ARC REGION, never the inverted array: outside-footprint pixels then
        # default to 'not arc' (masked out), the safe side of the stamp edge.
        band_arc = arc_region if is_draw else make_masks.reproject_mask_bool(
            arc_region, src_wcs, WCS(band_hdr).celestial,
            (band_hdr['NAXIS2'], band_hdr['NAXIS1']))
        mask_path = os.path.join(cutout_dir, f'{prefix}_{ARC_SUFFIX}.fits')
        write_arc_mask(band_arc, band_ps, mask_path)
        save_arc_qc_png(cutout_dir, prefix, band_arc, band_ps,
                        os.path.join(cutout_dir, f'{prefix}_{ARC_SUFFIX}.png'), lens, filt,
                        stretch=stretch, vmin_percent=vmin_percent,
                        vmax_percent=vmax_percent, asinh_a=asinh_a)
        n_arc = int(np.asarray(band_arc).sum())
        print(f"  wrote {mask_path}  ({'drawn' if is_draw else f'reprojected<-{display_filt}'}, "
              f"{n_arc}px in the arc region, everything else masked out)")
        entry = {
            'prefix': prefix,
            'drizzle_pass': 'cr' if prefix == 'cutout_cr' else 'nocrrej',
            'pixel_scale_arcsec': round(band_ps, 6),
            'display': source_label,
            # Which reduction the stamp was drawn on, from its own header (cutout_paths.py).
            'bcfill': bool(band_hdr.get('BCFILL', False)),
            'source': 'drawn' if is_draw else f'reprojected_from_{display_filt}',
            'n_arc_px': n_arc,
            'polarity': 'saved True=excluded (al.Mask2D); ~mask = arc region',
        }
        if is_draw:
            entry['brush_width'] = round(brush_width, 6)
            entry['brush_radius_px'] = start_radius
            entry['n_positions_shown'] = len(positions)
            entry['subtract_radial'] = bool(subtract_radial)
        info_json.update(ARC_MASKS_JSON, sample, lens, filt, entry)
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
                        'band per lens (f814W>f606W>f555W>f160W>...). Only that band gets an '
                        'arc mask unless --broadcast is passed')
    p.add_argument('--broadcast', action=argparse.BooleanOptionalAction, default=False,
                   help="ALSO write the arc mask to the lens's other bands, WCS-reprojected "
                        'onto each grid. Default OFF: run again with --filt for the next band. '
                        'An arc is often detected in one filter and not another, and its usable '
                        'extent differs with each band\'s depth and PSF, so the drawn region is '
                        'not in practice band-independent')
    p.add_argument('--pass', dest='drizzle_pass', choices=['auto', 'cr', 'nocrrej'],
                   default='auto',
                   help="which cutout pass to mask, matching make_masks.py/make_cutouts.py: "
                        "'auto' (default) prefers cutout_cr_*, falling back to cutout_* where "
                        "no CR pass exists (F160W)")
    p.add_argument('--size', type=float, default=cutout_paths.DEFAULT_SIZE,
                   help=f'cutout tree to draw arc masks for (default '
                        f'{cutout_paths.DEFAULT_SIZE:g}", data/cutouts/ -- the git-tracked '
                        'tree; like every hand-drawn product these are non-regenerable)')
    p.add_argument('--force', action='store_true', default=False,
                   help='redraw an arc mask that already exists (default: skip it)')
    p.add_argument('--brush-radius', type=int, default=4,
                   help='starting brush radius in PIXELS (default 4 -- finer than make_masks.py, '
                        "since arcs are thin). Resize live in the GUI with '='/'-'")
    p.add_argument('--brush-width', type=float, default=None,
                   help='explicit brush width as a FRACTION of image height; overrides '
                        '--brush-radius when given')
    p.add_argument('--display', choices=['snr', 'sci'], default='snr',
                   help="image to scribble on: 'snr' (default) = signal/noise, which tracks "
                        "statistical significance given the non-uniform noise map; 'sci' = raw "
                        'signal. Choice does not affect the saved mask')
    p.add_argument('--stretch', choices=list(make_masks._STRETCHES), default='asinh',
                   help='intensity stretch for the DISPLAYED image only (default asinh)')
    p.add_argument('--asinh-a', type=float, default=0.1,
                   help='asinh softening (smaller = stronger stretch; --stretch asinh only)')
    p.add_argument('--vmin-percent', type=float, default=5.0,
                   help='lower clip percentile for the display stretch (default 5.0)')
    p.add_argument('--vmax-percent', type=float, default=99.5,
                   help='upper clip percentile for the display stretch (default 99.5)')
    p.add_argument('--subtract-radial', action=argparse.BooleanOptionalAction, default=True,
                   help="subtract the deflector's azimuthally-averaged radial profile from the "
                        'DISPLAYED image (DEFAULT ON here -- the arcs are this mask\'s subject '
                        'and are usually invisible under the galaxy envelope without it). '
                        'Display only; --no-subtract-radial draws against the raw image')
    p.add_argument('--side-by-side', action=argparse.BooleanOptionalAction, default=True,
                   help='show BOTH views at once -- radial-subtracted left, as-observed right '
                        '-- and accept scribbles on either panel, combined onto the one mask. '
                        'The subtracted panel is where the arcs are visible; the as-observed '
                        'one is where their real extent against the galaxy envelope is '
                        '(default on)')
    p.add_argument('--show-positions', action=argparse.BooleanOptionalAction, default=True,
                   help='ring the images already marked by make_positions.py in the display, '
                        'as a guide for where the multiple images are (default on; silently '
                        'inactive for a band with no positions file)')
    p.add_argument('--ring-radius', type=int, default=6,
                   help='radius in pixels of those position rings (default 6)')
    a = p.parse_args()

    root = cutout_paths.cutouts_root(ws_path, a.size)

    lens_filts = {}
    for lens, filt, cutout_dir in make_masks.discover_targets(root, a.sample, a.lens, None):
        lens_filts.setdefault(lens, {})[filt] = cutout_dir

    if not lens_filts:
        raise SystemExit(f"no cutouts found under {root} for sample {a.sample} matching "
                         f"lens={a.lens!r}")

    print(f"{len(lens_filts)} lens(es) to process")
    made = skipped = 0
    for lens in sorted(lens_filts):
        if process_lens_arc_mask(lens, lens_filts[lens], a.sample, a.filt, a.drizzle_pass,
                                 a.force, a.broadcast, a.brush_radius, a.brush_width,
                                 a.display, a.stretch, a.vmin_percent, a.vmax_percent,
                                 a.asinh_a, a.subtract_radial, a.side_by_side,
                                 a.show_positions, a.ring_radius):
            made += 1
        else:
            skipped += 1

    print(f"\nDone: {made} lens(es) arc-masked, {skipped} skipped")


if __name__ == '__main__':
    main()
