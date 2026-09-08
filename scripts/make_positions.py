#!/usr/bin/env python
"""
Interactive GUI for marking multiply-imaged (lensed) source positions, per lens.

For each lens in a sample this launches PyAutoLens's `al.Clicker` GUI over the lens's
best-available cutout -- the same tool as
autolens_workspace:scripts/imaging/data_preparation/gui/positions.py -- so you can
DOUBLE-CLICK each lensed image of the source (2 for a double, 4 for a quad, ...). Each
click snaps to the brightest pixel within `--search-box-size` and is recorded in
arcsec relative to the stamp centre; the collected clicks are saved as an
`al.Grid2DIrregular` in `cutout_[cr_]positions.json` (via `al.output_to_json`), the file a
modelling script loads to build a positions likelihood penalty (`al.PositionsLH`) that
rejects mass models mapping the images too far apart in the source plane.

WHY ONE BAND PER LENS, THEN BROADCAST. Image positions are band-INDEPENDENT sky
coordinates: every band is cut about the same shared centre (`make_cutouts.py
--center-band f814W`) and pinned to the same output WCS geometry (`final_rot=0`, tangent
point at the lens position), so a position measured in ARCSEC relative to the stamp centre
is identical in every band -- even where the pixel scales differ (f160W 0.06" vs 0.05" for
the optical bands), because arcsec, not pixels, is the shared frame PyAutoLens works in. So
you mark ONCE, on the highest-S/N band available (f814W > f606W > f555W > f160W > ...), and
the same Grid2DIrregular is broadcast into every band's cutout dir under that band's own
pass prefix, so any per-band modelling script finds a local positions.json next to its
sci/noise/psf/mask. `--filt` forces which band you mark on.

Trees mirror make_masks.py exactly: `--size` selects the cutout tree (default 12",
data/cutouts/, the git-tracked one -- and, like hand-drawn masks, these hand-marked
positions are a NON-regenerable product, so they belong in the tracked tree), and
`--variant {auto,bcfill,standard}` picks the reduction shown (bcfill and standard share
crop geometry EXACTLY, so positions marked on one are valid on the other; auto prefers
bcfill as the cleaner image). Each band's file is written into its PRIORITY tree only
(bcfill where a bcfill cutout exists, else standard) -- one file serves whichever reduction
is modelled, same convention as make_masks.py.

This is a manual, one-lens-at-a-time tool (not a batch driver): each lens blocks on its own
Tk window until you close it. A lens that already has positions is skipped so a run resumes
across a sample; --force re-marks.

Usage:
    uv run python scripts/make_positions.py --sample slacs_gold            # best band per lens
    uv run python scripts/make_positions.py --lens J0008-0004              # one lens, best band
    uv run python scripts/make_positions.py --lens J0008-0004 --filt f606W # force the band
    uv run python scripts/make_positions.py --sample slacs_gold --force    # re-mark everything
    uv run python scripts/make_positions.py --lens J0008-0004 --no-subtract-radial
        # opt OUT of the default radial-profile subtraction and click on the raw image
    uv run python scripts/make_positions.py --lens J0008-0004 --mask-center 0.8 --vmax-value 10
        # hide the central 0.8" of deflector light and saturate above S/N 10, so the arcs
        # are the brightest thing on screen
    (all three are DISPLAY ONLY -- the click snap uses true flux, saved positions are identical)
"""

import argparse
import os
import sys
import warnings

import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt

# autonerves prints a workspace-version-mismatch UserWarning on import in this repo -- see
# make_masks.py; silence it so it doesn't repeat per lens.
warnings.filterwarnings('ignore', category=UserWarning, module='autonerves')

import autolens as al

ws_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mast_target_names
import info_json
import cutout_paths
# Reuse make_masks' display/tree/discovery helpers so the two GUIs stay in lock-step
# (same pass-prefix rule, same variant-sharing, same stretched display).
import make_masks

POSITIONS_JSON = os.path.join(ws_path, 'info', 'lens_positions.json')

# "Best band to mark on" is shared with make_masks so the two tools pick the same primary
# band per lens (f814W>f606W>f555W>f160W>...).
pick_display_filt = make_masks.pick_display_filt


def save_overlay_png(sci_native, pixel_scales, positions, out_path,
                     stretch='asinh', vmin_percent=5.0, vmax_percent=99.5, asinh_a=0.1,
                     vmax_value=None):
    """QC overlay: the (stretched) science image with the marked positions scattered on
    top, in the same arcsec frame the Clicker used. Written per band so you can confirm the
    shared positions land on the lensed images in EVERY band, not just the marked one.

    `vmax_value` (--vmax-value) saturates the bright core here too, so the arcs the
    positions sit on are as visible in the QC PNG as they were in the GUI. The central
    blank (--mask-center) is deliberately NOT applied: the QC image should still show the
    deflector, so you can judge the positions against the whole system.
    """
    disp = make_masks.stretched_display(sci_native, pixel_scales, stretch=stretch,
                                         vmin_percent=vmin_percent, vmax_percent=vmax_percent,
                                         asinh_a=asinh_a, vmax_value=vmax_value)
    n_y, n_x = disp.shape_native
    hw = int(n_x / 2) * pixel_scales
    ext = [-hw, hw, -hw, hw]
    fig = plt.figure(figsize=(7, 7))
    plt.imshow(np.asarray(disp.native), cmap='jet', extent=ext)
    if len(positions) > 0:
        ys = [p[0] for p in positions]
        xs = [p[1] for p in positions]
        plt.scatter(xs, ys, marker='+', c='white', s=180, linewidths=1.5)
        for i, (y, x) in enumerate(positions):
            plt.annotate(str(i), (x, y), color='white', fontsize=9,
                         xytext=(4, 4), textcoords='offset points')
    plt.xlabel('arcsec'); plt.ylabel('arcsec')
    plt.tight_layout()
    fig.savefig(out_path, dpi=110)
    plt.close(fig)


def mark_positions_gui(sci_native, pixel_scales, lens, filt, search_box_size,
                       display='snr', display_base=None, stretch='asinh',
                       vmin_percent=5.0, vmax_percent=99.5, asinh_a=0.1,
                       mask_center=0.0, vmax_value=None, subtract_radial=True):
    """Run the Clicker GUI once and return the clicked positions as a list of (y, x) arcsec.

    `sci_native` (REAL flux) is handed to al.Clicker so its brightest-pixel snap uses true
    counts; `display_base` (raw signal, or S/N -- see load_display_base) is what we STRETCH
    and imshow, so bright cores and faint arcs are both visible while clicking. Display and
    Clicker arrays share shape + orientation + extent, so what you click is what gets snapped.
    Figure construction mirrors the canonical workspace positions.py (jet, arcsec extent,
    onclick connected directly) -- only the imshow'd array is swapped for the stretched copy.

    Two DISPLAY-ONLY levers make the source arcs readable under the deflector light, which
    otherwise sets the colour scale on its own:
      * `mask_center` (arcsec) blanks a disc at the stamp centre -- the deflector, since
        every stamp is recentred on it -- and drops those pixels from the percentile
        statistics, so the remaining pixels (the arcs) get the full colour range. The hidden
        region is outlined so you can see what is behind the blank.
      * `vmax_value` clips the display at an absolute value in the base array's units (S/N
        by default), saturating everything brighter.
      * `subtract_radial` removes the deflector's azimuthally-averaged radial profile,
        which is the strongest of the three -- it reaches the images buried INSIDE the
        galaxy envelope, where a central blank would hide them too.
    Neither touches `sci_native`, so al.Clicker still snaps on true flux. Do note the snap
    is flux-based: a click placed inside/adjacent to the blanked core can still snap onto a
    core pixel within `search_box_size`.
    """
    clicker = al.Clicker(image=sci_native, pixel_scales=pixel_scales,
                         search_box_size=search_box_size)

    base = display_base if display_base is not None else sci_native
    if subtract_radial:
        base = make_masks.radial_median_subtract(base)
    exclude = make_masks.central_disc(np.asarray(base).shape, pixel_scales, mask_center)
    disp = make_masks.stretched_display(base, pixel_scales, stretch=stretch,
                                        vmin_percent=vmin_percent, vmax_percent=vmax_percent,
                                        asinh_a=asinh_a, exclude=exclude,
                                        vmax_value=vmax_value)
    n_y, n_x = disp.shape_native
    hw = int(n_x / 2) * pixel_scales
    ext = [-hw, hw, -hw, hw]

    fig = plt.figure(figsize=(12, 12))
    ax = plt.gca()
    plt.imshow(np.asarray(disp.native), cmap='jet', extent=ext)
    plt.colorbar()
    if exclude is not None:
        # outline the blanked deflector so it is obvious what the black disc is hiding
        ax.add_patch(plt.Circle((0.0, 0.0), mask_center, fill=False, color='white',
                                lw=1.0, ls='--'))
    plt.title(f"{lens}  {filt}  --  DOUBLE-CLICK each lensed image, then close the window")
    plt.xlabel('arcsec'); plt.ylabel('arcsec')
    cid = fig.canvas.mpl_connect('button_press_event', clicker.onclick)
    plt.show()
    fig.canvas.mpl_disconnect(cid)
    plt.close(fig)

    return list(clicker.click_list)  # each entry (y_arcsec, x_arcsec)


def process_lens(lens, filt_dirs, sample, requested_filt, drizzle_pass, force,
                 search_box_size, display, stretch, vmin_percent, vmax_percent, asinh_a,
                 mask_center=0.0, vmax_value=None, subtract_radial=True):
    """Mark positions once for one lens and broadcast the result to every band.

    `filt_dirs` maps filt -> write_dirs (an ordered list of (variant, cutout_dir), bcfill
    before standard) for that lens, from make_masks.discover_targets. Returns True if
    positions were written.
    """
    filts = sorted(filt_dirs)
    display_filt = pick_display_filt(filts, requested_filt)
    if display_filt is None:
        print(f"{lens}: requested --filt {requested_filt} not present "
              f"(have {', '.join(filts)}), skipping")
        return False

    display_dir = filt_dirs[display_filt][0][1]          # priority-tree dir of marked band
    display_prefix = make_masks.find_prefix(display_dir, drizzle_pass)
    if display_prefix is None:
        print(f"{lens} {display_filt}: no cutout sci for --pass {drizzle_pass}, skipping")
        return False

    # Skip if the marked band already carries positions (one file per lens), unless --force.
    marked_json = os.path.join(display_dir, f'{display_prefix}_positions.json')
    if os.path.exists(marked_json) and not force:
        print(f"{lens}: positions already exist ({display_filt}), skipping (--force to re-mark)")
        return False

    sci_path = os.path.join(display_dir, f'{display_prefix}_sci.fits')
    with fits.open(sci_path) as hdul:
        sci_hdr = hdul[0].header
    pixel_scales = make_masks.pixel_scale_from_header(sci_hdr)
    sci_native = al.Array2D.from_fits(file_path=sci_path, pixel_scales=pixel_scales).native
    base, source_label = make_masks.load_display_base(display_dir, display_prefix,
                                                      display, pixel_scales)

    print(f"\n{lens}  marking on {display_filt}  [{display_prefix} pass]  "
          f"pixel_scale={pixel_scales:.4f}\"/pix")
    print(f"  {sci_path}")
    print(f"  display: {source_label}, {stretch} stretch  |  snap search box "
          f"= {search_box_size}px")
    if mask_center and mask_center > 0:
        print(f"  central {mask_center:g}\" blanked in the DISPLAY (deflector light hidden "
              f"and dropped from the stretch percentiles); snap still uses true flux")
    if vmax_value is not None:
        print(f"  display saturated above {vmax_value:g} ({source_label})")
    if subtract_radial:
        print("  deflector's radial-median profile subtracted from the DISPLAY "
              "(finding aid only -- expect a quadrupole residual on an elliptical)")
    print("  DOUBLE-CLICK each lensed image of the source (2 for a double, 4 for a quad).")
    print("  Each click snaps to the brightest nearby pixel. Close the window when done.")

    positions = mark_positions_gui(
        sci_native, pixel_scales, lens, display_filt, search_box_size,
        display=display, display_base=base, stretch=stretch,
        vmin_percent=vmin_percent, vmax_percent=vmax_percent, asinh_a=asinh_a,
        mask_center=mask_center, vmax_value=vmax_value,
        subtract_radial=subtract_radial)

    if len(positions) == 0:
        print(f"  no positions clicked -- nothing written for {lens}")
        return False
    if len(positions) < 2:
        print(f"  WARNING: only {len(positions)} position clicked -- a lensed source has "
              f">=2 images. Re-run with --force if this was a mis-click.")

    grid = al.Grid2DIrregular(values=positions)

    # Broadcast the SAME arcsec Grid2DIrregular to every band, ONE file per band into that
    # band's PRIORITY tree only (filt_dirs[filt][0] -- bcfill where a bcfill cutout exists,
    # else standard), matching make_masks.py exactly: bcfill and standard share crop geometry,
    # so one file serves whichever reduction is modelled, and downstream reads the priority
    # tree. Each under that dir's own pass prefix (f160W has no CR pass -> cutout_positions.json).
    written = []
    for filt in filts:
        variant, cutout_dir = filt_dirs[filt][0]         # priority tree for this band
        prefix = make_masks.find_prefix(cutout_dir, drizzle_pass)
        if prefix is None:
            continue
        json_path = os.path.join(cutout_dir, f'{prefix}_positions.json')
        al.output_to_json(obj=grid, file_path=json_path)
        # QC overlay in each band, in that band's own arcsec frame.
        band_sci = os.path.join(cutout_dir, f'{prefix}_sci.fits')
        band_ps = make_masks.pixel_scale_from_header(fits.getheader(band_sci))
        band_native = al.Array2D.from_fits(file_path=band_sci, pixel_scales=band_ps).native
        save_overlay_png(band_native, band_ps, positions,
                         os.path.join(cutout_dir, f'{prefix}_positions.png'),
                         stretch=stretch, vmin_percent=vmin_percent,
                         vmax_percent=vmax_percent, asinh_a=asinh_a,
                         vmax_value=vmax_value)
        written.append({'filt': filt, 'variant': variant or 'standard', 'prefix': prefix})
        print(f"  wrote {json_path}")

    info_json.update(POSITIONS_JSON, sample, lens, display_filt, {
        'n_positions': len(positions),
        'positions_arcsec': [[round(y, 6), round(x, 6)] for (y, x) in positions],
        'marked_filt': display_filt,
        'marked_prefix': display_prefix,
        'drizzle_pass': 'cr' if display_prefix == 'cutout_cr' else 'nocrrej',
        'pixel_scale_arcsec': round(pixel_scales, 6),
        'search_box_size': search_box_size,
        'display': source_label,
        'display_mask_center_arcsec': mask_center or None,
        'display_vmax_value': vmax_value,
        'display_subtract_radial': bool(subtract_radial),
        'broadcast_to': written,
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
                   help='force which band you mark on (e.g. f606W); default the best '
                        'available band per lens (f814W>f606W>f555W>f160W>...). The marked '
                        'positions are broadcast to every band regardless')
    p.add_argument('--pass', dest='drizzle_pass', choices=['auto', 'cr', 'nocrrej'],
                   default='auto',
                   help="which cutout pass to mark on, matching make_cutouts.py/make_masks.py: "
                        "'auto' (default) prefers cutout_cr_*, falling back to cutout_* where "
                        "no CR pass exists (F160W)")
    p.add_argument('--size', type=float, default=cutout_paths.DEFAULT_SIZE,
                   help=f'cutout tree to mark positions in (default {cutout_paths.DEFAULT_SIZE:g}", '
                        'data/cutouts/ -- the git-tracked tree; hand-marked positions are '
                        'non-regenerable, like masks, so they belong here)')
    p.add_argument('--variant', choices=['auto', 'bcfill', 'standard'], default='auto',
                   help="which reduction's cutouts to show: 'auto' (default) marks on the "
                        "bcfill sci where it exists (cleaner image; shared geometry), else "
                        "standard; 'bcfill'/'standard' restrict to one tree. Positions are "
                        "broadcast to every variant dir regardless")
    p.add_argument('--force', action='store_true', default=False,
                   help='re-mark a lens that already has positions (default: skip it)')
    p.add_argument('--search-box-size', type=int, default=5,
                   help='half-width in pixels of the brightest-pixel snap around each click '
                        '(al.Clicker search_box_size; default 5). Set 0-ish to take the click '
                        'as-is')
    p.add_argument('--display', choices=['snr', 'sci'], default='snr',
                   help="image to click on: 'snr' (default) = signal/noise (tracks "
                        "significance given the non-uniform noise map); 'sci' = raw signal. "
                        "Snap always uses true flux; choice does not affect saved positions")
    p.add_argument('--stretch', choices=list(make_masks._STRETCHES), default='asinh',
                   help='intensity stretch for the DISPLAYED image only (asinh default shows '
                        'bright cores and faint arcs together)')
    p.add_argument('--asinh-a', type=float, default=0.1,
                   help='asinh softening (smaller = stronger stretch; --stretch asinh only)')
    p.add_argument('--vmin-percent', type=float, default=5.0,
                   help='lower clip percentile for the display stretch (default 5.0)')
    p.add_argument('--vmax-percent', type=float, default=99.5,
                   help='upper clip percentile for the display stretch (default 99.5)')
    p.add_argument('--mask-center', type=float, default=0.0, metavar='ARCSEC',
                   help='blank a disc of this radius (arcsec) at the stamp centre in the '
                        'DISPLAYED image, to hide the deflector light that otherwise drowns '
                        'the source arcs (default 0 = off; ~0.5-1.0 is typical). The hidden '
                        'pixels are also dropped from the stretch percentiles, so the arcs '
                        'get the full colour range. Display only -- the click snap still '
                        'uses true flux and the saved positions are unaffected')
    p.add_argument('--vmax-value', type=float, default=None,
                   help='saturate the DISPLAYED image above this absolute value, in the '
                        'units of --display (S/N by default, e.g. 10; raw flux with '
                        '--display sci). Overrides --vmax-percent, and is applied to the QC '
                        'overlay PNG as well. Display only')
    p.add_argument('--subtract-radial', action=argparse.BooleanOptionalAction, default=True,
                   help="subtract the deflector's azimuthally-averaged radial profile from "
                        'the DISPLAYED image (DEFAULT ON -- it is what makes the images '
                        'clickable at all on most lenses). The strongest of the three '
                        'arc-finding levers: unlike --mask-center it reveals images buried '
                        'inside the galaxy envelope rather than hiding them. Expect a '
                        'quadrupole residual (the deflector is elliptical, not circular) '
                        'and treat it as a finding aid, not photometry. Display only -- pass '
                        '--no-subtract-radial to click on the raw image instead')
    a = p.parse_args()

    # Same variant->tree resolution as make_masks (bcfill first for display preference).
    if a.variant == 'bcfill':
        variants = ['bcfill']
    elif a.variant == 'standard':
        variants = ['']
    else:  # auto
        variants = ['bcfill', '']
    trees = [(v, cutout_paths.cutouts_root(ws_path, a.size, variant=v)) for v in variants]

    # discover_targets yields per (lens, filt); regroup to filt -> write_dirs per lens so we
    # mark once and broadcast across the lens's bands.
    lens_filts = {}
    for lens, filt, _display_dir, write_dirs in make_masks.discover_targets(
            trees, a.sample, a.lens, None):   # filt=None: consider all bands, pick best below
        lens_filts.setdefault(lens, {})[filt] = write_dirs

    if not lens_filts:
        roots = ', '.join(r for _, r in trees)
        raise SystemExit(f"no cutouts found under [{roots}] for sample {a.sample} matching "
                         f"lens={a.lens!r}")

    print(f"{len(lens_filts)} lens(es) to process")
    made = skipped = 0
    for lens in sorted(lens_filts):
        if process_lens(lens, lens_filts[lens], a.sample, a.filt, a.drizzle_pass, a.force,
                        a.search_box_size, a.display, a.stretch,
                        a.vmin_percent, a.vmax_percent, a.asinh_a,
                        a.mask_center, a.vmax_value, a.subtract_radial):
            made += 1
        else:
            skipped += 1

    print(f"\nDone: {made} lens(es) marked, {skipped} skipped")


if __name__ == '__main__':
    main()
