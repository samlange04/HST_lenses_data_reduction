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

WHAT IS ON SCREEN. Two panels by default -- the deflector's radial profile subtracted on
the LEFT, the band as observed on the RIGHT -- because each answers a different question: the
subtracted panel is where an image buried in the galaxy envelope shows up at all, the
as-observed panel is where you judge whether that blob is really there or is a subtraction
artefact. They are separate axes over the same arcsec extent, so DOUBLE-CLICK EITHER and the
snap is identical. Anything this band's `cutout_[cr_]mask.fits` already calls a contaminant is
blanked (and outlined) on both, and dropped from the stretch and the radial profile: that
verdict is already in, and a masked neighbour otherwise owns the colour scale, invites a
mis-click, and stamps a dark ring across the arcs at its own radius. All of it is display
only -- the snap runs on true flux and the saved positions do not depend on any of it.

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
positions are a NON-regenerable product, so they belong in the tracked tree). There is one
cutout dir per band (cutout_paths.py), so each band's file lands beside the sci it was
marked on, with nothing to choose and no second copy to keep in step.

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
    uv run python scripts/make_positions.py --lens J0008-0004 --no-side-by-side
        # one subtracted panel instead of the default subtracted-LEFT/as-observed-RIGHT pair
    uv run python scripts/make_positions.py --sample slacs_gold --overlays-only
        # rebuild the gitignored QC PNGs from the tracked positions JSONs -- no GUI
    uv run python scripts/make_positions.py --lens J0008-0004 --no-hide-masked
        # show the contaminants this band's mask already excludes, instead of blanking them
    (all three are DISPLAY ONLY -- the click snap uses true flux, saved positions are identical)
"""

import argparse
import os
import sys
import warnings

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
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
# (same pass-prefix rule, same tree discovery, same stretched display).
import make_masks

POSITIONS_JSON = os.path.join(ws_path, 'info', 'lens_positions.json')

# "Best band to mark on" is shared with make_masks so the two tools pick the same primary
# band per lens (f814W>f606W>f555W>f160W>...).
pick_display_filt = make_masks.pick_display_filt


def arcsec_to_pixel(y, x, shape, pixel_scales):
    """(y, x) arcsec from the stamp centre -> (row, col). Row 0 is at +y (autoarray's
    native convention), which is the orientation trap recorded in AGENTS.md."""
    n_y, n_x = shape
    return (n_y - 1) / 2.0 - y / pixel_scales, (n_x - 1) / 2.0 + x / pixel_scales


def pixel_to_arcsec(row, col, shape, pixel_scales):
    """Inverse of arcsec_to_pixel."""
    n_y, n_x = shape
    return ((n_y - 1) / 2.0 - row) * pixel_scales, (col - (n_x - 1) / 2.0) * pixel_scales


def reproject_positions(positions, src_hdr, dst_hdr):
    """Marked-band arcsec -> target-band arcsec, THROUGH SKY COORDINATES.

    WHY NOT JUST COPY THE ARCSEC. The old broadcast wrote the same (y, x) into every band,
    on the reasoning that all stamps share a centre and an output WCS, so an arcsec offset
    is band-independent. That is *nearly* true and not exactly true: the stamps' tangent
    points differ by a measured 21-56 mas on the f160W pairs (make_cutouts recentres per
    band, and 0.06"/px f160W cannot land on the same sky point as 0.05"/px ACS anyway), so
    the same arcsec offset is a DIFFERENT sky position in each band -- up to a whole f160W
    pixel. Going pixel -> sky -> pixel removes that term exactly, for any pixel scale,
    stamp size or tangent point, and costs one WCS round trip per position.

    What this CANNOT fix is a band whose pixels are misregistered against the marked band's
    (the same physical object sitting at different sky coords in the two WCS). That is an
    astrometry problem upstream, not a broadcast problem, so process_lens measures it
    separately and warns rather than silently absorbing it.
    """
    ws, wd = WCS(src_hdr), WCS(dst_hdr)
    src_shape = (src_hdr['NAXIS2'], src_hdr['NAXIS1'])
    dst_shape = (dst_hdr['NAXIS2'], dst_hdr['NAXIS1'])
    src_ps = make_masks.pixel_scale_from_header(src_hdr)
    dst_ps = make_masks.pixel_scale_from_header(dst_hdr)
    out = []
    for (y, x) in positions:
        row, col = arcsec_to_pixel(y, x, src_shape, src_ps)
        sky = ws.pixel_to_world(col, row)
        dcol, drow = wd.world_to_pixel(sky)
        out.append(pixel_to_arcsec(float(drow), float(dcol), dst_shape, dst_ps))
    return out


def deflector_sky(cutout_dir, prefix):
    """Sky position of the deflector in this band, by flux centroid at the stamp centre.

    Used only to CHECK the band tie: every stamp is recentred on the deflector, so the same
    galaxy should land on the same sky coordinate in every band. Where it does not, the
    bands are misregistered and a broadcast position will miss the lensed images in that
    band by the same amount, however it is computed.
    """
    from scipy import ndimage
    hdr = fits.getheader(os.path.join(cutout_dir, f'{prefix}_sci.fits'))
    sci = np.asarray(fits.getdata(os.path.join(cutout_dir, f'{prefix}_sci.fits')), float)
    ps = make_masks.pixel_scale_from_header(hdr)
    n_y, n_x = sci.shape
    cy, cx = (n_y - 1) / 2.0, (n_x - 1) / 2.0
    h = int(round(1.2 / ps))
    sub = sci[int(cy)-h:int(cy)+h+1, int(cx)-h:int(cx)+h+1]
    iy, ix = np.unravel_index(np.argmax(ndimage.gaussian_filter(sub, 1.0)), sub.shape)
    iy += int(cy) - h; ix += int(cx) - h
    b = int(round(0.6 / ps))
    box = sci[iy-b:iy+b+1, ix-b:ix+b+1].astype(float)
    w = np.clip(box - np.median(box), 0, None)
    if w.sum() <= 0:
        return None, ps
    yy, xx = np.mgrid[iy-b:iy+b+1, ix-b:ix+b+1]
    return WCS(hdr).pixel_to_world((xx*w).sum()/w.sum(), (yy*w).sum()/w.sum()), ps


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
                       mask_center=0.0, vmax_value=None, subtract_radial=True,
                       side_by_side=True, contaminants=None):
    """Run the Clicker GUI once and return the clicked positions as a list of (y, x) arcsec.

    `sci_native` (REAL flux) is handed to al.Clicker so its brightest-pixel snap uses true
    counts; `display_base` (raw signal, or S/N -- see load_display_base) is what we STRETCH
    and imshow, so bright cores and faint arcs are both visible while clicking.

    PANELS ARE REAL AXES HERE, NOT A COMPOSITED ARRAY -- the opposite of make_masks.py, and
    deliberately. The Scribbler reads one array back, so there the panels must be hstacked
    into a single image and the scribble folded up afterwards (make_masks.fold_panels).
    al.Clicker takes only `event.xdata/ydata`, so a click needs no folding at all if every
    panel is its own axes over the SAME arcsec extent: a double-click anywhere then already
    carries panel-local sky coordinates, and is forwarded to `clicker.onclick` untouched.
    Nothing is translated, so what you click is what gets snapped, on either panel, and the
    single-panel path is the pre-existing one unchanged.

    `side_by_side` (default True, only meaningful with subtract_radial) shows both views at
    once -- subtracted LEFT, as-observed RIGHT -- because each answers a different question:
    the subtracted panel is where a faint image inside the galaxy envelope is visible at all,
    the as-observed panel is where you judge whether that blob is really there. CLICK EITHER;
    they are the same sky. Each is stretched independently, their dynamic ranges differing by
    orders of magnitude.

    Three DISPLAY-ONLY levers make the source arcs readable under the deflector light, which
    otherwise sets the colour scale on its own:
      * `contaminants` (bool array, True = masked) blanks everything this band's
        `cutout_[cr_]mask.fits` already calls a contaminant, and drops it from the stretch
        percentiles. A neighbour bright enough to be masked is also bright enough to own the
        colour scale and to be mis-clicked as an image, and it is the one region whose verdict
        is already in. Outlined, so a blank is never mistaken for empty sky.
      * `mask_center` (arcsec) blanks a disc at the stamp centre -- the deflector, since
        every stamp is recentred on it -- and drops those pixels from the percentile
        statistics, so the remaining pixels (the arcs) get the full colour range. The hidden
        region is outlined so you can see what is behind the blank.
      * `vmax_value` clips the display at an absolute value in the base array's units (S/N
        by default), saturating everything brighter.
      * `subtract_radial` removes the deflector's azimuthally-averaged radial profile,
        which is the strongest of the three -- it reaches the images buried INSIDE the
        galaxy envelope, where a central blank would hide them too. The blanked regions are
        excluded from that profile too, so a masked neighbour cannot stamp a dark ring
        across the arcs at its own radius.
    None of them touches `sci_native`, so al.Clicker still snaps on true flux and the saved
    positions are identical whatever is on screen. Do note the snap is flux-based: a click
    placed inside/adjacent to a blanked region can still snap onto a hidden pixel within
    `search_box_size` -- process_lens checks the result and warns if one did.
    """
    clicker = al.Clicker(image=sci_native, pixel_scales=pixel_scales,
                         search_box_size=search_box_size)

    base = display_base if display_base is not None else sci_native
    # Everything hidden from the display, and from the stretch and radial statistics with it.
    exclude = make_masks.central_disc(np.asarray(base).shape, pixel_scales, mask_center)
    if contaminants is not None:
        cont = np.asarray(contaminants, dtype=bool)
        exclude = cont if exclude is None else (exclude | cont)

    panels = []
    if subtract_radial:
        panels.append(('radial-subtracted',
                       make_masks.radial_median_subtract(base, exclude=exclude)))
    if not subtract_radial or side_by_side:
        panels.append(('as-observed', base))

    n_y, n_x = np.asarray(base).shape
    hw = int(n_x / 2) * pixel_scales
    ext = [-hw, hw, -hw, hw]

    fig, axes = plt.subplots(1, len(panels), figsize=(12 * len(panels), 12), squeeze=False)
    axes = list(axes[0])
    for ax, (name, values) in zip(axes, panels):
        disp = np.asarray(make_masks.stretched_display(
            values, pixel_scales, stretch=stretch, vmin_percent=vmin_percent,
            vmax_percent=vmax_percent, asinh_a=asinh_a, exclude=exclude,
            vmax_value=vmax_value).native)
        if contaminants is not None and cont.any():
            # Outline ON the masked pixels (mask_boundary is an INNER boundary), so the
            # outline survives the interior being blanked to the colourmap floor.
            disp[make_masks.mask_boundary(cont)] = 1.0
        im = ax.imshow(disp, cmap='jet', extent=ext)
        fig.colorbar(im, ax=ax, fraction=0.046)
        if mask_center and mask_center > 0:
            ax.add_patch(plt.Circle((0.0, 0.0), mask_center, fill=False, color='white',
                                    lw=1.0, ls='--'))
        ax.set_title(name if len(panels) > 1 else '')
        ax.set_xlabel('arcsec')
        ax.set_ylabel('arcsec')

    where = ' on either panel' if len(panels) > 1 else ''
    fig.suptitle(f"{lens}  {filt}  --  DOUBLE-CLICK each lensed image{where}, "
                 f"then close the window")

    axes_set = set(axes)
    moved = []

    def onclick(event):
        # Every panel is the same sky on the same arcsec extent, so a click on any of them
        # is already in the single-band frame Clicker expects. Clicks outside the panels
        # (margins, colourbars) are dropped rather than snapped to a nonsense pixel.
        if event.inaxes not in axes_set:
            return
        before = len(clicker.click_list)
        clicker.onclick(event)
        if len(clicker.click_list) > before:
            # How far the brightest-pixel snap moved the click. Printed because the snap is
            # the one step between what you aimed at and what gets saved, and a large jump
            # usually means the box reached a brighter neighbour, not a better centre.
            sy, sx = clicker.click_list[-1]
            d_as = float(np.hypot(sy - event.ydata, sx - event.xdata))
            moved.append(d_as)
            note = ('   <-- LARGE, check it'
                    if d_as > 1.5 * search_box_size * pixel_scales else '')
            print(f"    snap moved {d_as/pixel_scales:4.1f} px ({d_as:.3f}\") "
                  f"-> ({sy:+.3f}, {sx:+.3f})\"{note}")

    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
    fig.canvas.mpl_disconnect(cid)
    plt.close(fig)

    if moved:
        print(f"  snap distances: max {max(moved)/pixel_scales:.1f} px, "
              f"mean {sum(moved)/len(moved)/pixel_scales:.1f} px "
              f"(--search-box-size {search_box_size})")
    return list(clicker.click_list), moved  # each entry (y_arcsec, x_arcsec)


def rebuild_overlays(lens, filt_dirs, drizzle_pass, stretch, vmin_percent,
                     vmax_percent, asinh_a, vmax_value=None):
    """Re-draw each band's `{prefix}_positions.png` from its tracked `_positions.json`.

    The PNG is a QC view of the JSON, not a product in its own right, so it is gitignored
    while the hand-marked JSON is tracked -- and that is only safe if it can be rebuilt
    WITHOUT the GUI, which is what this does. No clicking, no `--force`, nothing written
    but the PNGs; a band with no positions JSON is skipped.
    """
    n = 0
    for filt in sorted(filt_dirs):
        cutout_dir = filt_dirs[filt]
        prefix = make_masks.find_prefix(cutout_dir, drizzle_pass)
        if prefix is None:
            continue
        json_path = os.path.join(cutout_dir, f'{prefix}_positions.json')
        if not os.path.exists(json_path):
            continue
        grid = al.from_json(file_path=json_path)
        positions = [(float(y), float(x)) for (y, x) in np.asarray(grid)]
        sci_path = os.path.join(cutout_dir, f'{prefix}_sci.fits')
        ps = make_masks.pixel_scale_from_header(fits.getheader(sci_path))
        sci_native = al.Array2D.from_fits(file_path=sci_path, pixel_scales=ps).native
        out_png = os.path.join(cutout_dir, f'{prefix}_positions.png')
        save_overlay_png(sci_native, ps, positions, out_png, stretch=stretch,
                         vmin_percent=vmin_percent, vmax_percent=vmax_percent,
                         asinh_a=asinh_a, vmax_value=vmax_value)
        print(f"  rebuilt {out_png}  ({len(positions)} positions)")
        n += 1
    if not n:
        print(f"{lens}: no positions JSON on any band -- nothing to rebuild")
    return n


def rebroadcast(lens, filt_dirs, sample, drizzle_pass, stretch, vmin_percent,
                vmax_percent, asinh_a, vmax_value=None):
    """Re-derive every band's positions from the MARKED band's, without re-clicking.

    The marked band's JSON is the hand-made measurement and is never touched; every other
    band is recomputed from it through `reproject_positions`. This exists because positions
    marked before 2026-09-17 were broadcast by copying arcsec verbatim, which carries the
    stamps' tangent-point offset (21-56 mas on the f160W pairs) into every non-marked band.
    Re-marking by hand would be the wrong fix -- the clicks were fine, the transfer was not.
    """
    # Provenance is {sample: {lens: {marked_filt: {...}}}} -- the KEY is the band that was
    # marked, and 'marked_filt' is a field repeating it inside. Normally exactly one key; if
    # a lens was re-marked on another band, take the highest-priority one, as the GUI would.
    prov = info_json.load(POSITIONS_JSON).get(sample, {}).get(lens, {})
    candidates = [k for k in prov if k in filt_dirs]
    marked_filt = pick_display_filt(sorted(candidates), None) if candidates else None
    if marked_filt is None:
        return 0                      # simply not marked yet -- nothing to say
    marked_dir = filt_dirs[marked_filt]
    marked_prefix = make_masks.find_prefix(marked_dir, drizzle_pass)
    marked_json = os.path.join(marked_dir, f'{marked_prefix}_positions.json')
    if not os.path.exists(marked_json):
        # Provenance says this lens WAS marked but the product is gone. The clicks are not
        # lost -- info/lens_positions.json keeps `positions_arcsec` -- but restoring them
        # here would be wrong: the skip-if-exists rule means a restored file silently
        # prevents the re-mark the deletion may have been for. Report, don't guess.
        n_rec = prov[marked_filt].get('n_positions')
        print(f"{lens}: provenance records {n_rec} position(s) marked on {marked_filt}, but "
              f"{os.path.relpath(marked_json, ws_path)} is MISSING.\n"
              f"    The clicks survive in {os.path.basename(POSITIONS_JSON)} "
              f"('positions_arcsec'). Re-mark with --force, or restore from there.")
        return 0
    positions = [(float(y), float(x))
                 for (y, x) in np.asarray(al.from_json(file_path=marked_json))]
    marked_hdr = fits.getheader(os.path.join(marked_dir, f'{marked_prefix}_sci.fits'))
    print(f"\n{lens}: re-deriving {len(positions)} position(s) from {marked_filt}")
    n = 0
    for filt in sorted(filt_dirs):
        if filt == marked_filt:
            continue
        cutout_dir = filt_dirs[filt]
        prefix = make_masks.find_prefix(cutout_dir, drizzle_pass)
        if prefix is None:
            continue
        band_sci = os.path.join(cutout_dir, f'{prefix}_sci.fits')
        band_hdr = fits.getheader(band_sci)
        band_ps = make_masks.pixel_scale_from_header(band_hdr)
        new = reproject_positions(positions, marked_hdr, band_hdr)
        json_path = os.path.join(cutout_dir, f'{prefix}_positions.json')
        old = None
        if os.path.exists(json_path):
            old = [(float(y), float(x))
                   for (y, x) in np.asarray(al.from_json(file_path=json_path))]
        al.output_to_json(obj=al.Grid2DIrregular(values=new), file_path=json_path)
        band_native = al.Array2D.from_fits(file_path=band_sci, pixel_scales=band_ps).native
        save_overlay_png(band_native, band_ps, new,
                         os.path.join(cutout_dir, f'{prefix}_positions.png'),
                         stretch=stretch, vmin_percent=vmin_percent,
                         vmax_percent=vmax_percent, asinh_a=asinh_a, vmax_value=vmax_value)
        moved = ('' if old is None or len(old) != len(new) else
                 f"  moved {1000*max(np.hypot(a-c, b-d) for (a, b), (c, d) in zip(new, old)):.1f} mas"
                 f" ({max(np.hypot(a-c, b-d) for (a, b), (c, d) in zip(new, old))/band_ps:.2f} px)")
        print(f"  {filt}: rewritten{moved}")
        n += 1
    return n


def process_lens(lens, filt_dirs, sample, requested_filt, drizzle_pass, force,
                 search_box_size, display, stretch, vmin_percent, vmax_percent, asinh_a,
                 mask_center=0.0, vmax_value=None, subtract_radial=True,
                 side_by_side=True, hide_masked=True):
    """Mark positions once for one lens and broadcast the result to every band.

    `filt_dirs` maps filt -> cutout_dir for that lens, from make_masks.discover_targets.
    Returns True if positions were written.
    """
    filts = sorted(filt_dirs)
    display_filt = pick_display_filt(filts, requested_filt)
    if display_filt is None:
        print(f"{lens}: requested --filt {requested_filt} not present "
              f"(have {', '.join(filts)}), skipping")
        return False

    display_dir = filt_dirs[display_filt]                # cutout dir of the marked band
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
    # This band's already-judged contaminants, loaded exactly as make_arc_masks.py does.
    contaminants = None
    if hide_masked:
        cont_path = os.path.join(display_dir, f'{display_prefix}_mask.fits')
        if os.path.exists(cont_path):
            cont = np.asarray(fits.getdata(cont_path), dtype=bool)
            if cont.any():
                contaminants = cont
                print(f"  {int(cont.sum())}px contaminant-masked on this band -- BLANKED in "
                      f"the display (outlined) and dropped from the stretch and radial "
                      f"profile; --no-hide-masked shows them")
        else:
            print(f"  NOTE: no {display_prefix}_mask.fits on this band -- nothing to hide "
                  f"(draw one with make_masks.py)")

    source_label += (', radial-subtracted + as-observed' if subtract_radial and side_by_side
                     else ', radial-subtracted' if subtract_radial else '')
    if contaminants is not None:
        source_label += ' (contaminants blanked)'

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
        if side_by_side:
            print("  side by side: radial-subtracted LEFT, as-observed RIGHT -- click "
                  "either, they are the same sky")
    print("  DOUBLE-CLICK each lensed image of the source (2 for a double, 4 for a quad).")
    print("  Each click snaps to the brightest nearby pixel. Close the window when done.")

    positions, moved = mark_positions_gui(
        sci_native, pixel_scales, lens, display_filt, search_box_size,
        display=display, display_base=base, stretch=stretch,
        vmin_percent=vmin_percent, vmax_percent=vmax_percent, asinh_a=asinh_a,
        mask_center=mask_center, vmax_value=vmax_value,
        subtract_radial=subtract_radial, side_by_side=side_by_side,
        contaminants=contaminants)

    if len(positions) == 0:
        print(f"  no positions clicked -- nothing written for {lens}")
        return False
    if len(positions) < 2:
        print(f"  WARNING: only {len(positions)} position clicked -- a lensed source has "
              f">=2 images. Re-run with --force if this was a mis-click.")

    # The snap uses true flux, so a click at the edge of a blanked region can still land
    # inside it -- the one way a display-only lever can reach the saved product. Check rather
    # than trust, and name the offenders: a position on a contaminant is a bad position.
    n_on_mask = 0
    if contaminants is not None:
        for (y, x) in positions:
            iy = int(round((np.asarray(contaminants).shape[0] - 1) / 2.0 - y / pixel_scales))
            ix = int(round((np.asarray(contaminants).shape[1] - 1) / 2.0 + x / pixel_scales))
            if (0 <= iy < contaminants.shape[0] and 0 <= ix < contaminants.shape[1]
                    and contaminants[iy, ix]):
                n_on_mask += 1
                print(f"  WARNING: position ({y:+.3f}, {x:+.3f})\" snapped ONTO a "
                      f"contaminant-masked pixel -- the snap uses true flux, so it can reach "
                      f"into a blanked region. Re-run with --force to re-mark it.")

    marked_hdr = fits.getheader(sci_path)
    marked_sky, _ = deflector_sky(display_dir, display_prefix)

    # Broadcast THROUGH SKY, not by copying arcsec (see reproject_positions). Each band gets
    # its own Grid2DIrregular in its own frame, one file per band into that band's cutout
    # dir, under that dir's own pass prefix (f160W has no CR pass -> cutout_positions.json).
    written = []
    for filt in filts:
        cutout_dir = filt_dirs[filt]
        prefix = make_masks.find_prefix(cutout_dir, drizzle_pass)
        if prefix is None:
            continue
        band_sci = os.path.join(cutout_dir, f'{prefix}_sci.fits')
        band_hdr = fits.getheader(band_sci)
        band_ps = make_masks.pixel_scale_from_header(band_hdr)
        if filt == display_filt:
            band_positions, shift_mas = positions, 0.0
        else:
            band_positions = reproject_positions(positions, marked_hdr, band_hdr)
            shift_mas = 1000.0 * max(np.hypot(by - y, bx - x)
                                     for (by, bx), (y, x) in zip(band_positions, positions))

        # Does this band's own image actually sit where the marked band says it does? A
        # broadcast cannot fix a misregistered band; it can only be honest about it.
        tie_mas = None
        if marked_sky is not None and filt != display_filt:
            band_sky, _ = deflector_sky(cutout_dir, prefix)
            if band_sky is not None:
                tie_mas = float(band_sky.separation(marked_sky).arcsec * 1000.0)
                if tie_mas / (band_ps * 1000.0) > 1.0:
                    print(f"  WARNING: {filt} is MISREGISTERED against {display_filt} by "
                          f"{tie_mas:.0f} mas = {tie_mas/(band_ps*1000.0):.2f} {filt} px "
                          f"(deflector centroid). The broadcast positions are correct on "
                          f"SKY but will sit that far off this band's images. Fix the "
                          f"astrometry, don't re-mark.")

        json_path = os.path.join(cutout_dir, f'{prefix}_positions.json')
        al.output_to_json(obj=al.Grid2DIrregular(values=band_positions), file_path=json_path)
        band_native = al.Array2D.from_fits(file_path=band_sci, pixel_scales=band_ps).native
        save_overlay_png(band_native, band_ps, band_positions,
                         os.path.join(cutout_dir, f'{prefix}_positions.png'),
                         stretch=stretch, vmin_percent=vmin_percent,
                         vmax_percent=vmax_percent, asinh_a=asinh_a,
                         vmax_value=vmax_value)
        written.append({'filt': filt, 'prefix': prefix,
                        'reprojection_shift_mas': round(shift_mas, 2),
                        'deflector_tie_mas': None if tie_mas is None else round(tie_mas, 1),
                        'deflector_tie_px': (None if tie_mas is None
                                             else round(tie_mas / (band_ps*1000.0), 3))})
        extra = '' if filt == display_filt else f'  (+{shift_mas:.0f} mas vs a raw copy)'
        print(f"  wrote {json_path}{extra}")

    info_json.update(POSITIONS_JSON, sample, lens, display_filt, {
        'n_positions': len(positions),
        'positions_arcsec': [[round(y, 6), round(x, 6)] for (y, x) in positions],
        'marked_filt': display_filt,
        'marked_prefix': display_prefix,
        'drizzle_pass': 'cr' if display_prefix == 'cutout_cr' else 'nocrrej',
        'pixel_scale_arcsec': round(pixel_scales, 6),
        'search_box_size': search_box_size,
        'snap_moved_px_max': (round(max(moved)/pixel_scales, 2) if moved else None),
        'snap_moved_px_mean': (round(sum(moved)/len(moved)/pixel_scales, 2) if moved else None),
        'broadcast': 'sky (WCS pixel->world->pixel per band)',
        'display': source_label,
        'display_mask_center_arcsec': mask_center or None,
        'display_vmax_value': vmax_value,
        'display_subtract_radial': bool(subtract_radial),
        'display_side_by_side': bool(subtract_radial and side_by_side),
        'display_hidden_contaminant_px': (int(contaminants.sum())
                                          if contaminants is not None else None),
        'n_positions_on_contaminant': n_on_mask,
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
                        'data/cutouts/ -- the git-tracked tree; the hand-marked '
                        '_positions.json is non-regenerable, like masks, so it belongs here. '
                        'Its _positions.png is a QC view of that JSON and IS regenerable '
                        '(--overlays-only), so it is gitignored)')
    p.add_argument('--force', action='store_true', default=False,
                   help='re-mark a lens that already has positions (default: skip it)')
    p.add_argument('--search-box-size', type=int, default=2,
                   help='half-width in pixels of the brightest-pixel snap around each click '
                        '(al.Clicker search_box_size; default 2, lowered from 5 on '
                        '2026-09-17). 5 let a click move up to ~0.25" on ACS -- far enough '
                        'to jump onto the deflector envelope or a neighbouring knot instead '
                        'of the image you aimed at. 2 still catches the peak of an image you '
                        'clicked squarely; use 1 to snap only to the immediate neighbours, '
                        'or 0 to take the click as-is. Each click prints how far it actually '
                        'moved, and the per-lens max is recorded in the provenance')
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
    p.add_argument('--rebroadcast', action='store_true', default=False,
                   help="re-derive every non-marked band's positions from the MARKED band's, "
                        'through the WCS, and exit -- no GUI. Use it on positions marked '
                        'before 2026-09-17, when the broadcast copied arcsec verbatim and so '
                        "carried each band's tangent-point offset (21-56 mas on the f160W "
                        'pairs) into the saved positions. The marked band is never touched')
    p.add_argument('--overlays-only', action='store_true', default=False,
                   help='rebuild each band\'s {prefix}_positions.png from the tracked '
                        '{prefix}_positions.json and exit -- no GUI, no clicking, nothing '
                        'else written. The PNGs are gitignored (the JSON is the product); '
                        'this is how you get them back in a fresh clone, or after changing '
                        'a display flag like --vmax-value')
    p.add_argument('--side-by-side', action=argparse.BooleanOptionalAction, default=True,
                   help='show the radial-subtracted and as-observed views SIDE BY SIDE '
                        '(default on; only meaningful with --subtract-radial). Each answers '
                        'a different question -- the subtracted panel is where a faint image '
                        'inside the galaxy envelope is visible at all, the as-observed panel '
                        'is where you judge whether it is really there. Double-click either; '
                        'both are the same sky on the same arcsec extent, so a click on '
                        'either snaps identically. --no-side-by-side gives the single '
                        'subtracted panel')
    p.add_argument('--hide-masked', action=argparse.BooleanOptionalAction, default=True,
                   help="BLANK this band's already-masked contaminants "
                        '(cutout_[cr_]mask.fits) in the DISPLAY, outlining what was hidden, '
                        'and drop them from both the stretch percentiles and the radial '
                        'profile (default on). A neighbour bright enough to be masked is also '
                        'bright enough to own the colour scale, to be mis-clicked as a lensed '
                        'image, and to stamp a dark ring across the arcs at its own radius. '
                        'Display only -- the snap still uses true flux, and a position that '
                        'lands on a masked pixel anyway is reported. --no-hide-masked shows '
                        'the band untouched; silently inactive where the band has no mask')
    p.add_argument('--subtract-radial', action=argparse.BooleanOptionalAction, default=True,
                   help="subtract the deflector's azimuthally-averaged radial profile from "
                        'the DISPLAYED image (DEFAULT ON -- it is what makes the images '
                        'clickable at all on most lenses). The strongest of the three '
                        'arc-finding levers: unlike --mask-center it reveals images buried '
                        'inside the galaxy envelope rather than hiding them. Expect a '
                        'quadrupole residual (the deflector is elliptical, not circular) '
                        'and treat it as a finding aid, not photometry. Display only -- pass '
                        '--no-subtract-radial to click on the raw image instead')
    p.add_argument('--include-ignored', action='store_true', default=False,
                   help='also offer bands whose cutout directory is GITIGNORED. Off by '
                        'default: such a band is not a science product (the three 420 s '
                        'SLACS SNAP f814W diagnostics are why the rule exists), and hand-marked '
                        'positions written there are non-regenerable and invisible to every '
                        'clone. This guard exists because it already happened -- J1538+5817 '
                        'was MARKED on its SNAP f814W band and broadcast from there across a '
                        '688 mas (13.8 px) misregistration; both copies were deleted 2026-09-22.')
    a = p.parse_args()

    root = cutout_paths.cutouts_root(ws_path, a.size)

    # discover_targets yields per (lens, filt); regroup to filt -> cutout_dir per lens so we
    # mark once and broadcast across the lens's bands.
    lens_filts = {}
    for lens, filt, cutout_dir in make_masks.discover_targets(
            root, a.sample, a.lens, None,     # filt=None: consider all bands, pick best below
            include_ignored=a.include_ignored):
        lens_filts.setdefault(lens, {})[filt] = cutout_dir

    if not lens_filts:
        raise SystemExit(f"no cutouts found under {root} for sample {a.sample} matching "
                         f"lens={a.lens!r}")

    print(f"{len(lens_filts)} lens(es) to process")

    if a.rebroadcast:
        total = 0
        for lens in sorted(lens_filts):
            total += rebroadcast(lens, lens_filts[lens], a.sample, a.drizzle_pass, a.stretch,
                                 a.vmin_percent, a.vmax_percent, a.asinh_a, a.vmax_value)
        print(f"\nDone: {total} band(s) re-derived through the WCS")
        return

    if a.overlays_only:
        total = 0
        for lens in sorted(lens_filts):
            print(f"\n{lens}")
            total += rebuild_overlays(lens, lens_filts[lens], a.drizzle_pass, a.stretch,
                                      a.vmin_percent, a.vmax_percent, a.asinh_a,
                                      a.vmax_value)
        print(f"\nDone: {total} QC overlay(s) rebuilt from the positions JSONs")
        return

    made = skipped = 0
    for lens in sorted(lens_filts):
        if process_lens(lens, lens_filts[lens], a.sample, a.filt, a.drizzle_pass, a.force,
                        a.search_box_size, a.display, a.stretch,
                        a.vmin_percent, a.vmax_percent, a.asinh_a,
                        a.mask_center, a.vmax_value, a.subtract_radial,
                        a.side_by_side, a.hide_masked):
            made += 1
        else:
            skipped += 1

    print(f"\nDone: {made} lens(es) marked, {skipped} skipped")


if __name__ == '__main__':
    main()
