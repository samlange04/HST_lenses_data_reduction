#!/usr/bin/env python
"""
Interactive GUI for marking multiply-imaged (lensed) source positions, per lens.

Merged 2026-11-28 across the HST and JWST repos: the file is identical in both (the repo
profile is the import seam below; the stamp layout is each repo's own `stamp_paths`). North
is DOWN and East LEFT in every panel (display_look.ORIGIN; the JWST side drew North-up until
the merge, and its positions GUI dropped the y-flip here).

For each lens in a sample this opens a window over the lens's
best-available stamp -- PyAutoLens's `al.Clicker`'s panel layout, with the repo's own snap --
so you can
DOUBLE-CLICK each lensed image of the source (2 for a double, 4 for a quad, ...). Each
click snaps to the brightest pixel within a `--snap-radius` DISC (snap_click) and is recorded in
arcsec relative to the stamp centre; the collected clicks are saved as an
`al.Grid2DIrregular` in `stamp12_positions.json` (via `al.output_to_json`), the file a
modelling script loads to build a positions likelihood penalty (`al.PositionsLH`) that
rejects mass models mapping the images too far apart in the source plane.

WHAT IS ON SCREEN. With --all-bands (default since 2026-11-17), EVERY band of the lens is
shown -- one column per band -- each with the band as observed on the top row and, below it,
the band with the deflector removed by an ISOPHOTE model (since 2026-11-01, smoothed on a
red-blue map: red is light the model lacks; `--subtract radial` restores the old radial-median
panel) -- because each answers a different question: the
subtracted panels are where an image buried in the galaxy envelope shows up at all, the
as-observed panels are where you judge whether that blob is really there or is a subtraction
artefact. They are separate axes over each band's own arcsec extent, so DOUBLE-CLICK ANY
PANEL: the snap runs on THAT band's true flux and the position is stored in the marked
band's frame through the two bands' WCS -- an image that is clearer in another band (say
f555W) no longer has to be marked through the marked band's glare. Anything a band's
`stamp12_mask.fits` already calls a contaminant is blanked (and outlined) in that band's
panels, and dropped from the stretch and the isophote fit: that verdict is already in, and a
masked neighbour otherwise owns the colour scale and invites a mis-click. --zoom (default)
crops the panels to the central 7" (--zoom-width); off automatically on multi_galaxy lenses.
All of it is display only -- the snap runs on true flux and the saved positions do not
depend on any of it.

WHY ONE BAND PER LENS, THEN BROADCAST. Image positions are band-INDEPENDENT sky
coordinates: every band is cut about the same shared centre (`make_stamps.py
--center-band`, the sample's first available centre band) and pinned to the same output WCS geometry (`final_rot=0`, tangent
point at the lens position), so a position measured in ARCSEC relative to the stamp centre
is identical in every band -- even where the pixel scales differ (f160W 0.06" vs 0.05" for
the optical bands), because arcsec, not pixels, is the shared frame PyAutoLens works in. So
you mark ONCE, on the highest-priority band available (f475X > f814W > f606W > f555W >
f160W > ...), and
the same Grid2DIrregular is broadcast into every band's stamp dir under that band's own
pass prefix, so any per-band modelling script finds a local positions.json next to its
sci/noise/psf/mask. `--filt` forces which band you mark on.

Trees mirror make_masks.py exactly: `--size` selects the stamp tree (default 12",
data/stamps/, the git-tracked one -- and, like hand-drawn masks, these hand-marked
positions are a NON-regenerable product, so they belong in the tracked tree). There is one
stamp dir per band (stamp_paths.py), so each band's file lands beside the sci it was
marked on, with nothing to choose and no second copy to keep in step.

This is a manual, one-lens-at-a-time tool (not a batch driver): each lens blocks on its own
Tk window until you close it. A lens that already has positions is skipped so a run resumes
across a sample; --force re-marks.

Usage:
    uv run python scripts/make_positions.py --sample slacs            # best band per lens
    uv run python scripts/make_positions.py --lens J0008-0004              # one lens, best band
    uv run python scripts/make_positions.py --lens J0008-0004 --filt f606W # force the band
    uv run python scripts/make_positions.py --sample slacs --force    # re-mark everything
    uv run python scripts/make_positions.py --lens J0008-0004 --no-subtract-radial
        # opt OUT of the default radial-profile subtraction and click on the raw image
    uv run python scripts/make_positions.py --lens J0008-0004 --mask-center 0.8 --vmax-value 10
        # hide the central 0.8" of deflector light and saturate above S/N 10, so the arcs
        # are the brightest thing on screen
    uv run python scripts/make_positions.py --lens J0008-0004 --no-side-by-side
        # one view (the subtracted one) per band instead of the default pair
    uv run python scripts/make_positions.py --sample slacs --overlays-only
        # rebuild the gitignored QC PNGs from the tracked positions JSONs -- no GUI
    uv run python scripts/make_positions.py --lens J0008-0004 --no-hide-masked
        # show the contaminants this band's mask already excludes, instead of blanking them
    uv run python scripts/make_positions.py --lens J0008-0004 --no-all-bands --no-zoom
        # the pre-2026-11-17 layout: the marked band alone, full stamp
    (all are DISPLAY ONLY -- the click snap uses true flux, saved positions are identical)
"""

import argparse
import os
import sys
import warnings

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib import pyplot as plt
from matplotlib import patheffects
from matplotlib.ticker import MaxNLocator

# autonerves prints a workspace-version-mismatch UserWarning on import in this repo -- see
# make_masks.py; silence it so it doesn't repeat per lens.
warnings.filterwarnings('ignore', category=UserWarning, module='autonerves')

import autolens as al

ws_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    import lens_samples as mast_target_names   # the JWST repo's profile module
except ImportError:
    import mast_target_names                   # the HST repo's
import info_json
import stamp_paths
import positions_io
import isophote_model
# Reuse make_masks' display/tree/discovery helpers so the two GUIs stay in lock-step
# (same pass-prefix rule, same tree discovery, same stretched display).
import make_masks

POSITIONS_JSON = os.path.join(ws_path, 'info', 'provenance', 'lens_positions.json')

# "Best band to mark on" is shared with make_masks so the two tools pick the same primary
# band per lens (f475X>f814W>f606W>f555W>f160W>...).
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
    points differ by a measured 21-56 mas on the f160W pairs (make_stamps recentres per
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


SNAP_RADIUS_DEFAULT = 1   # px: the click snap's disc radius (brightest pixel, true flux)

TIE_EXCLUDE_ARCSEC = 2.5   # field sources only: lens, ring and arcs live inside this radius
TIE_APERTURE_ARCSEC = 0.24 # centroid window, the same SKY size in every band
TIE_MIN_SOURCES = 3


def snap_click(sci_native, pixel_scales, y, x, radius):
    """The click at (y, x) arcsec moved to the centre of the brightest pixel (true flux) within
    `radius` px of the clicked pixel, a DISC, so no direction is favoured; radius 0 keeps the
    clicked pixel. Returns (y, x) arcsec.

    Replaces al.Clicker's snap (the JWST repo, 2026-11-02; adopted by the merge 2026-11-28),
    which rounded the click to a pixel CORNER on these even-sided stamps before picking a
    pixel, and searched an asymmetric box (-b .. b-1 px), so a click could move further one
    way than the other."""
    img = np.asarray(sci_native, dtype=float)
    row, col = arcsec_to_pixel(y, x, img.shape, pixel_scales)
    r0, c0 = int(round(row)), int(round(col))
    best, br, bc = -np.inf, r0, c0
    rad = int(max(radius, 0))
    for dr in range(-rad, rad + 1):
        for dc in range(-rad, rad + 1):
            r_, c_ = r0 + dr, c0 + dc
            if dr * dr + dc * dc > radius * radius or not (0 <= r_ < img.shape[0] and 0 <= c_ < img.shape[1]):
                continue
            v = img[r_, c_]
            if np.isfinite(v) and v > best:
                best, br, bc = v, r_, c_
    return pixel_to_arcsec(br, bc, img.shape, pixel_scales)


def _windowed_centroid(img, x, y, r_px, n_iter=5):
    """Iterated flux centroid in a circle of r_px about (x, y); None if it walks off."""
    ny, nx = img.shape
    b = int(np.ceil(r_px))
    for _ in range(n_iter):
        ix, iy = int(round(x)), int(round(y))
        if ix - b < 0 or iy - b < 0 or ix + b >= nx or iy + b >= ny:
            return None
        yy, xx = np.mgrid[iy-b:iy+b+1, ix-b:ix+b+1]
        box = img[iy-b:iy+b+1, ix-b:ix+b+1]
        w = np.clip(box - np.median(box), 0, None) * (np.hypot(xx - x, yy - y) <= r_px)
        if not np.isfinite(w).all() or w.sum() <= 0:
            return None
        x, y = (xx*w).sum()/w.sum(), (yy*w).sum()/w.sum()
    return x, y


def band_tie(marked_dir, marked_prefix, band_dir, band_prefix):
    """How far this band's pixels sit from the marked band's on SKY, from FIELD sources.

    Compact sources outside TIE_EXCLUDE_ARCSEC are detected on the coarser of the two
    stamps, centroided in the same sky aperture in both, and the median (dRA*, dDec) offset
    is returned as {'mas', 'scatter_mas', 'n', 'method': 'field sources'}.

    A 12" stamp of a shallow field can hold no compact source at all (the SPT stamps: 0-1).
    Then, and only then, it falls back to the lens core: both bands centroided in the same
    TIE_APERTURE_ARCSEC window, iterated from the SAME sky point (the marked stamp's centre,
    which make_stamps puts on the deflector). Starting both from one sky point, not from each
    band's own brightest pixel, is what keeps a ring that outshines its lens in one band from
    being taken for the lens; it still assumes the core is the local light peak in both
    bands, so the result is labelled 'lens core (fallback)'. 'mas' is None if that fails too.

    WHY NOT THE DEFLECTOR. Centroiding "the brightest thing near the stamp centre" in each
    band assumes that thing is the same object in both. On B1938+666 (2026-10-04) it was not:
    the stamps were 0.52" off the lens and in F277W the Einstein ring outshines the galaxy,
    so the check found the lens core in F115W and the top of the ring in F277W and reported
    a 281 mas misregistration between bands that field sources tie to 10 mas. (The HST repo
    hit the same failure on its SL2S F475X and GALLERY F438W bands, where the arcs outshine
    the deflector.) Any lens whose arcs rival its light in one band fails the same way.
    """
    from scipy import ndimage

    def load(d, prefix):
        path = stamp_paths.find_stamp(d, prefix, 'sci')
        hdr = fits.getheader(path)
        img = np.nan_to_num(np.asarray(fits.getdata(path), float))
        return img, WCS(hdr), make_masks.pixel_scale_from_header(hdr)

    stamps = [load(marked_dir, marked_prefix), load(band_dir, band_prefix)]
    det, det_wcs, det_ps = max(stamps, key=lambda t: t[2])
    ny, nx = det.shape
    med = np.median(det)
    sig = 1.4826 * np.median(np.abs(det - med))
    lab, _ = ndimage.label(ndimage.gaussian_filter(det, 1.5) > med + 8.0 * sig)
    r_det = TIE_APERTURE_ARCSEC / det_ps
    offsets = []
    for sl in ndimage.find_objects(lab):
        cy = (sl[0].start + sl[0].stop - 1) / 2.0
        cx = (sl[1].start + sl[1].stop - 1) / 2.0
        size = (sl[0].stop - sl[0].start) * (sl[1].stop - sl[1].start)
        if (np.hypot(cx - (nx - 1) / 2.0, cy - (ny - 1) / 2.0) * det_ps < TIE_EXCLUDE_ARCSEC
                or size * det_ps**2 > 0.5):            # extended: no well-defined centroid
            continue
        sky = []
        for img, w, ps in stamps:
            if not sky:
                guess = (cx, cy) if img is det else tuple(float(v) for v in w.world_to_pixel(
                    det_wcs.pixel_to_world(cx, cy)))
            else:
                guess = tuple(float(v) for v in w.world_to_pixel(sky[0]))
            c = _windowed_centroid(img, guess[0], guess[1], TIE_APERTURE_ARCSEC / ps)
            if c is None:
                break
            sky.append(w.pixel_to_world(*c))
        if len(sky) == 2:
            dra, ddec = sky[0].spherical_offsets_to(sky[1])
            offsets.append((dra.arcsec * 1e3, ddec.arcsec * 1e3))
    if len(offsets) >= TIE_MIN_SOURCES:
        off = np.asarray(offsets)
        med_off = np.median(off, axis=0)
        scatter = 1.4826 * np.median(np.hypot(*(off - med_off).T))
        return {'mas': float(np.hypot(*med_off)), 'scatter_mas': float(scatter),
                'n': len(offsets), 'method': 'field sources'}

    img0, w0, ps0 = stamps[0]
    start = w0.pixel_to_world((img0.shape[1] - 1) / 2.0, (img0.shape[0] - 1) / 2.0)
    core = []
    for img, w, ps in stamps:
        gx, gy = (float(v) for v in w.world_to_pixel(start))
        c = _windowed_centroid(img, gx, gy, TIE_APERTURE_ARCSEC / ps)
        if c is None:
            return {'mas': None, 'scatter_mas': None, 'n': len(offsets), 'method': None}
        core.append(w.pixel_to_world(*c))
    return {'mas': float(core[0].separation(core[1]).arcsec * 1e3), 'scatter_mas': None,
            'n': len(offsets), 'method': 'lens core (fallback)'}


# One colour per position set (source 1, 2, ...), each drawn with a black edge so it reads
# on the jet display as well as on a light or dark patch of it.
STROKE = [patheffects.withStroke(linewidth=4, foreground='black')]
MARK_STROKE = [patheffects.withStroke(linewidth=2.5, foreground='black')]  # the GUI's click crosses
# The as-observed panels (GUI and PNG) use the dataset QC subplots' look since 2026-11-02: inferno
# with 0.5-99.5 percentile limits (PercentileInterval(99)) and AsinhStretch(0.1).
DISPLAY_CMAP = 'inferno'
SET_COLOURS = ['white', 'magenta', 'lime', 'orange', 'cyan', 'yellow', 'red', 'deepskyblue', 'pink']


def save_overlay_png(sci_native, pixel_scales, sets, out_path,
                     stretch='asinh', vmin_percent=0.5, vmax_percent=99.5, asinh_a=0.1,
                     vmax_value=None, subtract='isophote'):
    """QC overlay: the marked positions over the band's AS-OBSERVED S/N (left) and its
    ISOPHOTE-SUBTRACTED view (right), the same two panels the GUI shows (standard since
    2026-11-01; `subtract=None` gives the as-observed panel alone). Written per band so you
    can confirm the shared positions land on the lensed images in EVERY band.

    The S/N base and the contaminant mask are read from the stamp directory `out_path` sits
    in (`{prefix}_noise.fits`, `{prefix}_mask.fits`), so every caller gets the same picture
    without passing them. Contaminants are blanked in the subtracted panel and outlined in
    both, and kept out of the isophote fit and the stretch. `vmax_value` (--vmax-value)
    saturates the bright core of the as-observed panel; the central blank (--mask-center) is
    deliberately NOT applied, so the deflector stays visible.
    """
    stamp_dir = os.path.dirname(out_path)
    prefix = os.path.basename(out_path)[:-len('_positions.png')]
    base = sci_native
    if stamp_paths.find_stamp(stamp_dir, prefix, 'noise') is not None:
        base, _ = make_masks.load_display_base(stamp_dir, prefix, 'snr', pixel_scales)
    mask_path = os.path.join(stamp_dir, f'{prefix}_mask.fits')
    cont = (np.asarray(fits.getdata(mask_path), dtype=bool) if os.path.exists(mask_path)
            else np.zeros(np.asarray(base).shape, bool))

    n_y, n_x = np.asarray(base).shape
    hw = int(n_x / 2) * pixel_scales
    ext = [-hw, hw, -hw, hw]
    panels = [('as-observed (S/N)', None)]
    if subtract == 'isophote':
        got = isophote_model.deflector_pixel(stamp_paths.find_stamp(stamp_dir, prefix, 'sci'))
        res, vmax_div = isophote_model.isophote_display(
            base, pixel_scales, exclude=cont if cont.any() else None,
            centre=got[:2] if got else None, refine_centre=got[2] if got else None)
        panels.append(('isophote-subtracted (smoothed)', (res, vmax_div)))
    fig, axes = plt.subplots(1, len(panels), figsize=(7 * len(panels), 7), squeeze=False)
    for ax, (name, sub) in zip(axes[0], panels):
        if sub is None:
            disp = make_masks.stretched_display(base, pixel_scales, stretch=stretch,
                                                vmin_percent=vmin_percent,
                                                vmax_percent=vmax_percent, asinh_a=asinh_a,
                                                exclude=cont if cont.any() else None,
                                                vmax_value=vmax_value)
            ax.imshow(np.asarray(disp.native), cmap=DISPLAY_CMAP, extent=ext)
        else:
            res, vmax_div = sub
            disp = np.array(res, dtype=float)
            disp[cont] = np.nan
            ax.imshow(disp, cmap='RdBu_r', vmin=-vmax_div, vmax=vmax_div, extent=ext)
        if cont.any():   # white on a black stroke over inferno, black over the red-blue residual
            cs = ax.contour(cont, levels=[0.5], colors='white' if sub is None else 'black',
                            linewidths=0.8, extent=ext, origin='upper')
            if sub is None:
                cs.set(path_effects=[patheffects.withStroke(linewidth=2.2, foreground='black')])
        for k, positions in enumerate(sets):
            if not positions:
                continue
            colour = SET_COLOURS[k % len(SET_COLOURS)]
            ax.scatter([p[1] for p in positions], [p[0] for p in positions], marker='+',
                       c=colour, s=80, linewidths=1.2, path_effects=MARK_STROKE, zorder=5)
            for i, (y, x) in enumerate(positions):
                label = str(i) if len(sets) == 1 else f's{k + 1}.{i}'
                ax.annotate(label, (x, y), color=colour, fontsize=9, xytext=(4, 4),
                            textcoords='offset points', path_effects=MARK_STROKE, zorder=5)
        ax.set_title(name if len(panels) > 1 else '', fontsize=10)
        ax.set_xlabel('arcsec'); ax.set_ylabel('arcsec')
    fig.tight_layout()
    fig.savefig(out_path, dpi=110)
    plt.close(fig)


def load_band_panel_info(lens, filt, stamp_dir, drizzle_pass, display, hide_masked):
    """One band's GUI data: the sci array the snap runs on, the display base, header/WCS,
    pixel scale, contaminant mask and deflector pixel. None if the band has no stamp for
    the pass. All panels of the all-bands GUI (and the single-band path, which is the
    all-bands machinery on one band) are built from this."""
    prefix = make_masks.find_prefix(stamp_dir, drizzle_pass)
    if prefix is None:
        print(f"{lens} {filt}: no stamp sci for --pass {drizzle_pass}, band not shown")
        return None
    sci_path = stamp_paths.find_stamp(stamp_dir, prefix, 'sci')
    hdr = fits.getheader(sci_path)
    ps = make_masks.pixel_scale_from_header(hdr)
    sci_native = al.Array2D.from_fits(file_path=sci_path, pixel_scales=ps).native
    noise_path = stamp_paths.find_stamp(stamp_dir, prefix, 'noise')
    noise = (np.asarray(al.Array2D.from_fits(file_path=noise_path,
                                             pixel_scales=ps).native, dtype=float)
             if noise_path else None)
    psf_path = os.path.join(stamp_dir, 'psf.fits')
    psf = (np.nan_to_num(np.asarray(fits.getdata(psf_path), dtype=float))
           if os.path.exists(psf_path) else None)
    base, source_label = make_masks.load_display_base(stamp_dir, prefix, display, ps)
    contaminants = None
    if hide_masked:
        cont_path = os.path.join(stamp_dir, f'{prefix}_mask.fits')
        if os.path.exists(cont_path):
            cont = np.asarray(fits.getdata(cont_path), dtype=bool)
            if cont.any():
                contaminants = cont
    return {'filt': filt, 'prefix': prefix, 'stamp_dir': stamp_dir, 'sci_path': sci_path,
            'hdr': hdr, 'ps': ps, 'sci': sci_native, 'noise': noise, 'psf': psf,
            'base': base,
            'source_label': source_label, 'contaminants': contaminants,
            # (col, row, refine) or None: the profile's deflector position, where known
            'deflector': isophote_model.deflector_pixel(sci_path)}


def mark_positions_gui(bands, marked_filt, lens, search_box_size,
                       display='snr', stretch='asinh',
                       vmin_percent=0.5, vmax_percent=99.5, asinh_a=0.1,
                       mask_center=0.0, vmax_value=None, subtract='isophote',
                       side_by_side=True, zoom_width=None, colour=True, prompt=None):
    """Run the Clicker GUI once and return the clicked position SETS, [[(y, x), ...], ...]
    arcsec IN THE MARKED BAND'S FRAME, one list per source, plus the per-position record of
    which band each click was made on.

    `bands` is a list of load_band_panel_info dicts -- EVERY band of the lens by default
    (--all-bands, 2026-11-17), so the band where an image is clearest can be the one it is
    clicked in (an arc invisible in f814W can be obvious in f555W). One row of panels per
    view (as-observed, and the subtracted view with --subtract), one COLUMN per band. A
    double-click on any panel snaps on THAT band's true flux (snap_click on that band's sci
    array), is converted into the marked band's frame through the two bands' WCS (pixel ->
    world -> pixel, the same hop the broadcast uses), stored there, and drawn back onto
    every panel.

    NUMBER KEYS SELECT THE SOURCE BEING MARKED (2026-11-01): the window opens on source 1;
    press 2 and every following double-click goes to source 2, press 1 to go back. Use a set
    per source plane (or per distinct source), because each set's images must map to their
    OWN source-plane point; the title names the active set and each set is drawn in its own
    colour. Empty sets are dropped, so pressing 3 and marking nothing changes nothing.
    **'z' (or Backspace) UNDOES the last click**, whichever source and band it was on
    (markers removed from every panel; 2026-11-17).

    PANELS ARE REAL AXES HERE, NOT A COMPOSITED ARRAY -- the opposite of the Scribbler, and
    deliberately. The Scribbler reads one array back, so there the panels are hstacked into a
    single image and the scribble folded up afterwards (al.Scribbler does that itself now:
    `autogalaxy.gui.display_util.composite_panels` / `fold_panels`).
    A click carries only `event.xdata/ydata`, so it needs no folding at all if every
    panel is its own axes over that band's arcsec extent: a double-click anywhere already
    carries panel-local sky coordinates and is snapped on that band's flux there.
    Nothing is translated within a band, so what you click is what gets snapped;
    the only conversion is the WCS hop from the clicked band to the marked band's frame.

    `zoom_width` (arcsec, from --zoom, default 7) crops every panel to the central square --
    the images of a galaxy-scale lens sit within ~3" of the deflector, and a full 12" (or
    20") stamp makes them tiny. Display only; the snap sees the full frame regardless.

    The subtracted panel is the deflector-removed view ('isophote' default, or 'radial'),
    shown per band; each view row answers a different question -- the subtracted row is where
    a faint image inside the galaxy envelope is visible at all, the as-observed row is where
    you judge whether that blob is really there. Each panel is stretched independently.

    Display-only levers (all per band): this band's contaminant mask is blanked (and
    outlined) and dropped from the stretch percentiles and the isophote fit; `mask_center`
    blanks a disc at the stamp centre; `vmax_value` saturates the as-observed panels. None
    of them touches the sci arrays, so the snaps use true flux and the saved positions do
    not depend on the display. Do note the snap is flux-based: a click placed
    inside/adjacent to a blanked region can still snap onto a hidden pixel within
    `search_box_size` -- process_lens checks the result and warns if one did.
    """
    # The snap is our own brightest-pixel disc (snap_click), not al.Clicker's: the Clicker
    # rounded the click to a pixel corner on these even-sided stamps and searched an
    # asymmetric box, so a click could move further one way than the other (2026-11-02, the
    # JWST repo's fix; adopted across both repos in the 2026-11-28 merge).

    # Rows are the views (as-observed, then the subtracted view), columns the bands.
    views = []
    if not subtract or side_by_side:
        views.append('as-observed')
    if subtract == 'isophote':
        views.append('isophote-subtracted (smoothed)')
    elif subtract == 'radial':
        views.append('radial-subtracted')

    n_rows = len(views)
    # The colour column (2026-11-17, --colour default on): the COLOUR-INDEX map
    # (log10(bluest/reddest vs the deflector's own ratio, PSF-matched) in the as-observed
    # row, and a channel-balanced Lupton RGB of the isophote residuals in the subtracted
    # row -- colour is the discriminant between (blue) lensed images and the (red)
    # deflector, which the single-band panels cannot show. The map is positive where a pixel
    # is BLUER than the deflector and is shown on RdBu (NOT RdBu_r -- the reversed map
    # printed the blue arcs RED, 2026-11-24): blue = bluer than the deflector, red = redder,
    # the deflector itself white. Clicks on either colour panel
    # snap on the MARKED band's flux; the panels are drawn on the marked band's arcsec frame.
    # A raw-flux Lupton RGB was tried first and dropped from THIS GUI: on a disk-dominated
    # lens the red deflector owns every pixel and dilutes the knots' colours away
    # (J0822+1828); the raw RGB is still the mask GUIs' companion window, where the question
    # is the raw colour of a contaminant.
    show_colour = colour and len(bands) >= 2
    col0 = 1 if show_colour else 0
    n_cols = len(bands) + col0
    # Cap the window at ~90% of the screen (Sam was having to resize every window,
    # 2026-11-17): a 4-band, 2-row grid at 6"/panel is 24"x13" before the colourbars.
    # The screen size comes from CoreGraphics (logical points = the pixels a figure window
    # is sized in). Do NOT probe via Tk: the default backend here is macosx, whose manager
    # has no Tk window (AttributeError -> the except leaves the cap off -> oversized GUI),
    # a probe figure survives plt.close as a blank 'Figure 1' (no event loop is running yet
    # to destroy it), and a standalone tkinter.Tk() root leaves a stray blank 'tk' window
    # because the macosx event loop never services Tk's destroy (all seen 2026-11-24).
    panel_in = 6.0
    try:
        if sys.platform == 'darwin':
            import ctypes
            cg = ctypes.CDLL('/System/Library/Frameworks/CoreGraphics.framework/CoreGraphics')

            class _CGRect(ctypes.Structure):
                _fields_ = [('x', ctypes.c_double), ('y', ctypes.c_double),
                            ('w', ctypes.c_double), ('h', ctypes.c_double)]

            cg.CGDisplayBounds.restype = _CGRect
            cg.CGDisplayBounds.argtypes = [ctypes.c_uint32]
            bnds = cg.CGDisplayBounds(cg.CGMainDisplayID())
            scr_w, scr_h = bnds.w, bnds.h
            dpi = float(plt.rcParams['figure.dpi'])
            scale = min(1.0, 0.92 * scr_w / (panel_in * n_cols * 1.15 * dpi),
                        0.86 * scr_h / ((panel_in * n_rows + 1) * dpi))
            panel_in = max(2.5, panel_in * scale)
    except Exception:
        pass
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(panel_in * n_cols * 1.15,
                                                      panel_in * n_rows + 1),
                             squeeze=False)
    # Many bands -> capped, small panels: scale every annotation down with the panel size
    # and thin the ticks out, or the tick and colourbar labels overlap (J1032+5322, 4 bands
    # + the colour column, 2026-11-24).
    fs = max(6.5, min(10.0, 1.7 * panel_in))          # tick and axis-label size
    title_fs = fs + 1.0
    ax_band = {}                       # axes -> band dict
    iso_res = {}                       # filt -> isophote residual array (for the colour panel)
    for ci, b in enumerate(bands):
        col = ci + col0
        base = b['base']
        cont = b['contaminants']
        exclude = make_masks.central_disc(np.asarray(base).shape, b['ps'], mask_center)
        if cont is not None:
            exclude = cont if exclude is None else (exclude | cont)

        n_y, n_x = np.asarray(base).shape
        hw = int(n_x / 2) * b['ps']
        ext = [-hw, hw, -hw, hw]
        for row, view in enumerate(views):
            ax = axes[row][col]
            ax_band[ax] = b
            if view == 'as-observed':
                disp = np.asarray(make_masks.stretched_display(
                    base, b['ps'], stretch=stretch, vmin_percent=vmin_percent,
                    vmax_percent=vmax_percent, asinh_a=asinh_a, exclude=exclude,
                    vmax_value=vmax_value).native)
                if cont is not None:
                    # Outline ON the masked pixels (mask_boundary is an INNER boundary), so
                    # the outline survives the interior being blanked to the colourmap floor.
                    disp[make_masks.mask_boundary(cont)] = 1.0
                im = ax.imshow(disp, cmap=DISPLAY_CMAP, extent=ext)
            elif view == 'isophote-subtracted (smoothed)':
                dp = b['deflector']
                res, vmax_div = isophote_model.isophote_display(
                    base, b['ps'], exclude=exclude,
                    centre=dp[:2] if dp else None, refine_centre=dp[2] if dp else None)
                iso_res[b['filt']] = np.array(res, dtype=float)
                disp = iso_res[b['filt']].copy()
                if cont is not None:
                    disp[cont] = np.nan
                im = ax.imshow(disp, cmap='RdBu_r', vmin=-vmax_div, vmax=vmax_div,
                               extent=ext)
                if cont is not None:
                    ax.contour(cont, levels=[0.5], colors='black', linewidths=0.8,
                               extent=ext, origin='upper')
            else:                        # radial-subtracted
                disp = np.asarray(make_masks.stretched_display(
                    make_masks.radial_median_subtract(base, exclude=exclude), b['ps'],
                    stretch=stretch, vmin_percent=vmin_percent, vmax_percent=vmax_percent,
                    asinh_a=asinh_a).native)
                im = ax.imshow(disp, cmap=DISPLAY_CMAP, extent=ext)
            cb = fig.colorbar(im, ax=ax, fraction=0.046)
            cb.locator = MaxNLocator(4)
            cb.update_ticks()
            cb.ax.tick_params(labelsize=fs)
            ax.tick_params(labelsize=fs)
            ax.xaxis.set_major_locator(MaxNLocator(4))
            ax.yaxis.set_major_locator(MaxNLocator(4))
            if mask_center and mask_center > 0:
                ax.add_patch(plt.Circle((0.0, 0.0), mask_center, fill=False, color='white',
                                        lw=1.0, ls='--'))
            if zoom_width:
                h = zoom_width / 2.0
                ax.set_xlim(-h, h)
                ax.set_ylim(-h, h)
            if row == 0:
                ax.set_title(b['filt'] + ('  [stored]' if b['filt'] == marked_filt else ''),
                             fontsize=title_fs)
            if col == 0:
                ax.set_ylabel(f'{view}\narcsec' if n_rows > 1 else 'arcsec', fontsize=fs)
            if row == n_rows - 1:
                ax.set_xlabel('arcsec', fontsize=fs)

    marked = next(b for b in bands if b['filt'] == marked_filt)
    if show_colour:
        n_y, n_x = np.asarray(marked['sci']).shape
        hw = int(n_x / 2) * marked['ps']
        ext_c = [-hw, hw, -hw, hw]
        rgb_bands = [{'filt': b['filt'], 'data': b['sci'], 'hdr': b['hdr'], 'ps': b['ps'],
                      'noise': b['noise'], 'psf': b['psf']} for b in bands]
        names = sorted((b['filt'] for b in bands), key=make_masks.band_pivot_wavelength)
        axc = axes[0][0]
        # As-observed row: the colour-index map (bluer/redder than the deflector). The
        # raw-flux RGB was tried first and is useless on a disk-dominated lens -- the red
        # deflector owns every pixel and dilutes the knots' colours away (J0822+1828,
        # 2026-11-17); it stays as the masks GUIs' companion, where the question is the
        # raw colour of a contaminant.
        cmap, _ = make_masks.colour_index_map(rgb_bands)
        if cmap is not None:
            imc = axc.imshow(cmap, extent=ext_c, cmap='RdBu', vmin=-0.3, vmax=0.3)
            cbc = fig.colorbar(imc, ax=axc, fraction=0.046,
                               label=f'log10({names[0]}/{names[-1]} vs deflector)')
            cbc.ax.tick_params(labelsize=fs)
            cbc.ax.yaxis.label.set_size(fs)
        axc.set_title(f'colour index ({names[0]}/{names[-1]})', fontsize=title_fs)
        axc.tick_params(labelsize=fs)
        axc.xaxis.set_major_locator(MaxNLocator(4))
        axc.yaxis.set_major_locator(MaxNLocator(4))
        ax_band[axc] = marked                # clicks here snap on the marked band's flux
        if n_rows == 1:
            axc.set_xlabel('arcsec', fontsize=fs)
        if n_rows > 1:
            axc.set_ylabel(f'{views[0]}\narcsec', fontsize=fs)
        if zoom_width:
            axc.set_xlim(-zoom_width / 2.0, zoom_width / 2.0)
            axc.set_ylim(-zoom_width / 2.0, zoom_width / 2.0)
        if n_rows > 1:
            axs = axes[1][0]
            if iso_res:
                rgb_sub, _ = make_masks.lupton_rgb(
                    [{'filt': b['filt'], 'data': np.clip(iso_res[b['filt']], 0, None),
                      'hdr': b['hdr'], 'ps': b['ps']} for b in bands if b['filt'] in iso_res],
                    balance=True)
                axs.imshow(rgb_sub, extent=ext_c)
                ax_band[axs] = marked
                ch = (f"B,G,R = {','.join(names)}" if len(names) == 3
                      else f"B..R = {names[0]}..{names[-1]}")
                axs.set_title(f'colour of the isophote residuals ({ch})', fontsize=title_fs)
            else:
                axs.axis('off')
            axs.set_xlabel('arcsec', fontsize=fs)
            axs.set_ylabel(f'{views[1]}\narcsec', fontsize=fs)
            axs.tick_params(labelsize=fs)
            axs.xaxis.set_major_locator(MaxNLocator(4))
            axs.yaxis.set_major_locator(MaxNLocator(4))
            if zoom_width and iso_res:
                axs.set_xlim(-zoom_width / 2.0, zoom_width / 2.0)
                axs.set_ylim(-zoom_width / 2.0, zoom_width / 2.0)

    # tight_layout AFTER every panel, colourbar and title exists; the rect leaves the top
    # ~10% for the 3-line suptitle drawn by title().
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.90))

    sets = [[]]                       # sets[k] = clicks of source k+1 (marked-band frame)
    clicked_filts = [[]]              # parallel per-position record of the clicked band
    active = [0]                      # boxed so the key handler can rebind it
    marked_hdr = next(b['hdr'] for b in bands if b['filt'] == marked_filt)
    history = []                      # (set index, clicked band, artists) per click, newest
                                    # last, for undo

    def title():
        colour = SET_COLOURS[active[0] % len(SET_COLOURS)]
        fig.suptitle((prompt or
                      f"{lens}  --  marking SOURCE {active[0] + 1} ({colour})\n"
                      f"DOUBLE-CLICK each of its images on ANY panel (snap: that band's "
                      f"flux; stored in the {marked_filt} frame)\n"
                      f"Keys 1-9 switch source;  z undoes the last click;  "
                      f"Esc (or close the window) when done")
                     + "   [North down, East left]",
                     fontsize=max(9.0, min(11.0, fs + 3.0)))
        fig.canvas.draw_idle()

    title()
    moved = []

    def onclick(event):
        b = ax_band.get(event.inaxes)
        if b is None or not event.dblclick:
            # Clicks outside the panels (margins, colourbars) are dropped rather than
            # snapped to a nonsense pixel; single clicks do nothing.
            return
        sy, sx = snap_click(b['sci'], b['ps'], event.ydata, event.xdata,
                            search_box_size)      # snapped, clicked band's arcsec frame
        if b['filt'] == marked_filt:
            my, mx = sy, sx
        else:
            (my, mx), = reproject_positions([(sy, sx)], b['hdr'], marked_hdr)
        sets[active[0]].append((my, mx))
        clicked_filts[active[0]].append(b['filt'])
        colour = SET_COLOURS[active[0] % len(SET_COLOURS)]
        artists = []
        draw_cols = ([0] if show_colour else []) + [ci + col0 for ci in range(len(bands))]
        for col_i in draw_cols:
            b2 = bands[col_i - col0]
            if b2['filt'] == marked_filt or col_i == 0 and show_colour:
                py, px = my, mx            # colour panels live on the marked band's frame
            else:
                (py, px), = reproject_positions([(my, mx)], marked_hdr, b2['hdr'])
            for ax in axes[:, col_i]:
                artists += ax.plot(px, py, '+', color=colour, ms=9, mew=1.5, zorder=5,
                                   path_effects=MARK_STROKE)
                artists.append(ax.annotate(f'{active[0] + 1}.{len(sets[active[0]]) - 1}',
                                           (px, py), color=colour, fontsize=10, xytext=(6, 6),
                                           textcoords='offset points', zorder=5,
                                           path_effects=MARK_STROKE))
        history.append((active[0], b['filt'], artists))
        fig.canvas.draw_idle()
        # How far the brightest-pixel snap moved the click. Printed because the snap is
        # the one step between what you aimed at and what gets saved, and a large jump
        # usually means the box reached a brighter neighbour, not a better centre.
        d_as = float(np.hypot(sy - event.ydata, sx - event.xdata))
        moved.append(d_as)                        # arcsec; process_lens divides by the
        note = ('   <-- LARGE, check it'          # marked band's pixel scale
                if d_as > (search_box_size + 0.75) * b['ps'] else '')
        band_note = '' if b['filt'] == marked_filt else f'  [{b["filt"]} panel]'
        print(f"    snap moved {d_as/b['ps']:4.1f} px ({d_as:.3f}\") "
              f"-> ({my:+.3f}, {mx:+.3f})\"{band_note}{note}")

    def undo():
        """Drop the last click, whichever source and band it was on: remove its markers
        from every panel and pop the set entry, the clicked-band record and the snap log."""
        if not history:
            print("    nothing to undo")
            return
        k, filt_, artists = history.pop()
        y_, x_ = sets[k].pop()
        clicked_filts[k].pop()
        moved.pop()
        for a_ in artists:
            a_.remove()
        print(f"    undone: source {k + 1} image {len(sets[k])} at ({y_:+.3f}, {x_:+.3f})\"")
        fig.canvas.draw_idle()

    def onkey(event):
        if event.key == 'escape':     # finish, as Esc does in the mask GUI
            plt.close(fig)
            return
        if event.key in ('z', 'backspace', 'ctrl+z', 'cmd+z'):   # 'z' undoes, as in the mask GUI
            undo()
            return
        if event.key is not None and event.key in '123456789' and len(event.key) == 1:
            active[0] = int(event.key) - 1
            while len(sets) <= active[0]:
                sets.append([])
                clicked_filts.append([])
            print(f"    -> now marking SOURCE {active[0] + 1}")
            title()

    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    kid = fig.canvas.mpl_connect('key_press_event', onkey)
    plt.show()
    fig.canvas.mpl_disconnect(cid)
    fig.canvas.mpl_disconnect(kid)
    plt.close(fig)
    make_masks._drain_tk()   # let Tk remove the window now, not after the per-band PNGs below

    if moved:
        marked_ps = next(b['ps'] for b in bands if b['filt'] == marked_filt)
        print(f"  snap distances: max {max(moved)/marked_ps:.1f} px, "
              f"mean {sum(moved)/len(moved)/marked_ps:.1f} px of the marked band "
              f"(--snap-radius {search_box_size})")
    keep = [i for i, s in enumerate(sets) if s]
    return ([sets[i] for i in keep],
            [clicked_filts[i] for i in keep],
            moved)   # each set a list of (y_arcsec, x_arcsec) in the marked band's frame


def background_renders(sample, lenses, vmax_value=None):
    """Draw the changed lenses' positions QC overlays (--overlays-only) in ONE detached
    background process, at the end of the run, so the GUI never waits on them (the overlay
    is an isophote fit per band, ~20 s a lens -- drawn inline, it left the closed GUI window
    frozen between lenses; make_masks.rerender_overviews explains the same problem). Log in
    data/run_logs/. The band overview does not show positions, so it is not re-rendered."""
    if not lenses:
        return
    import shlex
    import subprocess
    import time
    plt.close('all')
    here = os.path.dirname(os.path.abspath(__file__))
    py = shlex.quote(sys.executable)
    steps = [f'{py} {shlex.quote(os.path.join(here, "make_positions.py"))} --sample {sample} '
             f'--lens {shlex.quote(l)} --overlays-only'
             + (f' --vmax-value {vmax_value}' if vmax_value is not None else '') for l in lenses]
    log = os.path.join(ws_path, 'data', 'run_logs', f'positions_renders_{sample}_{time.strftime("%Y%m%d_%H%M%S")}.log')
    os.makedirs(os.path.dirname(log), exist_ok=True)
    with open(log, 'w') as fh:
        subprocess.Popen(['/bin/sh', '-c', '; '.join(steps) + '; echo "=== renders done"'],
                         stdout=fh, stderr=subprocess.STDOUT, start_new_session=True,
                         env=dict(os.environ, MPLBACKEND='Agg'))
    print(f"\npositions overlays of {len(lenses)} lens(es) rendering IN THE BACKGROUND "
          f"(~20 s each); log: {os.path.relpath(log, ws_path)}")


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
        stamp_dir = filt_dirs[filt]
        prefix = make_masks.find_prefix(stamp_dir, drizzle_pass)
        if prefix is None:
            continue
        json_path = os.path.join(stamp_dir, f'{prefix}_positions.json')
        if not os.path.exists(json_path):
            continue
        sets = positions_io.read_sets(json_path)
        sci_path = stamp_paths.find_stamp(stamp_dir, prefix, 'sci')
        ps = make_masks.pixel_scale_from_header(fits.getheader(sci_path))
        sci_native = al.Array2D.from_fits(file_path=sci_path, pixel_scales=ps).native
        out_png = os.path.join(stamp_dir, f'{prefix}_positions.png')
        save_overlay_png(sci_native, ps, sets, out_png, stretch=stretch,
                         vmin_percent=vmin_percent, vmax_percent=vmax_percent,
                         asinh_a=asinh_a, vmax_value=vmax_value)
        print(f"  rebuilt {out_png}  ({len(positions_io.flatten(sets))} positions"
              f"{f', {len(sets)} sets' if len(sets) > 1 else ''})")
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
    marked_filt = pick_display_filt(sorted(candidates), None, lens=lens, sample=sample) if candidates else None
    if marked_filt is None:
        return 0                      # simply not marked yet -- nothing to say
    marked_dir = filt_dirs[marked_filt]
    marked_prefix = make_masks.find_prefix(marked_dir, drizzle_pass)
    marked_json = os.path.join(marked_dir, f'{marked_prefix}_positions.json')
    if not os.path.exists(marked_json):
        # Provenance says this lens WAS marked but the product is gone. The clicks are not
        # lost -- info/provenance/lens_positions.json keeps `positions_arcsec` -- but restoring them
        # here would be wrong: the skip-if-exists rule means a restored file silently
        # prevents the re-mark the deletion may have been for. Report, don't guess.
        n_rec = prov[marked_filt].get('n_positions')
        print(f"{lens}: provenance records {n_rec} position(s) marked on {marked_filt}, but "
              f"{os.path.relpath(marked_json, ws_path)} is MISSING.\n"
              f"    The clicks survive in {os.path.basename(POSITIONS_JSON)} "
              f"('positions_arcsec'). Re-mark with --force, or restore from there.")
        return 0
    sets = positions_io.read_sets(marked_json)
    marked_hdr = fits.getheader(stamp_paths.find_stamp(marked_dir, marked_prefix, 'sci'))
    print(f"\n{lens}: re-deriving {len(positions_io.flatten(sets))} position(s) "
          f"({len(sets)} set(s)) from {marked_filt}")
    n = 0
    for filt in sorted(filt_dirs):
        if filt == marked_filt:
            continue
        stamp_dir = filt_dirs[filt]
        prefix = make_masks.find_prefix(stamp_dir, drizzle_pass)
        if prefix is None:
            continue
        band_sci = stamp_paths.find_stamp(stamp_dir, prefix, 'sci')
        band_hdr = fits.getheader(band_sci)
        band_ps = make_masks.pixel_scale_from_header(band_hdr)
        new_sets = [reproject_positions(s_, marked_hdr, band_hdr) for s_ in sets]
        new = positions_io.flatten(new_sets)
        json_path = os.path.join(stamp_dir, f'{prefix}_positions.json')
        old = None
        if os.path.exists(json_path):
            old = positions_io.flatten(positions_io.read_sets(json_path))
        positions_io.write_sets(new_sets, json_path)
        band_native = al.Array2D.from_fits(file_path=band_sci, pixel_scales=band_ps).native
        save_overlay_png(band_native, band_ps, new_sets,
                         os.path.join(stamp_dir, f'{prefix}_positions.png'),
                         stretch=stretch, vmin_percent=vmin_percent,
                         vmax_percent=vmax_percent, asinh_a=asinh_a, vmax_value=vmax_value)
        moved = ('' if old is None or len(old) != len(new) or not new else
                 f"  moved {1000*max(np.hypot(a-c, b-d) for (a, b), (c, d) in zip(new, old)):.1f} mas"
                 f" ({max(np.hypot(a-c, b-d) for (a, b), (c, d) in zip(new, old))/band_ps:.2f} px)")
        print(f"  {filt}: rewritten{moved}")
        n += 1
    return n


def process_lens(lens, filt_dirs, sample, requested_filt, drizzle_pass, force,
                 search_box_size, display, stretch, vmin_percent, vmax_percent, asinh_a,
                 mask_center=0.0, vmax_value=None, subtract='isophote',
                 side_by_side=True, hide_masked=True, all_bands=True, zoom=True,
                 zoom_width=7.0, colour=True):
    """Mark positions once for one lens and broadcast the result to every band.

    `filt_dirs` maps filt -> stamp_dir for that lens, from make_masks.discover_targets.
    With `all_bands` (default) every band is loaded and shown -- one column per band -- and
    a click on any panel snaps on that band's flux and is stored in the MARKED band's frame
    through the WCS. Returns True if positions were written.
    """
    filts = sorted(filt_dirs)
    display_filt = pick_display_filt(filts, requested_filt, lens=lens, sample=sample)
    if display_filt is None:
        print(f"{lens}: requested --filt {requested_filt} not present "
              f"(have {', '.join(filts)}), skipping")
        return False

    display_dir = filt_dirs[display_filt]                # stamp dir of the marked band
    display_prefix = make_masks.find_prefix(display_dir, drizzle_pass)
    if display_prefix is None:
        print(f"{lens} {display_filt}: no stamp sci for --pass {drizzle_pass}, skipping")
        return False

    # Skip if the marked band already carries positions (one file per lens), unless --force.
    marked_json = os.path.join(display_dir, f'{display_prefix}_positions.json')
    if os.path.exists(marked_json) and not force:
        print(f"{lens}: positions already exist ({display_filt}), skipping (--force to re-mark)")
        return False

    # Load the panels: every band with --all-bands, else the marked one. Columns are in
    # strict pivot-WAVELENGTH order with the marked band in place (tagged [stored] in its
    # title) -- Sam's call over the old marked-band-first ordering (2026-11-24).
    shown = sorted(filts, key=make_masks.band_pivot_wavelength) if all_bands else [display_filt]
    bands = []
    for filt in shown:
        info = load_band_panel_info(lens, filt, filt_dirs[filt], drizzle_pass, display,
                                    hide_masked)
        if info is not None:
            bands.append(info)
    marked = next(b for b in bands if b['filt'] == display_filt)  # verified present above
    pixel_scales = marked['ps']
    contaminants = marked['contaminants']

    zoom_arcsec = None
    if zoom:
        if mast_target_names.category(lens, sample) == 'multi_galaxy':
            print(f"  NOTE: multi_galaxy lens -- zoom left OFF (the group needs the full "
                  f"stamp); --zoom-width applies to the rest")
        else:
            zoom_arcsec = zoom_width

    print(f"\n{lens}  marking on {display_filt}  [{display_prefix} pass]  "
          f"pixel_scale={pixel_scales:.4f}\"/pix")
    print(f"  {marked['sci_path']}")
    if all_bands and len(bands) > 1:
        print(f"  all-bands display: {', '.join(b['filt'] for b in bands)} -- click any "
              f"panel; the snap runs on THAT band's flux and the position is stored in "
              f"the {display_filt} frame through the WCS")
    for b in bands:
        if b['contaminants'] is not None:
            print(f"  {b['filt']}: {int(b['contaminants'].sum())}px contaminant-masked -- "
                  f"BLANKED in the display (outlined) and dropped from the stretch and "
                  f"profile; --no-hide-masked shows them")
        elif hide_masked:
            print(f"  NOTE: no {b['prefix']}_mask.fits on {b['filt']} -- nothing to hide "
                  f"(draw one with make_masks.py)")

    source_label = marked['source_label']
    source_label += (f', as-observed + {subtract}-subtracted' if subtract and side_by_side
                     else f', {subtract}-subtracted' if subtract else '')
    if contaminants is not None:
        source_label += ' (contaminants blanked)'

    print(f"  display: {source_label}, {stretch} stretch  |  snap disc radius "
          f"= {search_box_size} px")
    if zoom_arcsec:
        print(f"  panels zoomed to the central {zoom_arcsec:g}\" (--no-zoom for the full "
              f"stamp); display only, the snap sees the whole frame")
    if mask_center and mask_center > 0:
        print(f"  central {mask_center:g}\" blanked in the DISPLAY (deflector light hidden "
              f"and dropped from the stretch percentiles); snap still uses true flux")
    if vmax_value is not None:
        print(f"  display saturated above {vmax_value:g} ({source_label})")
    if subtract == 'isophote':
        print("  deflector's ISOPHOTE model subtracted from the DISPLAY, smoothed, red = light "
              "the model lacks (finding aid only)")
    elif subtract == 'radial':
        print("  deflector's radial-median profile subtracted from the DISPLAY "
              "(finding aid only -- expect a quadrupole residual on an elliptical)")
    print("  DOUBLE-CLICK each lensed image of the source (2 for a double, 4 for a quad).")
    print("  Each click snaps to the brightest nearby pixel. Esc (or close the window) when done.")

    sets, clicked_filts, moved = mark_positions_gui(
        bands, display_filt, lens, search_box_size,
        display=display, stretch=stretch,
        vmin_percent=vmin_percent, vmax_percent=vmax_percent, asinh_a=asinh_a,
        mask_center=mask_center, vmax_value=vmax_value,
        subtract=subtract, side_by_side=side_by_side, zoom_width=zoom_arcsec,
        colour=colour)

    positions = positions_io.flatten(sets)
    if len(positions) == 0:
        print(f"  no positions clicked -- nothing written for {lens}")
        return False
    for k, s_ in enumerate(sets):
        if len(s_) < 2:
            print(f"  WARNING: source {k + 1} has only {len(s_)} position -- a lensed source "
                  f"has >=2 images. Re-run with --force if this was a mis-click.")
    if len(sets) > 1:
        print(f"  {len(sets)} position sets: " + ', '.join(
            f'source {k + 1} = {len(s_)}' for k, s_ in enumerate(sets)))

    # The snap uses true flux, so a click at the edge of a blanked region can still land
    # inside it -- the one way a display-only lever can reach the saved product. Check rather
    # than trust, and name the offenders: a position on a contaminant is a bad position.
    # With all-bands the relevant mask is the CLICKED band's (that is the panel the click
    # was aimed in), so each position is tested against the mask of the band it came from.
    n_on_mask = 0
    band_by_filt = {b['filt']: b for b in bands}
    for (y, x), cf in zip(positions, positions_io.flatten(clicked_filts)):
        b = band_by_filt.get(cf, marked)
        cont = b['contaminants']
        if cont is None:
            continue
        if b['filt'] == display_filt:
            py, px = y, x
        else:
            (py, px), = reproject_positions([(y, x)], marked['hdr'], b['hdr'])
        iy = int(round((np.asarray(cont).shape[0] - 1) / 2.0 - py / b['ps']))
        ix = int(round((np.asarray(cont).shape[1] - 1) / 2.0 + px / b['ps']))
        if (0 <= iy < cont.shape[0] and 0 <= ix < cont.shape[1] and cont[iy, ix]):
            n_on_mask += 1
            print(f"  WARNING: position ({y:+.3f}, {x:+.3f})\" (clicked on {cf}) snapped "
                  f"ONTO a contaminant-masked pixel in {b['filt']} -- the snap uses true "
                  f"flux, so it can reach into a blanked region. Re-run with --force to "
                  f"re-mark it.")

    sci_path = marked['sci_path']
    marked_hdr = fits.getheader(sci_path)

    # Broadcast THROUGH SKY, not by copying arcsec (see reproject_positions). Each band gets
    # its own Grid2DIrregular in its own frame, one file per band into that band's stamp
    # dir, under that dir's own pass prefix (f160W has no CR pass -> stamp12_positions.json).
    written = []
    for filt in filts:
        stamp_dir = filt_dirs[filt]
        prefix = make_masks.find_prefix(stamp_dir, drizzle_pass)
        if prefix is None:
            continue
        band_sci = stamp_paths.find_stamp(stamp_dir, prefix, 'sci')
        band_hdr = fits.getheader(band_sci)
        band_ps = make_masks.pixel_scale_from_header(band_hdr)
        if filt == display_filt:
            band_sets, shift_mas = sets, 0.0
        else:
            band_sets = [reproject_positions(s_, marked_hdr, band_hdr) for s_ in sets]
            shift_mas = 1000.0 * max(np.hypot(by - y, bx - x)
                                     for (by, bx), (y, x) in zip(positions_io.flatten(band_sets),
                                                                 positions))

        # Does this band's own image actually sit where the marked band says it does? A
        # broadcast cannot fix a misregistered band; it can only be honest about it. The
        # check is FIELD SOURCES (band_tie), not a deflector centroid: a deflector estimate
        # reads arc light as an offset wherever the arcs rival the lens light.
        tie = None
        if filt != display_filt:
            tie = band_tie(display_dir, display_prefix, stamp_dir, prefix)
            if tie['mas'] is None:
                print(f"  {filt}: band tie unmeasured (too few field sources and the core "
                      f"centroid failed)")
            elif tie['mas'] / (band_ps * 1000.0) > 1.0:
                print(f"  WARNING: {filt} is MISREGISTERED against {display_filt} by "
                      f"{tie['mas']:.0f} mas = {tie['mas']/(band_ps*1000.0):.2f} {filt} px "
                      f"({tie['method']}, n={tie['n']}). The broadcast positions are correct "
                      f"on SKY but will sit that far off this band's images. Fix the "
                      f"astrometry, don't re-mark.")

        json_path = os.path.join(stamp_dir, f'{prefix}_positions.json')
        positions_io.write_sets(band_sets, json_path)
        # The QC overlay PNGs are NOT drawn here: an isophote fit per band between one GUI
        # and the next froze the window ~20 s a lens. They are rebuilt from the JSONs just
        # written, in one detached process at the end of the run (background_renders).
        written.append({'filt': filt, 'prefix': prefix,
                        'reprojection_shift_mas': round(shift_mas, 2),
                        'band_tie_mas': (None if not tie else tie['mas']),
                        'band_tie_px': (None if not tie or tie['mas'] is None
                                        else round(tie['mas'] / (band_ps * 1000.0), 3)),
                        'band_tie_method': (None if not tie else tie['method']),
                        'band_tie_n_sources': (None if not tie else tie['n'])})
        extra = '' if filt == display_filt else f'  (+{shift_mas:.0f} mas vs a raw copy)'
        print(f"  wrote {json_path}{extra}")

    info_json.update(POSITIONS_JSON, sample, lens, display_filt, {
        'n_positions': len(positions),
        # Every position, flattened, so one-set records read exactly as before; a multi-set
        # lens also carries the split (the JSON file holds one Grid2DIrregular per set).
        'positions_arcsec': [[round(y, 6), round(x, 6)] for (y, x) in positions],
        **({'n_sets': len(sets),
            'positions_arcsec_by_set': [[[round(y, 6), round(x, 6)] for (y, x) in s_]
                                        for s_ in sets]} if len(sets) > 1 else {}),
        'marked_filt': display_filt,
        # Which band's panel each position was CLICKED in (the marked band's frame is what
        # is stored; a click on another band snapped on that band's flux and hopped here
        # through the WCS). Same flattened order as positions_arcsec.
        'clicked_filts': positions_io.flatten(clicked_filts),
        'display_bands': [b['filt'] for b in bands],
        'marked_prefix': display_prefix,
        'drizzle_pass': make_masks.stamp_pass_of(display_dir, display_prefix) or 'single',
        'pixel_scale_arcsec': round(pixel_scales, 6),
        'snap_radius_px': search_box_size,
        'snap_moved_px_max': (round(max(moved)/pixel_scales, 2) if moved else None),
        'snap_moved_px_mean': (round(sum(moved)/len(moved)/pixel_scales, 2) if moved else None),
        'broadcast': 'sky (WCS pixel->world->pixel per band)',
        'display': source_label,
        'display_mask_center_arcsec': mask_center or None,
        'display_vmax_value': vmax_value,
        'display_subtract': subtract or None,
        'display_subtract_radial': subtract == 'radial',
        'display_side_by_side': bool(subtract and side_by_side),
        'display_hidden_contaminant_px': (int(contaminants.sum())
                                          if contaminants is not None else None),
        'n_positions_on_contaminant': n_on_mask,
        'broadcast_to': written,
    })
    # One record per lens, keyed on the MARKED band: a re-mark on a different band than the
    # record was made under (e.g. an old deliberately-empty record on f814W, re-marked on the
    # primary band by the all-bands GUI) would otherwise leave the stale record in place,
    # describing a product the broadcast just overwrote.
    for dropped in info_json.prune_lens_keys(POSITIONS_JSON, sample, lens, display_filt):
        print(f"  NOTE: dropped the stale {dropped} record from "
              f"{os.path.basename(POSITIONS_JSON)} (superseded by the {display_filt} mark)")
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
                   help='force which band you mark on (e.g. f606W); default the best '
                        'available band per lens (f475X>f814W>f606W>f555W>f160W>...). The marked '
                        'positions are broadcast to every band regardless')
    p.add_argument('--pass', dest='drizzle_pass', choices=['auto', 'cr', 'nocrrej'],
                   default='auto',
                   help="which stamp pass to mark on, matching make_stamps.py/make_masks.py: "
                        "'auto' (default) uses the band's science stamp (stamp<S>_*: the CR pass, "
                        "or the no-CR pass where none exists, F160W); 'cr'/'nocrrej' ask for one pass")
    p.add_argument('--size', type=float, default=None,
                   help='stamp size in arcsec (default: each lens\'s modelling size, its stamp_size in '
                        'info/samples/lens_samples.json, else 12"). Only modelling-size stamps are '
                        'tracked, so another size is skipped unless --include-ignored')
    p.add_argument('--force', action='store_true', default=False,
                   help='re-mark a lens that already has positions (default: skip it)')
    p.add_argument('--snap-radius', '--search-box-size', dest='snap_radius', type=int,
                   default=SNAP_RADIUS_DEFAULT,
                   help='radius in pixels of the brightest-pixel snap DISC around each '
                        'double-click (snap_click; default 1). The snap was al.Clicker\'s '
                        'box until 2026-11-28: it rounded the click to a pixel corner on '
                        'these even-sided stamps and searched an asymmetric box, so a click '
                        'could move further one way than the other. 1 snaps only to the '
                        'immediate neighbours; use 2 to catch the peak from a coarse click, '
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
    p.add_argument('--vmin-percent', type=float, default=0.5,
                   help='lower clip percentile for the display stretch (default 0.5, the '
                        'dataset QC subplots\' PercentileInterval(99))')
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
    p.add_argument('--all-bands', action=argparse.BooleanOptionalAction, default=True,
                   help='show EVERY band of the lens, one column per band, each with its '
                        'as-observed and subtracted views (default on, 2026-11-17), and let '
                        'a double-click on ANY panel snap on that band\'s true flux; the '
                        'position is stored in the marked band\'s frame through the two '
                        'bands\' WCS (the same hop the broadcast uses). This is what lets '
                        'you mark an image in the band where it is actually visible '
                        '(e.g. f555W) instead of only the highest-priority one. '
                        '--no-all-bands shows the marked band alone, the pre-2026-11-17 '
                        'behaviour. Display-only difference; the saved positions are the '
                        'same sky either way')
    p.add_argument('--zoom', action=argparse.BooleanOptionalAction, default=True,
                   help='crop every GUI panel to the central --zoom-width arcsec square '
                        '(default on). The images of a galaxy-scale lens sit within ~3" of '
                        'the deflector, so a full 12" (or 20") stamp makes them tiny on '
                        'screen. Display only -- the snap sees the whole frame regardless. '
                        'Automatically OFF on multi_galaxy lenses (the group needs the '
                        'full stamp), with a printed note')
    p.add_argument('--zoom-width', type=float, default=7.0, metavar='ARCSEC',
                   help='side of the zoomed display square in arcsec (default 7.0; '
                        '--zoom only)')
    p.add_argument('--colour', action=argparse.BooleanOptionalAction, default=True,
                   help='show a Lupton RGB colour panel of all the lens\'s bands as the '
                        'left column -- the raw fluxes in the as-observed row and the '
                        'isophote residuals in the subtracted row (default on; 2026-11-17). '
                        'Colour is the discriminant between blue lensed images and the red '
                        'deflector. Clickable like any panel (snaps on the marked band\'s '
                        'flux). Inactive on a single-band lens. --no-colour disables')
    p.add_argument('--hide-masked', action=argparse.BooleanOptionalAction, default=True,
                   help="BLANK this band's already-masked contaminants "
                        '(stamp12_mask.fits) in the DISPLAY, outlining what was hidden, '
                        'and drop them from both the stretch percentiles and the radial '
                        'profile (default on). A neighbour bright enough to be masked is also '
                        'bright enough to own the colour scale, to be mis-clicked as a lensed '
                        'image, and to stamp a dark ring across the arcs at its own radius. '
                        'Display only -- the snap still uses true flux, and a position that '
                        'lands on a masked pixel anyway is reported. --no-hide-masked shows '
                        'the band untouched; silently inactive where the band has no mask')
    p.add_argument('--subtract', choices=['isophote', 'radial', 'none'], default=None,
                   help='how the deflector is removed in the subtracted panel: isophote '
                        '(DEFAULT since 2026-11-01: elliptical isophotes free in ellipticity '
                        'and PA, shown smoothed on a red-blue map), radial (the old '
                        'radial-median profile) or none (the raw image only). Display only')
    p.add_argument('--subtract-radial', action=argparse.BooleanOptionalAction, default=None,
                   help="subtract the deflector's azimuthally-averaged radial profile from "
                        'the DISPLAYED image (DEFAULT ON -- it is what makes the images '
                        'clickable at all on most lenses). The strongest of the three '
                        'arc-finding levers: unlike --mask-center it reveals images buried '
                        'inside the galaxy envelope rather than hiding them. Expect a '
                        'quadrupole residual (the deflector is elliptical, not circular) '
                        'and treat it as a finding aid, not photometry. Display only. Kept for '
                        'old command lines: --subtract-radial = --subtract radial, '
                        '--no-subtract-radial = --subtract none')
    p.add_argument('--include-ignored', action='store_true', default=False,
                   help='also offer bands whose stamp directory is GITIGNORED. Off by '
                        'default: such a band is not a science product (the three 420 s '
                        'SLACS SNAP f814W diagnostics, since deleted, are why the rule exists), and hand-marked '
                        'positions written there are non-regenerable and invisible to every '
                        'clone. This guard exists because it already happened -- J1538+5817 '
                        'was MARKED on its SNAP f814W band and broadcast from there across a '
                        '688 mas (13.8 px) misregistration; both copies were deleted 2026-09-22.')
    a = p.parse_args()
    mast_target_names.resolve_sample(p, a, default=mast_target_names.DEFAULT_SAMPLE)
    sample = a.sample

    make_masks.STAMP_SIZE = a.size
    root = stamp_paths.stamps_root(ws_path)

    # discover_targets yields per (lens, filt); regroup to filt -> stamp_dir per lens so we
    # mark once and broadcast across the lens's bands.
    lens_filts = {}
    for lens, filt, stamp_dir in make_masks.discover_targets(
            root, sample, a.lens, None,       # filt=None: consider all bands, pick best below
            include_ignored=a.include_ignored):
        lens_filts.setdefault(lens, {})[filt] = stamp_dir

    if not lens_filts:
        raise SystemExit(f"no stamps found under {root} for sample {sample} matching "
                         f"lens={a.lens!r}")

    print(f"{len(lens_filts)} lens(es) to process")

    if a.rebroadcast:
        total = 0
        for lens in sorted(lens_filts):
            total += rebroadcast(lens, lens_filts[lens], sample, a.drizzle_pass, a.stretch,
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

    # --subtract wins; the old --(no-)subtract-radial still maps onto it; default isophote.
    subtract = (a.subtract if a.subtract is not None
                else 'radial' if a.subtract_radial is True
                else 'none' if a.subtract_radial is False else 'isophote')
    subtract = None if subtract == 'none' else subtract
    made = skipped = 0
    changed = []
    for lens in sorted(lens_filts):
        if process_lens(lens, lens_filts[lens], sample, a.filt, a.drizzle_pass, a.force,
                        a.snap_radius, a.display, a.stretch,
                        a.vmin_percent, a.vmax_percent, a.asinh_a,
                        a.mask_center, a.vmax_value, subtract,
                        a.side_by_side, a.hide_masked, a.all_bands, a.zoom, a.zoom_width,
                        a.colour):
            made += 1
            changed.append(lens)
        else:
            skipped += 1

    print(f"\nDone: {made} lens(es) marked, {skipped} skipped")
    background_renders(sample, changed, a.vmax_value)


if __name__ == '__main__':
    main()
