#!/usr/bin/env python
"""
Build a per-lens, per-filter PSF from the drizzled mosaics for strong-lens modelling.

Empirical route (default): detect isolated field stars in the full drizzled science
frame, extract postage stamps, and build an oversampled effective PSF (ePSF) with
photutils (DAOStarFinder -> extract_stars -> EPSFBuilder). This is a production port of
scripts/old_notebooks/load_data.ipynb, with every manual step of that notebook replaced
by automatic quality cuts plus a per-lens override file:

  * the hand-typed NaN rectangles that blanked the lens galaxy / bad regions  ->  an
    automatic mask (deflector circle from info/slacs_coords.py, zero-weight pixels from
    the drizzle WHT map, and a frame-edge margin);
  * the frame-max detection threshold ->  a background-sigma threshold (threshold_scale x
    the sigma-clipped background std);
  * the hand-deletion of stars 5..48 ->  isolation / sharpness / roundness / bounds cuts,
    keep the brightest N, plus an optional info/psf_stars.json entry (include / exclude
    coords, or per-lens parameter overrides) for the few fields the automatics get wrong.

Model route (fallback): when a field has too few usable stars (< --min-stars) or with
--method model, evaluate an STScI STDPSF library grid at the deflector position via
scripts/psf_models.py. Covers ACS/WFC, WFC3/IR, WFPC2/WF3 and WFC3/UVIS (BELLS GALLERY).

Outputs (data/psf/<sample>/<lens>/<filt>/):
  * psf_kernel.fits  -- image-scale, odd, unit-sum: a drop-in al.Kernel2D at the band
                        pixel scale (built by block_reduce of the in-memory oversampled
                        ePSF, which is not itself written to disk -- nothing reads it back).
  * psf.png          -- QA panel: selected-star montage, the kernel, and a radial profile.

For a MODEL-tier build (method_used starts with 'model'), this auto-chains into
make_psf_inject.run_injection(..., promote=True): the drizzle-broadened injected kernel
becomes the canonical psf_kernel.fits / cutout_[cr_]psf.fits, and the analytic model this
function just wrote is moved aside to psf_kernel_analytic.fits / cutout_[cr_]psf_analytic.fits
(comparison only). See make_psf_inject.py's docstring. If injection fails (e.g. the
data/drizzle_files/ inputs were cleared), it's a graceful degradation: a warning prints and
the analytic model stays canonical.

Records the build in info/lens_psf.json (method, n_stars, fwhm_pix, oversample,
kernel_size); no data for a filter records null and exits 0, matching the other scripts.

Usage:
    conda run -n stenv python scripts/make_psf.py --lens J0252+0039 --filt f814W
"""

import argparse
import glob
import json
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.modeling import fitting, models
from astropy.nddata import NDData, block_reduce
from astropy.stats import sigma_clipped_stats
from astropy.table import Table
from astropy.visualization import (AsinhStretch, ImageNormalize, LinearStretch,
                                   LogStretch, PercentileInterval)
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
import astropy.units as u
from photutils.background import Background2D, MedianBackground
from photutils.centroids import centroid_com
from photutils.detection import DAOStarFinder
from photutils.psf import EPSFBuilder, extract_stars

ws_path = '/Users/samlange/Code/HST_lenses_data_reduction'
sys.path.insert(0, os.path.join(ws_path, 'info'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mast_target_names
# Reuse the canonical product locator and the coordinate resolver rather than
# re-implementing the instrument-varying prefix/suffix globbing or the per-sample coords
# table lookup (SLACS in slacs_coords.py, BELLS GALLERY in gallery_coords.py). Importing
# make_cutouts is safe: its work is under __main__.
from make_cutouts import find_products, _has_products, catalogue_coord_for
import psf_models
import info_json


# ── Per-instrument defaults ─────────────────────────────────────────────────────
# Keyed on a normalised instrument label (see instrument_key). WFPC2's drizzled primary
# header records DETECTOR=PC (the aperture, not the chip -- the lens is on WF3), so WFPC2
# is keyed on INSTRUME alone. Every value here is overridable by info/psf_stars.json and
# then by an explicit CLI flag; the drizzled pixel scale is read from the WCS, not stored.
_BASE = dict(
    threshold_scale=5.0,     # detection threshold = this x sigma-clipped background std
    min_stars=3,             # fewer usable stars than this -> model fallback
    kernel_size=31,          # odd size of the delivered image-scale kernel (pixels)
    lens_mask_radius=5.0,    # arcsec circle around the deflector, masked out of detection
    maxiters=15,             # EPSFBuilder iterations
    fwhm_tol_lo=0.75,        # per-candidate fitted FWHM must exceed this x psf_fwhm
    fwhm_tol_hi=1.4,         # ...and stay below this x psf_fwhm (rejects galaxies / CRs)
    max_ellip=0.15,          # reject elongated objects (edge-on galaxies, streaks)
    min_snr=30.0,            # absolute peak-S/N floor for a PSF-grade star. A 5-sigma
                             # DAO detection is not a PSF star: on star-poor fields the
                             # only round detections are ~5-sigma noise blobs (measured
                             # peak-S/N 4-6 on J0252 F606W) that build a pure-noise ePSF.
                             # This absolute gate drops them so the field goes to the
                             # model, while real stars (peak-S/N hundreds) are kept.
    flux_floor_frac=0.05,    # drop stars fainter than this x the brightest kept star:
                             # EPSFBuilder normalises each star by its flux, so a faint
                             # star's background noise is amplified into the ePSF wings.
    pedestal_bad=3.0e-3,     # empirical wing pedestal (median of the outer annulus / peak)
                             # above this -> poor build, fall back to the model. Chosen to
                             # pass every ACS/F555W empirical build (worst ~2.6e-4) and the
                             # good F160W builds (<=1.2e-3) while dropping bad ones (J0936
                             # 6.6e-3). See wing_stats / the hybrid gate in main().
    scatter_bad=3.0e-3,      # empirical wing scatter (std of the outer annulus / peak)
                             # above this -> poor build. Worst clean ACS ~2.0e-3; catches
                             # the noisy F160W builds (J0936 7.8e-3, J0946 6.6e-3).
)
# `fwhm` is the DAOStarFinder detection kernel; `psf_fwhm` is the true stellar PSF FWHM
# (pixels, at the drizzled scale) that the per-candidate shape cut selects around --
# distinct because detection tolerates a coarser kernel than the star/galaxy separation.
# oversample=2: SLACS fields yield only ~10-20 usable stars, so a 4x-oversampled ePSF
# has too many pixels per constraining star and comes out noisy. 2x is well-matched to
# these (near-)critically-sampled drizzled grids and much more robust; raise per lens via
# --oversample or info/psf_stars.json where a field is unusually star-rich.
_INSTR = {
    'ACS/WFC': dict(fwhm=2.5, psf_fwhm=2.0, oversample=2, star_size=51, max_stars=25,
                    min_sep=25.0, sharplo=0.3, sharphi=1.0, roundlo=-0.7, roundhi=0.7),
    'WFC3/IR': dict(fwhm=2.8, psf_fwhm=2.6, oversample=2, star_size=41, max_stars=25,
                    min_sep=20.0, sharplo=0.2, sharphi=1.0, roundlo=-0.7, roundhi=0.7),
    # star_size=35 (smaller than ACS): WFPC2 fields are star-poor and EPSFBuilder diverges
    # on a 51px stamp with only ~3 stars but converges on 35px (verified J0252+0039 F606W).
    'WFPC2':   dict(fwhm=2.5, psf_fwhm=2.0, oversample=2, star_size=35, max_stars=25,
                    min_sep=25.0, sharplo=0.3, sharphi=1.0, roundlo=-0.7, roundhi=0.7),
    # WFC3/UVIS (BELLS GALLERY): a two-CCD optical detector like ACS/WFC, drizzled at the
    # native 0.0396"/px. The measured drizzled stellar FWHM is ~2.4px (F606W), so psf_fwhm
    # is 2.4 (the [0.75,1.4]x shape window then brackets the ~2.2-2.6px star locus with
    # margin) and the DAO detection kernel 2.8. The wide UVIS field (~162") is star-rich, so
    # the empirical build is the norm. flux_floor_frac is lowered to 0.02 (from the 0.05
    # base): these fields routinely contain one very bright (near-saturated) star whose flux
    # sets the 5%-floor high enough to discard a dozen genuinely good, high-S/N mid-range
    # stars -- the "brightest star collapses the flux floor" trap. 0.02 keeps them.
    'WFC3/UVIS': dict(fwhm=2.8, psf_fwhm=2.4, oversample=2, star_size=51, max_stars=25,
                      min_sep=25.0, sharplo=0.3, sharphi=1.0, roundlo=-0.7, roundhi=0.7,
                      flux_floor_frac=0.02),
}


def instrument_key(hdr):
    """Normalised instrument label for the defaults table.

    WFPC2's drizzled primary header carries DETECTOR=PC (the full-field aperture centres
    the target on WF3; DETECTOR names the aperture, not the chip), so key WFPC2 on
    INSTRUME alone. ACS/WFC3 use INSTRUME/DETECTOR.
    """
    inst = hdr.get('INSTRUME', '').strip()
    det = hdr.get('DETECTOR', '').strip()
    if inst == 'WFPC2':
        return 'WFPC2'
    return f'{inst}/{det}'


def resolve_params(inst_key, json_entry, cli):
    """Merge parameter sources: instrument defaults < JSON override < explicit CLI flag."""
    if inst_key not in _INSTR:
        raise KeyError(f'no PSF defaults for instrument {inst_key!r}; add it to _INSTR')
    params = {**_BASE, **_INSTR[inst_key]}
    # JSON override: any key in params that the per-lens entry sets.
    for k in list(params):
        if k in json_entry:
            params[k] = json_entry[k]
    # CLI override: argparse leaves tunables as None unless the user passed them.
    for k in list(params):
        v = getattr(cli, k, None)
        if v is not None:
            params[k] = v
    return params


# ── Masking and star selection ──────────────────────────────────────────────────
def build_mask(shape, wht_data, wcs, catalogue_coord, scale_arcsec, lens_mask_radius,
               edge_margin):
    """Boolean mask of pixels excluded from star detection.

    Masks (a) a circle of `lens_mask_radius` arcsec around the deflector so the lens
    galaxy is never mistaken for a star, (b) zero-weight pixels (chip gaps / no coverage)
    from the drizzle WHT map, and (c) a frame-edge margin so extracted stamps stay in
    bounds. Replaces the notebook's hand-typed NaN rectangles.
    """
    ny, nx = shape
    mask = np.zeros(shape, dtype=bool)

    # (a) deflector circle
    px, py = wcs.world_to_pixel(catalogue_coord)
    r = lens_mask_radius / scale_arcsec
    yy, xx = np.ogrid[:ny, :nx]
    mask |= (xx - px) ** 2 + (yy - py) ** 2 <= r ** 2

    # (b) no-coverage pixels
    if wht_data is not None:
        mask |= ~np.isfinite(wht_data) | (wht_data <= 0)

    # (c) edge margin
    m = int(edge_margin)
    if m > 0:
        mask[:m, :] = True
        mask[-m:, :] = True
        mask[:, :m] = True
        mask[:, -m:] = True

    return mask


def _in_any_box(x, y, boxes):
    for b in boxes:
        x0, x1, y0, y1 = b
        if x0 <= x <= x1 and y0 <= y <= y1:
            return True
    return False


def _fit_stamp_shape(data, x, y, psf_fwhm):
    """Fit a 2D Gaussian to a small window around (x, y); return (fwhm, ellip, centered).

    Measured on a window ~6 PSF-widths across: wide enough that a resolved galaxy's
    extended light inflates its fitted FWHM (and so fails the upper bound), which an
    11-pixel window around a compact nucleus does not. `centered` guards against the fit
    locking onto an off-centre companion. Returns None if the fit fails or the window
    falls off the frame.
    """
    win = max(15, int(round(6 * psf_fwhm)) | 1)   # odd
    half = win // 2
    ix, iy = int(round(x)), int(round(y))
    ny, nx = data.shape
    if not (half <= ix < nx - half and half <= iy < ny - half):
        return None
    stamp = data[iy - half:iy + half + 1, ix - half:ix + half + 1].astype(float)
    stamp = stamp - np.median(stamp)          # local background
    yy, xx = np.mgrid[:win, :win]
    g0 = models.Gaussian2D(amplitude=max(stamp[half, half], 1e-6),
                           x_mean=half, y_mean=half,
                           x_stddev=psf_fwhm / 2.355, y_stddev=psf_fwhm / 2.355)
    try:
        fit = fitting.TRFLSQFitter()(g0, xx, yy, stamp, maxiter=200)
    except Exception:
        return None
    sx, sy = abs(fit.x_stddev.value), abs(fit.y_stddev.value)
    if not (np.isfinite(sx) and np.isfinite(sy)) or (sx + sy) <= 0:
        return None
    fwhm = 2.0 * np.sqrt(2.0 * np.log(2.0)) * np.sqrt(sx * sy)
    ellip = abs(sx - sy) / (sx + sy)
    centered = (abs(fit.x_mean.value - half) < 2.0 and abs(fit.y_mean.value - half) < 2.0)
    return fwhm, ellip, centered


def select_stars(data, mask, params, overrides, star_size):
    """Detect and filter isolated field stars; return an (x, y) table plus a QA dict.

    Detection uses DAOStarFinder with a background-sigma threshold and per-instrument
    sharpness/roundness/min_separation cuts. Because SLACS fields are galaxy-rich and the
    brightest objects are usually galaxies (flux-ranking alone picks them), each candidate
    that clears the bounds/mask/exclusion filters is then Gaussian-fit and kept only if its
    FWHM sits near the true stellar PSF FWHM (`psf_fwhm`) with low ellipticity and a
    centred fit -- this is what actually separates stars from galaxies, diffraction
    streaks and cosmic rays. Survivors are ranked by flux and the brightest `max_stars`
    kept.

    `overrides` (from info/psf_stars.json) is applied last: `exclude` points/boxes remove
    detected stars, `include` points force extra positions in.
    """
    hsize = (star_size - 1) // 2
    mean, median, std = sigma_clipped_stats(data, mask=mask, sigma=3.0)
    threshold = params['threshold_scale'] * std
    psf_fwhm = params['psf_fwhm']

    finder = DAOStarFinder(
        threshold=threshold, fwhm=params['fwhm'],
        sharplo=params['sharplo'], sharphi=params['sharphi'],
        roundlo=params['roundlo'], roundhi=params['roundhi'],
        min_separation=params['min_sep'], exclude_border=True)
    sources = finder(data - median, mask=mask)

    qa = {'n_detected': 0 if sources is None else len(sources)}
    if sources is None or len(sources) == 0:
        rows = []
    else:
        ny, nx = data.shape
        exclude_pts = overrides.get('exclude', [])
        exclude_boxes = [b for b in exclude_pts if len(b) == 4]
        exclude_points = [b for b in exclude_pts if len(b) == 2]
        tol = params['fwhm'] * 2.0
        lo, hi = params['fwhm_tol_lo'] * psf_fwhm, params['fwhm_tol_hi'] * psf_fwhm
        min_peak = params['min_snr'] * std

        rows = []
        for s in sources:
            x, y = float(s['xcentroid']), float(s['ycentroid'])
            ix, iy = int(round(x)), int(round(y))
            # bounds: full stamp must be inside the frame
            if not (hsize < x < nx - 1 - hsize and hsize < y < ny - 1 - hsize):
                continue
            # stamp must not overlap the mask (lens galaxy / chip gap / edge)
            sub = mask[iy - hsize:iy + hsize + 1, ix - hsize:ix + hsize + 1]
            if sub.any():
                continue
            # per-lens exclusions
            if _in_any_box(x, y, exclude_boxes):
                continue
            if any((x - ex) ** 2 + (y - ey) ** 2 <= tol ** 2 for ex, ey in exclude_points):
                continue
            # absolute peak-S/N: reject faint noise blobs before the shape fit
            peak = float(data[iy - 2:iy + 3, ix - 2:ix + 3].max())
            if peak < min_peak:
                continue
            # star/galaxy separation by fitted shape
            shape = _fit_stamp_shape(data, x, y, psf_fwhm)
            if shape is None:
                continue
            fwhm_i, ellip_i, centered = shape
            if not (lo <= fwhm_i <= hi and ellip_i <= params['max_ellip'] and centered):
                continue
            rows.append((x, y, float(s['flux'])))

        rows.sort(key=lambda r: r[2], reverse=True)
        if rows:
            floor = params['flux_floor_frac'] * rows[0][2]
            rows = [r for r in rows if r[2] >= floor]
        rows = rows[:int(params['max_stars'])]

    qa['n_auto'] = len(rows)

    # Forced include positions (bypass detection, still bounds-checked).
    ny, nx = data.shape
    n_incl = 0
    for pt in overrides.get('include', []):
        x, y = float(pt[0]), float(pt[1])
        if hsize < x < nx - 1 - hsize and hsize < y < ny - 1 - hsize:
            rows.append((x, y, np.nan))
            n_incl += 1
    qa['n_included'] = n_incl

    tbl = Table()
    tbl['x'] = [r[0] for r in rows]
    tbl['y'] = [r[1] for r in rows]
    return tbl, qa


# ── ePSF build and kernel extraction ────────────────────────────────────────────
def build_epsf(data, stars_tbl, star_size, oversample, maxiters):
    """Extract star stamps and build an oversampled ePSF. Returns (epsf, fitted_stars)."""
    nddata = NDData(data=data)
    stars = extract_stars(nddata, stars_tbl, size=star_size)
    builder = EPSFBuilder(oversampling=oversample, maxiters=maxiters,
                          progress_bar=False)
    epsf, fitted_stars = builder(stars)
    return epsf, stars, fitted_stars


def _crop_centered(arr, out, cx, cy):
    """Crop `arr` to (out, out) centred on integer (cx, cy), zero-padding if needed."""
    half = out // 2
    ny, nx = arr.shape
    y0, x0 = int(round(cy)) - half, int(round(cx)) - half
    result = np.zeros((out, out), dtype=float)
    sy0, sx0 = max(0, y0), max(0, x0)
    sy1, sx1 = min(ny, y0 + out), min(nx, x0 + out)
    result[sy0 - y0:sy1 - y0, sx0 - x0:sx1 - x0] = arr[sy0:sy1, sx0:sx1]
    return result


def validate_epsf(epsf_data, oversample):
    """True if the ePSF looks like a real PSF (finite, positive core peaked near centre).

    With too few or faint stars EPSFBuilder can diverge and drift the peak to a corner
    (seen on a 3-star WFPC2 field). Catching that here lets the caller fall back to the
    model rather than crop an empty region into a garbage kernel.
    """
    arr = np.asarray(epsf_data, dtype=float)
    if not np.all(np.isfinite(arr)) or arr.max() <= 0:
        return False
    ny, nx = arr.shape
    py, px = np.unravel_index(np.argmax(arr), arr.shape)
    # the peak must land within a few image pixels of the centre
    tol = max(4 * oversample, oversample + 2)
    if not (abs(px - nx // 2) <= tol and abs(py - ny // 2) <= tol):
        return False
    # the core must stand well above the ePSF's own outskirt noise (a pure-noise ePSF
    # from junk detections has peak/std of only a few); measure std on the outer frame.
    r = max(ny, nx) // 4
    edge = np.concatenate([arr[:r].ravel(), arr[-r:].ravel(),
                           arr[:, :r].ravel(), arr[:, -r:].ravel()])
    std = np.std(edge)
    return std <= 0 or arr.max() / std >= 15.0


def oversampled_to_kernel(epsf_data, oversample, kernel_size):
    """Bin the oversampled ePSF down to an image-scale, unit-sum, odd kernel.

    Centres on the ePSF core, crops to kernel_size x oversample so the block reduction
    lands an integer number of output pixels, sums each oversample x oversample block
    (flux-conserving for a convolution kernel), then renormalises to unit sum. Falls back
    from a centroid to the peak, then to the geometric centre, if a step is degenerate.
    """
    arr = np.asarray(epsf_data, dtype=float)
    ny, nx = arr.shape
    # locate the core: centroid of a window around the peak, guarding non-finite results
    py, px = np.unravel_index(np.argmax(arr), arr.shape)
    win = max(oversample * 5, 9)
    y0, x0 = max(0, py - win // 2), max(0, px - win // 2)
    cy0, cx0 = centroid_com(arr[y0:y0 + win, x0:x0 + win])
    cx, cy = x0 + cx0, y0 + cy0
    if not (np.isfinite(cx) and np.isfinite(cy)):
        cx, cy = float(px), float(py)

    n = kernel_size * oversample
    for centre in ((cx, cy), (float(px), float(py)), (nx / 2.0, ny / 2.0)):
        sub = _crop_centered(arr, n, *centre)
        kernel = block_reduce(sub, oversample, func=np.sum)
        if kernel.sum() > 0:
            return kernel / kernel.sum()
    raise ValueError('ePSF integrates to <= 0 at every candidate centre; '
                     'cannot normalise a kernel from it')


def trim_kernel_to_amplitude(kernel, threshold=1e-3, min_half=3):
    """Crop a full image-scale kernel to the radius where the azimuthally-averaged PSF
    drops below `threshold` x peak, then renormalise to unit sum.

    This is the amplitude-based extent that matters for a *convolution* kernel: it keeps
    the PSF out to where it still spreads flux above ~this fraction of the peak. It is
    deliberately NOT an enclosed-energy cut -- EE is area-weighted and, being an integral
    of the noisy empirical wing, under-sizes (a 95% EE cut truncates while the PSF is still
    ~1% of peak; see the info/ PSF handbooks and the psf.png log panel). Band-adaptive by
    construction: broad PSFs (F160W, F606W) keep more pixels than the sharp ACS core.
    Returns (odd, centred, unit-sum kernel, its size).
    """
    arr = np.asarray(kernel, float)
    ny, nx = arr.shape
    cy, cx = ny // 2, nx // 2
    yy, xx = np.mgrid[:ny, :nx]
    rint = np.round(np.hypot(xx - cx, yy - cy)).astype(int)
    peak = arr[cy, cx] if arr[cy, cx] > 0 else float(arr.max())
    prof = np.array([arr[rint == k].mean() if (rint == k).any() else 0.0
                     for k in range(rint.max() + 1)])          # azimuthal mean per radius
    below = np.where(prof < threshold * peak)[0]
    half = int(below[0]) if len(below) else int(rint.max())
    half = min(max(min_half, half), cy, cx)                    # never exceed the full kernel
    sub = arr[cy - half:cy + half + 1, cx - half:cx + half + 1]
    s = sub.sum()
    return (sub / s if s > 0 else sub), 2 * half + 1


def wing_stats(kernel, inner_frac=0.75):
    """(pedestal_frac, scatter_frac): median and std of the kernel's outer annulus over the
    central peak. The annulus is radius in [inner_frac*half, half] of the full kernel; a
    healthy ePSF wing sits near zero there. Gates poor empirical builds and quantifies the
    residual background the pedestal subtraction removes.
    """
    arr = np.asarray(kernel, float)
    ny, nx = arr.shape
    cy, cx = ny // 2, nx // 2
    yy, xx = np.mgrid[:ny, :nx]
    rint = np.round(np.hypot(xx - cx, yy - cy)).astype(int)
    peak = arr[cy, cx] if arr[cy, cx] > 0 else float(arr.max())
    half = min(cy, cx)
    ann = (rint >= int(inner_frac * half)) & (rint <= half)
    if not ann.any() or peak <= 0:
        return 0.0, 0.0
    return float(np.median(arr[ann]) / peak), float(arr[ann].std() / peak)


def subtract_pedestal(kernel):
    """Remove the ePSF's residual flat background (median of the outer annulus) and
    renormalise to unit sum. Returns (kernel, pedestal_frac).

    EPSFBuilder leaves a small DC floor in the ePSF wings. Left in, a ~1e-3-of-peak pedestal
    across the kernel becomes several percent of the (renormalised) flux as a spurious
    uniform background -- and it stops trim_kernel_to_amplitude from ever crossing the
    amplitude threshold, so the trimmed kernel caps at the full size. This is a near no-op
    for clean builds (ACS ~1e-4) and matters for oversampled bands (F160W empirical ~1e-3,
    the WFPC2 F606W DB ePSF ~1e-3). Applied to every kernel regardless of method.
    """
    arr = np.asarray(kernel, float)
    ny, nx = arr.shape
    cy, cx = ny // 2, nx // 2
    yy, xx = np.mgrid[:ny, :nx]
    rint = np.round(np.hypot(xx - cx, yy - cy)).astype(int)
    peak = arr[cy, cx] if arr[cy, cx] > 0 else float(arr.max())
    half = min(cy, cx)
    ann = (rint >= int(0.75 * half)) & (rint <= half)
    ped = float(np.median(arr[ann])) if ann.any() else 0.0
    out = arr - ped
    s = out.sum()
    if s > 0:
        out = out / s
    return out, (float(ped / peak) if peak > 0 else 0.0)


def representative_input_cd(inst_key, lens, filt, sample):
    """2x2 CD (deg/detector-pixel) of a contributing exposure, in the detector frame a model
    ePSF lives in -- so a detector-frame model PSF can be rotated into the North-up drizzle
    frame (psf_models._northup_M). Empirical ePSFs need no CD (built from the North-up mosaic).

    WFPC2: the extracted WF3 file under data/drizzle_files (single WF3 chip; carries the
    correct per-visit roll for split-visit products, keyed by the suffixed `filt`). WFC3/IR:
    the calibrated FLT SCI. ACS/WFC: the calibrated FLC first SCI (fallback only -- the
    focus-diverse path reads per-exposure CDs itself). None if none is found.
    """
    def cd_of(hdr):
        if 'CD1_1' not in hdr:
            return None
        return np.array([[hdr['CD1_1'], hdr['CD1_2']],
                         [hdr['CD2_1'], hdr['CD2_2']]], float)

    def first_sci_cd(path):
        with fits.open(path) as h:
            for hd in h:
                if hd.header.get('EXTNAME') == 'SCI':
                    cd = cd_of(hd.header)
                    if cd is not None:
                        return cd
        return None

    if inst_key == 'WFPC2':
        dd = os.path.join(ws_path, 'data', 'drizzle_files', sample, lens, filt)
        cands = [f for f in sorted(glob.glob(os.path.join(dd, 'wf3_*flt.fits')))
                 if all(t not in os.path.basename(f) for t in ('ivm', 'd2im', 'hlet'))]
        for f in cands:
            cd = first_sci_cd(f)
            if cd is not None:
                return cd
        return None

    pat = {'WFC3/IR': '*flt.fits', 'ACS/WFC': '*flc.fits',
           'WFC3/UVIS': '*flc.fits'}.get(inst_key)
    if pat is None:
        return None
    cd_dir = os.path.join(ws_path, 'data', 'calibrated', sample, lens, filt)
    for f in sorted(glob.glob(os.path.join(cd_dir, pat))):
        cd = first_sci_cd(f)
        if cd is not None:
            return cd
    return None


def measure_fwhm(kernel):
    """Approximate FWHM (pixels) of a kernel via a 2D Gaussian fit; None on failure."""
    try:
        ny, nx = kernel.shape
        yy, xx = np.mgrid[:ny, :nx]
        g0 = models.Gaussian2D(amplitude=kernel.max(), x_mean=nx // 2, y_mean=ny // 2,
                               x_stddev=1.5, y_stddev=1.5)
        fit = fitting.TRFLSQFitter()(g0, xx, yy, kernel)
        s = 0.5 * (abs(fit.x_stddev.value) + abs(fit.y_stddev.value))
        return float(2.0 * np.sqrt(2.0 * np.log(2.0)) * s)
    except Exception as exc:
        print(f'  WARNING: FWHM fit failed ({exc}); recording null')
        return None


# ── Outputs ─────────────────────────────────────────────────────────────────────
def write_fits(data, header, path):
    fits.PrimaryHDU(data=data.astype(np.float32), header=header).writeto(path,
                                                                         overwrite=True)
    print(f'  wrote {path}  shape={data.shape}')


def plot_psf(kernel, epsf_data, star_stamps, scale_arcsec, fwhm_pix, output_path, title):
    """QA panel: selected-star montage, the image-scale kernel on linear *and* log
    stretch (core vs faint wings), and its radial profile."""
    fig = plt.figure(figsize=(20, 5.5))
    gs = fig.add_gridspec(1, 4, width_ratios=[1.2, 1, 1, 1])

    # (1) montage of the extracted star stamps
    ax0 = fig.add_subplot(gs[0])
    n = len(star_stamps)
    if n:
        ncol = int(np.ceil(np.sqrt(n)))
        nrow = int(np.ceil(n / ncol))
        h, w = star_stamps[0].shape
        montage = np.full((nrow * h, ncol * w), np.nan)
        for i, st in enumerate(star_stamps):
            r, c = divmod(i, ncol)
            montage[r * h:(r + 1) * h, c * w:(c + 1) * w] = st
        finite = montage[np.isfinite(montage)]
        vmin, vmax = PercentileInterval(99.0).get_limits(finite)
        norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=AsinhStretch(0.1))
        ax0.imshow(montage, norm=norm, origin='upper', cmap='inferno')
    ax0.set_title(f'selected stars (n={n})')
    ax0.set_xticks([]); ax0.set_yticks([])

    # (2) the delivered image-scale kernel -- linear stretch (shows the core)
    ax1 = fig.add_subplot(gs[1])
    peak = float(kernel.max())
    norm_lin = ImageNormalize(vmin=0.0, vmax=peak, stretch=LinearStretch())
    im1 = ax1.imshow(kernel, norm=norm_lin, origin='upper', cmap='inferno')
    ttl = 'psf_kernel (linear)'
    if fwhm_pix is not None:
        ttl += f'  FWHM {fwhm_pix:.2f}px={fwhm_pix * scale_arcsec:.3f}"'
    ax1.set_title(ttl)
    ax1.set_xlabel('pixels')
    fig.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)

    # (3) same kernel -- log stretch (shows the faint wings)
    ax_log = fig.add_subplot(gs[2])
    norm_log = ImageNormalize(vmin=peak * 1e-4, vmax=peak, stretch=LogStretch())
    im_log = ax_log.imshow(kernel, norm=norm_log, origin='upper', cmap='inferno')
    ax_log.set_title('psf_kernel (log, 1e-4..1)')
    ax_log.set_xlabel('pixels')
    fig.colorbar(im_log, ax=ax_log, fraction=0.046, pad=0.04)

    # (4) radial profile of the kernel
    ax2 = fig.add_subplot(gs[3])
    ny, nx = kernel.shape
    yy, xx = np.mgrid[:ny, :nx]
    rr = np.sqrt((xx - nx // 2) ** 2 + (yy - ny // 2) ** 2).ravel()
    order = np.argsort(rr)
    ax2.plot(rr[order] * scale_arcsec, kernel.ravel()[order], '.', ms=3, color='C1')
    ax2.set_yscale('symlog', linthresh=kernel.max() * 1e-3)
    ax2.set_xlabel('radius (arcsec)')
    ax2.set_ylabel('kernel value')
    ax2.set_title('radial profile')
    ax2.grid(alpha=0.3)

    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(output_path, dpi=120, bbox_inches='tight')
    plt.close(fig)
    print(f'  wrote {output_path}')


def lens_products_rootnames(sample, lens, filt_key):
    """Exposure rootnames that reached the drizzle for (lens, filt), from lens_products.json.

    Used by the ACS focus-diverse ePSF model to retrieve a focus-matched ePSF per exposure.
    Returns [] if the file or entry is absent.
    """
    data = info_json.load(os.path.join(ws_path, 'info', 'lens_products.json'))
    val = data.get(sample, {}).get(lens, {}).get(filt_key)
    return list(val) if val else []


# ── Main ────────────────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--lens', default='J0252+0039')
    p.add_argument('--filt', default='f814W')
    p.add_argument('--sample', default=mast_target_names.DEFAULT_SAMPLE,
                   help='sample subdirectory of data/drizzled/ and data/psf/. Defined in '
                        f'info/lens_samples.json (default {mast_target_names.DEFAULT_SAMPLE})')
    p.add_argument('--pass', dest='drizzle_pass',
                   choices=['auto', 'cr', 'nocrrej'], default='auto',
                   help="drizzle pass to detect stars in. 'auto' (default) prefers the "
                        "CR-rejected mosaic (no cosmic rays to mistake for stars) and "
                        "falls back to no-CR (e.g. WFC3/IR F160W, no CR pass).")
    p.add_argument('--method', choices=['auto', 'empirical', 'model'], default='auto',
                   help="'auto' (default) builds an empirical ePSF and falls back to the "
                        "STDPSF model only when usable stars < --min-stars; 'empirical' "
                        "forces the star build (fails if too few); 'model' forces STDPSF.")
    p.add_argument('--trim-threshold', dest='trim_threshold', type=float, default=1e-3,
                   help="amplitude fraction of peak at which the trimmed modelling kernel "
                        "(data/cutouts/<...>_psf.fits) is cut. Default 1e-3 -- an "
                        "amplitude criterion, not enclosed-energy (which under-sizes).")
    # Tunables: None here means 'unset' so instrument defaults / JSON overrides win; a
    # value on the command line takes precedence over both. See resolve_params().
    p.add_argument('--threshold-scale', dest='threshold_scale', type=float, default=None)
    p.add_argument('--fwhm', type=float, default=None)
    p.add_argument('--max-stars', dest='max_stars', type=int, default=None)
    p.add_argument('--min-sep', dest='min_sep', type=float, default=None)
    p.add_argument('--min-stars', dest='min_stars', type=int, default=None)
    p.add_argument('--lens-mask-radius', dest='lens_mask_radius', type=float, default=None,
                   help='arcsec circle around the deflector masked out of detection')
    p.add_argument('--star-size', dest='star_size', type=int, default=None,
                   help='side (pixels) of the extracted star stamps')
    p.add_argument('--oversample', type=int, default=None,
                   help='ePSF build oversampling factor')
    p.add_argument('--kernel-size', dest='kernel_size', type=int, default=None,
                   help='odd side (pixels) of the delivered image-scale kernel')
    p.add_argument('--output', default=None,
                   help='output dir, default data/psf/<sample>/<lens>/<filt>')
    a = p.parse_args()

    drizzled_dir = os.path.join(ws_path, 'data', 'drizzled', a.sample, a.lens, a.filt)
    output_dir = a.output or os.path.join(ws_path, 'data', 'psf', a.sample, a.lens, a.filt)
    psf_json = os.path.join(ws_path, 'info', 'lens_psf.json')

    # No-data outcome: matches the drizzle/cutout scripts -- record null, exit 0.
    if not os.path.isdir(drizzled_dir) or not _has_products(drizzled_dir):
        print(f'=== NO DATA: {a.lens} {a.filt} (no drizzled products in {drizzled_dir})')
        info_json.update(psf_json, a.sample, a.lens, a.filt, None)
        sys.exit(0)

    os.makedirs(output_dir, exist_ok=True)

    # Per-lens overrides (loaded before the pass resolution: an override can redirect
    # which drizzle pass the empirical build reads).
    overrides_all = info_json.load(os.path.join(ws_path, 'info', 'psf_stars.json'))
    overrides = overrides_all.get(a.sample, {}).get(a.lens, {}).get(a.filt, {})

    # Output pass -- names the cutout (cutout_cr_* vs cutout_*) so the PSF matches the
    # science stamp downstream reads. auto -> cr where a CR pass exists.
    if a.drizzle_pass == 'auto':
        has_cr = bool(glob.glob(os.path.join(drizzled_dir, '*_cr_*_sci.fits')))
        drizzle_pass = 'cr' if has_cr else 'nocrrej'
    else:
        drizzle_pass = a.drizzle_pass

    # Star-detection pass -- the mosaic the *empirical* ePSF is actually built from.
    # The LACosmic CR pass flags sharp field-star cores as cosmic rays and masks them,
    # punching a hole in the ePSF core and dropping stars below the SNR/shape gates
    # (stark on star-poor fields, e.g. J0008-0004: 3 usable stars + a core hole from the
    # CR pass vs 11 clean stars from no-CR). So build the empirical ePSF from the
    # least-CR-rejected mosaic available: prefer a no-CR ('*_nocrrej_*') pass when it is
    # on disk (all ACS once drizzled with --nocrrej; WFC3/IR is nocrrej-named already),
    # else the output pass (WFPC2 keeps its cr pass -- no no-CR exists). The point-source
    # PSF shape is pass-independent, so a no-CR-built PSF stays correct for the CR science
    # stamp, whose name still follows the output pass above. A per-lens "psf_star_pass"
    # override forces a specific pass (e.g. "cr" to opt back out).
    if 'psf_star_pass' in overrides:
        star_pass = overrides['psf_star_pass']
    else:
        has_nocr = bool(glob.glob(os.path.join(drizzled_dir, '*_nocrrej_*_sci.fits')))
        star_pass = 'nocrrej' if has_nocr else drizzle_pass
    sci_file, wht_file = find_products(drizzled_dir, star_pass)
    if star_pass == drizzle_pass:
        print(f'{a.lens} {a.filt}  [{drizzle_pass} pass]')
    else:
        print(f'{a.lens} {a.filt}  [build from {star_pass} pass, name as {drizzle_pass}]')
    print(f'  sci: {os.path.basename(sci_file)}')

    with fits.open(sci_file) as hdul:
        sci_hdr = hdul[0].header
        sci_data = np.array(hdul[0].data, dtype=float)
    with fits.open(wht_file) as hdul:
        wht_data = np.array(hdul[0].data, dtype=float)
    wcs = WCS(sci_hdr)

    inst_key = instrument_key(sci_hdr)
    scale = proj_plane_pixel_scales(wcs.celestial)[0] * 3600.0
    print(f'  instrument: {inst_key}   pixel scale: {scale:.4f}"/pix')

    # Resolved parameters.
    if a.method == 'auto' and 'method' in overrides:
        method = overrides['method']
    else:
        method = a.method
    params = resolve_params(inst_key, overrides, a)
    star_size = int(params['star_size'])
    oversample = int(params['oversample'])
    kernel_size = int(params['kernel_size'])
    if kernel_size % 2 == 0:
        kernel_size += 1

    ra, dec = catalogue_coord_for(a.lens)   # SLACS or BELLS GALLERY, whichever holds it
    catalogue_coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')

    # ── Background subtraction and star selection ──────────────────────────────
    mask = build_mask(sci_data.shape, wht_data, wcs, catalogue_coord, scale,
                      params['lens_mask_radius'], edge_margin=(star_size - 1) // 2)
    box = max(32, min(sci_data.shape) // 40)
    try:
        bkg = Background2D(sci_data, box, coverage_mask=mask, exclude_percentile=90.0,
                          bkg_estimator=MedianBackground()).background
        data = sci_data - bkg
    except Exception as exc:
        print(f'  WARNING: Background2D failed ({exc}); using median subtraction')
        _, med, _ = sigma_clipped_stats(sci_data, mask=mask, sigma=3.0)
        data = sci_data - med

    stars_tbl, qa = select_stars(data, mask, params, overrides, star_size)
    n_stars = len(stars_tbl)
    print(f'  stars: {qa["n_detected"]} detected -> {qa["n_auto"]} auto-kept'
          f' + {qa["n_included"]} forced = {n_stars}')

    too_few = n_stars < params['min_stars']
    use_model = (method == 'model') or (method == 'auto' and too_few)
    if method == 'empirical' and too_few:
        raise SystemExit(
            f'ERROR: only {n_stars} usable stars for {a.lens} {a.filt} '
            f'(need >= {params["min_stars"]}). Loosen cuts, add include stars to '
            f'info/psf_stars.json, or use --method model.')

    star_stamps = []
    method_used = None
    epsf_data = None

    # Empirical build first (unless forced to model). Validate it: a diverged ePSF (peak
    # drifted to a corner on a star-poor field) falls back to the model under --method auto.
    if not use_model:
        print(f'  method: EMPIRICAL  oversample={oversample} star_size={star_size}')
        epsf, stars, fitted = build_epsf(data, stars_tbl, star_size, oversample,
                                         params['maxiters'])
        cand = np.asarray(epsf.data, dtype=float)
        if validate_epsf(cand, oversample):
            # Hybrid quality gate: a validated ePSF can still have a noisy / pedestalled wing
            # (star-poor oversampled fields, esp. F160W). Measure the wing and, under auto,
            # fall back to the model rather than ship a poor empirical kernel.
            k_prov = oversampled_to_kernel(cand, oversample, star_size)
            ped_frac, scat_frac = wing_stats(k_prov)
            poor = ped_frac > params['pedestal_bad'] or scat_frac > params['scatter_bad']
            if poor and method == 'auto':
                print(f'  WARNING: empirical ePSF wings poor (pedestal={ped_frac:.1e}, '
                      f'scatter={scat_frac:.1e} > {params["pedestal_bad"]:.0e}); '
                      f'falling back to the model')
                use_model = True
            else:
                if poor:
                    print(f'  WARNING: keeping forced-empirical ePSF despite poor wings '
                          f'(pedestal={ped_frac:.1e}, scatter={scat_frac:.1e})')
                epsf_data = cand
                star_stamps = [np.asarray(s.data, dtype=float) for s in stars]
                method_used = 'empirical'
        elif method == 'empirical':
            raise SystemExit(
                f'ERROR: the empirical ePSF for {a.lens} {a.filt} diverged ({n_stars} '
                f'stars). Add cleaner include stars to info/psf_stars.json or use '
                f'--method model.')
        else:
            print(f'  WARNING: empirical ePSF diverged with {n_stars} stars; '
                  f'falling back to the STDPSF model')
            use_model = True

    if method_used is None:
        tag = 'forced' if method == 'model' else \
              (f'{n_stars} stars < {params["min_stars"]}' if too_few else 'empirical failed')
        # Detector CD for rotating the (detector-frame) model PSF into the North-up drizzle
        # frame. ACS focus-diverse reads per-exposure CDs itself, so this is only used by the
        # WFPC2 DB and STDPSF paths below.
        cd_det = representative_input_cd(inst_key, a.lens, a.filt, a.sample)
        if cd_det is None and inst_key != 'ACS/WFC':
            print('  WARNING: no input CD found; model PSF will not be rotated to North-up')
        # ACS/WFC: prefer the focus-diverse, observation-matched ePSF (native F555W/F814W,
        # focus-corrected, at the lens position) over the static STDPSF. STDPSF stays the
        # fallback-of-the-fallback if a retrieval fails. Other instruments use STDPSF.
        if inst_key == 'ACS/WFC':
            try:
                rootnames = lens_products_rootnames(a.sample, a.lens, a.filt)
                calibrated_dir = os.path.join(ws_path, 'data', 'calibrated', a.sample,
                                              a.lens, a.filt)
                print(f'  method: MODEL (ACS focus-diverse ePSF)  [{tag}]')
                epsf_data = psf_models.acs_focus_diverse_psf(
                    a.lens, a.filt, rootnames, calibrated_dir, catalogue_coord,
                    oversample=oversample, size=star_size, out_scale=scale)
                method_used = 'model_acs_fdpsf'
            except Exception as exc:
                print(f'  WARNING: focus-diverse ePSF unavailable ({exc}); '
                      f'falling back to STDPSF')
        # WFPC2 F606W: prefer a native-filter ePSF built from the MAST PSF database
        # (Dauphin ISR WFC3 2021-12) over the STDPSF F555W-substituted proxy. STDPSF
        # stays the fallback-of-the-fallback if the DB query/build fails.
        elif inst_key == 'WFPC2' and a.filt.split('_')[0].upper() == 'F606W':
            try:
                print(f'  method: MODEL (WFPC2 F606W MAST PSF DB)  [{tag}]')
                epsf_data = psf_models.wfpc2_f606w_db_epsf(
                    oversample=oversample, size=star_size, out_scale=scale,
                    cd_detector=cd_det)
                method_used = 'model_wfpc2_psfdb'
            except Exception as exc:
                print(f'  WARNING: WFPC2 F606W PSF-DB unavailable ({exc}); '
                      f'falling back to STDPSF')
        if method_used is None:
            print(f'  method: MODEL (STDPSF)  [{tag}]')
            epsf_data = psf_models.model_psf(inst_key, a.filt, oversample=oversample,
                                             size=star_size, out_scale=scale,
                                             cd_detector=cd_det)
            method_used = 'model'

    # Full kernel (whole ePSF footprint binned to image scale) is the archival product in
    # data/psf/. The trimmed, amplitude-sized kernel for modelling goes next to the science
    # stamp in data/cutouts/ (see trim_kernel_to_amplitude / the --trim-threshold flag).
    kernel = oversampled_to_kernel(epsf_data, oversample, star_size)
    kernel, pedestal_frac = subtract_pedestal(kernel)
    if abs(pedestal_frac) >= 1e-4:
        print(f'  pedestal removed: {pedestal_frac:.2e} of peak (flat ePSF-wing floor)')
    trimmed, trim_size = trim_kernel_to_amplitude(kernel, a.trim_threshold)
    fwhm_pix = measure_fwhm(kernel)
    if fwhm_pix is not None:
        print(f'  kernel FWHM: {fwhm_pix:.2f} px = {fwhm_pix * scale:.3f}"')
    print(f'  full kernel {kernel.shape[0]}px -> trimmed {trim_size}px '
          f'(amplitude<{a.trim_threshold:g} of peak)')

    # ── Write the full kernel + ePSF + QA plot to data/psf/ ─────────────────────
    khdr = fits.Header()
    khdr['PSFMETH'] = (method_used, 'empirical ePSF or STDPSF model')
    khdr['PSFNSTAR'] = (n_stars, 'stars used (0 for model)')
    khdr['PSFOVSMP'] = (oversample, 'ePSF build oversampling factor')
    khdr['PSFPXSCL'] = (round(scale, 5), 'kernel pixel scale (arcsec)')
    khdr['PSFFWHM'] = (round(fwhm_pix, 4) if fwhm_pix is not None else 0.0,
                       'fitted kernel FWHM (pixels)')
    khdr['PSFLENS'] = (a.lens, 'lens')
    khdr['PSFFILT'] = (a.filt, 'filter')
    khdr['PSFKIND'] = ('full', 'full ePSF footprint (trimmed copy in data/cutouts/)')
    khdr['PSFPED'] = (round(pedestal_frac, 6), 'ePSF-wing pedestal removed (fraction of peak)')
    write_fits(kernel, khdr, os.path.join(output_dir, 'psf_kernel.fits'))

    plot_psf(kernel, epsf_data, star_stamps, scale, fwhm_pix,
             os.path.join(output_dir, 'psf.png'),
             title=f'{a.lens}  {a.filt}  [{method_used}]  ({n_stars} stars)')

    # ── Write the trimmed modelling kernel to data/cutouts/ (pass-matched prefix) ─
    cutouts_dir = os.path.join(ws_path, 'data', 'cutouts', a.sample, a.lens, a.filt)
    os.makedirs(cutouts_dir, exist_ok=True)
    prefix = 'cutout_cr' if drizzle_pass == 'cr' else 'cutout'
    thdr = khdr.copy()
    thdr['PSFKIND'] = ('trimmed', 'amplitude-trimmed modelling kernel')
    thdr['PSFTRIM'] = (a.trim_threshold, 'azimuthal PSF < this x peak sets the radius')
    thdr['PSFPASS'] = (drizzle_pass, 'drizzle pass the PSF matches')
    if star_pass != drizzle_pass:
        thdr['PSFSTARP'] = (star_pass, 'drizzle pass the empirical ePSF was built from')
    write_fits(trimmed, thdr, os.path.join(cutouts_dir, f'{prefix}_psf.fits'))

    info_json.update(psf_json, a.sample, a.lens, a.filt, {
        'method': method_used,
        'n_stars': n_stars,
        'fwhm_pix': round(fwhm_pix, 4) if fwhm_pix is not None else None,
        'oversample': oversample,
        'kernel_size': int(kernel.shape[0]),
        'cutout_kernel_size': int(trim_size),
        'trim_threshold': a.trim_threshold,
        'pedestal_frac': round(pedestal_frac, 6),
    })
    print('  done')

    # Model tier lacks drizzle broadening (see make_psf_inject.py docstring) -- auto-chain
    # the injection + promotion so this single command already leaves the drizzle-broadened
    # kernel as the canonical psf_kernel.fits / cutout_[cr_]psf.fits. Empirical builds are
    # already the drizzled PSF (cut from the mosaic), so this is skipped for them. A failure
    # here (e.g. data/drizzle_files/ was cleared) is a graceful degradation: the analytic
    # model this function just wrote stays canonical.
    if method_used.startswith('model'):
        try:
            import make_psf_inject
            print(f'  model-tier build ({method_used}) -> running drizzle-broadened '
                  f'injection + promoting to canonical ...')
            make_psf_inject.run_injection(a.lens, a.filt, a.sample,
                                          drizzle_pass=drizzle_pass,
                                          trim_threshold=a.trim_threshold, promote=True)
        except Exception as exc:
            print(f'  WARNING: injection/promotion failed ({exc}); canonical PSF stays '
                  f'the analytic model ({method_used})')


if __name__ == '__main__':
    main()
