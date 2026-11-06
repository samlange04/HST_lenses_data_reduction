#!/usr/bin/env python
"""
Build a per-lens, per-filter PSF from the drizzled mosaics for strong-lens modelling.

The canonical PSF (--method auto, the default, and model) is the model PSF injected into
every science frame at the lens and drizzled exactly as the science was (make_psf_inject.py:
the mean of a grid of injected stars, PSFMETH 'inject_*'), wherever a model exists: ACS/WFC
the focus-diverse ePSF (else STDPSF), WFC3/IR and the gallery UVIS bands STDPSF, WFPC2 the
MAST PSF-database ePSF. Since 2026-11-03 (AGENTS.md *The 2026-11-03 PSF audit*): a kernel
built from bright field stars is not the PSF a faint extended source sees in a drizzled
image (their cores' ERR weights vary frame to frame with sub-pixel phase), and the injected
grid mean is the one kernel that recovers an injected Sersic galaxy unbiased.

The empirical ePSF from the field stars is built beside it as the CHECK
(psf_kernel_empirical.fits, `empirical_check` in the record), and SHIPS only with
--method empirical, where the band has no model (WFC3/UVIS F475X/F600LP), or if the
injection fails: isolated field stars in the full drizzled science frame, postage stamps,
an oversampled effective PSF (ePSF) with photutils (DAOStarFinder -> extract_stars ->
EPSFBuilder), point-sampled to the image scale (its samples are already pixel-integrated),
with the measured sky residual removed. This is a production port of an early exploratory
notebook (not tracked), with every manual step of it replaced by automatic quality cuts plus
a per-lens override file:

  * the hand-typed NaN rectangles that blanked the lens galaxy / bad regions  ->  an
    automatic mask (deflector circle from info/coords/lens_coords.py, zero-weight pixels from
    the drizzle WHT map, and a frame-edge margin);
  * the frame-max detection threshold ->  a background-sigma threshold (threshold_scale x
    the sigma-clipped background std);
  * the hand-deletion of stars 5..48 ->  isolation / sharpness / roundness / bounds cuts,
    keep the brightest N, plus an optional info/overrides/psf_stars.json entry (include / exclude
    coords, or per-lens parameter overrides) for the few fields the automatics get wrong.

Model route (fallback): when a field has too few usable stars (< --min-stars), a
diverged ePSF or wings failing the pedestal/scatter gate, or with --method model, build
a model PSF at the deflector position via scripts/psf_models.py (ACS: focus-diverse ePSF,
else STDPSF; WFPC2 F606W: the MAST PSF-database ePSF, else STDPSF; WFC3/IR and the
gallery UVIS bands: STDPSF) and auto-chain make_psf_inject.py so the canonical kernel is
drizzle-broadened. SL2S F475X/F600LP have no model grid; make_psf_median.py covers them.

Outputs (data/psf/<sample>/<band_mode>/<category>/<lens>/<filt>/):
  * psf_kernel.fits  -- image-scale, odd, unit-sum: a drop-in al.Kernel2D at the band
                        pixel scale (an ePSF is point-sampled from the in-memory oversampled
                        build, which is not itself written to disk -- nothing reads it back).
  * psf_kernel_empirical.fits -- on a model-tier product, the field-star ePSF (the check).
  * psf_kernel_phase_std.fits -- on an injected product, the sub-pixel-phase scatter.
  * psf.png          -- QA panel: selected-star montage, the kernel, and a radial profile.

The trimmed, modelling-ready kernel (psf.fits, plus its _err/_analytic
siblings) goes to stamp_paths.psf_stamp_dir(): the band's stamp dir, beside the
sci/noise/mask a fit consumes. One kernel per band, and NOT keyed on stamp size -- it is
trimmed by amplitude, a property of the band, so a --size variant stamp pairs with this
same kernel in the default-size tree.

For a MODEL-tier build (method_used starts with 'model'), this auto-chains into
make_psf_inject.run_injection(..., promote=True): the drizzle-broadened injected kernel
becomes the canonical psf_kernel.fits / psf.fits, and the analytic model this
function just wrote is moved aside to psf_kernel_analytic.fits / psf_analytic.fits
(comparison only). See make_psf_inject.py's docstring. If injection fails (e.g. the
data/drizzle_files/ inputs were cleared), it falls back to
make_psf_inject.analytic_broadened_fallback() -- a cheap analytic drop-box broadening of the
analytic kernel, promoted the same way -- and only if THAT also fails does the un-broadened
analytic model stay canonical.

Records the build in info/provenance/lens_psf.json (method, n_stars, fwhm_pix, oversample,
kernel_size); no data for a filter records null and exits 0, matching the other scripts.

Usage:
    uv run python scripts/make_psf.py --lens J0252+0039 --filt f814W
"""

import argparse
import glob
import json
import os
import shutil
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

ws_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# HST_INFO_DIR: per-shard private info dir for parallel kernel rebuilds
# (concurrent info_json.update calls tear shared JSONs); unset = the checkout's info/.
_INFO_ROOT = os.environ.get('HST_INFO_DIR') or os.path.join(ws_path, 'info')
sys.path.insert(0, os.path.join(_INFO_ROOT, 'coords'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mast_target_names
# Reuse the canonical product locator and the coordinate resolver rather than
# re-implementing the instrument-varying prefix/suffix globbing or the per-sample coords
# lookup (info/coords/lens_coords.py, every sample's table merged). Importing
# make_stamps is safe: its work is under __main__.
from make_stamps import find_products, _has_products, catalogue_coord_for
import psf_models
import info_json
import stamp_paths


# ── Per-instrument defaults ─────────────────────────────────────────────────────
# Keyed on a normalised instrument label (see instrument_key). WFPC2's drizzled primary
# header records DETECTOR=PC (the aperture, not the chip -- the lens is on WF3), so WFPC2
# is keyed on INSTRUME alone. Every value here is overridable by info/overrides/psf_stars.json and
# then by an explicit CLI flag; the drizzled pixel scale is read from the WCS, not stored.
_BASE = dict(
    threshold_scale=5.0,     # detection threshold = this x sigma-clipped background std
    min_stars=3,             # fewer usable stars than this -> model fallback
    kernel_size=31,          # odd size of the delivered image-scale kernel (pixels)
    lens_mask_radius=5.0,    # arcsec circle around the deflector, masked out of detection
    maxiters=15,             # EPSFBuilder iterations
    fwhm_tol_lo=0.75,        # per-candidate fitted FWHM must exceed this x psf_fwhm
    fwhm_tol_hi=1.2,         # ...and stay below this x psf_fwhm (rejects galaxies / CRs). With
                             # psf_fwhm the RESAMPLED MODEL's _fit_stamp_shape FWHM (2026-11-03;
                             # was a looser guess x1.4): genuine stars fit 0.99-1.07x their
                             # injected twins, compact galaxies 1.24-1.35x, and the old bounds
                             # (ACS 2.8 px, UVIS 3.36 px against models of 1.92 and 2.11 px) let
                             # 4 of 25 and 6 of 25 selected "stars" in at 1.15-1.35x
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
    crowd_sep_frac=1.0,      # candidate-pair rejection radius, as a fraction of star_size
                             # (the extract_stars() stamp width). DAOStarFinder's own
                             # min_separation (below, ~0.5x star_size) only keeps it from
                             # double-detecting one blended peak; it does not know the
                             # downstream stamp size, so two accepted detections can still
                             # sit close enough for their star_size stamps to overlap and
                             # contaminate each other's wings. 1.0 = stamps never overlap.
)
# `fwhm` is the DAOStarFinder detection kernel; `psf_fwhm` is the true stellar PSF FWHM
# (pixels, at the drizzled scale) that the per-candidate shape cut selects around --
# distinct because detection tolerates a coarser kernel than the star/galaxy separation.
# oversample=2: SLACS fields yield only ~10-20 usable stars, so a 4x-oversampled ePSF
# has too many pixels per constraining star and comes out noisy. 2x is well-matched to
# these (near-)critically-sampled drizzled grids and much more robust; raise per lens via
# --oversample or info/overrides/psf_stars.json where a field is unusually star-rich.
# wing_radius (output px): where the HYBRID kernel switches from the injected model (core)
# to the field-star ePSF (wing) -- see hybrid_kernel. Set from the closure tests of
# 2026-11-03: the ERR-weight artefact that deforms a star-built ePSF's core has deposited all
# its displaced flux inside 13 px on ACS/UVIS (closure EE reaches 1.000 by r = 14) and ~13 px
# on F160W, and the model ePSFs hold < 0.1% of their flux beyond 12 px (ACS/UVIS), so beyond
# these radii the ePSF carries only real light and the model none. None: no hybrid (WFPC2
# fields have no usable stars).
_INSTR = {
    'ACS/WFC': dict(fwhm=2.5, psf_fwhm=1.95, oversample=2, star_size=51, max_stars=25, wing_radius=12,
                    min_sep=25.0, sharplo=0.3, sharphi=1.0, roundlo=-0.7, roundhi=0.7),
    'WFC3/IR': dict(fwhm=2.8, psf_fwhm=3.55, oversample=2, star_size=41, max_stars=25, wing_radius=13,
                    # 3.55 px is the real F160W stars' gate FWHM (3.46 for the uncalibrated
                    # STDPSF twins); the x1.2 window then spans 2.7-4.3 px
                    min_sep=20.0, sharplo=0.2, sharphi=1.0, roundlo=-0.7, roundhi=0.7),
    # star_size=35 (smaller than ACS): WFPC2 fields are star-poor and EPSFBuilder diverges
    # on a 51px stamp with only ~3 stars but converges on 35px (verified J0252+0039 F606W).
    'WFPC2':   dict(fwhm=2.5, psf_fwhm=4.5, oversample=2, star_size=35, max_stars=25, wing_radius=None,
                    min_sep=25.0, sharplo=0.3, sharphi=1.0, roundlo=-0.7, roundhi=0.7),
    # WFC3/UVIS (BELLS GALLERY): a two-CCD optical detector like ACS/WFC, drizzled at the
    # native 0.0396"/px. psf_fwhm is the resampled STDPSF's gate FWHM (2.11 px on the F475X
    # twins, 2026-11-03; the [0.75,1.2]x shape window brackets the 2.0-2.3px star locus) and
    # the DAO detection kernel 2.8. The wide UVIS field (~162") is star-rich, so
    # the empirical build is the norm. flux_floor_frac is lowered to 0.02 (from the 0.05
    # base): these fields routinely contain one very bright (near-saturated) star whose flux
    # sets the 5%-floor high enough to discard a dozen genuinely good, high-S/N mid-range
    # stars -- the "brightest star collapses the flux floor" trap. 0.02 keeps them.
    'WFC3/UVIS': dict(fwhm=2.8, psf_fwhm=2.1, oversample=2, star_size=51, max_stars=25, wing_radius=12,
                      min_sep=25.0, sharplo=0.3, sharphi=1.0, roundlo=-0.7, roundhi=0.7,
                      flux_floor_frac=0.02),
}


# Detector full well in the mosaic's native per-frame units (ACS/UVIS FLC: electrons; WFPC2
# c0m: DN, 12-bit ADC). A candidate whose drizzled peak exceeds SAT_PEAK_FRAC x full well /
# the LONGEST frame exposure is saturated in at least one frame (DQ 256 is drizzled in as
# "good", so nothing upstream removes it) and must not enter the ePSF. 0.5 allows for the
# drizzle-lowered peak and the flat top: a star just at full well in the longest frame lands
# at ~0.7-0.9 of sat/dexp in the mosaic. Losing an unsaturated star in the top half of the
# dynamic range is harmless (the ePSF is flux-normalised); admitting a saturated one is not.
# WFC3/IR: up-the-ramp fitting handles saturation, no gate. Found 2026-09-29 on the pixfrac
# 1.0 gallery rebuild: J0742+3341 f814W's 195 e-/s saturated star had failed the shape gate
# by accident at pixfrac 0.7, passed it at 1.0, set the 2% flux floor (8 -> 4 stars) and
# went into the ePSF. Overridable via info/overrides/psf_stars.json `satlevel` / `max_peak`.
_SATLEVEL = {'ACS/WFC': 84700.0, 'WFC3/UVIS': 70000.0, 'WFPC2': 4095.0, 'WFC3/IR': None}
SAT_PEAK_FRAC = 0.5


def saturation_peak(inst_key, sci_hdr, overrides):
    """Mosaic peak (per-second units) above which a candidate is saturated, or None."""
    if 'max_peak' in overrides:
        return overrides['max_peak']
    sat = overrides.get('satlevel', _SATLEVEL.get(inst_key))
    dexp = [v for k, v in sci_hdr.items() if k.startswith('D') and k.endswith('DEXP')]
    if not sat or not dexp:
        return None
    return SAT_PEAK_FRAC * float(sat) / max(float(v) for v in dexp)


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


def reject_crowded(rows, min_separation):
    """Greedy-by-flux crowding cut: drop any candidate whose extract_stars() stamp would
    overlap an already-kept (brighter) one's stamp.

    `rows` must already be flux-sorted, brightest first (as select_stars sorts them), so
    the brightest of any close pair is always the one kept -- a fainter neighbour is
    dropped in favour of it rather than either being silently mixed into the ePSF. Distance
    is Euclidean, `min_separation` is normally `star_size` (see crowd_sep_frac): a Euclidean
    cut at the full stamp width is a conservative (not exact-square) non-overlap test, which
    is the right side to err on for a "don't contaminate the wings" gate.
    """
    kept = []
    for x, y, flux in rows:
        if all((x - kx) ** 2 + (y - ky) ** 2 >= min_separation ** 2 for kx, ky, _ in kept):
            kept.append((x, y, flux))
    return kept


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

    `overrides` (from info/overrides/psf_stars.json) is applied last: `exclude` points/boxes remove
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

    qa = {'n_detected': 0 if sources is None else len(sources), 'n_crowded_rejected': 0,
          'n_saturated_rejected': 0}
    max_peak = params.get('max_peak')
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
            # saturated in at least one frame (see _SATLEVEL): never a PSF star
            if max_peak is not None and peak > max_peak:
                qa['n_saturated_rejected'] += 1
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
        n_before_crowd = len(rows)
        crowd_sep = params['crowd_sep_frac'] * star_size
        rows = reject_crowded(rows, crowd_sep)
        qa['n_crowded_rejected'] = n_before_crowd - len(rows)
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


def oversampled_to_kernel(epsf_data, oversample, kernel_size, integrated=True):
    """Reduce an oversampled PSF to an image-scale, unit-sum, odd kernel centred on its core.

    `integrated=True` (the default; a photutils ePSF, or a model ePSF resampled by
    psf_models._resample_centered): the samples are ALREADY pixel-integrated -- an ePSF's
    value at an offset is the flux a pixel centred there collects (Anderson & King 2000) --
    so the image-scale kernel is the ePSF POINT-SAMPLED at integer-pixel offsets, every
    `oversample`-th sample through the core sample. Block-summing it (this function until
    2026-11-03, and the JWST repo's until 2026-10-31) integrates over the pixel a second
    time: on a noiseless synthetic closure test with this repo's own build it lowered the
    ACS/WFC peak by 12.5%, broadened the radial FWHM by 15% and, because a block of an even
    number of samples has no centre sample, shifted every kernel by +0.25 px in both axes
    (the "0.26 px median kernel offset" once recorded as a measurement floor). Point
    sampling reproduces the truth to < 1% in peak, width and encircled energy, and to
    < 0.01 px (AGENTS.md *PSF generation*; the JWST repo found the same, its AGENTS.md
    *PSF*, finding 1).

    `integrated=False` (a raw optical PSF that has not been integrated over a pixel):
    sum each oversample x oversample block. No model this repo uses is of that kind (STDPSF,
    the focus-diverse ePSF and the MAST PSF-database ePSF are all Anderson ePSFs).

    The core is the centroid of a window around the peak, falling back to the peak, then
    to the geometric centre, if a step is degenerate.
    """
    arr = np.asarray(epsf_data, dtype=float)
    ny, nx = arr.shape
    # locate the core: centroid of a window around the peak, guarding non-finite results
    py, px = np.unravel_index(np.argmax(arr), arr.shape)
    win = max(oversample * 5, 9)
    y0, x0 = max(0, py - win // 2), max(0, px - win // 2)
    # photutils centroid_com returns (x, y) -- COLUMN first. Unpacking it as (y, x) transposes
    # the sub-pixel offset, so the crop below recentres by (dx, dy) instead of (dy, dx) and
    # leaves the kernel off its own array centre by an antisymmetric residual (dy = -dx).
    # That was the bug behind kernels sitting up to 1.7px off centre; measured repo-wide the
    # residual offsets were anticorrelated at -0.76, which is that transposition's signature.
    cx0, cy0 = centroid_com(arr[y0:y0 + win, x0:x0 + win])
    cx, cy = x0 + cx0, y0 + cy0
    if not (np.isfinite(cx) and np.isfinite(cy)):
        cx, cy = float(px), float(py)

    for centre in ((cx, cy), (float(px), float(py)), (nx / 2.0, ny / 2.0)):
        if integrated:
            half = kernel_size // 2
            ic, jc = int(round(centre[1])), int(round(centre[0]))
            pad = half * oversample + 1
            big = np.pad(arr, pad)
            idx = np.arange(-half, half + 1) * oversample
            kernel = big[np.ix_(ic + pad + idx, jc + pad + idx)]
        else:
            sub = _crop_centered(arr, kernel_size * oversample, *centre)
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


SKY_ANNULUS = (1.5, 2.5)   # x the kernel half-width: where sky_residual_level reads the sky


def sky_residual_level(data, stars_tbl, star_size, annulus=SKY_ANNULUS):
    """The flat background the ePSF inherits from its stars, per pixel, in units of the
    kernel's unit sum: the median over the PSF stars of (sky level in an annulus OUTSIDE the
    star cutout) / (the cutout's sum, which is the flux EPSFBuilder normalises each star by).
    Measured on the background-subtracted mosaic, so it is Background2D's residual under the
    stars, not the sky itself. The annulus is 1.5-2.5 x the kernel half-width (ACS: 1.9-3.1"),
    where the stellar wing is ~1e-5 of the peak (ACS/WFC encircled-energy curve) against the
    1e-4 the ePSF's own outer annulus reads. Returns 0.0 with no usable star.
    """
    half = star_size // 2
    r_in, r_out = annulus[0] * half, annulus[1] * half
    w = int(np.ceil(r_out)) + 1
    ny, nx = data.shape
    yy, xx = np.mgrid[-w:w + 1, -w:w + 1]
    ann = (np.hypot(xx, yy) >= r_in) & (np.hypot(xx, yy) <= r_out)
    levels = []
    for x, y in zip(stars_tbl['x'], stars_tbl['y']):
        ix, iy = int(round(float(x))), int(round(float(y)))
        if not (w <= ix < nx - w and w <= iy < ny - w):
            continue
        big = data[iy - w:iy + w + 1, ix - w:ix + w + 1]
        cut = data[iy - half:iy + half + 1, ix - half:ix + half + 1]
        flux = float(np.nansum(cut))
        if flux <= 0 or not np.isfinite(big[ann]).any():
            continue
        levels.append(float(np.nanmedian(big[ann])) / flux)
    return float(np.median(levels)) if levels else 0.0


def subtract_pedestal(kernel, level=None):
    """Remove a flat background from the kernel and renormalise to unit sum. Returns
    (kernel, pedestal_frac), the pedestal as a fraction of the peak.

    EPSFBuilder leaves a flat floor in the ePSF: the residual background under its stars.
    Left in, a ~1e-4-of-peak floor across a 51 px kernel is several percent of the
    (renormalised) flux as a spurious uniform background, and it stops
    trim_kernel_to_amplitude from ever crossing the amplitude threshold.

    `level` (per pixel, in units of the unit-sum kernel; from sky_residual_level) is the
    floor to remove. Without it the median of the kernel's own outer annulus (0.75-1 x the
    half-width) is removed, as until 2026-11-03 -- but that annulus is not empty sky: on
    ACS/WFC it reads 8-13e-5 of the peak where the sky residual outside the cutout is 2-3e-5
    and the published encircled-energy curve puts the real wing at ~4e-5, so the old rule
    removed a median 4% of the flux on ACS F814W (6% UVIS, 9% F160W), most of it real wing,
    and renormalised the core up by as much (AGENTS.md *PSF generation*; the JWST repo found
    the same, its finding 5). Noiseless model and injected kernels are not passed through
    this at all since 2026-11-03.
    """
    arr = np.asarray(kernel, float)
    ny, nx = arr.shape
    cy, cx = ny // 2, nx // 2
    s0 = arr.sum()
    arr = arr / s0 if s0 > 0 else arr
    peak = arr[cy, cx] if arr[cy, cx] > 0 else float(arr.max())
    if level is None:
        yy, xx = np.mgrid[:ny, :nx]
        rint = np.round(np.hypot(xx - cx, yy - cy)).astype(int)
        half = min(cy, cx)
        ann = (rint >= int(0.75 * half)) & (rint <= half)
        ped = float(np.median(arr[ann])) if ann.any() else 0.0
    else:
        ped = float(level)
    out = arr - ped
    s = out.sum()
    if s > 0:
        out = out / s
    return out, (float(ped / peak) if peak > 0 else 0.0)


# ── Empirical PSF uncertainty (bootstrap / jackknife over the star sample) ────────
# The empirical ePSF is built from a handful of field stars (often 3-20, sometimes as few as
# 3), which makes the kernel a point estimate with real, un-propagated uncertainty. We
# estimate a per-pixel PSF error map by resampling the star sample and rebuilding the ePSF:
# bootstrap-with-replacement where there are enough stars, leave-one-out jackknife where there
# are too few for bootstrap to be non-degenerate. Empirical tier only -- the model/injected
# tiers have their own (deferred) error budget. See the PSF uncertainty note in AGENTS.md.
JACKKNIFE_MAX_STARS = 6   # < this, bootstrap draws collapse onto 1-2 stars; use jackknife


def _kernel_from_epsf(epsf_data, oversample, star_size, level=None, pedestal=True):
    """One ensemble member: oversampled ePSF -> the same full image-scale, unit-sum kernel
    the primary build produces (oversampled_to_kernel, then subtract_pedestal(level) for an
    empirical member; `pedestal=False` for a noiseless model member), so members are directly
    comparable to psf_kernel.fits."""
    kernel = oversampled_to_kernel(epsf_data, oversample, star_size)
    if pedestal:
        kernel, _ = subtract_pedestal(kernel, level)
    return kernel


def build_epsf_ensemble(data, stars_tbl, star_size, oversample, maxiters, n_boot, rng=None,
                        level=None):
    """Resample the star sample and rebuild the ePSF per draw for a PSF-uncertainty estimate.

    len(stars_tbl) >= JACKKNIFE_MAX_STARS: bootstrap-with-replacement, `n_boot` draws.
    Fewer: leave-one-out jackknife (n draws), since bootstrap draws degenerate at tiny N.
    Each draw runs the exact primary build (build_epsf) and is validate_epsf-gated; a
    diverged/failed draw is discarded. Returns (kernels, method, n_valid) where `kernels` is a
    list of full image-scale unit-sum kernels (the same post-processing the primary gets).
    """
    rng = rng if rng is not None else np.random.default_rng(0)
    n = len(stars_tbl)
    if n < JACKKNIFE_MAX_STARS:
        method = 'jackknife'
        index_sets = [[j for j in range(n) if j != i] for i in range(n)]
    else:
        method = 'bootstrap'
        index_sets = [rng.integers(0, n, size=n).tolist() for _ in range(n_boot)]

    kernels = []
    for idx in index_sets:
        if not idx:
            continue
        try:
            epsf, _, _ = build_epsf(data, stars_tbl[idx], star_size, oversample, maxiters)
            cand = np.asarray(epsf.data, dtype=float)
            if validate_epsf(cand, oversample):
                kernels.append(_kernel_from_epsf(cand, oversample, star_size, level))
        except Exception:
            continue                       # a degenerate resample is dropped, not fatal
    return kernels, method, len(kernels)


def _fwhm_spread(kernels, method):
    """Standard error of the per-member kernel FWHM, matching the ensemble estimator
    (jackknife leave-one-out vs bootstrap replicate std). None if < 2 finite FWHMs."""
    fw = np.array([f for f in (measure_fwhm(k) for k in kernels) if f is not None],
                  dtype=float)
    if len(fw) < 2:
        return None
    if method == 'jackknife':
        return float(np.sqrt((len(fw) - 1) / len(fw) * np.sum((fw - fw.mean()) ** 2)))
    return float(fw.std(ddof=1))


def psf_error_map(kernels, method, trim_half):
    """Per-pixel PSF uncertainty from an ensemble of unit-sum kernels.

    Members are co-registered (each re-centred by oversampled_to_kernel) and unit-sum, so the
    spread is already in the kernel's own amplitude units -- NOT renormalised. Bootstrap: the
    standard error is the replicate std, sqrt(SS/(B-1)). Jackknife: the leave-one-out standard
    error sqrt((n-1)/n . SS). Returns (err_full, err_trim, err_frac, fwhm_err): err_trim is
    err_full cropped to the primary kernel's trim window (same `trim_half`, not a fresh
    amplitude cut) so it matches psf.fits pixel-for-pixel; err_frac =
    sqrt(sum(err_full**2)) (integrated relative PSF uncertainty; the kernels sum to 1).
    """
    stack = np.stack(kernels, axis=0)
    m = stack.shape[0]
    ss = np.sum((stack - stack.mean(axis=0)) ** 2, axis=0)
    if method == 'jackknife':
        err_full = np.sqrt((m - 1) / m * ss)
    else:
        err_full = np.sqrt(ss / (m - 1))
    ny, nx = err_full.shape
    cy, cx = ny // 2, nx // 2
    err_trim = err_full[cy - trim_half:cy + trim_half + 1, cx - trim_half:cx + trim_half + 1]
    err_frac = float(np.sqrt(np.sum(err_full ** 2)))
    return err_full, err_trim, err_frac, _fwhm_spread(kernels, method)


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
        dd = os.path.join(ws_path, 'data', 'drizzle_files', stamp_paths.lens_rel(sample, lens), filt)
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
    cd_dir = os.path.join(ws_path, 'data', 'calibrated', stamp_paths.lens_rel(sample, lens), filt)
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
    data = info_json.load(os.path.join(_INFO_ROOT, 'provenance', 'lens_products.json'))
    val = data.get(sample, {}).get(lens, {}).get(filt_key)
    return list(val) if val else []


# ── Main ────────────────────────────────────────────────────────────────────────
HYBRID_TAPER = 2.0   # px: the model->ePSF switch is a linear ramp over wing_radius +- this


def hybrid_kernel(model, emp, wing_radius, taper=HYBRID_TAPER):
    """Model core, empirical wing (Sam's call, 2026-11-03). Inside `wing_radius` the injected
    model (the faint-source PSF at the lens: unbiased core, but no wing beyond ~0.6" and a
    6-48% wing deficit); outside it the field-star ePSF (real wing, but a core deformed by the
    stars' ERR weights, which the closure tests show is confined inside wing_radius). The
    weight w(r) ramps linearly from 1 to 0 over wing_radius +- taper. The model is unit-
    normalised over its own weighted footprint and scaled by (1 - f_wing), f_wing being the
    ePSF's flux in the wing, so the kernel sums to 1 with the ePSF's own wing fraction.
    Returns (kernel, w, f_wing) on the larger of the two grids (the smaller zero-padded)."""
    n = max(model.shape[0], emp.shape[0])
    m = _crop_centered(np.asarray(model, float), n, model.shape[1] // 2, model.shape[0] // 2)
    e = _crop_centered(np.asarray(emp, float), n, emp.shape[1] // 2, emp.shape[0] // 2)
    m, e = m / m.sum(), e / e.sum()
    yy, xx = np.mgrid[:n, :n]
    r = np.hypot(xx - n // 2, yy - n // 2)
    w = np.clip((wing_radius + taper - r) / (2 * taper), 0.0, 1.0)
    f_wing = float(((1 - w) * e).sum())
    core = w * m
    k = core / core.sum() * (1 - f_wing) + (1 - w) * e
    return k / k.sum(), w, f_wing


EE_RADII = (1, 2, 3, 5)


def compare_kernels(kernel, ref):
    """The empirical check on a model kernel: FWHM (measure_fwhm, px), encircled energy within
    EE_RADII px of both, and the core L1 sum |k - ref| within 5 px, on unit-sum kernels cropped
    to the smaller one (ported from the JWST repo, 2026-11-03)."""
    n = min(kernel.shape[0], ref.shape[0])
    k = _crop_centered(kernel, n, kernel.shape[1] // 2, kernel.shape[0] // 2)
    r = _crop_centered(ref, n, ref.shape[1] // 2, ref.shape[0] // 2)
    k, r = k / k.sum(), r / r.sum()
    yy, xx = np.mgrid[:n, :n]
    rad = np.hypot(xx - n // 2, yy - n // 2)
    fw, fr = measure_fwhm(k), measure_fwhm(r)
    return dict(fwhm_pix=round(fw, 4) if fw else None, fwhm_pix_model=round(fr, 4) if fr else None,
                ee=[round(float(k[rad <= R].sum()), 4) for R in EE_RADII],
                ee_model=[round(float(r[rad <= R].sum()), 4) for R in EE_RADII],
                ee_radii_px=list(EE_RADII), l1_core=round(float(np.abs(k - r)[rad <= 5].sum()), 4),
                peak_ratio=round(float(k.max() / r.max()), 4))


def prepare_field(sample, lens, filt, cli=None, drizzle_pass='auto', quiet=False):
    """The band as the empirical build sees it: the mosaic of the star pass, background-
    subtracted, masked, with the PSF stars selected -- everything main() needs before a
    kernel is built, and what psf_star_check.py needs to put a twin at every star. Returns a
    dict (sci_hdr, data, wht, wcs, inst_key, scale, params, star_size, oversample,
    catalogue_coord, stars_tbl, qa, star_pass, drizzle_pass, sci_file, overrides, method,
    drizzled_dir), or None when the band has no drizzled products."""
    class _NoCli:
        def __getattr__(self, k):
            return None
    a = cli if cli is not None else _NoCli()
    drizzled_dir = os.path.join(ws_path, 'data', 'drizzled', stamp_paths.lens_rel(sample, lens), filt)
    if not os.path.isdir(drizzled_dir) or not _has_products(drizzled_dir):
        return None
    _print = (lambda *x, **k: None) if quiet else print
    # Per-lens overrides (loaded before the pass resolution: an override can redirect
    # which drizzle pass the empirical build reads).
    overrides_all = info_json.load(os.path.join(_INFO_ROOT, 'overrides', 'psf_stars.json'))
    overrides = overrides_all.get(sample, {}).get(lens, {}).get(filt, {})

    # Output pass -- the science stamp's pass, recorded in the kernel header (PSFPASS).
    # auto -> cr where a CR pass exists. The kernel file is psf.fits either way.
    if drizzle_pass == 'auto':
        has_cr = bool(glob.glob(os.path.join(drizzled_dir, '*_cr_*_sci.fits')))
        drizzle_pass = 'cr' if has_cr else 'nocrrej'

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
        _print(f'{lens} {filt}  [{drizzle_pass} pass]')
    else:
        _print(f'{lens} {filt}  [build from {star_pass} pass, name as {drizzle_pass}]')
    _print(f'  sci: {os.path.basename(sci_file)}')

    with fits.open(sci_file) as hdul:
        sci_hdr = hdul[0].header
        sci_data = np.array(hdul[0].data, dtype=float)
    with fits.open(wht_file) as hdul:
        wht_data = np.array(hdul[0].data, dtype=float)
    wcs = WCS(sci_hdr)

    inst_key = instrument_key(sci_hdr)
    scale = proj_plane_pixel_scales(wcs.celestial)[0] * 3600.0
    _print(f'  instrument: {inst_key}   pixel scale: {scale:.4f}"/pix')

    # Resolved parameters.
    method = overrides.get('method', 'auto') if (getattr(a, 'method', None) in (None, 'auto')) else a.method
    params = resolve_params(inst_key, overrides, a)
    params['max_peak'] = saturation_peak(inst_key, sci_hdr, overrides)
    if params['max_peak'] is not None:
        _print(f'  saturation gate: mosaic peak > {params["max_peak"]:.1f} rejected')
    star_size = int(params['star_size'])
    oversample = int(params['oversample'])
    kernel_size = int(params['kernel_size'])
    if kernel_size % 2 == 0:
        kernel_size += 1

    ra, dec = catalogue_coord_for(lens)   # any sample, via info/coords/lens_coords.py
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
        _print(f'  WARNING: Background2D failed ({exc}); using median subtraction')
        _, med, _ = sigma_clipped_stats(sci_data, mask=mask, sigma=3.0)
        data = sci_data - med

    stars_tbl, qa = select_stars(data, mask, params, overrides, star_size)
    return dict(sci_hdr=sci_hdr, data=data, wht=wht_data, wcs=wcs, inst_key=inst_key, scale=scale,
                params=params, star_size=star_size, oversample=oversample,
                catalogue_coord=catalogue_coord, stars_tbl=stars_tbl, qa=qa, star_pass=star_pass,
                drizzle_pass=drizzle_pass, sci_file=sci_file, overrides=overrides, method=method,
                drizzled_dir=drizzled_dir, mask=mask)



def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--lens', default='J0252+0039')
    p.add_argument('--filt', default='f814W')
    p.add_argument('--sample', default=mast_target_names.DEFAULT_SAMPLE,
                   help='sample subdirectory of data/drizzled/ and data/psf/. Defined in '
                        f'info/samples/lens_samples.json (default {mast_target_names.DEFAULT_SAMPLE})')
    p.add_argument('--pass', dest='drizzle_pass',
                   choices=['auto', 'cr', 'nocrrej'], default='auto',
                   help="drizzle pass the kernel is NAMED for (psf.fits vs "
                        "psf.fits). 'auto' (default) matches the science stamp. Stars "
                        "are detected in the no-CR mosaic whenever one exists (LACosmic eats "
                        "star cores), else in this pass; override with psf_star_pass in "
                        "info/overrides/psf_stars.json.")
    p.add_argument('--method', choices=['auto', 'empirical', 'model'], default='auto',
                   help="'auto' (default) builds an empirical ePSF and falls back to the "
                        "model tier (focus-diverse / PSF-DB / STDPSF, then injected) when "
                        "usable stars < --min-stars, the ePSF diverges or its wings fail the "
                        "quality gate; 'empirical' forces the star build (fails if too few); "
                        "'model' forces the model tier.")
    p.add_argument('--trim-threshold', dest='trim_threshold', type=float, default=1e-3,
                   help="amplitude fraction of peak at which the trimmed modelling kernel "
                        "(<psf stamp dir>/<...>_psf.fits) is cut. Default 1e-3 -- an "
                        "amplitude criterion, not enclosed-energy (which under-sizes).")
    p.add_argument('--n-boot', dest='n_boot', type=int, default=None,
                   help='bootstrap resamples for the empirical PSF error map (default 100; '
                        'leave-one-out jackknife is used instead when < %d stars). Empirical '
                        'tier only.' % JACKKNIFE_MAX_STARS)
    p.add_argument('--no-psf-err', dest='no_psf_err', action='store_true', default=False,
                   help='skip the PSF error map (psf_kernel_err.fits and '
                        'psf_err.fits) for any tier -- empirical (star '
                        'resample), ACS focus-diverse (per-exposure), or WFPC2 MAST-DB '
                        '(bootstrap over DB stars). STDPSF has none regardless. The '
                        'point-estimate kernel is unchanged either way.')
    # Tunables: None here means 'unset' so instrument defaults / JSON overrides win; a
    # value on the command line takes precedence over both. See resolve_params().
    p.add_argument('--threshold-scale', dest='threshold_scale', type=float, default=None)
    p.add_argument('--fwhm', type=float, default=None)
    p.add_argument('--max-stars', dest='max_stars', type=int, default=None)
    p.add_argument('--min-sep', dest='min_sep', type=float, default=None)
    p.add_argument('--crowd-sep-frac', dest='crowd_sep_frac', type=float, default=None,
                   help='candidate-pair rejection radius as a fraction of star_size (the '
                        'extract_stars() stamp width); default 1.0 so two kept stars never '
                        'share stamp footprint')
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
                   help='output dir, default data/psf/<sample>/<band_mode>/<category>/<lens>/<filt>')
    a = p.parse_args()

    drizzled_dir = os.path.join(ws_path, 'data', 'drizzled', stamp_paths.lens_rel(a.sample, a.lens), a.filt)
    output_dir = a.output or os.path.join(ws_path, 'data', 'psf', stamp_paths.lens_rel(a.sample, a.lens), a.filt)
    psf_json = os.path.join(_INFO_ROOT, 'provenance', 'lens_psf.json')

    # No-data outcome: matches the drizzle/stamp scripts -- record null, exit 0.
    if not os.path.isdir(drizzled_dir) or not _has_products(drizzled_dir):
        print(f'=== NO DATA: {a.lens} {a.filt} (no drizzled products in {drizzled_dir})')
        info_json.update(psf_json, a.sample, a.lens, a.filt, None)
        sys.exit(0)

    os.makedirs(output_dir, exist_ok=True)

    fld = prepare_field(a.sample, a.lens, a.filt, cli=a, drizzle_pass=a.drizzle_pass)
    sci_hdr, data, wht_data, wcs = fld['sci_hdr'], fld['data'], fld['wht'], fld['wcs']
    inst_key, scale, params = fld['inst_key'], fld['scale'], fld['params']
    star_size, oversample = fld['star_size'], fld['oversample']
    catalogue_coord, stars_tbl, qa = fld['catalogue_coord'], fld['stars_tbl'], fld['qa']
    star_pass, drizzle_pass, overrides, method = fld['star_pass'], fld['drizzle_pass'], fld['overrides'], fld['method']
    n_stars = len(stars_tbl)
    print(f'  stars: {qa["n_detected"]} detected -> {qa["n_auto"]} auto-kept'
          f' + {qa["n_included"]} forced = {n_stars}'
          f' ({qa["n_crowded_rejected"]} dropped as crowded, '
          f'{qa["n_saturated_rejected"]} saturated)')

    too_few = n_stars < params['min_stars']
    if method == 'empirical' and too_few:
        raise SystemExit(
            f'ERROR: only {n_stars} usable stars for {a.lens} {a.filt} '
            f'(need >= {params["min_stars"]}). Loosen cuts, add include stars to '
            f'info/overrides/psf_stars.json, or use --method model.')

    # ── The empirical ePSF from the field stars ──────────────────────────────────────
    # Shipped with --method empirical, or when the band has no model (WFC3/UVIS F475X and
    # F600LP); otherwise built beside the model as its CHECK (written as
    # psf_kernel_empirical.fits, compared in the record's `empirical_check`). Since 2026-11-03
    # the model is canonical wherever one exists: a kernel built from bright field stars is
    # not the PSF a faint, extended source sees in a drizzled image. Their cores carry ERR
    # weights ~100x below the sky's that differ frame to frame with sub-pixel phase, so the
    # resample takes a star's core from the frames where it is sharpest; the ePSF then
    # registers the phases on top. Against the model injected at the lens on the lens's own
    # frames and weights (make_psf_inject), a point-sampled ePSF's core is 11-21% too peaked
    # (ACS, UVIS, IR closure tests), and a Sersic lens galaxy drizzled through the real frames
    # is recovered by the injected kernel with its true size and index where every star-built
    # kernel biases R_e by 10-35% (AGENTS.md *PSF generation*; the JWST repo's finding 4).
    emp = None
    star_stamps, epsf_data, method_used = [], None, None
    if not too_few and method != 'model':
        print(f'  empirical ePSF: {n_stars} stars, oversample={oversample} star_size={star_size}')
        epsf, stars, fitted = build_epsf(data, stars_tbl, star_size, oversample,
                                         params['maxiters'])
        cand = np.asarray(epsf.data, dtype=float)
        if validate_epsf(cand, oversample):
            k_emp = oversampled_to_kernel(cand, oversample, star_size)
            ped_frac, scat_frac = wing_stats(k_emp)
            poor = ped_frac > params['pedestal_bad'] or scat_frac > params['scatter_bad']
            if poor:
                print(f'  WARNING: empirical ePSF wings poor (pedestal={ped_frac:.1e}, '
                      f'scatter={scat_frac:.1e} > {params["pedestal_bad"]:.0e})')
            sky_level = sky_residual_level(data, stars_tbl, star_size)
            k_emp, ped_emp = subtract_pedestal(k_emp, sky_level)
            emp = dict(kernel=k_emp, epsf=cand, pedestal_frac=ped_emp, poor=poor,
                       sky_level=sky_level, annulus_frac=ped_frac,
                       star_stamps=[np.asarray(s.data, dtype=float) for s in stars])
        else:
            print(f'  WARNING: the empirical ePSF diverged ({n_stars} stars)')
    if method == 'empirical' and (emp is None or emp['poor']):
        if emp is None:
            raise SystemExit(
                f'ERROR: the empirical ePSF for {a.lens} {a.filt} diverged ({n_stars} '
                f'stars). Add cleaner include stars to info/overrides/psf_stars.json or use '
                f'--method model.')
        print('  WARNING: keeping the forced-empirical ePSF despite poor wings')

    # Model-tier PSF error ensemble (ACS focus-diverse: per-exposure; WFPC2 DB: bootstrap
    # over the shared archival stars). None for STDPSF, which has no natural per-lens
    # ensemble. Populated below, consumed by the psf-err block after the kernel is built.
    model_ensemble = model_err_method = err_source = None

    if method != 'empirical':
        tag = 'forced' if method == 'model' else \
              (f'{n_stars} stars < {params["min_stars"]}' if too_few else 'canonical; ePSF is the check')
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
                calibrated_dir = os.path.join(ws_path, 'data', 'calibrated', stamp_paths.lens_rel(a.sample, a.lens), a.filt)
                print(f'  method: MODEL (ACS focus-diverse ePSF)  [{tag}]')
                epsf_data, model_ensemble, model_err_method = psf_models.acs_focus_diverse_psf(
                    a.lens, a.filt, rootnames, calibrated_dir, catalogue_coord,
                    oversample=oversample, size=star_size, out_scale=scale,
                    return_ensemble=True)
                method_used = 'model_acs_fdpsf'
                err_source = 'exposures'
            except Exception as exc:
                print(f'  WARNING: focus-diverse ePSF unavailable ({exc}); '
                      f'falling back to STDPSF')
        # WFPC2 F606W/F814W: prefer a native-filter ePSF built from the MAST PSF database
        # (Dauphin ISR WFC3 2021-12) over STDPSF (F555W-substituted for F606W; exact for
        # F814W but with no error calibration). STDPSF stays the fallback-of-the-fallback
        # if the DB query/build fails.
        elif (inst_key == 'WFPC2'
              and a.filt.split('_')[0].upper() in psf_models.WFPC2_DB_FILTERS):
            db_filt = a.filt.split('_')[0].upper()
            try:
                print(f'  method: MODEL (WFPC2 {db_filt} MAST PSF DB)  [{tag}]')
                epsf_data, model_ensemble, model_err_method = psf_models.wfpc2_f606w_db_epsf(
                    oversample=oversample, size=star_size, out_scale=scale,
                    cd_detector=cd_det, return_ensemble=True, filt=db_filt)
                method_used = 'model_wfpc2_psfdb'
                err_source = 'db_stars'
            except Exception as exc:
                print(f'  WARNING: WFPC2 {db_filt} PSF-DB unavailable ({exc}); '
                      f'falling back to STDPSF')
        if method_used is None:
            try:
                epsf_data = psf_models.model_psf(inst_key, a.filt, oversample=oversample,
                                                 size=star_size, out_scale=scale,
                                                 cd_detector=cd_det)
                method_used = 'model'
                print(f'  method: MODEL (STDPSF)  [{tag}]')
            except Exception as exc:
                print(f'  no model PSF for {inst_key} {a.filt} ({exc}); '
                      f'the empirical ePSF ships')


    if method_used is None:
        # no model for this band (or --method empirical): the ePSF ships
        if emp is None:
            raise SystemExit(f'ERROR: {a.lens} {a.filt} has no model PSF and no usable '
                             f'empirical ePSF ({n_stars} stars).')
        epsf_data, star_stamps, method_used = emp['epsf'], emp['star_stamps'], 'empirical'
        print(f'  method: EMPIRICAL{"" if method == "empirical" else " (no model for this band)"}')

    # Full kernel (whole ePSF footprint binned to image scale) is the archival product in
    # data/psf/. The trimmed, amplitude-sized kernel for modelling goes next to the science
    # stamp (see trim_kernel_to_amplitude / the --trim-threshold flag).
    # Pedestal: only the EMPIRICAL kernel carries one (the residual background under its
    # stars, measured outside their cutouts by sky_residual_level); a model kernel is
    # noiseless and is left alone. See subtract_pedestal for the old annulus rule it replaces.
    sky_level = pedestal_frac = 0.0
    if method_used == 'empirical':
        kernel, pedestal_frac, sky_level = emp['kernel'], emp['pedestal_frac'], emp['sky_level']
        print(f'  pedestal removed: {pedestal_frac:.2e} of peak (sky residual outside the '
              f'star cutouts; the ePSF annulus itself reads {emp["annulus_frac"]:.2e}, mostly wing)')
    else:
        kernel = oversampled_to_kernel(epsf_data, oversample, star_size)
    trimmed, trim_size = trim_kernel_to_amplitude(kernel, a.trim_threshold)
    fwhm_pix = measure_fwhm(kernel)
    if fwhm_pix is not None:
        print(f'  kernel FWHM: {fwhm_pix:.2f} px = {fwhm_pix * scale:.3f}"')
    print(f'  full kernel {kernel.shape[0]}px -> trimmed {trim_size}px '
          f'(amplitude<{a.trim_threshold:g} of peak)')

    # ── PSF uncertainty: per-pixel error map from a tier-appropriate ensemble ───────
    # Empirical: resample the lens's own field stars and rebuild the ePSF (bootstrap, or
    # jackknife when star-poor). Model tier (ACS focus-diverse / WFPC2 MAST-DB): reuse the
    # ensemble psf_models already built (per-exposure or bootstrap-over-DB-stars, see
    # psf_models._reduce_ensemble / _wfpc2_f606w_db_ensemble) via the SAME psf_error_map()
    # statistics, just converted to kernels first. STDPSF has no natural per-lens ensemble
    # and is skipped (err stays null) -- see the PSF uncertainty note in AGENTS.md. The
    # point-estimate kernel above is untouched either way; this only adds parallel
    # *_err.fits products.
    err_full = err_trim = err_frac = fwhm_err = err_method = n_boot_valid = None
    do_psf_err = not a.no_psf_err and overrides.get('psf_err', True)
    trim_half = (trim_size - 1) // 2
    kernels = None
    if do_psf_err and method_used == 'empirical' and n_stars >= 2:
        n_boot = a.n_boot if a.n_boot is not None else int(overrides.get('n_boot', 100))
        mode = 'jackknife' if n_stars < JACKKNIFE_MAX_STARS else f'bootstrap x{n_boot}'
        print(f'  PSF uncertainty: resampling {n_stars} stars ({mode}) ...')
        kernels, err_method, n_boot_valid = build_epsf_ensemble(
            data, stars_tbl, star_size, oversample, params['maxiters'], n_boot,
            level=sky_level)
        err_source = 'stars'
    elif (do_psf_err and method_used in ('model_acs_fdpsf', 'model_wfpc2_psfdb')
          and model_ensemble):
        print(f'  PSF uncertainty: {len(model_ensemble)} {err_source} ensemble members '
              f'({model_err_method}) ...')
        kernels = [_kernel_from_epsf(m, oversample, star_size, pedestal=False)
                   for m in model_ensemble]
        err_method, n_boot_valid = model_err_method, len(kernels)

    if kernels is not None:
        if n_boot_valid >= 2:
            err_full, err_trim, err_frac, fwhm_err = psf_error_map(
                kernels, err_method, trim_half)
            print(f'    {err_method}: {n_boot_valid} valid members, PSF error fraction '
                  f'{err_frac:.2e}'
                  + (f', FWHM err {fwhm_err:.3f}px' if fwhm_err is not None else ''))
        else:
            print(f'    WARNING: only {n_boot_valid} valid ensemble member(s); '
                  f'no PSF error map written')
            err_method = err_source = None

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
    khdr['PSFKIND'] = ('full', 'full ePSF footprint (trimmed copy in the psf stamp dir)')
    khdr['PSFPED'] = (round(pedestal_frac, 6), 'ePSF-wing pedestal removed (fraction of peak)')
    write_fits(kernel, khdr, os.path.join(output_dir, 'psf_kernel.fits'))

    plot_psf(kernel, epsf_data, star_stamps, scale, fwhm_pix,
             os.path.join(output_dir, 'psf.png'),
             title=f'{a.lens}  {a.filt}  [{method_used}]  ({n_stars} stars)')

    # ── Write the trimmed modelling kernel (psf.fits) into the band's stamp dir ─────
    # stamp_paths.psf_stamp_dir pins it to the default-size tree (the kernel is not
    # size-keyed), beside the sci/noise/mask a fit consumes.
    stamps_dir = stamp_paths.psf_stamp_dir(ws_path, a.sample, a.lens, a.filt)
    os.makedirs(stamps_dir, exist_ok=True)
    thdr = khdr.copy()
    thdr['PSFKIND'] = ('trimmed', 'amplitude-trimmed modelling kernel')
    thdr['PSFTRIM'] = (a.trim_threshold, 'azimuthal PSF < this x peak sets the radius')
    thdr['PSFPASS'] = (drizzle_pass, 'drizzle pass the PSF matches')
    if star_pass != drizzle_pass:
        thdr['PSFSTARP'] = (star_pass, 'drizzle pass the empirical ePSF was built from')
    write_fits(trimmed, thdr, os.path.join(stamps_dir, f'psf.fits'))

    # ── Write the PSF error maps (full to data/psf/, trimmed next to the stamp) ─
    if err_full is not None:
        for base, kind, arr, dest in (
                (khdr, 'full_err', err_full, os.path.join(output_dir, 'psf_kernel_err.fits')),
                (thdr, 'trimmed_err', err_trim,
                 os.path.join(stamps_dir, f'psf_err.fits'))):
            ehdr = base.copy()
            ehdr['PSFKIND'] = (kind, 'per-pixel PSF std (ensemble)')
            ehdr['PSFERR'] = (err_method, 'bootstrap or jackknife ensemble method')
            ehdr['PSFESRC'] = (err_source, 'ensemble source: stars/exposures/db_stars')
            ehdr['PSFNVALD'] = (n_boot_valid, 'valid ensemble members')
            ehdr['PSFEFRAC'] = (round(err_frac, 8), 'integrated PSF error fraction')
            write_fits(arr, ehdr, dest)

    info_json.update(psf_json, a.sample, a.lens, a.filt, {
        'method': method_used,
        'n_stars': n_stars,
        'n_crowded_rejected': qa['n_crowded_rejected'],
        'fwhm_pix': round(fwhm_pix, 4) if fwhm_pix is not None else None,
        'oversample': oversample,
        'kernel_size': int(kernel.shape[0]),
        'trimmed_kernel_size': int(trim_size),
        'trim_threshold': a.trim_threshold,
        'pedestal_frac': round(pedestal_frac, 6),
        'sky_level': round(sky_level, 9) if method_used == 'empirical' else None,
        'kernel_sampling': 'point',
        'fwhm_pix_err': round(fwhm_err, 4) if fwhm_err is not None else None,
        'psf_err_frac': round(err_frac, 8) if err_frac is not None else None,
        'err_method': err_method,
        'err_source': err_source,
        'n_boot_valid': n_boot_valid,
    })
    print('  done')

    # Model tier lacks drizzle broadening (see make_psf_inject.py docstring) -- auto-chain
    # the injection + promotion so this single command already leaves the drizzle-broadened
    # kernel as the canonical psf_kernel.fits / psf.fits. Empirical builds are
    # already the drizzled PSF (cut from the mosaic), so this is skipped for them. A failed
    # injection (e.g. data/drizzle_files/ was cleared) degrades to the cheap analytic
    # drop-box broadening fallback (psf_models.analytic_drop_broaden) rather than silently
    # leaving the sharp, un-broadened analytic model canonical -- and only if that fallback
    # itself fails does the analytic model stay canonical.
    if method_used.startswith('model'):
        try:
            import make_psf_inject
            print(f'  model-tier build ({method_used}) -> running drizzle-broadened '
                  f'injection + promoting to canonical ...')
            make_psf_inject.run_injection(a.lens, a.filt, a.sample,
                                          drizzle_pass=drizzle_pass,
                                          trim_threshold=a.trim_threshold, promote=True)
        except Exception as exc:
            print(f'  WARNING: injection failed ({exc}); trying the cheap analytic '
                  f'drizzle-broadening fallback instead ...')
            try:
                record = make_psf_inject.analytic_broadened_fallback(
                    output_dir, stamps_dir, sci_hdr, method_used, drizzle_pass,
                    trim_threshold=a.trim_threshold)
                info_json.update(psf_json, a.sample, a.lens, a.filt, dict(record))
                print(f'  promoted: broadened_{method_used} is now the canonical PSF '
                      f'(analytic model kept as *_analytic)')
            except Exception as exc2:
                print(f'  WARNING: analytic broadening fallback also failed ({exc2}); '
                      f'canonical PSF stays the un-broadened analytic model ({method_used})')

    # ── The empirical check on a model-tier product ──────────────────────────────────
    # The field-star ePSF (where there were enough stars) beside the canonical kernel now on
    # disk, as psf_kernel_empirical.fits plus `empirical_check` in the record. Expect its
    # core 10-20% more peaked than the model's (the ERR-weight effect, see above) and a
    # fuller wing: the model ePSFs are ~0.6" wide, the real star has light beyond that.
    if method_used.startswith('model') and emp is not None:
        canon = fits.getdata(os.path.join(output_dir, 'psf_kernel.fits')).astype(float)
        check = compare_kernels(emp['kernel'], canon)
        check.update(n_stars=n_stars, pedestal_frac=round(emp['pedestal_frac'], 6),
                     poor_wings=bool(emp['poor']))
        vhdr = fits.Header()
        vhdr['PSFMETH'] = ('empirical', 'field-star ePSF: the check on the shipped model')
        vhdr['PSFNSTAR'] = (n_stars, 'stars used')
        vhdr['PSFPXSCL'] = (round(scale, 5), 'kernel pixel scale (arcsec)')
        vhdr['PSFPED'] = (round(emp['pedestal_frac'], 6), 'sky residual removed (fraction of peak)')
        write_fits(emp['kernel'], vhdr, os.path.join(output_dir, 'psf_kernel_empirical.fits'))
        info_json.update(psf_json, a.sample, a.lens, a.filt,
                         dict(info_json.load(psf_json)[a.sample][a.lens][a.filt], empirical_check=check))
        print(f'  check: empirical ({n_stars} stars) FWHM {check["fwhm_pix"]} px vs model '
              f'{check["fwhm_pix_model"]} px, peak ratio {check["peak_ratio"]:.3f}, core L1 {check["l1_core"]:.3f}')
    elif method_used.startswith('model'):
        f = os.path.join(output_dir, 'psf_kernel_empirical.fits')
        if os.path.exists(f):
            os.remove(f)                       # no check this time; do not leave a stale one

    # ── The HYBRID kernel: model core + empirical wing (canonical where both exist) ────
    # The canonical file at this point is the injected (or analytically broadened) model.
    # With a usable ePSF and a wing_radius for the camera, the shipped kernel becomes the
    # hybrid; the pure model is kept as psf_kernel_model.fits, the blend weight as
    # psf_kernel_blend.fits, and the ePSF's bootstrap/jackknife error map as
    # psf_kernel_empirical_err.fits (make_psf_err_injected.py combines them).
    wing_r = params.get('wing_radius')
    for stale in ('psf_kernel_model.fits', 'psf_kernel_blend.fits', 'psf_kernel_empirical_err.fits'):
        if os.path.exists(os.path.join(output_dir, stale)):
            os.remove(os.path.join(output_dir, stale))
    if method_used.startswith('model') and emp is not None and not emp['poor'] and wing_r:
        rec = info_json.load(psf_json)[a.sample][a.lens][a.filt]
        model_method = rec['method']
        canon_path = os.path.join(output_dir, 'psf_kernel.fits')
        with fits.open(canon_path) as hdul:
            model_kernel = np.asarray(hdul[0].data, float)
            model_hdr = hdul[0].header.copy()
        shutil.copy(canon_path, os.path.join(output_dir, 'psf_kernel_model.fits'))
        hyb, w_blend, f_wing = hybrid_kernel(model_kernel, emp['kernel'], float(wing_r))
        trimmed_h, trim_h = trim_kernel_to_amplitude(hyb, a.trim_threshold)
        fwhm_h = measure_fwhm(hyb)
        hhdr = model_hdr.copy()
        hhdr['PSFMETH'] = (f'hybrid_{model_method}', 'model core + empirical wing')
        hhdr['PSFMODEL'] = (model_method, 'the model the core comes from (psf_kernel_model.fits)')
        hhdr['PSFNSTAR'] = (n_stars, 'stars in the empirical wing')
        hhdr['PSFWINGR'] = (float(wing_r), 'model->ePSF switch radius (px), +-%g px ramp' % HYBRID_TAPER)
        hhdr['PSFWINGF'] = (round(f_wing, 6), 'flux fraction in the empirical wing')
        hhdr['PSFFWHM'] = (round(fwhm_h, 4) if fwhm_h else 0.0, 'fitted kernel FWHM (pixels)')
        hhdr['PSFKIND'] = ('full', 'full kernel (trimmed copy in the psf stamp dir)')
        write_fits(hyb, hhdr, canon_path)
        bhdr = fits.Header()
        bhdr['PSFKIND'] = ('blend', 'hybrid weight of the model core (1) vs the ePSF wing (0)')
        write_fits(w_blend, bhdr, os.path.join(output_dir, 'psf_kernel_blend.fits'))
        plot_psf(hyb, emp['epsf'], emp['star_stamps'], scale, fwhm_h, os.path.join(output_dir, 'psf.png'),
                 title=f'{a.lens}  {a.filt}  [hybrid: {model_method} core + ePSF wing ({n_stars} stars)]')
        thdr = hhdr.copy()
        thdr['PSFKIND'] = ('trimmed', 'amplitude-trimmed modelling kernel')
        thdr['PSFTRIM'] = (a.trim_threshold, 'azimuthal PSF < this x peak sets the radius')
        thdr['PSFPASS'] = (drizzle_pass, 'drizzle pass the PSF matches')
        if star_pass != drizzle_pass:
            thdr['PSFSTARP'] = (star_pass, 'drizzle pass the ePSF wing was built from')
        write_fits(trimmed_h, thdr, os.path.join(stamps_dir, 'psf.fits'))
        for stale in ('psf_kernel_err.fits', os.path.join(stamps_dir, 'psf_err.fits')):
            f = stale if os.path.isabs(stale) else os.path.join(output_dir, stale)
            if os.path.exists(f):
                os.remove(f)             # described the pure model; make_psf_err_injected rewrites it
        # the ePSF wing's own uncertainty: the star-resample ensemble, as for an empirical build
        emp_err_frac = None
        if do_psf_err and n_stars >= 2:
            n_boot = a.n_boot if a.n_boot is not None else int(overrides.get('n_boot', 100))
            print(f'  ePSF wing uncertainty: resampling {n_stars} stars '
                  f'({"jackknife" if n_stars < JACKKNIFE_MAX_STARS else f"bootstrap x{n_boot}"}) ...')
            ks, emeth, nvalid = build_epsf_ensemble(data, stars_tbl, star_size, oversample, params['maxiters'],
                                                    n_boot, level=emp['sky_level'])
            if nvalid >= 2:
                e_full, _, emp_err_frac, _ = psf_error_map(ks, emeth, (trim_h - 1) // 2)
                ehdr = fits.Header()
                ehdr['PSFKIND'] = ('empirical_err', 'per-pixel std of the ePSF (star resample)')
                ehdr['PSFERR'] = (emeth, 'bootstrap or jackknife')
                ehdr['PSFNVALD'] = (nvalid, 'valid ensemble members')
                write_fits(e_full, ehdr, os.path.join(output_dir, 'psf_kernel_empirical_err.fits'))
        rec.update(method=f'hybrid_{model_method}', model_method=model_method,
                   fwhm_pix=round(fwhm_h, 4) if fwhm_h else None, kernel_size=int(hyb.shape[0]),
                   trimmed_kernel_size=int(trim_h), wing_radius_px=float(wing_r), wing_frac=round(f_wing, 6),
                   wing_n_stars=n_stars, wing_err_frac=round(emp_err_frac, 8) if emp_err_frac else None)
        info_json.update(psf_json, a.sample, a.lens, a.filt, rec)
        print(f'  HYBRID shipped: {model_method} core inside {wing_r:g} px, ePSF wing beyond '
              f'({f_wing:.3f} of the flux; {n_stars} stars); FWHM {fwhm_h:.2f} px, trimmed {trim_h} px')


if __name__ == '__main__':
    main()
