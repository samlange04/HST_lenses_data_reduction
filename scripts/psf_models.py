#!/usr/bin/env python
"""
Model-PSF fallback for make_psf.py: evaluate an STScI STDPSF library grid when a field
has too few field stars for an empirical ePSF (SLACS WFPC2/WF3 fields in particular are
star-poor -- Anderson & King 2000 build the ePSF from rich globular-cluster fields, which
these are not, so a library model is the right fallback there).

The STDPSFs (Anderson's empirical library; Dauphin et al. 2021, the WFPC2/WFC3 PSF
Database ISR) are 4x-supersampled 101x101 ePSF grids, read natively by photutils
(GriddedPSFModel.read(..., format='stdpsf')). Two facts from the ISR shape this module:

  * The library does NOT carry every filter. Neither our WFPC2 F606W nor our ACS F555W
    has an exact grid; the nearest available band is substituted (WFPC2 F606W -> F555W,
    ACS F555W -> F606W), with a warning. The optical PSF varies slowly enough with
    wavelength over ~50 nm for a fallback model.
  * WFPC2 grids hold a 3x3 array of fiducial PSFs per chip, so a chip must be selected
    (WF3 = detector 3); ACS/WFC holds two chips (default chip 1). WFC3/IR is a single
    detector.

The delivered PSF is evaluated at the detector grid centre -- the drizzled-mosaic pixel
coordinate is not in the chip frame, and the across-chip spatial variation is a secondary
effect for a fallback -- then resampled from the detector's native pixel scale to the
drizzled output scale (`out_scale`).

Caveat (Anderson 2016, WFC3/IR ISR): the STDPSF is the *detector-frame* ePSF. The true
drizzled-frame PSF differs -- AstroDrizzle resampling broadens point sources -- and the
rigorous route is to inject artificial stars into the flt and drizzle them. This module's
resample omits that broadening, so the model kernel runs slightly sharp. It is a fallback
for star-poor fields; the empirical ePSF (built from stars in the drizzled image itself)
is the true drizzled PSF and is always preferred when enough stars exist.

Grids are cached under data/reference_files/stdpsf/ and downloaded once from STScI.
"""

import glob
import json
import os
import re
import urllib.request

import numpy as np
from photutils.psf import GriddedPSFModel

ws_path = '/Users/samlange/Code/HST_lenses_data_reduction'
_CACHE = os.path.join(ws_path, 'data', 'reference_files', 'stdpsf')
_STDPSF_BASE = 'https://www.stsci.edu/~jayander/HST1PASS/LIB/PSFs/STDPSFs'

# Cache for the ACS/WFC focus-diverse ePSF grids (acstools.focus_diverse_epsfs).
_ACS_FD_CACHE = os.path.join(ws_path, 'data', 'reference_files', 'acs_fdpsf')
# STDPBF focus-diverse grids are supersampled 4x relative to the detector pixel scale.
_ACS_FD_SUPERSAMPLE = 4

# Native detector pixel scale (arcsec) each STDPSF grid is defined on. The delivered
# kernel is resampled from this to the drizzled output scale, or the model comes out the
# wrong size (e.g. WFC3/IR 0.1283" vs the 0.06" drizzled grid).
_DET_SCALE = {'ACS/WFC': 0.05, 'WFC3/IR': 0.1283, 'WFPC2': 0.0996, 'WFC3/UVIS': 0.0396}

# STDPSF subdir, WFPC2/ACS chip to select, and the filters actually published for each
# instrument (verified against the STScI listing). Missing filters fall back to the
# nearest by pivot wavelength (below).
_LIB = {
    'ACS/WFC': dict(subdir='ACSWFC', detector_id=1,
                    filters=['F435W', 'F475W', 'F606W', 'F775W', 'F814W']),
    'WFC3/IR': dict(subdir='WFC3IR', detector_id=None,
                    filters=['F098M', 'F105W', 'F110W', 'F125W', 'F127M', 'F140W',
                             'F153M', 'F160W']),
    'WFPC2':   dict(subdir='WFPC2', detector_id=3,   # WF3 -- the lens chip
                    filters=['F555W', 'F658N', 'F675W', 'F814W']),
    # WFC3/UVIS is a two-CCD optical detector; its STDPSF grid stacks both chips, so a
    # detector_id (1 or 2, like ACS/WFC) must be selected -- chip 1 by default, the same
    # secondary-effect choice as ACS. Every gallery band (F225W/F275W/F438W/F606W/F814W)
    # has an EXACT published grid, so no nearest-filter substitution ever fires here.
    'WFC3/UVIS': dict(subdir='WFC3UV', detector_id=1,
                      filters=['F225W', 'F275W', 'F336W', 'F390W', 'F410M', 'F438W',
                               'F467M', 'F555W', 'F606W', 'F621M', 'F775W', 'F814W',
                               'F850L']),
}

# Filter pivot wavelengths (nm), enough to pick a nearest neighbour.
_PIVOT = {
    'F225W': 236, 'F275W': 271, 'F336W': 335, 'F390W': 392, 'F410M': 411,
    'F435W': 431, 'F438W': 433, 'F467M': 468, 'F475W': 477, 'F555W': 539,
    'F606W': 589, 'F621M': 622, 'F625W': 632,
    'F658N': 658, 'F675W': 673, 'F775W': 765, 'F814W': 802, 'F850L': 917,
    'F098M': 986, 'F105W': 1055, 'F110W': 1153, 'F125W': 1248, 'F127M': 1274,
    'F140W': 1392, 'F153M': 1531, 'F160W': 1537,
}

_GRID_CACHE = {}


# ── Detector-frame -> North-up rotation ──────────────────────────────────────────
# Every model ePSF here (STDPSF, the WFPC2 F606W DB build, the ACS focus-diverse ePSF)
# is in the *detector* frame: its x/y axes are the exposure's detector axes, so its
# diffraction spikes and asymmetric wing structure are rotated by the exposure roll
# (ORIENTAT, up to ~105 deg for SLACS) relative to the North-up drizzled science image
# (final_rot=0.0). Resampling to the output scale WITHOUT this rotation leaves the model
# PSF misoriented against the data. The empirical ePSF is immune -- it is built from the
# North-up mosaic itself. We rotate via the exposure CD matrix rather than an ORIENTAT
# angle so rotation, parity (detector-to-sky handedness) and scale are all handled at once.

def _northup_M(cd_det, out_scale, oversample, pix_per_src):
    """2x2 map: output-subpixel index offset [i(col); j(row)] -> source-pixel offset
    [dcol; drow], taking a detector-frame ePSF into the North-up output grid.

    `cd_det` is the exposure CD (deg per detector pixel; its columns are how detector
    x,y map to sky). The output grid is North-up at `out_scale` arcsec/output-pixel
    (CD_out = diag(-s, s), matching every drizzled product: CD1_1<0, CD2_2>0). `pix_per_src`
    is source pixels per detector pixel (the ePSF supersampling; 1 when sampling a
    GriddedPSFModel directly in detector pixels). Reduces to an isotropic scale when the
    detector frame is itself North-up.
    """
    s = out_scale / 3600.0
    cd_out = np.array([[-s, 0.0], [0.0, s]])
    return pix_per_src * (np.linalg.inv(np.asarray(cd_det, float)) @ cd_out) / oversample


def _grid_from_M(cx, cy, M, n):
    """Sampling coords (xx=col, yy=row) for an n x n stamp centred at (cx, cy) under map M."""
    offs = np.arange(n) - (n - 1) / 2.0
    I, J = np.meshgrid(offs, offs)          # I: col(x) offset, J: row(y) offset
    xx = cx + M[0, 0] * I + M[0, 1] * J
    yy = cy + M[1, 0] * I + M[1, 1] * J
    return xx, yy


def _base_filter(filt):
    """Strip a split-visit suffix (f606W_v1 -> f606W) for filter-library lookups.

    Multi-visit lenses (J0728, J0822) are keyed per visit at the product-directory
    level, but the optics are identical across visits, so the model PSF is the base
    filter's. Without this the STDPSF path KeyErrors on the raw directory name.
    """
    return re.sub(r'_v\d+$', '', filt)


def _resolve_filter(inst_key, filt):
    """Exact STDPSF filter if published, else the nearest by pivot wavelength."""
    if inst_key not in _LIB:
        raise KeyError(f'no STDPSF library mapping for instrument {inst_key!r}')
    fu = _base_filter(filt).upper()
    avail = _LIB[inst_key]['filters']
    if fu in avail:
        return fu, False
    if fu not in _PIVOT:
        raise KeyError(f'unknown filter {filt!r}; add its pivot wavelength to _PIVOT')
    nearest = min(avail, key=lambda f: abs(_PIVOT[f] - _PIVOT[fu]))
    return nearest, True


def _grid_path(inst_key, stdpsf_filt):
    subdir = _LIB[inst_key]['subdir']
    fname = f'STDPSF_{subdir}_{stdpsf_filt}.fits'
    os.makedirs(_CACHE, exist_ok=True)
    local = os.path.join(_CACHE, fname)
    if os.path.exists(local):
        return local
    url = f'{_STDPSF_BASE}/{subdir}/{fname}'
    print(f'  fetching STDPSF grid: {url}')
    try:
        tmp = local + '.part'
        urllib.request.urlretrieve(url, tmp)
        os.replace(tmp, local)
    except Exception as exc:
        raise SystemExit(
            f'ERROR: could not obtain the STDPSF model grid {fname}.\n'
            f'       Tried: {url}\n       Reason: {exc}\n'
            f'       Download it by hand into {_CACHE}/ and re-run, or build an\n'
            f'       empirical PSF instead (add include stars to info/psf_stars.json).')
    return local


def _load_grid(inst_key, stdpsf_filt):
    ck = (inst_key, stdpsf_filt)
    if ck not in _GRID_CACHE:
        path = _grid_path(inst_key, stdpsf_filt)
        _GRID_CACHE[ck] = GriddedPSFModel.read(
            path, format='stdpsf', detector_id=_LIB[inst_key]['detector_id'])
    return _GRID_CACHE[ck]


def model_psf(inst_key, filt, oversample, size, out_scale, cd_detector=None):
    """Oversampled model-PSF stamp for (inst_key, filt) at the detector grid centre.

    Returns an array on a grid `oversample` times finer than the drizzled output scale
    (`out_scale`, arcsec), spanning `size` output pixels -- the same convention as the
    empirical ePSF, so make_psf.oversampled_to_kernel() can bin either to the kernel.

    When `cd_detector` (the exposure CD, deg/detector-pixel) is given, the sampling grid is
    rotated so the detector-frame STDPSF lands in the North-up drizzle frame; without it the
    stamp is only scaled (legacy behaviour), leaving it misoriented by the exposure roll.
    """
    if inst_key not in _DET_SCALE:
        raise KeyError(f'no detector pixel scale recorded for {inst_key!r}')
    stdpsf_filt, substituted = _resolve_filter(inst_key, filt)
    if substituted:
        print(f'  NOTE: no STDPSF for {inst_key} {filt}; using nearest band '
              f'{stdpsf_filt} ({_PIVOT[_base_filter(filt).upper()]}->'
              f'{_PIVOT[stdpsf_filt]} nm)')
    grid = _load_grid(inst_key, stdpsf_filt)

    # Evaluate at the centre of the fiducial-PSF grid (representative; the drizzled-mosaic
    # coordinate is not in the chip frame, and across-chip variation is second order here).
    xp = np.array([p[0] for p in grid.grid_xypos], float)
    yp = np.array([p[1] for p in grid.grid_xypos], float)
    xc, yc = float(xp.mean()), float(yp.mean())

    n = size * oversample
    if cd_detector is not None:
        # grid.evaluate samples in detector pixels, so pix_per_src = 1.
        M = _northup_M(cd_detector, out_scale, oversample, pix_per_src=1.0)
        xx, yy = _grid_from_M(xc, yc, M, n)
    else:
        print('  NOTE: no detector CD available; STDPSF left in the detector frame '
              '(not rotated to North-up)')
        step = (out_scale / _DET_SCALE[inst_key]) / oversample   # detector px per sub-sample
        offs = (np.arange(n) - (n - 1) / 2.0) * step
        xx, yy = np.meshgrid(xc + offs, yc + offs)

    stamp = np.asarray(grid.evaluate(xx, yy, flux=1.0, x_0=xc, y_0=yc), dtype=float)
    if stamp.sum() > 0:
        stamp /= stamp.sum()
    return stamp


# ── Cheap analytic drizzle-broadening (fast fallback for make_psf_inject.py) ─────
# make_psf_inject.py's artificial-star injection + re-drizzle is the rigorous way to give
# a model-tier kernel the resampling broadening AstroDrizzle puts on a real point source
# (Anderson 2016) -- but it needs the persisted data/drizzle_files/ inputs and a real
# re-drizzle per lens. When that isn't available (inputs cleared, injection failed for any
# other reason), this gives a cheap analytic stand-in: convolve the already-built,
# already-North-up analytic kernel with a box the width the drizzle drop projects onto the
# OUTPUT grid -- `pixfrac * native_scale / out_scale` output pixels -- in Fourier space
# (exact for a fractional box width; PyAutoReduce's psf/frame_combine._drop_convolve does
# the same thing to a detector-frame kernel before its own North-up resample). Doing it
# post-resample, on an already-North-up array, needs no detector CD/Jacobian machinery at
# all -- the real drizzle drop is a square in DETECTOR pixels, which after rotation to
# North-up is a rotated square in output-pixel space, not axis-aligned unless the roll is
# ~0; approximating it as axis-aligned here is exactly the "cheap" trade this function is
# for. It is not a substitute for the real injection where that's available.

def drop_convolve_box(kernel, box_width_pix):
    """Convolve `kernel` with a `box_width_pix`-wide box (axis-aligned, same width both
    axes), done in Fourier space so a fractional width is exact: the box's transform is
    `sinc(box_width_pix * f)` per axis. `box_width_pix <= 0` is a no-op (returns `kernel`
    unchanged) -- native-scale drizzle with no oversampling has nothing to broaden.
    """
    if box_width_pix <= 0:
        return np.asarray(kernel, dtype=float)
    kernel = np.asarray(kernel, dtype=float)
    ny, nx = kernel.shape
    fy = np.fft.fftfreq(ny)[:, None]
    fx = np.fft.rfftfreq(nx)[None, :]
    transfer = np.sinc(box_width_pix * fy) * np.sinc(box_width_pix * fx)
    return np.fft.irfft2(np.fft.rfft2(kernel) * transfer, s=kernel.shape)


def analytic_drop_broaden(kernel, pixfrac, native_scale, out_scale, n_frames=1):
    """Cheap analytic approximation to AstroDrizzle's resampling broadening of a point
    source, applied to an already-North-up, output-pixel-scale kernel (e.g. the analytic
    model kernel make_psf.py just wrote, before injection promotion).

    `n_frames` (contributing dithered exposures) matters: a single frame's drop-box is the
    full `pixfrac * native_scale / out_scale` box, but with N frames at different sub-pixel
    dither phases the box averages down -- the same "interlacing limit" (r -> 1, R -> 1 as
    dithering gets dense) casertano_r assumes in make_cutouts.py -- so the effective box is
    scaled by `1/sqrt(n_frames)`. This is calibrated, not guessed: on J0008-0004 F606W
    (pixfrac=1, native/out=1.99, 4 frames) the n=1 box overshoots the real injected FWHM
    (4.35 -> 4.67px analytic vs the true re-drizzled 4.42px); n=4 lands at 4.43px, matching
    to 0.01px. `n_frames=1` (default) is the conservative no-dither-correction box, for a
    caller that doesn't know or care to pass the frame count.

    Returns `(broadened, box_width_pix)`: the drop-convolved, unit-sum-renormalised kernel
    (clipped to non-negative first -- the Fourier convolution can ring slightly negative at
    the sub-percent level in the wings) and the box width actually used (output pixels,
    post dither-scaling), for the caller to log. See the section docstring above for why
    this is only a cheap stand-in, not the real re-drizzle injection.
    """
    box_width = (pixfrac * native_scale / out_scale) / np.sqrt(max(1, n_frames))
    broadened = drop_convolve_box(kernel, box_width)
    broadened = np.clip(broadened, 0.0, None)
    total = broadened.sum()
    if total > 0:
        broadened = broadened / total
    return broadened, box_width


# ── WFPC2 F606W native ePSF from the MAST PSF database ───────────────────────────
# The STDPSF library carries no WFPC2 F606W grid, so model_psf() substitutes WFPC2
# F555W (right chip, wrong filter) -- the least-verified product in the pipeline. The
# MAST PSF database (Dauphin et al., ISR WFC3 2021-12) instead holds ~140k good,
# unsaturated WF3 F606W point-source cutouts. We build ONE native-F606W WF3 ePSF from
# the best-qfit stars near the lens position -- every WFPC2 F606W lens puts its target
# at the same WF3 spot (~435,424), so a single shared model serves all of them -- and
# cache it. Still a detector-frame ePSF (omits AstroDrizzle broadening), so the
# lens's-own-field empirical build is preferred where stars exist; this is the model
# FALLBACK, slotted ABOVE the STDPSF F555W proxy.

_WFPC2_F606W_DB_CACHE = os.path.join(ws_path, 'data', 'reference_files', 'wfpc2_f606w_psfdb')
_WF3_LENS_XY = (435, 424)     # lens galaxy position on WF3 (chip 3); see CLAUDE.md
_F606W_DB = dict(
    qfit_max=0.05,        # low qfit == good template fit; ~8.6k WF3 F606W stars qualify
    radius=200,           # x_cal/y_cal box half-width around the lens WF3 position
    n_download=300,       # best-qfit candidates to fetch (many are edge-contaminated)
    build_oversample=4,   # EPSFBuilder oversampling relative to the WF3 detector scale
    inner_half=13,        # crop to a 27x27 inner box -> drop edge neighbours / warm pixels
    max_center_off=2.0,   # reject a cutout whose windowed centroid is >2px off centre
    min_stars=15,         # need at least this many clean stars, else fall back to STDPSF
    maxiters=12,
)


def wfpc2_f606w_db_epsf(oversample, size, out_scale, cd_detector=None,
                        force_rebuild=False, return_ensemble=False):
    """Native WF3 F606W ePSF (MAST PSF database) on the make_psf output grid.

    Returns an array `oversample` times finer than `out_scale`, spanning `size` output
    pixels -- the same convention as model_psf(). Builds+caches a single shared *detector-
    frame* ePSF the first time; later calls reload it. The shared ePSF is rotated per lens
    into the North-up output frame via `cd_detector` (that lens's WF3 exposure CD) at resample
    time -- so one cached build serves every roll. Raises RuntimeError if a usable ePSF can't
    be built so make_psf falls back to the STDPSF F555W proxy.

    `return_ensemble=True` additionally returns (ensemble, method): a list of the shared
    build's bootstrap/jackknife-over-DB-stars ensemble members (see
    `_wfpc2_f606w_db_ensemble`), each resampled into this lens's North-up frame the same way
    as the point estimate, plus the resampling method ('bootstrap' or 'jackknife') --
    together in the exact (members, method) convention make_psf.psf_error_map expects. An
    ensemble-build failure degrades to an empty list (point estimate is unaffected).
    """
    data, samp = _wfpc2_f606w_db_build(force_rebuild=force_rebuild)
    point = _resample_centered(data, samp, out_scale, oversample, size,
                               cd_detector=cd_detector, det_scale=_DET_SCALE['WFPC2'])
    if not return_ensemble:
        return point
    try:
        members, msamp, method = _wfpc2_f606w_db_ensemble(force_rebuild=force_rebuild)
        ensemble = [_resample_centered(m, msamp, out_scale, oversample, size,
                                       cd_detector=cd_detector, det_scale=_DET_SCALE['WFPC2'])
                   for m in members]
    except Exception as exc:
        print(f'  WARNING: WFPC2 F606W DB error ensemble unavailable ({exc})')
        ensemble, method = [], None
    return point, ensemble, method


def _wfpc2_f606w_collect_stars(force_rebuild=False):
    """Query + download (cached) the WF3 F606W MAST PSF-DB star cutouts as EPSFStars.

    Shared by the point-estimate build and the bootstrap/jackknife error ensemble so both
    read the same star sample and re-use the same on-disk cutout cache (no repeat downloads).
    """
    import mast_api_psf
    from astropy.io import fits
    from photutils.psf import EPSFStar
    from photutils.centroids import centroid_com

    cfg = _F606W_DB
    lx, ly = _WF3_LENS_XY
    r = cfg['radius']
    cols = ['id', 'rootname', 'filter_1', 'chip', 'x_cal', 'y_cal', 'qfit',
            'n_sat_pixels', 'subarray']
    params = {'filter_1': ['F606W'], 'chip': ['3'], 'n_sat_pixels': ['0'],
              'qfit': [{'min': 0.0, 'max': cfg['qfit_max']}],
              'x_cal': [{'min': lx - r, 'max': lx + r}],
              'y_cal': [{'min': ly - r, 'max': ly + r}]}
    obs = mast_api_psf.mast_query_psf_database(
        'WFPC2', mast_api_psf.set_filters(params), columns=cols)
    if len(obs) == 0:
        raise RuntimeError('MAST PSF DB returned no WF3 F606W stars')
    obs.sort('qfit')
    obs = obs[:cfg['n_download']]
    uris = mast_api_psf.make_dataURIs(obs, 'WFPC2', file_suffix=['c0m'])

    cutdir = os.path.join(_WFPC2_F606W_DB_CACHE, 'cutouts')
    os.makedirs(cutdir, exist_ok=True)
    half = cfg['inner_half']
    stars = []
    for uri in uris:
        fn = os.path.join(cutdir, uri.split('/')[-1])
        if not os.path.exists(fn) or os.path.getsize(fn) < 1000:
            try:
                mast_api_psf.download_request_file(uri, fn)
            except Exception:
                continue
        try:
            d = np.nan_to_num(np.asarray(fits.getdata(fn), float))
        except Exception:
            continue
        cy, cx = (np.array(d.shape) - 1) // 2
        sub = d[cy - half:cy + half + 1, cx - half:cx + half + 1]
        if sub.shape != (2 * half + 1, 2 * half + 1):
            continue
        sub = sub - np.median(d)                 # per-stamp sky
        c = centroid_com(np.clip(sub, 0, None))  # windowed centroid drops edge contaminants
        if not np.all(np.isfinite(c)) or np.hypot(c[0] - half, c[1] - half) > cfg['max_center_off']:
            continue
        peak = sub.max()
        if peak <= 0:
            continue
        stars.append(EPSFStar(sub / peak, cutout_center=(float(c[0]), float(c[1]))))
    if len(stars) < cfg['min_stars']:
        raise RuntimeError(f'only {len(stars)} usable WF3 F606W DB stars '
                           f'(<{cfg["min_stars"]})')
    return stars


def _wfpc2_f606w_db_build(force_rebuild=False):
    """(oversampled ePSF array, arcsec/pixel sampling) for WF3 F606W, cached to FITS."""
    from astropy.io import fits
    from photutils.psf import EPSFStars, EPSFBuilder

    cfg = _F606W_DB
    samp = _DET_SCALE['WFPC2'] / cfg['build_oversample']
    os.makedirs(_WFPC2_F606W_DB_CACHE, exist_ok=True)
    cache = os.path.join(_WFPC2_F606W_DB_CACHE, 'wf3_f606w_epsf.fits')
    if os.path.exists(cache) and not force_rebuild:
        return np.asarray(fits.getdata(cache), float), samp

    stars = _wfpc2_f606w_collect_stars(force_rebuild=force_rebuild)
    builder = EPSFBuilder(oversampling=cfg['build_oversample'], maxiters=cfg['maxiters'],
                          recentering_maxiters=8, progress_bar=False)
    epsf, _ = builder(EPSFStars(stars))
    data = np.nan_to_num(np.asarray(epsf.data, float))

    hdr = fits.Header()
    hdr['PSFSRC'] = ('MAST_PSF_DB', 'WFPC2 F606W PSF DB (Dauphin ISR WFC3 2021-12)')
    hdr['NSTARS'] = (len(stars), 'WF3 F606W star cutouts used')
    hdr['QFITMAX'] = (cfg['qfit_max'], 'max qfit selected')
    hdr['BUILDOVR'] = (cfg['build_oversample'], 'EPSFBuilder oversampling vs WF3 scale')
    hdr['PIXSCALE'] = (round(samp, 5), 'ePSF arcsec/pixel')
    fits.writeto(cache, data, hdr, overwrite=True)
    print(f'    built WF3 F606W MAST-PSF-DB ePSF from {len(stars)} stars -> {cache}')
    return data, samp


# JACKKNIFE_MAX_STARS mirrors make_psf.JACKKNIFE_MAX_STARS (duplicated, not imported --
# psf_models must stay importable standalone / without a make_psf -> psf_models cycle).
_JACKKNIFE_MAX_STARS = 6


def _reduce_ensemble(items, n_boot=100, rng=None):
    """Turn a list of raw per-item arrays (e.g. one North-up kernel per contributing
    exposure) into the "reduced-member" convention make_psf.psf_error_map expects: each
    returned member is itself a MEAN over a resample of `items` -- leave-one-out when there
    are too few items for bootstrap to be non-degenerate (< _JACKKNIFE_MAX_STARS), else
    bootstrap-with-replacement -- exactly mirroring make_psf.build_epsf_ensemble's contract,
    just resampling precomputed items instead of rebuilding from a star table. Returns
    (members, method).
    """
    rng = rng if rng is not None else np.random.default_rng(0)
    n = len(items)
    if n < _JACKKNIFE_MAX_STARS:
        method = 'jackknife'
        index_sets = [[j for j in range(n) if j != i] for i in range(n)]
    else:
        method = 'bootstrap'
        index_sets = [rng.integers(0, n, size=n).tolist() for _ in range(n_boot)]
    members = []
    for idx in index_sets:
        if not idx:
            continue
        m = np.mean([items[j] for j in idx], axis=0)
        if m.sum() > 0:
            m = m / m.sum()
        members.append(m)
    return members, method


def _wfpc2_f606w_db_ensemble(force_rebuild=False, n_boot=100, rng=None):
    """Bootstrap (or jackknife if star-poor) ensemble of the detector-frame, oversampled
    WF3 F606W DB ePSF, over the shared DB star sample -- cached to FITS (one 3D cube) like
    the point estimate. Returns (list_of_arrays, arcsec/pixel sampling, method).

    This is the WFPC2 F606W analogue of make_psf.build_epsf_ensemble (which resamples a
    lens's own field stars): here the star sample is the ~150 shared archival DB stars, so
    the ensemble -- and its cache -- is built ONCE and reused by every WFPC2 F606W lens,
    each of which only re-resamples/rotates it into its own North-up frame at call time.
    """
    from astropy.io import fits
    from photutils.psf import EPSFStars, EPSFBuilder

    cfg = _F606W_DB
    samp = _DET_SCALE['WFPC2'] / cfg['build_oversample']
    cache = os.path.join(_WFPC2_F606W_DB_CACHE, 'wf3_f606w_epsf_ensemble.fits')
    if os.path.exists(cache) and not force_rebuild:
        cube = np.asarray(fits.getdata(cache), float)
        method = fits.getheader(cache).get('ENSMETH', 'bootstrap')
        return [cube[i] for i in range(cube.shape[0])], samp, method

    stars = _wfpc2_f606w_collect_stars(force_rebuild=force_rebuild)
    n = len(stars)
    rng = rng if rng is not None else np.random.default_rng(0)
    if n < _JACKKNIFE_MAX_STARS:
        method = 'jackknife'
        index_sets = [[j for j in range(n) if j != i] for i in range(n)]
    else:
        method = 'bootstrap'
        index_sets = [rng.integers(0, n, size=n).tolist() for _ in range(n_boot)]

    builder = EPSFBuilder(oversampling=cfg['build_oversample'], maxiters=cfg['maxiters'],
                          recentering_maxiters=8, progress_bar=False)
    members = []
    for idx in index_sets:
        if not idx:
            continue
        try:
            epsf, _ = builder(EPSFStars([stars[j] for j in idx]))
            d = np.nan_to_num(np.asarray(epsf.data, float))
            if np.all(np.isfinite(d)) and d.sum() > 0:
                members.append(d)
        except Exception:
            continue                       # a degenerate resample is dropped, not fatal
    if len(members) < 2:
        raise RuntimeError(f'only {len(members)} valid WF3 F606W DB ensemble member(s)')

    cube = np.stack(members, axis=0)
    hdr = fits.Header()
    hdr['PSFSRC'] = ('MAST_PSF_DB_ENS', 'WF3 F606W PSF DB error ensemble')
    hdr['ENSMETH'] = (method, 'bootstrap or jackknife over DB stars')
    hdr['NMEMBER'] = (len(members), 'valid ensemble members')
    hdr['NSTARS'] = (n, 'DB stars in the shared sample')
    hdr['PIXSCALE'] = (round(samp, 5), 'ePSF arcsec/pixel')
    fits.writeto(cache, cube, hdr, overwrite=True)
    print(f'    built WF3 F606W DB error ensemble: {len(members)} members '
          f'({method}) -> {cache}')
    return members, samp, method


# ── ACS/WFC focus-diverse ePSF (preferred ACS model over the static STDPSF) ──────
# The STDPSF library carries no ACS/WFC F555W grid and is one static, grid-centre,
# fixed-focus PSF. acstools.focus_diverse_epsfs instead serves an *observation-matched*,
# focus-corrected empirical ePSF (Bellini et al., ACS ISR 2018-08 / 2023-06) keyed to each
# exposure rootname -- native F555W and F814W (no nearest-filter substitution), matched to
# the HST focus (breathing) at the time of each exposure, and interpolable to the lens's
# actual detector position. The drizzled PSF is the exposure-average, so we retrieve the
# ePSF per contributing exposure, interpolate each to the lens position on that exposure's
# detector, average, then resample the 4x-supersampled detector grid to the drizzled
# output scale -- the same output convention as model_psf(). It stays a *fallback*: like
# any library PSF it is the detector-frame ePSF and omits AstroDrizzle broadening (Anderson
# 2016), so an empirical ePSF built from stars in the drizzled image is still preferred
# when the field has enough. Raises on any unavailability so make_psf can fall back to the
# STDPSF model_psf().

def _fd_detector_position(rootname, calibrated_dir, catalogue_coord):
    """(x, y, chip, cd) of the lens on an ACS/WFC exposure, or None if unresolved.

    Reads the FLC's per-chip SCI WCS and returns the 1-indexed detector coordinate on
    whichever chip contains the lens, naming it 'WFC1' (CCDCHIP=1) or 'WFC2' (CCDCHIP=2)
    as acstools.interp_epsf expects, plus that chip's 2x2 CD matrix (deg/detector-pixel)
    for the North-up rotation. None -> caller uses the grid centre and no rotation.
    """
    from astropy.io import fits
    from astropy.wcs import WCS
    flc = os.path.join(calibrated_dir, f'{rootname}_flc.fits')
    if not os.path.isfile(flc):
        return None
    try:
        with fits.open(flc) as hdul:
            for hdu in hdul:
                if hdu.header.get('EXTNAME') != 'SCI':
                    continue
                h = hdu.header
                x, y = WCS(h, hdul).world_to_pixel(catalogue_coord)
                if 0.5 < x < 4096.5 and 0.5 < y < 2048.5:
                    chip = 'WFC1' if int(h.get('CCDCHIP', 1)) == 1 else 'WFC2'
                    cd = np.array([[h['CD1_1'], h['CD1_2']],
                                   [h['CD2_1'], h['CD2_2']]], float)
                    return float(x), float(y), chip, cd
    except Exception:
        return None
    return None


def _resample_centered(src, src_scale, out_scale, oversample, size,
                       cd_detector=None, det_scale=None):
    """Cubic-resample `src` (centred on its centroid) onto the make_psf output grid.

    `src` is sampled at `src_scale` arcsec/pixel; the returned (size*oversample)^2 stamp is
    at `out_scale/oversample` arcsec/pixel -- matching model_psf()'s convention so
    oversampled_to_kernel() can bin it. Unit-sum normalised.

    When `cd_detector` (exposure CD, deg/detector-pixel) and `det_scale` (arcsec/detector-
    pixel) are given, the detector-frame `src` is rotated into the North-up output frame;
    otherwise it is only scaled (legacy behaviour, leaving the roll misorientation in place).
    """
    from photutils.centroids import centroid_com
    from scipy.ndimage import map_coordinates
    src = np.asarray(src, dtype=float)
    ny, nx = src.shape
    # centre on the ePSF core (fall back to the geometric centre if the centroid is bad)
    py, px = np.unravel_index(np.argmax(src), src.shape)
    win = 15
    y0, x0 = max(0, py - win // 2), max(0, px - win // 2)
    cy, cx = centroid_com(src[y0:y0 + win, x0:x0 + win])
    cx, cy = x0 + cx, y0 + cy
    if not (np.isfinite(cx) and np.isfinite(cy)):
        cx, cy = nx / 2.0, ny / 2.0

    n = size * oversample
    if cd_detector is not None and det_scale is not None:
        # src is supersampled relative to the detector by det_scale/src_scale.
        M = _northup_M(cd_detector, out_scale, oversample,
                       pix_per_src=det_scale / src_scale)
        xx, yy = _grid_from_M(cx, cy, M, n)
    else:
        step = (out_scale / oversample) / src_scale     # source pixels per output subpixel
        offs = (np.arange(n) - (n - 1) / 2.0) * step
        xx, yy = np.meshgrid(cx + offs, cy + offs)
    stamp = map_coordinates(src, [yy.ravel(), xx.ravel()], order=3,
                            mode='constant', cval=0.0).reshape(n, n)
    stamp = np.clip(stamp, 0.0, None)
    if stamp.sum() > 0:
        stamp /= stamp.sum()
    return stamp


def acs_focus_diverse_psf(lens, filt, rootnames, calibrated_dir, catalogue_coord,
                          oversample, size, out_scale, return_ensemble=False):
    """Focus-diverse ACS/WFC ePSF model stamp, averaged over the contributing exposures.

    `rootnames` are the exposures that reached the drizzle (from info/lens_products.json);
    `calibrated_dir` holds their FLCs (for the per-exposure detector position). Returns an
    oversampled stamp on the same grid as model_psf(). Raises RuntimeError if no exposure's
    focus-diverse ePSF could be retrieved, so make_psf falls back to the STDPSF model.

    `return_ensemble=True` additionally returns (ensemble, method): leave-one-exposure-out
    (or bootstrap, if there are ever >= _JACKKNIFE_MAX_STARS exposures) reduced members over
    the per-exposure North-up kernels that went into the average -- the natural error
    ensemble for this tier, since each exposure already carries its own HST focus/breathing
    and that's the only thing this model averages over. `_reduce_ensemble` puts these in the
    same (members, method) convention make_psf.psf_error_map expects, so no downstream
    special-casing is needed vs. the empirical tier's own star-resample ensemble.
    """
    from acstools.focus_diverse_epsfs import psf_retriever, interp_epsf
    from astropy.io import fits

    if not rootnames:
        raise RuntimeError('no exposure rootnames in info/lens_products.json for '
                           f'{lens} {filt}; cannot build a focus-diverse ePSF')
    os.makedirs(_ACS_FD_CACHE, exist_ok=True)
    src_scale = _DET_SCALE['ACS/WFC'] / _ACS_FD_SUPERSAMPLE   # arcsec per supersampled pixel

    # Each exposure's ePSF is detector-frame at that exposure's roll. The drizzled PSF is
    # the exposure-average IN THE NORTH-UP FRAME, so we rotate every exposure to North-up
    # (via its own chip CD) BEFORE averaging -- correct even if the exposures span rolls.
    northup = []
    for root in rootnames:
        try:
            path = psf_retriever(root, _ACS_FD_CACHE)
        except Exception as exc:
            print(f'    focus-diverse ePSF retrieval failed for {root}: {exc}')
            continue
        grid = fits.getdata(path, ext=0)
        pos = _fd_detector_position(root, calibrated_dir, catalogue_coord)
        if pos is None:
            x, y, chip, cd = 2048.0, 1024.0, 'WFC1', None   # grid centre, no rotation
            print(f'    {root}: FLC position/CD unavailable, using grid centre unrotated')
        else:
            x, y, chip, cd = pos
        P = interp_epsf(grid, int(round(x)), int(round(y)), chip)
        if P is None or not np.all(np.isfinite(P)) or np.asarray(P).sum() <= 0:
            print(f'    {root}: interp_epsf returned no usable ePSF; skipping')
            continue
        P = np.asarray(P, dtype=float)
        northup.append(_resample_centered(P / P.sum(), src_scale, out_scale, oversample,
                                          size, cd_detector=cd,
                                          det_scale=_DET_SCALE['ACS/WFC']))

    if not northup:
        raise RuntimeError(f'no focus-diverse ePSF could be retrieved for {lens} {filt}')

    mean_epsf = np.mean(northup, axis=0)
    print(f'    focus-diverse ePSF from {len(northup)}/{len(rootnames)} exposures')
    if mean_epsf.sum() > 0:
        mean_epsf = mean_epsf / mean_epsf.sum()
    if not return_ensemble:
        return mean_epsf
    if len(northup) < 2:
        print(f'    only {len(northup)} exposure(s); no PSF error ensemble')
        return mean_epsf, [], None
    ensemble, method = _reduce_ensemble(northup)
    return mean_epsf, ensemble, method
