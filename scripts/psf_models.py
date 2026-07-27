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
_DET_SCALE = {'ACS/WFC': 0.05, 'WFC3/IR': 0.1283, 'WFPC2': 0.0996}

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
}

# Filter pivot wavelengths (nm), enough to pick a nearest neighbour.
_PIVOT = {
    'F435W': 431, 'F475W': 477, 'F555W': 539, 'F606W': 589, 'F625W': 632,
    'F658N': 658, 'F675W': 673, 'F775W': 765, 'F814W': 802,
    'F098M': 986, 'F105W': 1055, 'F110W': 1153, 'F125W': 1248, 'F127M': 1274,
    'F140W': 1392, 'F153M': 1531, 'F160W': 1537,
}

_GRID_CACHE = {}


def _resolve_filter(inst_key, filt):
    """Exact STDPSF filter if published, else the nearest by pivot wavelength."""
    if inst_key not in _LIB:
        raise KeyError(f'no STDPSF library mapping for instrument {inst_key!r}')
    fu = filt.upper()
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


def model_psf(inst_key, filt, oversample, size, out_scale):
    """Oversampled model-PSF stamp for (inst_key, filt) at the detector grid centre.

    Returns an array on a grid `oversample` times finer than the drizzled output scale
    (`out_scale`, arcsec), spanning `size` output pixels -- the same convention as the
    empirical ePSF, so make_psf.oversampled_to_kernel() can bin either to the kernel.
    """
    if inst_key not in _DET_SCALE:
        raise KeyError(f'no detector pixel scale recorded for {inst_key!r}')
    stdpsf_filt, substituted = _resolve_filter(inst_key, filt)
    if substituted:
        print(f'  NOTE: no STDPSF for {inst_key} {filt}; using nearest band '
              f'{stdpsf_filt} ({_PIVOT[filt.upper()]}->{_PIVOT[stdpsf_filt]} nm)')
    grid = _load_grid(inst_key, stdpsf_filt)

    # Evaluate at the centre of the fiducial-PSF grid (representative; the drizzled-mosaic
    # coordinate is not in the chip frame, and across-chip variation is second order here).
    xp = np.array([p[0] for p in grid.grid_xypos], float)
    yp = np.array([p[1] for p in grid.grid_xypos], float)
    xc, yc = float(xp.mean()), float(yp.mean())

    step = (out_scale / _DET_SCALE[inst_key]) / oversample   # detector px per sub-sample
    n = size * oversample
    offs = (np.arange(n) - (n - 1) / 2.0) * step
    xx, yy = np.meshgrid(xc + offs, yc + offs)

    stamp = np.asarray(grid.evaluate(xx, yy, flux=1.0, x_0=xc, y_0=yc), dtype=float)
    if stamp.sum() > 0:
        stamp /= stamp.sum()
    return stamp


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
    """(x, y, chip) of the lens on an ACS/WFC exposure, or None if it can't be resolved.

    Reads the FLC's per-chip SCI WCS and returns the 1-indexed detector coordinate on
    whichever chip contains the lens, naming it 'WFC1' (CCDCHIP=1) or 'WFC2' (CCDCHIP=2)
    as acstools.interp_epsf expects. None -> caller uses the grid centre.
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
                x, y = WCS(hdu.header, hdul).world_to_pixel(catalogue_coord)
                if 0.5 < x < 4096.5 and 0.5 < y < 2048.5:
                    chip = 'WFC1' if int(hdu.header.get('CCDCHIP', 1)) == 1 else 'WFC2'
                    return float(x), float(y), chip
    except Exception:
        return None
    return None


def _resample_centered(src, src_scale, out_scale, oversample, size):
    """Bilinear-resample `src` (centred on its centroid) onto the make_psf output grid.

    `src` is sampled at `src_scale` arcsec/pixel; the returned (size*oversample)^2 stamp is
    at `out_scale/oversample` arcsec/pixel -- matching model_psf()'s convention so
    oversampled_to_kernel() can bin it. Unit-sum normalised.
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

    step = (out_scale / oversample) / src_scale     # source pixels per output subpixel
    n = size * oversample
    offs = (np.arange(n) - (n - 1) / 2.0) * step
    xx, yy = np.meshgrid(cx + offs, cy + offs)
    stamp = map_coordinates(src, [yy.ravel(), xx.ravel()], order=1,
                            mode='constant', cval=0.0).reshape(n, n)
    stamp = np.clip(stamp, 0.0, None)
    if stamp.sum() > 0:
        stamp /= stamp.sum()
    return stamp


def acs_focus_diverse_psf(lens, filt, rootnames, calibrated_dir, catalogue_coord,
                          oversample, size, out_scale):
    """Focus-diverse ACS/WFC ePSF model stamp, averaged over the contributing exposures.

    `rootnames` are the exposures that reached the drizzle (from info/lens_products.json);
    `calibrated_dir` holds their FLCs (for the per-exposure detector position). Returns an
    oversampled stamp on the same grid as model_psf(). Raises RuntimeError if no exposure's
    focus-diverse ePSF could be retrieved, so make_psf falls back to the STDPSF model.
    """
    from acstools.focus_diverse_epsfs import psf_retriever, interp_epsf
    from astropy.io import fits

    if not rootnames:
        raise RuntimeError('no exposure rootnames in info/lens_products.json for '
                           f'{lens} {filt}; cannot build a focus-diverse ePSF')
    os.makedirs(_ACS_FD_CACHE, exist_ok=True)

    supersampled = []
    for root in rootnames:
        try:
            path = psf_retriever(root, _ACS_FD_CACHE)
        except Exception as exc:
            print(f'    focus-diverse ePSF retrieval failed for {root}: {exc}')
            continue
        grid = fits.getdata(path, ext=0)
        pos = _fd_detector_position(root, calibrated_dir, catalogue_coord)
        if pos is None:
            x, y, chip = 2048.0, 1024.0, 'WFC1'   # grid centre when the FLC is unavailable
            print(f'    {root}: FLC position unavailable, using grid centre')
        else:
            x, y, chip = pos
        P = interp_epsf(grid, int(round(x)), int(round(y)), chip)
        if P is None or not np.all(np.isfinite(P)) or np.asarray(P).sum() <= 0:
            print(f'    {root}: interp_epsf returned no usable ePSF; skipping')
            continue
        P = np.asarray(P, dtype=float)
        supersampled.append(P / P.sum())

    if not supersampled:
        raise RuntimeError(f'no focus-diverse ePSF could be retrieved for {lens} {filt}')

    mean_epsf = np.mean(supersampled, axis=0)
    src_scale = _DET_SCALE['ACS/WFC'] / _ACS_FD_SUPERSAMPLE   # arcsec per supersampled pixel
    print(f'    focus-diverse ePSF from {len(supersampled)}/{len(rootnames)} exposures')
    return _resample_centered(mean_epsf, src_scale, out_scale, oversample, size)
