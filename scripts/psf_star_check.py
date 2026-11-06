#!/usr/bin/env python
"""
The injected PSF model against the real PSF stars, one by one: check it, calibrate its
width, and measure its error. Ported from the JWST repo's psf_star_check.py (2026-11-03).

    uv run python scripts/psf_star_check.py --sample slacs --lens J0728+3835 --filt f160W
    uv run python scripts/psf_star_check.py --calibrate          # the extra blur per camera+filter
    uv run python scripts/psf_star_check.py --error-maps         # model-vs-star error maps
    bash scripts/run_psf_star_check_all.sh slacs                 # every band of a sample

Band mode. The band's PSF stars are make_psf.prepare_field's selection, the ones the
empirical check and the hybrid wing use. The model (the same source make_psf_inject uses:
ACS focus-diverse ePSF or STDPSF, WFC3/IR and UVIS STDPSF, WFPC2 the PSF-database ePSF) is
injected noiselessly at each star's own sky position into every staged science frame and
drizzled as the science was (make_psf_inject.injected_mosaic). Each real star therefore has
a twin at the same sub-pixel position, on the same frames, with the same ERR weights, and
nothing else differs -- a like-for-like test of the model that no ePSF build can give. Per
star it records the Gaussian FWHM of the real star and its twin and their ratio, the peak
fractions, the core amplitude fit real = a * twin + b, and `extended` (real more than
EXTENDED_RATIO broader than its twin: a compact galaxy or a blend). Crops go to
data/psf/<...>/star_check.npz, the summary to info/provenance/lens_psf_starcheck.json.

--calibrate. Pools every band's crops per camera+filter (stars with real/twin FWHM within
EXTENDED_RATIO either way and peak S/N >= CAL_MIN_SNR). Per star, the extra Gaussian blur
of the twin that minimises the core chi^2 of real = a * (twin * G) + b; the filter's value is
the median over its stars. It is converted to DETECTOR pixels and written to
info/overrides/psf_model_calibration.json (`blur_sigma_det_px`), which
psf_models.model_blur_sigma applies to every later injection and analytic model. A blur on
the resampled star stands in for one on the detector only approximately, so iterate: re-run
the bands (the twins now carry the blur), --calibrate again; with two or more passes the next
width is the secant zero of residual^2 against width^2 (logged in `iterations`). Stop when the
residual is a few mas. Measured before this existed: ACS/WFC focus-diverse 1.007, STDPSF
F814W 0.994, WFC3/IR F160W 1.067 (real/twin FWHM): only F160W needs it.

--error-maps. Per camera+filter, over the same stars, the per-pixel rms residual of the
flux-normalised real star about its twin beyond its own noise (sigma-clipped mean of
r^2 - e^2) goes to data/psf/_model_error/<INST>_<FILT>.fits; make_psf_err_injected.py ships
it as the model term of the injected / hybrid kernel's error map (err_method
model_vs_stars). Per filter, not per lens, at the stars' sub-pixel phases.
"""
import argparse
import contextlib
import io
import json
import os
import sys
import time
import warnings

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
ws_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

import make_psf
import make_psf_inject
import psf_models
import info_json
import mast_target_names
import stamp_paths

EXTENDED_RATIO = 1.15       # real/twin FWHM above this: not a point source
CAL_MIN_SNR = 100.0         # peak S/N of a star used for calibration and error maps
CROP = 31                   # star crops (output px) for the record; the calibration uses r_fit
SIGMA_GRID = np.arange(0.0, 1.501, 0.05)   # extra blur grid, OUTPUT px
STARCHECK_JSON = os.environ.get('PSF_STARCHECK_JSON') or os.path.join(ws_path, 'info', 'provenance', 'lens_psf_starcheck.json')
MODEL_ERROR_DIR = os.path.join(ws_path, 'data', 'psf', '_model_error')


def err_key(inst_key, filt):
    return f'{inst_key.replace("/", "")}_{psf_models._base_filter(filt).upper()}'


def core_fit(r, t, e, rad, r_fit):
    """Weighted least squares real = a * twin + b within r_fit px: (a, b, chi2, n)."""
    m = (rad <= r_fit) & np.isfinite(r) & (e > 0)
    w = 1.0 / e[m] ** 2
    A = np.column_stack([t[m], np.ones(m.sum())])
    coef = np.linalg.lstsq(A * np.sqrt(w)[:, None], r[m] * np.sqrt(w), rcond=None)[0]
    chi2 = float(np.sum(w * (r[m] - A @ coef) ** 2))
    return float(coef[0]), float(coef[1]), chi2, int(m.sum())


def band(sample, lens, filt, blur_sigma_det=None):
    from astropy.io import fits
    t0 = time.time()
    fld = make_psf.prepare_field(sample, lens, filt, quiet=True)
    if fld is None:
        print(f'{lens} {filt}: no mosaic')
        return
    tbl, sc, inst_key = fld['stars_tbl'], fld['scale'], fld['inst_key']
    rec = dict(n_stars=len(tbl), scale=sc, instrument=inst_key, star_pass=fld['star_pass'])
    if len(tbl) == 0:
        print(f'{lens} {filt}: no PSF stars')
        info_json.update(STARCHECK_JSON, sample, lens, filt, rec)
        return
    x, y = np.asarray(tbl['x'], float), np.asarray(tbl['y'], float)
    coords = [fld['wcs'].pixel_to_world(float(xi), float(yi)) for xi, yi in zip(x, y)]
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            twin, twcs, inj = make_psf_inject.injected_mosaic(sample, lens, filt, coords, blur_sigma_det=blur_sigma_det)
    except KeyError as e:
        # No model source for this filter (STDPSF/focus-diverse have no f475X or f600LP
        # grid, and every such band is empirical anyway): record a deliberate skip so the
        # calibration excludes the band instead of dying.
        if 'unknown filter' not in str(e):
            raise
        rec.update(model=None, skipped=f'no model PSF source: {e}', seconds=round(time.time() - t0, 1))
        info_json.update(STARCHECK_JSON, sample, lens, filt, rec)
        print(f'{lens} {filt}: skipped ({e})')
        return
    # the real star's noise: the star pass's weight map (1/sqrt(WHT) is the per-pixel sigma)
    err = np.sqrt(1.0 / np.clip(fld['wht'], 1e-30, None))
    err[~np.isfinite(fld['wht']) | (fld['wht'] <= 0)] = 0.0
    C = CROP
    yy, xx = np.mgrid[:C, :C]
    rad = np.hypot(xx - C // 2, yy - C // 2)
    r_fit = 3.0 * float(fld['params']['psf_fwhm'])
    B = max(15, int(round(6 * fld['params']['psf_fwhm'])) | 1)
    R, T, E, stars = [], [], [], []
    for xi, yi, c in zip(x, y, coords):
        cx, cy = int(round(xi)), int(round(yi))
        tx, ty = (float(v) for v in twcs.world_to_pixel(c))
        # the twin mosaic is the science lattice up to an integer crop offset: the star's
        # sub-pixel phase must agree, its pixel index may not
        dx, dy = tx - xi, ty - yi
        if max(abs(dx - round(dx)), abs(dy - round(dy))) > 0.05:
            print(f'  WARNING: twin lattice off by ({dx - round(dx):+.3f}, {dy - round(dy):+.3f}) px at star ({xi:.1f}, {yi:.1f})')
        r = make_psf._crop_centered(fld['data'], C, cx, cy)
        t = make_psf._crop_centered(twin, C, cx + int(round(dx)), cy + int(round(dy)))
        e = make_psf._crop_centered(err, C, cx, cy)
        if t.sum() <= 0:
            continue
        a, b, chi2, npx = core_fit(r, t, e, rad, r_fit)
        rb, tb = (make_psf._crop_centered(v, B, C // 2, C // 2) for v in (r - b, t))
        fr = make_psf.measure_fwhm(rb / rb.sum()) if rb.sum() > 0 else None
        ft = make_psf.measure_fwhm(tb / tb.sum())
        ratio = fr / ft if fr and ft else None
        peak_snr = float(np.max((r - b)[rad <= 1.5] / np.maximum(e[rad <= 1.5], 1e-30)))
        stars.append(dict(x=round(float(xi), 2), y=round(float(yi), 2), flux=round(a * float(t.sum()), 4),
                          peak_snr=round(peak_snr, 1), fwhm_real=round(fr * sc, 5) if fr else None,
                          fwhm_twin=round(ft * sc, 5) if ft else None,
                          ratio=round(ratio, 4) if ratio else None,
                          peak_frac_real=round(float((r - b).max() / max((r - b)[rad <= C // 2].sum(), 1e-30)), 4),
                          peak_frac_twin=round(float(t.max() / t[rad <= C // 2].sum()), 4),
                          chi2_per_px=round(chi2 / max(npx - 2, 1), 3),
                          extended=bool(ratio is not None and ratio > EXTENDED_RATIO)))
        R.append(r); T.append(t); E.append(e)
    out = os.path.join(ws_path, 'data', 'psf', stamp_paths.lens_rel(sample, lens), filt)
    os.makedirs(out, exist_ok=True)
    np.savez_compressed(os.path.join(out, 'star_check.npz'), real=np.array(R, np.float32), twin=np.array(T, np.float32),
                        err=np.array(E, np.float32), ratio=np.array([s['ratio'] or np.nan for s in stars]),
                        peak_snr=np.array([s['peak_snr'] for s in stars]),
                        extended=np.array([s['extended'] for s in stars]),
                        scale=sc, blur_sigma_det=inj['blur_sigma_det_px'], r_fit=r_fit,
                        det_scale=make_psf_inject._DRIZ[inst_key]['det_scale'])
    pts = [s['ratio'] for s in stars if s['ratio'] and not s['extended']]
    rec.update(stars=stars, n_extended=sum(s['extended'] for s in stars), n_point=len(pts),
               median_ratio=round(float(np.median(pts)), 4) if pts else None,
               model=inj['method'], blur_sigma_det_px=inj['blur_sigma_det_px'], seconds=round(time.time() - t0, 1))
    info_json.update(STARCHECK_JSON, sample, lens, filt, rec)
    print(f'{lens} {filt}: {len(stars)} stars ({rec["n_extended"]} extended), real/twin FWHM median '
          f'{rec["median_ratio"]} ({inj["method"]}, blur {inj["blur_sigma_det_px"]:.3f} det px), {rec["seconds"]:.0f}s')


def load_checks(snr_min=CAL_MIN_SNR):
    """{inst_filt key: [(r, t, e, scale, det_scale, blur_used, r_fit), ...]} over every band's
    saved crops, point-like stars at or above snr_min."""
    rec = info_json.load(STARCHECK_JSON)
    out = {}
    for sample, lenses in rec.items():
        for lens, filts in lenses.items():
            for filt, r in filts.items():
                if not r or not r.get('stars'):
                    continue
                f = os.path.join(ws_path, 'data', 'psf', stamp_paths.lens_rel(sample, lens), filt, 'star_check.npz')
                if not os.path.isfile(f):
                    continue
                z = np.load(f)
                key = err_key(r['instrument'], filt)
                for i in range(len(z['ratio'])):
                    ratio = z['ratio'][i]
                    if not np.isfinite(ratio) or z['peak_snr'][i] < snr_min or not (1 / EXTENDED_RATIO <= ratio <= EXTENDED_RATIO):
                        continue
                    out.setdefault(key, []).append((z['real'][i].astype(float), z['twin'][i].astype(float),
                                                    z['err'][i].astype(float), float(z['scale']), float(z['det_scale']),
                                                    float(z['blur_sigma_det']), float(z['r_fit'])))
    return out


def _load_cal():
    if os.path.isfile(psf_models.CALIBRATION_JSON):
        with open(psf_models.CALIBRATION_JSON) as f:
            return json.load(f)
    return {}


def calibrate():
    cal = _load_cal()
    cal.setdefault('blur_sigma_det_px', {})
    hist = cal.setdefault('iterations', {})
    fits_ = cal.setdefault('calibration_fit', {})
    for key, stars in sorted(load_checks().items()):
        best = []
        for r, t, e, sc, det, used, r_fit in stars:
            C = r.shape[0]
            yy, xx = np.mgrid[:C, :C]
            rad = np.hypot(xx - C // 2, yy - C // 2)
            chi = np.array([core_fit(r, psf_models.gaussian_blur(t, sg), e, rad, r_fit)[2] for sg in SIGMA_GRID])
            k = int(np.argmin(chi))
            b = SIGMA_GRID[k]
            if 0 < k < len(chi) - 1:
                c0, c1, c2 = chi[k - 1], chi[k], chi[k + 1]
                b += 0.5 * (c0 - c2) / (c0 - 2 * c1 + c2) * (SIGMA_GRID[1] - SIGMA_GRID[0])
            best.append(b * sc / det)          # output px -> detector px
        best = np.array(best)
        med = float(np.median(best))
        used = stars[0][5]
        inst, filt = key.split('_')
        inst_key = {'ACSWFC': 'ACS/WFC', 'WFC3IR': 'WFC3/IR', 'WFC3UVIS': 'WFC3/UVIS', 'WFPC2': 'WFPC2'}[inst]
        h = [it for it in hist.get(key, []) if abs(it['blur_used'] - used) > 1e-6]
        h.append(dict(blur_used=round(used, 4), residual_det_px=round(med, 4), n_stars=len(best)))
        hist[key] = h
        if len(h) >= 2:
            (u0, e0), (u1, e1) = [(it['blur_used'] ** 2, it['residual_det_px'] ** 2) for it in h[-2:]]
            new = float(np.sqrt(max(u1 - e1 * (u1 - u0) / (e1 - e0), 0.0))) if e1 != e0 else float(np.sqrt(u1 + e1))
        else:
            new = float(np.hypot(used, med))
        cal['blur_sigma_det_px'].setdefault(inst_key, {})[filt] = round(new, 4)
        fits_[key] = dict(n_stars=len(best), residual_det_px=round(med, 4),
                          residual_16_84_det_px=[round(float(v), 4) for v in np.percentile(best, [16, 84])],
                          blur_used=round(used, 4), blur_next=round(new, 4))
        print(f'{key}: {len(best)} stars, residual blur {med:.3f} detector px (16-84% '
              f'{np.percentile(best, 16):.3f}-{np.percentile(best, 84):.3f}) at {used:.3f} -> {new:.3f}')
    cal['_how'] = ('psf_star_check.py --calibrate: each real PSF star against its noiseless twin injected at the '
                   f'same position into the same frames (point-like stars, peak S/N >= {CAL_MIN_SNR:g}); per star the '
                   'extra Gaussian blur of the twin that minimises the core chi^2, the median over the stars, in '
                   'detector px; iterated (secant) because a blur on the resampled star stands in for one on the detector.')
    with open(psf_models.CALIBRATION_JSON, 'w') as f:
        json.dump(cal, f, indent=4)
    print(f'wrote {psf_models.CALIBRATION_JSON}')


def error_maps():
    from astropy.io import fits
    from astropy.stats import sigma_clip
    cal = _load_cal()
    os.makedirs(MODEL_ERROR_DIR, exist_ok=True)
    errs = cal.setdefault('model_error', {})
    for key, stars in sorted(load_checks().items()):
        res, var = [], []
        for r, t, e, sc, det, used, r_fit in stars:
            C = r.shape[0]
            yy, xx = np.mgrid[:C, :C]
            rad = np.hypot(xx - C // 2, yy - C // 2)
            a, b, _, _ = core_fit(r, t, e, rad, r_fit)
            tot = t.sum()
            res.append(((r - b) / a - t) / tot)
            var.append((e / (a * tot)) ** 2)
        res, var = np.array(res), np.array(var)
        sq = sigma_clip(res ** 2 - var, sigma=3.0, axis=0)
        em = np.sqrt(np.clip(np.ma.mean(sq, axis=0).filled(0.0), 0.0, None))
        hdr = fits.Header()
        hdr['PSFEKEY'] = (key, 'camera_filter')
        hdr['PSFERR'] = ('model_vs_stars', 'rms real-star residual about its injected twin')
        hdr['PSFNSTAR'] = (len(stars), 'stars (point-like, peak S/N >= %g)' % CAL_MIN_SNR)
        hdr['PSFBLUR'] = (stars[0][5], 'blur the twins carried (detector px)')
        fits.PrimaryHDU(em.astype(np.float32), hdr).writeto(os.path.join(MODEL_ERROR_DIR, f'{key}.fits'), overwrite=True)
        errs[key] = dict(n_stars=len(stars), err_frac=round(float(np.sqrt(np.sum(em ** 2))), 5))
        print(f'{key}: {len(stars)} stars, model error fraction {errs[key]["err_frac"]:.4f}')
    with open(psf_models.CALIBRATION_JSON, 'w') as f:
        json.dump(cal, f, indent=4)
    print(f'wrote {psf_models.CALIBRATION_JSON} and {MODEL_ERROR_DIR}/<INST>_<FILT>.fits')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--lens')
    p.add_argument('--filt')
    p.add_argument('--sample', default=mast_target_names.DEFAULT_SAMPLE)
    p.add_argument('--blur', type=float, default=None,
                   help='extra blur of the twins in detector px (default: the calibrated one; 0 for none)')
    p.add_argument('--calibrate', action='store_true')
    p.add_argument('--error-maps', action='store_true')
    a = p.parse_args()
    warnings.simplefilter('ignore')
    if a.calibrate:
        return calibrate()
    if a.error_maps:
        return error_maps()
    if not (a.lens and a.filt):
        p.error('give --lens and --filt, or --calibrate / --error-maps')
    band(a.sample, a.lens, a.filt, a.blur)


if __name__ == '__main__':
    sys.exit(main())
