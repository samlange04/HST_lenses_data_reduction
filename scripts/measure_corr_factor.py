"""
Measure the noise-map inflation factor for every cutout -> info/noise_corr_factors.json.

Drizzle correlates neighbouring output pixels but conserves noise: for white input noise
the sigma-normalised ACF has C(0) < 1, positive neighbours, and sum C = 1 over the drop
width, so 1/sqrt(WHT) is already right for any patch >= the drop and drizzle itself needs
no factor. What the measured sum C > 1 catches is upstream of drizzle: an ERR-array deficit
(~10% rms in ACS/UVIS FLCs, white) plus a few percent of faint unmasked structure. See
AGENTS.md *Drizzle correlated noise* and local/literature/drizzle_noise.md. This measures,
per stamp, on blank sky in its own mosaic:

  * the ACF of the sigma-normalised sky (acf_ratio) -- THE factor: sqrt(sum C) over the
    drop window |dx|,|dy| <= round(pixfrac * native/output) px (the smallest window over
    which drizzle's redistribution has summed back to 1), less the far floor. This is what
    make_cutouts.py applies.
  * the blank-sky block-sum ratio at 0.05-2" (block_ratio) -- a diagnostic only. It never
    plateaus: past the kernel a weak, environment-dependent long-range floor (faint
    wings, sky/flat residuals) keeps adding up, so it has no single right answer. See
    AGENTS.md *Drizzle correlated noise* before changing which quantity is used.

Per stamp:
  * The stamp's own mosaic (same pass, same reduction) is checked to still reproduce the
    stamp's sci pixels exactly; one that doesn't is measured anyway and flagged `stale`.
  * Region: +-HALF arcsec about the stamp centre, minus r < LENS_R, minus detected sources
    (photutils segmentation, dilated), minus WHT <= 0. Run at two detection thresholds;
    3 sigma is the one used, 1.5 sigma is the sensitivity check.

Writes every measurement to diagnostics/corr_factor/ and prints the per-(instrument,
filter) medians in the shape of info/noise_corr_factors.json (which it does NOT write --
copy the values in by hand after reading them).

    uv run python scripts/measure_corr_factor.py [--sample S ...] [--jobs 8]
"""
import argparse
import glob
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.convolution import Gaussian2DKernel, convolve
from scipy import ndimage

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cutout_paths                                  # noqa: E402
from make_cutouts import find_products              # noqa: E402

WS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HALF = 30.0                              # arcsec half-width of the measurement region
LENS_R = 6.0                             # arcsec radius around the lens left out
BLOCKS_ARCSEC = (0.05, 0.25, 0.5, 1.0, 1.5, 2.0)
THRESHOLDS = (1.5, 3.0)                  # detection thresholds (x background rms)
N_BOOT = 300


def robust_std(x):
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def source_mask(img, good, nsig):
    from photutils.background import Background2D, MedianBackground
    from photutils.segmentation import detect_sources
    data = np.where(good, img, 0.0)
    bkg = Background2D(data, 64, filter_size=3, mask=~good,
                       bkg_estimator=MedianBackground())
    sm = convolve(data - bkg.background, Gaussian2DKernel(1.5), mask=~good)
    seg = detect_sources(sm, nsig * bkg.background_rms, npixels=8, mask=~good)
    if seg is None:
        return np.zeros_like(good)
    return ndimage.binary_dilation(seg.data > 0, iterations=7)


def block_ratio(img, var, usable, n, rng):
    ny, nx = img.shape[0] // n, img.shape[1] // n
    if ny < 1 or nx < 3:
        return None
    sl = (slice(0, ny * n), slice(0, nx * n))
    S = img[sl].reshape(ny, n, nx, n).sum(axis=(1, 3))
    V = var[sl].reshape(ny, n, nx, n).sum(axis=(1, 3))
    ok = usable[sl].reshape(ny, n, nx, n).all(axis=(1, 3))
    # Pairs one block apart horizontally (a gap of one block width, so the pair shares no
    # drizzle-correlated pixels across their touching edge).
    pair = ok[:, :-2] & ok[:, 2:]
    z = ((S[:, :-2] - S[:, 2:]) / np.sqrt(V[:, :-2] + V[:, 2:]))[pair]
    if z.size < 20:
        return None
    r = robust_std(z)
    boot = [robust_std(rng.choice(z, z.size)) for _ in range(N_BOOT)]
    return dict(ratio=float(r), err=float(np.std(boot)), npairs=int(z.size))


ACF_MAXLAG = 5                            # kernel windows L = 0..ACF_MAXLAG px
FLOOR_LAGS = (7, 12)                      # Chebyshev-lag annulus that defines the floor
HIPASS = 51                               # px median filter removed before the ACF


def acf_ratio(img, var, usable):
    """Noise-map inflation factor from the sigma-normalised blank-sky ACF.

    C(d) = <u(x) u(x+d)> with u = (sky - local median) / sigma_map, so C(0) is the
    per-pixel (emp/map)^2 and the integrated factor over a window |dx|,|dy| <= L is
    R(L) = sqrt(sum C). Three parts: drizzle redistribution (C(0) < 1, neighbours > 0,
    summing to exactly 1 by the drop width -- Casertano's R is 1/sqrt of its C(0)); a white
    ERR deficit that scales every lag alike; and a weak long-range floor (~0.01-0.05 per
    lag from faint sources / sky and flat residuals) that never stops adding up -- which is
    why a block-sum ratio does not plateau. The floor is the mean C over the FLOOR_LAGS
    annulus; R_sub(L) subtracts it from every lag in the window. R_sub at the drop window
    is the --corr-factor quantity."""
    from numpy.fft import rfft2, irfft2
    filled = np.where(usable, img, np.median(img[usable]))
    u = np.where(usable, (img - ndimage.median_filter(filled, size=HIPASS)) /
                 np.sqrt(np.where(usable, var, 1.0)), 0.0)
    m = usable.astype(np.float64)
    sh = [2 * n for n in u.shape]
    fu, fm = rfft2(u, sh), rfft2(m, sh)
    cu, cm = irfft2(fu * np.conj(fu), sh), irfft2(fm * np.conj(fm), sh)
    L = FLOOR_LAGS[1]
    idx = np.r_[-L:0, 0:L + 1]
    C = cu[np.ix_(idx, idx)] / np.maximum(cm[np.ix_(idx, idx)], 1)   # lag 0 at [L, L]
    cheb = np.maximum(*np.abs(np.mgrid[-L:L + 1, -L:L + 1]))
    floor = float(C[(cheb >= FLOOR_LAGS[0]) & (cheb <= FLOOR_LAGS[1])].mean())
    raw, sub = {}, {}
    for l in range(ACF_MAXLAG + 1):
        w = cheb <= l
        raw[str(l)] = float(np.sqrt(max(C[w].sum(), 0)))
        sub[str(l)] = float(np.sqrt(max((C[w] - floor).sum(), 0)))
    return dict(c0=float(C[L, L]), c1=float((C[L, L + 1] + C[L + 1, L]) / 2),
                floor=floor, R_raw=raw, R_sub=sub)


def measure(job):
    sample, lens, filt, stamp = job
    out = dict(sample=sample, lens=lens, filt=filt, stamp=os.path.basename(stamp))
    try:
        prefix = 'cutout_cr' if os.path.basename(stamp).startswith('cutout_cr') else 'cutout'
        variant = cutout_paths.stamp_variant(stamp)
        ddir = os.path.join(cutout_paths.drizzled_root(WS, variant), sample, lens, filt)
        sci_f, wht_f = find_products(ddir, 'cr' if prefix == 'cutout_cr' else 'nocrrej')
        stamp_sci = fits.getdata(stamp).astype(np.float64)
        nh = fits.getheader(stamp.replace('_sci.fits', '_noise.fits'))
        k = float(nh.get('NOISEK', 1.0))
        sh = fits.getheader(stamp)
        with fits.open(sci_f, memmap=True) as h:
            hdr = h[0].header
            wcs = WCS(hdr)
            # Stamp pixel (0,0) in the mosaic, via the sky (Cutout2D keeps whole pixels).
            sky0 = WCS(sh).celestial.pixel_to_world(0, 0)
            x0, y0 = (int(round(float(v))) for v in wcs.world_to_pixel(sky0))
            ny, nx = stamp_sci.shape
            mos_patch = np.asarray(h[0].data[y0:y0 + ny, x0:x0 + nx], dtype=np.float64)
            out['stale'] = not (mos_patch.shape == stamp_sci.shape and
                                np.array_equal(mos_patch, stamp_sci, equal_nan=True))
            scale = proj_plane_pixel_scales(wcs.celestial)[0] * 3600.0
            cx, cy = x0 + (nx - 1) / 2.0, y0 + (ny - 1) / 2.0
            hw = int(HALF / scale)
            Y0, Y1 = max(0, int(cy) - hw), min(hdr['NAXIS2'], int(cy) + hw)
            X0, X1 = max(0, int(cx) - hw), min(hdr['NAXIS1'], int(cx) + hw)
            img = np.asarray(h[0].data[Y0:Y1, X0:X1], dtype=np.float64)
        with fits.open(wht_f, memmap=True) as h:
            wht = np.asarray(h[0].data[Y0:Y1, X0:X1], dtype=np.float64)
        good = np.isfinite(img) & np.isfinite(wht) & (wht > 0)
        yy, xx = np.mgrid[Y0:Y1, X0:X1]
        good &= np.hypot(xx - cx, yy - cy) * scale > LENS_R
        var = np.where(good, k * k / np.where(good, wht, 1.0), 0.0)
        img = np.where(good, img, 0.0)
        out.update(scale=float(scale), instrume=str(hdr.get('INSTRUME', '')).strip(),
                   detector=str(hdr.get('DETECTOR', '')).strip(),
                   casr=nh.get('CASR'), ndrizim=hdr.get('NDRIZIM'),
                   pixfrac=hdr.get('D001PIXF'), native_scale=hdr.get('D001ISCL'),
                   region_arcsec=[(X1 - X0) * scale,
                                                               (Y1 - Y0) * scale])
        rng = np.random.default_rng(0)
        for nsig in THRESHOLDS:
            usable = good & ~source_mask(img, good, nsig)
            res = {}
            for b in BLOCKS_ARCSEC:
                n = max(1, int(round(b / scale)))
                r = block_ratio(img, var, usable, n, rng)
                if r is not None:
                    res[f'{n * scale:.2f}'] = dict(npix=n, **r)
            out[f'thr{nsig:g}'] = dict(frac_usable=float(usable.sum() / good.sum()),
                                       blocks=res, acf=acf_ratio(img, var, usable))
    except Exception as e:                                   # noqa: BLE001
        out['error'] = f'{type(e).__name__}: {e}'
    return out


def kernel_window(r):
    """Drizzle kernel half-width in output px: the drop is pixfrac * native/output px wide,
    so beyond round(that) lag the kernel overlap is negligible (ACS/UVIS 1, WFPC2/F160W 2)."""
    return max(1, int(round(float(r['pixfrac']) * float(r['native_scale']) / r['scale'])))


def summarise(results, thr='thr3'):
    groups = {}
    for r in results:
        if 'error' in r:
            continue
        inst = 'WFPC2' if r['instrume'] == 'WFPC2' else f"{r['instrume']}/{r['detector']}"
        groups.setdefault((inst, r['filt'].split('_')[0]), []).append(r)
    print(f'\nper-(instrument, filter) factor, {thr} mask, kernel-window ACF (median, 16-84%):')
    for (inst, band), rs in sorted(groups.items()):
        v = [r[thr]['acf']['R_sub'][str(kernel_window(r))] for r in rs]
        print(f'  {inst:10s} {band:6s} n={len(rs):3d}  L={kernel_window(rs[0])}  '
              f'{max(1.0, np.median(v)):.2f}  [{np.percentile(v, 16):.2f}-'
              f'{np.percentile(v, 84):.2f}]' + ('  (floored at 1.0)' if np.median(v) < 1 else ''))


def jobs_for(samples):
    root = cutout_paths.cutouts_root(WS)
    for sample in samples:
        for d in sorted(glob.glob(os.path.join(root, sample, '*', '*'))):
            for prefix in cutout_paths.PREFIXES:
                stamps = list(cutout_paths.list_stamps(d, prefix).values())
                if stamps:
                    lens, filt = d.split(os.sep)[-2:]
                    yield sample, lens, filt, stamps[0]
                    break


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--sample', nargs='+', default=['slacs_gold', 'slacs_other', 'gallery'])
    p.add_argument('--lens', default=None, help='restrict to one lens (testing)')
    p.add_argument('--jobs', type=int, default=6)
    p.add_argument('--out', default=os.path.join(WS, 'diagnostics', 'corr_factor',
                                                 'corr_factor_measurements.json'))
    a = p.parse_args()
    jobs = [j for j in jobs_for(a.sample) if a.lens is None or j[1] == a.lens]
    print(f'{len(jobs)} stamps')
    with ProcessPoolExecutor(a.jobs) as ex:
        results = list(ex.map(measure, jobs))
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, 'w') as f:
        json.dump(results, f, indent=1)
    for r in results:
        if 'error' in r:
            print(f"ERROR {r['sample']} {r['lens']} {r['filt']}: {r['error']}")
    print(f'wrote {a.out}')
    summarise(results)


if __name__ == '__main__':
    main()
