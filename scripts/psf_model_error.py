"""Calibrate the STDPSF model-tier PSF error against empirical drizzled truth.

READ-ONLY QC: reads the kernels already in `data/psf/`, rebuilds the analytic model
in memory, and prints a comparison. Writes nothing under `data/` and only the
optional `--out` JSON elsewhere. It never touches a canonical product.

**Why this exists.** The model tier's only uncertainty story is an ensemble scatter
(bootstrap over stars for WFPC2's MAST-DB tier, leave-one-exposure-out for ACS
focus-diverse) -- and STDPSF has no natural ensemble at all, so it carries none. But
ensemble scatter is the *standard error of the mean ePSF*, which shrinks as
1/sqrt(N), while the model's error against a specific observation is set by
focus/breathing mismatch and does not shrink at all. Measuring the model against an
empirical build on the SAME lens+filter -- where the empirical ePSF is the true
drizzled PSF of that field -- gives the model-vs-truth quantity a likelihood actually
wants. Measured 2026-08-02: the model error runs a median 2x the bootstrap scatter on
well-constrained builds and up to 7x. See memory: stdpsf_model_error_measured.

Targets default to every empirical build whose model tier would be STDPSF: F160W
(WFC3/IR, any sample) and all of `gallery` (WFC3/UVIS). ACS/WFC is excluded because
its model tier is the focus-diverse ePSF, not STDPSF.

The metric of record is `resid` = sqrt(sum((mod-emp)^2)) on aligned, unit-sum kernels
-- the same construction as `psf_err_frac` in info/lens_psf.json, so the two are
directly comparable. `corrected` removes the empirical truth's own bootstrap noise in
quadrature, since a star-poor empirical build is itself uncertain.

**Three traps this script exists to get right:**

 1. `make_psf.measure_fwhm`'s 2D-Gaussian fit is unreliable on sharp, near-critically-
    sampled kernels -- it returns a LARGER FWHM for a HIGHER-peak kernel on the
    unbroadened analytic model. `fwhm_radial()` below uses a half-max crossing on the
    azimuthal radial profile instead. (This is not a live pipeline bug: the pipeline
    only ever fits the broader drizzled/injected kernels, and every recorded
    `fwhm_pix` was checked and is sane.)
 2. Raw residuals are swamped by centring: model-vs-empirical centroid offsets reach
    ~1.9px. The model is shifted sub-pixel to minimise the residual before it is
    measured, so `resid` is a SHAPE error and the fitted shift is reported separately.
 3. The L2 residual is 94-97% core-dominated, so it barely sees the wing deficit --
    the model's enclosed wing flux runs 6-48% low. `wing%` is reported separately
    because for lens modelling the extended arc flux is where that bites.

Note this script broadens the analytic model with `psf_models.analytic_drop_broaden`
as a stand-in for the real injection promotion. That is right in the MEDIAN but
unreliable per-lens (individual lenses moved up to 2.7x, both directions, when checked
against real injection). Quote group medians from this script, not per-lens values;
for a per-lens number use psf_model_error_injected.py against a real injected build.

    uv run python scripts/psf_model_error.py
    uv run python scripts/psf_model_error.py --sample gallery --filt f606W
    uv run python scripts/psf_model_error.py --out /tmp/psf_model_error.json
"""
import argparse
import glob
import json
import os
import warnings

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from scipy.ndimage import shift as ndshift
from scipy.optimize import minimize

import make_psf
import psf_models

ws_path = '/Users/samlange/Code/HST_lenses_data_reduction'


def unit(k):
    """Unit-sum copy of `k` (a no-op guard on an all-zero kernel)."""
    k = np.asarray(k, dtype=float)
    s = k.sum()
    return k / s if s else k


def radial_profile(k):
    """Azimuthally-averaged profile: (radii, mean value) in 1px annuli about the centre."""
    ny, nx = k.shape
    cy, cx = (ny - 1) / 2.0, (nx - 1) / 2.0
    yy, xx = np.mgrid[:ny, :nx]
    r = np.hypot(xx - cx, yy - cy).ravel()
    v = k.ravel()
    cent, prof = [], []
    for i in range(int(np.floor(min(cy, cx)))):
        m = (r >= i) & (r < i + 1)
        if m.any():
            cent.append(i + 0.5)
            prof.append(float(v[m].mean()))
    return np.array(cent), np.array(prof)


def fwhm_radial(k):
    """FWHM (px) from the half-max crossing of the radial profile; NaN if unresolved.

    Robust where make_psf.measure_fwhm's 2D-Gaussian fit is not -- see trap 1 in the
    module docstring. Returns NaN when the half-max radius falls inside the first
    annulus (the kernel is too sharp for this estimator to resolve at 1px binning).
    """
    rr, pp = radial_profile(k)
    half = k.max() / 2.0
    below = np.where(pp < half)[0]
    if not len(below) or below[0] == 0:
        return float('nan')
    i = below[0]
    r_half = rr[i - 1] + (half - pp[i - 1]) * (rr[i] - rr[i - 1]) / (pp[i] - pp[i - 1])
    return float(2.0 * r_half)


def align_resid(mod, emp):
    """Shift `mod` sub-pixel to minimise ||mod-emp||_2. Returns (resid, dy, dx, aligned).

    Without this the residual measures centring, not shape -- see trap 2 above.
    """
    def cost(p):
        m = unit(np.clip(ndshift(mod, p, order=3, mode='constant', cval=0.0), 0, None))
        return float(((m - emp) ** 2).sum())

    best = minimize(cost, [0.0, 0.0], method='Nelder-Mead',
                    options=dict(xatol=1e-3, fatol=1e-12, maxiter=400))
    dy, dx = best.x
    aligned = unit(np.clip(ndshift(mod, [dy, dx], order=3, mode='constant', cval=0.0),
                           0, None))
    return float(np.sqrt(max(best.fun, 0.0))), float(dy), float(dx), aligned


def split_residuals(mod_aligned, emp, core_radius=3.0):
    """(resid_core, resid_wing, wing_flux_fraction) about the kernel centre."""
    ny, nx = emp.shape
    cy, cx = (ny - 1) / 2.0, (nx - 1) / 2.0
    yy, xx = np.mgrid[:ny, :nx]
    r = np.hypot(xx - cx, yy - cy)
    core, wing = r < core_radius, r >= core_radius
    d = mod_aligned - emp
    rc = float(np.sqrt((d[core] ** 2).sum()))
    rw = float(np.sqrt((d[wing] ** 2).sum()))
    denom = emp[wing].sum()
    wf = float((mod_aligned[wing].sum() - denom) / denom) if denom else float('nan')
    return rc, rw, wf


def is_stdpsf_tier(sample, filt):
    """True where the model fallback for this product would be STDPSF.

    F160W (WFC3/IR) and every gallery band (WFC3/UVIS). ACS/WFC goes to the
    focus-diverse ePSF and WFPC2 F606W to the MAST PSF DB, so neither is STDPSF.
    """
    return filt.lower() == 'f160w' or sample == 'gallery'


def collect_targets(samples=None, filts=None):
    """Empirical builds on the STDPSF tier, as (sample, lens, filt, json_entry)."""
    psf_json = json.load(open(os.path.join(ws_path, 'info', 'lens_psf.json')))
    out = []
    for sample, lenses in psf_json.items():
        if samples and sample not in samples:
            continue
        for lens, per_filt in (lenses or {}).items():
            for filt, entry in (per_filt or {}).items():
                if not entry or entry.get('method') != 'empirical':
                    continue
                if filts and filt not in filts:
                    continue
                if is_stdpsf_tier(sample, filt):
                    out.append((sample, lens, filt, entry))
    return sorted(out)


def build_model_kernel(sample, lens, filt, star_size):
    """Rebuild the analytic STDPSF kernel + its drizzle-broadened form, in memory.

    Mirrors make_psf.py's model branch exactly (same oversample/out_scale, rotated to
    North-up through the same representative_input_cd), then applies the analytic
    stand-in for the injection promotion. Returns (sharp, broadened, meta).
    """
    scis = sorted(glob.glob(os.path.join(ws_path, 'data', 'drizzled', sample, lens,
                                         filt, '*_sci.fits')))
    if not scis:
        raise FileNotFoundError(f'no drizzled sci for {sample}/{lens}/{filt}')
    # Either pass carries the same drizzle geometry keywords; prefer no-CR when present.
    sci = next((s for s in scis if 'nocrrej' in os.path.basename(s)), scis[0])
    hdr = fits.getheader(sci)

    inst_key = make_psf.instrument_key(hdr)
    scale = proj_plane_pixel_scales(WCS(hdr).celestial)[0] * 3600.0
    oversample = int(make_psf.resolve_params(inst_key, {}, None)['oversample'])

    cd = make_psf.representative_input_cd(inst_key, lens, filt, sample)
    eps = psf_models.model_psf(inst_key, filt, oversample=oversample, size=star_size,
                               out_scale=scale, cd_detector=cd)
    sharp = make_psf.oversampled_to_kernel(eps, oversample, star_size)
    sharp, _ = make_psf.subtract_pedestal(sharp)

    # NDRIZIM counts drizzled IMAGES; ACS/WFC and WFC3/UVIS FLCs are 2-chip MEFs whose
    # chips share a dither phase, so halve them -- same correction make_psf_inject uses.
    n_chips = 2 if inst_key in ('ACS/WFC', 'WFC3/UVIS') else 1
    n_frames = max(1, int(hdr.get('NDRIZIM', n_chips)) // n_chips)
    broad, box = psf_models.analytic_drop_broaden(
        sharp, float(hdr['D001PIXF']), float(hdr['D001ISCL']), float(hdr['D001SCAL']),
        n_frames=n_frames)
    broad, _ = make_psf.subtract_pedestal(broad)
    return unit(sharp), unit(broad), dict(inst=inst_key, scale=scale,
                                          n_frames=n_frames, box=box, cd=cd is not None)


def summarise(rows):
    """Print per-group medians. Group medians are the trustworthy output of this
    script -- per-lens analytic-broadened values are not (see the module docstring)."""
    print()
    print(f'{"group":22s}{"n":>3s}{"resid med":>11s}{"range":>19s}'
          f'{"boot med":>10s}{"corrected":>11s}{"wing% med":>11s}')
    keys = []
    for r in rows:
        k = 'F160W (WFC3/IR)' if r['filt'].lower() == 'f160w' else \
            f'{r["sample"]} {r["filt"]}'
        if k not in keys:
            keys.append(k)
    for k in keys:
        g = [r for r in rows
             if (('F160W' in k and r['filt'].lower() == 'f160w') or
                 (k == f'{r["sample"]} {r["filt"]}' and r['filt'].lower() != 'f160w'))]
        d = np.array([r['resid'] for r in g])
        b = np.array([r['boot'] or 0.0 for r in g])
        c = np.sqrt(np.clip(d * d - b * b, 0, None))
        w = 100 * np.array([r['wing_frac'] for r in g])
        print(f'{k:22s}{len(g):3d}{np.median(d):11.4f}'
              f'{f"{d.min():.4f}-{d.max():.4f}":>19s}'
              f'{np.median(b):10.4f}{np.median(c):11.4f}{np.median(w):+11.1f}')

    good = [r for r in rows if (r['n_stars'] or 0) >= 6]
    if good:
        d = np.median([r['resid'] for r in good])
        b = np.median([r['boot'] or 0.0 for r in good])
        print(f'\nwell-constrained truth (n_stars>=6, n={len(good)}): '
              f'resid={d:.4f}  boot={b:.4f}  ratio={d / b if b else float("nan"):.1f}x')
    frac = np.median([r['resid_core'] / r['resid'] for r in rows if r['resid']])
    print(f'core share of residual: {frac:.3f}  '
          f'(so the L2 number barely sees the wing deficit)')


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--sample', action='append',
                   help='restrict to this sample (repeatable; default: all)')
    p.add_argument('--filt', action='append',
                   help='restrict to this filter (repeatable; default: all STDPSF-tier)')
    p.add_argument('--core-radius', type=float, default=3.0,
                   help='core/wing split radius in pixels (default 3)')
    p.add_argument('--out', help='also write the per-product rows to this JSON path')
    a = p.parse_args()

    targets = collect_targets(a.sample, a.filt)
    if not targets:
        print('no matching empirical STDPSF-tier products')
        return

    print(f'{"sample":12s}{"lens":12s}{"filt":7s}{"nst":>4s}{"resid":>8s}{"boot":>8s}'
          f'{"corr":>8s}{"wing%":>8s}{"shift(y,x)":>16s}')
    rows = []
    for sample, lens, filt, entry in targets:
        kpath = os.path.join(ws_path, 'data', 'psf', sample, lens, filt,
                             'psf_kernel.fits')
        if not os.path.isfile(kpath):
            print(f'  SKIP {sample}/{lens}/{filt}: no psf_kernel.fits')
            continue
        emp = unit(fits.getdata(kpath).astype(float))
        try:
            sharp, broad, meta = build_model_kernel(sample, lens, filt, emp.shape[0])
        except Exception as exc:
            print(f'  SKIP {sample}/{lens}/{filt}: {exc}')
            continue

        resid, dy, dx, aligned = align_resid(broad, emp)
        resid_sharp, _, _, _ = align_resid(sharp, emp)
        rc, rw, wf = split_residuals(aligned, emp, a.core_radius)
        boot = entry.get('psf_err_frac') or 0.0
        corr = float(np.sqrt(max(resid ** 2 - boot ** 2, 0.0)))

        rows.append(dict(sample=sample, lens=lens, filt=filt, inst=meta['inst'],
                         n_stars=entry.get('n_stars'), boot=boot, resid=resid,
                         corrected=corr, resid_sharp=resid_sharp, resid_core=rc,
                         resid_wing=rw, wing_frac=wf, dy=dy, dx=dx,
                         fwhm_emp=fwhm_radial(emp), fwhm_model=fwhm_radial(broad),
                         n_frames=meta['n_frames'], box=meta['box'],
                         rotated=meta['cd'], scale=meta['scale']))
        print(f'{sample:12s}{lens:12s}{filt:7s}{str(entry.get("n_stars")):>4s}'
              f'{resid:8.4f}{boot:8.4f}{corr:8.4f}{100 * wf:+8.1f}'
              f'{f"({dy:+.2f},{dx:+.2f})":>16s}')

    if rows:
        summarise(rows)
    if a.out and rows:
        with open(a.out, 'w') as fh:
            json.dump(rows, fh, indent=1)
        print(f'\nwrote {a.out}  ({len(rows)} products)')


if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    main()
