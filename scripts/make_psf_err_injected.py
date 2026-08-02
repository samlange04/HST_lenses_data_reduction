"""Per-pixel error maps for the drizzle-broadened INJECTED kernels.

The injected kernel is canonical for every model-tier product (see
make_psf_inject.py), but until now it carried no uncertainty of its own: the
model-tier ensemble maps describe the pre-broadening analytic kernel and are parked
under `*_analytic_err.fits`, and the STDPSF tier had nothing at all. This writes
`psf_kernel_err.fits` (full) and `cutout_[cr_]psf_err.fits` (trimmed) describing the
kernel that is actually canonical, in the same single-HDU, non-renormalised,
`sqrt(sum(err**2))`-scalar convention make_psf.py uses for the empirical tier -- so
any `[0]` / `al.Kernel2D` reader keeps working.

**Two sources, because the tiers differ in what is measurable. The header and
info/lens_psf.json always say which one a given map came from -- they are NOT
interchangeable.**

`ensemble` (WFPC2 MAST-DB, ACS focus-diverse) -- the existing analytic ensemble map
    propagated through the same drizzle broadening the kernel got. This is a LOWER
    BOUND, flagged `PSFEBND=True`: an ensemble bootstrap estimates the standard error
    of the MEAN ePSF, which shrinks as 1/sqrt(N), while the model's error against a
    specific observation is set by focus/breathing mismatch and does not shrink at
    all. It is the only thing measurable for WFPC2 F606W, whose fields are too
    star-poor for any empirical build to exist to compare against.
    Propagation assumes the ensemble perturbations are spatially correlated on scales
    >= the drop box, so the std transforms like the signal (box (x) sigma) rather than
    in quadrature. That holds for focus/star-sampling perturbations, which are smooth.

`calibrated` (STDPSF: F160W, gallery UVIS) -- built from the measured model-vs-truth
    residual, the quantity `ensemble` cannot see. `scripts/psf_model_error.py` /
    `psf_model_error_injected.py` measure the model against empirical builds on OTHER
    lenses of the same instrument+filter (the model-tier lenses have no empirical
    build themselves -- that is why they fell back to the model), so this is
    necessarily a group-level transfer by filter, not a per-lens measurement. The
    RADIAL SHAPE comes from the measured residual profile and the MAGNITUDE from the
    group's median quadrature-corrected residual; a flat map would be badly wrong,
    since the residual is ~95% core-dominated while the wings carry a 6-48% flux
    deficit. Measured 2026-08-02: ~4x the bootstrap scatter it replaces.

Empirical-tier products are never touched -- their `psf_kernel_err.fits` already
describes their canonical kernel. Only products whose info/lens_psf.json `method`
starts with `inject` are considered.

    uv run python scripts/make_psf_err_injected.py --dry-run
    uv run python scripts/make_psf_err_injected.py
    uv run python scripts/make_psf_err_injected.py --sample gallery --filt f606W
"""
import argparse
import glob
import os
import warnings

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

import info_json
import make_psf
import psf_models
from psf_model_error import (align_resid, build_model_kernel, collect_targets,
                             split_residuals, unit, ws_path)

PSF_JSON = os.path.join(ws_path, 'info', 'lens_psf.json')


# ── Target selection ────────────────────────────────────────────────────────────
def injected_products(samples=None, filts=None):
    """(sample, lens, filt, entry) for every product whose canonical kernel is injected."""
    data = info_json.load(PSF_JSON)
    out = []
    for sample, lenses in data.items():
        if samples and sample not in samples:
            continue
        for lens, per_filt in (lenses or {}).items():
            for filt, entry in (per_filt or {}).items():
                if not entry or not str(entry.get('method', '')).startswith('inject'):
                    continue
                if filts and filt not in filts:
                    continue
                out.append((sample, lens, filt, entry))
    return sorted(out)


def calibration_group(sample, filt):
    """Group key a product's calibrated error map is transferred from."""
    return 'f160W' if filt.lower() == 'f160w' else f'{sample}:{filt}'


def cutout_prefix(cutouts_dir):
    """'cutout_cr' where a CR-pass kernel exists, else 'cutout' (F160W has no CR pass)."""
    return ('cutout_cr'
            if os.path.isfile(os.path.join(cutouts_dir, 'cutout_cr_psf.fits'))
            else 'cutout')


# ── Source 1: propagate the analytic ensemble through the broadening ────────────
def drizzle_box(sample, lens, filt):
    """(box_width_pix, n_frames) of the drop the science drizzle projects onto the grid."""
    scis = sorted(glob.glob(os.path.join(ws_path, 'data', 'drizzled', sample, lens,
                                         filt, '*_sci.fits')))
    if not scis:
        raise FileNotFoundError(f'no drizzled sci for {sample}/{lens}/{filt}')
    hdr = fits.getheader(next((s for s in scis if 'nocrrej' in os.path.basename(s)),
                              scis[0]))
    inst_key = make_psf.instrument_key(hdr)
    n_chips = 2 if inst_key in ('ACS/WFC', 'WFC3/UVIS') else 1
    n_frames = max(1, int(hdr.get('NDRIZIM', n_chips)) // n_chips)
    box = (float(hdr['D001PIXF']) * float(hdr['D001ISCL']) / float(hdr['D001SCAL'])
           ) / np.sqrt(n_frames)
    return float(box), n_frames


def ensemble_err(psf_dir, sample, lens, filt):
    """Broaden the analytic ensemble error map onto the injected kernel's resolution.

    Not renormalised: drop_convolve_box preserves the sum (the box transfer is 1 at
    zero frequency), which is exactly the renormalisation the broadened kernel gets,
    so the map stays in the kernel's own amplitude units like every other *_err.fits.
    """
    src = os.path.join(psf_dir, 'psf_kernel_analytic_err.fits')
    if not os.path.isfile(src):
        return None, None
    err = np.asarray(fits.getdata(src), dtype=float)
    box, n_frames = drizzle_box(sample, lens, filt)
    out = np.clip(psf_models.drop_convolve_box(err, box), 0.0, None)
    return out, dict(box=box, n_frames=n_frames)


# ── Source 2: calibrated model-vs-empirical residual ────────────────────────────
def _residual_profile(emp, mod):
    """(radii, per-pixel RMS residual by radius) after sub-pixel alignment."""
    _, _, _, aligned = align_resid(mod, emp)
    d = aligned - emp
    ny, nx = d.shape
    cy, cx = (ny - 1) / 2.0, (nx - 1) / 2.0
    yy, xx = np.mgrid[:ny, :nx]
    r = np.hypot(xx - cx, yy - cy).ravel()
    dd = (d ** 2).ravel()
    cent, prof = [], []
    for i in range(int(np.floor(min(cy, cx)))):
        m = (r >= i) & (r < i + 1)
        if m.any():
            cent.append(i + 0.5)
            prof.append(float(np.sqrt(dd[m].mean())))
    return np.array(cent), np.array(prof), aligned, emp


def build_calibrations(verbose=True):
    """{group: dict(radii, rms, err_frac, n_lens, source)} from empirical-truth lenses.

    Derived at run time from the products on disk rather than hardcoded, so the maps
    track the kernels. Uses the REAL injected build (psf_kernel_injected.fits) where
    one exists -- currently F160W -- and the analytic drop-broadened model elsewhere,
    which is right in the median but unreliable per-lens (hence group medians only).
    """
    groups = {}
    for sample, lens, filt, entry in collect_targets():
        d = os.path.join(ws_path, 'data', 'psf', sample, lens, filt)
        emp_path = os.path.join(d, 'psf_kernel.fits')
        if not os.path.isfile(emp_path):
            continue
        emp = unit(fits.getdata(emp_path).astype(float))

        inj_path = os.path.join(d, 'psf_kernel_injected.fits')
        if os.path.isfile(inj_path):
            mod = unit(fits.getdata(inj_path).astype(float))
            src = 'injected'
        else:
            try:
                _, mod, _ = build_model_kernel(sample, lens, filt, emp.shape[0])
            except Exception:
                continue
            src = 'analytic'

        radii, rms, aligned, _ = _residual_profile(emp, mod)
        resid = float(np.sqrt(((aligned - emp) ** 2).sum()))
        boot = entry.get('psf_err_frac') or 0.0
        corrected = float(np.sqrt(max(resid ** 2 - boot ** 2, 0.0)))
        _, _, wing_frac = split_residuals(aligned, emp)

        g = calibration_group(sample, filt)
        groups.setdefault(g, dict(radii=radii, rms=[], corrected=[], wing=[],
                                  sources=set()))
        if len(rms) != len(groups[g]['radii']):
            continue                      # different kernel size; skip rather than pad
        groups[g]['rms'].append(rms)
        groups[g]['corrected'].append(corrected)
        groups[g]['wing'].append(wing_frac)
        groups[g]['sources'].add(src)

    out = {}
    for g, v in groups.items():
        if not v['rms']:
            continue
        out[g] = dict(radii=v['radii'],
                      rms=np.sqrt(np.mean(np.stack(v['rms']) ** 2, axis=0)),
                      err_frac=float(np.median(v['corrected'])),
                      wing_frac=float(np.median(v['wing'])),
                      n_lens=len(v['rms']),
                      source='+'.join(sorted(v['sources'])))
        if verbose:
            c = out[g]
            print(f'  calibration {g:22s} n={c["n_lens"]:2d} '
                  f'err_frac={c["err_frac"]:.4f} wing={100 * c["wing_frac"]:+.1f}% '
                  f'[{c["source"]}]')
    return out


def calibrated_err(calib, shape):
    """Error map on `shape` with the calibration's radial shape, scaled to its err_frac."""
    ny, nx = shape
    cy, cx = (ny - 1) / 2.0, (nx - 1) / 2.0
    yy, xx = np.mgrid[:ny, :nx]
    r = np.hypot(xx - cx, yy - cy)
    err = np.interp(r, calib['radii'], calib['rms'],
                    left=calib['rms'][0], right=calib['rms'][-1])
    total = float(np.sqrt((err ** 2).sum()))
    if total > 0:
        err = err * (calib['err_frac'] / total)
    return err


# ── Write ───────────────────────────────────────────────────────────────────────
def write_maps(psf_dir, cutouts_dir, prefix, err_full, trim_half, meta, dry_run):
    """Write the full + trimmed error maps; returns (err_frac, core, wing, trim_size)."""
    ny, nx = err_full.shape
    cy, cx = ny // 2, nx // 2
    err_trim = err_full[cy - trim_half:cy + trim_half + 1,
                        cx - trim_half:cx + trim_half + 1]
    err_frac = float(np.sqrt((err_full ** 2).sum()))

    yy, xx = np.mgrid[:ny, :nx]
    r = np.hypot(xx - (ny - 1) / 2.0, yy - (nx - 1) / 2.0)
    core = float(np.sqrt((err_full[r < 3.0] ** 2).sum()))
    wing = float(np.sqrt((err_full[r >= 3.0] ** 2).sum()))

    if not dry_run:
        for kind, arr, dest in (
                ('full_err', err_full, os.path.join(psf_dir, 'psf_kernel_err.fits')),
                ('trimmed_err', err_trim,
                 os.path.join(cutouts_dir, f'{prefix}_psf_err.fits'))):
            hdr = fits.Header()
            hdr['PSFKIND'] = (kind, 'per-pixel PSF std for the injected kernel')
            hdr['PSFERR'] = (meta['err_method'], 'error-map construction')
            hdr['PSFESRC'] = (meta['err_source'], 'ensemble / empirical-calibration source')
            hdr['PSFEBND'] = (meta['lower_bound'],
                              'True = ensemble scatter only, a LOWER BOUND')
            hdr['PSFEFRAC'] = (round(err_frac, 8), 'integrated PSF error fraction')
            hdr['PSFECORE'] = (round(core, 8), 'error fraction inside r=3px')
            hdr['PSFEWING'] = (round(wing, 8), 'error fraction outside r=3px')
            hdr['PSFINJ'] = (True, 'describes the injection-drizzled kernel')
            if meta.get('calib_group'):
                hdr['PSFECAL'] = (meta['calib_group'], 'calibration group transferred from')
                hdr['PSFECALN'] = (meta['calib_n'], 'lenses in the calibration')
            make_psf.write_fits(arr, hdr, dest)
    return err_frac, core, wing, int(err_trim.shape[0])


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--sample', action='append', help='restrict to this sample (repeatable)')
    p.add_argument('--filt', action='append', help='restrict to this filter (repeatable)')
    p.add_argument('--dry-run', action='store_true',
                   help='report what would be written, touch nothing')
    a = p.parse_args()

    print('Building calibrations from empirical-truth lenses ...')
    calibs = build_calibrations()
    if not calibs:
        print('  (none available -- STDPSF-tier products will be skipped)')

    targets = injected_products(a.sample, a.filt)
    print(f'\n{len(targets)} injected product(s)\n')
    print(f'{"sample":12s}{"lens":12s}{"filt":10s}{"source":12s}{"err_frac":>10s}'
          f'{"core":>9s}{"wing":>9s}{"trim":>6s}')

    n_written, n_skipped, counts = 0, [], {}
    for sample, lens, filt, entry in targets:
        psf_dir = os.path.join(ws_path, 'data', 'psf', sample, lens, filt)
        cutouts_dir = os.path.join(ws_path, 'data', 'cutouts', sample, lens, filt)
        kernel_path = os.path.join(psf_dir, 'psf_kernel.fits')
        if not os.path.isfile(kernel_path):
            n_skipped.append(f'{sample}/{lens}/{filt} (no psf_kernel.fits)')
            continue
        kernel = np.asarray(fits.getdata(kernel_path), dtype=float)

        err_full, meta = None, None
        ens, ens_meta = ensemble_err(psf_dir, sample, lens, filt)
        if ens is not None:
            err_full = ens
            meta = dict(err_method='ensemble_broadened',
                        err_source=f'analytic_ensemble/box={ens_meta["box"]:.3f}',
                        lower_bound=True)
        else:
            g = calibration_group(sample, filt)
            calib = calibs.get(g)
            if calib is None:
                n_skipped.append(f'{sample}/{lens}/{filt} (no calibration for {g})')
                continue
            err_full = calibrated_err(calib, kernel.shape)
            meta = dict(err_method='calibrated_vs_empirical',
                        err_source=f'empirical_calibration/{calib["source"]}',
                        lower_bound=False, calib_group=g, calib_n=calib['n_lens'])

        trim_size = int(entry.get('cutout_kernel_size') or kernel.shape[0])
        trim_half = (trim_size - 1) // 2
        os.makedirs(cutouts_dir, exist_ok=True)
        err_frac, core, wing, tsz = write_maps(psf_dir, cutouts_dir,
                                               cutout_prefix(cutouts_dir), err_full,
                                               trim_half, meta, a.dry_run)

        if not a.dry_run:
            rec = dict(entry)
            rec.update({'psf_err_frac': round(err_frac, 8),
                        'psf_err_frac_core': round(core, 8),
                        'psf_err_frac_wing': round(wing, 8),
                        'err_method': meta['err_method'],
                        'err_source': meta['err_source'],
                        'err_lower_bound': meta['lower_bound']})
            info_json.update(PSF_JSON, sample, lens, filt, rec)

        tag = 'ensemble' if meta['lower_bound'] else 'calibrated'
        counts[tag] = counts.get(tag, 0) + 1
        n_written += 1
        print(f'{sample:12s}{lens:12s}{filt:10s}{tag:12s}{err_frac:10.4f}'
              f'{core:9.4f}{wing:9.4f}{tsz:6d}')

    print(f'\n{n_written} written' + (' (dry run -- nothing touched)' if a.dry_run else '')
          + ''.join(f', {v} {k}' for k, v in sorted(counts.items())))
    for s in n_skipped:
        print(f'  skipped: {s}')


if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    main()
