"""Model-tier PSF error measured against empirical truth using the REAL injected kernel.

READ-ONLY QC. The companion `psf_model_error.py` approximates the injection promotion
with `psf_models.analytic_drop_broaden`, because a real injection needs a re-drizzle
per lens. That stand-in is right in the MEDIAN but unreliable per-lens (checked
2026-08-02 on the five F160W lenses that have both: individual residuals moved by up
to 2.7x, in both directions, while the group median only moved 0.032 -> 0.037). This
script does the honest version -- comparing `psf_kernel_injected.fits` against the
empirical `psf_kernel.fits` on the same lens+filter -- for products where a real
injected build exists.

**Producing the injected builds first.** On an EMPIRICAL primary, injection is
non-destructive: `make_psf_inject.run_injection(lens, filt, sample)` auto-detects
`promote=False` and writes only the parallel `*_injected` names, leaving the canonical
`psf_kernel.fits` / `cutout_[cr_]psf.fits` untouched (verified by mtime). That is what
`run_psf_inject_all.sh --all` does. Requires `data/drizzle_files/<sample>/<lens>/<filt>/`
to still hold the science drizzle inputs -- re-run the band's drizzle first if cleared.

Measured F160W result (5 lenses, 2026-08-02): model error median **0.036** against a
bootstrap ensemble scatter of 0.0086 -- a ~4x underestimate -- with the model's
enclosed wing flux a median 9.9% low (-6.4 to -13.4%) and the residual 94% core-
dominated. The injected model is simultaneously slightly too broad in the core and too
faint in the wings, so one scalar cannot represent this error. J0728+3835 reproduced
the earlier prototype's 3.93px injected FWHM exactly. See memory:
stdpsf_model_error_measured.

    uv run python scripts/psf_model_error_injected.py
    uv run python scripts/psf_model_error_injected.py --sample slacs_gold --filt f160W
"""
import argparse
import json
import os
import warnings

import numpy as np
from astropy.io import fits

from psf_model_error import (align_resid, collect_targets, fwhm_radial,
                             split_residuals, unit, ws_path)


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
    rows, missing = [], []

    print(f'{"sample":12s}{"lens":12s}{"filt":7s}{"nst":>4s}{"resid":>8s}{"boot":>8s}'
          f'{"corr":>8s}{"wing%":>8s}{"fwhm_e":>8s}{"fwhm_m":>8s}')
    for sample, lens, filt, entry in targets:
        d = os.path.join(ws_path, 'data', 'psf', sample, lens, filt)
        emp_path = os.path.join(d, 'psf_kernel.fits')
        inj_path = os.path.join(d, 'psf_kernel_injected.fits')
        if not (os.path.isfile(emp_path) and os.path.isfile(inj_path)):
            missing.append(f'{sample}/{lens}/{filt}')
            continue

        emp = unit(fits.getdata(emp_path).astype(float))
        mod = unit(fits.getdata(inj_path).astype(float))
        resid, dy, dx, aligned = align_resid(mod, emp)
        rc, rw, wf = split_residuals(aligned, emp, a.core_radius)
        boot = entry.get('psf_err_frac') or 0.0
        corr = float(np.sqrt(max(resid ** 2 - boot ** 2, 0.0)))

        rows.append(dict(sample=sample, lens=lens, filt=filt,
                         n_stars=entry.get('n_stars'), boot=boot, resid=resid,
                         corrected=corr, resid_core=rc, resid_wing=rw, wing_frac=wf,
                         dy=dy, dx=dx, fwhm_emp=fwhm_radial(emp),
                         fwhm_model=fwhm_radial(mod)))
        print(f'{sample:12s}{lens:12s}{filt:7s}{str(entry.get("n_stars")):>4s}'
              f'{resid:8.4f}{boot:8.4f}{corr:8.4f}{100 * wf:+8.1f}'
              f'{fwhm_radial(emp):8.2f}{fwhm_radial(mod):8.2f}')

    if not rows:
        print('\nno products have BOTH psf_kernel.fits and psf_kernel_injected.fits.')
        print('Build the injected comparison products first (non-destructive on an '
              'empirical primary):  bash scripts/run_psf_inject_all.sh --all')
        return

    d = np.array([r['resid'] for r in rows])
    b = np.array([r['boot'] for r in rows])
    c = np.array([r['corrected'] for r in rows])
    w = 100 * np.array([r['wing_frac'] for r in rows])
    print(f'\nn={len(rows)}  resid median={np.median(d):.4f} '
          f'[{d.min():.4f}-{d.max():.4f}]  corrected median={np.median(c):.4f}')
    print(f'bootstrap median={np.median(b):.4f}  -> model error is '
          f'{np.median(d) / np.median(b) if np.median(b) else float("nan"):.1f}x '
          f'the ensemble scatter')
    print(f'wing flux deficit: median {np.median(w):+.1f}% '
          f'[{w.min():+.1f}% .. {w.max():+.1f}%]')
    print(f'core share of residual: '
          f'{np.median([r["resid_core"] / r["resid"] for r in rows if r["resid"]]):.3f}')
    if missing:
        print(f'\nno injected build for {len(missing)} product(s): '
              f'{", ".join(missing[:6])}{" ..." if len(missing) > 6 else ""}')

    if a.out:
        with open(a.out, 'w') as fh:
            json.dump(rows, fh, indent=1)
        print(f'\nwrote {a.out}  ({len(rows)} products)')


if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    main()
