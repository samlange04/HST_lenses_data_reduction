#!/usr/bin/env python
"""
Build a PSF for a star-poor product from the median of other fields' empirical PSFs.

For a band with no STDPSF model grid (WFC3/UVIS F475X, F600LP), a field with too few usable
stars for make_psf.py's empirical ePSF has no model fallback. Chosen 2026-10-04 for
those products: take the per-pixel MEDIAN of the empirical kernels already built for the
same filter in other lenses (same instrument pixel scale), as the drizzle products are all
on one output grid per band. Used for SL2S J0233-0438 F475X (0 usable stars; 29 F475X
donors) and J1358+5459 F600LP (2 stars; 11 F600LP donors).

Outputs mirror make_psf.py exactly, so everything downstream reads them unchanged:
    data/psf/<sample>/<lens>/<filt>/psf_kernel.fits, psf_kernel_err.fits, psf.png
    <psf cutout dir>/cutout_[cr_]psf.fits, cutout_[cr_]psf_err.fits
and info/lens_psf.json gets method 'empirical_median' with the donor list.

The error map is the donor-to-donor per-pixel std (NOT the std of the median): the target's
true PSF is one more draw from that field-to-field spread (focus/breathing, position on the
chip, drizzle geometry), so the scatter of single fields is its honest uncertainty.

    uv run python scripts/make_psf_median.py --lens J0233-0438 --filt f475X --sample sl2s_dinos_x

**Single-donor mode** (``--donor SAMPLE/LENS``, method ``empirical_donor``): the kernel is
that one lens's empirical PSF instead of the median, and the error map is still the
donor-pool scatter. Used for BELLS J0830+5116 F814W (2026-10-05). Its focus-diverse ePSF
does not exist on STScI's service, and only 2 field stars are usable. Fitted to those 2 stars,
J0747+4448's empirical PSF (observed the same day) beat all 22 other BELLS empirical PSFs, the
BELLS median and the STDPSF+injection kernel on the bright star: chi2/dof 156 against 192-345,
224 and 251. It ranked second on the faint star.

    uv run python scripts/make_psf_median.py --lens J0830+5116 --filt f814W --sample bells \
        --donor bells/J0747+4448 --pool-samples bells
"""
import argparse
import glob
import os
import sys

import numpy as np
from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ws_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ws_path, 'info'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import info_json
import cutout_paths
# make_psf's work is under main(), so importing it only brings in the helpers.
from make_psf import trim_kernel_to_amplitude, subtract_pedestal, measure_fwhm, write_fits


def donor_kernels(filt, target_lens, scale, psf_json, pool_samples=None):
    """[(sample, lens, unit-sum full kernel)] for every OTHER empirical build of `filt`
    at pixel scale `scale`, across all samples (or only `pool_samples`)."""
    rec = info_json.load(psf_json)
    out = []
    for path in sorted(glob.glob(os.path.join(ws_path, 'data', 'psf', '*', '*', filt,
                                              'psf_kernel.fits'))):
        sample, lens = path.split(os.sep)[-4:-2]
        if lens == target_lens or (pool_samples and sample not in pool_samples):
            continue
        r = (rec.get(sample, {}).get(lens) or {}).get(filt) or {}
        if r.get('method') != 'empirical':
            continue
        with fits.open(path) as h:
            k, hdr = h[0].data.astype(float), h[0].header
        if abs(float(hdr.get('PSFPXSCL', -1)) - scale) > 1e-4:
            continue
        out.append((sample, lens, k / k.sum()))
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--lens', required=True)
    p.add_argument('--filt', required=True)
    p.add_argument('--sample', required=True)
    p.add_argument('--trim-threshold', dest='trim_threshold', type=float, default=1e-3)
    p.add_argument('--min-donors', type=int, default=5)
    p.add_argument('--donor', default=None, metavar='SAMPLE/LENS',
                   help='use this one empirical kernel instead of the median (error map is '
                        'still the donor-pool scatter)')
    p.add_argument('--pool-samples', nargs='+', default=None,
                   help='restrict the donor pool to these samples (default: all)')
    a = p.parse_args()

    psf_json = os.path.join(ws_path, 'info', 'lens_psf.json')
    # Pixel scale and drizzle pass come from the target's own cutout, so the kernel matches it.
    stamp_dir = cutout_paths.psf_cutout_dir(ws_path, a.sample, a.lens, a.filt)
    sci = (glob.glob(os.path.join(stamp_dir, 'cutout_cr_sci.fits'))
           or glob.glob(os.path.join(stamp_dir, 'cutout_sci.fits')))
    if not sci:
        sys.exit(f'no cutout for {a.lens} {a.filt} in {stamp_dir}; run make_cutouts.py first')
    prefix = 'cutout_cr' if os.path.basename(sci[0]).startswith('cutout_cr') else 'cutout'
    drizzle_pass = 'cr' if prefix == 'cutout_cr' else 'nocrrej'
    from astropy.wcs import WCS
    from astropy.wcs.utils import proj_plane_pixel_scales
    scale = float(proj_plane_pixel_scales(WCS(fits.getheader(sci[0])).celestial)[0] * 3600.0)
    scale = round(scale, 4)

    donors = donor_kernels(a.filt, a.lens, scale, psf_json, a.pool_samples)
    print(f'{a.lens} {a.filt}: {len(donors)} empirical donor kernels at {scale}"/px')
    if len(donors) < a.min_donors:
        sys.exit(f'only {len(donors)} donors (< {a.min_donors}); not building a median PSF')
    shapes = {k.shape for _, _, k in donors}
    if len(shapes) != 1:
        sys.exit(f'donor kernels differ in shape {shapes}; refusing to combine')
    stack = np.stack([k for _, _, k in donors])

    if a.donor:
        pick = [k for s_, l_, k in donors if f'{s_}/{l_}' == a.donor]
        if not pick:
            sys.exit(f'--donor {a.donor} is not an empirical {a.filt} kernel in the pool')
        kernel = pick[0].copy()
        method = 'empirical_donor'
    else:
        kernel = np.median(stack, axis=0)
        method = 'empirical_median'
    kernel /= kernel.sum()
    kernel, pedestal_frac = subtract_pedestal(kernel)
    trimmed, trim_size = trim_kernel_to_amplitude(kernel, a.trim_threshold)
    trim_half = trim_size // 2
    fwhm_pix = measure_fwhm(kernel)

    err_full = stack.std(axis=0, ddof=1)
    cy, cx = err_full.shape[0] // 2, err_full.shape[1] // 2
    err_trim = err_full[cy - trim_half:cy + trim_half + 1, cx - trim_half:cx + trim_half + 1]
    err_frac = float(np.sqrt(np.sum(err_full ** 2)))
    fwhms = [f for f in (measure_fwhm(k) for _, _, k in donors) if f is not None]
    fwhm_err = float(np.std(fwhms, ddof=1)) if len(fwhms) > 1 else None
    print(f'  {method} kernel FWHM {fwhm_pix:.3f} px (donor spread {fwhm_err:.3f}), '
          f'trimmed {kernel.shape[0]} -> {trim_size} px, PSF error fraction {err_frac:.2e}')

    output_dir = os.path.join(ws_path, 'data', 'psf', a.sample, a.lens, a.filt)
    os.makedirs(output_dir, exist_ok=True)
    khdr = fits.Header()
    khdr['PSFMETH'] = (method, 'other fields\' empirical PSF(s): median or one donor')
    if a.donor:
        khdr['PSFDONOR'] = (a.donor, 'donor whose empirical PSF this is')
    khdr['PSFNSTAR'] = (0, 'stars used in THIS field (none)')
    khdr['PSFNDON'] = (len(donors), 'donor empirical kernels median-combined')
    khdr['PSFPXSCL'] = (scale, 'kernel pixel scale (arcsec)')
    khdr['PSFFWHM'] = (round(fwhm_pix, 4) if fwhm_pix is not None else 0.0,
                       'fitted kernel FWHM (pixels)')
    khdr['PSFLENS'] = (a.lens, 'lens')
    khdr['PSFFILT'] = (a.filt, 'filter')
    khdr['PSFKIND'] = ('full', 'full footprint (trimmed copy in the psf cutout dir)')
    khdr['PSFPED'] = (round(pedestal_frac, 6), 'wing pedestal removed (fraction of peak)')
    write_fits(kernel, khdr, os.path.join(output_dir, 'psf_kernel.fits'))

    thdr = khdr.copy()
    thdr['PSFKIND'] = ('trimmed', 'amplitude-trimmed modelling kernel')
    thdr['PSFTRIM'] = (a.trim_threshold, 'azimuthal PSF < this x peak sets the radius')
    thdr['PSFPASS'] = (drizzle_pass, 'drizzle pass the PSF matches')
    write_fits(trimmed, thdr, os.path.join(stamp_dir, f'{prefix}_psf.fits'))

    for base, kind, arr, dest in (
            (khdr, 'full_err', err_full, os.path.join(output_dir, 'psf_kernel_err.fits')),
            (thdr, 'trimmed_err', err_trim, os.path.join(stamp_dir, f'{prefix}_psf_err.fits'))):
        ehdr = base.copy()
        ehdr['PSFKIND'] = (kind, 'per-pixel donor-to-donor std')
        ehdr['PSFERR'] = ('donor_scatter', 'std of the single-field donor kernels')
        ehdr['PSFESRC'] = ('donors', 'ensemble source')
        ehdr['PSFNVALD'] = (len(donors), 'donor kernels')
        ehdr['PSFEFRAC'] = (round(err_frac, 8), 'integrated PSF error fraction')
        write_fits(arr, ehdr, dest)

    fig, ax = plt.subplots(1, 3, figsize=(12, 4))
    for axx, img, t in ((ax[0], kernel, 'kernel'), (ax[1], err_full, 'donor std'),
                        (ax[2], err_full / np.maximum(kernel, 1e-12), 'std / kernel')):
        im = axx.imshow(np.log10(np.clip(img, 1e-6, None)) if t != 'std / kernel'
                        else np.clip(img, 0, 2), origin='lower')
        axx.set_title(t)
        plt.colorbar(im, ax=axx, fraction=0.046)
    fig.suptitle(f'{a.lens}  {a.filt}  [{method}{" " + a.donor if a.donor else ""}]  '
                 f'({len(donors)} donors in the scatter pool)')
    fig.savefig(os.path.join(output_dir, 'psf.png'), dpi=80, bbox_inches='tight')
    plt.close(fig)

    info_json.update(psf_json, a.sample, a.lens, a.filt, {
        'method': method,
        **({'donor': a.donor} if a.donor else {}),
        'n_stars': 0,
        'n_donors': len(donors),
        'donors': [f'{s}/{l}' for s, l, _ in donors],
        'fwhm_pix': round(fwhm_pix, 4) if fwhm_pix is not None else None,
        'oversample': None,
        'kernel_size': int(kernel.shape[0]),
        'cutout_kernel_size': int(trim_size),
        'trim_threshold': a.trim_threshold,
        'pedestal_frac': round(pedestal_frac, 6),
        'fwhm_pix_err': round(fwhm_err, 4) if fwhm_err is not None else None,
        'psf_err_frac': round(err_frac, 8),
        'err_method': 'donor_scatter',
        'err_source': 'donors',
        'n_boot_valid': len(donors),
    })
    print('  done')


if __name__ == '__main__':
    main()
