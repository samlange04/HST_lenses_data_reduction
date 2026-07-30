#!/usr/bin/env python
"""
Tile every lens's PSF kernel into per-filter-group QC mosaics.

Sibling to make_mosaics.py: same 5-wide grid, same filter groups (mosaic_groups.py,
shared between the two scripts). Reads the trimmed, modelling-ready kernels already
written by make_psf.py
(data/cutouts/<sample>/<lens>/<filt>/cutout[_cr]_psf.fits) - nothing is rebuilt.

Kernels are unit-sum normalised by construction (make_psf.trim_kernel_to_amplitude), so
raw peak amplitude reflects kernel *size* (a broader/larger-footprint PSF has a lower
peak for the same total flux) as much as PSF sharpness. Each panel is instead
peak-normalised before display, so panels are comparable regardless of trim size. The
display uses a log stretch over 1e-4..1 (pooled_log_norm), matching the 'log' wing
panel in each lens's psf.png rather than make_mosaics.py's pooled asinh - the PSF wings
are what matters here, and they read better on a log stretch. Each panel label is
tagged 'emp' (empirical ePSF,
cut from the drizzled mosaic - the closest thing to ground truth) or 'mod' (STDPSF /
focus-diverse / MAST PSF DB - a detector-frame model resampled to North-up), from the
PSFMETH keyword make_psf.py stamps on every kernel.

Only slacs_gold has PSF products today (info/lens_psf.json has just that sample key).
slacs_other uses the same three instruments as slacs_gold, so `bash
scripts/run_psf_all.sh slacs_other` should build its kernels with no code changes.
gallery is WFC3/UVIS, which make_psf.py/psf_models.py has no instrument support for yet
(no detector scale, no STDPSF grid config) - this script needs no changes when that
support lands, it will just start finding cutout_psf.fits files under
data/cutouts/gallery/*/<filt>/ once they exist (see mosaic_groups.SAMPLE_GROUPS for
gallery's 5 filters, already listed there).

Usage:
    uv run python scripts/make_psf_mosaics.py --sample slacs_gold
"""

import argparse
import glob
import os

import numpy as np
from astropy.io import fits
from astropy.visualization import ImageNormalize, LogStretch

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mast_target_names
import mosaic_groups
from make_mosaics import short_filt, plot_mosaic

ws_path = '/Users/samlange/Code/HST_lenses_data_reduction'


def pooled_log_norm(arrays):
    """Log stretch over 1e-4..1, matching the 'log' wing panel in each lens's psf.png
    (make_psf.plot_psf: LogStretch, vmin=peak*1e-4, vmax=peak). Every panel here is
    already peak-normalised to 1 (build_group divides by its own peak), so a fixed
    1e-4..1 log norm reproduces that per-kernel view across the whole mosaic. `arrays`
    is unused - the range is fixed by the peak-normalisation, not pooled from pixels."""
    return ImageNormalize(vmin=1e-4, vmax=1.0, stretch=LogStretch())


def find_psf_path(filt_dir):
    """Prefer the CR-pass kernel, matching make_mosaics.find_cutout_pair's precedence."""
    for prefix in ('cutout_cr', 'cutout'):
        path = os.path.join(filt_dir, f'{prefix}_psf.fits')
        if os.path.exists(path):
            return path
    return None


def method_tag(method):
    """'emp' for the empirical ePSF (cut from the North-up mosaic, the drizzled truth);
    'mod' for any model tier (STDPSF, ACS focus-diverse, WFPC2 MAST-DB, injected)."""
    return 'emp' if str(method).startswith('empirical') else 'mod'


def build_group(cutouts_dir, precedence):
    """Entries for one mosaic group - see make_mosaics.build_group for the precedence
    convention. `group` (the per-panel colourbar-split key) is only set for multi-filter
    groups, same reasoning as make_mosaics.py."""
    entries = []
    for lens_dir in sorted(glob.glob(os.path.join(cutouts_dir, '*'))):
        lens = os.path.basename(lens_dir)
        for filt in precedence:
            path = find_psf_path(os.path.join(lens_dir, filt))
            if path is None:
                continue
            with fits.open(path) as hdul:
                data = hdul[0].data.astype(np.float64)
                method = hdul[0].header.get('PSFMETH', '?')
            peak = float(np.nanmax(data))
            if peak > 0:
                data = data / peak
            tag = method_tag(method)
            if len(precedence) > 1:
                short = short_filt(filt)
                label, group = f'{lens} [{short}] {tag}', short.split('_')[0]
            else:
                label, group = f'{lens} {tag}', None
            entries.append((lens, label, data, group))
            break
    return entries


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--sample', default=mast_target_names.DEFAULT_SAMPLE,
                   help='sample subdirectory under data/cutouts/ to mosaic. Defined in '
                        f'info/lens_samples.json (default {mast_target_names.DEFAULT_SAMPLE})')
    a = p.parse_args()

    cutouts_dir = os.path.join(ws_path, 'data', 'cutouts', a.sample)
    out_dir = os.path.join(ws_path, 'data', 'mosaics', a.sample)
    os.makedirs(out_dir, exist_ok=True)

    groups = mosaic_groups.groups_for_sample(a.sample, cutouts_dir)
    for group_name, precedence in groups.items():
        entries_raw = build_group(cutouts_dir, precedence)
        if not entries_raw:
            print(f"{group_name}: no PSF kernels found under {cutouts_dir}, skipping")
            continue

        entries = [{'label': label} for _, label, _, _ in entries_raw]
        arrays = [data for _, _, data, _ in entries_raw]
        raw_groups = [group for _, _, _, group in entries_raw]
        split_by = raw_groups if any(g is not None for g in raw_groups) else None
        print(f"{group_name}: {len(entries)} PSF kernels")

        plot_mosaic(entries, arrays, 'label',
                    os.path.join(out_dir, f'{group_name}_psf.png'),
                    'PSF (peak-normalised, log)', split_by=split_by,
                    norm_fn=pooled_log_norm)


if __name__ == '__main__':
    main()
