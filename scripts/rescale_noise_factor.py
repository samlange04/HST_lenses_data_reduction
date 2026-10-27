#!/usr/bin/env python
"""Bring stamp noise maps in line with info/overrides/noise_corr_factors.json without re-cutting.

make_stamps.py writes noise = NOISEK * NOISECOR / sqrt(WHT) and records NOISECOR/NOISECRS in
the noise header. When the JSON's factor for a stamp changes, multiplying the noise map by
new/old is exact, and the science pixels, pixel grid, masks and positions are untouched. A
re-cut, by contrast, would re-read the mosaic. A stamp with the PSF error folded in (PSFERRAP)
is NOT linear in the factor and is skipped with a message; re-cut it with make_stamps.py
--center-from instead.

    uv run python scripts/rescale_noise_factor.py --sample slacs_gold bells    # dry run
    uv run python scripts/rescale_noise_factor.py --sample slacs_gold bells --apply
"""
import argparse
import glob
import os
import sys

import numpy as np
from astropy.io import fits

WS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(WS, 'info', 'coords'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stamp_paths                                   # noqa: E402
import mast_target_names                              # noqa: E402
from make_stamps import resolve_corr_factor          # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--sample', nargs='+', default=mast_target_names.samples())
    p.add_argument('--apply', action='store_true', help='write (default: dry run)')
    a = p.parse_args()
    root = stamp_paths.stamps_root(WS)
    n_change = n_skip = 0
    for sample in a.sample:
        for nf in sorted(glob.glob(os.path.join(root, sample, '*', '*', 'stamp*_noise.fits'))):
            lens, filt = nf.split(os.sep)[-3:-1]
            sci_hdr = fits.getheader(nf.replace('_noise.fits', '_sci.fits'))
            with fits.open(nf, mode='update' if a.apply else 'readonly') as h:
                hdr = h[0].header
                old = float(hdr.get('NOISECOR', 1.0))
                new, src = resolve_corr_factor(None, sci_hdr, sample, lens, filt, WS)
                if abs(new - old) < 1e-9:
                    continue
                tag = f'{sample} {lens} {filt} {os.path.basename(nf)}: {old:g} -> {new:g} ({src})'
                if hdr.get('PSFERRAP'):
                    print(f'SKIP (PSF error folded; re-cut instead) {tag}')
                    n_skip += 1
                    continue
                print(('rescaled ' if a.apply else 'would rescale ') + tag)
                n_change += 1
                if a.apply:
                    h[0].data = (np.asarray(h[0].data, dtype=np.float64) * (new / old)
                                 ).astype(h[0].data.dtype)
                    hdr['NOISECOR'] = (new, 'noise-map inflation applied (see NOISECRS)')
                    hdr['NOISECRS'] = (src, 'NOISECOR source: json/lens/cli/none')
                    hdr['NOISERSC'] = (f'{old:g}->{new:g}', 'rescaled in place by rescale_noise_factor.py')
                    h.flush()
    print(f'{n_change} {"rescaled" if a.apply else "to rescale"}, {n_skip} skipped (PSF error folded)')


if __name__ == '__main__':
    main()
