#!/usr/bin/env python
"""
Correct the WFPC2 F606W absolute astrometry to the GAIA frame.

WFPC2 SLACS frames carry only a GSC 2.4.0 (`GSC240`) astrometric solution, which has
~0.3-1" absolute error. ACS (F814W) and WFC3/IR (F160W) carry GAIA eDR3 / GSC242
solutions accurate to <0.02" and agree with each other to ~0.01". So the WFPC2 F606W
mosaic sits ~0.5-0.9" off from the other bands (measured 0.66" on J0252+0039).

This script registers each F606W mosaic to its ACS F814W counterpart using the
deflector as the tie point: it measures the deflector light-centroid (iterative
windowed centroid, robust to the roughly symmetric lensed ring) in both bands and
shifts the F606W CRVAL so the two coincide. Verified on J0252+0039: F606W went from
0.66" to 0.009" of F814W, and independently to 0.009" of F160W (F160W is not used in
the fit, so it is a clean check that F606W landed on the GAIA frame).

Run AFTER both the ACS and WFPC2 drizzles exist, BEFORE make_cutouts. Idempotent:
re-running re-measures the residual and applies ~0.

    uv run python scripts/align_wfpc2_to_acs.py --lens J0252+0039
    uv run python scripts/align_wfpc2_to_acs.py --all

--target / --ref generalise the same tie to any band pair, because the band that needs
correcting is not *always* F606W. The assumption in the paragraphs above -- that ACS and
WFC3/IR always carry a GAIA-grade solution -- is a property of the delivered WCS, not of
the instrument, and it does fail: J1016+3859's ACS F814W (the force_copy COPY visit) came
down with a bare `-GSC240` fit and sits 0.60" off its own F160W, which is what makes it
the reference-band trap for make_cutouts' shared --center-band. Check WCSNAME before
assuming which band is the truth; the defaults (target F606W, ref F814W) are unchanged
and remain right for every other SLACS lens.

    uv run python scripts/align_wfpc2_to_acs.py --lens J1016+3859 --sample slacs_other \
        --target f814W --ref f160W
"""
import argparse
import glob
import os
import sys
import warnings
warnings.filterwarnings('ignore')
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.nddata import Cutout2D
from scipy.ndimage import median_filter
from photutils.centroids import centroid_com

WS = '/Users/samlange/Code/HST_lenses_data_reduction'
sys.path.insert(0, os.path.join(WS, 'info'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from slacs_coords import slacs_coords
import mast_target_names

MAX_SHIFT = 1.5   # arcsec; a larger apparent offset means the tie failed -> skip

# Product-filename prefix per band directory, so --target/--ref can name any band. The
# prefix is set by the drizzle script that wrote the band, not by the filter, which is why
# f814W and f555W share one (both come from drizzle_acs_wfc.py) and the split-visit
# f606W_v2 shares WFPC2's.
BAND_PREFIX = {
    'f606W':    'wfpc2_wf3',
    'f606W_v2': 'wfpc2_wf3',
    'f814W':    'acs_wfc_flc',
    'f555W':    'acs_wfc_flc',
    'f160W':    'wfc3_ir_flt',
}


def band_prefix(filt):
    """Product prefix for a band directory, tolerating any f606W_v<N> split-visit key."""
    if filt in BAND_PREFIX:
        return BAND_PREFIX[filt]
    base = filt.split('_')[0]
    if base in BAND_PREFIX:
        return BAND_PREFIX[base]
    raise KeyError(f'no product prefix known for band {filt!r}; add it to BAND_PREFIX')


def stable_centroid(path, cat, box=2.5, niter=8, tol=0.003):
    """Iterative windowed light-centroid of the deflector, as a SkyCoord.

    The window re-centres on the running centroid each pass, so it converges to the
    galaxy light centroid independent of the starting guess and of CRVAL. The lensed
    ring is roughly symmetric, so it does not bias the centroid.
    """
    h = fits.open(path)[0]
    w = WCS(h.header)
    pos = cat
    for _ in range(niter):
        c = Cutout2D(np.nan_to_num(h.data), pos, u.Quantity((box, box), u.arcsec), wcs=w)
        d = np.clip(median_filter(c.data, 3) - np.percentile(c.data, 30), 0, None)
        if d.sum() <= 0:
            return None
        cx, cy = centroid_com(d)
        newpos = c.wcs.pixel_to_world(cx, cy)
        if pos.separation(newpos).arcsec < tol:
            return newpos
        pos = newpos
    return pos


def find_product(lens, filt, prefix, sample=mast_target_names.DEFAULT_SAMPLE):
    """The CR-pass sci product for a band (falls back to no-CR)."""
    d = os.path.join(WS, 'data', 'drizzled', sample, lens, filt)
    for tag in ('cr', 'nocrrej'):
        hit = glob.glob(os.path.join(d, f'{prefix}_{tag}_*_sci.fits'))
        if hit:
            return hit[0]
    return None


def align_lens(lens, f606_dir='f606W', sample=mast_target_names.DEFAULT_SAMPLE,
               target=None, ref_filt='f814W', check_filt='f160W'):
    """Shift `target`'s CRVAL so its deflector centroid coincides with `ref_filt`'s.

    `target` defaults to `f606_dir` (the F606W band directory), keeping the original
    F606W->F814W behaviour. `check_filt` is measured but never fitted, so it stays an
    independent check that the target landed on the right frame; it is skipped when it
    is one of the two bands in the tie.
    """
    target = target or f606_dir
    cat = SkyCoord(*slacs_coords[lens], unit=(u.hourangle, u.deg))
    ref = find_product(lens, ref_filt, band_prefix(ref_filt), sample)
    tgt = find_product(lens, target, band_prefix(target), sample)
    if ref is None:
        print(f'{lens}: no {ref_filt} product to align against — skip')
        return
    if tgt is None:
        print(f'{lens}: no {target} product — skip')
        return
    cref = stable_centroid(ref, cat)
    ctgt = stable_centroid(tgt, cat)
    if cref is None or ctgt is None:
        print(f'{lens}: deflector centroid failed — skip')
        return
    sep = cref.separation(ctgt).arcsec
    if sep > MAX_SHIFT:
        print(f'{lens}: measured offset {sep:.2f}" > {MAX_SHIFT}" — refusing (tie looks wrong)')
        return
    dra = cref.ra.deg - ctgt.ra.deg
    ddec = cref.dec.deg - ctgt.dec.deg
    # apply to every product of the target band (both passes, sci + wht)
    d = os.path.join(WS, 'data', 'drizzled', sample, lens, target)
    n = 0
    for fn in glob.glob(os.path.join(d, f'{band_prefix(target)}_*.fits')):
        with fits.open(fn, mode='update') as h:
            h[0].header['CRVAL1'] += dra
            h[0].header['CRVAL2'] += ddec
            h[0].header['GSC240FX'] = (True, f'CRVAL shifted to {ref_filt} deflector frame')
            h[0].header['ASTROREF'] = (ref_filt, 'band this product was astrometrically tied to')
            h.flush()
        n += 1
    # verify against a band that took no part in the fit
    chk = ''
    if check_filt not in (target, ref_filt):
        fchk = find_product(lens, check_filt, band_prefix(check_filt), sample)
        if fchk:
            cchk = stable_centroid(fchk, cat)
            cnew = stable_centroid(tgt, cat)
            if cchk is not None and cnew is not None:
                chk = (f'  ({target}->{check_filt} check: '
                       f'{cnew.separation(cchk).arcsec:.4f}")')
    print(f'{lens}: {target} offset {sep:.3f}" from {ref_filt} -> corrected {n} files '
          f'(dRA={dra*3600:+.3f}", dDec={ddec*3600:+.3f}"){chk}')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--lens')
    p.add_argument('--all', action='store_true')
    p.add_argument('--f606-dir', default='f606W',
                   help='F606W band subdirectory (e.g. f606W_v2 for a split lens\'s second visit)')
    p.add_argument('--target', default=None,
                   help='band whose CRVAL is shifted (default: --f606-dir, i.e. F606W). Use '
                        'this when the mis-fit band is not F606W -- check WCSNAME first')
    p.add_argument('--ref', dest='ref_filt', default='f814W',
                   help='band providing the reference astrometric frame (default f814W). '
                        'Set to f160W where ACS itself carries only a GSC240 fit')
    p.add_argument('--sample', default=mast_target_names.DEFAULT_SAMPLE,
                   help='sample subdirectory of data/drizzled/. Defined in '
                        f'info/lens_samples.json (default {mast_target_names.DEFAULT_SAMPLE})')
    a = p.parse_args()
    if a.all:
        # Glob the products rather than the sample list: only a lens that actually has a
        # drizzled F606W product can be tied, and split-visit lenses have no bare f606W.
        lenses = sorted(
            os.path.basename(os.path.dirname(d))
            for d in glob.glob(os.path.join(WS, 'data', 'drizzled', a.sample, '*', 'f606W'))
        )
        for lens in lenses:
            align_lens(lens, sample=a.sample, ref_filt=a.ref_filt)
    elif a.lens:
        align_lens(a.lens, f606_dir=a.f606_dir, sample=a.sample,
                   target=a.target, ref_filt=a.ref_filt)
    else:
        p.error('give --lens LENS or --all')
