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

    conda run -n stenv python scripts/align_wfpc2_to_acs.py --lens J0252+0039
    conda run -n stenv python scripts/align_wfpc2_to_acs.py --all
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

WS = '/Users/samlange/Code/data_reduction'
sys.path.insert(0, os.path.join(WS, 'info'))
from slacs_coords import slacs_coords

MAX_SHIFT = 1.5   # arcsec; a larger apparent offset means the tie failed -> skip


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


def find_product(lens, filt, prefix):
    """The CR-pass sci product for a band (falls back to no-CR)."""
    d = os.path.join(WS, 'data', 'drizzled', 'slacs', lens, filt)
    for tag in ('cr', 'nocrrej'):
        hit = glob.glob(os.path.join(d, f'{prefix}_{tag}_*_sci.fits'))
        if hit:
            return hit[0]
    return None


def align_lens(lens):
    cat = SkyCoord(*slacs_coords[lens], unit=(u.hourangle, u.deg))
    ref = find_product(lens, 'f814W', 'acs_wfc_flc')      # GAIA-accurate reference
    f606 = find_product(lens, 'f606W', 'wfpc2_wf3')
    if ref is None:
        print(f'{lens}: no ACS F814W product to align against — skip')
        return
    if f606 is None:
        print(f'{lens}: no WFPC2 F606W product — skip')
        return
    cref = stable_centroid(ref, cat)
    c606 = stable_centroid(f606, cat)
    if cref is None or c606 is None:
        print(f'{lens}: deflector centroid failed — skip')
        return
    sep = cref.separation(c606).arcsec
    if sep > MAX_SHIFT:
        print(f'{lens}: measured offset {sep:.2f}" > {MAX_SHIFT}" — refusing (tie looks wrong)')
        return
    dra = cref.ra.deg - c606.ra.deg
    ddec = cref.dec.deg - c606.dec.deg
    # apply to every F606W product (both passes, sci + wht)
    d = os.path.join(WS, 'data', 'drizzled', 'slacs', lens, 'f606W')
    n = 0
    for fn in glob.glob(os.path.join(d, 'wfpc2_wf3_*_drw_*.fits')):
        with fits.open(fn, mode='update') as h:
            h[0].header['CRVAL1'] += dra
            h[0].header['CRVAL2'] += ddec
            h[0].header['GSC240FX'] = (True, 'CRVAL shifted to GAIA/ACS deflector frame')
            h.flush()
        n += 1
    # verify independently against F160W if present
    f160 = find_product(lens, 'f160W', 'wfc3_ir_flt')
    chk = ''
    if f160:
        c160 = stable_centroid(f160, cat)
        cnew = stable_centroid(f606, cat)
        if c160 is not None and cnew is not None:
            chk = f'  (F606W->F160W check: {cnew.separation(c160).arcsec:.4f}")'
    print(f'{lens}: F606W offset {sep:.3f}" -> corrected {n} files '
          f'(dRA={dra*3600:+.3f}", dDec={ddec*3600:+.3f}"){chk}')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--lens')
    p.add_argument('--all', action='store_true')
    a = p.parse_args()
    if a.all:
        lenses = sorted(
            os.path.basename(os.path.dirname(d))
            for d in glob.glob(os.path.join(WS, 'data', 'drizzled', 'slacs', '*', 'f606W'))
        )
        for lens in lenses:
            align_lens(lens)
    elif a.lens:
        align_lens(a.lens)
    else:
        p.error('give --lens LENS or --all')
