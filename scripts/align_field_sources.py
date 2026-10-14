#!/usr/bin/env python
"""
Tie a band's astrometry to a reference band with FIELD SOURCES, then (optionally) re-cut its
stamp and carry its hand-drawn mask with the data.

Why (2026-10-14): `align_wfpc2_to_acs.py` ties bands with a deflector-light centroid. On SL2S
that leaves tens to hundreds of mas between bands (up to 231 mas, J0214 F814W), because the
reference is often F475X, where blue arcs pull the centroid. Transferred mask outlines then sit
pixels off the data. This script instead builds an INDEPENDENT catalogue of compact sources in each
full drizzled mosaic (S/N>15, each centroided on its own image with a self-re-centring window),
pairs them as nearest sky neighbours within 1.5", and takes the MODE of the pair offsets (a tight
cluster over a flat background of chance pairs), refined to a median within 60 mas of it.

**Two failure modes, both hit while writing this -- don't reintroduce them:**
- Predicting a source's position in the reference and centroiding there with a single-pass window
  follows the PREDICTION, not the source, so the measured offset ignores the WCS being tested. The
  first version reported 94 mas on J0214 where the truth is 231, and could not see its own
  correction.
- Matching without requiring a detection in the reference pairs blank sky with sky and drags the
  median to ~0.
- A FIXED tight match radius (0.3") drops every true pair once the offset exceeds it: J1406 F606W
  (0.35") came out as 31 +/- 111 mas of pure noise. Search wide, then take the mode.
Validate any change by shifting CRVAL by a known amount: the measurement must move by exactly that.

Apply (`--apply`) shifts CRVAL1/2 of every drizzled product of the target band, so that band's
sky positions move onto the reference frame. The original CRVALs are kept in FTORIG1/2, and the
tie is stamped in FTREF/FTDE/FTDN/FTNSRC and recorded in info/lens_field_ties.json. Re-running
re-measures the residual, which should be ~0 (idempotent); a tie is skipped if it is smaller than
`--min-shift` or under 3 sigma.

Re-cut (`--recut`) re-cuts the target band's 12" stamp about the REFERENCE stamp's central-pixel
sky position (`make_cutouts.py --center-from`), so the reference band's stamp, and every mask
drawn on it, is untouched. A mask that already exists on the target band was drawn on that band's
DATA. It therefore moves with the data, not with the sky: the new stamp is an integer-pixel
window of the same mosaic, so the shift is recovered exactly by matching the old and new sci
arrays, and the mask is rolled by it. Pixels shifted in from outside the old stamp become
unmasked (False) and are counted. Provenance goes to info/lens_masks.json under `field_tie_shift`.

    uv run python scripts/align_field_sources.py --sample sl2s_dinos_a --lens J0217-0513 --ref f475X           # measure only
    uv run python scripts/align_field_sources.py --sample sl2s_dinos_a --lens J0217-0513 --ref f475X --apply --recut

Run AFTER the drizzle and any `align_wfpc2_to_acs.py` tie, BEFORE make_cutouts. A re-drizzle
resets CRVAL, so re-run this afterwards.
"""
import argparse
import datetime
import glob
import os
import subprocess
import sys

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy import ndimage as ndi

WS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(WS, 'scripts'))
import info_json  # noqa: E402
from cutout_paths import find_prefix  # noqa: E402

TIES_JSON = os.path.join(WS, 'info', 'lens_field_ties.json')
MASKS_JSON = os.path.join(WS, 'info', 'lens_masks.json')


def mosaic(sample, lens, band):
    """CR-pass sci mosaic of a band (no-CR when that is all there is, e.g. WFC3/IR)."""
    d = os.path.join(WS, 'data', 'drizzled', sample, lens, band)
    for tag in ('_cr_', '_nocrrej_'):
        hit = sorted(p for p in glob.glob(os.path.join(d, '*_sci.fits')) if tag in os.path.basename(p))
        if hit:
            return hit[0]
    return None


def _snr(d):
    sm = ndi.gaussian_filter(d, 1.5)
    nz = sm[sm != 0]
    return (sm - np.median(nz)) / (1.4826 * np.median(np.abs(nz - np.median(nz))))


def _centroid(d, x, y, rad, niter=10, tol=0.01):
    """Iterative windowed centroid that RE-CENTRES on itself, so it converges onto the source
    and not onto the starting guess. (A single-pass window follows the guess -- that bug made
    the first version of this script blind to the very shift it was measuring.)"""
    for _ in range(niter):
        x0, y0 = int(round(x)), int(round(y))
        if not (rad + 1 < x0 < d.shape[1] - rad - 2 and rad + 1 < y0 < d.shape[0] - rad - 2):
            return None
        cut = d[y0 - rad:y0 + rad + 1, x0 - rad:x0 + rad + 1]
        c = np.clip(cut - np.median(cut), 0, None)
        if c.sum() <= 0:
            return None
        yy, xx = np.mgrid[-rad:rad + 1, -rad:rad + 1]
        nx, ny = x0 + (c * xx).sum() / c.sum(), y0 + (c * yy).sum() / c.sum()
        if np.hypot(nx - x, ny - y) < tol:
            return nx, ny
        x, y = nx, ny
    return x, y


def _catalogue(path, snr_min):
    """Compact sources detected and centroided on THIS image alone, as (ra, dec) arrays."""
    h = fits.getheader(path)
    d = np.nan_to_num(fits.getdata(path).astype(float))
    w = WCS(h).celestial
    sc = abs(w.pixel_scale_matrix[1, 1]) * 3600
    snr = _snr(d)
    pk = (snr == ndi.maximum_filter(snr, size=max(5, int(0.4 / sc)))) & (snr > snr_min)
    rad = max(2, int(round(0.1 / sc)))
    xy = []
    for y, x in zip(*np.nonzero(pk)):
        c = _centroid(d, x, y, rad)
        if c is not None and np.hypot(c[0] - x, c[1] - y) < rad:
            xy.append(c)
    xy = np.array(xy)
    sky = w.pixel_to_world(xy[:, 0], xy[:, 1])
    return sky.ra.deg, sky.dec.deg


def measure(ref_path, tgt_path, snr_min=15.0, search=1.5, tight=0.06):
    """Median (dE, dN) offset in mas of TARGET sky positions relative to REF.

    Sources are detected and centroided INDEPENDENTLY in each image. Every target source is
    paired with its nearest reference source within `search` arcsec, which is wide on purpose:
    a fixed tight match radius silently drops all true pairs once the offset exceeds it (0.3"
    missed J1406's 0.35" entirely and returned noise). The true offset is the MODE of the pair
    offsets -- a tight cluster above a flat background of chance pairs -- found as the peak of a
    smoothed 2-D histogram, then refined by keeping pairs within `tight` arcsec of it.
    Returns (dE, dN, err, n_kept, n_paired) or None."""
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    ra_r, de_r = _catalogue(ref_path, snr_min)
    ra_t, de_t = _catalogue(tgt_path, snr_min)
    cr = SkyCoord(ra_r * u.deg, de_r * u.deg)
    ct = SkyCoord(ra_t * u.deg, de_t * u.deg)
    idx, sep, _ = ct.match_to_catalog_sky(cr)
    ok = sep.arcsec < search
    if ok.sum() < 5:
        return None
    dE = (ra_t[ok] - ra_r[idx[ok]]) * np.cos(np.radians(de_r[idx[ok]])) * 3.6e6
    dN = (de_t[ok] - de_r[idx[ok]]) * 3.6e6
    lim, step = search * 1000, 20.0
    edges = np.arange(-lim, lim + step, step)
    H, _, _ = np.histogram2d(dE, dN, bins=[edges, edges])
    H = ndi.gaussian_filter(H, 1.5)
    ie, jn = np.unravel_index(np.argmax(H), H.shape)
    med = np.array([edges[ie] + step / 2, edges[jn] + step / 2])
    a = np.c_[dE, dN]
    for _ in range(6):
        keep = np.hypot(*(a - med).T) < tight * 1000
        if keep.sum() < 5:
            return None
        med = np.median(a[keep], 0)
    mad = 1.4826 * np.median(np.abs(a[keep] - med), 0)
    err = np.hypot(*mad) / np.sqrt(keep.sum()) * 1.2533
    return float(med[0]), float(med[1]), float(err), int(keep.sum()), int(ok.sum())


def apply_tie(sample, lens, band, dE, dN, ref, n_src, err):
    """Shift CRVAL of every drizzled product of `band` by -(dE, dN), keeping the originals."""
    d = os.path.join(WS, 'data', 'drizzled', sample, lens, band)
    n = 0
    for fn in sorted(glob.glob(os.path.join(d, '*.fits'))):
        with fits.open(fn, mode='update') as h:
            hd = h[0].header
            if 'CRVAL1' not in hd:
                continue
            if 'FTORIG1' not in hd:
                hd['FTORIG1'] = (hd['CRVAL1'], 'CRVAL1 before the field-source tie')
                hd['FTORIG2'] = (hd['CRVAL2'], 'CRVAL2 before the field-source tie')
            hd['CRVAL1'] -= dE / 3.6e6 / np.cos(np.radians(hd['CRVAL2']))
            hd['CRVAL2'] -= dN / 3.6e6
            hd['FTREF'] = (ref, 'band this product is field-source tied to')
            hd['FTDE'] = (hd.get('FTDE', 0.0) + dE, '[mas] cumulative East shift removed')
            hd['FTDN'] = (hd.get('FTDN', 0.0) + dN, '[mas] cumulative North shift removed')
            hd['FTNSRC'] = (n_src, 'field sources in the last tie')
            hd['ASTROREF'] = (ref, 'band this product was astrometrically tied to')
            h.flush()
        n += 1
    return n


def integer_shift(old, new, maxs=8):
    """(dy, dx) with new[y, x] == old[y + dy, x + dx] over the overlap (exact pixel match)."""
    best = None
    for dy in range(-maxs, maxs + 1):
        for dx in range(-maxs, maxs + 1):
            ys, ye = max(0, -dy), min(new.shape[0], old.shape[0] - dy)
            xs, xe = max(0, -dx), min(new.shape[1], old.shape[1] - dx)
            a = new[ys:ye, xs:xe]
            b = old[ys + dy:ye + dy, xs + dx:xe + dx]
            frac = np.mean(a == b)
            if best is None or frac > best[0]:
                best = (frac, dy, dx)
    return best


def roll_mask(mask, dy, dx):
    """new[y, x] = mask[y + dy, x + dx]. Pixels with no counterpart (the strip shifted in at the
    stamp edge) copy the nearest edge pixel, so a mask that runs off the stamp keeps running off
    it instead of leaving an unmasked sliver."""
    ny, nx = mask.shape
    pad = max(abs(dy), abs(dx))
    big = np.pad(mask, pad, mode='edge')
    return big[pad + dy:pad + dy + ny, pad + dx:pad + dx + nx]


def reset_tie(sample, lens, band):
    """Restore every product of `band` to its pre-tie CRVAL (FTORIG1/2) and drop the FT* cards."""
    d = os.path.join(WS, 'data', 'drizzled', sample, lens, band)
    n = 0
    for fn in sorted(glob.glob(os.path.join(d, '*.fits'))):
        with fits.open(fn, mode='update') as h:
            hd = h[0].header
            if 'FTORIG1' not in hd:
                continue
            hd['CRVAL1'], hd['CRVAL2'] = hd['FTORIG1'], hd['FTORIG2']
            for k in ('FTORIG1', 'FTORIG2', 'FTREF', 'FTDE', 'FTDN', 'FTNSRC'):
                if k in hd:
                    del hd[k]
            h.flush()
        n += 1
    return n


def recut(sample, lens, band, ref):
    cdir = os.path.join(WS, 'data', 'cutouts', sample, lens, band)
    rdir = os.path.join(WS, 'data', 'cutouts', sample, lens, ref)
    prefix, rprefix = find_prefix(cdir), find_prefix(rdir)
    old_sci = fits.getdata(os.path.join(cdir, f'{prefix}_sci.fits')).astype(float)
    mpath = os.path.join(cdir, f'{prefix}_mask.fits')
    old_mask = fits.getdata(mpath).astype(bool) if os.path.exists(mpath) else None
    subprocess.run(['uv', 'run', 'python', os.path.join(WS, 'scripts', 'make_cutouts.py'),
                    '--sample', sample, '--lens', lens, '--filt', band,
                    '--center-from', os.path.join(rdir, f'{rprefix}_sci.fits')],
                   check=True, cwd=WS, stdout=subprocess.DEVNULL)
    new_sci = fits.getdata(os.path.join(cdir, f'{prefix}_sci.fits')).astype(float)
    frac, dy, dx = integer_shift(old_sci, new_sci)
    msg = f'    re-cut: data moved by (drow, dcol)=({-dy:+d},{-dx:+d}) px in the new stamp (exact match {frac*100:.1f}%)'
    if old_mask is None:
        print(msg + '; no mask on this band')
        return
    if frac < 0.95:
        print(msg + '\n    !! old/new sci do not match by an integer shift -- mask NOT moved; check by hand')
        return
    new_mask = roll_mask(old_mask, dy, dx)
    lost = int(old_mask.sum() - new_mask.sum())   # net change; edge strips are filled from the edge
    hl = fits.open(mpath)
    hl[0].data = new_mask.astype(hl[0].data.dtype)
    hl[0].header['HISTORY'] = (f'{datetime.date.today()} field-source re-tie: mask rolled with the data '
                               f'by (drow,dcol)=({-dy:+d},{-dx:+d}) px')
    hl.writeto(mpath, overwrite=True)
    e = info_json.entry(MASKS_JSON, sample, lens, band) or {}
    e['field_tie_shift'] = {'date': str(datetime.date.today()), 'drow_px': -dy, 'dcol_px': -dx,
                            'net_masked_px_change': -lost, 'ref': ref}
    e['n_excluded_px'] = int(new_mask.sum())
    info_json.update(MASKS_JSON, sample, lens, band, dict(sorted(e.items())))
    print(msg + f'; mask rolled with it (net masked-px change {-lost:+d})')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--sample', required=True)
    p.add_argument('--lens', required=True)
    p.add_argument('--ref', default='f475X')
    p.add_argument('--target', nargs='*', help='bands to tie (default: every other band with a mosaic)')
    p.add_argument('--apply', action='store_true')
    p.add_argument('--recut', action='store_true', help='re-cut the tied band stamps and move their masks')
    p.add_argument('--min-shift', type=float, default=10.0, help='[mas] skip ties smaller than this')
    p.add_argument('--reset', action='store_true', help='restore the pre-tie CRVAL of the target bands and exit')
    a = p.parse_args()
    ref_path = mosaic(a.sample, a.lens, a.ref)
    if ref_path is None:
        sys.exit(f'{a.lens}: no {a.ref} mosaic')
    bands = a.target or sorted(b for b in os.listdir(os.path.join(WS, 'data', 'drizzled', a.sample, a.lens))
                               if b != a.ref and mosaic(a.sample, a.lens, b))
    if a.reset:
        for b in bands:
            print(f'{a.lens} {b}: reset {reset_tie(a.sample, a.lens, b)} files to their pre-tie CRVAL')
        return
    for b in bands:
        tp = mosaic(a.sample, a.lens, b)
        r = measure(ref_path, tp)
        if r is None:
            print(f'{a.lens} {b}: too few common field sources -- not tied')
            continue
        dE, dN, err, nk, nall = r
        d = np.hypot(dE, dN)
        do = a.apply and d >= a.min_shift and d >= 3 * err
        print(f'{a.lens} {b} vs {a.ref}: dE={dE:+.0f} dN={dN:+.0f} |d|={d:.0f} mas (+/-{err:.0f}, {nk}/{nall} sources)'
              + ('  -> APPLY' if do else ('  (below threshold, not applied)' if a.apply else '')))
        if do:
            n = apply_tie(a.sample, a.lens, b, dE, dN, a.ref, nk, err)
            r2 = measure(ref_path, tp)
            resid = f'{np.hypot(r2[0], r2[1]):.0f} mas' if r2 else 'n/a'
            rec = info_json.entry(TIES_JSON, a.sample, a.lens, b) or {}
            hist = rec.get('history', [])
            hist.append({'date': str(datetime.date.today()), 'ref': a.ref, 'dE_mas': round(dE, 1), 'dN_mas': round(dN, 1),
                         'err_mas': round(err, 1), 'n_sources': nk, 'files': n, 'residual_after': resid})
            info_json.update(TIES_JSON, a.sample, a.lens, b, {'ref': a.ref, 'history': hist})
            print(f'    shifted CRVAL of {n} files; residual after: {resid}')
            if a.recut and os.path.isdir(os.path.join(WS, 'data', 'cutouts', a.sample, a.lens, b)):
                recut(a.sample, a.lens, b, a.ref)


if __name__ == '__main__':
    main()
