#!/usr/bin/env python
"""
Measure the WCS offset between the VISITS of a WFPC2 band, from compact sources in the
distortion-corrected input frames, to decide whether the visits may be drizzled together.

Why (2026-10-29): roll is not the criterion. Two visits can be combined only when their
delivered WCSs agree, and that depends on the guide-star solution, not on PA_V3. J1318-0313's
visits are 0.97 deg apart in roll and 0.37" apart on the sky; J1719+2939's are 1.7 deg apart and
agree to 0.01" (both Gaia-fitted); J1142+1001's share a roll to 0.13 deg and are 0.21" apart.
**Combine only below ~0.03"**, the frame-to-frame scatter the WFPC2 core-registration audit
accepted; otherwise drizzle per visit (`run_wfpc2_wf3.sh` `split_visits_for`).

Method: the frames the drizzle actually used (`data/drizzle_files/<sample>/<band_mode>/<category>/<lens>/<band>*/
wf3_*_flt.fits`, after `updatewcs`, so the SIP distortion is in the WCS) are grouped by visit
(the first six characters of the rootname). In each frame, compact sources are found with
DAOStarFinder (8 sigma, CR-like sharpness rejected) away from the vignetted chip border and
put on the sky with that frame's own WCS. Every source pair within 1.5" between two frames gives
an offset; the MODE of those offsets (a tight cluster over chance pairs) is refined to a median
within 0.12". A cosmic ray is in one frame only, so it never forms a cluster. Each visit pair
is measured on the frame pair with the most matches, and a within-visit pair is printed as the
control (it should read ~0.01"). Validated against the deflector-centroid tie
(`align_wfpc2_to_acs.py`) on the 2026-10-29 split lenses: J1318-0313 0.368" vs 0.366",
J1306+0600 0.484" vs 0.477", J1436+3640 0.258" vs 0.255". Read-only; writes nothing.

Usage:
    uv run python scripts/measure_visit_offset.py --sample slacs --lens J1142+1001
    uv run python scripts/measure_visit_offset.py --sample slacs --lens J0822+2652 --filt f606W
"""
import argparse
import glob
import itertools
import os
import sys
import warnings

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS, FITSFixedWarning
from photutils.detection import DAOStarFinder

warnings.simplefilter('ignore', FITSFixedWarning)
ws_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ws_path, 'scripts'))
import mast_target_names  # noqa: E402
import stamp_paths  # noqa: E402

BORDER = (60, 20)        # px trimmed from the low / high edges (WF3's vignetted L-shaped border)
MATCH_RADIUS = 1.5       # arcsec: search wide, then take the mode (see align_field_sources.py)
COMBINE_MAX = 0.03       # arcsec: the largest visit-to-visit offset that may be drizzled together


def frame_sources(path):
    """Sky positions (deg) of compact sources in a frame's WF3 SCI extension."""
    with fits.open(path) as h:
        sci = [e for e in h if e.header.get('EXTNAME') == 'SCI']
        ext = next((e for e in sci if e.header.get('DETECTOR') == 3), sci[0])
        data = np.nan_to_num(ext.data.astype(float))
        wcs = WCS(ext.header, h)
    lo, hi = BORDER
    d = data[lo:-hi, lo:-hi]
    _, med, std = sigma_clipped_stats(d, sigma=3)
    tab = DAOStarFinder(fwhm=1.6, threshold=8 * std, sharplo=0.3, sharphi=0.9,
                        roundlo=-0.6, roundhi=0.6)(d - med)
    if tab is None:
        return np.empty((0, 2))
    ra, dec = wcs.all_pix2world(tab['xcentroid'] + lo, tab['ycentroid'] + lo, 0)
    return np.c_[ra, dec]


def pair_offset(a, b):
    """(dE, dN) arcsec of frame b relative to frame a, number of matched pairs, scatter."""
    if len(a) == 0 or len(b) == 0:
        return None
    cosd = np.cos(np.deg2rad(np.median(a[:, 1])))
    d = ((b[None, :, :] - a[:, None, :]) * [cosd * 3600, 3600]).reshape(-1, 2)
    d = d[(np.abs(d) < MATCH_RADIUS).all(axis=1)]
    if len(d) < 3:
        return None
    hist, xe, ye = np.histogram2d(d[:, 0], d[:, 1], bins=60,
                                  range=[[-MATCH_RADIUS, MATCH_RADIUS]] * 2)
    i, j = np.unravel_index(hist.argmax(), hist.shape)
    peak = np.array([(xe[i] + xe[i + 1]) / 2, (ye[j] + ye[j + 1]) / 2])
    for _ in range(5):
        sel = np.hypot(*(d - peak).T) < 0.12
        peak = np.median(d[sel], axis=0)
    return peak, int(sel.sum()), float(np.std(np.hypot(*(d[sel] - peak).T)))


def best_pair(frames_a, frames_b, cat):
    """The frame pair with the most matched sources, and its offset."""
    best = None
    for fa, fb in itertools.product(frames_a, frames_b):
        if fa == fb:
            continue
        r = pair_offset(cat[fa], cat[fb])
        if r and (best is None or r[1] > best[1][1]):
            best = ((fa, fb), r)
    return best


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--sample', default=mast_target_names.DEFAULT_SAMPLE)
    p.add_argument('--lens', required=True)
    p.add_argument('--filt', default='f606W',
                   help='band; every <filt>* work dir is read, so a split lens contributes all its visits')
    a = p.parse_args()

    frames = sorted(glob.glob(os.path.join(ws_path, 'data', 'drizzle_files', stamp_paths.lens_rel(a.sample, a.lens),
                                           a.filt + '*', 'wf3_*_flt.fits')))
    if not frames:
        sys.exit(f'no wf3_*_flt.fits for {a.sample}/{a.lens}/{a.filt}* in data/drizzle_files/ '
                 f'(re-run the drizzle; the work dir is gitignored)')
    visits = {}
    for f in frames:
        visits.setdefault(os.path.basename(f)[4:10], []).append(f)
    cat = {f: frame_sources(f) for f in frames}
    print(f'{a.lens} {a.filt}: visits {", ".join(f"{v} ({len(fs)} frames)" for v, fs in visits.items())}')

    for v, fs in visits.items():
        r = best_pair(fs, fs, cat) if len(fs) > 1 else None
        if r:
            (dx, dy), n, sc = r[1]
            print(f'  control {v} within visit      {np.hypot(dx, dy):.3f}"  n={n}  scatter {sc:.3f}"')
    for (va, fa), (vb, fb) in itertools.combinations(visits.items(), 2):
        r = best_pair(fa, fb, cat)
        if r is None:
            print(f'  {va} vs {vb}: too few matched sources')
            continue
        (dx, dy), n, sc = r[1]
        off = np.hypot(dx, dy)
        verdict = 'may combine' if off < COMBINE_MAX else 'drizzle per visit'
        print(f'  {va} vs {vb}  offset {off:.3f}" (dE {dx:+.3f}, dN {dy:+.3f})  n={n}  '
              f'scatter {sc:.3f}"  -> {verdict}')


if __name__ == '__main__':
    main()
