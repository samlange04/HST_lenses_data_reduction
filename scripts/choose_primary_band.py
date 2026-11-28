#!/usr/bin/env python
"""
Merged 2026-11-28 across the HST and JWST repos: the file is identical in both (the repo
profile is the import seam below).

Choose each lens's PRIMARY band: the one where a lensed image's position is best defined. The
mask and arc-mask GUIs draw it first (every other band starts from its mask), positions are
clicked on it, and it is the band to model first.

    uv run python scripts/choose_primary_band.py                    # every sample
    uv run python scripts/choose_primary_band.py --sample swells [--lens LENS]

Ported from JWST_lenses_data_reduction 2026-11-10. A fixed order cannot do this. The HST order
(f475X > f814W > f606W > ...) is a sensible default for SLACS and gallery, where the deep red
band is always the best exposed, but it need not hold per lens elsewhere: swells disks are
bright and extended in every band, and the arc's S/N depends on the source's colour and on how
well the band's PSF resolves the Einstein ring.

Per band, on the modelling-size stamp:
  * the lens light is the azimuthal median profile of the stamp, interpolated in log-log
    across the ring annulus |r - theta_E| < w (w = max(0.25 theta_E, the PSF FWHM, 0.15"));
  * arc S/N = the SIGNED sum of the residual in that annulus over its noise in quadrature,
    with the band's contaminant mask applied if it has one. Signed, so that the quadrupole an
    elliptical lens leaves about a circular profile cancels;
  * score = arc S/N x min(1, theta_E / (R0 x FWHM)), R0 = 5: a clicked position's precision
    goes as FWHM / S/N while the ring is under-resolved; once theta_E spans R0 PSF widths,
    more resolution no longer helps, and S/N alone decides.
The suggestion is the best score. None if no band reaches MIN_ARC_SNR: the arcs are then too
faint against the error of the lens-light profile.

So the measurement (`best_band_arc_snr`) is a guide, not a verdict. The band chosen by eye with
choose_best_band.py is `best_band_visual` in info/samples/lens_samples.json and always wins.
mast_target_names.primary_band() resolves visual > arc_snr > none (the GUIs then fall back to
their fixed order, make_masks._BAND_PRIORITY).

REVIEWED vs NOT. Every record carries `reviewed` (bool): an arc_snr suggestion is born
unreviewed, and choose_best_band.py marks it reviewed when the reviewer accepts the suggested
band (accepting a DIFFERENT band writes best_band_visual instead -- source 'visual', inherently
reviewed). record() preserves a reviewed flag only while the resolved primary band is unchanged;
a re-measurement that moves the suggestion drops it back to unreviewed. `reviewed_date` /
`review_note` say when and by what.

theta_E comes from info/catalogues/lens_einstein_radii.json (measured values only); a lens with
a null theta_E (no published lens model) is skipped -- the annulus cannot be placed.

Writes:
  * info/provenance/lens_primary_band.json: per lens, best_band_arc_snr, best_band_visual, the
    resolved primary_band and its source, theta_E, and every band's arc S/N, FWHM and score;
  * data/stamps/<sample>/<band_mode>/<category>/<lens>/primary_band.json: the resolved choice,
    so it travels with the lens's stamps.
"""
import argparse
import glob
import json
import os
import sys
import warnings

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import info_json
import make_band_overview
import stamp_paths
try:
    import lens_samples as mast_target_names   # the JWST repo's profile module
except ImportError:
    import mast_target_names                   # the HST repo's

ws_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

R0 = 5.0               # theta_E >= R0 PSF FWHMs: resolution no longer limits the position
MIN_ARC_SNR = 15.0     # below this in every band the measurement is the lens-profile error
# Fallback PSF FWHM (arcsec) when the band has no lens_psf.json record. The per-filter table
# (the JWST repo's, from the JDox nominal values) is tried first, then the per-instrument one
# (the HST repo's: its filters recur across cameras, so key on the instrument where the filter
# is ambiguous -- a WFPC2 f814W must not read the ACS value). HST values are the drizzled
# kernels' typical FWHM: ACS/WFC + WFC3/UVIS ~2 px at 0.05/0.0396", WFPC2/WF3 ~4.5 px at
# 0.05", WFC3/IR ~3.45 px at 0.06".
JDOX_FWHM = {'F070W': 0.030, 'F090W': 0.034, 'F115W': 0.040, 'F150W': 0.050, 'F200W': 0.066,
             'F277W': 0.092, 'F356W': 0.116, 'F444W': 0.145}
FALLBACK_FWHM = {'ACS': 0.10, 'WFC3_UVIS': 0.085, 'WFPC2': 0.23, 'WFC3_IR': 0.21,
                 'NIRCAM_SW': 0.045, 'NIRCAM_LW': 0.12}


def _instrument_key(hdr):
    inst = str(hdr.get('INSTRUME', '')).strip().upper()
    det = str(hdr.get('DETECTOR', '')).strip().upper()
    if inst == 'WFPC2':
        return 'WFPC2'
    if inst == 'ACS':
        return 'ACS'
    if inst == 'WFC3':
        return 'WFC3_IR' if det == 'IR' else 'WFC3_UVIS'
    if inst == 'NIRCAM':
        from make_band_overview import FILTER_UM
        um = FILTER_UM.get(str(hdr.get('FILTER', '')).strip(), 99.0)
        return 'NIRCAM_SW' if um < 2.4 else 'NIRCAM_LW'
    return 'ACS'


def arc_snr(sci, noise, scale, theta_e, fwhm, mask=None):
    """Signed arc S/N in the ring annulus after subtracting the interpolated lens profile."""
    n = sci.shape[0]
    yy, xx = np.mgrid[:n, :n]
    r = np.hypot(xx - (n - 1) / 2, yy - (n - 1) / 2) * scale
    good = np.isfinite(sci) & (noise > 0) & (noise < 1e7)   # 1e8 is the no-coverage sentinel
    if mask is not None:
        good &= ~mask
    w = max(0.25 * theta_e, fwhm, 0.15)
    edges = np.arange(0.0, r.max() + scale, scale)
    rc = 0.5 * (edges[1:] + edges[:-1])
    prof = np.full(len(rc), np.nan)
    for k, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
        sel = good & (r >= a) & (r < b)
        if sel.any():
            prof[k] = np.median(sci[sel])
    ok = np.isfinite(prof) & (prof > 0) & (np.abs(rc - theta_e) >= w)
    if ok.sum() < 3:
        return None
    lens = np.exp(np.interp(np.log(np.maximum(r, scale / 2)), np.log(rc[ok]), np.log(prof[ok])))
    ann = good & (np.abs(r - theta_e) < w)
    if not ann.any():
        return None
    return float((sci - lens)[ann].sum() / np.sqrt((noise[ann] ** 2).sum()))


def measure(sample, lens):
    """Per-band arc S/N and resolution-discounted score for `lens`, and the suggested band."""
    te = mast_target_names.theta_e(lens, sample)
    psf = info_json.load(info_json.path(ws_path, 'lens_psf.json')).get(sample, {}).get(lens, {})
    prefix = stamp_paths.stamp_prefix(stamp_paths.modelling_size(lens, sample))
    lens_dir = stamp_paths.lens_dir(ws_path, 'stamps', sample, lens)
    scores = {}
    for d in sorted(glob.glob(os.path.join(lens_dir, '*'))):
        if not os.path.isdir(d) or not os.path.basename(d).startswith(('f', 'F')):
            continue
        filt = os.path.basename(d)
        sci_f = stamp_paths.find_stamp(d, prefix, 'sci')
        noise_f = stamp_paths.find_stamp(d, prefix, 'noise')
        if not sci_f or not noise_f or te is None:
            continue
        hdr = fits.getheader(sci_f)
        scale = float(proj_plane_pixel_scales(WCS(hdr).celestial)[0]) * 3600.0
        sci = fits.getdata(sci_f).astype(float)
        noise = fits.getdata(noise_f).astype(float)
        mask_f = os.path.join(d, f'{prefix}_mask.fits')
        mask = fits.getdata(mask_f).astype(bool) if os.path.exists(mask_f) else None
        rec_psf = psf.get(filt) or psf.get(filt.split('_')[0]) or {}
        # lens_psf.json records fwhm_pix in the HST repo, fwhm_arcsec in the JWST one
        fwhm = ((rec_psf.get('fwhm_pix') or 0) * scale
                or rec_psf.get('fwhm_arcsec')
                or JDOX_FWHM.get(filt.upper().split('_')[0], 0)
                or FALLBACK_FWHM[_instrument_key(hdr)])
        snr = arc_snr(sci, noise, scale, te, fwhm, mask)
        res = min(1.0, te / (R0 * fwhm))
        scores[filt] = dict(arc_snr=None if snr is None else round(snr, 2),
                            psf_fwhm_arcsec=round(fwhm, 4),
                            resolution_factor=round(res, 3),
                            score=None if snr is None else round(snr * res, 2),
                            contaminant_mask=mask is not None)
    valid = {f: v for f, v in scores.items() if v['score'] is not None}
    best = max(valid, key=lambda f: valid[f]['score']) if valid else None
    measurable = bool(valid) and max(v['arc_snr'] for v in valid.values()) >= MIN_ARC_SNR
    return dict(theta_e=te, scores=scores, best_band_arc_snr=best if measurable else None,
                measurable=measurable, min_arc_snr=MIN_ARC_SNR, r0=R0)


def record(sample, lens, rerender=True):
    """Measure the lens, resolve its primary band (visual > arc_snr), and write both records
    (lens_primary_band.json, the lens stamps dir's primary_band.json) and, with `rerender`,
    the band overview. Returns the record, or None without stamps or theta_E.

    The `reviewed` flag survives a re-measurement only while the resolved primary band is
    unchanged; a moved suggestion drops back to unreviewed."""
    rec = measure(sample, lens)
    if not rec['scores']:
        return None
    out_json = info_json.path(ws_path, 'lens_primary_band.json')
    existing = info_json.load(out_json).get(sample, {}).get(lens, {})   # BEFORE the writes
    for k, v in rec.items():                   # the suggestion first, so primary_band() sees it
        info_json.update(out_json, sample, lens, k, v)
    band, source = mast_target_names.primary_band(lens, sample)
    reviewed = (bool(existing.get('reviewed')) and existing.get('primary_band') == band
                and existing.get('primary_band_source') == source
                and existing.get('best_band_arc_snr') == rec['best_band_arc_snr'])
    info_json.update(out_json, sample, lens, 'primary_band', band)
    info_json.update(out_json, sample, lens, 'primary_band_source', source)
    info_json.update(out_json, sample, lens, 'reviewed', reviewed)
    visual = mast_target_names.entry(lens, sample).get('best_band_visual')
    info_json.update(out_json, sample, lens, 'best_band_visual', visual)
    rec.update(primary_band=band, primary_band_source=source, reviewed=reviewed)
    stamp_rec = dict(lens=lens, sample=sample, primary_band=band, primary_band_source=source,
                     reviewed=reviewed,
                     best_band_visual=visual,
                     best_band_arc_snr=rec['best_band_arc_snr'], theta_e=rec['theta_e'],
                     scores=rec['scores'],
                     how='best_band_arc_snr: choose_primary_band.py (measured); '
                         'best_band_visual: choose_best_band.py (by eye). reviewed: the '
                         'arc_snr choice has been confirmed by eye (visual choices are '
                         'reviewed by construction).')
    lens_dir = stamp_paths.lens_dir(ws_path, 'stamps', sample, lens)
    with open(os.path.join(lens_dir, 'primary_band.json'), 'w') as f:
        json.dump(stamp_rec, f, indent=4)
    if rerender:
        make_band_overview.rerender(sample, lens)          # its star marks the primary band
    return rec


def set_reviewed(sample, lens, reviewed=True, note=None):
    """Mark the lens's arc_snr primary-band choice reviewed (choose_best_band.py accepting the
    suggested band). Writes lens_primary_band.json and the per-lens primary_band.json."""
    import datetime
    out_json = info_json.path(ws_path, 'lens_primary_band.json')
    rec = info_json.load(out_json).get(sample, {}).get(lens, {})
    if not rec:
        print(f'{lens}: no primary-band record to mark reviewed -- run choose_primary_band.py first')
        return
    info_json.update(out_json, sample, lens, 'reviewed', bool(reviewed))
    info_json.update(out_json, sample, lens, 'reviewed_date', datetime.date.today().isoformat())
    if note:
        info_json.update(out_json, sample, lens, 'review_note', note)
    lens_dir = stamp_paths.lens_dir(ws_path, 'stamps', sample, lens)
    p = os.path.join(lens_dir, 'primary_band.json')
    try:
        stamp_rec = json.load(open(p))
    except FileNotFoundError:
        return
    stamp_rec['reviewed'] = bool(reviewed)
    stamp_rec['reviewed_date'] = datetime.date.today().isoformat()
    if note:
        stamp_rec['review_note'] = note
    with open(p, 'w') as f:
        json.dump(stamp_rec, f, indent=4)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--sample', type=mast_target_names.sample_arg, default=None,
                   help='default: every sample')
    p.add_argument('--lens', default=None)
    a = p.parse_args()
    mast_target_names.resolve_sample(p, a)          # validates the lens, fills a.sample
    warnings.simplefilter('ignore')
    samples = [a.sample] if a.sample else mast_target_names.samples()
    for sample in samples:
        for lens in mast_target_names.lenses(sample):
            if a.lens and lens != a.lens:
                continue
            rec = record(sample, lens)
            if rec is None:
                print(f'{lens}: no stamps (or no theta_e), skipped')
                continue
            band, source = rec['primary_band'], rec['primary_band_source']
            row = '  '.join(f'{f} {v["arc_snr"]}/{v["score"]}' for f, v in rec['scores'].items())
            flag = ('' if rec['measurable'] or source == 'visual'
                    else f'  [no band reaches arc S/N {MIN_ARC_SNR:g}: set primary_band by eye]')
            print(f'{lens:18s} theta_E {rec["theta_e"]:.2f}"  arc-S/N best {rec["best_band_arc_snr"]}  -> primary {band} ({source})'
                  f'{flag}\n    arc S/N/score: {row}')
    print(f"wrote {info_json.path(ws_path, 'lens_primary_band.json')} and "
          f"data/stamps/<sample>/<band_mode>/<category>/<lens>/primary_band.json")


if __name__ == '__main__':
    sys.exit(main())
