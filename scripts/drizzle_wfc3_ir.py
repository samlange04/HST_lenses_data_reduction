#!/usr/bin/env python
"""
Drizzle HST WFC3/IR FLT images into a combined mosaic.
Downloads FLT files from MAST into data/calibrated/ if not already present.
Writes final sci/wht FITS to data/drizzled/, intermediates to data/drizzle_files/.
Produces a no-CR-rejection drizzle by default. Pass --cr to also run the
CR-rejected pass.

FLT files are the recommended input for WFC3/IR: the IR detector is not a CCD
so there is no CTE effect and no FLC (CTE-corrected) variant is produced.
"""

import argparse
import json
import os
import glob
import shutil
import subprocess
import sys
import numpy as np
from scipy.ndimage import median_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import LogStretch, ImageNormalize
from astroquery.mast import Observations
from drizzlepac import tweakreg, astrodrizzle
from stwcs.updatewcs import updatewcs

# Route large FITS output writes through mmap+memcpy (vm_fault path) instead of
# fwrite/cluster_write copyin, to dodge the macOS U-state write-path lost-wakeup.
# Must run before AstroDrizzle. Import from this script's own dir regardless of cwd.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mmap_fits_write
mmap_fits_write.install()
import mast_target_names

# ── Configuration ──────────────────────────────────────────────────────────────
_p = argparse.ArgumentParser()
_p.add_argument('--lens',        default='J0008-0004')
_p.add_argument('--filt',        default='f160W')
_p.add_argument('--sample',      default='slacs')
_p.add_argument('--align',       default='mast', choices=['mast', 'tweakreg'],
                help="'mast' (default) trusts the GSC242/GAIAeDR3-fitted WCS in the "
                     "delivered files and runs neither updatewcs nor TweakReg; "
                     "'tweakreg' restores the old re-solve, which erases the dither.")
_p.add_argument('--cr',          action='store_true', default=False)
_p.add_argument('--_subprocess', action='store_true', default=False, help=argparse.SUPPRESS)
# Drizzle output weight type. 'ERR' (default) makes the WHT extension a full
# inverse-variance map (source Poisson + sky + read + dark), so 1/sqrt(WHT) is a
# CALIBRATED per-pixel noise map -- what make_cutouts uses. AstroDrizzle's own
# default 'EXP' is only an effective-exposure-time map (uncalibrated, missing
# source shot noise). See DrizzlePac Handbook pp.103,139 and Bayer et al. 2023.
_p.add_argument('--wht-type',    default='ERR', choices=['ERR', 'IVM', 'EXP'])
_p.add_argument('--cr-method',   default='drizcr', choices=['lacosmic', 'drizcr'],
                help="How the --cr pass rejects cosmic rays. Neither is recommended for "
                     "WFC3/IR -- see the note above run_lacosmic(). 'drizcr' (default) is "
                     "AstroDrizzle's median/blot route; 'lacosmic' flags per frame with "
                     "astroscrappy into DQ 4096, and destroys point sources on this detector.")
_p.add_argument('--lacosmic-sigclip', type=float, default=4.5)
_p.add_argument('--lacosmic-objlim',  type=float, default=5.0)
_p.add_argument('--dq-refine', type=float, default=3.0, metavar='SIGMA',
                help="Un-flag DQ 8/16/32 pixels that are not actually deviant in their "
                     "own exposure, at this sigma (default 3.0; 0 disables). The dark "
                     "reference file flags these over a whole anneal cycle, so they are "
                     "conservative: measured on J0841+3824 only 35-43%% of them deviate "
                     "from a local median by >3 sigma, against 1.2%% of unflagged pixels. "
                     "Masking all of them costs 26%% of the mosaic's coverage -- see "
                     "refine_dq_flags().")
_a = _p.parse_args()

lens           = _a.lens
sample         = _a.sample
filt           = _a.filt
do_cr          = _a.cr
_is_subprocess = _a._subprocess

ws_path     = '/Users/samlange/Code/data_reduction'
data_path   = os.path.join(ws_path, 'data', 'calibrated', sample, lens, filt)
output_path = os.path.join(ws_path, 'data', 'drizzled', sample, lens, filt)
work_path   = os.path.join(ws_path, 'data', 'drizzle_files', sample, lens, filt)
ref_path    = os.path.join(ws_path, 'data', 'reference_files')

# ── Common output WCS across filters (orientation + centre) ────────────────────
# Each band was drizzled at its native, differing roll (e.g. J0252+0039: F606W 136deg,
# F814W -75deg, F160W -144deg) with its own tangent point, so filters of one lens did
# not overlay. Pin every band to North-up (final_rot=0) with the tangent point
# (final_ra/dec) at the lens position, so all filters share orientation and WCS centre
# and the lens sits ~centred. Pixel scale still differs by instrument (intentional);
# only rotation and centre are unified. Unknown lens -> native WCS (empty dict).
sys.path.insert(0, os.path.join(ws_path, 'info'))
try:
    from slacs_coords import slacs_coords as _slacs_coords
    from astropy.coordinates import SkyCoord as _SkyCoord
    import astropy.units as _u
    _lc = _SkyCoord(*_slacs_coords[lens], unit=(_u.hourangle, _u.deg))
    _common_wcs = dict(final_rot=0.0, final_ra=float(_lc.ra.deg), final_dec=float(_lc.dec.deg))
    print(f'=== Common WCS: North-up, centre ({_lc.ra.deg:.5f}, {_lc.dec.deg:.5f}) ===')
except (KeyError, ImportError):
    _common_wcs = {}
    print('=== Common WCS: lens not in slacs_coords -> native drizzle WCS ===')


# ── Info JSON paths ────────────────────────────────────────────────────────────
filt_key             = filt
json_path            = os.path.join(ws_path, 'info', 'lens_products.json')
exptime_json_path    = os.path.join(ws_path, 'info', 'lens_exptime.json')
instrument_json_path = os.path.join(ws_path, 'info', 'lens_instrument.json')

def _update_info_json(path, lens, filt_key, value):
    try:
        with open(path) as _f:
            _data = json.load(_f)
    except FileNotFoundError:
        _data = {}
    _entry = _data.setdefault(lens, {})
    _entry[filt_key] = value
    # Keep filters ordered within each lens entry, as well as lenses across the file.
    _data[lens] = dict(sorted(_entry.items()))
    with open(path, 'w') as _f:
        json.dump(dict(sorted(_data.items())), _f, indent=4)

# ── Output sampling ────────────────────────────────────────────────────────────
# WFC3/IR native is 0.1283"/px. The WFC3-IR-DITHER-BOX-MIN pattern gives 4 distinct
# sub-pixel phases on both axes, so the data support a sub-native grid; drizzling at
# native throws that sampling away. 0.06"/px is just under half-native and is the
# standard WFC3/IR choice (CANDELS, 3D-HST). F160W stays on 0.06" while F606W/F814W
# are on 0.05" -- do NOT re-drizzle to 0.05" to grid-match: PyAutoLens ingests each
# band at its native scale and pixel-matches at the modelling stage (verified in the
# pyauto multi-wavelength API, which even fits sub-pixel inter-band grid offsets), and
# 0.05" only worsens the weight non-uniformity (8-19% under-covered across the 13
# lenses; J1430/J1029/J0841 worst) with no resolution gain.
#
# pixfrac 1.0 (not 0.8): with the noise now on calibrated ERR weight maps, the goal is
# a *uniform, low-correlation* noise map for the likelihood, not the sharpest PSF.
# pixfrac 1.0 fills the drizzle grid so adjacent-pixel noise correlation collapses.
# Measured on J0252+0039 (0.06", correctly registered frames):
#   pixfrac   noise texture   adjacent-pixel RMS
#   0.8       5.4%            7.3%
#   1.0       3.3%            2.4%   <- chosen
# The PSF is marginally softer at 1.0, but PyAutoLens fits the PSF explicitly, so the
# uniform low-covariance noise map is the better trade for the strong-lens modelling.
# (Earlier pixfrac tuning preferred 0.8 on stacked FWHM; that predates ERR weighting
# and the shift to prioritising noise-map covariance.)
IR_OUT_SCALE, IR_OUT_PIXFRAC = 0.06, 1.0

# ── DQ bits treated as good ────────────────────────────────────────────────────
# Only 512 (IR flat-field "blob"). Everything else flagged is rejected -- in
# particular 16 (hot) and 32 (unstable), which are the bright defects: measured on
# J0841+3824 their median |residual| is 1.9-3.3x the clean-pixel MAD with a 90th
# percentile of 16-62x, and 25-38% of them exceed 5 MAD. Because the no-CR pass runs
# no cross-frame rejection at all, any bit kept as "good" is drizzled as-is at each
# dither position, so one detector pixel becomes 4 separate replicas in the output
# (the "each spot turned into 4" report). 512 is kept because blob pixels are
# photometrically indistinguishable from clean ones (median |residual| 0.70 MAD vs
# 0.67) while masking them would raise zero-coverage sky pixels from 0.19% to 0.77%.
# Bit 64 (warm) is never actually set in these FLTs, so it is not listed.
#
# NEVER write this as '' or None. astropy's interpret_bit_flags maps both to None,
# which AstroDrizzle logs as `bits : None` and which disables DQ masking *entirely* --
# it silently keeps every flagged pixel, and also voids any CR flag written into DQ
# 4096 by driz_cr or LACosmic. The assert below is the guard against that.
from astropy.nddata.bitmask import interpret_bit_flags as _interpret_bit_flags
_DQ_GOOD = '512'
assert _interpret_bit_flags(_DQ_GOOD) is not None, \
    "_DQ_GOOD must name real bits: '' and None disable DQ masking entirely"

CR_BIT = 4096          # DRIZ_CR; not in _DQ_GOOD, so flagged pixels are excluded

# DQ bits whose flags are advisory rather than measured per-exposure, and which
# refine_dq_flags() is allowed to clear when the pixel behaves normally in the frame
# being masked: 8 (deviant zero-read), 16 (hot), 32 (unstable response).
# NOT included, deliberately: 4 (bad detector pixel -- a permanent defect, and 68% of
# them are deviant anyway), 256 (saturated -- a real ceiling, not a noise statement),
# 2048, and 512 (already kept as good).
_SOFT_DQ = 8 | 16 | 32

_num_cores = 1

# ── Helpers ────────────────────────────────────────────────────────────────────
def refine_dq_flags(files, sigma=3.0, drop_stale_cr=True):
    """Clear _SOFT_DQ flags on pixels that are not deviant in their own exposure.

    Why this exists. WFC3/IR hot/unstable flags are inherited from the dark reference
    file, which characterises a pixel over a whole anneal cycle; a pixel that misbehaved
    once is flagged in every exposure of the cycle. 1.96% of input pixels were flagged,
    and the 4-point dither lands each one at a different sky position, so the noise map
    is speckled with their coverage deficits.

    Note the sky area lost is set by the INPUT pixel size (0.1283"), not the output
    grid. Each masked input pixel renders across (0.1283/0.06)**2 = 4.57 output pixels
    at 0.06", which is why these read as visible blobs rather than the single pixels
    ACS produces -- but that is a rendering property. Coarsening the output grid does
    not recover the sky: measured, 0.06" -> 0.08" cuts the affected area only 18%
    (2.69 -> 2.20 arcsec**2). Do not reach for the output scale to fix this.

    But the flags are far more conservative than the data justify. Fraction of pixels
    whose residual from a 5x5 local median exceeds 3 MAD, in the same exposure that
    flags them (unflagged pixels sit at 1.2% for scale):

        DQ 16 hot 37%    DQ 32 unstable 43%    DQ 8 deviant-zero-read 35%
        DQ 4 bad detector px 68%   DQ 512 blob 1.9% (kept as good)

    So ~60% of what is masked is indistinguishable from clean sky in the frame where it
    is being discarded. Clearing those takes the masked fraction 1.96% -> 0.83% and the
    degraded area from ~26% to ~15%, at no cost in resolution.

    This tests each exposure separately, so a pixel stays masked in the frames where it
    actually misbehaves. It errs conservative on real sources: a flagged pixel sitting
    on a steep PSF has a large local-median residual and keeps its flag.

    `drop_stale_cr` also clears DQ 4096, which is not a calibration flag at all -- it is
    written by driz_cr/LACosmic and persists in the file afterwards, so a rejected
    experiment's flags would otherwise leak into later runs.
    """
    if sigma <= 0:
        print('DQ refinement disabled (--dq-refine 0)')
        return
    print(f'\n=== Refining DQ flags (bits {_SOFT_DQ}, {sigma} sigma) ===')
    for f in files:
        with fits.open(f, mode='update') as h:
            sci = h['SCI', 1].data.astype(np.float64)
            dq = h['DQ', 1].data
            resid = sci - median_filter(sci, 5)
            clean = (dq == 0)
            mad = 1.4826 * np.median(np.abs(resid[clean] - np.median(resid[clean])))
            benign = np.abs(resid) <= sigma * mad

            before = int(((dq & ~_interpret_bit_flags(_DQ_GOOD)) > 0).sum())
            cleared = (dq & _SOFT_DQ) * benign          # bits to drop, per pixel
            dq[:] = dq & ~cleared
            if drop_stale_cr:
                dq[:] = dq & ~CR_BIT
            after = int(((dq & ~_interpret_bit_flags(_DQ_GOOD)) > 0).sum())

            h['DQ', 1].header['DQREFINE'] = (
                sigma, 'sigma for un-flagging non-deviant DQ 8/16/32')
            print(f'  {os.path.basename(f)}: masked {before:,} -> {after:,} px '
                  f'({100 * before / dq.size:.2f}% -> {100 * after / dq.size:.2f}%), '
                  f'MAD={mad:.4f}')


def crop_to_coverage(sci_file, wht_file):
    """Trim sci and wht to the bounding box of wht>0, updating CRPIX in-place."""
    with fits.open(wht_file) as h:
        wht = h[0].data
        rows = np.any(wht > 0, axis=1)
        cols = np.any(wht > 0, axis=0)
        y0, y1 = np.where(rows)[0][[0, -1]]
        x0, x1 = np.where(cols)[0][[0, -1]]
    for fname in (sci_file, wht_file):
        with fits.open(fname, mode='update') as h:
            h[0].data = h[0].data[y0:y1+1, x0:x1+1]
            h[0].header['CRPIX1'] -= x0
            h[0].header['CRPIX2'] -= y0
            h.flush()
    print(f'  Cropped to [{y0}:{y1+1}, {x0}:{x1+1}]  ({x1-x0+1} x {y1-y0+1} px)')

def make_log_norm(data, wht):
    covered = data[wht > 0]
    covered = covered[np.isfinite(covered)]
    vmin = max(np.percentile(covered, 10), 1e-4)
    vmax = np.percentile(covered, 99.9)
    return ImageNormalize(vmin=vmin, vmax=vmax, stretch=LogStretch())

# WFC3/IR LACosmic, mirroring run_lacosmic() in drizzle_acs_wfc.py. Differences:
# one SCI extension (not two), and the FLT is in ELECTRONS/S, so gain=EXPTIME is what
# converts to the electrons astroscrappy's noise model expects.
#
# DO NOT make this the default for WFC3/IR. Measured on J0841+3824: it zeroes 59 of the
# 121 pixels around the field star (the drizzled star goes to WHT=0, SNR=0), because at
# 0.1283"/px the IR PSF is ~1 px FWHM and Laplacian edge detection cannot tell a real
# point source from a cosmic ray. objlim does not rescue it -- at 5/15/20/30/50/100 the
# star is still clipped (23/6/6/6/5/6 px), while the flagged fraction of the unflagged
# >5-sigma outliers it is meant to catch never exceeds 12%. driz_cr is milder but also
# point-source-lossy here (star peak -10%, deflector F(1") -1.7%). The F160W science
# product is therefore the no-CR pass; correct DQ masking (_DQ_GOOD) does the work.
_IR_RDNOISE, _IR_SATLEVEL = 12.0, 78000.0

def run_lacosmic(files, sigclip, objlim):
    """Flag cosmic rays and unflagged outliers into DQ bit 4096 of each FLT, in place."""
    import astroscrappy
    total = 0
    for fname in files:
        with fits.open(fname, mode='update') as hdul:
            exptime = hdul[0].header['EXPTIME']
            sci = hdul['SCI', 1].data.astype(np.float32)
            dq = hdul['DQ', 1].data
            dq &= ~CR_BIT              # clear any previous CR flags
            # genuine defects, so LACosmic does not key off them
            bad = (dq & (256 | 2048)) > 0
            mask, _ = astroscrappy.detect_cosmics(
                sci, inmask=bad, sigclip=sigclip, sigfrac=0.3, objlim=objlim,
                gain=exptime, readnoise=_IR_RDNOISE, satlevel=_IR_SATLEVEL,
                niter=4, sepmed=True, cleantype='medmask', fsmode='median')
            dq[mask] |= CR_BIT
            total += int(mask.sum())
            hdul['DQ', 1].data = dq
    print(f'  LACosmic flagged {total} pixels across {len(files)} frames')

# ── Subprocess mode: run only the no-CR drizzle pass with pre-aligned files ───
# Launched by the main process after the CR pass to get a clean memory slate.
if _is_subprocess:
    os.environ['CRDS_SERVER_URL'] = 'https://hst-crds.stsci.edu'
    os.environ['CRDS_PATH']       = ref_path
    os.environ['iref']            = os.path.join(ref_path, 'references', 'hst', 'wfc3') + os.sep
    os.chdir(work_path)
    flt_files = sorted(glob.glob('*flt.fits'))
    print('\n=== AstroDrizzle (no CR rejection) ===')
    astrodrizzle.AstroDrizzle(flt_files,
                               output='wfc3_ir_flt_nocrrej',
                               preserve=False, build=False, context=False,
                               skysub=True, skymethod='localmin',
                               driz_sep_wcs=True, driz_sep_scale=0.1283,
                               driz_sep_bits=_DQ_GOOD, driz_sep_fillval=-1,
                               median=False, blot=False, driz_cr=False,
                               resetbits=CR_BIT,
                               final_fillval=None, final_bits=_DQ_GOOD,
                               final_wcs=True, final_scale=IR_OUT_SCALE,
                               final_pixfrac=IR_OUT_PIXFRAC,
                               final_wht_type=_a.wht_type,
                               **_common_wcs,
                               num_cores=_num_cores)
    for f in glob.glob('*ask.fits'):
        os.remove(f)
    print('\n=== Cropping (no CR) ===')
    crop_to_coverage('wfc3_ir_flt_nocrrej_drz_sci.fits', 'wfc3_ir_flt_nocrrej_drz_wht.fits')
    sys.exit(0)

# ── Download FLT files from MAST (skip if already present) ────────────────────
with open(json_path) as _f:
    lens_products = json.load(_f)

if glob.glob(os.path.join(data_path, '*flt.fits')):
    print(f'=== MAST download: FLT files already present in {data_path}, skipping ===')
else:
    print(f'=== MAST download: querying {lens} {filt_key} WFC3/IR ===')
    os.makedirs(data_path, exist_ok=True)
    try:
        obs_table = None
        for _pat in mast_target_names.target_patterns(lens):
            obs_table = Observations.query_criteria(
                target_name=_pat,
                obs_collection='HST',
                instrument_name='WFC3/IR',
                filters=[filt_key.upper()],
            )
            if len(obs_table) > 0:
                print(f'  MAST target {_pat}: {len(obs_table)} observations')
                break
            print(f'  MAST target {_pat}: no observations')
        _non_copy = [t for t in obs_table['target_name'] if 'COPY' not in t.upper()]
        if _non_copy:
            obs_table = obs_table[np.array(['COPY' not in t.upper() for t in obs_table['target_name']])]
            print(f'  Using {len(obs_table)} non-COPY observations')
        else:
            print(f'  No non-COPY observations found, using COPY data')
        products = Observations.get_product_list(obs_table)
        Observations.download_products(
            products,
            download_dir=data_path,
            productSubGroupDescription=['FLT'],
            project=['CALWF3'],
        )
        for flt_file in glob.glob(os.path.join(data_path, 'mastDownload', '**', '*flt.fits'),
                                   recursive=True):
            shutil.move(flt_file, data_path)
        shutil.rmtree(os.path.join(data_path, 'mastDownload'), ignore_errors=True)

        print(f'  Downloaded '
              f'{len(glob.glob(os.path.join(data_path, "*flt.fits")))} exposures')
    except Exception as e:
        print(f'  MAST query failed: {e}')

if not glob.glob(os.path.join(data_path, '*flt.fits')):
    _update_info_json(exptime_json_path,    lens, filt_key, None)
    _update_info_json(instrument_json_path, lens, filt_key, None)
    sys.exit(f'No files found for {lens} {filt_key} — check target name and filter')

# Save instrument from first FLT header
with fits.open(sorted(glob.glob(os.path.join(data_path, '*flt.fits')))[0]) as _h:
    _instrume = _h[0].header['INSTRUME'].strip()
    _detector = _h[0].header.get('DETECTOR', 'IR').strip()
_update_info_json(instrument_json_path, lens, filt_key, f'{_instrume}/{_detector}')

# ── Provenance ────────────────────────────────────────────────────────────────
# Record the frames that actually reach the drizzle, not everything the download left
# in data/calibrated/, and refresh on every run -- this used to sit inside the download
# block, so a re-run on already-present files never updated it. EXPTIME=0 frames are
# excluded because AstroDrizzle drops them (that mismatch was real on four ACS
# entries; no WFC3/IR lens currently has one, but the rule is the same).
_obs_ids = sorted(
    os.path.basename(f).replace('_flt.fits', '')
    for f in glob.glob(os.path.join(data_path, '*flt.fits'))
    if fits.getheader(f)['EXPTIME'] > 0
)
lens_products.setdefault(lens, {})[filt_key] = _obs_ids
lens_products[lens] = dict(sorted(lens_products[lens].items()))
with open(json_path, 'w') as _f:
    json.dump(dict(sorted(lens_products.items())), _f, indent=4)

# Skip drizzle if final products already exist
_skip_sentinel = 'wfc3_ir_flt_cr_drz_sci.fits' if do_cr else 'wfc3_ir_flt_nocrrej_drz_sci.fits'
if os.path.exists(os.path.join(output_path, _skip_sentinel)):
    print(f'=== {lens} {filt_key}: drizzled products already exist, skipping ===')
    sys.exit(0)

for p in (output_path, work_path):
    if os.path.exists(p):
        shutil.rmtree(p)
    os.makedirs(p)

# ── Capture stdout/stderr to log file ─────────────────────────────────────────
_log_file = open(os.path.join(work_path, 'run.log'), 'w')
_orig_out = sys.stdout
_orig_err = sys.stderr

class _Tee:
    def __init__(self, primary, secondary):
        self._p, self._s = primary, secondary
    def write(self, data):
        try:
            self._p.write(data)
        except Exception:
            pass
        self._s.write(data)
    def flush(self):
        try:
            self._p.flush()
        except Exception:
            pass
        self._s.flush()
    def __getattr__(self, name):
        return getattr(self._p, name)

sys.stdout = _Tee(_orig_out, _log_file)
sys.stderr = _Tee(_orig_err, _log_file)

# ── CRDS / reference files ─────────────────────────────────────────────────────
os.environ['CRDS_SERVER_URL'] = 'https://hst-crds.stsci.edu'
os.environ['CRDS_PATH']       = ref_path
os.environ['iref']            = os.path.join(
    ref_path, 'references', 'hst', 'wfc3') + os.sep

# ── Copy FLT files to work directory and work there ───────────────────────────
for f in glob.glob(os.path.join(data_path, '*flt.fits')):
    shutil.copy(f, work_path)

os.chdir(work_path)
print(f'Working directory: {work_path}')

# Refine the inherited DQ flags on the *copies*, never on data/calibrated/.
refine_dq_flags(sorted(glob.glob('*flt.fits')), sigma=_a.dq_refine)

# ── Download reference files ──────────────────────────────────────────────────
# Always run bestrefs against the actual input files: it is idempotent (CRDS only
# fetches missing refs) and cheap when they are present. A "skip if the ref dir is
# non-empty" guard was wrong across filters — one filter's run leaves another filter's
# refs unfetched, then updatewcs fails with a missing reference file.
print('\n=== CRDS bestrefs ===')
os.system('crds bestrefs --files *flt.fits --sync-references=1 --update-bestrefs')

# ── Alignment ─────────────────────────────────────────────────────────────────
# DO NOT re-solve the WCS. MAST delivers these files already fitted to GSC 2.4.2 or
# GAIA eDR3, and their *relative* astrometry across a dither sequence is far better
# than anything TweakReg derives here. TweakReg aligns every frame onto the FIRST
# frame, so for dithered exposures it measures the dither itself as an error and
# removes it, leaving AstroDrizzle to stack dithered frames as if they shared a
# pointing. Measured as the rms frame-to-frame scatter of the WCS error (a common
# offset is harmless; scatter is what smears the stack):
#
#   dataset                    as delivered   after updatewcs   after TweakReg
#   WFC3/IR  J0728+3835           0.05 px         0.01 px        0.89 px (max 3.26)
#   WFPC2 1-visit J0029-0055      0.19 px         0.38 px        1.54 px (max 4.69)
#   WFPC2 2-visit J0822+2652      0.77 px         0.90 px        1.55 px (max 3.90)
#   ACS      J0330-0020           0.82 px*        0.82 px*       3.61 px
#     (*ACS measured as spread of |err|, same conclusion)
#
# Even the two-visit WFPC2 case -- different guide stars, 14 deg roll difference --
# is already registered to 0.77 px by the delivered WCS, and TweakReg doubles the
# error. The TweakReg threshold/searchrad tuning kept below was fixing failures in a
# step that should not run at all; it is retained only for --align tweakreg.
flt_files = sorted(glob.glob('*flt.fits'))

if _a.align == 'tweakreg':
    print('\n=== updatewcs ===')
    updatewcs('*flt.fits', use_db=False)
    print('\n=== TweakReg ===')
    tweakreg.TweakReg(flt_files,
                      updatehdr=True,
                      clean=True,
                      reusename=True,
                      interactive=False,
                      conv_width=2.5,
                      # 200 was inherited from the WFPC2/PC tuning and is meaningless
                      # for WFC3/IR, whose FLTs are in ELECTRONS/S: the sky sits at
                      # ~0.7 e/s and the 99.9th percentile at ~8 e/s, so a 200 e/s cut
                      # left only 4-6 objects per image, below minobj=7, and TweakReg
                      # aligned nothing (no shiftfile written). Measured on J0029-0055:
                      # every threshold from 50 down to 1.5 yields the same 4/4 match and
                      # the same 4.61 px solution, so the fit is insensitive here; 20
                      # gives ~40-50 objects/image, clear of both the minobj floor and
                      # the noise.
                      threshold=20.0,
                      ylimit=0.2,
                      shiftfile=True,
                      outshifts='shift_flt.txt',
                      searchrad=1,
                      tolerance=3,
                      minobj=7)

    with open('shift_flt.txt') as f:
        for i, line in enumerate(f, 1):
            if 'nan' in line:
                raise ValueError(f'nan in shift_flt.txt line {i} — TweakReg alignment failed')
else:
    _names = set()
    for _f in sorted(glob.glob('*flt.fits')):
        try:
            _names.add(fits.getval(_f, 'WCSNAME', ext=('SCI', 1)))
        except Exception:
            _names.add('<none>')
    print(f'\n=== Alignment: using MAST WCS as delivered ({", ".join(sorted(_names))}) ===')
    if not all(('GSC' in n or 'GAIA' in n or 'HSC' in n) for n in _names):
        print('  WARNING: a frame lacks a fitted (GSC/GAIA/HSC) WCS. Relative alignment is unverified —')
        print('           check the stacked PSF, or re-run with --align tweakreg.')

# ── AstroDrizzle pass 1: with CR rejection ────────────────────────────────────
if do_cr:
    # Masking DQ alone cannot clean this data: ~8400 px/frame carry DQ==0 yet exceed
    # 5 MAD (1052 clumps on J0841+3824), and with no cross-frame rejection each one
    # drizzles into 4 dither replicas. Only a rejection step reaches them.
    #
    # LACosmic (default) flags per frame into DQ 4096 before drizzling, then combines
    # a plain weighted mean -- same route as ACS and WFPC2, where driz_cr was found to
    # eat the deflector core. resetbits=0 is mandatory: AstroDrizzle defaults it to
    # 4096 and would clear the very bit the mask was written into.
    if _a.cr_method == 'lacosmic':
        print('\n=== LACosmic CR masking ===')
        run_lacosmic(flt_files, _a.lacosmic_sigclip, _a.lacosmic_objlim)
        _cr_kw = dict(median=False, blot=False, driz_cr=False, resetbits=0)
    else:
        _cr_kw = dict(median=True, blot=True, driz_cr=True,
                      driz_cr_snr='3.5 3.0', driz_cr_scale='1.2 0.7')
    print(f'\n=== AstroDrizzle (with CR rejection, {_a.cr_method}) ===')
    astrodrizzle.AstroDrizzle(flt_files,
                               output='wfc3_ir_flt_cr',
                               preserve=False, build=False, context=False,
                               skysub=True, skymethod='localmin',
                               driz_sep_wcs=True, driz_sep_scale=0.1283,
                               driz_sep_bits=_DQ_GOOD, driz_sep_fillval=-1,
                               **_cr_kw,
                               final_fillval=None, final_bits=_DQ_GOOD,
                               final_wcs=True, final_scale=IR_OUT_SCALE,
                               final_pixfrac=IR_OUT_PIXFRAC,
                               final_wht_type=_a.wht_type,
                               **_common_wcs,
                               num_cores=_num_cores)
    for f in glob.glob('*ask.fits'):
        os.remove(f)

    print('\n=== Cropping (CR) ===')
    crop_to_coverage('wfc3_ir_flt_cr_drz_sci.fits', 'wfc3_ir_flt_cr_drz_wht.fits')

    # ── AstroDrizzle pass 2: no CR rejection in subprocess to free memory ─────
    print('\n=== AstroDrizzle (no CR rejection) — launching subprocess ===')
    _result = subprocess.run(
        [sys.executable, os.path.abspath(__file__),
         '--lens', lens, '--filt', filt, '--sample', sample,
         '--wht-type', _a.wht_type, '--_subprocess'],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
    )
    sys.stdout.write(_result.stdout)
    sys.stdout.flush()
    if _result.returncode != 0:
        raise RuntimeError(f'no-CR subprocess exited with code {_result.returncode}')

else:
    # ── AstroDrizzle pass 2: no CR rejection ──────────────────────────────────
    print('\n=== AstroDrizzle (no CR rejection) ===')
    astrodrizzle.AstroDrizzle(flt_files,
                               output='wfc3_ir_flt_nocrrej',
                               preserve=False, build=False, context=False,
                               skysub=True, skymethod='localmin',
                               driz_sep_wcs=True, driz_sep_scale=0.1283,
                               driz_sep_bits=_DQ_GOOD, driz_sep_fillval=-1,
                               median=False, blot=False, driz_cr=False,
                               resetbits=CR_BIT,
                               final_fillval=None, final_bits=_DQ_GOOD,
                               final_wcs=True, final_scale=IR_OUT_SCALE,
                               final_pixfrac=IR_OUT_PIXFRAC,
                               final_wht_type=_a.wht_type,
                               **_common_wcs,
                               num_cores=_num_cores)
    for f in glob.glob('*ask.fits'):
        os.remove(f)

    print('\n=== Cropping (no CR) ===')
    crop_to_coverage('wfc3_ir_flt_nocrrej_drz_sci.fits', 'wfc3_ir_flt_nocrrej_drz_wht.fits')

print('\n=== Exposure times ===')
exptime = fits.getheader('wfc3_ir_flt_nocrrej_drz_sci.fits')['EXPTIME']
print(f'  No CR rejection: {exptime:.1f} s')
if do_cr:
    print(f'  CR rejected: {fits.getheader("wfc3_ir_flt_cr_drz_sci.fits")["EXPTIME"]:.1f} s')
_update_info_json(exptime_json_path, lens, filt_key, exptime)

# ── Copy final sci/wht to output_path ────────────────────────────────────────
print('\n=== Copying final products to output directory ===')
_copy = ['wfc3_ir_flt_nocrrej_drz_sci.fits', 'wfc3_ir_flt_nocrrej_drz_wht.fits']
if do_cr:
    _copy = ['wfc3_ir_flt_cr_drz_sci.fits', 'wfc3_ir_flt_cr_drz_wht.fits'] + _copy
for fname in _copy:
    shutil.copy(fname, output_path)
    print(f'  {fname}')

# ── Plots ──────────────────────────────────────────────────────────────────────
print('\n=== Saving plots ===')

single_sci_files = sorted(glob.glob('*_single_sci.fits'))
single_wht_files = sorted(glob.glob('*_single_wht.fits'))
if single_sci_files:
    n = len(single_sci_files)
    fig, axes = plt.subplots(1, n, figsize=(8 * n, 8))
    if n == 1:
        axes = [axes]
    for ax, sci_f, wht_f in zip(axes, single_sci_files, single_wht_files):
        data     = fits.getdata(sci_f)
        wht_data = fits.getdata(wht_f)
        ax.imshow(data, norm=make_log_norm(data, wht_data), cmap='gray', origin='lower')
        ax.set_title(os.path.basename(sci_f).replace('_single_sci.fits', ''))
    fig.tight_layout()
    fig.savefig(os.path.join(work_path, 'single_sci.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print('  single_sci.png')

def save_drizzled_png(sci_file, wht_file, out_png, title=''):
    sci = fits.getdata(sci_file)
    wht = fits.getdata(wht_file)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
    ax1.imshow(sci, norm=make_log_norm(sci, wht), cmap='gray', origin='lower')
    ax1.set_title(f'SCI  {title}')
    ax2.imshow(wht, cmap='gray', vmin=0, vmax=wht.max(), origin='lower')
    ax2.set_title('WHT')
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches='tight')
    plt.close(fig)

if do_cr:
    save_drizzled_png('wfc3_ir_flt_cr_drz_sci.fits', 'wfc3_ir_flt_cr_drz_wht.fits',
                      os.path.join(output_path, 'drizzled_cr.png'), title='(CR rejection)')
    print('  drizzled_cr.png')

save_drizzled_png('wfc3_ir_flt_nocrrej_drz_sci.fits', 'wfc3_ir_flt_nocrrej_drz_wht.fits',
                  os.path.join(output_path, 'drizzled_nocrrej.png'), title='(no CR rejection)')
print('  drizzled_nocrrej.png')

if do_cr:
    sci_cr      = fits.getdata('wfc3_ir_flt_cr_drz_sci.fits')
    sci_nocrrej = fits.getdata('wfc3_ir_flt_nocrrej_drz_sci.fits')
    wht         = fits.getdata('wfc3_ir_flt_cr_drz_wht.fits')
    residual    = sci_cr - sci_nocrrej

    combined = np.concatenate([sci_nocrrej[wht > 0], sci_cr[wht > 0]])
    combined = combined[np.isfinite(combined)]
    vmin      = max(np.percentile(combined, 10), 1e-4)
    vmax_sci  = np.percentile(combined, 99.9)
    shared_norm = ImageNormalize(vmin=vmin, vmax=vmax_sci, stretch=LogStretch())

    covered_res = residual[wht > 0]
    covered_res = covered_res[np.isfinite(covered_res)]
    vmax_res    = np.percentile(np.abs(covered_res), 99.5)

    fig, axes = plt.subplots(1, 3, figsize=(24, 8))
    axes[0].imshow(sci_nocrrej, norm=shared_norm, cmap='gray', origin='lower')
    axes[0].set_title('No CR rejection')
    axes[1].imshow(sci_cr, norm=shared_norm, cmap='gray', origin='lower')
    axes[1].set_title('CR rejection (snr=3.5/3.0, scale=1.2/0.7)')
    im = axes[2].imshow(residual, cmap='RdBu_r', vmin=-vmax_res, vmax=vmax_res, origin='lower')
    axes[2].set_title('CR rejected minus no CR rejection')
    plt.colorbar(im, ax=axes[2], label='counts/s')
    fig.tight_layout()
    fig.savefig(os.path.join(work_path, 'drizzled_diff.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print('  drizzled_diff.png')

print(f'\nDone. Output in: {output_path}')

sys.stdout = _orig_out
sys.stderr = _orig_err
_log_file.close()
