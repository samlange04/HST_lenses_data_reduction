#!/usr/bin/env python
"""
Drizzle HST WFC3/UVIS FLC images into a combined mosaic (BELLS GALLERY sample).
Downloads FLC files from MAST into data/calibrated/ if not already present.
Writes final sci/wht FITS to data/drizzled/, intermediates to data/drizzle_files/.
Runs the LACosmic CR-rejected pass by default (the science product); the no-CR
comparison pass is opt-in via --nocrrej.

WFC3/UVIS is a two-CCD optical detector (a pair of 4096x2051 chips, ~0.0396"/px,
full well ~70,000 e-; WFC3 Data Handbook 2024 s1.2), so it is the direct analog of
ACS/WFC and this script is modeled closely on drizzle_acs_wfc.py: same two-chip
LACosmic masking, same "trust the delivered GSC242/GAIAeDR3 WCS" alignment default,
same ERR-weighting. FLC (CTE-corrected, CALWF3) is the recommended input, as for ACS.

The gallery is WFC3/UVIS multi-band (F225W, F275W, F438W, F606W, F814W); this script
is filter-agnostic (--filt), F606W is the primary band. The alignment and CR reasoning
is inherited from the ACS sibling (same optical-CCD dither/CR physics) and has NOT yet
been independently audited on gallery data -- see the comments below and re-measure
(stacked stellar FWHM, core-registration scatter) before trusting a batch campaign.
"""

import argparse
import json
import os
import glob
import shutil
import signal
import subprocess
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
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
import info_json

# WFC3/UVIS native plate scale, UVIS1 x-axis (WFC3 DHB 2024 s4.1.1: 0.0396"/px x,
# 0.0393 y on chip 1; 0.0398 both on chip 2). Drizzle to the native scale -- like
# ACS at 0.05", one masked input pixel maps to ~one output pixel, and PyAutoLens
# ingests each band at its native scale so nothing is gained by oversampling here.
FINAL_SCALE = 0.0396

# Drizzle drop size. Chosen by a pixfrac scan on J1110+3649 F606W (12 exposures) measuring
# drizzle correlated noise against weight-map uniformity on a fixed blank-sky box:
#     pixfrac  lag-1 corr  integrated inflation(8x8)  wht-CoV
#       1.0       0.29            1.81                 0.073
#       0.8       0.22            1.65                 0.078
#       0.7       0.18            1.58                 0.081
#       0.6       0.15            1.52                 0.084
# Going 1.0->0.7 cuts the neighbour correlation markedly (the noise map stops looking
# "blurry") while weight uniformity barely moves (median weight flat -> no coverage holes
# on this densely-dithered stack). 0.7 keeps margin for sparser-dither gallery lenses;
# below it the gains shrink and hole risk grows. Native-scale regime, so this is the
# opposite lever from the oversampled F606W/F160W bands (which chose pixfrac 1.0).
# → memory: drizzle_correlated_noise. User-settled 2026-07-28.
FINAL_PIXFRAC = 0.7

# DQ bits treated as GOOD so they are not drizzled out (and not mistaken for CRs by
# LACosmic). WFC3/UVIS DQ definitions (WFC3 DHB 2024 Table 2.4, p.38):
#   16  stable hot pixel      64  warm pixel (obsolete flag)   256  full-well saturation
# This mirrors the ACS/WFC choice '256,64,16' (the bit meanings coincide). The STScI
# UVIS drizzle examples (uvis_time_dependent_photometry notebook) use '80' (=16+64,
# i.e. drop saturation from "good"); revisit if saturated cores prove a problem.
_DQ_GOOD = '256,64,16'

# ── Configuration ──────────────────────────────────────────────────────────────
_p = argparse.ArgumentParser()
_p.add_argument('--lens',        default='J1110+3649')
_p.add_argument('--filt',        default='f606W')
_p.add_argument('--sample',      default='gallery',
                help='sample the lens belongs to; sets the <sample> level of every '
                     'data/ path. Defined in info/lens_samples.json (default gallery, '
                     'the WFC3/UVIS BELLS GALLERY sample)')
# CR rejection is the science default: the LACosmic CR pass is the product the rest of
# the pipeline consumes. Pass --no-cr to skip it. The no-CR ("nocrrej") pass is opt-in
# for comparison only via --nocrrej. (Same convention as drizzle_acs_wfc.py.)
_p.add_argument('--cr',      action=argparse.BooleanOptionalAction, default=True)
_p.add_argument('--nocrrej', action=argparse.BooleanOptionalAction, default=False,
                help='also produce the no-CR comparison drizzle (off by default)')
_p.add_argument('--cr-method',   default='lacosmic', choices=['lacosmic', 'drizcr'])
# Alignment source. 'mast' trusts the GSC242/GAIAeDR3-fitted WCS in the delivered FLCs
# and runs neither updatewcs nor TweakReg; 'tweakreg' is the old re-solve, which erases
# the dither on these data. See the block comment above the alignment step.
#
# No per-lens gallery override has been measured yet (unlike ACS's J1213+6708). An
# explicit --align on the command line always wins over this table.
ALIGN_OVERRIDES = {}
_p.add_argument('--align',       default=None, choices=['mast', 'tweakreg'])
# Total-exposure-time gate (frames that reach the drizzle, i.e. EXPTIME>0 only).
# Below BLOCK_EXPTIME no product is written (same outcome as no MAST data); between
# BLOCK and WARN the drizzle proceeds but is flagged.
WARN_EXPTIME, BLOCK_EXPTIME = 1200.0, 500.0
# Drizzle output weight type. 'ERR' (default) makes the WHT extension a full
# inverse-variance map (source Poisson + sky + read + dark), so 1/sqrt(WHT) is a
# CALIBRATED per-pixel noise map -- what make_cutouts uses. The UVIS FLC ERR array is
# in electrons (like ACS FLC), so 1/sqrt(WHT) is a sigma map with K=1 (no per-frame
# EXPTIME rescale, unlike WFC3/IR whose ERR is in electrons/s). See noise_map_err_units.
_p.add_argument('--wht-type',    default='ERR', choices=['ERR', 'IVM', 'EXP'])
_p.add_argument('--lacosmic-sigclip', type=float, default=4.5)
_p.add_argument('--lacosmic-objlim',  type=float, default=5.0)
# driz_cr tuning, only used by --cr-method drizcr (see run_lacosmic note for why
# lacosmic is the default). Loosened defaults for the sharp UVIS core, as for ACS.
_p.add_argument('--driz-cr-snr',   default='3.5 3.0')
_p.add_argument('--driz-cr-scale', default='1.2 0.7')
_p.add_argument('--combine-type',  default=None)
_p.add_argument('--_subprocess', action='store_true', default=False, help=argparse.SUPPRESS)
_a = _p.parse_args()

if _a.align is None:
    _a.align = ALIGN_OVERRIDES.get((_a.lens, _a.filt), 'mast')
    print(f'=== align: {_a.align} (default for this lens/filter) ===')

lens           = _a.lens
sample         = _a.sample
filt           = _a.filt
do_cr          = _a.cr
do_nocrrej     = _a.nocrrej
_is_subprocess = _a._subprocess
if not do_cr and not do_nocrrej:
    _p.error('nothing to do: --no-cr given without --nocrrej (no drizzle pass requested)')

ws_path     = '/Users/samlange/Code/HST_lenses_data_reduction'
data_path   = os.path.join(ws_path, 'data', 'calibrated', sample, lens, filt)
output_path = os.path.join(ws_path, 'data', 'drizzled', sample, lens, filt)
# DRIZZLE_WORK_ROOT lets us relocate the AstroDrizzle working dir (where the large
# output FITS is written) onto a different volume/backing store -- used to isolate the
# APFS-write U-state hang (e.g. point it at a RAM disk). Falls back to the in-repo path.
_work_root  = os.environ.get('DRIZZLE_WORK_ROOT') or os.path.join(ws_path, 'data', 'drizzle_files')
work_path   = os.path.join(_work_root, sample, lens, filt)
ref_path    = os.path.join(ws_path, 'data', 'reference_files')

# ── Common output WCS across filters (orientation + centre) ────────────────────
# Pin every band to North-up (final_rot=0) with the tangent point (final_ra/dec) at the
# lens position, so all five gallery filters share orientation and WCS centre and the
# lens sits ~centred. gallery lenses live in info/gallery_coords.py (SLACS lenses live
# in slacs_coords.py -- this script is the WFC3/UVIS gallery driver, so it reads the
# gallery table). Unknown lens -> native drizzle WCS (empty dict), with a warning.
sys.path.insert(0, os.path.join(ws_path, 'info'))
try:
    from gallery_coords import gallery_coords as _coords
    from astropy.coordinates import SkyCoord as _SkyCoord
    import astropy.units as _u
    _lc = _SkyCoord(*_coords[lens], unit=(_u.hourangle, _u.deg))
    _common_wcs = dict(final_rot=0.0, final_ra=float(_lc.ra.deg), final_dec=float(_lc.dec.deg))
    print(f'=== Common WCS: North-up, centre ({_lc.ra.deg:.5f}, {_lc.dec.deg:.5f}) ===')
except (KeyError, ImportError):
    _common_wcs = {}
    _lc = None
    print('=== Common WCS: lens not in gallery_coords -> native drizzle WCS ===')


# ── Info JSON paths ────────────────────────────────────────────────────────────
filt_key             = filt
json_path            = os.path.join(ws_path, 'info', 'lens_products.json')
exptime_json_path    = os.path.join(ws_path, 'info', 'lens_exptime.json')
instrument_json_path = os.path.join(ws_path, 'info', 'lens_instrument.json')

_num_cores = 1

# ── Helpers ────────────────────────────────────────────────────────────────────
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

# ── Subprocess mode: run only the no-CR drizzle pass with pre-aligned files ───
# Launched by the main process after the CR pass to get a clean memory slate.
if _is_subprocess:
    os.environ['CRDS_SERVER_URL'] = 'https://hst-crds.stsci.edu'
    os.environ['CRDS_PATH']       = ref_path
    os.environ['iref']            = os.path.join(ref_path, 'references', 'hst', 'wfc3') + os.sep
    os.chdir(work_path)
    flc_files = sorted(glob.glob('*flc.fits'))
    print('\n=== AstroDrizzle (no CR rejection) ===')
    astrodrizzle.AstroDrizzle(flc_files,
                               output='wfc3_uvis_flc_nocrrej',
                               preserve=False, build=False, context=False,
                               skysub=True, skymethod='localmin',
                               driz_sep_wcs=True, driz_sep_scale=FINAL_SCALE,
                               driz_sep_bits=_DQ_GOOD, driz_sep_fillval=-1,
                               median=False, blot=False, driz_cr=False,
                               # Must clear bit 4096: when --cr ran first it left the
                               # LACosmic mask in the input DQ, and inheriting it here
                               # would make the "no CR rejection" product silently
                               # CR-rejected. AstroDrizzle's default, pinned explicitly.
                               resetbits=4096,
                               final_fillval=None, final_bits=_DQ_GOOD,
                               final_wcs=True, final_scale=FINAL_SCALE, final_pixfrac=FINAL_PIXFRAC,
                               final_wht_type=_a.wht_type,
                               **_common_wcs,
                               num_cores=_num_cores)
    for f in glob.glob('*ask.fits'):
        os.remove(f)
    print('\n=== Cropping (no CR) ===')
    crop_to_coverage('wfc3_uvis_flc_nocrrej_drc_sci.fits', 'wfc3_uvis_flc_nocrrej_drc_wht.fits')
    sys.exit(0)

# ── Download FLC files from MAST (skip if already present) ────────────────────
# Distinguishes "MAST has nothing for this lens+filter" (an ordinary outcome -- every
# lens in a sample is tried on every run, most have no data in most bands) from "the
# download broke" (a real failure). Conflating them is what would let a network error
# be silently recorded as `null` = no data. See the no-data exit below.
_mast_empty = False
_mast_error = None

if glob.glob(os.path.join(data_path, '*flc.fits')):
    print(f'=== MAST download: FLC files already present in {data_path}, skipping ===')
else:
    print(f'=== MAST download: querying {lens} {filt_key} WFC3/UVIS ===')
    os.makedirs(data_path, exist_ok=True)
    try:
        obs_table = None
        for _pat in mast_target_names.target_patterns(lens, sample):
            obs_table = Observations.query_criteria(
                target_name=_pat,
                obs_collection='HST',
                instrument_name='WFC3/UVIS',
                filters=[filt_key.upper()],
            )
            if len(obs_table) > 0:
                print(f'  MAST target {_pat}: {len(obs_table)} observations')
                break
            print(f'  MAST target {_pat}: no observations')
        if obs_table is None or len(obs_table) == 0:
            # Not an error: this lens simply has no WFC3/UVIS data in this filter.
            _mast_empty = True
            raise mast_target_names.NoMastData()
        _copy_mask     = np.array(['COPY' in t.upper() for t in obs_table['target_name']])
        _non_copy      = [t for t in obs_table['target_name'] if 'COPY' not in t.upper()]
        if mast_target_names.force_copy(lens) and _copy_mask.any():
            obs_table = obs_table[_copy_mask]
            print(f'  Forcing {len(obs_table)} COPY observations for {lens} (non-COPY unusable)')
        elif _non_copy:
            obs_table = obs_table[~_copy_mask]
            print(f'  Using {len(obs_table)} non-COPY observations')
        else:
            print(f'  No non-COPY observations found, using COPY data')
        products = Observations.get_product_list(obs_table)
        Observations.download_products(
            products,
            download_dir=data_path,
            productSubGroupDescription=['FLC'],
            project=['CALWF3'],
        )
        for flc_file in glob.glob(os.path.join(data_path, 'mastDownload', '**', '*flc.fits'),
                                   recursive=True):
            shutil.move(flc_file, data_path)
        shutil.rmtree(os.path.join(data_path, 'mastDownload'), ignore_errors=True)

        print(f'  Downloaded '
              f'{len(glob.glob(os.path.join(data_path, "*flc.fits")))} exposures')
    except mast_target_names.NoMastData:
        print(f'  No WFC3/UVIS {filt_key.upper()} observations for {lens} on MAST')
    except Exception as e:
        _mast_error = e
        print(f'  MAST query failed: {e}')

if not glob.glob(os.path.join(data_path, '*flc.fits')):
    info_json.update(exptime_json_path,    sample, lens, filt_key, None)
    info_json.update(instrument_json_path, sample, lens, filt_key, None)
    if _mast_empty:
        # Ordinary outcome, not a failure: exit 0 so a batch runner sweeping the whole
        # sample records "no data" and moves on instead of counting it as an error.
        print(f'=== NO DATA: {lens} has no WFC3/UVIS {filt_key.upper()} on MAST '
              f'(recorded as null) ===')
        sys.exit(0)
    if _mast_error is not None:
        sys.exit(f'MAST download failed for {lens} {filt_key}: {_mast_error}')
    sys.exit(f'No files found for {lens} {filt_key} — MAST listed observations but no '
             f'FLC files landed in {data_path}')

# Save instrument from first FLC header
with fits.open(sorted(glob.glob(os.path.join(data_path, '*flc.fits')))[0]) as _h:
    _instrume  = _h[0].header['INSTRUME'].strip()
    _detector  = _h[0].header.get('DETECTOR', 'UVIS').strip()
info_json.update(instrument_json_path, sample, lens, filt_key, f'{_instrume}/{_detector}')

# ── Provenance ────────────────────────────────────────────────────────────────
# Record the frames that actually reach the drizzle, not everything the download left
# in data/calibrated/: AstroDrizzle silently drops EXPTIME=0 frames. EXPTIME=0 frames
# are only filtered out of the record here, not deleted. Note when auditing: WFC3/UVIS
# FLCs are 2-chip MEFs, so the product's NDRIZIM is 2 x the number of exposures listed.
_frames = [
    (os.path.basename(f).replace('_flc.fits', ''), fits.getheader(f)['EXPTIME'])
    for f in glob.glob(os.path.join(data_path, '*flc.fits'))
]
_obs_ids = sorted(rootname for rootname, exp in _frames if exp > 0)
_total_exptime = sum(exp for _, exp in _frames if exp > 0)

if _total_exptime < BLOCK_EXPTIME:
    info_json.update(exptime_json_path,    sample, lens, filt_key, None)
    info_json.update(instrument_json_path, sample, lens, filt_key, None)
    print(f'=== BLOCKED (exptime): {lens} {filt_key} total exptime {_total_exptime:.1f}s '
          f'< {BLOCK_EXPTIME:.0f}s minimum (recorded as null) ===')
    sys.exit(0)
if _total_exptime < WARN_EXPTIME:
    print(f'  EXPTIME WARNING: {lens} {filt_key} total exptime {_total_exptime:.1f}s '
          f'< {WARN_EXPTIME:.0f}s (proceeding)')

info_json.update(json_path, sample, lens, filt_key, _obs_ids)

# Skip drizzle if final products already exist
_skip_sentinel = 'wfc3_uvis_flc_cr_drc_sci.fits' if do_cr else 'wfc3_uvis_flc_nocrrej_drc_sci.fits'
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

def _sig_handler(signum, frame):
    _log_file.write(f'\n[SIGNAL] Received signal {signum} ({signal.Signals(signum).name})\n')
    _log_file.flush()
    signal.signal(signum, signal.SIG_DFL)
    os.kill(os.getpid(), signum)

for _s in (signal.SIGHUP, signal.SIGTERM, signal.SIGINT, signal.SIGPIPE):
    signal.signal(_s, _sig_handler)

# ── CRDS / reference files ─────────────────────────────────────────────────────
os.environ['CRDS_SERVER_URL'] = 'https://hst-crds.stsci.edu'
os.environ['CRDS_PATH']       = ref_path
os.environ['iref']            = os.path.join(
    ref_path, 'references', 'hst', 'wfc3') + os.sep

# ── Copy FLC files to work directory and work there ───────────────────────────
for f in glob.glob(os.path.join(data_path, '*flc.fits')):
    shutil.copy(f, work_path)

os.chdir(work_path)
print(f'Working directory: {work_path}')

# ── Download reference files ──────────────────────────────────────────────────
# Always run bestrefs against the actual input files: it is idempotent (CRDS only
# fetches refs that are missing from the local cache) and cheap when already present.
print('\n=== CRDS bestrefs ===')
os.system('crds bestrefs --files *flc.fits --sync-references=1 --update-bestrefs')

# ── Alignment ─────────────────────────────────────────────────────────────────
# DO NOT re-solve the WCS by default. MAST delivers WFC3/UVIS FLCs already fitted to
# GSC 2.4.2 / GAIA eDR3 (WCSNAME *-FIT_REL_GSC242 / *-GAIAeDR3), whose *relative*
# astrometry across a dither sequence is good to ~0.5 px. TweakReg aligns every frame
# onto the first, so on dithered exposures it measures the dither as an error and
# removes it -- smearing point sources and splitting lensed arcs into offset copies.
# This was measured extensively on ACS/WFC (see drizzle_acs_wfc.py) and the physics is
# identical for UVIS. It has NOT been re-measured on gallery data yet -- if a stacked
# PSF looks off, compare --align mast vs tweakreg (stacked stellar FWHM) before trusting
# a batch run. --align tweakreg restores the re-solve for that comparison.
flc_files = sorted(glob.glob('*flc.fits'))

if _a.align == 'tweakreg':
    print('\n=== updatewcs ===')
    updatewcs('*flc.fits', use_db=False)
    print('\n=== TweakReg ===')
    tweakreg.TweakReg(flc_files,
                      updatehdr=True,
                      clean=True,
                      reusename=True,
                      interactive=False,
                      conv_width=3.5,
                      threshold=200.0,
                      ylimit=0.2,
                      shiftfile=True,
                      outshifts='shift_flc.txt',
                      searchrad=1,
                      tolerance=3,
                      minobj=7)
    with open('shift_flc.txt') as f:
        for i, line in enumerate(f, 1):
            if 'nan' in line:
                raise ValueError(f'nan in shift_flc.txt line {i} — TweakReg alignment failed')
else:
    _wcsnames = {fits.getval(f, 'WCSNAME', ext=('SCI', 1)) for f in flc_files}
    print(f'\n=== Alignment: using MAST WCS as delivered ({", ".join(sorted(_wcsnames))}) ===')
    if not all(('GSC' in n or 'GAIA' in n or 'HSC' in n) for n in _wcsnames):
        print('  WARNING: a frame carries no fitted (GSC/GAIA/HSC) WCS. Relative '
              'alignment is unverified — check the stacked PSF, or re-run with '
              '--align tweakreg.')

# ── CR masking: LACosmic, per frame ───────────────────────────────────────────
# AstroDrizzle's driz_cr detects cosmic rays by comparing each frame against a blotted
# median of the stack; on a steep PSF core that reference is systematically low, so the
# core's residual reads as a cosmic ray (measured on ACS: it destroyed ~37% of the
# deflector flux). LACosmic works one frame at a time with an explicit object-protection
# term (objlim), so there is no stacked reference to bias, and the final drizzle is then
# a plain weighted mean (median/blot/driz_cr all off). Same route and defaults as ACS;
# UVIS shares the optical-CCD/CR physics. --cr-method drizcr restores the old route.
CR_BIT = 4096          # DRIZ_CR; excluded because it is not in _DQ_GOOD
# UVIS FLC is in electrons (CALWF3 applies the gain), so gain=1 for astroscrappy.
# readnoise ~3.1 e- (per-amp READNSEA-D, WFC3 DHB 2024); full well ~70,000 e- (s1.2).
_UVIS_RDNOISE, _UVIS_SATLEVEL = 3.1, 70000.0


def run_lacosmic(files, sigclip, objlim):
    """Flag cosmic rays into DQ bit 4096 of each FLC, in place."""
    import astroscrappy
    total = 0
    for fname in files:
        with fits.open(fname, mode='update') as hdul:
            for ext in range(1, 3):        # WFC3/UVIS has two SCI chips
                sci = hdul['SCI', ext].data.astype(np.float32)
                dq = hdul['DQ', ext].data
                dq &= ~CR_BIT              # clear any previous CR flags
                # genuine defects, so LACosmic does not key off them (UVIS DQ:
                # 4 bad pixel, 8 deviant bias, 128 bad bias ref, 512 bad flat/blob)
                bad = (dq & (4 | 8 | 128 | 512)) > 0
                mask, _ = astroscrappy.detect_cosmics(
                    sci, inmask=bad, sigclip=sigclip, sigfrac=0.3, objlim=objlim,
                    gain=1.0, readnoise=_UVIS_RDNOISE, satlevel=_UVIS_SATLEVEL,
                    niter=4, sepmed=True, cleantype='medmask', fsmode='median')
                dq[mask] |= CR_BIT
                total += int(mask.sum())
                hdul['DQ', ext].data = dq
    print(f'  LACosmic flagged {total} pixels across {len(files)} frames')


# ── AstroDrizzle pass 1: with CR rejection ────────────────────────────────────
if do_cr:
    if _a.cr_method == 'lacosmic':
        print('\n=== LACosmic CR masking ===')
        run_lacosmic(flc_files, _a.lacosmic_sigclip, _a.lacosmic_objlim)
        print('\n=== AstroDrizzle (LACosmic-masked, plain weighted mean) ===')
        astrodrizzle.AstroDrizzle(flc_files,
                                   output='wfc3_uvis_flc_cr',
                                   preserve=False, build=False, context=False,
                                   skysub=True, skymethod='localmin',
                                   driz_sep_wcs=True, driz_sep_scale=FINAL_SCALE,
                                   driz_sep_bits=_DQ_GOOD, driz_sep_fillval=-1,
                                   median=False, blot=False, driz_cr=False,
                                   # resetbits defaults to 4096 and would clear the very
                                   # bit LACosmic just wrote, silently reverting to an
                                   # unmasked drizzle that still looks plausible.
                                   resetbits=0,
                                   final_fillval=None, final_bits=_DQ_GOOD,
                                   final_wcs=True, final_scale=FINAL_SCALE, final_pixfrac=FINAL_PIXFRAC,
                                   final_wht_type=_a.wht_type,
                                   **_common_wcs,
                                   num_cores=_num_cores)
    else:
        print('\n=== AstroDrizzle (with driz_cr CR rejection) ===')
        _combine_kw = {} if _a.combine_type is None else {'combine_type': _a.combine_type}
        astrodrizzle.AstroDrizzle(flc_files,
                                   **_combine_kw,
                                   output='wfc3_uvis_flc_cr',
                                   preserve=False, build=False, context=False,
                                   skysub=True, skymethod='localmin',
                                   driz_sep_wcs=True, driz_sep_scale=FINAL_SCALE,
                                   driz_sep_bits=_DQ_GOOD, driz_sep_fillval=-1,
                                   median=True, blot=True, driz_cr=True,
                                   driz_cr_snr=_a.driz_cr_snr, driz_cr_scale=_a.driz_cr_scale,
                                   final_fillval=None, final_bits=_DQ_GOOD,
                                   final_wcs=True, final_scale=FINAL_SCALE, final_pixfrac=FINAL_PIXFRAC,
                                   final_wht_type=_a.wht_type,
                                   **_common_wcs,
                                   num_cores=_num_cores)
    for f in glob.glob('*ask.fits'):
        os.remove(f)

    print('\n=== Cropping (CR) ===')
    crop_to_coverage('wfc3_uvis_flc_cr_drc_sci.fits', 'wfc3_uvis_flc_cr_drc_wht.fits')

    # ── AstroDrizzle pass 2: no CR rejection in subprocess to free memory ─────
    if do_nocrrej:
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
    astrodrizzle.AstroDrizzle(flc_files,
                               output='wfc3_uvis_flc_nocrrej',
                               preserve=False, build=False, context=False,
                               skysub=True, skymethod='localmin',
                               driz_sep_wcs=True, driz_sep_scale=FINAL_SCALE,
                               driz_sep_bits=_DQ_GOOD, driz_sep_fillval=-1,
                               median=False, blot=False, driz_cr=False,
                               resetbits=4096,
                               final_fillval=None, final_bits=_DQ_GOOD,
                               final_wcs=True, final_scale=FINAL_SCALE, final_pixfrac=FINAL_PIXFRAC,
                               final_wht_type=_a.wht_type,
                               **_common_wcs,
                               num_cores=_num_cores)
    for f in glob.glob('*ask.fits'):
        os.remove(f)

    print('\n=== Cropping (no CR) ===')
    crop_to_coverage('wfc3_uvis_flc_nocrrej_drc_sci.fits', 'wfc3_uvis_flc_nocrrej_drc_wht.fits')

print('\n=== Exposure times ===')
# The CR pass is the primary product; record its EXPTIME (identical between passes,
# but read from whichever product actually exists).
_primary_sci = 'wfc3_uvis_flc_cr_drc_sci.fits' if do_cr else 'wfc3_uvis_flc_nocrrej_drc_sci.fits'
if do_cr:
    print(f'  CR rejected: {fits.getheader("wfc3_uvis_flc_cr_drc_sci.fits")["EXPTIME"]:.1f} s')
if do_nocrrej:
    print(f'  No CR rejection: {fits.getheader("wfc3_uvis_flc_nocrrej_drc_sci.fits")["EXPTIME"]:.1f} s')
exptime = fits.getheader(_primary_sci)['EXPTIME']
info_json.update(exptime_json_path, sample, lens, filt_key, exptime)

# ── Copy final sci/wht to output_path ────────────────────────────────────────
print('\n=== Copying final products to output directory ===')
_copy = []
if do_cr:
    _copy += ['wfc3_uvis_flc_cr_drc_sci.fits', 'wfc3_uvis_flc_cr_drc_wht.fits']
if do_nocrrej:
    _copy += ['wfc3_uvis_flc_nocrrej_drc_sci.fits', 'wfc3_uvis_flc_nocrrej_drc_wht.fits']
for fname in _copy:
    shutil.copy(fname, output_path)
    print(f'  {fname}')

# ── Plots ──────────────────────────────────────────────────────────────────────
print('\n=== Saving plots ===')

single_sci_files = sorted(glob.glob('*_single_sci.fits'))
single_wht_files = sorted(glob.glob('*_single_wht.fits'))

# ── Registration QC ──────────────────────────────────────────────────────────────
# A per-visit WCS solution can land far off even though --align mast trusts every
# delivered WCS as correct (e.g. an old GSC240 fit vs a GAIA-tied fit on another visit
# of the same lens). Catch it by checking, for every contributing frame, whether its
# own single-drizzle image actually has flux at the catalogue position -- a real
# mismatch shows up as one frame's peak there being far below its siblings'. (Modeled
# on the ACS QC block.)
if _lc is not None and len(single_sci_files) >= 2:
    print('\n=== Registration QC ===')
    _qc_box = 30  # +-1.2" around the catalogue position at 0.0396"/px
    _hdr0 = fits.getheader(single_sci_files[0])
    _xt, _yt = WCS(_hdr0).all_world2pix(_lc.ra.deg, _lc.dec.deg, 0)
    _xt, _yt = int(round(float(_xt))), int(round(float(_yt)))
    print(f'  checking {len(single_sci_files)} contributing frames at catalogue position...')
    _peaks = []
    for _f in single_sci_files:
        _data = fits.getdata(_f)
        _y0, _y1 = max(0, _yt - _qc_box), min(_data.shape[0], _yt + _qc_box)
        _x0, _x1 = max(0, _xt - _qc_box), min(_data.shape[1], _xt + _qc_box)
        _sub = _data[_y0:_y1, _x0:_x1]
        _peaks.append((os.path.basename(_f).replace('_single_sci.fits', ''),
                        float(np.nanmax(_sub)) if _sub.size else float('nan')))
    _med = float(np.nanmedian([_p for _, _p in _peaks])) if _peaks else float('nan')
    if not (_med > 0):
        print('  (zero/undefined median peak -- nothing to compare)')
    else:
        _bad = False
        for _name, _p in _peaks:
            _ratio = _p / _med if np.isfinite(_p) else float('nan')
            _flag = np.isfinite(_ratio) and _ratio < 0.35
            _bad = _bad or _flag
            print(f'  {_name}: peak={_p:.3g} ({_ratio:.2f}x median)' + ('  <-- WARNING' if _flag else ''))
        if _bad:
            print('  WARNING: possible cross-visit WCS mismatch -- inspect before '
                  'trusting this product (see CLAUDE.md "WCS alignment")')
        else:
            print('  OK - all frames consistent')

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
    save_drizzled_png('wfc3_uvis_flc_cr_drc_sci.fits', 'wfc3_uvis_flc_cr_drc_wht.fits',
                      os.path.join(output_path, 'drizzled_cr.png'), title='(CR rejection)')
    print('  drizzled_cr.png')

if do_nocrrej:
    save_drizzled_png('wfc3_uvis_flc_nocrrej_drc_sci.fits', 'wfc3_uvis_flc_nocrrej_drc_wht.fits',
                      os.path.join(output_path, 'drizzled_nocrrej.png'), title='(no CR rejection)')
    print('  drizzled_nocrrej.png')

if do_cr and do_nocrrej:
    sci_cr      = fits.getdata('wfc3_uvis_flc_cr_drc_sci.fits')
    sci_nocrrej = fits.getdata('wfc3_uvis_flc_nocrrej_drc_sci.fits')
    wht         = fits.getdata('wfc3_uvis_flc_cr_drc_wht.fits')
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
    axes[1].set_title('CR rejection (LACosmic)')
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
