#!/usr/bin/env python
"""
Drizzle HST WFPC2 WF3-chip FLT images into a combined mosaic.
Downloads FLT files from MAST into data/calibrated/ if not already present.
Writes final sci/wht FITS to data/drizzled/, intermediates to data/drizzle_files/.
Produces both a CR-rejected and a no-CR-rejection drizzle for comparison.

WHY WF3 AND NOT THE PC: for every SLACS lens with WFPC2 F606W data, the lens
galaxy falls on WF3 (extension 3), never on the PC. `DETECTOR = PC` in the
primary header — and hence the "WFPC2/PC" label on MAST — describes the aperture,
not the chip the target lands on; the standard WFPC2 full-field aperture centres
the target on WF3. Verified on pixel data: 327-774 sigma detections at the
WCS-predicted WF3 position while the PC chip centre is empty background.
The superseded PC version of this script produced blank-sky mosaics.
"""

import argparse
import json
import os
import glob
import shutil
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
from drizzlepac.wfpc2Data import wfpc2_to_flt
from stwcs.updatewcs import updatewcs

# Resolve MAST target names (some lenses are archived under GAL-* not SDSS<LENS>).
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mast_target_names
import info_json

# Route large FITS output writes through mmap+memcpy (vm_fault path) instead of
# fwrite/cluster_write copyin, to dodge the macOS U-state write-path lost-wakeup.
# Must run before AstroDrizzle. Output SIZE is not a reliable predictor of the wedge
# (NICMOS hung on a sub-MB write), so this is wired in here despite the modest WF3 output.
import mmap_fits_write
mmap_fits_write.install()

# ── Configuration ──────────────────────────────────────────────────────────────
_p = argparse.ArgumentParser()
_p.add_argument('--lens',   default='J0008-0004')
_p.add_argument('--filt',   default='f606W')
_p.add_argument('--sample', default=mast_target_names.DEFAULT_SAMPLE,
                help='sample the lens belongs to; sets the <sample> level of every '
                     'data/ path. Defined in info/lens_samples.json '
                     f'(default {mast_target_names.DEFAULT_SAMPLE})')
# WFPC2 runs updatewcs but NOT TweakReg. It is the one instrument here that cannot
# skip updatewcs -- the NPOL/D2IM distortion arrays are required and the WF3 chip
# extraction does not carry them -- so --align mast means updatewcs(use_db=True),
# which restores the GSC240 fit, and then stops. ACS/WFC3 skip updatewcs entirely.
# Do not "unify" the three scripts on this point.
#
# DEFAULT REVERSED 2026-07-26: was 'tweakreg'. The per-lens core-registration audit
# (info/wfpc2_alignment.json) put ALL 22 lenses on 'mast', but the default was never
# changed to match, so any direct run -- and the batch runner, which passed no
# --align at all -- silently used the mode the audit rejected. TweakReg aligns every
# frame onto the first, which for a dithered single-visit sequence measures the
# dither itself as error and removes it: frames scatter by ~0.7" and the deflector
# core splits into ~4 offset knots (visible on J0252+0039), against ~0.02" for MAST.
# The older FWHM comparison that favoured TweakReg (J0029-0055, 0.309" vs 0.388")
# was misleading -- it centroided a single extended galaxy, which rewards TweakReg's
# internal self-consistency even while the deflector is being split.
_p.add_argument('--align',   default='mast', choices=['mast', 'tweakreg'],
                help="'mast' (default) runs updatewcs and trusts the delivered WCS; "
                     "'tweakreg' additionally re-solves with TweakReg, which the "
                     "per-lens audit rejected for every lens -- comparison runs only.")
# Drizzle output weight type. 'IVM' (default HERE, unlike the other three drizzle
# scripts which default to 'ERR') makes the WHT extension a full inverse-variance
# map (source Poisson + sky + read + dark), so 1/sqrt(WHT) is a CALIBRATED
# per-pixel noise map -- what make_cutouts uses. AstroDrizzle's own default 'EXP'
# is only an effective-exposure-time map (uncalibrated, missing source shot
# noise). See DrizzlePac Handbook pp.103,139 and Bayer et al. 2023.
#
# WFPC2-SPECIFIC: 'ERR' does NOT work here. drizzlepac.wfpc2Data.WFPC2InputImage
# hardcodes self.errExt = None unconditionally (standard WFPC2 pipeline products
# never carry an ERR extension), so imageObject.buildERRmask() always takes the
# "WFPC2 not supported" branch and silently falls back to exposure-time-only
# weighting -- confirmed firing in every run.log ("No ERR weighting will be
# applied ... WFPC2 data is not supported by this weighting type"). Measured
# effect: the resulting noise map is flat, core/sky ratio 1.000 on J0330-0020 and
# J1213+6708, vs ~3.5 for ACS F814W on the same lenses (source Poisson correctly
# elevating the core). build_ivm_files() below builds a genuine per-frame IVM and
# DrizzlePac *does* honor a user-supplied IVM file. 'ERR' is kept as a choice only
# for comparison/debugging; do not use it expecting calibrated noise.
#
# DO NOT build the IVM from the file's own ERR array. It is not a real error array:
# ERR == sqrt(SCI) exactly (an identity over 100.00% of good pixels, not a fit),
# i.e. Poisson applied to DN as if DN were electrons -- no gain conversion, no
# read-noise term. It overstates the true noise by ~2.1x at sky level. The model
# used instead is var_DN = SCI/gain + floor^2, with the floor measured per frame.
_p.add_argument('--wht-type',    default='IVM', choices=['ERR', 'IVM', 'EXP'])
# CR-rejection method for the CR pass. 'lacosmic' (default) masks cosmic rays per
# frame with astroscrappy; 'drizcr' is the old AstroDrizzle median/blot/driz_cr route.
# WFPC2's driz_cr punches holes and masks real core pixels exactly like it did on ACS
# (biased blotted-median reference on a steep core), which spikes the noise map; the
# per-frame LACosmic route has no stacked reference to bias. See run_lacosmic below.
_p.add_argument('--cr-method',   default='lacosmic', choices=['lacosmic', 'drizcr'])
_p.add_argument('--lacosmic-sigclip', type=float, default=4.5)
_p.add_argument('--lacosmic-objlim',  type=float, default=5.0)
# The LACosmic CR pass is the science product for WFPC2 (LACosmic preserves the core, so
# the CR mosaic is what make_cutouts --pass auto cuts). The no-CR ("nocrrej") pass used to
# run unconditionally as a comparison; it is now opt-in via --nocrrej and off by default —
# it doubles the drizzle time and nothing downstream reads it.
_p.add_argument('--nocrrej', action=argparse.BooleanOptionalAction, default=False,
                help='also produce the no-CR comparison drizzle (off by default)')
# Single-visit drizzling for multi-visit lenses. --pa restricts the drizzle to frames
# within 1 deg of that PA_V3 (one guide-star solution -> clean MAST registration);
# --out-suffix tags the output/work dir (e.g. f606W -> f606W_v2) so the two visits do
# not overwrite each other -- empty (default) for the primary, longer-exptime visit,
# same as any other filter. Used for J0728+3835 / J0822+2652, whose two visits have a
# 14-15 deg roll split that TweakReg cannot co-register below ~0.3": drizzle each visit
# separately and let the modelling (PyAutoLens DatasetModel.grid_offset) fit the offset.
_p.add_argument('--pa',          type=float, default=None)
_p.add_argument('--out-suffix',  default='')
_a = _p.parse_args()

lens       = _a.lens
sample     = _a.sample
filt       = _a.filt
do_nocrrej = _a.nocrrej

# ── Output sampling ────────────────────────────────────────────────────────────
# WF3 native scale is 0.0996"/px. The WFPC2-BOX 4-point pattern (spacing 0.559")
# dithers by half a pixel in both axes, which supports drizzling to a ~2x finer
# grid; 0.05"/px also matches the ACS F814W mosaics, putting both bands on a
# common grid. A lens without sub-pixel phase coverage in *both* axes must stay at
# native scale or it aliases, so the scale is chosen per lens from the actual
# exposure WCSs (see dither_phase_counts) rather than from a hardcoded list —
# the set of exposures depends on which MAST instrument labels are queried, so a
# fixed table goes stale as soon as that changes.
WF3_NATIVE_SCALE = 0.0996

# Frames below this are pointing/CR check shots, not science: the WFC-labelled
# frames on J0728+3835, J0822+2652 and J1142+1001 are 0.5 s. This also catches the
# EXPTIME=0 frames that the old "prefer non-COPY" rule existed to avoid.
MIN_EXPTIME = 10.0

# Total-exposure-time gate, checked on the surviving (post-MIN_EXPTIME, post--pa)
# frame set. slacs_other runs generally shorter total exposures than slacs_gold;
# below BLOCK_EXPTIME no product is written (same outcome as no MAST data), between
# BLOCK and WARN the drizzle proceeds but is flagged.
WARN_EXPTIME, BLOCK_EXPTIME = 1200.0, 500.0

# pixfrac 1.0 (not 0.8): F606W keeps its dither-supported 0.05" scale (WFPC2-BOX is a
# 4-point HALF-pixel dither, purpose-built for 2x oversampling of the 0.0996" native --
# do NOT coarsen to 0.06"; that would waste the sub-pixel sampling the dither encodes).
# But with the noise now on calibrated ERR weight maps, the drizzle goal is a *uniform,
# low-correlation* noise map for the likelihood, not the sharpest PSF. pixfrac 1.0 fills
# the grid so adjacent-pixel noise correlation collapses. Measured on J0252+0039 (0.05",
# correctly MAST-registered frames):
#   pixfrac   noise texture   adjacent-pixel RMS
#   0.8       3.9%            5.3%
#   1.0       2.1%            2.1%   <- chosen
# 0.06"/1.0 was fractionally smoother (1.9%) but sacrifices resolution for ~0.2% and is
# not the dither-matched scale, so it is rejected. Matches the WFC3/IR F160W choice
# (0.06"/1.0). The old pixfrac-0.8 tuning (stacked FWHM on J0029) predates ERR weighting
# and the shift to prioritising noise-map covariance; PyAutoLens fits the PSF explicitly.
DEFAULT_SCALE, DEFAULT_PIXFRAC = 0.05, 1.0


def dither_phase_counts(flt_files, ext=3, ref_pix=(400.0, 400.0)):
    """Number of distinct sub-pixel dither phases per axis, in WF3 pixels.

    POSTARG1/2 are zero for these datasets — WFPC2 pattern offsets are applied by
    the spacecraft and land in the WCS, not in POSTARG — so measure the offsets by
    mapping a fixed sky position through each exposure's WF3 WCS. Two or more
    distinct phases on an axis means that axis is sub-pixel sampled and can be
    drizzled to a finer grid; one phase means it cannot.
    """
    with fits.open(flt_files[0]) as hl0:
        ref = WCS(hl0[ext].header, hl0).pixel_to_world(*ref_pix)
    fx, fy = set(), set()
    for fname in flt_files:
        with fits.open(fname) as hl:
            x, y = WCS(hl[ext].header, hl).world_to_pixel(ref)
        fx.add(round((x - ref_pix[0]) % 1, 1))
        fy.add(round((y - ref_pix[1]) % 1, 1))
    # phase 1.0 is phase 0.0
    norm = lambda s: len({0.0 if v in (0.0, 1.0) else v for v in s})
    return norm(fx), norm(fy)

ws_path     = '/Users/samlange/Code/HST_lenses_data_reduction'
# data_path (calibrated source) always uses the base filter; only the output/work dirs
# take --out-suffix, so both visits share the one download but write to f606W / _v2.
data_path   = os.path.join(ws_path, 'data', 'calibrated', sample, lens, filt)
output_path = os.path.join(ws_path, 'data', 'drizzled', sample, lens, filt + _a.out_suffix)
work_path   = os.path.join(ws_path, 'data', 'drizzle_files', sample, lens, filt + _a.out_suffix)
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
filt_key             = filt  # e.g. 'f606W' — the MAST filter name, used for querying
# Tracking-JSON key. It mirrors the product directory, so a split-visit run
# (--out-suffix _v2) is recorded under 'f606W_v2', not 'f606W' -- only the primary
# (longer-exptime, empty-suffix) visit takes the bare filter key. Writing every visit
# under the bare filter name is what once left J0728+3835 claiming a combined
# 6600s/6-obsid f606W product that has never existed on disk (the real product is
# now 'f606W', 4400s/4 frames — its 2-frame visit is dropped for want of dither phase).
# That error was invisible precisely because the key was present and plausible.
product_key          = filt + _a.out_suffix
json_path            = os.path.join(ws_path, 'info', 'lens_products.json')
exptime_json_path    = os.path.join(ws_path, 'info', 'lens_exptime.json')
instrument_json_path = os.path.join(ws_path, 'info', 'lens_instrument.json')

# ── Download FLT files from MAST (skip if already present) ────────────────────
# Distinguishes "MAST has nothing for this lens+filter" (an ordinary outcome -- every
# lens in a sample is tried on every run, most have no data in most bands) from "the
# download broke" (a real failure). Conflating them is what would let a network error
# be silently recorded as `null` = no data. See the no-data exit below.
_mast_empty = False
_mast_error = None
_all_too_short = False

if glob.glob(os.path.join(data_path, 'u*flt.fits')):
    print(f'=== MAST download: FLT files already present in {data_path}, skipping ===')
else:
    print(f'=== MAST download: querying {lens} {filt_key} WFPC2 (PC/WFC apertures) ===')
    os.makedirs(data_path, exist_ok=True)
    try:
        obs_table = None
        for _pat in mast_target_names.target_patterns(lens, sample):
            obs_table = Observations.query_criteria(
                target_name=_pat,
                obs_collection='HST',
                instrument_name=['WFPC2/PC', 'WFPC2/WFC'],
                filters=[filt_key.upper()],
            )
            if len(obs_table) > 0:
                print(f'  MAST target {_pat}: {len(obs_table)} observations')
                break
            print(f'  MAST target {_pat}: no observations')
        if obs_table is None or len(obs_table) == 0:
            # Not an error: this lens simply has no WFPC2 data in this filter. 16 of the
            # 38 slacs_gold lenses are in exactly this position.
            _mast_empty = True
            raise mast_target_names.NoMastData()
        # Keep COPY *and* non-COPY. For WFPC2 SLACS these are genuine repeat visits
        # months apart, not archive duplicates of the same photons (J0728+3835:
        # 2 x 1100s on 2007-09-14 plus 4 x 1100s "-COPY" on 2007-11-05; distinct
        # t_min for every frame), so preferring non-COPY threw away two thirds of
        # the usable data and left too few dither phases to reach 0.05"/px.
        # Unusable frames are rejected below on exposure time instead, which also
        # covers the EXPTIME=0 case the non-COPY preference was originally for.
        print(f'  Using all {len(obs_table)} observations (COPY and non-COPY)')
        products = Observations.get_product_list(obs_table)
        # C0M/C1M as well as FLT: MAST only generated FLT products for the
        # observations it labelled WFPC2/PC. The WFPC2/WFC ones carry just the raw
        # GEIS-style C0M (science) + C1M (DQ) pair, so without these the extra
        # exposures are invisible — J1218+0830 has 4 x 1100s frames but only the
        # 2 PC-labelled ones ship an FLT. They are converted to FLT below.
        Observations.download_products(
            products,
            download_dir=data_path,
            productSubGroupDescription=['FLT', 'C0M', 'C1M'],
            project=['CALWFPC2'],
        )
        for pat in ('*flt.fits', '*c0m.fits', '*c1m.fits'):
            for f in glob.glob(os.path.join(data_path, 'mastDownload', '**', pat),
                               recursive=True):
                shutil.move(f, data_path)
        shutil.rmtree(os.path.join(data_path, 'mastDownload'), ignore_errors=True)

        # Build FLT files for any exposure MAST shipped only as C0M/C1M.
        for _c0m in sorted(glob.glob(os.path.join(data_path, '*c0m.fits'))):
            _flt = _c0m.replace('c0m', 'flt')
            if os.path.exists(_flt):
                continue
            if not os.path.exists(_c0m.replace('c0m', 'c1m')):
                print(f'  {os.path.basename(_c0m)}: no matching c1m, skipping conversion')
                continue
            wfpc2_to_flt(_c0m)
            print(f'  converted {os.path.basename(_c0m)} -> {os.path.basename(_flt)}')

        # Drop check shots / zero-exposure frames before they reach the drizzle.
        _n_before = len(glob.glob(os.path.join(data_path, 'u*flt.fits')))
        for _f in sorted(glob.glob(os.path.join(data_path, 'u*flt.fits'))):
            _exp = fits.getheader(_f)['EXPTIME']
            if _exp < MIN_EXPTIME:
                os.remove(_f)
                print(f'  rejected {os.path.basename(_f)}: EXPTIME={_exp}s < {MIN_EXPTIME}s')
        _n_after = len(glob.glob(os.path.join(data_path, 'u*flt.fits')))
        # Every frame was a check shot. Also an ordinary outcome, not a broken download:
        # some lenses' only WFPC2 F606W frames are 0.5s pointing verifications.
        if _n_before and not _n_after:
            _all_too_short = True

        print(f'  Downloaded {_n_after} exposures')
    except mast_target_names.NoMastData:
        print(f'  No WFPC2 {filt_key.upper()} observations for {lens} on MAST')
    except Exception as e:
        _mast_error = e
        print(f'  MAST query failed: {e}')

if not glob.glob(os.path.join(data_path, 'u*flt.fits')):
    info_json.update(exptime_json_path,    sample, lens, product_key, None)
    info_json.update(instrument_json_path, sample, lens, product_key, None)
    if _mast_empty or _all_too_short:
        # Ordinary outcome, not a failure: exit 0 so a batch runner sweeping the whole
        # sample records "no data" and moves on instead of counting it as an error.
        _why = (f'no WFPC2 {filt_key.upper()} on MAST' if _mast_empty
                else f'every frame is shorter than MIN_EXPTIME={MIN_EXPTIME}s')
        print(f'=== NO DATA: {lens} — {_why} (recorded as null) ===')
        sys.exit(0)
    if _mast_error is not None:
        sys.exit(f'MAST download failed for {lens} {filt_key}: {_mast_error}')
    sys.exit(f'No files found for {lens} {filt_key} — MAST listed observations but no '
             f'FLT files landed in {data_path}')

# Save instrument from first FLT header
with fits.open(sorted(glob.glob(os.path.join(data_path, 'u*flt.fits')))[0]) as _h:
    _instrume = _h[0].header['INSTRUME'].strip()
info_json.update(instrument_json_path, sample, lens, product_key, f'{_instrume}/WF3')

# ── Choose output scale from the actual sub-pixel dither coverage ─────────────
_inputs = sorted(glob.glob(os.path.join(data_path, 'u*flt.fits')))
if _a.pa is not None:
    # keep only frames within 1 deg of the requested visit roll (--pa single-visit mode)
    _inputs = [f for f in _inputs
               if abs(float(fits.getheader(f).get('PA_V3', 1e9)) - _a.pa) < 1.0]
    print(f'=== Single-visit: {len(_inputs)} frames within 1deg of PA_V3={_a.pa} ===')

_total_exptime = sum(fits.getheader(f)['EXPTIME'] for f in _inputs)
if _total_exptime < BLOCK_EXPTIME:
    info_json.update(exptime_json_path,    sample, lens, product_key, None)
    info_json.update(instrument_json_path, sample, lens, product_key, None)
    print(f'=== BLOCKED (exptime): {lens} {product_key} total exptime '
          f'{_total_exptime:.1f}s < {BLOCK_EXPTIME:.0f}s minimum (recorded as null) ===')
    sys.exit(0)
if _total_exptime < WARN_EXPTIME:
    print(f'  EXPTIME WARNING: {lens} {product_key} total exptime {_total_exptime:.1f}s '
          f'< {WARN_EXPTIME:.0f}s (proceeding)')

# Record provenance from the frames that actually reach the drizzle, not from
# everything the MAST download left in data/calibrated/. Those differ whenever --pa
# selects one visit out of two, and recording the download instead is what made
# lens_products.json list all 6 J0728+3835 obsids against a 4-frame product. Written
# here rather than inside the download block so that a re-run on already-downloaded
# data still refreshes it; a lens that exits below for want of dither phase writes
# nothing, which is correct — no product exists for it.
info_json.update(json_path, sample, lens, product_key, sorted(
    os.path.basename(f).replace('_flt.fits', '') for f in _inputs
))

_nx, _ny = dither_phase_counts(_inputs)
print(f'=== Dither sampling: {len(_inputs)} exposures, '
      f'{_nx} x-phases / {_ny} y-phases ===')
# Reaching 0.05"/px needs sub-pixel phase coverage on both axes. Without it the
# only honest option is the native 0.0996"/px grid, and a native-scale mosaic is
# not wanted here — skip the lens rather than emit a product at the wrong scale.
if min(_nx, _ny) < 2:
    sys.exit(f'{lens} {filt_key}: only {_nx} x-phase(s) / {_ny} y-phase(s) across '
             f'{len(_inputs)} exposures — cannot drizzle to {DEFAULT_SCALE}"/px '
             f'without aliasing. Skipping (no product written).')
out_scale, out_pixfrac = DEFAULT_SCALE, DEFAULT_PIXFRAC
print(f'=== Output grid: final_scale={out_scale}"/px, pixfrac={out_pixfrac} ===')

# Skip drizzle if final products already exist
if os.path.exists(os.path.join(output_path, 'wfpc2_wf3_cr_drw_sci.fits')):
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
os.environ['uref']            = os.path.join(
    ref_path, 'references', 'hst', 'wfpc2') + os.sep

# ── Copy FLT files to work directory and work there ───────────────────────────
# Copy only the selected frames (_inputs is PA-filtered in single-visit mode), so
# every downstream step that globs u*flt.fits in the work dir sees just this visit.
for f in _inputs:
    shutil.copy(f, work_path)

os.chdir(work_path)
print(f'Working directory: {work_path}')

# ── Download reference files ──────────────────────────────────────────────────
# Always run bestrefs against the actual input files: it is idempotent (CRDS only
# fetches missing refs) and cheap when they are present. A "skip if the ref dir is
# non-empty" guard was wrong across filters — one filter's run leaves another filter's
# refs unfetched, then updatewcs fails with a missing reference file.
print('\n=== CRDS bestrefs ===')
os.system('crds bestrefs --files u*flt.fits --sync-references=1 --update-bestrefs')

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
# updatewcs runs for WFPC2 in BOTH modes, unlike ACS/WFC3 where it is skipped. It is
# not optional here: it attaches the non-polynomial (NPOLFILE) and detector-to-image
# (D2IMFILE) distortion arrays, and the WF3 chip extraction below does not carry them
# over, so without it AstroDrizzle stops and prompts for the missing DGEO correction.
# It costs a little relative accuracy (spread 0.19 -> 0.38 px on J0029-0055) but that
# is far cheaper than TweakReg (-> 1.54 px).
# use_db matters here. use_db=False strips the delivered astrometric fit
# (IDC_ta81040lu-FIT_IMG_GSC242 -> bare IDC_ta81040lu), which leaves relative
# alignment resting on the IDCTAB alone -- that is why dropping TweakReg without
# this made J0029-0055 worse (FWHM 0.309" -> 0.388"). use_db=True pulls a fitted
# solution back from the astrometry database (IDC_ta81040lu-GSC240).
print('\n=== updatewcs ===')
updatewcs('u*flt.fits', use_db=(_a.align == 'mast'))

if _a.align == 'tweakreg':
    # ── TweakReg: align exposures ─────────────────────────────────────────────
    print('\n=== TweakReg ===')
    tweakreg.TweakReg(sorted(glob.glob('u*flt.fits')),
                      updatehdr=True,
                      clean=True,
                      reusename=True,
                      interactive=False,
                      conv_width=3.0,
                      # 200 was a PC-chip value and is too strict once WF3 sources
                      # dominate the catalogue: on J0822+2652 it left only the very
                      # brightest objects, so 1 of 6 frames went nan and another matched
                      # spuriously at (40, -65) px. Measured on that lens, 100 and 50 both
                      # match 6/6 with a coherent ~9.5 px visit offset; 20 and below let
                      # false matches back in (max |shift| 52-69 px). 100 is the
                      # conservative end of the working range.
                      threshold=100.0,
                      ylimit=1,
                      shiftfile=True,
                      outshifts='shift_flt.txt',
                      # 1" was enough when every exposure came from one visit. Combining
                      # the COPY and non-COPY visits means matching across a roll
                      # difference (J0822+2652: PA_V3 101.85 vs 87.92) and a separate
                      # guide-star solution, which needs shifts up to ~0.85" — right at
                      # a 1" radius, where 2 of 6 frames failed to match and TweakReg
                      # wrote nan shifts. 3" is still tight against the 80" WF3 field.
                      searchrad=3,
                      tolerance=3,
                      minobj=7)

    with open('shift_flt.txt') as f:
        for i, line in enumerate(f, 1):
            if 'nan' in line:
                raise ValueError(f'nan in shift_flt.txt line {i} — TweakReg alignment failed')
else:
    _names = set()
    for _f in sorted(glob.glob('u*flt.fits')):
        try:
            _names.add(fits.getval(_f, 'WCSNAME', ext=('SCI', 1)))
        except Exception:
            _names.add('<none>')
    print(f'\n=== Alignment: using MAST WCS as delivered ({", ".join(sorted(_names))}) ===')
    if not all(('GSC' in n or 'GAIA' in n or 'HSC' in n) for n in _names):
        print('  WARNING: a frame lacks a fitted (GSC/GAIA/HSC) WCS. Relative alignment is unverified —')
        print('           check the stacked PSF, or re-run with --align tweakreg.')

# ── Extract WF3 chip only ──────────────────────────────────────────────────────
print('\n=== Masking PC / WF2 / WF4 chips ===')

def extract_wf3_chip(flt_files):
    """Extract only the WF3 chip (SCI,3 / DQ,3 / ERR,3) from each flt file.
    Prefix with 'wf3_' so the filename still ends in _flt.fits for DrizzlePac."""
    wf3_files = []
    for fname in flt_files:
        out = 'wf3_' + os.path.basename(fname)
        with fits.open(fname) as hdul:
            new_hdul = fits.HDUList([hdul[0].copy()])
            new_hdul.append(hdul['SCI', 3].copy())
            new_hdul.append(hdul['DQ',  3].copy())
            new_hdul.append(hdul['ERR', 3].copy())
            for ext in hdul:
                if ext.name in ('D2IMARR', 'WCSCORR'):
                    new_hdul.append(ext.copy())
            # DrizzlePac indexes chips positionally as (SCI, 1..numchips), so the
            # lone remaining extension must be EXTVER=1 or WFPC2InputImage raises
            # KeyError: Extension ('SCI', 1) not found. This does NOT mislabel the
            # chip: detnum comes from the DETECTOR keyword (stwcs instruments.py
            # set_chip), which stays 3, so the WF3 gain/readnoise row and the "WF3"
            # detector name are still the ones used.
            for ext in new_hdul[1:4]:
                ext.header['EXTVER'] = 1
            new_hdul.writeto(out, overwrite=True)
        print(f'  {fname} -> {out}')
        wf3_files.append(out)
    return wf3_files

flt_wf3_files = extract_wf3_chip(sorted(glob.glob('u*flt.fits')))

# ── Build per-frame IVM files from the real ERR array ──────────────────────────
# See the --wht-type help above: DrizzlePac's WFPC2 driver hardcodes errExt=None
# and never reads ERR for WFPC2, no matter what's in the file. IVM = 1/ERR^2 is
# the inverse-variance-map DrizzlePac *does* support as a user-supplied file, fed
# via a two-column @-file (irafglob's ivmlist convention: 'atfile_ivm' reads the
# second whitespace-separated word per line). ERR itself doesn't change between
# the CR and no-CR passes (only DQ does), so these are built once and reused.
WF3_READNOISE = 5.2    # electrons, WFPC2 WF3 chip (WFPC2 Instrument Handbook Table 4.2)


def measure_noise_floor(sci, dq, gain):
    """Measure the frame's own additive variance floor, in DN^2.

    The Poisson term of the noise model is physics (counts/gain), but the additive
    term is not just read noise: dark current and flat-field error land there too,
    and on these frames the total measures ~7.4 e against a nominal RN of 5.2 e. So
    it is measured per frame rather than assumed.

    The estimator is the MAD of the second difference along x,
    0.5*(x[i-1] + x[i+1]) - x[i], which has variance 1.5*sigma^2. Two properties
    matter: it cancels linear gradients (so the galaxy and any sky slope do not
    inflate it), and it involves no selection on pixel *value* -- clipping pixels
    at +-N MAD and then re-measuring the MAD of the survivors biases the width low
    by ~4%, which is how an earlier version of this measurement overstated the
    error in the delivered ERR array as 2.28x when it is 2.11x.
    """
    ok = (dq[:, :-2] == 0) & (dq[:, 1:-1] == 0) & (dq[:, 2:] == 0)
    d = (0.5 * (sci[:, :-2] + sci[:, 2:]) - sci[:, 1:-1])[ok]
    sky = float(np.median(sci[dq == 0]))
    var = (1.4826 * np.median(np.abs(d - np.median(d)))) ** 2 / 1.5
    floor2 = var - sky / gain
    if not np.isfinite(floor2) or floor2 <= 0:
        # Measured sky variance below the Poisson floor is unphysical; fall back to
        # nominal read noise rather than silently producing an over-confident IVM.
        floor2 = (WF3_READNOISE / gain) ** 2
        print(f'    WARNING: measured floor non-positive, using nominal RN '
              f'{WF3_READNOISE} e -> {floor2:.3f} DN^2')
    return floor2, sky, var


# See the --wht-type help above for why a user-supplied IVM file is the only route
# that works for WFPC2 (DrizzlePac's WFPC2 driver hardcodes errExt=None). It is fed
# via a two-column @-file (irafglob's ivmlist convention: 'atfile_ivm' reads the
# second whitespace-separated word per line). The model depends only on SCI and DQ,
# and DQ changes between the CR and no-CR passes only in bit 4096, which does not
# enter the variance -- so these are built once and reused.
#
# The IVM is built from a NOISE MODEL, not from the file's ERR array. The ERR array
# delivered in these WFPC2 FLTs is exactly sqrt(SCI) -- verified as an identity, not
# a fit: ERR/sqrt(SCI) == 1.000000 for 100.00% of good pixels. That is Poisson
# statistics applied to DATA NUMBERS as though they were electrons, so it omits the
# gain conversion entirely and carries no read-noise term at all. It overstates the
# true noise by ~2.1x at sky level (see the ERR-array section in CLAUDE.md).
def build_ivm_files(wf3_files):
    """Write '<sci> -> IVM,1' files with IVM = 1/(SCI/gain + floor^2), in DN^-2."""
    ivm_files = []
    for fname in wf3_files:
        with fits.open(fname) as hdul:
            sci  = hdul['SCI', 1].data.astype(np.float64)
            dq   = hdul['DQ', 1].data
            gain = float(hdul[0].header.get('ATODGAIN', 7.0))

        floor2, sky, skyvar = measure_noise_floor(sci, dq, gain)
        # Poisson on total collected counts (sky + source; bias and dark are already
        # subtracted by CALWFPC2, so their variance sits in the floor). Negative
        # pixels are noise excursions about zero, not negative variance.
        var = np.clip(sci, 0.0, None) / gain + floor2
        ivm = np.zeros(sci.shape, dtype=np.float32)
        good = np.isfinite(var) & (var > 0)
        ivm[good] = (1.0 / var[good]).astype(np.float32)

        out = 'ivm_' + os.path.basename(fname)
        ivm_hdu = fits.ImageHDU(data=ivm, name='IVM')
        ivm_hdu.header['EXTVER']   = 1
        ivm_hdu.header['IVMGAIN']  = (gain, 'e/DN used for the Poisson term')
        ivm_hdu.header['IVMFLOOR'] = (floor2, 'DN^2 additive variance floor (measured)')
        ivm_hdu.header['IVMSKY']   = (sky, 'DN median sky level of the frame')
        fits.HDUList([fits.PrimaryHDU(), ivm_hdu]).writeto(out, overwrite=True)
        ivm_files.append(out)

        print(f'  {os.path.basename(fname)}: sky {sky:7.2f} DN, gain {gain:.2f} e/DN, '
              f'measured sky var {skyvar:6.3f} = Poisson {sky / gain:6.3f} + floor '
              f'{floor2:5.3f} DN^2 (floor {np.sqrt(floor2) * gain:.1f} e)')
    return ivm_files

if _a.wht_type == 'IVM':
    print('\n=== Building per-frame IVM files (1/ERR^2) — WFPC2 ERR weighting is unsupported by DrizzlePac ===')
    _ivm_files = build_ivm_files(flt_wf3_files)
    _ivm_assoc_path = 'wf3_ivm_association.lst'
    with open(_ivm_assoc_path, 'w') as _f:
        for _sci_f, _ivm_f in zip(flt_wf3_files, _ivm_files):
            _f.write(f'{_sci_f} {_ivm_f}\n')
    drizzle_input = '@' + _ivm_assoc_path
    print(f'  Built {len(_ivm_files)} IVM files -> {_ivm_assoc_path}')
else:
    drizzle_input = flt_wf3_files

_num_cores = 1

# ── Helpers ────────────────────────────────────────────────────────────────────
def coverage_bbox(wht_file):
    """Bounding box (y0, y1, x0, x1) of wht>0 in a weight map."""
    with fits.open(wht_file) as h:
        wht = h[0].data
        rows = np.any(wht > 0, axis=1)
        cols = np.any(wht > 0, axis=0)
        y0, y1 = np.where(rows)[0][[0, -1]]
        x0, x1 = np.where(cols)[0][[0, -1]]
    return int(y0), int(y1), int(x0), int(x1)


def crop_to_bbox(sci_file, wht_file, bbox):
    """Trim sci and wht to a given (y0, y1, x0, x1) bbox, updating CRPIX in-place."""
    y0, y1, x0, x1 = bbox
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

# ── CR rejection: LACosmic (default), same rationale as drizzle_acs_wfc.py ─────
# driz_cr compares each frame to a blotted median of the stack; on a steep PSF core
# that reference is biased low, so real core pixels read as CRs (measured -37% core
# flux on ACS, and it spikes the WFPC2 noise map with masked patches + a CR trail).
# LACosmic works one frame at a time with an object-protection term (objlim), so no
# stacked reference biases the core, and the final drizzle is a plain weighted mean.
CR_BIT = 4096          # not in final_bits='8,1024', so flagged CRs are excluded
# WF3_READNOISE is defined above build_ivm_files, which runs first and needs it.


def run_lacosmic(files, sigclip, objlim):
    """Flag cosmic rays into DQ bit 4096 of each WF3 FLT, in place.

    WF3 FLT is in DN (BUNIT=COUNTS): gain=ATODGAIN converts to electrons internally,
    satlevel=SATURATE*gain is the full well in electrons, readnoise is the WF3 value.
    The extracted chip is the lone (SCI,1)/(DQ,1) pair. Saturated pixels (DQ bit 8)
    are protected so the bright core is not mistaken for cosmic rays.
    """
    import astroscrappy
    total = 0
    for fname in files:
        with fits.open(fname, mode='update') as hdul:
            gain  = float(hdul[0].header.get('ATODGAIN', 7.0))
            satdn = float(hdul[0].header.get('SATURATE', 4095))
            sci = hdul['SCI', 1].data.astype(np.float32)
            dq  = hdul['DQ', 1].data
            dq &= ~CR_BIT                      # clear any previous CR flags
            bad = (dq & 8) > 0                 # protect full-well saturated pixels
            mask, _ = astroscrappy.detect_cosmics(
                sci, inmask=bad, sigclip=sigclip, sigfrac=0.3, objlim=objlim,
                gain=gain, readnoise=WF3_READNOISE, satlevel=satdn * gain,
                niter=4, sepmed=True, cleantype='medmask', fsmode='median')
            dq[mask] |= CR_BIT
            total += int(mask.sum())
            hdul['DQ', 1].data = dq
    print(f'  LACosmic flagged {total} pixels across {len(files)} frames')


# ── AstroDrizzle pass 1: with CR rejection ────────────────────────────────────
if _a.cr_method == 'lacosmic':
    print('\n=== LACosmic CR masking ===')
    run_lacosmic(flt_wf3_files, _a.lacosmic_sigclip, _a.lacosmic_objlim)
    print('\n=== AstroDrizzle (LACosmic-masked, plain weighted mean) ===')
    astrodrizzle.AstroDrizzle(drizzle_input,
                               output='wfpc2_wf3_cr',
                               preserve=False, build=False, context=False,
                               skysub=True, skymethod='localmin',
                               driz_sep_wcs=True, driz_sep_scale=WF3_NATIVE_SCALE,
                               driz_sep_bits='8,1024', driz_sep_fillval=-1,
                               median=False, blot=False, driz_cr=False,
                               # resetbits defaults to 4096 and would clear the very
                               # bit LACosmic just wrote, silently reverting to an
                               # unmasked drizzle that still looks plausible.
                               resetbits=0,
                               final_fillval=None, final_bits='8,1024',
                               final_wcs=True, final_scale=out_scale,
                               final_pixfrac=out_pixfrac,
                               final_wht_type=_a.wht_type,
                               **_common_wcs,
                               num_cores=_num_cores)
else:
    print('\n=== AstroDrizzle (with driz_cr CR rejection) ===')
    astrodrizzle.AstroDrizzle(drizzle_input,
                               output='wfpc2_wf3_cr',
                               preserve=False, build=False, context=False,
                               skysub=True, skymethod='localmin',
                               driz_sep_wcs=True, driz_sep_scale=WF3_NATIVE_SCALE,
                               driz_sep_bits='8,1024', driz_sep_fillval=-1,
                               median=True, blot=True, driz_cr=True,
                               driz_cr_snr='15.0 10.0', driz_cr_scale='1.5 1.0',
                               final_fillval=None, final_bits='8,1024',
                               final_wcs=True, final_scale=out_scale,
                               final_pixfrac=out_pixfrac,
                               final_wht_type=_a.wht_type,
                               **_common_wcs,
                               num_cores=_num_cores)
for f in glob.glob('*ask.fits'):
    os.remove(f)

# ── AstroDrizzle pass 2: no CR rejection (opt-in via --nocrrej) ───────────────
if do_nocrrej:
    print('\n=== AstroDrizzle (no CR rejection) ===')
    astrodrizzle.AstroDrizzle(drizzle_input,
                               output='wfpc2_wf3_nocrrej',
                               preserve=False, build=False, context=False,
                               skysub=True, skymethod='localmin',
                               driz_sep_wcs=True, driz_sep_scale=WF3_NATIVE_SCALE,
                               driz_sep_bits='8,1024', driz_sep_fillval=-1,
                               median=False, blot=False, driz_cr=False,
                               # this pass runs on the same FLTs after the LACosmic pass
                               # wrote DQ bit 4096; reset it so the no-CR mosaic is genuinely
                               # un-masked and not silently CR-rejected.
                               resetbits=4096,
                               final_fillval=None, final_bits='8,1024',
                               final_wcs=True, final_scale=out_scale,
                               final_pixfrac=out_pixfrac,
                               final_wht_type=_a.wht_type,
                               **_common_wcs,
                               num_cores=_num_cores)
    for f in glob.glob('*ask.fits'):
        os.remove(f)

if do_nocrrej:
    print('\n=== Cropping (shared bbox for both passes) ===')
    # Crop BOTH passes to the union of their wht>0 boxes. LACosmic masks different edge
    # pixels in the CR pass than the no-CR pass, so per-pass cropping left the two products
    # differing by ~1 column and broke the residual diagnostic (and mis-registers the passes
    # pixel-for-pixel). The union box keeps every covered pixel and gives identical shapes.
    _bb_cr    = coverage_bbox('wfpc2_wf3_cr_drw_wht.fits')
    _bb_nocr  = coverage_bbox('wfpc2_wf3_nocrrej_drw_wht.fits')
    _bbox = (min(_bb_cr[0], _bb_nocr[0]), max(_bb_cr[1], _bb_nocr[1]),
             min(_bb_cr[2], _bb_nocr[2]), max(_bb_cr[3], _bb_nocr[3]))
    crop_to_bbox('wfpc2_wf3_cr_drw_sci.fits',     'wfpc2_wf3_cr_drw_wht.fits',     _bbox)
    crop_to_bbox('wfpc2_wf3_nocrrej_drw_sci.fits', 'wfpc2_wf3_nocrrej_drw_wht.fits', _bbox)
else:
    print('\n=== Cropping (CR pass) ===')
    _bbox = coverage_bbox('wfpc2_wf3_cr_drw_wht.fits')
    crop_to_bbox('wfpc2_wf3_cr_drw_sci.fits', 'wfpc2_wf3_cr_drw_wht.fits', _bbox)

print('\n=== Exposure times ===')
_passes = [('CR rejected', 'wfpc2_wf3_cr_drw_sci.fits')]
if do_nocrrej:
    _passes.append(('No CR rejection', 'wfpc2_wf3_nocrrej_drw_sci.fits'))
for label, fname in _passes:
    print(f'  {label}: {fits.getheader(fname)["EXPTIME"]:.1f} s')
# Record from the CR pass (the science product); EXPTIME is identical between passes.
exptime = fits.getheader('wfpc2_wf3_cr_drw_sci.fits')['EXPTIME']
info_json.update(exptime_json_path, sample, lens, product_key, exptime)

# ── Copy final sci/wht to output_path ────────────────────────────────────────
# Stamp the noise model into the products. Nothing else distinguishes a calibrated
# WFPC2 weight map from the two earlier uncalibrated generations (exptime-only, and
# IVM built from the bogus ERR=sqrt(SCI) array), and they are not separable by
# inspection -- so make_cutouts.py keys its "do not use for a likelihood" warning on
# this keyword rather than on the instrument name.
_ivm_model = ('SCI/gain+floor^2' if _a.wht_type == 'IVM' else _a.wht_type)
print('\n=== Copying final products to output directory ===')
_copy_files = ['wfpc2_wf3_cr_drw_sci.fits', 'wfpc2_wf3_cr_drw_wht.fits']
if do_nocrrej:
    _copy_files += ['wfpc2_wf3_nocrrej_drw_sci.fits', 'wfpc2_wf3_nocrrej_drw_wht.fits']
for fname in _copy_files:
    with fits.open(fname, mode='update') as _h:
        _h[0].header['IVMMODEL'] = (_ivm_model, 'WFPC2 noise model behind the WHT map')
        # AstroDrizzle writes count RATES (D001OUUN='cps') but leaves BUNIT at the
        # value inherited from the input FLT, which for WFPC2 is 'COUNTS' (DN). Left
        # alone, the product claims counts while holding DN/s — an EXPTIME-sized
        # (4400x) error for anything that reads BUNIT to set units. ACS/WFC3 are not
        # affected: their FLT/FLC BUNIT already says ELECTRONS/S or is rewritten.
        # Only the SCI files: the WHT map is correctly 'UNITLESS'.
        if fname.endswith('_sci.fits'):
            _h[0].header['BUNIT'] = ('COUNTS/S', 'DN per second (drizzle final_units=cps)')
    shutil.copy(fname, output_path)
    print(f'  {fname}')

# ── Plots ──────────────────────────────────────────────────────────────────────
print('\n=== Saving plots ===')

# Individual single-drizzle frames
single_sci_files = sorted(glob.glob('wf3_*_single_sci.fits'))
single_wht_files = sorted(glob.glob('wf3_*_single_wht.fits'))
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

save_drizzled_png('wfpc2_wf3_cr_drw_sci.fits', 'wfpc2_wf3_cr_drw_wht.fits',
                  os.path.join(output_path, 'drizzled_cr.png'), title='(CR rejection)')
print('  drizzled_cr.png')

if do_nocrrej:
    save_drizzled_png('wfpc2_wf3_nocrrej_drw_sci.fits', 'wfpc2_wf3_nocrrej_drw_wht.fits',
                      os.path.join(output_path, 'drizzled_nocrrej.png'), title='(no CR rejection)')
    print('  drizzled_nocrrej.png')

    # 3-panel comparison with shared SCI normalization
    sci_cr      = fits.getdata('wfpc2_wf3_cr_drw_sci.fits')
    sci_nocrrej = fits.getdata('wfpc2_wf3_nocrrej_drw_sci.fits')
    wht         = fits.getdata('wfpc2_wf3_cr_drw_wht.fits')
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
    axes[1].set_title('CR rejection (snr=15/10, scale=1.5/1.0)')
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
