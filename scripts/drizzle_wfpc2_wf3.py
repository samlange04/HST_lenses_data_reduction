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
_p.add_argument('--sample', default='slacs')
# WFPC2 KEEPS TweakReg, unlike ACS and WFC3/IR. It is the one instrument here that
# cannot skip updatewcs -- the NPOL/D2IM distortion arrays are required and the WF3
# chip extraction does not carry them -- and updatewcs strips the delivered
# astrometric fit (IDC_ta81040lu-FIT_IMG_GSC242 -> bare IDC_ta81040lu; use_db=True
# only restores an older GSC240 fit, which measured identically). With the fit gone,
# TweakReg is the only source of relative alignment. Measured on J0029-0055
# (stacked stellar FWHM, 32 vs 9 stars):
#     updatewcs + TweakReg     0.309"   <- default
#     updatewcs, no TweakReg   0.388"
# ACS/WFC3 are the opposite: they skip updatewcs, keep their FIT_REL solutions, and
# TweakReg only degrades them. Do not "unify" the three scripts on this point.
_p.add_argument('--align',   default='tweakreg', choices=['mast', 'tweakreg'],
                help="'tweakreg' (default for WFPC2) runs updatewcs + TweakReg; "
                     "'mast' skips TweakReg and trusts the delivered WCS, which is "
                     "correct for ACS/WFC3 but measurably worse here.")
# Drizzle output weight type. 'ERR' (default) makes the WHT extension a full
# inverse-variance map (source Poisson + sky + read + dark), so 1/sqrt(WHT) is a
# CALIBRATED per-pixel noise map -- what make_cutouts uses. AstroDrizzle's own
# default 'EXP' is only an effective-exposure-time map (uncalibrated, missing
# source shot noise). See DrizzlePac Handbook pp.103,139 and Bayer et al. 2023.
_p.add_argument('--wht-type',    default='ERR', choices=['ERR', 'IVM', 'EXP'])
_a = _p.parse_args()

lens   = _a.lens
sample = _a.sample
filt   = _a.filt

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
# pixfrac 0.8 measured on J0029-0055 (4 exposures): it is both the sharpest and
# well covered. 0.7 leaves genuine weight holes (std/mean 0.44, 8.4% of interior
# pixels below half-median weight) which degrade the PSF rather than improve it;
# 1.0 is the most uniform (std/mean 0.09) but softer.
#   pixfrac   std/mean   holes    median stellar FWHM
#   0.7       0.435      8.4%     4.38 px / 0.219"
#   0.8       0.282      1.3%     3.75 px / 0.188"   <- chosen
#   1.0       0.093      1.1%     4.05 px / 0.203"
DEFAULT_SCALE, DEFAULT_PIXFRAC = 0.05, 0.8


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
filt_key             = filt  # e.g. 'f606W'
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

# ── Download FLT files from MAST (skip if already present) ────────────────────
with open(json_path) as _f:
    lens_products = json.load(_f)

if glob.glob(os.path.join(data_path, 'u*flt.fits')):
    print(f'=== MAST download: FLT files already present in {data_path}, skipping ===')
else:
    print(f'=== MAST download: querying {lens} {filt_key} WFPC2 (PC/WFC apertures) ===')
    os.makedirs(data_path, exist_ok=True)
    try:
        obs_table = None
        for _pat in mast_target_names.target_patterns(lens):
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
        for _f in sorted(glob.glob(os.path.join(data_path, 'u*flt.fits'))):
            _exp = fits.getheader(_f)['EXPTIME']
            if _exp < MIN_EXPTIME:
                os.remove(_f)
                print(f'  rejected {os.path.basename(_f)}: EXPTIME={_exp}s < {MIN_EXPTIME}s')

        obs_ids = sorted(
            os.path.basename(f).replace('_flt.fits', '')
            for f in glob.glob(os.path.join(data_path, 'u*flt.fits'))
        )
        if lens not in lens_products:
            lens_products[lens] = {}
        lens_products[lens][filt_key] = obs_ids
        lens_products[lens] = dict(sorted(lens_products[lens].items()))
        with open(json_path, 'w') as _f:
            json.dump(dict(sorted(lens_products.items())), _f, indent=4)
        print(f'  Downloaded {len(obs_ids)} exposures, updated lens_products.json')
    except Exception as e:
        print(f'  MAST query failed: {e}')

if not glob.glob(os.path.join(data_path, 'u*flt.fits')):
    _update_info_json(exptime_json_path,    lens, filt_key, None)
    _update_info_json(instrument_json_path, lens, filt_key, None)
    sys.exit(f'No files found for {lens} {filt_key} — check target name and filter')

# Save instrument from first FLT header
with fits.open(sorted(glob.glob(os.path.join(data_path, 'u*flt.fits')))[0]) as _h:
    _instrume = _h[0].header['INSTRUME'].strip()
_update_info_json(instrument_json_path, lens, filt_key, f'{_instrume}/WF3')

# ── Choose output scale from the actual sub-pixel dither coverage ─────────────
_inputs = sorted(glob.glob(os.path.join(data_path, 'u*flt.fits')))
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
for f in glob.glob(os.path.join(data_path, 'u*flt.fits')):
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

# ── AstroDrizzle pass 1: with CR rejection ────────────────────────────────────
print('\n=== AstroDrizzle (with CR rejection) ===')
astrodrizzle.AstroDrizzle(flt_wf3_files,
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

print('\n=== Cropping (CR) ===')
crop_to_coverage('wfpc2_wf3_cr_drw_sci.fits', 'wfpc2_wf3_cr_drw_wht.fits')

# ── AstroDrizzle pass 2: no CR rejection ──────────────────────────────────────
print('\n=== AstroDrizzle (no CR rejection) ===')
astrodrizzle.AstroDrizzle(flt_wf3_files,
                           output='wfpc2_wf3_nocrrej',
                           preserve=False, build=False, context=False,
                           skysub=True, skymethod='localmin',
                           driz_sep_wcs=True, driz_sep_scale=WF3_NATIVE_SCALE,
                           driz_sep_bits='8,1024', driz_sep_fillval=-1,
                           median=False, blot=False, driz_cr=False,
                           final_fillval=None, final_bits='8,1024',
                           final_wcs=True, final_scale=out_scale,
                           final_pixfrac=out_pixfrac,
                           final_wht_type=_a.wht_type,
                           **_common_wcs,
                           num_cores=_num_cores)
for f in glob.glob('*ask.fits'):
    os.remove(f)

print('\n=== Cropping (no CR) ===')
crop_to_coverage('wfpc2_wf3_nocrrej_drw_sci.fits', 'wfpc2_wf3_nocrrej_drw_wht.fits')

print('\n=== Exposure times ===')
for label, fname in (('CR rejected', 'wfpc2_wf3_cr_drw_sci.fits'),
                     ('No CR rejection', 'wfpc2_wf3_nocrrej_drw_sci.fits')):
    exptime = fits.getheader(fname)['EXPTIME']
    print(f'  {label}: {exptime:.1f} s')
_update_info_json(exptime_json_path, lens, filt_key, exptime)

# ── Copy final sci/wht to output_path ────────────────────────────────────────
print('\n=== Copying final products to output directory ===')
for fname in ('wfpc2_wf3_cr_drw_sci.fits',     'wfpc2_wf3_cr_drw_wht.fits',
              'wfpc2_wf3_nocrrej_drw_sci.fits', 'wfpc2_wf3_nocrrej_drw_wht.fits'):
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
