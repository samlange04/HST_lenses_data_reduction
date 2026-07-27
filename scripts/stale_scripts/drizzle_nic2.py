#!/usr/bin/env python
"""
Drizzle HST NICMOS/NIC2 CAL images into a combined mosaic.
Downloads CAL files from MAST into data/calibrated/ if not already present.
Writes final sci/wht FITS to data/drizzled/, intermediates to data/drizzle_files/.
Produces a no-CR-rejection drizzle by default. Pass --cr to also run the
CR-rejected pass.

CAL files (_cal.fits) are the recommended input for NICMOS: they are the
per-exposure calibrated images produced by CALNICA. Unlike CCD instruments,
NICMOS uses HgCdTe detectors with up-the-ramp sampling, so there is no CTE
correction and no FLC equivalent. DQ bits are kept at 0 (only clean pixels)
because NICMOS DQ flags are less standardized than those of ACS/WFC3.
"""

import argparse
import json
import os
import glob
import shutil
import subprocess
import sys

# ── DEPRIORITISED — do not run ────────────────────────────────────────────────
# NIC2's field of view is far too small to be useful here (258x256 px at 0.0756" is
# ~19" across, against ~139" for WFC3/IR), and the pipeline may be unsound. All
# NICMOS data was deleted on 2026-07-21 — drizzled products, working dirs, 108
# *cal.fits exposures, run logs and the NICMOS CRDS cache, 472 MB — and the 24
# affected lenses now carry f160W: null in all three tracking JSONs.
#
# The guard exists because that deletion is silently reversible: this script
# re-downloads from MAST and re-fetches the CRDS refs automatically, so an
# accidental run brings all 472 MB back AND repopulates the 24 null entries,
# overwriting the record that they were dropped on purpose. The script itself is
# kept deliberately (F160W is answered from WFC3/IR), hence an override rather than
# deletion. Raised, not sys.exit()'d, so the failure is loud and non-zero.
if not os.environ.get('ALLOW_NICMOS'):
    raise NotImplementedError(
        'drizzle_nic2.py is deprioritised: NIC2 has too small a field of view to be '
        'useful and all NICMOS data was deliberately deleted (2026-07-21). Running '
        'this re-downloads ~472 MB from MAST and overwrites the f160W: null entries '
        'that record the deletion. Answer F160W questions from WFC3/IR '
        '(drizzle_wfc3_ir.py) instead. Set ALLOW_NICMOS=1 to override.'
    )

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import LogStretch, ImageNormalize
from astroquery.mast import Observations
from drizzlepac import tweakreg, astrodrizzle
from stwcs.updatewcs import updatewcs

# Resolve MAST target names (some lenses are archived under GAL-* not SDSS<LENS>).
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mast_target_names

# Route large FITS output writes through mmap+memcpy (vm_fault path) instead of
# fwrite/cluster_write copyin, to dodge the macOS U-state write-path lost-wakeup.
# Must run before AstroDrizzle. NICMOS outputs are small but still triggered the hang
# in testing, so this is required here too (not just for ACS/WFC3).
import mmap_fits_write
mmap_fits_write.install()

# ── Configuration ──────────────────────────────────────────────────────────────
_p = argparse.ArgumentParser()
_p.add_argument('--lens',        default='J0008-0004')
_p.add_argument('--filt',        default='f160W')
_p.add_argument('--sample',      default='slacs')
_p.add_argument('--cr',          action='store_true', default=False)
_p.add_argument('--_subprocess', action='store_true', default=False, help=argparse.SUPPRESS)
_a = _p.parse_args()

lens           = _a.lens
sample         = _a.sample
filt           = _a.filt
do_cr          = _a.cr
_is_subprocess = _a._subprocess

ws_path     = '/Users/samlange/Code/HST_lenses_data_reduction'
data_path   = os.path.join(ws_path, 'data', 'calibrated', sample, lens, filt)
output_path = os.path.join(ws_path, 'data', 'drizzled', sample, lens, filt)
work_path   = os.path.join(ws_path, 'data', 'drizzle_files', sample, lens, filt)
ref_path    = os.path.join(ws_path, 'data', 'reference_files')

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
    os.environ['nref']            = os.path.join(ref_path, 'references', 'hst', 'nicmos') + os.sep
    os.chdir(work_path)
    cal_files = sorted(glob.glob('*cal.fits'))
    print('\n=== AstroDrizzle (no CR rejection) ===')
    astrodrizzle.AstroDrizzle(cal_files,
                               output='nic2_cal_nocrrej',
                               preserve=False, build=False, context=False,
                               skysub=True, skymethod='localmin',
                               driz_sep_wcs=True, driz_sep_scale=0.0756,
                               driz_sep_bits='2,4,8', driz_sep_fillval=-1,
                               median=False, blot=False, driz_cr=False,
                               final_fillval=None, final_bits='2,4,8',
                               final_wcs=True, final_scale=0.0756,
                               num_cores=_num_cores)
    for f in glob.glob('*ask.fits'):
        os.remove(f)
    print('\n=== Cropping (no CR) ===')
    crop_to_coverage('nic2_cal_nocrrej_drz_sci.fits', 'nic2_cal_nocrrej_drz_wht.fits')
    sys.exit(0)

# ── Download CAL files from MAST (skip if already present) ────────────────────
with open(json_path) as _f:
    lens_products = json.load(_f)

if glob.glob(os.path.join(data_path, '*cal.fits')):
    print(f'=== MAST download: CAL files already present in {data_path}, skipping ===')
else:
    print(f'=== MAST download: querying {lens} {filt_key} NICMOS/NIC2 ===')
    os.makedirs(data_path, exist_ok=True)
    try:
        obs_table = None
        for _pat in mast_target_names.target_patterns(lens):
            obs_table = Observations.query_criteria(
                target_name=_pat,
                obs_collection='HST',
                instrument_name='NICMOS/NIC2',
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
            productSubGroupDescription=['CAL'],
            project=['CALNIC'],
        )
        for cal_file in glob.glob(os.path.join(data_path, 'mastDownload', '**', '*cal.fits'),
                                   recursive=True):
            shutil.move(cal_file, data_path)
        shutil.rmtree(os.path.join(data_path, 'mastDownload'), ignore_errors=True)

        obs_ids = sorted(
            os.path.basename(f).replace('_cal.fits', '')
            for f in glob.glob(os.path.join(data_path, '*cal.fits'))
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

if not glob.glob(os.path.join(data_path, '*cal.fits')):
    _update_info_json(exptime_json_path,    lens, filt_key, None)
    _update_info_json(instrument_json_path, lens, filt_key, None)
    sys.exit(f'No files found for {lens} {filt_key} — check target name and filter')

# Save instrument from first CAL header
with fits.open(sorted(glob.glob(os.path.join(data_path, '*cal.fits')))[0]) as _h:
    _instrume = _h[0].header['INSTRUME'].strip()
    _camera   = str(_h[0].header.get('CAMERA', 2))
_update_info_json(instrument_json_path, lens, filt_key, f'{_instrume}/NIC{_camera}')

# Skip drizzle if final products already exist
_skip_sentinel = 'nic2_cal_cr_drz_sci.fits' if do_cr else 'nic2_cal_nocrrej_drz_sci.fits'
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
os.environ['nref']            = os.path.join(
    ref_path, 'references', 'hst', 'nicmos') + os.sep

# ── Copy CAL files to work directory and work there ───────────────────────────
for f in glob.glob(os.path.join(data_path, '*cal.fits')):
    shutil.copy(f, work_path)

os.chdir(work_path)
print(f'Working directory: {work_path}')

# ── Download reference files ──────────────────────────────────────────────────
# Always run bestrefs against the actual input files: it is idempotent (CRDS only
# fetches missing refs) and cheap when they are present. A "skip if the ref dir is
# non-empty" guard was wrong across filters — one filter's run leaves another filter's
# refs unfetched, then updatewcs fails with a missing reference file.
print('\n=== CRDS bestrefs ===')
os.system('crds bestrefs --files *cal.fits --sync-references=1 --update-bestrefs')

# ── Update WCS ────────────────────────────────────────────────────────────────
print('\n=== updatewcs ===')
updatewcs('*cal.fits', use_db=False)

# ── TweakReg: align exposures ─────────────────────────────────────────────────
cal_files = sorted(glob.glob('*cal.fits'))

# Default source-finding / matching params, tuned for the small (256x256), faint
# (~0.03 counts/s), widely-dithered (~59 px) NIC2 F160W fields. The ACS defaults fail:
#  - threshold=200 finds only ~2 sources -> 2-point convex hull crashes TweakReg's
#    SphericalPolygon.from_radec ("Polygon made of too few points"); threshold=4 -> ~19.
#  - searchrad=1/tolerance=3 cross-match only 1 source -> nan shifts; searchrad=5/
#    tolerance=8 matches ~17 of ~20 sources on a typical field.
_tweak_params = dict(conv_width=2.5, threshold=4.0, searchrad=5, tolerance=8, minobj=5)
# Per-lens overrides for fields where the sample defaults still fail (see run 2026-07-18):
_TWEAK_OVERRIDES = {
    'J0252+0039': {'tolerance': 5},    # tol=8 -> stimage.xyxymatch "output coordinates exceeded allocation"
    'J1213+6708': {'threshold': 3.0},  # thr=4 -> one exposure has only 6 sources -> <minobj matches -> nan
}
_tweak_params.update(_TWEAK_OVERRIDES.get(lens, {}))

print('\n=== TweakReg ===')
if lens in _TWEAK_OVERRIDES:
    print(f'  Using per-lens TweakReg overrides for {lens}: {_TWEAK_OVERRIDES[lens]}')
tweakreg.TweakReg(cal_files,
                  updatehdr=True,
                  clean=True,
                  reusename=True,
                  interactive=False,
                  ylimit=0.2,
                  shiftfile=True,
                  outshifts='shift_cal.txt',
                  **_tweak_params)

with open('shift_cal.txt') as f:
    for i, line in enumerate(f, 1):
        if 'nan' in line:
            raise ValueError(f'nan in shift_cal.txt line {i} — TweakReg alignment failed')

# ── AstroDrizzle pass 1: with CR rejection ────────────────────────────────────
if do_cr:
    # DQ bits 2,4,8: uncertain linearity/dark/flat corrections — calibration
    # imperfections that are acceptable to include. All other flags are kept active:
    # 1=telemetry error, 16=grot, 32=defective pixel, 64=saturated, 128=missing
    # data, 256=bad cal pixel, 512=CR, 1024=source, 2048=0th-read signal,
    # 4096=CR from MultiDrizzle. (NICMOS Data Handbook, Table 2.3)
    print('\n=== AstroDrizzle (with CR rejection) ===')
    astrodrizzle.AstroDrizzle(cal_files,
                               output='nic2_cal_cr',
                               preserve=False, build=False, context=False,
                               skysub=True, skymethod='localmin',
                               driz_sep_wcs=True, driz_sep_scale=0.0756,
                               driz_sep_bits='2,4,8', driz_sep_fillval=-1,
                               median=True, blot=True, driz_cr=True,
                               driz_cr_snr='4.0 3.5', driz_cr_scale='1.2 0.7',
                               final_fillval=None, final_bits='2,4,8',
                               final_wcs=True, final_scale=0.0756,
                               num_cores=_num_cores)
    for f in glob.glob('*ask.fits'):
        os.remove(f)

    print('\n=== Cropping (CR) ===')
    crop_to_coverage('nic2_cal_cr_drz_sci.fits', 'nic2_cal_cr_drz_wht.fits')

    # ── AstroDrizzle pass 2: no CR rejection in subprocess to free memory ─────
    print('\n=== AstroDrizzle (no CR rejection) — launching subprocess ===')
    _result = subprocess.run(
        [sys.executable, os.path.abspath(__file__),
         '--lens', lens, '--filt', filt, '--sample', sample, '--_subprocess'],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
    )
    sys.stdout.write(_result.stdout)
    sys.stdout.flush()
    if _result.returncode != 0:
        raise RuntimeError(f'no-CR subprocess exited with code {_result.returncode}')

else:
    # ── AstroDrizzle pass 2: no CR rejection ──────────────────────────────────
    print('\n=== AstroDrizzle (no CR rejection) ===')
    astrodrizzle.AstroDrizzle(cal_files,
                               output='nic2_cal_nocrrej',
                               preserve=False, build=False, context=False,
                               skysub=True, skymethod='localmin',
                               driz_sep_wcs=True, driz_sep_scale=0.0756,
                               driz_sep_bits='2,4,8', driz_sep_fillval=-1,
                               median=False, blot=False, driz_cr=False,
                               final_fillval=None, final_bits='2,4,8',
                               final_wcs=True, final_scale=0.0756,
                               num_cores=_num_cores)
    for f in glob.glob('*ask.fits'):
        os.remove(f)

    print('\n=== Cropping (no CR) ===')
    crop_to_coverage('nic2_cal_nocrrej_drz_sci.fits', 'nic2_cal_nocrrej_drz_wht.fits')

print('\n=== Exposure times ===')
exptime = fits.getheader('nic2_cal_nocrrej_drz_sci.fits')['EXPTIME']
print(f'  No CR rejection: {exptime:.1f} s')
if do_cr:
    print(f'  CR rejected: {fits.getheader("nic2_cal_cr_drz_sci.fits")["EXPTIME"]:.1f} s')
_update_info_json(exptime_json_path, lens, filt_key, exptime)

# ── Copy final sci/wht to output_path ────────────────────────────────────────
print('\n=== Copying final products to output directory ===')
_copy = ['nic2_cal_nocrrej_drz_sci.fits', 'nic2_cal_nocrrej_drz_wht.fits']
if do_cr:
    _copy = ['nic2_cal_cr_drz_sci.fits', 'nic2_cal_cr_drz_wht.fits'] + _copy
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
    save_drizzled_png('nic2_cal_cr_drz_sci.fits', 'nic2_cal_cr_drz_wht.fits',
                      os.path.join(output_path, 'drizzled_cr.png'), title='(CR rejection)')
    print('  drizzled_cr.png')

save_drizzled_png('nic2_cal_nocrrej_drz_sci.fits', 'nic2_cal_nocrrej_drz_wht.fits',
                  os.path.join(output_path, 'drizzled_nocrrej.png'), title='(no CR rejection)')
print('  drizzled_nocrrej.png')

if do_cr:
    sci_cr      = fits.getdata('nic2_cal_cr_drz_sci.fits')
    sci_nocrrej = fits.getdata('nic2_cal_nocrrej_drz_sci.fits')
    wht         = fits.getdata('nic2_cal_cr_drz_wht.fits')
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
    axes[1].set_title('CR rejection (snr=4.0/3.5, scale=1.2/0.7)')
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
