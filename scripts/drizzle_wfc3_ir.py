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
_p.add_argument('--cr',          action='store_true', default=False)
_p.add_argument('--_subprocess', action='store_true', default=False, help=argparse.SUPPRESS)
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
# standard WFC3/IR choice (CANDELS, 3D-HST). Measured on J0728+3835 (4 exposures):
#   scale   pixfrac   wht std/mean   holes
#   0.1283  1.0       0.044          0.00%
#   0.0800  0.8       0.144          0.03%
#   0.0650  0.8       0.176          0.13%
#   0.0600  0.8       0.188          0.20%   <- chosen
#   0.0600  0.7       0.288          2.45%
# pixfrac 0.7 is clearly worse than 0.8 at either scale.
# Note this leaves F160W on a 0.06" grid while F606W/F814W are on 0.05".
IR_OUT_SCALE, IR_OUT_PIXFRAC = 0.06, 0.8

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
    flt_files = sorted(glob.glob('*flt.fits'))
    print('\n=== AstroDrizzle (no CR rejection) ===')
    astrodrizzle.AstroDrizzle(flt_files,
                               output='wfc3_ir_flt_nocrrej',
                               preserve=False, build=False, context=False,
                               skysub=True, skymethod='localmin',
                               driz_sep_wcs=True, driz_sep_scale=0.1283,
                               driz_sep_bits='64,512', driz_sep_fillval=-1,
                               median=False, blot=False, driz_cr=False,
                               final_fillval=None, final_bits='64,512',
                               final_wcs=True, final_scale=IR_OUT_SCALE,
                               final_pixfrac=IR_OUT_PIXFRAC,
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

        obs_ids = sorted(
            os.path.basename(f).replace('_flt.fits', '')
            for f in glob.glob(os.path.join(data_path, '*flt.fits'))
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

if not glob.glob(os.path.join(data_path, '*flt.fits')):
    _update_info_json(exptime_json_path,    lens, filt_key, None)
    _update_info_json(instrument_json_path, lens, filt_key, None)
    sys.exit(f'No files found for {lens} {filt_key} — check target name and filter')

# Save instrument from first FLT header
with fits.open(sorted(glob.glob(os.path.join(data_path, '*flt.fits')))[0]) as _h:
    _instrume = _h[0].header['INSTRUME'].strip()
    _detector = _h[0].header.get('DETECTOR', 'IR').strip()
_update_info_json(instrument_json_path, lens, filt_key, f'{_instrume}/{_detector}')

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

# ── Download reference files ──────────────────────────────────────────────────
# Always run bestrefs against the actual input files: it is idempotent (CRDS only
# fetches missing refs) and cheap when they are present. A "skip if the ref dir is
# non-empty" guard was wrong across filters — one filter's run leaves another filter's
# refs unfetched, then updatewcs fails with a missing reference file.
print('\n=== CRDS bestrefs ===')
os.system('crds bestrefs --files *flt.fits --sync-references=1 --update-bestrefs')

# ── Update WCS ────────────────────────────────────────────────────────────────
print('\n=== updatewcs ===')
updatewcs('*flt.fits', use_db=False)

# ── TweakReg: align exposures ─────────────────────────────────────────────────
flt_files = sorted(glob.glob('*flt.fits'))

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

# ── AstroDrizzle pass 1: with CR rejection ────────────────────────────────────
if do_cr:
    # DQ bits 64,512: warm pixels and blob pixels (flat-field artifacts dithered over).
    print('\n=== AstroDrizzle (with CR rejection) ===')
    astrodrizzle.AstroDrizzle(flt_files,
                               output='wfc3_ir_flt_cr',
                               preserve=False, build=False, context=False,
                               skysub=True, skymethod='localmin',
                               driz_sep_wcs=True, driz_sep_scale=0.1283,
                               driz_sep_bits='64,512', driz_sep_fillval=-1,
                               median=True, blot=True, driz_cr=True,
                               driz_cr_snr='3.5 3.0', driz_cr_scale='1.2 0.7',
                               final_fillval=None, final_bits='64,512',
                               final_wcs=True, final_scale=IR_OUT_SCALE,
                               final_pixfrac=IR_OUT_PIXFRAC,
                               num_cores=_num_cores)
    for f in glob.glob('*ask.fits'):
        os.remove(f)

    print('\n=== Cropping (CR) ===')
    crop_to_coverage('wfc3_ir_flt_cr_drz_sci.fits', 'wfc3_ir_flt_cr_drz_wht.fits')

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
    astrodrizzle.AstroDrizzle(flt_files,
                               output='wfc3_ir_flt_nocrrej',
                               preserve=False, build=False, context=False,
                               skysub=True, skymethod='localmin',
                               driz_sep_wcs=True, driz_sep_scale=0.1283,
                               driz_sep_bits='64,512', driz_sep_fillval=-1,
                               median=False, blot=False, driz_cr=False,
                               final_fillval=None, final_bits='64,512',
                               final_wcs=True, final_scale=IR_OUT_SCALE,
                               final_pixfrac=IR_OUT_PIXFRAC,
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
