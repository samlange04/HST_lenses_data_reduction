#!/usr/bin/env python
"""
Drizzle HST WFPC2 PC-chip FLT images into a combined mosaic.
Downloads FLT files from MAST into data/calibrated/ if not already present.
Writes final sci/wht FITS to data/drizzled/, intermediates to data/drizzle_files/.
Produces both a CR-rejected and a no-CR-rejection drizzle for comparison.
"""

import argparse
import json
import os
import glob
import shutil
import sys

# ── SUPERSEDED — do not run ───────────────────────────────────────────────────
# This script extracts the PC chip, but the SLACS lenses all fall on WF3; it
# produced 22 mosaics of blank sky ~79" from the target. Use drizzle_wfpc2_wf3.py.
# The guard matters because this script rmtree's data/drizzled/<lens>/<filt>/
# before drizzling, so running it would delete the good WF3 products.
# Raised rather than sys.exit()'d so that an accidental run fails loudly with a
# traceback and a non-zero status that a batch runner cannot mistake for a clean skip.
if not os.environ.get('ALLOW_SUPERSEDED_WFPC2_PC'):
    raise NotImplementedError(
        'drizzle_wfpc2_pc.py is superseded: the SLACS lenses sit on WF3, not the PC, '
        'so this script writes blank-sky mosaics and deletes the WF3 products on its '
        'way there. Use scripts/drizzle_wfpc2_wf3.py instead. '
        'Set ALLOW_SUPERSEDED_WFPC2_PC=1 to override.'
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
# Must run before AstroDrizzle. Output SIZE is not a reliable predictor of the wedge
# (NICMOS hung on a sub-MB write), so this is wired in here despite the small PC output.
import mmap_fits_write
mmap_fits_write.install()

# ── Configuration ──────────────────────────────────────────────────────────────
_p = argparse.ArgumentParser()
_p.add_argument('--lens',   default='J0008-0004')
_p.add_argument('--filt',   default='f606W')
_p.add_argument('--sample', default='slacs')
_a = _p.parse_args()

lens   = _a.lens
sample = _a.sample
filt   = _a.filt

ws_path     = '/Users/samlange/Code/HST_lenses_data_reduction'
data_path   = os.path.join(ws_path, 'data', 'calibrated', sample, lens, filt)
output_path = os.path.join(ws_path, 'data', 'drizzled', sample, lens, filt)
work_path   = os.path.join(ws_path, 'data', 'drizzle_files', sample, lens, filt)
ref_path    = os.path.join(ws_path, 'data', 'reference_files')

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
    print(f'=== MAST download: querying {lens} {filt_key} WFPC2/PC ===')
    os.makedirs(data_path, exist_ok=True)
    try:
        obs_table = None
        for _pat in mast_target_names.target_patterns(lens):
            obs_table = Observations.query_criteria(
                target_name=_pat,
                obs_collection='HST',
                instrument_name='WFPC2/PC',
                filters=[filt_key.upper()],
            )
            if len(obs_table) > 0:
                print(f'  MAST target {_pat}: {len(obs_table)} observations')
                break
            print(f'  MAST target {_pat}: no observations')
        # Prefer non-COPY observations; fall back to COPY only if nothing else exists
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
            project=['CALWFPC2'],
        )
        for flt_file in glob.glob(os.path.join(data_path, 'mastDownload', '**', '*flt.fits'),
                                   recursive=True):
            shutil.move(flt_file, data_path)
        shutil.rmtree(os.path.join(data_path, 'mastDownload'), ignore_errors=True)

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
_update_info_json(instrument_json_path, lens, filt_key, f'{_instrume}/PC')

# Skip drizzle if final products already exist
if os.path.exists(os.path.join(output_path, 'wfpc2_flt_cr_drw_sci.fits')):
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

# ── Update WCS ────────────────────────────────────────────────────────────────
print('\n=== updatewcs ===')
updatewcs('u*flt.fits', use_db=False)

# ── TweakReg: align exposures ─────────────────────────────────────────────────
print('\n=== TweakReg ===')
tweakreg.TweakReg(sorted(glob.glob('u*flt.fits')),
                  updatehdr=True,
                  clean=True,
                  reusename=True,
                  interactive=False,
                  conv_width=3.0,
                  threshold=200.0,
                  ylimit=1,
                  shiftfile=True,
                  outshifts='shift_flt.txt',
                  searchrad=1,
                  tolerance=3,
                  minobj=7)

with open('shift_flt.txt') as f:
    for i, line in enumerate(f, 1):
        if 'nan' in line:
            raise ValueError(f'nan in shift_flt.txt line {i} — TweakReg alignment failed')

# ── Extract PC chip only ───────────────────────────────────────────────────────
print('\n=== Masking WF chips ===')

def extract_pc_chip(flt_files):
    """Extract only the PC chip (SCI,1 / DQ,1 / ERR,1) from each flt file.
    Prefix with 'pc_' so the filename still ends in _flt.fits for DrizzlePac."""
    pc_files = []
    for fname in flt_files:
        out = 'pc_' + os.path.basename(fname)
        with fits.open(fname) as hdul:
            new_hdul = fits.HDUList([hdul[0].copy()])
            new_hdul.append(hdul['SCI', 1].copy())
            new_hdul.append(hdul['DQ',  1].copy())
            new_hdul.append(hdul['ERR', 1].copy())
            for ext in hdul:
                if ext.name in ('D2IMARR', 'WCSCORR'):
                    new_hdul.append(ext.copy())
            new_hdul.writeto(out, overwrite=True)
        print(f'  {fname} -> {out}')
        pc_files.append(out)
    return pc_files

flt_pc_files = extract_pc_chip(sorted(glob.glob('u*flt.fits')))

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
astrodrizzle.AstroDrizzle(flt_pc_files,
                           output='wfpc2_flt_cr',
                           preserve=False, build=False, context=False,
                           skysub=True, skymethod='localmin',
                           driz_sep_wcs=True, driz_sep_scale=0.0455,
                           driz_sep_bits='8,1024', driz_sep_fillval=-1,
                           median=True, blot=True, driz_cr=True,
                           driz_cr_snr='15.0 10.0', driz_cr_scale='1.5 1.0',
                           final_fillval=None, final_bits='8,1024',
                           final_wcs=True, final_scale=0.0455,
                           num_cores=_num_cores)
for f in glob.glob('*ask.fits'):
    os.remove(f)

print('\n=== Cropping (CR) ===')
crop_to_coverage('wfpc2_flt_cr_drw_sci.fits', 'wfpc2_flt_cr_drw_wht.fits')

# ── AstroDrizzle pass 2: no CR rejection ──────────────────────────────────────
print('\n=== AstroDrizzle (no CR rejection) ===')
astrodrizzle.AstroDrizzle(flt_pc_files,
                           output='wfpc2_flt_nocrrej',
                           preserve=False, build=False, context=False,
                           skysub=True, skymethod='localmin',
                           driz_sep_wcs=True, driz_sep_scale=0.0455,
                           driz_sep_bits='8,1024', driz_sep_fillval=-1,
                           median=False, blot=False, driz_cr=False,
                           final_fillval=None, final_bits='8,1024',
                           final_wcs=True, final_scale=0.0455,
                           num_cores=_num_cores)
for f in glob.glob('*ask.fits'):
    os.remove(f)

print('\n=== Cropping (no CR) ===')
crop_to_coverage('wfpc2_flt_nocrrej_drw_sci.fits', 'wfpc2_flt_nocrrej_drw_wht.fits')

print('\n=== Exposure times ===')
for label, fname in (('CR rejected', 'wfpc2_flt_cr_drw_sci.fits'),
                     ('No CR rejection', 'wfpc2_flt_nocrrej_drw_sci.fits')):
    exptime = fits.getheader(fname)['EXPTIME']
    print(f'  {label}: {exptime:.1f} s')
_update_info_json(exptime_json_path, lens, filt_key, exptime)

# ── Copy final sci/wht to output_path ────────────────────────────────────────
print('\n=== Copying final products to output directory ===')
for fname in ('wfpc2_flt_cr_drw_sci.fits',     'wfpc2_flt_cr_drw_wht.fits',
              'wfpc2_flt_nocrrej_drw_sci.fits', 'wfpc2_flt_nocrrej_drw_wht.fits'):
    shutil.copy(fname, output_path)
    print(f'  {fname}')

# ── Plots ──────────────────────────────────────────────────────────────────────
print('\n=== Saving plots ===')

# Individual single-drizzle frames
single_sci_files = sorted(glob.glob('pc_*_single_sci.fits'))
single_wht_files = sorted(glob.glob('pc_*_single_wht.fits'))
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

save_drizzled_png('wfpc2_flt_cr_drw_sci.fits', 'wfpc2_flt_cr_drw_wht.fits',
                  os.path.join(output_path, 'drizzled_cr.png'), title='(CR rejection)')
print('  drizzled_cr.png')

save_drizzled_png('wfpc2_flt_nocrrej_drw_sci.fits', 'wfpc2_flt_nocrrej_drw_wht.fits',
                  os.path.join(output_path, 'drizzled_nocrrej.png'), title='(no CR rejection)')
print('  drizzled_nocrrej.png')

# 3-panel comparison with shared SCI normalization
sci_cr      = fits.getdata('wfpc2_flt_cr_drw_sci.fits')
sci_nocrrej = fits.getdata('wfpc2_flt_nocrrej_drw_sci.fits')
wht         = fits.getdata('wfpc2_flt_cr_drw_wht.fits')
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
