#!/usr/bin/env python
"""
Drizzle-broadened model PSF by artificial-star injection (the Anderson 2016 route).

The model-PSF tier in make_psf.py (STDPSF, the ACS focus-diverse ePSF, the WFPC2 F606W
MAST-DB ePSF) is a *detector-frame* ePSF: psf_models.py resamples and rotates it to the
North-up output scale analytically, which reproduces the orientation and pixel scale but
NOT the extra blur AstroDrizzle's resampling puts on a point source. So the analytic model
kernel runs slightly sharp (e.g. detector-frame F160W ~3.15px FWHM vs the drizzle-broadened
empirical build ~3.8px). The rigorous fix, noted in the WFC3/IR ISR (Anderson 2016) and in
the psf_models docstring, is to inject the model PSF as artificial stars into the individual
exposures and drizzle them exactly as the science frames were drizzled; the drizzled star
then carries the resampling broadening for free, plus the correct North-up orientation and
the exposure-average weighting -- all produced by the real drizzle rather than emulated.

This is a *parallel* product for the model-PSF tier only (empirical builds already carry the
broadening -- they are cut from the drizzled mosaic itself). It writes alongside make_psf's
outputs and changes nothing downstream:

  data/psf/<sample>/<lens>/<filt>/psf_kernel_injected.fits   full image-scale kernel
  data/psf/<sample>/<lens>/<filt>/psf_injected.png           QA panel
  data/cutouts/<sample>/<lens>/<filt>/cutout_[cr_]psf_injected.fits   trimmed modelling kernel
  info/lens_psf_injected.json                                {sample:{lens:{filt:{...}}}}

How it works, and why it is faithful:
  * It reuses the already-prepared inputs the science drizzle consumed, which persist in
    data/drizzle_files/<sample>/<lens>/<filt>/ -- ACS `<root>_flc.fits`, WFC3/IR
    `<root>_flt.fits`, WFPC2 extracted `wf3_<root>_flt.fits` + per-frame IVM + the two-column
    `@`-association. Reusing them inherits every instrument-specific prep step (WF3 chip
    extraction + EXTVER renumbering, distortion arrays, updatewcs, the IVM weighting) that
    the drizzle scripts do -- there is no need to reconstruct any of it here.
  * For each input frame it ZEROES the SCI and adds the detector-frame model PSF at the lens
    position on that frame (via the frame's own WCS), keeping ERR / DQ / IVM untouched. The
    weight arrays are unchanged, so the per-frame drizzle weighting is identical to science;
    the injected point source is clean (no CRs), so no CR pass is needed.
  * It re-drizzles onto the SAME output grid the science product used -- final_rot=0,
    final_ra/dec at the lens, and the per-instrument final_scale / final_pixfrac / final_bits
    / final_wht_type read straight from the drizzle scripts. The drizzled star is then the
    drizzle-broadened, North-up PSF at the modelling pixel scale.

Requires the persisted data/drizzle_files/ inputs (present after a normal drizzle run; the
scripts only delete them on an explicit force re-run). If they are gone, re-run the band's
drizzle first. Records `null` + exits 0 on no data / no model tier, like the other scripts.

Usage:
    conda run -n stenv python scripts/make_psf_inject.py --lens J0252+0039 --filt f814W
    conda run -n stenv python scripts/make_psf_inject.py --lens J0728+3835 --filt f606W_v2
"""

import argparse
import glob
import json
import os
import shutil
import sys

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
import astropy.units as u
from photutils.centroids import centroid_com
from scipy.ndimage import map_coordinates

ws_path = '/Users/samlange/Code/HST_lenses_data_reduction'
sys.path.insert(0, os.path.join(ws_path, 'info'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# macOS AstroDrizzle write-hang workaround -- must be installed before AstroDrizzle runs
# (see CLAUDE.md / mmap_fits_write). No-op off macOS.
import mmap_fits_write
mmap_fits_write.install()
from drizzlepac import astrodrizzle

from slacs_coords import slacs_coords
import mast_target_names
import psf_models
# Reuse make_psf's kernel post-processing + QA + JSON helpers verbatim so the injected
# product is treated identically to the analytic one (importing make_psf is safe: its work
# is under __main__). instrument_key / representative field sizes also come from there.
import make_psf
from make_cutouts import find_products, _has_products
import info_json


# ── Per-instrument drizzle configuration ─────────────────────────────────────────
# Every value is lifted directly from the matching drizzle script so the injected re-drizzle
# reproduces the science output grid exactly. `input` selects how inputs are staged;
# `det_scale` is the detector pixel scale the model PSF is rendered at. resetbits matches the
# pass we emulate: 0 keeps the science CR mask (ACS/WFPC2 CR pass), 4096 clears stale CR
# flags on F160W (which has no CR pass). final_pixfrac None -> omit (ACS keeps the default).
_DRIZ = {
    'ACS/WFC': dict(input='flc', suffix='drc', det_scale=0.05,
                    driz_sep_scale=0.05, driz_sep_bits='256,64,16',
                    final_scale=0.05, final_pixfrac=None, final_bits='256,64,16',
                    wht_type='ERR', resetbits=0),
    'WFC3/IR': dict(input='flt', suffix='drz', det_scale=0.1283,
                    driz_sep_scale=0.1283, driz_sep_bits='512',
                    final_scale=0.06, final_pixfrac=1.0, final_bits='512',
                    wht_type='ERR', resetbits=4096),
    'WFPC2':   dict(input='wf3', suffix='drw', det_scale=0.0996,
                    driz_sep_scale=0.0996, driz_sep_bits='8,1024',
                    final_scale=0.05, final_pixfrac=1.0, final_bits='8,1024',
                    wht_type='IVM', resetbits=0),
}

_DET_STAMP = 51          # side (detector px) of the rendered model-PSF stamp injected per frame
_INJECT_FLUX = 1.0e4     # arbitrary positive total flux per injected star (kernel is renormalised)


def _set_ref_env():
    """Point CRDS / instrument ref-dir env vars at the cached reference files, so AstroDrizzle
    resolves the distortion references (IDCTAB etc.) at drizzle time -- same as the drizzle
    scripts. The files are already synced from the science run; no bestrefs needed here."""
    ref_path = os.path.join(ws_path, 'data', 'reference_files')
    os.environ['CRDS_SERVER_URL'] = 'https://hst-crds.stsci.edu'
    os.environ['CRDS_PATH'] = ref_path
    for inst in ('acs', 'wfc3', 'wfpc2'):
        var = {'acs': 'jref', 'wfc3': 'iref', 'wfpc2': 'uref'}[inst]
        os.environ[var] = os.path.join(ref_path, 'references', 'hst', inst) + os.sep


# ── Staging the prepared inputs from data/drizzle_files/ ──────────────────────────
def stage_inputs(inst_key, src_dir, work_dir):
    """Copy the science-prepared inputs into a fresh work dir and return
    (drizzle_input, sci_paths, ivm_of).

    `sci_paths` are the copied SCI files to inject into; `drizzle_input` is what AstroDrizzle
    is called with (a file list for ACS/WFC3, the copied `@`-association for WFPC2); `ivm_of`
    maps a WFPC2 sci path to its IVM file (empty for ACS/WFC3). Raises FileNotFoundError with
    a clear message if the prepared inputs are missing (drizzle_files was cleared)."""
    os.makedirs(work_dir, exist_ok=True)
    cfg = _DRIZ[inst_key]

    if cfg['input'] in ('flc', 'flt'):
        pat = '*_flc.fits' if cfg['input'] == 'flc' else '*_flt.fits'
        srcs = sorted(glob.glob(os.path.join(src_dir, pat)))
        # Drop EXPTIME=0 frames -- exactly what the science ACS/WFC3 drizzle does (AstroDrizzle
        # silently drops them; lens_products.json records only the survivors). A dud frame left
        # in drizzle_files (e.g. a zero-exposure exposure MAST returned) would otherwise be
        # staged, and for ACS its focus-diverse ePSF lookup fails -> the whole lens wrongly
        # falls back to STDPSF (this is exactly what happened to J1213+6708 f814W).
        kept = []
        for s in srcs:
            try:
                et = float(fits.getheader(s).get('EXPTIME', 0) or 0)
            except Exception:
                et = 0.0
            if et > 0:
                kept.append(s)
            else:
                print(f'  skipping EXPTIME=0 frame {os.path.basename(s)} '
                      f'(dropped by the science drizzle)')
        srcs = kept
    else:  # wf3: read the science IVM association for the exact (sci, ivm) pairs
        assoc = os.path.join(src_dir, 'wf3_ivm_association.lst')
        if not os.path.isfile(assoc):
            raise FileNotFoundError(
                f'no wf3_ivm_association.lst in {src_dir}; re-run the WFPC2 drizzle for this '
                f'lens (data/drizzle_files inputs are required by the injection route).')
        pairs = []
        with open(assoc) as f:
            for line in f:
                parts = line.split()
                if len(parts) == 2:
                    pairs.append(parts)
        sci_paths, ivm_of = [], {}
        for sci_name, ivm_name in pairs:
            s = shutil.copy(os.path.join(src_dir, sci_name), work_dir)
            shutil.copy(os.path.join(src_dir, ivm_name), work_dir)
            sci_paths.append(s)
            ivm_of[s] = os.path.basename(ivm_name)
        assoc_out = os.path.join(work_dir, 'wf3_ivm_association.lst')
        with open(assoc_out, 'w') as f:
            for s in sci_paths:
                f.write(f'{os.path.basename(s)} {ivm_of[s]}\n')
        return '@wf3_ivm_association.lst', sci_paths, ivm_of

    if not srcs:
        raise FileNotFoundError(
            f'no {cfg["input"]} inputs in {src_dir}; re-run the drizzle for this lens '
            f'(data/drizzle_files inputs are required by the injection route).')
    sci_paths = [shutil.copy(s, work_dir) for s in srcs]
    return [os.path.basename(s) for s in sci_paths], sci_paths, {}


# ── Rendering the detector-frame model PSF onto a frame ───────────────────────────
def _render_from_grid(grid, x, y, flux, S=_DET_STAMP):
    """(stamp, iy0, ix0): a photutils GriddedPSFModel (STDPSF) rendered at subpixel (x, y)
    on the detector, unit-sum x flux. grid.evaluate handles the subpixel phase directly."""
    half = S // 2
    ix, iy = int(round(x)), int(round(y))
    xs = np.arange(ix - half, ix + half + 1)
    ys = np.arange(iy - half, iy + half + 1)
    xx, yy = np.meshgrid(xs, ys)
    stamp = np.asarray(grid.evaluate(xx, yy, flux=1.0, x_0=x, y_0=y), float)
    s = stamp.sum()
    return (stamp / s * flux if s > 0 else stamp), iy - half, ix - half


def _render_from_array(P, ov, x, y, flux, S=_DET_STAMP):
    """(stamp, iy0, ix0): a supersampled detector-frame ePSF array `P` (ov source pixels per
    detector pixel) rendered at subpixel (x, y) on the detector, unit-sum x flux. Samples P
    with cubic interpolation; the detector and ePSF-source frames share axes (both detector),
    so the map is a pure scale by `ov` plus the subpixel offset."""
    P = np.asarray(P, float)
    ny, nx = P.shape
    py, px = np.unravel_index(np.argmax(P), P.shape)
    w = 15
    y0, x0 = max(0, py - w // 2), max(0, px - w // 2)
    cy, cx = centroid_com(P[y0:y0 + w, x0:x0 + w])
    cx, cy = x0 + cx, y0 + cy
    if not (np.isfinite(cx) and np.isfinite(cy)):
        cx, cy = nx / 2.0, ny / 2.0

    half = S // 2
    ix, iy = int(round(x)), int(round(y))
    offs = np.arange(-half, half + 1)
    XX, YY = np.meshgrid(offs, offs)
    # detector offset from (ix, iy) -> source-array coord: cx + (detcol - x)*ov
    srcx = cx + (XX + (ix - x)) * ov
    srcy = cy + (YY + (iy - y)) * ov
    stamp = map_coordinates(P, [srcy.ravel(), srcx.ravel()], order=3,
                            mode='constant', cval=0.0).reshape(S, S)
    stamp = np.clip(stamp, 0.0, None)
    s = stamp.sum()
    return (stamp / s * flux if s > 0 else stamp), iy - half, ix - half


def build_renderers(inst_key, filt, sample, lens, sci_paths, src_dir, catalogue_coord):
    """Return {sci_path: renderer} where renderer(sci_hdr_wcs_shape) -> (stamp, iy0, ix0, ext),
    plus a short method tag. Mirrors make_psf's model-tier selection so the injected PSF is
    built from the SAME model source the analytic kernel would use:

      ACS/WFC  -> focus-diverse ePSF per exposure (native F555W/F814W, focus/position matched),
                  falling back to STDPSF for ALL frames if any retrieval fails.
      WFC3/IR  -> exact-filter STDPSF grid (F160W has a real grid).
      WFPC2    -> the shared native WF3 F606W MAST-DB ePSF, falling back to STDPSF F555W.

    The drizzle itself does the exposure-average and North-up rotation, so each frame simply
    injects its own detector-frame model at its own lens position."""
    def grid_renderer(grid):
        def r(x, y, ext):
            stamp, iy0, ix0 = _render_from_grid(grid, x, y, _INJECT_FLUX)
            return stamp, iy0, ix0, ext
        return r

    def array_renderer(P, ov):
        def r(x, y, ext):
            stamp, iy0, ix0 = _render_from_array(P, ov, x, y, _INJECT_FLUX)
            return stamp, iy0, ix0, ext
        return r

    if inst_key == 'ACS/WFC':
        from acstools.focus_diverse_epsfs import psf_retriever, interp_epsf
        os.makedirs(psf_models._ACS_FD_CACHE, exist_ok=True)
        per_root, ok = {}, True
        for sci in sci_paths:
            root = os.path.basename(sci).split('_flc')[0]
            pos = psf_models._fd_detector_position(root, src_dir, catalogue_coord)
            if pos is None:
                print(f'    {root}: lens position/CD unavailable on FLC; FD unusable')
                ok = False
                break
            x, y, chip, _cd = pos
            try:
                grid = fits.getdata(psf_retriever(root, psf_models._ACS_FD_CACHE), ext=0)
            except Exception as exc:
                print(f'    {root}: focus-diverse retrieval failed ({exc})')
                ok = False
                break
            P = interp_epsf(grid, int(round(x)), int(round(y)), chip)
            if P is None or not np.all(np.isfinite(P)) or np.asarray(P).sum() <= 0:
                print(f'    {root}: interp_epsf returned no usable ePSF')
                ok = False
                break
            per_root[sci] = array_renderer(np.asarray(P, float),
                                           psf_models._ACS_FD_SUPERSAMPLE)
        if ok:
            return per_root, 'inject_acs_fdpsf'
        print('    falling back to STDPSF for all ACS frames')
        stdpsf_filt, _ = psf_models._resolve_filter(inst_key, filt)
        grid = psf_models._load_grid(inst_key, stdpsf_filt)
        return {s: grid_renderer(grid) for s in sci_paths}, 'inject_stdpsf'

    if inst_key == 'WFPC2' and psf_models._base_filter(filt).upper() == 'F606W':
        try:
            P, _samp = psf_models._wfpc2_f606w_db_build()
            r = array_renderer(np.asarray(P, float), psf_models._F606W_DB['build_oversample'])
            return {s: r for s in sci_paths}, 'inject_wfpc2_psfdb'
        except Exception as exc:
            print(f'    WFPC2 F606W PSF-DB unavailable ({exc}); falling back to STDPSF F555W')

    stdpsf_filt, sub = psf_models._resolve_filter(inst_key, filt)
    if sub:
        print(f'    NOTE: no STDPSF for {inst_key} {filt}; using nearest band {stdpsf_filt}')
    grid = psf_models._load_grid(inst_key, stdpsf_filt)
    return {s: grid_renderer(grid) for s in sci_paths}, 'inject_stdpsf'


def inject_frame(sci_path, renderer, catalogue_coord):
    """Zero every SCI extension of `sci_path` (in place) and add the rendered model PSF at the
    lens position on the chip that contains it. Multi-chip (ACS/WFC): the lens is on one chip;
    both chips are zeroed, the star added only to the containing chip. Returns True if a star
    was injected (lens on-frame), False otherwise."""
    injected = False
    with fits.open(sci_path, mode='update') as hdul:
        sci_hdus = [h for h in hdul if h.header.get('EXTNAME') == 'SCI']
        for h in sci_hdus:
            h.data = np.zeros_like(h.data, dtype=np.float32)
        for h in sci_hdus:
            ny, nx = h.data.shape
            try:
                x, y = WCS(h.header, hdul).world_to_pixel(catalogue_coord)
            except Exception:
                continue
            if not (0 <= x < nx and 0 <= y < ny):
                continue
            ext = h.header.get('EXTVER', 1)
            stamp, iy0, ix0, _ext = renderer(float(x), float(y), ext)
            S = stamp.shape[0]
            # clip the stamp to the frame
            dy0, dx0 = max(0, -iy0), max(0, -ix0)
            gy0, gx0 = max(0, iy0), max(0, ix0)
            gy1, gx1 = min(ny, iy0 + S), min(nx, ix0 + S)
            if gy1 <= gy0 or gx1 <= gx0:
                continue
            h.data[gy0:gy1, gx0:gx1] += stamp[dy0:dy0 + (gy1 - gy0),
                                              dx0:dx0 + (gx1 - gx0)].astype(np.float32)
            injected = True
            break  # lens is on a single chip
    return injected


# ── Drizzle and kernel extraction ────────────────────────────────────────────────
def run_injection_drizzle(inst_key, drizzle_input, catalogue_coord):
    """Run AstroDrizzle on the injected inputs in the current working dir, onto the science
    output grid. Returns the drizzled SCI path."""
    cfg = _DRIZ[inst_key]
    kw = dict(output='psf_inject', preserve=False, build=False, context=False,
              skysub=True, skymethod='localmin',
              driz_sep_wcs=True, driz_sep_scale=cfg['driz_sep_scale'],
              driz_sep_bits=cfg['driz_sep_bits'], driz_sep_fillval=-1,
              median=False, blot=False, driz_cr=False,
              resetbits=cfg['resetbits'],
              final_fillval=None, final_bits=cfg['final_bits'],
              final_wcs=True, final_scale=cfg['final_scale'],
              final_wht_type=cfg['wht_type'],
              final_rot=0.0, final_ra=float(catalogue_coord.ra.deg),
              final_dec=float(catalogue_coord.dec.deg),
              num_cores=1)
    if cfg['final_pixfrac'] is not None:
        kw['final_pixfrac'] = cfg['final_pixfrac']
    print('\n=== AstroDrizzle (injected artificial star) ===')
    astrodrizzle.AstroDrizzle(drizzle_input, **kw)
    for f in glob.glob('*ask.fits'):
        os.remove(f)
    out = sorted(glob.glob('psf_inject*_sci.fits'))
    if not out:
        raise RuntimeError('injection drizzle produced no psf_inject*_sci.fits')
    return out[0]


def extract_kernel(out_sci, catalogue_coord, size):
    """Cut a centred, unit-sum, odd `size`x`size` image-scale kernel from the drizzled star.

    The drizzled star is already at the science output pixel scale and North-up, so it IS the
    drizzle-broadened PSF -- we only need to centre it on the central pixel. Locate the star at
    the lens sky position, refine with a windowed centroid, then cubic-resample onto a grid
    centred on the centroid (a sub-pixel shift far smaller than the drizzle kernel). Returns
    (kernel, raw_stamp) where raw_stamp is the un-recentred cut for QA."""
    with fits.open(out_sci) as hdul:
        # AstroDrizzle fills uncovered pixels with NaN (final_fillval=None); zero them so
        # argmax / centroid / cubic sampling don't latch onto NaN outside the coverage.
        data = np.nan_to_num(np.asarray(hdul[0].data, float))
        wcs = WCS(hdul[0].header)
    px, py = wcs.world_to_pixel(catalogue_coord)
    ny, nx = data.shape
    ix, iy = int(round(float(px))), int(round(float(py)))

    # windowed centroid around the local peak near the lens position
    w = max(size, 21)
    y0, x0 = max(0, iy - w // 2), max(0, ix - w // 2)
    sub = data[y0:min(ny, y0 + w), x0:min(nx, x0 + w)]
    py2, px2 = np.unravel_index(np.argmax(sub), sub.shape)
    cw = 11
    yy0, xx0 = max(0, py2 - cw // 2), max(0, px2 - cw // 2)
    cy, cx = centroid_com(sub[yy0:min(sub.shape[0], yy0 + cw),
                              xx0:min(sub.shape[1], xx0 + cw)])
    cy, cx = y0 + yy0 + cy, x0 + xx0 + cx
    if not (np.isfinite(cx) and np.isfinite(cy)):
        cy, cx = float(iy), float(ix)

    half = size // 2
    raw = make_psf._crop_centered(data, size, cx, cy)

    offs = np.arange(size) - half
    XX, YY = np.meshgrid(offs, offs)
    kern = map_coordinates(data, [(cy + YY).ravel(), (cx + XX).ravel()], order=3,
                           mode='constant', cval=0.0).reshape(size, size)
    kern = np.clip(kern, 0.0, None)
    s = kern.sum()
    if s > 0:
        kern = kern / s
    return kern, raw


# ── Main ──────────────────────────────────────────────────────────────────────────
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--lens', default='J0252+0039')
    p.add_argument('--filt', default='f814W')
    p.add_argument('--sample', default=mast_target_names.DEFAULT_SAMPLE,
                   help=f'sample subdirectory (default {mast_target_names.DEFAULT_SAMPLE})')
    p.add_argument('--pass', dest='drizzle_pass', choices=['auto', 'cr', 'nocrrej'],
                   default='auto', help="science pass to match for the trimmed cutout prefix "
                                        "and star recentring (default auto: CR if present).")
    p.add_argument('--trim-threshold', dest='trim_threshold', type=float, default=1e-3,
                   help='amplitude fraction of peak at which the trimmed modelling kernel is '
                        'cut (matches make_psf; an amplitude criterion, not enclosed-energy).')
    p.add_argument('--keep-work', action='store_true',
                   help='keep the injection working dir (data/drizzle_files/.../psf_inject/) '
                        'for inspection instead of removing it.')
    a = p.parse_args()

    drizzled_dir = os.path.join(ws_path, 'data', 'drizzled', a.sample, a.lens, a.filt)
    src_dir = os.path.join(ws_path, 'data', 'drizzle_files', a.sample, a.lens, a.filt)
    psf_dir = os.path.join(ws_path, 'data', 'psf', a.sample, a.lens, a.filt)
    cutouts_dir = os.path.join(ws_path, 'data', 'cutouts', a.sample, a.lens, a.filt)
    json_path = os.path.join(ws_path, 'info', 'lens_psf_injected.json')

    # No-data outcome: matches the other scripts -- record null, exit 0.
    if not os.path.isdir(drizzled_dir) or not _has_products(drizzled_dir):
        print(f'=== NO DATA: {a.lens} {a.filt} (no drizzled products in {drizzled_dir})')
        info_json.update(json_path, a.sample, a.lens, a.filt, None)
        sys.exit(0)

    # Resolve the science pass (for the cutout prefix + recentring reference).
    if a.drizzle_pass == 'auto':
        has_cr = bool(glob.glob(os.path.join(drizzled_dir, '*_cr_*_sci.fits')))
        drizzle_pass = 'cr' if has_cr else 'nocrrej'
    else:
        drizzle_pass = a.drizzle_pass
    sci_file, _wht = find_products(drizzled_dir, drizzle_pass)
    with fits.open(sci_file) as hdul:
        sci_hdr = hdul[0].header
    inst_key = make_psf.instrument_key(sci_hdr)
    scale = proj_plane_pixel_scales(WCS(sci_hdr).celestial)[0] * 3600.0
    if inst_key not in _DRIZ:
        raise KeyError(f'no injection drizzle config for instrument {inst_key!r}')
    size = int(make_psf._INSTR[inst_key]['star_size'])
    print(f'{a.lens} {a.filt}  [{inst_key}, {drizzle_pass} pass, {scale:.4f}"/px, kernel {size}px]')

    if a.lens not in slacs_coords:
        raise KeyError(f'{a.lens} not in slacs_coords (info/slacs_coords.py)')
    ra, dec = slacs_coords[a.lens]
    catalogue_coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')

    _set_ref_env()
    work_dir = os.path.join(src_dir, 'psf_inject')
    if os.path.isdir(work_dir):
        shutil.rmtree(work_dir)

    cwd = os.getcwd()
    try:
        print('  staging science-prepared inputs ...')
        drizzle_input, sci_paths, _ivm = stage_inputs(inst_key, src_dir, work_dir)
        print(f'  {len(sci_paths)} frames')

        print('  building model-PSF renderers ...')
        renderers, method = build_renderers(inst_key, a.filt, a.sample, a.lens,
                                            sci_paths, src_dir, catalogue_coord)
        print(f'  model source: {method}')

        n_inj = 0
        for sci in sci_paths:
            if inject_frame(sci, renderers[sci], catalogue_coord):
                n_inj += 1
        print(f'  injected the model PSF into {n_inj}/{len(sci_paths)} frames')
        if n_inj == 0:
            raise RuntimeError('lens fell off every frame; nothing injected')

        os.chdir(work_dir)
        out_sci = run_injection_drizzle(inst_key, drizzle_input, catalogue_coord)
        kernel, raw = extract_kernel(out_sci, catalogue_coord, size)
    finally:
        os.chdir(cwd)

    # Post-process identically to make_psf's analytic kernel.
    kernel, pedestal_frac = make_psf.subtract_pedestal(kernel)
    if abs(pedestal_frac) >= 1e-4:
        print(f'  pedestal removed: {pedestal_frac:.2e} of peak')
    ped_frac, scat_frac = make_psf.wing_stats(kernel)
    trimmed, trim_size = make_psf.trim_kernel_to_amplitude(kernel, a.trim_threshold)
    fwhm_pix = make_psf.measure_fwhm(kernel)
    if fwhm_pix is not None:
        print(f'  kernel FWHM: {fwhm_pix:.2f} px = {fwhm_pix * scale:.3f}"')
    print(f'  full kernel {kernel.shape[0]}px -> trimmed {trim_size}px '
          f'(amplitude<{a.trim_threshold:g} of peak)')

    # ── Write products ──────────────────────────────────────────────────────────
    os.makedirs(psf_dir, exist_ok=True)
    khdr = fits.Header()
    khdr['PSFMETH'] = (method, 'injection-drizzled model PSF (Anderson 2016)')
    khdr['PSFOVSMP'] = (1, 'injected kernel is at image scale (no oversampling)')
    khdr['PSFPXSCL'] = (round(scale, 5), 'kernel pixel scale (arcsec)')
    khdr['PSFFWHM'] = (round(fwhm_pix, 4) if fwhm_pix is not None else 0.0,
                       'fitted kernel FWHM (pixels)')
    khdr['PSFLENS'] = (a.lens, 'lens')
    khdr['PSFFILT'] = (a.filt, 'filter')
    khdr['PSFKIND'] = ('full', 'full kernel (trimmed copy in data/cutouts/)')
    khdr['PSFPED'] = (round(pedestal_frac, 6), 'ePSF-wing pedestal removed (fraction of peak)')
    khdr['PSFINJ'] = (True, 'built by artificial-star injection + re-drizzle')
    make_psf.write_fits(kernel, khdr, os.path.join(psf_dir, 'psf_kernel_injected.fits'))

    # QA panel: reuse make_psf.plot_psf; the "star montage" shows the raw drizzled star.
    raw_norm = raw / raw.max() if raw.max() > 0 else raw
    make_psf.plot_psf(kernel, None, [raw_norm], scale, fwhm_pix,
                      os.path.join(psf_dir, 'psf_injected.png'),
                      title=f'{a.lens}  {a.filt}  [{method}]  (injected)')

    os.makedirs(cutouts_dir, exist_ok=True)
    prefix = 'cutout_cr' if drizzle_pass == 'cr' else 'cutout'
    thdr = khdr.copy()
    thdr['PSFKIND'] = ('trimmed', 'amplitude-trimmed modelling kernel')
    thdr['PSFTRIM'] = (a.trim_threshold, 'azimuthal PSF < this x peak sets the radius')
    thdr['PSFPASS'] = (drizzle_pass, 'drizzle pass the PSF matches')
    make_psf.write_fits(trimmed, thdr,
                        os.path.join(cutouts_dir, f'{prefix}_psf_injected.fits'))

    info_json.update(json_path, a.sample, a.lens, a.filt, {
        'method': method,
        'fwhm_pix': round(fwhm_pix, 4) if fwhm_pix is not None else None,
        'kernel_size': int(kernel.shape[0]),
        'cutout_kernel_size': int(trim_size),
        'trim_threshold': a.trim_threshold,
        'pedestal_frac': round(pedestal_frac, 6),
        'wing_scatter': round(scat_frac, 6),
        'n_frames': n_inj,
    })

    if not a.keep_work:
        shutil.rmtree(work_dir, ignore_errors=True)
    print('  done')


if __name__ == '__main__':
    main()
