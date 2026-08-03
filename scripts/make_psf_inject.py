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

For the model-PSF tier, this is now the CANONICAL downstream product -- empirical builds
already carry the broadening (cut from the drizzled mosaic itself), but a model-tier build
does not, and the injected kernel is strictly the more correct one when it is available. So
when the product's primary build (per info/lens_psf.json) is on the model tier, this PROMOTES
the injected kernel: the existing pre-broadening analytic-model files are moved aside to
*_analytic (comparison only), and the injected build takes over the canonical names:

  data/psf/<sample>/<lens>/<filt>/psf_kernel.fits             canonical (now the injected build)
  data/psf/<sample>/<lens>/<filt>/psf_kernel_analytic.fits    the superseded analytic model
  data/psf/<sample>/<lens>/<filt>/psf.png / psf_analytic.png  QA panels, same split
  data/cutouts/<sample>/<lens>/<filt>/cutout_[cr_]psf.fits            canonical (injected)
  data/cutouts/<sample>/<lens>/<filt>/cutout_[cr_]psf_analytic.fits   superseded analytic
  info/lens_psf.json                                          method becomes inject_*
  info/lens_psf_injected.json                                 {sample:{lens:{filt:{...}}}}

make_psf.py auto-chains into this (run_injection(..., promote=True)) right after building a
model-tier product, so a single `make_psf.py --lens X --filt Y` run already leaves the
promoted, drizzle-broadened kernel as the canonical file -- this script rarely needs to be
run by hand for that case. It IS still run standalone for: (1) the one lens missing its
injected build, (2) re-promoting after an injection-code change without rebuilding the
analytic model, and (3) `--all` mode, which also runs on EMPIRICAL primaries purely as a
validation comparison (injected-model FWHM should approach the empirical truth) -- there,
promotion does not apply and the parallel `*_injected`-suffixed names are used instead, as
before.

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
    uv run python scripts/make_psf_inject.py --lens J0252+0039 --filt f814W
    uv run python scripts/make_psf_inject.py --lens J0822+2652 --filt f606W_v2
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

import mast_target_names
import psf_models
# Reuse make_psf's kernel post-processing + QA + JSON helpers verbatim so the injected
# product is treated identically to the analytic one (importing make_psf is safe: its work
# is under __main__). instrument_key / representative field sizes also come from there.
import make_psf
from make_cutouts import find_products, _has_products, catalogue_coord_for
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
    'WFC3/UVIS': dict(input='flc', suffix='drc', det_scale=0.0396,
                      driz_sep_scale=0.0396, driz_sep_bits='256,64,16',
                      final_scale=0.0396, final_pixfrac=0.7, final_bits='256,64,16',
                      wht_type='ERR', resetbits=0),
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


# ── Promotion: make the injected (drizzle-broadened) kernel the canonical product ──
def _is_model_tier(sample, lens, filt):
    """True if the CURRENT primary build recorded in info/lens_psf.json for this product is
    on the model tier -- either not yet promoted ('model...') or already promoted
    ('inject...') by a previous run of this function. False for 'empirical' or no record."""
    psf_json = os.path.join(ws_path, 'info', 'lens_psf.json')
    rec = info_json.load(psf_json).get(sample, {}).get(lens, {}).get(filt)
    method = str((rec or {}).get('method', ''))
    return method.startswith('model') or method.startswith('inject')


def _promote(psf_dir, cutouts_dir, prefix):
    """Move the current canonical model-tier files aside to *_analytic (once -- a file
    already carrying PSFINJ=True or PSFBROAD=True is itself a previously-promoted injected
    or analytically-broadened build, not the sharp analytic model, so it is left alone
    rather than clobbering the real analytic backup with a non-analytic file).

    Also moves aside any analytic-tier PSF error map (psf_kernel_err.fits /
    cutout_[cr_]psf_err.fits -- see make_psf.py's ACS focus-diverse / WFPC2 MAST-DB
    ensembles) to *_analytic_err, rather than leaving it under the canonical name where it
    would describe the just-superseded analytic kernel, not the injected one that replaces
    it. The injected kernel has no uncertainty estimate of its own yet, so after promotion
    there is deliberately no canonical psf_kernel_err.fits -- not a wrong one.
    """
    canonical_kernel = os.path.join(psf_dir, 'psf_kernel.fits')
    analytic_kernel = os.path.join(psf_dir, 'psf_kernel_analytic.fits')
    if os.path.isfile(canonical_kernel):
        hdr = fits.getheader(canonical_kernel)
        already_promoted = bool(hdr.get('PSFINJ', False)) or bool(hdr.get('PSFBROAD', False))
        if not already_promoted:
            shutil.move(canonical_kernel, analytic_kernel)
            png = os.path.join(psf_dir, 'psf.png')
            if os.path.isfile(png):
                shutil.move(png, os.path.join(psf_dir, 'psf_analytic.png'))
            cut = os.path.join(cutouts_dir, f'{prefix}_psf.fits')
            if os.path.isfile(cut):
                shutil.move(cut, os.path.join(cutouts_dir, f'{prefix}_psf_analytic.fits'))
            kernel_err = os.path.join(psf_dir, 'psf_kernel_err.fits')
            if os.path.isfile(kernel_err):
                shutil.move(kernel_err, os.path.join(psf_dir, 'psf_kernel_analytic_err.fits'))
            cut_err = os.path.join(cutouts_dir, f'{prefix}_psf_err.fits')
            if os.path.isfile(cut_err):
                shutil.move(cut_err,
                           os.path.join(cutouts_dir, f'{prefix}_psf_analytic_err.fits'))


def analytic_broadened_fallback(psf_dir, cutouts_dir, prefix, sci_hdr, method,
                                drizzle_pass, trim_threshold=1e-3):
    """Cheap stand-in for a FAILED run_injection(): approximate AstroDrizzle's resampling
    broadening analytically (psf_models.analytic_drop_broaden -- a pixfrac-wide box
    convolution, sized from the science drizzle's own D001PIXF/D001ISCL/D001SCAL and the
    number of contributing dithered frames) instead of the real re-drizzle injection. This
    is make_psf.py's fallback when
    make_psf_inject.run_injection() raises (e.g. data/drizzle_files/ was cleared): rather
    than silently leaving the sharp, un-broadened analytic model canonical, promote the
    analytically-broadened one -- an approximation, but strictly closer to the true
    drizzled PSF than the un-broadened model, same logic as why the real injected kernel is
    promoted when it succeeds. Idempotent like _promote(): a canonical file already carrying
    PSFBROAD=True (from a previous call) is treated as already-promoted and re-broadened in
    place from the untouched *_analytic backup, not re-derived from itself.

    Distinguished from a real injected build by PSFINJ=False / PSFBROAD=True and a
    'broadened_<method>' PSFMETH/info-json method, so nothing downstream mistakes it for
    the rigorous re-drizzled kernel. Also prints a sanity-check FWHM comparison against the
    un-broadened analytic model -- the real injected FWHM, on lenses where both exist, is
    always the one to trust over this approximation.

    Returns the record dict written to info/lens_psf.json (mirroring run_injection's
    promoted-record shape, minus the injection-only fields), or raises if there is no
    canonical kernel to broaden.
    """
    canonical_kernel = os.path.join(psf_dir, 'psf_kernel.fits')
    if not os.path.isfile(canonical_kernel):
        raise FileNotFoundError(f'{canonical_kernel} missing; nothing to broaden')
    if fits.getheader(canonical_kernel).get('PSFINJ', False):
        raise RuntimeError(f'{canonical_kernel} is already a real injected build; '
                           'refusing to overwrite it with an analytic approximation')

    pixfrac = float(sci_hdr['D001PIXF'])
    native_scale = float(sci_hdr['D001ISCL'])
    out_scale = float(sci_hdr['D001SCAL'])
    # NDRIZIM is the drizzled IMAGE count, not exposure count: ACS/WFC and WFC3/UVIS FLCs
    # are 2-chip MEFs (NDRIZIM = 2x exposures; see CLAUDE.md's tracking-JSON section), and
    # both chips of one exposure share the same dither phase, so they must not be counted
    # as two independent phase samples. WFC3/IR and WFPC2 (single chip) need no correction.
    inst_key = make_psf.instrument_key(sci_hdr)
    n_chips = 2 if inst_key in ('ACS/WFC', 'WFC3/UVIS') else 1
    n_frames = max(1, int(sci_hdr.get('NDRIZIM', n_chips)) // n_chips)

    _promote(psf_dir, cutouts_dir, prefix)   # moves the sharp analytic files aside (no-op
                                              # if already broadened by a prior call)

    analytic_kernel_path = os.path.join(psf_dir, 'psf_kernel_analytic.fits')
    with fits.open(analytic_kernel_path) as hdul:
        analytic_kernel = np.asarray(hdul[0].data, dtype=float)
        analytic_hdr = hdul[0].header.copy()

    broadened, box_width = psf_models.analytic_drop_broaden(
        analytic_kernel, pixfrac, native_scale, out_scale, n_frames=n_frames)
    broadened, pedestal_frac = make_psf.subtract_pedestal(broadened)
    fwhm_before = make_psf.measure_fwhm(analytic_kernel)
    fwhm_after = make_psf.measure_fwhm(broadened)
    print(f'  analytic drizzle-broadening fallback: box {box_width:.2f} output-px '
         f'(pixfrac={pixfrac:g}, native={native_scale:.4f}"/px, out={out_scale:.4f}"/px, '
         f'{n_frames} frames)')
    if fwhm_before is not None and fwhm_after is not None:
        print(f'  sanity check: FWHM {fwhm_before:.2f} -> {fwhm_after:.2f} px '
             f'(analytic model -> analytically broadened; trust a real injected FWHM '
             f'over this where one exists)')

    trimmed, trim_size = make_psf.trim_kernel_to_amplitude(broadened, trim_threshold)

    khdr = analytic_hdr
    khdr['PSFMETH'] = (f'broadened_{method}', 'analytic drop-box broadened model PSF')
    khdr['PSFINJ'] = (False, 'NOT built by artificial-star injection + re-drizzle')
    khdr['PSFBROAD'] = (True, 'cheap analytic pixfrac drop-box broadening (fallback)')
    khdr['PSFBOXPX'] = (round(box_width, 4), 'drop-box width applied, output pixels')
    khdr['PSFPED'] = (round(pedestal_frac, 6), 'pedestal removed after broadening (frac of peak)')
    make_psf.write_fits(broadened, khdr, canonical_kernel)

    raw_norm = analytic_kernel / analytic_kernel.max() if analytic_kernel.max() > 0 else analytic_kernel
    make_psf.plot_psf(broadened, None, [raw_norm], out_scale, fwhm_after,
                      os.path.join(psf_dir, 'psf.png'),
                      title=f'{khdr.get("PSFLENS", "")}  {khdr.get("PSFFILT", "")}  '
                            f'[broadened_{method}]  (analytic fallback)')

    thdr = khdr.copy()
    thdr['PSFKIND'] = ('trimmed', 'amplitude-trimmed modelling kernel')
    thdr['PSFTRIM'] = (trim_threshold, 'azimuthal PSF < this x peak sets the radius')
    thdr['PSFPASS'] = (drizzle_pass, 'drizzle pass the PSF matches')
    make_psf.write_fits(trimmed, thdr, os.path.join(cutouts_dir, f'{prefix}_psf.fits'))

    return {
        'method': f'broadened_{method}',
        'fwhm_pix': round(fwhm_after, 4) if fwhm_after is not None else None,
        'kernel_size': int(broadened.shape[0]),
        'cutout_kernel_size': int(trim_size),
        'trim_threshold': trim_threshold,
        'pedestal_frac': round(pedestal_frac, 6),
        'box_width_px': round(box_width, 4),
    }


# ── Callable entry point (also used by make_psf.py to auto-chain model-tier builds) ─
def run_injection(lens, filt, sample=None, drizzle_pass='auto', trim_threshold=1e-3,
                  keep_work=False, promote=None):
    """Build the drizzle-broadened injected PSF for (sample, lens, filt) and, when the
    product's primary build is on the model tier, PROMOTE it to the canonical
    psf_kernel.fits / cutout_[cr_]psf.fits (the pre-broadening analytic model is moved
    aside to *_analytic for comparison). Returns the record dict written to
    info/lens_psf_injected.json, or None on a no-data outcome.

    `promote`: None (default) auto-detects from info/lens_psf.json (True for a model-tier
    primary, False -- i.e. a parallel, non-canonical '*_injected' product -- for an
    empirical primary, matching the --all validation mode). Pass True/False to override,
    e.g. make_psf.py knows definitively it just built a model-tier product."""
    sample = sample or mast_target_names.DEFAULT_SAMPLE
    drizzled_dir = os.path.join(ws_path, 'data', 'drizzled', sample, lens, filt)
    src_dir = os.path.join(ws_path, 'data', 'drizzle_files', sample, lens, filt)
    psf_dir = os.path.join(ws_path, 'data', 'psf', sample, lens, filt)
    cutouts_dir = os.path.join(ws_path, 'data', 'cutouts', sample, lens, filt)
    json_path = os.path.join(ws_path, 'info', 'lens_psf_injected.json')

    # No-data outcome: matches the other scripts -- record null, exit 0.
    if not os.path.isdir(drizzled_dir) or not _has_products(drizzled_dir):
        print(f'=== NO DATA: {lens} {filt} (no drizzled products in {drizzled_dir})')
        info_json.update(json_path, sample, lens, filt, None)
        return None

    # Resolve the science pass (for the cutout prefix + recentring reference).
    if drizzle_pass == 'auto':
        has_cr = bool(glob.glob(os.path.join(drizzled_dir, '*_cr_*_sci.fits')))
        drizzle_pass = 'cr' if has_cr else 'nocrrej'
    sci_file, _wht = find_products(drizzled_dir, drizzle_pass)
    with fits.open(sci_file) as hdul:
        sci_hdr = hdul[0].header
    inst_key = make_psf.instrument_key(sci_hdr)
    scale = proj_plane_pixel_scales(WCS(sci_hdr).celestial)[0] * 3600.0
    if inst_key not in _DRIZ:
        raise KeyError(f'no injection drizzle config for instrument {inst_key!r}')
    size = int(make_psf._INSTR[inst_key]['star_size'])
    print(f'{lens} {filt}  [{inst_key}, {drizzle_pass} pass, {scale:.4f}"/px, kernel {size}px]')

    ra, dec = catalogue_coord_for(lens)  # slacs_coords or gallery_coords
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
        renderers, method = build_renderers(inst_key, filt, sample, lens,
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
    trimmed, trim_size = make_psf.trim_kernel_to_amplitude(kernel, trim_threshold)
    fwhm_pix = make_psf.measure_fwhm(kernel)
    if fwhm_pix is not None:
        print(f'  kernel FWHM: {fwhm_pix:.2f} px = {fwhm_pix * scale:.3f}"')
    print(f'  full kernel {kernel.shape[0]}px -> trimmed {trim_size}px '
          f'(amplitude<{trim_threshold:g} of peak)')

    if promote is None:
        promote = _is_model_tier(sample, lens, filt)

    os.makedirs(psf_dir, exist_ok=True)
    os.makedirs(cutouts_dir, exist_ok=True)
    prefix = 'cutout_cr' if drizzle_pass == 'cr' else 'cutout'

    khdr = fits.Header()
    khdr['PSFMETH'] = (method, 'injection-drizzled model PSF (Anderson 2016)')
    khdr['PSFOVSMP'] = (1, 'injected kernel is at image scale (no oversampling)')
    khdr['PSFPXSCL'] = (round(scale, 5), 'kernel pixel scale (arcsec)')
    khdr['PSFFWHM'] = (round(fwhm_pix, 4) if fwhm_pix is not None else 0.0,
                       'fitted kernel FWHM (pixels)')
    khdr['PSFLENS'] = (lens, 'lens')
    khdr['PSFFILT'] = (filt, 'filter')
    khdr['PSFKIND'] = ('full', 'full kernel (trimmed copy in data/cutouts/)')
    khdr['PSFPED'] = (round(pedestal_frac, 6), 'ePSF-wing pedestal removed (fraction of peak)')
    khdr['PSFINJ'] = (True, 'built by artificial-star injection + re-drizzle')

    if promote:
        # Move the existing analytic-model files aside, then write the injected build
        # under the CANONICAL names -- this is what downstream (and every other script)
        # reads. No separate "_injected"-suffixed files: promotion means there is exactly
        # one canonical kernel again, just now the drizzle-broadened one.
        _promote(psf_dir, cutouts_dir, prefix)
        kernel_path = os.path.join(psf_dir, 'psf_kernel.fits')
        png_path = os.path.join(psf_dir, 'psf.png')
        cutout_path = os.path.join(cutouts_dir, f'{prefix}_psf.fits')
        title = f'{lens}  {filt}  [{method}]'
    else:
        # Empirical primary (--all validation mode): keep this as a parallel, clearly-named
        # comparison product; nothing canonical changes.
        kernel_path = os.path.join(psf_dir, 'psf_kernel_injected.fits')
        png_path = os.path.join(psf_dir, 'psf_injected.png')
        cutout_path = os.path.join(cutouts_dir, f'{prefix}_psf_injected.fits')
        title = f'{lens}  {filt}  [{method}]  (injected)'

    make_psf.write_fits(kernel, khdr, kernel_path)

    # QA panel: reuse make_psf.plot_psf; the "star montage" shows the raw drizzled star.
    raw_norm = raw / raw.max() if raw.max() > 0 else raw
    make_psf.plot_psf(kernel, None, [raw_norm], scale, fwhm_pix, png_path, title=title)

    thdr = khdr.copy()
    thdr['PSFKIND'] = ('trimmed', 'amplitude-trimmed modelling kernel')
    thdr['PSFTRIM'] = (trim_threshold, 'azimuthal PSF < this x peak sets the radius')
    thdr['PSFPASS'] = (drizzle_pass, 'drizzle pass the PSF matches')
    make_psf.write_fits(trimmed, thdr, cutout_path)

    record = {
        'method': method,
        'fwhm_pix': round(fwhm_pix, 4) if fwhm_pix is not None else None,
        'kernel_size': int(kernel.shape[0]),
        'cutout_kernel_size': int(trim_size),
        'trim_threshold': trim_threshold,
        'pedestal_frac': round(pedestal_frac, 6),
        'wing_scatter': round(scat_frac, 6),
        'n_frames': n_inj,
    }
    info_json.update(json_path, sample, lens, filt, record)

    if promote:
        # The canonical lens_psf.json entry must describe what's now actually in
        # cutout_[cr_]psf.fits: the injected build, not the superseded analytic one.
        psf_json = os.path.join(ws_path, 'info', 'lens_psf.json')
        info_json.update(psf_json, sample, lens, filt, dict(record))
        print(f'  promoted: {method} is now the canonical PSF '
              f'(analytic model kept as *_analytic)')

    if not keep_work:
        shutil.rmtree(work_dir, ignore_errors=True)
    print('  done')
    return record


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
    run_injection(a.lens, a.filt, a.sample, drizzle_pass=a.drizzle_pass,
                 trim_threshold=a.trim_threshold, keep_work=a.keep_work)


if __name__ == '__main__':
    main()
