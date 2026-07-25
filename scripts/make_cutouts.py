#!/usr/bin/env python
"""
Cut science + noise-map postage stamps out of the drizzled mosaics.

Reads the no-CR-rejection drizzle products from data/drizzled/<sample>/<lens>/<filt>/,
centres a cutout on the lens coordinate from info/slacs_coords.py, then recentres
onto the brightest pixel found within the central 100x100 pixels of that first
cutout and re-cuts so the stamp is exactly centred on the lens galaxy.

The noise map is derived from the drizzle weight map as 1/sqrt(weight), with
zero-weight pixels mapped to a large value so they are excluded from any fit.

Writes cutout_sci.fits, cutout_noise.fits and a 3-panel cutout.png
(signal / noise / signal-to-noise) for visual inspection.

Usage:
    conda run -n stenv python scripts/make_cutouts.py --lens J0008-0004 --filt f814W
"""

import argparse
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.nddata.utils import NoOverlapError
from astropy.wcs import WCS
from astropy.visualization import AsinhStretch, ImageNormalize, PercentileInterval
from astropy.wcs.utils import proj_plane_pixel_scales
from scipy.ndimage import median_filter
import astropy.units as u

ws_path = '/Users/samlange/Code/data_reduction'
sys.path.insert(0, os.path.join(ws_path, 'info'))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from slacs_coords import slacs_coords
import mast_target_names


# Detectors whose calibrated ERR array is already in ELECTRONS/S rather than counts.
# This single fact decides whether 1/sqrt(WHT) is a calibrated sigma map or is wrong
# by a factor of the per-frame exposure time -- see weight_to_sigma_scale().
_ERR_IN_RATE_UNITS = {('WFC3', 'IR')}


def weight_to_sigma_scale(sci_hdr):
    """
    Return (K, note) such that the calibrated per-pixel sigma is K / sqrt(WHT).

    `1/sqrt(WHT)` is only a sigma map when WHT is a true inverse-variance map, and
    with `final_wht_type=ERR` that depends on the *units of the input ERR array*.
    DrizzlePac computes the per-frame weight as

        weight = (EXPTIME / ERR)**2                  # imageObject.buildERRmask

    For ACS FLC (SCI and ERR both in ELECTRONS) `EXPTIME/ERR` is exactly
    `1/sigma_rate`, so WHT is a genuine inverse variance of the ELECTRONS/S output
    and K = 1. For WFC3/IR FLT, SCI and ERR are **already ELECTRONS/S**, so the same
    expression evaluates to `EXPTIME/sigma_rate`: every weight is inflated by
    EXPTIME**2 and 1/sqrt(WHT) comes out a factor EXPTIME too small. K = per-frame
    EXPTIME undoes exactly that.

    Verified on J0841+3824 with a blank-sky block-sum test (scatter of integrated
    flux in NxN blank blocks against what the noise map predicts, which converges to
    the truth once the block exceeds the drizzle correlation length):

        F814W (ACS, K=1)          ratio 1.04 at 0.24", 1.24 at 1.44"  -> correct
        F160W (WFC3/IR, uncorr.)  ratio 477 at 0.24",  725 at 1.44"   -> ~700x low

    and 700 / 599.23 = 1.17, i.e. once K = EXPTIME is applied the only thing left is
    the same drizzle correlated-noise factor ACS shows independently (1.24). The
    `D001WTSC = 1/scale**4` term does *not* enter: it cancels against the finer
    output grid, which is why ACS (native scale, WTSC 1.0) and WFC3/IR (WTSC 20.9)
    share one formula.

    WFPC2 F606W also takes K = 1, and that is now correct rather than a placeholder.
    It reaches the same place by a different route: DrizzlePac's IVM branch leaves the
    supplied IVM unscaled but sets `wt_scl = exptime**2/scale**4` (against the ERR
    branch's `1/scale**4`), so the exptime**2 the ERR mask carries internally is
    supplied by wt_scl instead and the two paths agree. With the IVM built as
    1/(SCI/gain + floor**2) in DN**-2, 1/sqrt(WHT) is a sigma in DN/s, matching the
    DN/s output. Measured on J0330-0020: block ratio 0.90 / 1.07 / 1.09.

    Products predating that fix are NOT calibrated and cannot be rescued by any K --
    the error is in the noise model, not the units. The caller warns on IVMMODEL.
    """
    inst = (sci_hdr.get('INSTRUME', '').strip(), sci_hdr.get('DETECTOR', '').strip())
    if inst not in _ERR_IN_RATE_UNITS:
        return 1.0, f'{"/".join(inst)}: ERR in counts, 1/sqrt(WHT) already calibrated'

    # Per-frame exposure times, not EXPTIME/NDRIZIM: the correction is only a single
    # constant when the frames are equal, so read them and refuse to guess if not.
    dexp = sorted({round(v, 3) for k, v in sci_hdr.items()
                   if k.startswith('D') and k.endswith('DEXP')})
    if not dexp:
        raise KeyError('no D00nDEXP keywords in the drizzle header; cannot determine '
                       'the per-frame exposure time needed to calibrate the noise map')
    if len(dexp) > 1:
        raise ValueError(
            f'unequal per-frame exposure times {dexp} for {"/".join(inst)}. The ERR '
            'units correction is a single constant only for equal-length frames; a '
            'mixed-exposure stack needs a per-pixel treatment that is not implemented.')

    return dexp[0], (f'{"/".join(inst)}: ERR in ELECTRONS/S, scaling by per-frame '
                     f'EXPTIME={dexp[0]:.3f}s')


def noise_map_via_weight_map_from(weight_map, scale=1.0):
    """
    Setup the noise-map from a weight map, which is a form of noise-map that comes via HST
    image-reduction and the software package MultiDrizzle.

    The noise in each pixel is computed as:

    sigma = scale / sqrt(weight_map).

    `scale` is the calibration constant from weight_to_sigma_scale() (times any
    correlated-noise inflation); it is 1.0 for detectors whose ERR array is in counts,
    which is the classic MultiDrizzle recipe.

    The weight map may contain zeros, in which case the variances are converted to large
    values to omit them from the analysis.

    Parameters
    ----------
    weight_map
        The weight-value of each pixel which is converted to a variance.
    scale
        Multiplicative constant converting 1/sqrt(weight) into a calibrated sigma.
    """
    np.seterr(divide="ignore")
    noise_map = scale / weight_map ** 0.5
    noise_map[noise_map > 1.0e8] = 1.0e8
    return noise_map


def find_products(drizzled_dir, drizzle_pass='nocrrej'):
    """
    Locate the sci/wht pair for a given drizzle pass in a drizzled output directory.

    `drizzle_pass` is 'nocrrej' (default) or 'cr'. The two patterns are disjoint:
    '*_cr_*' does not match '..._nocrrej_drc_...' since that has no '_cr_' substring.

    The filename prefix and drizzle suffix vary by instrument (_drw_ for WFPC2 FLT,
    _drc_ for ACS FLC, _drz_ for WFC3/IR and NICMOS), so glob rather than hardcode.
    """
    sci = glob.glob(os.path.join(drizzled_dir, f'*_{drizzle_pass}_*_sci.fits'))
    wht = glob.glob(os.path.join(drizzled_dir, f'*_{drizzle_pass}_*_wht.fits'))

    if len(sci) != 1 or len(wht) != 1:
        raise FileNotFoundError(
            f"Expected exactly one {drizzle_pass} sci and wht file in {drizzled_dir}, "
            f"found sci={sci}, wht={wht}"
        )
    return sci[0], wht[0]


def brightest_pixel_position(cutout, box=100, median_size=5):
    """
    Return the (x, y) pixel position of the brightest pixel within the central
    `box` x `box` pixel region of a Cutout2D, in the cutout's own pixel frame.

    Restricting to the central region stops a bright neighbour or a field star
    near the stamp edge from stealing the recentre.

    These are the no-CR-rejection mosaics, so cosmic rays are present by construction
    and a raw brightest-pixel search will happily lock onto a CR spike instead of the
    galaxy (J0252+0039: raw peak is a CR 48 pix off-centre, galaxy at 9 pix).

    Peak-finding therefore runs on a median-filtered copy. A median is the right tool
    here and Gaussian smoothing is not: smoothing conserves flux, so a compact 30-count
    CR still outranks a 1.3-count galaxy peak at any sigma (J0008-0004 fails this way).
    A rank filter removes any feature narrower than its window outright. Only the search
    is filtered — the returned position indexes the unmodified data.
    """
    data = np.nan_to_num(cutout.data)
    ny, nx = data.shape

    half = box // 2
    ylo, yhi = max(0, ny // 2 - half), min(ny, ny // 2 + half)
    xlo, xhi = max(0, nx // 2 - half), min(nx, nx // 2 + half)

    sub = data[ylo:yhi, xlo:xhi]
    search = median_filter(sub, median_size) if median_size > 1 else sub
    iy, ix = np.unravel_index(np.argmax(search), search.shape)

    # A peak pinned to the edge of the search box usually means the true maximum lies
    # outside it — i.e. we locked onto a neighbour rather than the lens. Worth flagging
    # rather than silently recentring onto the wrong object.
    edge = min(ix, iy, sub.shape[1] - 1 - ix, sub.shape[0] - 1 - iy)
    if edge <= 2:
        print(f"  WARNING: peak sits {edge} pix from the edge of the {box}x{box} search box "
              f"- likely a neighbouring source, not the lens. Check this cutout.")

    return xlo + ix, ylo + iy


def _has_products(drizzled_dir):
    return bool(glob.glob(os.path.join(drizzled_dir, '*_sci.fits')))


def find_peak_coord(drizzled_dir, catalogue_coord, size, box, median_size):
    """Deflector brightest-pixel sky position from a band's mosaic.

    Uses the CR-rejected pass when present (no cosmic rays to lock onto), else the
    no-CR pass. Returns (SkyCoord, path_used). May raise NoOverlapError if the
    catalogue position is off the mosaic.
    """
    try:
        cen_sci, _ = find_products(drizzled_dir, 'cr')
    except FileNotFoundError:
        cen_sci, _ = find_products(drizzled_dir, 'nocrrej')
    with fits.open(cen_sci) as hdul:
        cut = make_cutout(hdul[0].data, WCS(hdul[0].header), catalogue_coord, size)
    x, y = brightest_pixel_position(cut, box=box, median_size=median_size)
    return cut.wcs.pixel_to_world(x, y), cen_sci


def make_cutout(data, wcs, position, size_arcsec):
    """Cut a square stamp of `size_arcsec` on a side, centred on `position`."""
    size = u.Quantity((size_arcsec, size_arcsec), u.arcsec)
    return Cutout2D(data, position=position, size=size, wcs=wcs)


def write_cutout(data, cutout_wcs, header, output_path):
    """Write cutout data to a single-extension FITS file with an updated WCS."""
    output_hdr = header.copy()
    output_hdr.update(cutout_wcs.to_header())

    fits.PrimaryHDU(data=data, header=output_hdr).writeto(output_path, overwrite=True)
    print(f"  wrote {output_path}  shape={data.shape}")


def plot_cutouts(sci_data, noise_data, output_path, title):
    """
    Write a 3-panel PNG of the cutout: signal, noise, and signal-to-noise.

    Zero-weight pixels carry the 1e8 sentinel from noise_map_via_weight_map_from; they
    are masked out for display so they don't flatten the colour scale of the noise panel.

    Panels use origin='upper' so the orientation matches PyAutoLens, which plots row 0
    at the top - with origin='lower' the same cutout comes out mirrored in y relative to
    the autolens plots and it is easy to misidentify which neighbour is which.

    Colour/stretch (inferno + asinh, astropy.visualization) matches scripts/make_mosaics.py:
    asinh is well-defined through zero and for negative values, so it displays the
    negative background-noise pixels directly instead of needing them masked to NaN
    (as a log stretch would), while still showing faint outskirts and bright cores
    together.
    """
    noise_display = np.where(noise_data >= 1.0e8, np.nan, noise_data)
    with np.errstate(divide='ignore', invalid='ignore'):
        snr = sci_data / noise_data

    panels = [
        ('signal',        sci_data),
        ('noise',         noise_display),
        ('signal / noise', snr),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))
    for ax, (label, data) in zip(axes, panels):
        finite = data[np.isfinite(data)]
        vmin, vmax = PercentileInterval(99.0).get_limits(finite)
        norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=AsinhStretch(0.1))
        im = ax.imshow(data, norm=norm, origin='upper', cmap='inferno')
        ax.set_title(label)
        ax.set_xlabel('pixels')
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(output_path, dpi=120, bbox_inches='tight')
    plt.close(fig)
    print(f"  wrote {output_path}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--lens',   default='J0008-0004')
    p.add_argument('--filt',   default='f814W')
    p.add_argument('--sample', default=mast_target_names.DEFAULT_SAMPLE,
                   help='sample subdirectory of data/drizzled/ and data/cutouts/. '
                        'Defined in info/lens_samples.json '
                        f'(default {mast_target_names.DEFAULT_SAMPLE})')
    p.add_argument('--size',   type=float, default=20.0,
                   help='cutout size in arcsec (square), default 20.0')
    p.add_argument('--box',    type=int, default=100,
                   help='central box size in pixels searched for the brightest pixel')
    p.add_argument('--median-size', type=int, default=5,
                   help='median filter window (pix) applied before peak-finding, to stop '
                        'cosmic rays in the no-CR data being mistaken for the galaxy; '
                        '1 disables')
    p.add_argument('--pass', dest='drizzle_pass',
                   choices=['auto', 'cr', 'nocrrej'], default='auto',
                   help="which drizzle pass to cut from. 'auto' (default) prefers the "
                        "CR-rejected pass and falls back to no-CR when no CR product "
                        "exists (e.g. WFC3/IR F160W, which has no CR pass). 'cr' and "
                        "'nocrrej' force one. The CR pass is the science default: ACS "
                        # '%%' not '%': argparse runs the help through %-interpolation,
                        # so a literal percent sign here crashes --help entirely.
                        "uses LACosmic masking that keeps ~99%% of the deflector core "
                        "while removing cosmic rays, so the old no-CR default put "
                        "CR-riddled stamps in front of the user.")
    p.add_argument('--cr', action='store_true', default=False,
                   help='deprecated alias for --pass cr')
    p.add_argument('--center-band', default='f814W',
                   help='band whose deflector peak defines the shared cutout centre for '
                        'ALL bands of a lens, so they co-register. Default f814W (highest '
                        'S/N, GAIA-accurate). Faint bands recentred on their own brightest '
                        'pixel mis-centre (F606W locks onto a ring knot ~0.45" off).')
    p.add_argument('--center-self', action='store_true', default=False,
                   help='recentre each band on its own brightest pixel (old behaviour), '
                        'not on the shared --center-band centre')
    p.add_argument('--corr-factor', type=float, default=1.0,
                   help='multiply the noise map by this factor to absorb drizzle '
                        'correlated noise (default 1.0 = off, leaving a pure per-pixel '
                        'sigma). Drizzling correlates neighbouring output pixels, so a '
                        'diagonal-covariance likelihood such as PyAutoLens understates '
                        'integrated-flux uncertainties by ~1.24 (ACS F814W) and ~1.17 '
                        '(WFC3/IR F160W), measured by a blank-sky block-sum test. Set '
                        'this per band to correct that; leave at 1.0 if the covariance '
                        'is handled at the modelling stage instead.')
    p.add_argument('--output', default=None,
                   help='output dir, default data/cutouts/<sample>/<lens>/<filt>')
    a = p.parse_args()

    drizzled_dir = os.path.join(ws_path, 'data', 'drizzled', a.sample, a.lens, a.filt)
    output_dir = a.output or os.path.join(ws_path, 'data', 'cutouts', a.sample, a.lens, a.filt)
    os.makedirs(output_dir, exist_ok=True)

    # Resolve which pass to cut from. --cr is a deprecated alias for --pass cr.
    requested = 'cr' if a.cr else a.drizzle_pass
    if requested == 'auto':
        has_cr = bool(glob.glob(os.path.join(drizzled_dir, '*_cr_*_sci.fits')))
        drizzle_pass = 'cr' if has_cr else 'nocrrej'
    else:
        drizzle_pass = requested

    # Distinct output names by pass so the two never clobber and can be compared in the
    # same directory: cutout_cr_* for CR, cutout_* for no-CR.
    prefix = 'cutout_cr' if drizzle_pass == 'cr' else 'cutout'

    sci_file, wht_file = find_products(drizzled_dir, drizzle_pass)
    print(f"{a.lens} {a.filt}  [{drizzle_pass} pass]")
    print(f"  sci: {os.path.basename(sci_file)}")
    print(f"  wht: {os.path.basename(wht_file)}")

    # Drizzled products from this pipeline are single-image files: data in the primary HDU.
    with fits.open(sci_file) as hdul:
        sci_hdr = hdul[0].header
        sci_data = hdul[0].data
    with fits.open(wht_file) as hdul:
        wht_hdr = hdul[0].header
        wht_data = hdul[0].data

    wcs = WCS(sci_hdr)

    if a.lens not in slacs_coords:
        raise KeyError(f"{a.lens} not in slacs_coords (info/slacs_coords.py)")
    ra, dec = slacs_coords[a.lens]
    catalogue_coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')
    print(f"  catalogue position: {ra} {dec}")

    # Shared recentre. Every band of a lens is centred on the SAME sky point, taken
    # from the --center-band mosaic (default F814W: highest S/N, GAIA-accurate, cleanest
    # deflector peak). Recentring each band on its own brightest pixel mis-centres the
    # faint bands -- F606W locks onto a lensed-ring knot ~0.45" from the deflector core
    # -- so the bands would not co-register. --center-self restores per-band recentring.
    #
    # Peak-finding uses the CR-rejected pass when present (no cosmic rays to lock onto),
    # else the no-CR pass; see find_peak_coord.
    center_dir = os.path.join(ws_path, 'data', 'drizzled', a.sample, a.lens, a.center_band)
    if a.center_self or a.filt == a.center_band or not _has_products(center_dir):
        peak_dir, peak_src = drizzled_dir, 'this band'
        if not a.center_self and a.filt != a.center_band:
            print(f"  WARNING: no {a.center_band} products for a shared centre; "
                  f"recentring on this band instead")
    else:
        peak_dir, peak_src = center_dir, a.center_band

    try:
        peak_coord, cen_used = find_peak_coord(peak_dir, catalogue_coord,
                                               a.size, a.box, a.median_size)
    except NoOverlapError:
        px, py = wcs.world_to_pixel(catalogue_coord)
        raise SystemExit(
            f"ERROR: the catalogue position for {a.lens} falls outside the mosaic.\n"
            f"       It maps to pixel ({px:.0f}, {py:.0f}) in a "
            f"{sci_hdr['NAXIS1']}x{sci_hdr['NAXIS2']} image.\n"
            f"       The mosaic WCS does not cover the target - check the drizzle output."
        )
    print(f"  centre from {peak_src}: {os.path.basename(cen_used)}")

    scale = proj_plane_pixel_scales(wcs.celestial)[0] * 3600.0
    offset = catalogue_coord.separation(peak_coord).arcsec
    print(f"  pixel scale: {scale:.4f}\"/pix")
    print(f"  recentred: offset {offset:.3f}\" ({offset / scale:.1f} pix) "
          f"from catalogue position")

    # Second pass: re-cut both sci and weight on the recentred position, so the stamp
    # keeps its full requested size rather than being shifted and trimmed.
    sci_cutout = make_cutout(sci_data, wcs, peak_coord, a.size)
    wht_cutout = make_cutout(wht_data, WCS(wht_hdr), peak_coord, a.size)

    # Calibrate 1/sqrt(WHT) into a real sigma map. This is a units correction, not a
    # tuning knob: see weight_to_sigma_scale().
    units_k, note = weight_to_sigma_scale(sci_hdr)
    print(f"  noise scale: {note}")
    # WFPC2 weight maps come in three generations and are not separable by inspection,
    # so key the warning on the IVMMODEL stamp written by drizzle_wfpc2_wf3.py.
    # Absent = pre-fix product: either exptime-only weighting (block-sum 5e-4) or an
    # IVM built from the file's ERR array, which is exactly sqrt(SCI) and so overstates
    # the noise ~2.1x (block-sum 0.46). Both need a re-drizzle, not a scale factor.
    if sci_hdr.get('INSTRUME', '').strip() == 'WFPC2':
        model = str(sci_hdr.get('IVMMODEL', '')).strip()
        if model == 'SCI/gain+floor^2':
            print(f"  WFPC2 noise model: {model} (calibrated; block-sum 0.90-1.09)")
        else:
            print("  WARNING: this WFPC2 product predates the noise-model fix "
                  f"(IVMMODEL={model or 'absent'}). Its noise map is NOT calibrated -- "
                  "block-sum 0.46 (old IVM) or 5e-4 (exptime-only) where 1.0 is "
                  "correct. Re-drizzle before using it for a likelihood.")
    if a.corr_factor != 1.0:
        print(f"  correlated-noise inflation: x{a.corr_factor:g}")
    noise_data = noise_map_via_weight_map_from(wht_cutout.data.astype(np.float64),
                                               scale=units_k * a.corr_factor)

    noise_hdr = wht_hdr.copy()
    noise_hdr['NOISEK'] = (units_k, 'ERR-units scale applied to 1/sqrt(WHT)')
    noise_hdr['NOISECOR'] = (a.corr_factor, 'drizzle correlated-noise inflation applied')

    write_cutout(sci_cutout.data, sci_cutout.wcs, sci_hdr,
                 os.path.join(output_dir, f'{prefix}_sci.fits'))
    write_cutout(noise_data, wht_cutout.wcs, noise_hdr,
                 os.path.join(output_dir, f'{prefix}_noise.fits'))

    plot_cutouts(sci_cutout.data, noise_data,
                 os.path.join(output_dir, f'{prefix}.png'),
                 title=f"{a.lens}  {a.filt}  [{drizzle_pass}]  ({a.size:g}\" cutout, "
                       f"recentred {offset:.2f}\" from catalogue)")


if __name__ == '__main__':
    main()
