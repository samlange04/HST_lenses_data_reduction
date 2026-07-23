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
from astropy.visualization import simple_norm
from astropy.wcs.utils import proj_plane_pixel_scales
from scipy.ndimage import median_filter
import astropy.units as u

ws_path = '/Users/samlange/Code/data_reduction'
sys.path.insert(0, os.path.join(ws_path, 'info'))
from slacs_coords import slacs_coords


def noise_map_via_weight_map_from(weight_map):
    """
    Setup the noise-map from a weight map, which is a form of noise-map that comes via HST
    image-reduction and the software package MultiDrizzle.

    The variance in each pixel is computed as:

    Variance = 1.0 / sqrt(weight_map).

    The weight map may contain zeros, in which case the variances are converted to large
    values to omit them from the analysis.

    Parameters
    ----------
    weight_map
        The weight-value of each pixel which is converted to a variance.
    """
    np.seterr(divide="ignore")
    noise_map = 1.0 / weight_map ** 0.5
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
    The signal panel is log-stretched for the same reason: it makes the faint envelope
    and any arcs visible rather than a single saturated core.
    """
    noise_display = np.where(noise_data >= 1.0e8, np.nan, noise_data)
    with np.errstate(divide='ignore', invalid='ignore'):
        snr = sci_data / noise_data

    panels = [
        ('signal',        sci_data,      'log'),
        ('noise',         noise_display, 'linear'),
        ('signal / noise', snr,          'linear'),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.5))
    for ax, (label, data, stretch) in zip(axes, panels):
        norm = simple_norm(data[np.isfinite(data)], stretch, percent=99.0)
        im = ax.imshow(data, norm=norm, origin='upper', cmap='viridis')
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
    p.add_argument('--sample', default='slacs')
    p.add_argument('--size',   type=float, default=10.0,
                   help='cutout size in arcsec (square), default 10.0')
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
                        "uses LACosmic masking that keeps ~99% of the deflector core "
                        "while removing cosmic rays, so the old no-CR default put "
                        "CR-riddled stamps in front of the user.")
    p.add_argument('--cr', action='store_true', default=False,
                   help='deprecated alias for --pass cr')
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

    # Peak-finding runs on the CR-rejected mosaic when one exists, even when the
    # science stamp is cut from the no-CR pass. Cosmic rays in the no-CR data
    # otherwise steal the recentre outright: at median_size=5 this mis-centred 3 of
    # 21 WFPC2 F606W lenses (J0008-0004 landed on a CR streak with the galaxy in the
    # stamp corner) and 3 of 38 ACS F814W lenses. Widening the median window far
    # enough to suppress the CRs (21 pix at 0.05"/px) starts biasing the peak of a
    # compact galaxy, so use data that has no CRs in it instead. The CR pass loses
    # core flux and is unfit for science, but the core is still the local maximum,
    # which is all that is needed to locate it.
    centring_data, centring_wcs = sci_data, wcs
    try:
        cen_sci, _ = find_products(drizzled_dir, 'cr')
    except FileNotFoundError:
        cen_sci = None
    if cen_sci and cen_sci != sci_file:
        with fits.open(cen_sci) as hdul:
            centring_data = hdul[0].data
            centring_wcs = WCS(hdul[0].header)
        print(f"  centring on: {os.path.basename(cen_sci)}")
    else:
        print("  centring on: science pass (no CR-rejected mosaic available)")

    # First pass: cut on the catalogue coordinate, then find the lens galaxy peak.
    try:
        first = make_cutout(centring_data, centring_wcs, catalogue_coord, a.size)
    except NoOverlapError:
        px, py = wcs.world_to_pixel(catalogue_coord)
        raise SystemExit(
            f"ERROR: the catalogue position for {a.lens} falls outside the {a.filt} mosaic.\n"
            f"       It maps to pixel ({px:.0f}, {py:.0f}) in a "
            f"{sci_hdr['NAXIS1']}x{sci_hdr['NAXIS2']} image.\n"
            f"       The mosaic WCS does not cover the target - check the drizzle output."
        )
    x, y = brightest_pixel_position(first, box=a.box, median_size=a.median_size)
    peak_coord = first.wcs.pixel_to_world(x, y)

    scale = proj_plane_pixel_scales(wcs.celestial)[0] * 3600.0
    offset = catalogue_coord.separation(peak_coord).arcsec
    print(f"  pixel scale: {scale:.4f}\"/pix")
    print(f"  recentred onto brightest pixel: offset {offset:.3f}\" "
          f"({offset / scale:.1f} pix) from catalogue position")

    # Second pass: re-cut both sci and weight on the recentred position, so the stamp
    # keeps its full requested size rather than being shifted and trimmed.
    sci_cutout = make_cutout(sci_data, wcs, peak_coord, a.size)
    wht_cutout = make_cutout(wht_data, WCS(wht_hdr), peak_coord, a.size)

    noise_data = noise_map_via_weight_map_from(wht_cutout.data.astype(np.float64))

    write_cutout(sci_cutout.data, sci_cutout.wcs, sci_hdr,
                 os.path.join(output_dir, f'{prefix}_sci.fits'))
    write_cutout(noise_data, wht_cutout.wcs, wht_hdr,
                 os.path.join(output_dir, f'{prefix}_noise.fits'))

    plot_cutouts(sci_cutout.data, noise_data,
                 os.path.join(output_dir, f'{prefix}.png'),
                 title=f"{a.lens}  {a.filt}  [{drizzle_pass}]  ({a.size:g}\" cutout, "
                       f"recentred {offset:.2f}\" from catalogue)")


if __name__ == '__main__':
    main()
