"""
Deflector light removed with an ISOPHOTE model, as a display for finding lensed images.

Why not the radial median (`make_masks.radial_median_subtract`): a circular profile leaves a
quadrupole on every elliptical deflector and manufactures features at a bright neighbour's
azimuth (AGENTS.md *Masks*, method lessons). Fitting elliptical isophotes whose ellipticity and
position angle vary with radius (photutils `Ellipse`, centre fixed on the deflector, masked
pixels excluded) removes the deflector far more cleanly, so arcs and counter-images inside the
envelope stand out. Adopted for the positions GUI on 2026-11-01, after the review sheets that
used it (diagnostics/feature_vetting/review_2026-11-01/) read far better than radial-subtracted.

It is a FINDING AID, like the radial median: it is not photometry, it absorbs a little of any
arc lying along an isophote, and lensed light biases the isophotes it crosses. The model is fitted
on a central crop (+-6", or the stamp if smaller), so on a large stamp the outskirts are left
unsubtracted. If photutils finds no fit for any of the starting shapes, the image is returned
unchanged (and a note is printed), so a display never fails outright.
"""
import os
import sys

import numpy as np
from scipy.ndimage import gaussian_filter

FIT_HALF_ARCSEC = 6.0
# (sma0 arcsec, ellipticity, position angle rad) tried in turn until photutils finds a fit
STARTS = [(0.5, 0.15, 0.5), (1.0, 0.3, 0.0), (0.3, 0.05, 1.5), (1.5, 0.4, 2.5)]


def isophote_subtract(image, pixel_scale, exclude=None, centre=None):
    """`image` minus its isophote model. `exclude` (bool, True = ignore) keeps masked pixels out
    of the fit; `centre` (col, row) is refined by a windowed centroid, default the stamp centre."""
    from photutils.isophote import EllipseGeometry, Ellipse, build_ellipse_model

    img = np.nan_to_num(np.asarray(image, dtype=float))
    msk = np.zeros(img.shape, bool) if exclude is None else np.asarray(exclude, dtype=bool)
    ny, nx = img.shape
    y, x = np.mgrid[:ny, :nx]
    cx, cy = centre if centre is not None else ((nx - 1) / 2.0, (ny - 1) / 2.0)
    for _ in range(5):
        sel = (np.hypot(x - cx, y - cy) < 0.3 / pixel_scale) & ~msk
        w = np.clip(img, 0, None) * sel
        if w.sum() <= 0:
            break
        cx, cy = (w * x).sum() / w.sum(), (w * y).sum() / w.sum()

    half = int(min(FIT_HALF_ARCSEC / pixel_scale, cx, cy, nx - 1 - cx, ny - 1 - cy))
    x0, y0 = int(round(cx)) - half, int(round(cy)) - half
    sub = img[y0:y0 + 2 * half + 1, x0:x0 + 2 * half + 1]
    smsk = msk[y0:y0 + 2 * half + 1, x0:x0 + 2 * half + 1]
    model = np.zeros_like(img)
    for sma0, eps, pa in STARTS:
        try:
            geometry = EllipseGeometry(x0=cx - x0, y0=cy - y0, sma=sma0 / pixel_scale, eps=eps, pa=pa)
            isolist = Ellipse(np.ma.masked_array(sub, mask=smsk), geometry).fit_image(
                sma0=sma0 / pixel_scale, minsma=1, maxsma=half - 3, step=0.12, fix_center=True,
                maxgerr=1.0)
            if len(isolist) < 3:
                continue
            model[y0:y0 + 2 * half + 1, x0:x0 + 2 * half + 1] = build_ellipse_model(
                sub.shape, isolist, fill=0.0)
            break
        except Exception:
            continue
    else:
        print('  NOTE: no isophote fit for any starting shape -- display left unsubtracted')
    return img - model


def deflector_pixel(sci_path):
    """(col, row) of the catalogue deflector on a stamp whose centre is NOT the deflector, else
    None (the stamp centre is then right). Only the group-centred `slacs_multi` stamps
    (`stamp_center` in lens_samples.json) qualify: centred on the group, an isophote fit started
    at the stamp centre finds no galaxy there and fails. The lens and sample are read from the
    path, data/stamps/<sample>/<lens>/<band>/<file>."""
    from astropy.coordinates import SkyCoord
    from astropy.io import fits
    from astropy.wcs import WCS
    import astropy.units as u
    scripts = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, scripts)
    sys.path.insert(0, os.path.join(os.path.dirname(scripts), 'info', 'coords'))
    import mast_target_names
    from lens_coords import lens_coords
    band_dir = os.path.dirname(os.path.abspath(sci_path))
    lens = os.path.basename(os.path.dirname(band_dir))
    sample = os.path.basename(os.path.dirname(os.path.dirname(band_dir)))
    try:
        if not mast_target_names.stamp_center(lens, sample) or lens not in lens_coords:
            return None
    except KeyError:
        return None
    wcs = WCS(fits.getheader(sci_path)).celestial
    c = SkyCoord(*lens_coords[lens], unit=(u.hourangle, u.deg), frame=wcs.wcs.radesys.lower() or 'icrs')
    return tuple(float(v) for v in wcs.world_to_pixel(c))


def isophote_display(image, pixel_scale, exclude=None, smooth_px=1.5, centre=None):
    """The positions-GUI panel: the isophote residual smoothed by `smooth_px` (a Gaussian sigma)
    and its symmetric display range, so positive (lensed) light shows red and over-subtraction
    blue on the diverging map. `centre` (col, row) starts the fit off-centre (`deflector_pixel`).
    Returns (array, vmax)."""
    res = gaussian_filter(isophote_subtract(image, pixel_scale, exclude=exclude, centre=centre),
                          smooth_px)
    keep = ~np.asarray(exclude, dtype=bool) if exclude is not None else np.ones(res.shape, bool)
    vmax = float(np.nanpercentile(np.abs(res[keep]), 99.5)) or 1.0
    return res, vmax
