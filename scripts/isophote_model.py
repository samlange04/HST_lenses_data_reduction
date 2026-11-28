"""
Deflector light removed with an ISOPHOTE model, as a display for finding lensed images.

Why not the radial median (`make_masks.radial_median_subtract`): a circular profile leaves a
quadrupole on every elliptical deflector and manufactures features at a bright neighbour's
azimuth (AGENTS.md *Masks*, method lessons). Fitting elliptical isophotes whose ellipticity and
position angle vary with radius (photutils `Ellipse`, centre fixed on the deflector, masked
pixels excluded) removes the deflector far more cleanly, so arcs and counter-images inside the
envelope stand out. Adopted for the positions GUI on 2026-11-01, after the review sheets that
used it read far better than radial-subtracted ones.

It is a FINDING AID, like the radial median: it is not photometry, it absorbs a little of any
arc lying along an isophote, and lensed light biases the isophotes it crosses. The model is fitted
on a central crop (+-6" by default, or `fit_half_arcsec`; the band overview fits only its zoom),
so on a large stamp the outskirts are left unsubtracted. If photutils finds no fit for any of
the starting shapes, the image is returned unchanged (and a note is printed), so a display never
fails outright.

Merged 2026-11-28 across the HST and JWST repos (the file is identical in both): the HST repo
brought `deflector_pixel` (off-centre deflectors) and the display cache, the JWST repo the
`refine_centre` / `fit_half_arcsec` / `step` / `maxit` knobs (band-overview speed).
"""
import os
import sys

import numpy as np
from scipy.ndimage import gaussian_filter

FIT_HALF_ARCSEC = 6.0
# (sma0 arcsec, ellipticity, position angle rad) tried in turn until photutils finds a fit
STARTS = [(0.5, 0.15, 0.5), (1.0, 0.3, 0.0), (0.3, 0.05, 1.5), (1.5, 0.4, 2.5)]


def isophote_subtract(image, pixel_scale, exclude=None, centre=None, refine_centre=None,
                      fit_half_arcsec=FIT_HALF_ARCSEC, step=0.12, maxit=50):
    """`image` minus its isophote model. `exclude` (bool, True = ignore) keeps masked pixels out
    of the fit; `centre` (col, row) is the fit's starting centre, default the stamp centre.
    `refine_centre` (default None) refines the start by an iterated 0.3" windowed centroid --
    but an explicit `centre` is a MEASURED position (the deflector, `deflector_pixel`), and the
    iteration creeps toward bright asymmetric light such as an arc or ring inside that window,
    so the default refines only the unmeasured stamp-centre case (refine_centre = centre is
    None). `fit_half_arcsec` sets the fitted crop. `step` (the isophotes' geometric spacing)
    and `maxit` (iterations per isophote) trade precision for time."""
    from photutils.isophote import EllipseGeometry, Ellipse, build_ellipse_model

    img = np.nan_to_num(np.asarray(image, dtype=float))
    msk = np.zeros(img.shape, bool) if exclude is None else np.asarray(exclude, dtype=bool)
    ny, nx = img.shape
    y, x = np.mgrid[:ny, :nx]
    cx, cy = centre if centre is not None else ((nx - 1) / 2.0, (ny - 1) / 2.0)
    if refine_centre is None:
        refine_centre = centre is None
    for _ in range(5 if refine_centre else 0):
        sel = (np.hypot(x - cx, y - cy) < 0.3 / pixel_scale) & ~msk
        w = np.clip(img, 0, None) * sel
        if w.sum() <= 0:
            break
        cx, cy = (w * x).sum() / w.sum(), (w * y).sum() / w.sum()

    half = int(min(fit_half_arcsec / pixel_scale, cx, cy, nx - 1 - cx, ny - 1 - cy))
    x0, y0 = int(round(cx)) - half, int(round(cy)) - half
    sub = img[y0:y0 + 2 * half + 1, x0:x0 + 2 * half + 1]
    smsk = msk[y0:y0 + 2 * half + 1, x0:x0 + 2 * half + 1]
    model = np.zeros_like(img)
    for sma0, eps, pa in STARTS:
        try:
            geometry = EllipseGeometry(x0=cx - x0, y0=cy - y0, sma=sma0 / pixel_scale, eps=eps, pa=pa)
            isolist = Ellipse(np.ma.masked_array(sub, mask=smsk), geometry).fit_image(
                sma0=sma0 / pixel_scale, minsma=1, maxsma=half - 3, step=step, fix_center=True,
                maxgerr=1.0, maxit=maxit)
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
    """(col, row, refine) of the deflector on a stamp whose centre is NOT exactly the
    deflector, else None (the isophote fit then refines from the stamp centre, as always).
    What qualifies is the profile's call (`deflector_radec`): in the HST repo only the
    group-centred multi_galaxy stamps (`stamp_center`) -- an isophote fit started at a group
    centre finds no galaxy there and fails -- with refine=True (the catalogue position is
    approximate); in the JWST repo the measured deflector position (`deflector_ra/dec`)
    whenever one exists, with refine=False (it is measured to <= 39 mas, and the refinement
    creeps toward arc light). The lens and sample are read from the path's band dir (lens
    names are unique across samples)."""
    from astropy.coordinates import SkyCoord
    from astropy.io import fits
    from astropy.wcs import WCS
    scripts = os.path.dirname(os.path.abspath(__file__))
    if scripts not in sys.path:
        sys.path.insert(0, scripts)
    try:
        import lens_samples as profile   # the JWST repo's profile module
    except ImportError:
        import mast_target_names as profile   # the HST repo's
    lens = os.path.basename(os.path.dirname(os.path.dirname(os.path.abspath(sci_path))))
    sample = profile.sample_of(lens)
    try:
        got = profile.deflector_radec(lens, sample)
    except (KeyError, SystemExit):
        return None
    if got is None:
        return None
    ra, dec, exact = got
    wcs = WCS(fits.getheader(sci_path)).celestial
    c = SkyCoord(ra, dec, unit='deg', frame=wcs.wcs.radesys.lower() or 'icrs')
    col, row = (float(v) for v in wcs.world_to_pixel(c))
    return col, row, not exact


_DISPLAY_CACHE = {}


def isophote_display(image, pixel_scale, exclude=None, smooth_px=1.5, centre=None,
                     refine_centre=None):
    """The positions-GUI panel: the isophote residual smoothed by `smooth_px` (a Gaussian sigma)
    and its symmetric display range, so positive (lensed) light shows red and over-subtraction
    blue on the diverging map. `centre` (col, row) starts the fit off-centre (`deflector_pixel`;
    pass its `refine` as `refine_centre`). Returns (array, vmax). Memoised on the inputs, so
    the positions GUI's fit is reused for the marked band's QC PNG written straight after it
    (each fit takes a few seconds)."""
    import hashlib
    img = np.asarray(image, dtype=float)
    key = hashlib.sha1(img.tobytes() + (b'' if exclude is None else np.asarray(exclude, bool).tobytes())
                       + repr((img.shape, pixel_scale, smooth_px, centre, refine_centre)).encode()).hexdigest()
    if key in _DISPLAY_CACHE:
        return _DISPLAY_CACHE[key]
    res = gaussian_filter(isophote_subtract(img, pixel_scale, exclude=exclude, centre=centre,
                                            refine_centre=refine_centre),
                          smooth_px)
    keep = ~np.asarray(exclude, dtype=bool) if exclude is not None else np.ones(res.shape, bool)
    vmax = float(np.nanpercentile(np.abs(res[keep]), 99.5)) or 1.0
    _DISPLAY_CACHE[key] = (res, vmax)
    return res, vmax
