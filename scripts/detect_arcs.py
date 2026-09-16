#!/usr/bin/env python
"""
Automatic arc detection: propose the lensed-source region for `make_arc_masks.py` to review.

WHAT THIS IS FOR. Drawing an arc mask from a blank canvas is the slow half of
`make_arc_masks.py`. This module finds the arcs automatically and hands the result over as a
PROPOSAL -- an outline the GUI draws over the image, which you then correct with the green
(add) / red (erase) brushes and accept, reject or skip, exactly like `make_masks.py`'s
cross-band proposals. Nothing here ever writes an arc mask: the detector proposes, the human
disposes (`make_arc_masks.py --propose-from detect`).

THE METHOD IS COLOUR, NOT BRIGHTNESS. A SLACS deflector is a red early-type; its lensed
source is a blue star-forming galaxy. That is the one discriminant that survives a bad
galaxy model, and it is how SLACS found these systems in the first place (Bolton+2006
difference imaging). So instead of hunting for residual flux in one band -- which finds the
deflector's own model error, every red neighbour, and the sky gradient -- we build a
BLUE-EXCESS image:

    be = B' - k * R'

where R' and B' are the red (f814W) and blue (f606W/f555W) stamps PSF-MATCHED to a common
resolution (each convolved with the OTHER band's kernel, so both end up with PSF_R * PSF_B),
and k is the deflector's own blue/red flux ratio measured ON EACH ISOPHOTE -- keyed on the red
band's own brightness, so it needs no galaxy model, no isophote shape, and no assumption that
the galaxy is smooth. By construction the deflector cancels everywhere, INCLUDING its colour
gradient, while anything bluer than it stands up positive. k is a ratio measured on the data,
so the bands' different units (WFPC2 DN/s vs ACS e/s -- see AGENTS.md "F606W BUNIT") cancel and
never need converting.

THIS TOOL DOES NOT RE-REGISTER THE BANDS, AND THAT IS DELIBERATE. `align_wfpc2_to_acs.py`
owns the astrometry, and its tie holds -- three independent measurements agree that f606W and
f814W are aligned to about a third of a pixel:
  * cross-matched FIELD SOURCES (43-175 per lens, the only check the deflector cannot bias):
    median vector offset 0.002-0.026" = 0.04-0.5px;
  * the deflector centroid by windowed centre-of-mass: 0.05-0.42px, and identical in the
    mosaic and in the cutout, so make_cutouts carries the WCS through correctly;
  * this file's own per-lens diagnostic: median 0.31px, max 0.69px, nothing over the flag.
So there is no offset to correct here, and the band offset is MEASURED AND REPORTED
(`band_offset_px`) rather than applied. Two wrong turns are recorded in the code where they
happened: a flux-weighted centroid that the arc walked off and "measured" 2.2px on a lens whose
bands agree to 0.04px, and a shift fitted by minimising the colour residual, which an imperfect
PSF match can always reward -- it ran to its limit on half the test lenses and did not reduce
the dipole it was meant to remove. A measurement that says the astrometry is broken is nearly
always the estimator: check it against something the deflector cannot bias before believing it.

WHAT ACTUALLY WAS BROKEN WAS THE PSF KERNELS. Convolving with an off-centre kernel translates
the image, and these kernels are not centred: 15 of the 90 slacs_gold cutout PSFs are off their
own array centre by >0.5px and 5 by >1px, worst ~1.7px, concentrated in the INJECTED builds
(`inject_wfpc2_psfdb` median 0.44px, `inject_acs_fdpsf` 1.16px, against 0.15-0.32px for
`inject_stdpsf` and `empirical`). In a kernel swap each band is moved by the OTHER band's
kernel offset, which is what the fictitious 2.2px was really made of. `centre_kernel` fixes it
here; anything else that convolves these kernels should do the same.

WHAT IS LEFT SHOWS UP AS A DIPOLE, AND IS FLAGGED RATHER THAN HIDDEN. A residual centring or
PSF mismatch is ANTISYMMETRIC about the deflector -- positive one side, negative the other --
where lensing never is: opposite an arc lies a counter-image or empty sky, never a deficit.
So `core_dipole` reports the antisymmetric fraction just inside the arc annulus, and any
component whose own mirror reads below MIRROR_DEFICIT_SNR is dropped as the bright half of a
dipole. On lenses with a strong dipole (core_dipole >~ 2, about a third of the sample) the
inner region is subtraction error and the detector proposes nothing there: those inner images
have to be drawn by hand.

WHERE IT LOOKS IS SET BY THE MEASURED EINSTEIN RADIUS. Candidates are kept only in an annulus
around the Auger+2009 (SLACS IX) theta_E for that lens -- `info/lens_einstein_radii.json`,
built from VizieR J/ApJ/705/1099/lenses by `--fetch-theta-e`. Every slacs_gold lens has a
measured value, and a lens in that table has a successful published lens model, so the arc is
known to exist and roughly where. This is the cheapest prior available and it kills the two
classic false positives at once: the deflector's central residual (inside the annulus) and
unrelated field galaxies (outside it).

THREE TRAPS THIS CODE IS WRITTEN AROUND, all of which produced wrong answers here before:
  1. An integrated-aperture S/N invents sources. Every component is judged on its PEAK
     PER-PIXEL S/N and its connected area above threshold, never on a flux sum.
  2. Compact-source detection is blind to a diffuse arc. Detection runs on a mildly smoothed
     map (`--smooth`, 0.07" by default -- an arc's width, not a point source's), so
     low-surface-brightness emission survives.
  3. The propagated variance map is NOT the noise of these images. Drizzle correlates
     neighbouring pixels and the PSF match correlates them further, so
     var -> conv(var, kernel^2) underestimates the real scatter by a factor of 3-4. Both S/N
     maps are therefore rescaled by the MAD width of their own empty-sky background
     (`alpha_*` in the provenance): a background pixel reads S/N ~ 1 by construction, so a
     threshold means what it says.
  Plus one data trap: a few noise maps flag zero-coverage pixels with 1e8 (12 products
  repo-wide, a handful of pixels each). Left in, that one value destroys the whole frame
  through the FFT. They are detected as sentinels, healed for the convolution, and excluded
  from the detection and from the background statistics.

Usage:
    # one-off: cache the Auger+2009 Einstein radii (needs network; already committed)
    uv run python scripts/detect_arcs.py --fetch-theta-e

    # test the detector over a sample: QC figure per lens + a summary JSON, no products
    uv run python scripts/detect_arcs.py --sample slacs_gold
    uv run python scripts/detect_arcs.py --lens J0330-0020 --snr 3.0    # tune one lens

    # then review/correct the proposals by hand, one lens at a time:
    uv run python scripts/make_arc_masks.py --lens J0330-0020 --propose-from detect
"""

import argparse
import json
import os
import sys
import warnings

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import (binary_closing, binary_dilation, binary_fill_holes,
                           binary_opening, gaussian_filter, label, map_coordinates)
from scipy.signal import fftconvolve

warnings.filterwarnings('ignore', category=UserWarning, module='autonerves')

ws_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mast_target_names
import cutout_paths

THETA_E_JSON = os.path.join(ws_path, 'info', 'lens_einstein_radii.json')

#: Bands bluer than the reference band, best first. The blue band carries the source; the red
#: band carries the deflector and sets the colour reference.
BLUE_PRIORITY = ['f606W', 'f555W', 'f606W_v2', 'f438W']
#: Reference (red, deflector-dominated) bands, best first. f160W is redder still but sits on a
#: different grid and a much broader PSF, so it is not used as the reference.
RED_PRIORITY = ['f814W', 'f555W']

#: Noise-map values this far above the frame's median mark zero-coverage pixels (the maps use
#: 1e8). One such pixel left in an FFT convolution ruins the entire frame.
SENTINEL_FACTOR = 1e3

#: Largest sub-pixel shift this tool will apply when re-registering the blue band on the red.
#: The bands are already tied by align_wfpc2_to_acs.py and measure within ~0.5px of each other
#: (field sources and deflector centroids agree); a larger measurement means a broken tie or a
#: fooled centroid, neither of which a detector should silently absorb.
MAX_REGISTRATION_PX = 0.8

#: A component whose 180-degree mirror about the deflector reads below this in calibrated S/N
#: is the positive half of a centring/PSF dipole, not a source: lensing puts a counter-image or
#: empty sky opposite an arc, never a deficit.
MIRROR_DEFICIT_SNR = -2.0

DEFAULTS = dict(snr=3.5, peak=4.0, min_area=15, smooth=0.07, r_in=0.35, r_out=2.0, grow=2)


# --------------------------------------------------------------------------- theta_E table

def fetch_theta_e(sample_lenses=None, timeout=60):
    """Fetch Auger+2009 (SLACS IX) Einstein radii from VizieR and cache them as arcsec.

    The catalogue column `RE` is in kpc, so it is converted with the angular diameter
    distance to the lens: theta_E["] = RE[kpc] / D_A(z_lens)[kpc] * 206265, on the flat
    (H0=70, Om=0.3) cosmology the paper assumes. Writes `info/lens_einstein_radii.json`.
    """
    import urllib.request
    from astropy.cosmology import FlatLambdaCDM
    url = ('https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/ApJ/705/1099/lenses'
           '&-out=SDSS,zlens,zsrc,sigma,RE&-out.max=500')
    with urllib.request.urlopen(url, timeout=timeout) as fh:
        text = fh.read().decode('utf-8', 'replace')
    cos = FlatLambdaCDM(H0=70, Om0=0.3)
    out = {}
    for line in text.splitlines():
        if not line.startswith('J'):
            continue
        cols = [c.strip() for c in line.split('\t')]
        if len(cols) < 5 or not cols[4] or not cols[1]:
            continue
        name, zl, zs, sig, re_kpc = cols[:5]
        d_a = cos.angular_diameter_distance(float(zl)).to('kpc').value
        out[name] = {
            'theta_e_arcsec': round(float(re_kpc) / d_a * 206265.0, 4),
            'RE_kpc': float(re_kpc), 'zlens': float(zl),
            'zsrc': float(zs) if zs else None, 'sigma_kms': int(sig) if sig else None,
            'source': 'Auger+2009 (SLACS IX), VizieR J/ApJ/705/1099/lenses',
        }
    with open(THETA_E_JSON, 'w') as fh:
        json.dump(dict(sorted(out.items())), fh, indent=1)
        fh.write('\n')
    print(f"wrote {THETA_E_JSON}: {len(out)} lenses with a measured Einstein radius")
    if sample_lenses:
        missing = [l for l in sample_lenses if l not in out]
        print(f"  of the {len(sample_lenses)} requested, {len(missing)} have no row"
              + (f": {', '.join(missing)}" if missing else ''))
    return out


def theta_e_table():
    """The cached Einstein-radius table (arcsec), or {} if it has not been fetched."""
    if not os.path.exists(THETA_E_JSON):
        return {}
    with open(THETA_E_JSON) as fh:
        return {k: v['theta_e_arcsec'] for k, v in json.load(fh).items()}


# ------------------------------------------------------------------------- band loading

def _sentinel_free(noise):
    """Return (variance, valid) with zero-coverage sentinel pixels healed and flagged.

    A handful of noise maps carry 1e8 where the drizzle weight is zero. Squared and pushed
    through an FFT convolution, that single value swamps the transform and returns garbage
    (verified: J0252+0039 f606W, 0.005% of pixels, made the whole blue-excess frame noise).
    The sentinels are replaced by the frame's median variance -- a value, not a hole, so the
    convolution stays well conditioned -- and returned as `valid=False` so nothing is ever
    detected there and they never enter the background statistics.
    """
    noise = np.asarray(noise, dtype=float)
    valid = np.isfinite(noise) & (noise > 0)
    if valid.any():
        med = float(np.median(noise[valid]))
        valid &= noise < SENTINEL_FACTOR * med
    else:
        med = 1.0
    var = np.where(valid, noise, med) ** 2
    return var, valid


def load_band(cutout_dir, prefix):
    """sci / variance / valid / psf / header for one band's cutout, PSF normalised to unit sum."""
    sci = fits.getdata(os.path.join(cutout_dir, f'{prefix}_sci.fits')).astype(float)
    hdr = fits.getheader(os.path.join(cutout_dir, f'{prefix}_sci.fits'))
    var, valid = _sentinel_free(fits.getdata(os.path.join(cutout_dir, f'{prefix}_noise.fits')))
    psf = centre_kernel(fits.getdata(os.path.join(cutout_dir, f'{prefix}_psf.fits')))
    mask_path = os.path.join(cutout_dir, f'{prefix}_mask.fits')
    cmask = (np.asarray(fits.getdata(mask_path), dtype=bool) if os.path.exists(mask_path)
             else np.zeros(sci.shape, dtype=bool))
    valid &= np.isfinite(sci)
    return dict(sci=np.where(valid, sci, 0.0), var=var, valid=valid, psf=psf, hdr=hdr,
                cmask=cmask)


def pick_bands(filt_dirs):
    """(red, blue) band names for the colour difference, or (red, None) if no blue band."""
    red = next((f for f in RED_PRIORITY if f in filt_dirs), None)
    blue = next((f for f in BLUE_PRIORITY if f in filt_dirs and f != red), None)
    return red, blue


def regrid(arr, src_hdr, dst_hdr, order=1):
    """Resample `arr` from the source cutout's WCS onto the target cutout's WCS grid.

    The same per-pixel sky-position transfer `make_masks.reproject_mask_bool` uses, in
    floating point: bilinear for data, order=0 for anything categorical. The optical bands
    share a grid to well under a pixel, but they are NOT index-aligned (a ~2px edge strip of
    each stamp falls outside the other), so an array copy would shear the colour difference
    by a fraction of a pixel across the frame -- which is exactly the scale a PSF-matched
    subtraction is sensitive to.
    """
    ny, nx = dst_hdr['NAXIS2'], dst_hdr['NAXIS1']
    yy, xx = np.mgrid[0:ny, 0:nx]
    sky = WCS(dst_hdr).celestial.pixel_to_world(xx.ravel(), yy.ravel())
    sx, sy = WCS(src_hdr).celestial.world_to_pixel(sky)
    vals = map_coordinates(np.nan_to_num(arr), np.vstack([sy, sx]), order=order,
                           mode='constant', cval=0.0)
    return vals.reshape(ny, nx)


# ------------------------------------------------------------------- geometry of the deflector

def centroid(sci, pixel_scales, r_arcsec=0.5, n_iter=5, start=None):
    """Deflector centroid: iterated centre-of-mass in a BOX with its local median removed.

    This is `align_wfpc2_to_acs.py`'s estimator, deliberately -- the two tools must not
    disagree about where the deflector is. A plain flux-weighted centroid over a circular
    aperture was tried here first and it is NOT safe on the blue band: with no local
    background removed it responds to every asymmetry in the frame, and iterating lets it
    WALK toward the arc. On J0029-0055 it read a 2.2px f606W-f814W offset that does not exist
    -- field sources there agree to 0.04px. A measurement that says the astrometry is broken
    is nearly always the estimator; check it against something the deflector cannot bias
    before believing it.
    """
    from photutils.centroids import centroid_com
    ny, nx = sci.shape
    cy, cx = start if start is not None else ((ny - 1) / 2.0, (nx - 1) / 2.0)
    half = max(int(round(r_arcsec / pixel_scales)), 3)
    for _ in range(n_iter):
        yi, xi = int(round(cy)), int(round(cx))
        if yi - half < 0 or xi - half < 0 or yi + half + 1 > ny or xi + half + 1 > nx:
            break
        box = np.nan_to_num(sci[yi - half:yi + half + 1, xi - half:xi + half + 1])
        box = np.clip(box - np.median(box), 0, None)
        if box.sum() <= 0:
            break
        dx, dy = centroid_com(box)
        if not (np.isfinite(dx) and np.isfinite(dy)):
            break
        cy, cx = yi - half + float(dy), xi - half + float(dx)
    return cy, cx


# ------------------------------------------------------------------------ the colour difference

def centre_kernel(psf):
    """Return `psf` shifted so its own centre sits on the array centre, renormalised.

    REQUIRED before the PSF swap, because these kernels are NOT centred: measured over the 90
    slacs_gold cutout PSFs, the median offset from the array centre is 0.2-0.5px by band, 15
    are off by >0.5px and 5 by >1px, the worst being f606W kernels at ~1.7px (J0822+2652,
    J2341+0000, J0029-0055) -- on several of them even the BRIGHTEST PIXEL is a whole pixel
    off centre, so this is not an estimator artefact.

    Convolving with an off-centre kernel TRANSLATES the image by that offset. In a kernel swap
    each band is therefore moved by the *other* band's kernel offset, and the difference
    image inherits the mismatch as a bright-plus-hole dipole at the core. Chasing that dipole
    with a centroid is what produced a fictitious '2.2px f606W astrometric offset' here; the
    astrometry was right and the kernels were off. Anything that convolves these kernels --
    not just this file -- should centre them first or know why it does not.
    """
    from photutils.centroids import centroid_com
    from scipy.ndimage import shift as ndshift
    k = np.nan_to_num(np.asarray(psf, dtype=float))
    ny, nx = k.shape
    py, px = np.unravel_index(np.argmax(k), k.shape)
    h = 3
    y0, x0 = max(py - h, 0), max(px - h, 0)
    box = np.clip(k[y0:py + h + 1, x0:px + h + 1], 0, None)
    if box.sum() <= 0:
        return k / k.sum()
    bx, by = centroid_com(box)
    dy, dx = (y0 + by) - (ny - 1) / 2.0, (x0 + bx) - (nx - 1) / 2.0
    if np.hypot(dy, dx) > 0.02 and np.isfinite(dy) and np.isfinite(dx):
        k = ndshift(k, (-dy, -dx), order=3, mode='constant', cval=0.0)
    total = k.sum()
    return k / total if total else k


def _convolve(data, var, kernel):
    """Convolve data with `kernel` and propagate the variance with `kernel**2`.

    The propagated variance is a LOWER BOUND on the truth -- it assumes independent pixels,
    and drizzled pixels are not -- which is why every S/N built from it is rescaled
    empirically afterwards (`calibrate`).
    """
    return (fftconvolve(np.nan_to_num(data), kernel, mode='same'),
            np.clip(fftconvolve(np.nan_to_num(var), kernel ** 2, mode='same'), 0.0, None))


def _mad(values):
    """Gaussian-equivalent sigma from the median absolute deviation (outlier-proof)."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if v.size < 50:
        return 1.0
    return float(max(1.4826 * np.median(np.abs(v - np.median(v))), 1e-6))


def colour_scale(red_m, blue_m, red_smooth, galaxy, sky_red):
    """k: the deflector's blue/red flux ratio as a function of the RED BAND'S BRIGHTNESS.

    LOG-spaced brightness bins, not equal-population ones. The core holds few pixels over a
    huge brightness range, so quantile bins collapse it into a single bin and force one colour
    on it -- which showed up as a broad NEGATIVE bowl over the middle of every lens, precisely
    where the inner arcs are. Log bins resolve the core; bins too sparse for a reliable median
    are merged into their neighbour rather than dropped, so the relation stays continuous.

    Bin assignment uses a SMOOTHED copy of the red band, so a pixel's own noise cannot scatter
    it into the wrong colour bin; the ratio itself is measured on the unsmoothed images.
    """
    from scipy.ndimage import median_filter as _median_filter
    k_map = np.zeros_like(red_m)
    if galaxy.sum() > 500:
        vals = red_smooth[galaxy]
        edges = np.geomspace(max(vals.min(), 2.0 * sky_red), vals.max(), 41)
        centres, ratios, pending = [], [], None
        for lo, hi in zip(edges[:-1], edges[1:]):
            lo = pending if pending is not None else lo
            sel = galaxy & (red_smooth >= lo) & (red_smooth < hi)
            if sel.sum() < 30:
                pending = lo                         # too few pixels: widen into the next bin
                continue
            pending = None
            centres.append(np.median(red_smooth[sel]))
            ratios.append(np.median(blue_m[sel]) / max(np.median(red_m[sel]), 1e-12))
        if len(centres) >= 3:
            ratios = _median_filter(np.asarray(ratios), size=3, mode='nearest')
            k_map = np.interp(red_smooth, np.asarray(centres), ratios,
                              left=ratios[0], right=ratios[-1])
    if not k_map.any():                              # degenerate frame: fall back to one ratio
        k_map = np.full_like(red_m, float(np.median(blue_m[galaxy]) /
                                          max(np.median(red_m[galaxy]), 1e-12))
                             if galaxy.any() else 1.0)
    return k_map


def colour_excess(filt_dirs, prefixes, lens, theta_e, smooth_arcsec=DEFAULTS['smooth']):
    """Build the blue-excess image and its calibrated S/N maps on the RED band's grid.

    Returns a dict with `be` (the blue-excess image), `snr` (smoothed, calibrated),
    `snr_raw` (unsmoothed, calibrated), the PSF-matched bands, the deflector geometry and
    the colour scaling k -- everything the detection step and the QC figure need.
    """
    red_f, blue_f = pick_bands(filt_dirs)
    if red_f is None or blue_f is None:
        raise ValueError(f'{lens}: need a red and a blue band, have {sorted(filt_dirs)}')
    red = load_band(filt_dirs[red_f], prefixes[red_f])
    blue = load_band(filt_dirs[blue_f], prefixes[blue_f])
    pixel_scales = float(abs(red['hdr']['CD1_1']) * 3600) if 'CD1_1' in red['hdr'] else None
    if pixel_scales is None:
        from astropy.wcs.utils import proj_plane_pixel_scales
        pixel_scales = float(proj_plane_pixel_scales(WCS(red['hdr']).celestial)[0] * 3600)

    # Blue onto the red grid, then PSF-match: each band convolved with the other's kernel
    # leaves both with the same effective PSF (PSF_R * PSF_B), so the difference is not a
    # sharp-minus-blurry dipole at every steep gradient -- above all at the galaxy core.
    b_sci = regrid(blue['sci'], blue['hdr'], red['hdr'])
    b_var = regrid(blue['var'], blue['hdr'], red['hdr'])
    b_valid = regrid(blue['valid'].astype(float), blue['hdr'], red['hdr'], order=0) > 0.5
    b_cmask = regrid(blue['cmask'].astype(float), blue['hdr'], red['hdr'], order=0) > 0.5
    red_m, red_var_m = _convolve(red['sci'], red['var'], blue['psf'])
    blue_m, blue_var_m = _convolve(b_sci, b_var, red['psf'])

    # SUB-PIXEL RE-REGISTRATION ON THE DEFLECTOR, and a hard limit on how far it may go.
    # The bands are already tied astrometrically by `align_wfpc2_to_acs.py` (F606W's CRVAL is
    # shifted onto F814W and stamped GSC240FX/ASTROREF), and that tie holds: cross-matched
    # FIELD SOURCES put f606W within 0.002-0.026" = 0.04-0.5px of f814W on every lens checked,
    # and the deflector centroids agree to 0.05-0.42px in the mosaics and in the cutouts
    # alike. So there is NO band offset to correct here -- only the last few tenths of a pixel,
    # which still matter because this subtraction is differential on a steep profile.
    #
    # Hence MAX_REGISTRATION_PX. Anything larger is not an astrometry problem this tool should
    # paper over: it is either a bad tie -- which belongs to align_wfpc2_to_acs.py, the tool
    # that owns the WCS -- or, far more often, a centroid being fooled. An earlier version
    # here shifted by whatever it measured and moved f606W by up to 2.2px on a lens whose
    # astrometry was fine, quietly dragging the arcs through the colour difference.
    cy, cx = centroid(red_m, pixel_scales)
    valid = red['valid'] & b_valid
    ny, nx = red_m.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    r_arcsec = np.hypot(yy - cy, xx - cx) * pixel_scales
    outskirts = r_arcsec > max(3.0 * theta_e, 2.5)
    sky_red = _mad(red_m[outskirts])
    red_smooth = gaussian_filter(red_m, 0.15 / pixel_scales, mode='nearest')

    # BAND OFFSET IS MEASURED AND REPORTED, NEVER APPLIED. `align_wfpc2_to_acs.py` owns the
    # astrometry and its tie holds -- see the module docstring for the three measurements. So
    # this tool does not re-register: it records what it sees, on the REGRIDDED BUT UNMATCHED
    # frames (a PSF-matched frame is not the sky, and reading a shift off one is how the 2.2px
    # ghost got in), and leaves any real offset to the tool that can fix it properly.
    #
    # Fitting a shift to minimise the colour residual was tried too and is worse than useless
    # here: with an imperfect PSF match a shift can always buy a smaller residual, so the fit
    # ran to its limit on 3 of 6 test lenses, disagreed with the field sources (0.89px where
    # they say 0.04px), and did not consistently reduce the core dipole it was meant to remove.
    # Both sides RAW: red['sci'] on its own grid against the regridded blue, neither
    # PSF-matched. Comparing a matched frame with an unmatched one measures the kernels as
    # much as the sky -- which is precisely the mistake that manufactured the 2.2px ghost.
    cy_r, cx_r = centroid(red['sci'], pixel_scales)
    cy_b, cx_b = centroid(b_sci, pixel_scales, start=(cy_r, cx_r))
    band_offset_px = (round(float(cy_b - cy_r), 3), round(float(cx_b - cx_r), 3))
    if np.hypot(*band_offset_px) > MAX_REGISTRATION_PX:
        print(f"  NOTE: {lens} {blue_f} sits {np.hypot(*band_offset_px):.2f}px from {red_f} by "
              f"deflector centroid -- check the tie against FIELD SOURCES, and if it is real "
              f"fix it in align_wfpc2_to_acs.py (this tool does not re-register).")

    # k: the DEFLECTOR's own blue/red flux ratio, keyed on the RED BAND'S OWN BRIGHTNESS
    # rather than on radius. An early-type is bluer outwards, so a single scalar k leaves that
    # colour gradient behind as a broad positive or negative disc lying exactly where the arcs
    # are -- k has to track the gradient. Binning by radius is the obvious way to do that and
    # it was tried first; it needs the isophote shape, and real ellipticals TWIST and grow
    # discs, so elliptical annuli of one fixed q and PA mis-scale the envelope and leave
    # signed lobes along the major axis (J0841+3824). Binning by brightness assumes NO shape
    # at all: the red image traces its own isophotes, whatever they do, and each brightness
    # bin's median blue/red ratio is the deflector's colour on that isophote. Arcs cannot bias
    # a bin's median -- a bin is a whole isophote's worth of pixels and the arc is a small part
    # of one. Measured against the radial version it cut the in-lens residual scatter on every
    # lens tried, and removed the bullseye rings a noisy per-annulus ratio printed round the
    # core.
    #
    # Bin assignment uses a SMOOTHED copy of the red band (0.15"), so a pixel's own noise
    # cannot scatter it into the wrong colour bin; the ratio itself is measured on the
    # unsmoothed images.
    galaxy = valid & ~(red['cmask'] | b_cmask) & (red_smooth > 2.0 * sky_red)
    k_map = colour_scale(red_m, blue_m, red_smooth, galaxy, sky_red)

    be = blue_m - k_map * red_m
    be_var = blue_var_m + k_map ** 2 * red_var_m

    # Empty sky for the background terms: outside any plausible arc, inside the frame, and
    # not a flagged pixel or a hand-masked contaminant.
    background = outskirts & valid & ~red['cmask'] & ~b_cmask
    if background.sum() < 500:                       # tiny stamp / heavily masked field
        background = (r_arcsec > max(2.0 * theta_e, 2.0)) & valid
    be = be - np.median(be[background])              # residual sky-term offset between bands

    sigma_px = smooth_arcsec / pixel_scales
    smoothed = gaussian_filter(be, sigma_px, mode='nearest')
    kern = np.zeros((25, 25))
    kern[12, 12] = 1.0
    kern = gaussian_filter(kern, sigma_px)
    den_smooth = np.sqrt(np.clip(fftconvolve(be_var, kern ** 2, mode='same'), 1e-30, None))
    den_raw = np.sqrt(np.clip(be_var, 1e-30, None))
    # Calibrate both maps on their own empty sky: propagated variance underestimates the real
    # scatter (drizzle + PSF matching correlate neighbouring pixels), typically by 3-4x on the
    # smoothed map. After this a blank-sky pixel reads S/N ~ 1, so a threshold is meaningful.
    alpha_smooth = _mad((smoothed / den_smooth)[background])
    alpha_raw = _mad((be / den_raw)[background])
    # CORE-DIPOLE DIAGNOSTIC, reported rather than corrected. A centring or PSF mismatch
    # leaves a residual that is ANTISYMMETRIC about the deflector -- positive one side,
    # negative the other -- while arcs and noise are not. Comparing `be` with its own
    # 180-degree rotation just inside the arc annulus gives a unitless number: ~1 or below is
    # a clean subtraction, well above 1 means the middle of that lens is model error and its
    # inner components deserve a harder look in the GUI.
    inner = (r_arcsec > 0.15 * theta_e) & (r_arcsec < 0.6 * theta_e) & valid
    rot = be[::-1, ::-1]
    core_dipole = round(float(np.median(np.abs(0.5 * (be - rot)[inner]))
                              / max(np.median(np.abs(0.5 * (be + rot)[inner])), 1e-30)), 2) \
        if inner.sum() > 50 else None
    return dict(
        lens=lens, red_filt=red_f, blue_filt=blue_f, pixel_scales=pixel_scales,
        theta_e=theta_e, cy=cy, cx=cx, be=be, k_map=k_map,
        snr=smoothed / den_smooth / alpha_smooth, snr_raw=be / den_raw / alpha_raw,
        red_matched=red_m, blue_matched=blue_m, valid=valid, r_arcsec=r_arcsec,
        cmask=red['cmask'] | b_cmask, hdr=red['hdr'],
        alpha_smooth=alpha_smooth, alpha_raw=alpha_raw, band_offset_px=band_offset_px,
        core_dipole=core_dipole,
        blue_snr_raw=blue_m / np.sqrt(np.clip(blue_var_m, 1e-30, None)) / max(
            _mad((blue_m / np.sqrt(np.clip(blue_var_m, 1e-30, None)))[background]), 1e-6),
    )


# ------------------------------------------------------------------------------- detection

def detect(maps, snr=DEFAULTS['snr'], peak=DEFAULTS['peak'], min_area=DEFAULTS['min_area'],
           r_in=DEFAULTS['r_in'], r_out=DEFAULTS['r_out'], grow=DEFAULTS['grow']):
    """Threshold the calibrated blue-excess S/N inside the theta_E annulus and keep the
    components that survive a per-pixel peak and a connected-area cut.

    Returns (region, components) with `region` a bool array on the red band's grid, True
    where the proposal says 'arc'. Per component it records the radius (in arcsec and in
    theta_E), the peak per-pixel S/N, the area, and how tangential it is -- the last as
    INFORMATION, not a cut: a knot or a demagnified counter-image is not elongated, and
    tangential elongation on its own is weak shear as often as it is an arc.
    """
    theta_e, ps = maps['theta_e'], maps['pixel_scales']
    annulus = ((maps['r_arcsec'] >= max(r_in * theta_e, 0.25))
               & (maps['r_arcsec'] <= r_out * theta_e))
    candidate = (maps['snr'] > snr) & annulus & maps['valid'] & ~maps['cmask']
    # 3x3 opening: a threshold crossing that does not survive the loss of its border was a
    # noise spike, not a resolved feature.
    candidate = binary_opening(candidate, np.ones((3, 3), dtype=bool))
    labels, n = label(candidate, structure=np.ones((3, 3)))
    keep = np.zeros_like(candidate)
    components = []
    for i in range(1, n + 1):
        sel = labels == i
        area = int(sel.sum())
        peak_snr = float(maps['snr_raw'][sel].max())
        if area < min_area or peak_snr < peak:
            continue
        # The blue excess can also go positive where the RED band is over-subtracted, which is
        # a model artefact and not a source. Require the feature to be there in the blue band
        # on its own too.
        if float(maps['blue_snr_raw'][sel].max()) < 2.0:
            continue
        ys, xs = np.nonzero(sel)
        # IS THERE A HOLE OPPOSITE? A centring or PSF mismatch leaves an ANTISYMMETRIC
        # residual: positive one side of the deflector, negative the other. Lensing does not
        # do that -- the mirror of an arc is either a counter-image (positive) or empty sky
        # (~0), never a deficit. So a component whose 180-degree mirror reads significantly
        # NEGATIVE is the bright half of a dipole, not a source, and is dropped. This is the
        # per-component form of the `core_dipole` diagnostic, and it is the cheapest guard
        # against the one artefact this method is prone to.
        my = np.clip(np.round(2 * maps['cy'] - ys).astype(int), 0, maps['snr'].shape[0] - 1)
        mx = np.clip(np.round(2 * maps['cx'] - xs).astype(int), 0, maps['snr'].shape[1] - 1)
        mirror_snr = float(np.median(maps['snr'][my, mx]))
        if mirror_snr < MIRROR_DEFICIT_SNR:
            continue
        w = np.clip(maps['snr'][sel], 0, None)
        rr = np.hypot(ys - maps['cy'], xs - maps['cx']) * ps
        r_mean = float((rr * w).sum() / w.sum())
        dx, dy = xs - maps['cx'], ys - maps['cy']
        dxm, dym = dx - dx.mean(), dy - dy.mean()
        m_xx = (w * dxm * dxm).sum() / w.sum()
        m_yy = (w * dym * dym).sum() / w.sum()
        m_xy = (w * dxm * dym).sum() / w.sum()
        major = 0.5 * np.arctan2(2 * m_xy, m_xx - m_yy)
        radial = np.arctan2(dy.mean(), dx.mean())
        offset = np.degrees(abs(((major - radial + np.pi / 2) % np.pi) - np.pi / 2))
        components.append(dict(
            area_px=area, peak_snr=round(peak_snr, 2),
            peak_snr_smoothed=round(float(maps['snr'][sel].max()), 2),
            r_arcsec=round(r_mean, 3), r_over_theta_e=round(r_mean / theta_e, 3),
            tangential_deg=round(float(90.0 - offset), 1),
            mirror_snr=round(mirror_snr, 2)))
        keep |= sel
    # Grow slightly and close: a threshold traces an arc's bright spine, and a mask wants its
    # wings too. Two pixels is 0.1" -- under the PSF, so it cannot invent extent.
    region = binary_fill_holes(binary_closing(
        binary_dilation(keep, np.ones((3, 3), dtype=bool), iterations=grow),
        np.ones((3, 3), dtype=bool)))
    region &= maps['valid'] & ~maps['cmask']
    components.sort(key=lambda c: -c['area_px'])
    return region, components


def propose(filt_dirs, prefixes, lens, theta_e=None, **kw):
    """Detect the arcs and return (region, components, maps) on the RED band's grid.

    This is the entry point `make_arc_masks.py --propose-from detect` calls; it reprojects the
    region onto whichever band is being drawn. Raises ValueError with a readable message when
    the lens cannot be run (no blue band, no measured theta_E).
    """
    if theta_e is None:
        theta_e = theta_e_table().get(lens)
    if theta_e is None:
        raise ValueError(f'{lens}: no measured Einstein radius in '
                         f'{os.path.relpath(THETA_E_JSON, ws_path)} '
                         f'(run --fetch-theta-e; not every sample is in Auger+2009)')
    det_kw = {k: kw.pop(k) for k in list(kw) if k in DEFAULTS and k != 'smooth'}
    maps = colour_excess(filt_dirs, prefixes, lens, theta_e,
                         smooth_arcsec=kw.pop('smooth', DEFAULTS['smooth']))
    region, components = detect(maps, **det_kw)
    return region, components, maps


# ----------------------------------------------------------------------------- QC figure

def qc_figure(maps, region, components, out_path):
    """Four panels: the two PSF-matched bands, the blue-excess image, and its S/N -- each with
    the proposed region outlined and the measured theta_E circled, so the proposal can be
    judged against the colour evidence it was made from."""
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    ps, theta_e = maps['pixel_scales'], maps['theta_e']
    ny, nx = maps['be'].shape
    half = nx * ps / 2.0
    extent = [-half, half, -half, half]
    cx_as = (maps['cx'] - (nx - 1) / 2.0) * ps
    cy_as = -(maps['cy'] - (ny - 1) / 2.0) * ps
    t = np.linspace(0, 2 * np.pi, 256)

    panels = [(maps['red_matched'], f"{maps['red_filt']} (PSF-matched)", None),
              (maps['blue_matched'], f"{maps['blue_filt']} (PSF-matched)", None),
              (maps['be'], f"blue excess: {maps['blue_filt']} - k x {maps['red_filt']}", None),
              (maps['snr'], 'blue-excess S/N (calibrated)', (-3, 8))]
    fig, axes = plt.subplots(1, 4, figsize=(20, 5.6))
    for ax, (arr, name, vlim) in zip(axes, panels):
        a = np.nan_to_num(np.asarray(arr, dtype=float))
        if vlim is None:
            lo, hi = np.percentile(a[maps['valid']], [5, 99.5])
            ax.imshow(np.arcsinh((a - lo) / max(hi - lo, 1e-12) / 0.1), cmap='jet',
                      extent=extent)
        else:
            ax.imshow(a, cmap='jet', extent=extent, vmin=vlim[0], vmax=vlim[1])
        if region.any():
            # origin='upper' is REQUIRED with an extent: imshow puts row 0 at the top and
            # contour puts it at the bottom, so without it the outline is drawn mirrored.
            ax.contour(region.astype(float), levels=[0.5], colors='white', linewidths=1.2,
                       extent=extent, origin='upper')
        ax.plot(cx_as + theta_e * np.cos(t), cy_as + theta_e * np.sin(t), 'w--', lw=0.7)
        ax.set_title(f"{maps['lens']}  {name}", fontsize=10)
        ax.set_xlabel('arcsec')
    off = maps['band_offset_px']
    summary = (f"theta_E={theta_e:.2f}\" (Auger+09)   "
               f"band offset {np.hypot(*off):.2f}px (measured, not applied)   "
               f"core dipole {maps['core_dipole']}   "
               f"{len(components)} component(s), {int(region.sum())}px proposed   "
               f"noise calib x{maps['alpha_smooth']:.1f}")
    if components:
        summary += '   |   ' + '; '.join(
            f"r={c['r_arcsec']:.2f}\"={c['r_over_theta_e']:.2f}th_E "
            f"A={c['area_px']}px pk={c['peak_snr']:.1f} tang={c['tangential_deg']:.0f}deg"
            for c in components[:4])
    fig.suptitle(summary, fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=95)
    plt.close(fig)


def contact_sheet(thumbs, out_path, columns=6):
    """One page of every lens's blue-excess S/N with its proposal outlined.

    The per-lens four-panel figures are where a single proposal is judged; this is where the
    SAMPLE is judged -- whether the detector is finding arcs where arcs are, at a glance,
    before anyone opens a GUI.
    """
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    rows = int(np.ceil(len(thumbs) / columns))
    fig, axes = plt.subplots(rows, columns, figsize=(2.6 * columns, 2.85 * rows))
    for ax, thumb in zip(np.ravel(axes), thumbs):
        snr, region, theta_e, ps, cy_as, cx_as, lens, n_comp = thumb
        half = snr.shape[1] * ps / 2.0
        extent = [-half, half, -half, half]
        ax.imshow(snr, cmap='jet', extent=extent, vmin=-3, vmax=8)
        if region.any():
            ax.contour(region.astype(float), levels=[0.5], colors='white', linewidths=0.9,
                       extent=extent, origin='upper')
        t = np.linspace(0, 2 * np.pi, 200)
        ax.plot(cx_as + theta_e * np.cos(t), cy_as + theta_e * np.sin(t), 'w--', lw=0.5)
        ax.set_xlim(-3.5, 3.5)
        ax.set_ylim(-3.5, 3.5)
        ax.set_title(f'{lens}  {n_comp}c {int(region.sum())}px', fontsize=7)
        ax.set_xticks([])
        ax.set_yticks([])
    for ax in np.ravel(axes)[len(thumbs):]:
        ax.axis('off')
    fig.suptitle('Automatic arc proposals: blue-excess S/N, proposed region outlined, '
                 'dashed circle = measured theta_E', fontsize=10)
    fig.tight_layout()
    fig.savefig(out_path, dpi=110)
    plt.close(fig)


# --------------------------------------------------------------------------------- CLI

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--sample', default=mast_target_names.DEFAULT_SAMPLE,
                   help=f'sample under data/cutouts/ (default {mast_target_names.DEFAULT_SAMPLE})')
    p.add_argument('--lens', default=None, help='restrict to one lens')
    p.add_argument('--size', type=float, default=cutout_paths.DEFAULT_SIZE,
                   help=f'cutout tree (default {cutout_paths.DEFAULT_SIZE:g}")')
    p.add_argument('--pass', dest='drizzle_pass', choices=['auto', 'cr', 'nocrrej'],
                   default='auto', help="cutout pass, as make_masks.py (default auto)")
    p.add_argument('--fetch-theta-e', action='store_true',
                   help='refresh info/lens_einstein_radii.json from VizieR and exit')
    p.add_argument('--out-dir', default=None,
                   help='where the QC figures go (default diagnostics/arc_detection/<sample>/)')
    p.add_argument('--no-figures', action='store_true', help='summary JSON only')
    for name, val in DEFAULTS.items():
        p.add_argument(f'--{name.replace("_", "-")}', type=type(val), default=val,
                       help=f'detection parameter (default {val})')
    a = p.parse_args()

    import make_masks                              # discovery + prefix helpers, one source
    root = cutout_paths.cutouts_root(ws_path, a.size)
    lens_filts, lens_prefixes = {}, {}
    for lens, filt, cutout_dir in make_masks.discover_targets(root, a.sample, a.lens, None):
        prefix = make_masks.find_prefix(cutout_dir, a.drizzle_pass)
        if prefix is None:
            continue
        lens_filts.setdefault(lens, {})[filt] = cutout_dir
        lens_prefixes.setdefault(lens, {})[filt] = prefix

    if a.fetch_theta_e:
        fetch_theta_e(sorted(lens_filts))
        return
    if not lens_filts:
        raise SystemExit(f'no cutouts under {root} for sample {a.sample}')

    out_dir = a.out_dir or os.path.join(ws_path, 'diagnostics', 'arc_detection', a.sample)
    os.makedirs(out_dir, exist_ok=True)
    kw = {k: getattr(a, k) for k in DEFAULTS}
    summary_path = os.path.join(out_dir, 'arc_detection_summary.json')
    # A --lens run updates that lens's entry and leaves the rest of the sample's alone: this
    # file is the record of a sweep, and re-tuning one lens must not erase the other 37.
    results = {}
    if os.path.exists(summary_path):
        with open(summary_path) as fh:
            results = json.load(fh)
    thumbs = []
    for lens in sorted(lens_filts):
        try:
            region, components, maps = propose(lens_filts[lens], lens_prefixes[lens], lens, **kw)
        except Exception as exc:
            print(f'{lens}: SKIPPED -- {exc}')
            results[lens] = {'error': str(exc)}
            continue
        n_px = int(region.sum())
        entry = dict(red_filt=maps['red_filt'], blue_filt=maps['blue_filt'],
                     theta_e_arcsec=round(maps['theta_e'], 3), n_proposed_px=n_px,
                     n_components=len(components),
                     area_arcsec2=round(n_px * maps['pixel_scales'] ** 2, 4),
                     noise_calibration=round(maps['alpha_smooth'], 3),
                     band_offset_px=list(maps['band_offset_px']),
                     core_dipole=maps['core_dipole'],
                     components=components, params=kw)
        results[lens] = entry
        if not a.no_figures:
            qc_figure(maps, region, components,
                      os.path.join(out_dir, f'{lens}_arc_detection.png'))
            ny, nx = maps['be'].shape
            thumbs.append((maps['snr'], region, maps['theta_e'], maps['pixel_scales'],
                           -(maps['cy'] - (ny - 1) / 2.0) * maps['pixel_scales'],
                           (maps['cx'] - (nx - 1) / 2.0) * maps['pixel_scales'],
                           lens, len(components)))
        print(f"{lens}: {maps['blue_filt']}-{maps['red_filt']}  theta_E="
              f"{maps['theta_e']:.2f}\"  {len(components)} comp, {n_px}px  "
              + '; '.join(f"r={c['r_over_theta_e']:.2f}th_E pk={c['peak_snr']:.0f} "
                          f"A={c['area_px']}" for c in components[:3]))

    with open(summary_path, 'w') as fh:
        json.dump(dict(sorted(results.items())), fh, indent=1)
        fh.write('\n')
    if len(thumbs) > 1:
        contact_sheet(thumbs, os.path.join(out_dir, 'contact_sheet.png'))
        print(f"contact sheet -> {os.path.relpath(out_dir, ws_path)}/contact_sheet.png")
    found = [l for l, r in results.items() if r.get('n_components')]
    print(f"\n{len(found)}/{len(results)} lenses with a detection; summary -> "
          f"{os.path.relpath(summary_path, ws_path)}")
    if not a.no_figures:
        print(f"QC figures -> {os.path.relpath(out_dir, ws_path)}/")


if __name__ == '__main__':
    main()
