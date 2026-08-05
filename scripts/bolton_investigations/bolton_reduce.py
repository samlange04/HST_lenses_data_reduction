#!/usr/bin/env python
"""Bolton et al. 2008 (SLACS V) style bilinear-interpolation reduction.

STANDALONE comparison script -- NOT part of the production pipeline. It reproduces the
methodological choices of the *original* SLACS imaging pipeline (Bolton et al. 2008,
arXiv:0805.1931) that make its noise maps free of the dead-column diagonal stripes this
repo's AstroDrizzle reduction correctly shows (see the `legacy-slacs-bolton-bilinear-no-
stripes` memory). The point is to demonstrate, on real data, *why* a legacy noise map on
the same lens looks clean -- it is a property of the reduction method, not a bug here.

Recipe (faithful to the paper, deliberately minimal):
  1. Start from calibrated ACS/WFC FLC frames (electrons, full distortion model).
  2. Mask cosmic rays (LACosmic / astroscrappy) + strongly-negative cold pixels ONLY.
     NO DQ bad-column / warm / hot masking -> no per-pixel weight deficit is ever made.
  3. Rectify each frame onto a common North-up 0.05"/px grid by BILINEAR INTERPOLATION
     (stwcs full-distortion WCS -> scipy map_coordinates order=1), NOT drizzle.
  4. Combine frames with a nan-aware sigma-clipped mean -- the dither+combine that
     removes residual defects without ever building a weight map.
  5. Build the noise map from the COMBINED COUNTS + a measured background RMS. There is
     no weight map, so a masked-column weight deficit cannot appear in it.

Applied to J1023+4230 F814W, whose bad detector column runs through the deflector core,
so the drizzle stripe is maximally visible for the side-by-side.

Outputs (all under output/bolton_investigations/):
  bolton_J1023_sci.fits    - Bolton-style science image (e/s)
  bolton_J1023_noise.fits  - Bolton-style noise map (e/s)
  bolton_J1023_compare.png - 4-panel comparison vs this repo's drizzle products
"""
import glob
import os
import sys
import warnings

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import sigma_clip, sigma_clipped_stats
from astropy.visualization import AsinhStretch, ImageNormalize, PercentileInterval
from scipy.ndimage import map_coordinates
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from stwcs.wcsutil import HSTWCS
import astroscrappy

# ── config ────────────────────────────────────────────────────────────────────
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
LENS = "J1023+4230"
FILT = "f814W"
CAL_DIR = f"{REPO}/data/calibrated/slacs_gold/{LENS}/{FILT}"
OUT_DIR = f"{REPO}/output/bolton_investigations"
os.makedirs(OUT_DIR, exist_ok=True)
SCALE = 0.05 / 3600.0          # deg/px, matches the ACS/WFC drizzle output scale
NPIX = 800                     # 40" output grid (central 20" == the repo's cutout)
READNOISE = 5.2                # e-, mean over the four WFC amps
SATLEVEL = 84000.0             # e-, WFC full-well
SIGCLIP, OBJLIM = 4.5, 5.0     # same LACosmic params the repo uses

# ── lens position (same source the drizzle scripts use) ─────────────────────────
sys.path.insert(0, f"{REPO}/info")
from slacs_coords import slacs_coords
from astropy.coordinates import SkyCoord
import astropy.units as u
lc = SkyCoord(*slacs_coords[LENS], unit=(u.hourangle, u.deg))
RA0, DEC0 = float(lc.ra.deg), float(lc.dec.deg)

# ── common North-up output WCS (0.05"/px, tangent at the lens) ───────────────────
out_wcs = WCS(naxis=2)
out_wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
out_wcs.wcs.crpix = [NPIX / 2 + 0.5, NPIX / 2 + 0.5]
out_wcs.wcs.crval = [RA0, DEC0]
out_wcs.wcs.cd = [[-SCALE, 0.0], [0.0, SCALE]]     # North-up, East-left

# output pixel grid -> sky (output WCS carries no distortion)
yy, xx = np.mgrid[0:NPIX, 0:NPIX]
ra, dec = out_wcs.all_pix2world(xx.ravel(), yy.ravel(), 0)

# ── rectify every FLC chip that overlaps the window ─────────────────────────────
frames = sorted(glob.glob(f"{CAL_DIR}/*flc.fits"))
print(f"{len(frames)} FLC frames for {LENS} {FILT}")
rect_stack, exptimes = [], []

for fp in frames:
    with fits.open(fp) as hdul:
        exptime = float(hdul[0].header["EXPTIME"])
        for ext in range(1, len(hdul)):
            hdr = hdul[ext].header
            if hdr.get("EXTNAME") != "SCI":
                continue
            wcs_in = HSTWCS(hdul, ext=ext)          # full SIP+NPOL+D2IM distortion
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                x_in, y_in = wcs_in.all_world2pix(ra, dec, 0)
            ny_in, nx_in = hdul[ext].data.shape
            inb = (x_in > 1) & (x_in < nx_in - 2) & (y_in > 1) & (y_in < ny_in - 2)
            if inb.sum() == 0:
                continue                            # this chip doesn't touch the window

            sci_e = hdul[ext].data.astype(np.float32)   # electrons (total)
            # (2) mask CRs (LACosmic) + strongly-negative cold pixels ONLY; no DQ masking
            crmask, _ = astroscrappy.detect_cosmics(
                sci_e, gain=1.0, readnoise=READNOISE, satlevel=SATLEVEL,
                sigclip=SIGCLIP, objlim=OBJLIM, niter=4, cleantype="medmask")
            cold = sci_e < -5.0 * READNOISE
            rate = (sci_e / exptime).astype(np.float32)     # electrons / s
            rate[crmask | cold] = np.nan

            # (3) bilinear rectification onto the output grid (order=1 == bilinear)
            samp = map_coordinates(
                rate, [y_in, x_in], order=1, mode="constant",
                cval=np.nan, prefilter=False)
            out = np.full(NPIX * NPIX, np.nan, np.float32)
            out[inb] = samp[inb]
            rect_stack.append(out.reshape(NPIX, NPIX))
            exptimes.append(exptime)
            print(f"  rectified {os.path.basename(fp)} ext {ext} "
                  f"({inb.sum()} px in window, {np.isnan(rate).mean()*100:.2f}% masked)")

rect = np.stack(rect_stack)                          # (Nchipframes, NPIX, NPIX)

# ── (4) combine: nan-aware sigma-clipped mean (dither+combine, no weight map) ────
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    clipped = sigma_clip(rect, sigma=3.0, maxiters=3, axis=0, masked=True)
    sci = np.ma.mean(clipped, axis=0).filled(np.nan).astype(np.float32)
    ncov = (~clipped.mask & np.isfinite(rect)).sum(axis=0).astype(np.float32)

# ── (5) noise map FROM THE COUNTS (no weight map exists) ─────────────────────────
# background RMS in rate units, measured from the combined image itself
_, _, bg_std = sigma_clipped_stats(sci[np.isfinite(sci)], sigma=3.0)
n_ref = np.nanmax(ncov)
t_eff = np.where(ncov > 0, ncov * np.median(exptimes), np.nan)   # per-pixel exposure
bg_term = (bg_std ** 2) * (n_ref / np.where(ncov > 0, ncov, np.nan))  # coverage-scaled
src_term = np.clip(sci, 0, None) / t_eff                          # source Poisson (rate)
noise = np.sqrt(bg_term + src_term).astype(np.float32)
sci_filled = np.where(np.isfinite(sci), sci, 0.0).astype(np.float32)
noise = np.where(np.isfinite(noise), noise, np.nanmax(noise[np.isfinite(noise)]))

# ── write FITS ──────────────────────────────────────────────────────────────────
hdr = out_wcs.to_header()
hdr["BUNIT"] = "ELECTRONS/S"
hdr["REDUCTN"] = ("bolton2008-bilinear", "reduction method (comparison, not pipeline)")
hdr["NFRAMES"] = (len(rect_stack), "rectified chip-frames combined")
hdr["BGRMS"] = (float(bg_std), "measured background RMS (e/s)")
fits.writeto(f"{OUT_DIR}/bolton_{LENS[:5]}_sci.fits", sci_filled, hdr, overwrite=True)
fits.writeto(f"{OUT_DIR}/bolton_{LENS[:5]}_noise.fits", noise, hdr, overwrite=True)
print(f"wrote bolton_{LENS[:5]}_sci.fits / _noise.fits  (bg_std={bg_std:.4g} e/s)")

# ── comparison figure vs this repo's drizzle cutouts ─────────────────────────────
cut = f"{REPO}/data/cutouts/slacs_gold/{LENS}/{FILT}"
drz_sci = fits.getdata(f"{cut}/cutout_cr_sci.fits")
drz_noise = fits.getdata(f"{cut}/cutout_cr_noise.fits")
# central 20" (400 px) of the 40" Bolton grid == same footprint as the repo cutout
c0 = (NPIX - drz_sci.shape[0]) // 2
c1 = c0 + drz_sci.shape[0]
bol_sci = sci_filled[c0:c1, c0:c1]
bol_noise = noise[c0:c1, c0:c1]


def sci_norm(a):
    return ImageNormalize(a, interval=PercentileInterval(99.3), stretch=AsinhStretch())


def noise_norm(a):
    lo, hi = np.nanpercentile(a, [2, 98])
    return dict(vmin=lo, vmax=hi)


fig, ax = plt.subplots(2, 2, figsize=(11, 11))
ax[0, 0].imshow(drz_sci, origin="lower", cmap="inferno", norm=sci_norm(drz_sci))
ax[0, 0].set_title("This repo: AstroDrizzle science (F814W)")
ax[0, 1].imshow(drz_noise, origin="lower", cmap="inferno", **noise_norm(drz_noise))
ax[0, 1].set_title("This repo: drizzle NOISE  — dead-column stripe visible")
ax[1, 0].imshow(bol_sci, origin="lower", cmap="inferno", norm=sci_norm(bol_sci))
ax[1, 0].set_title("Bolton bilinear science (this script)")
ax[1, 1].imshow(bol_noise, origin="lower", cmap="inferno", **noise_norm(bol_noise))
ax[1, 1].set_title("Bolton bilinear NOISE — no stripe (count-derived, no weight map)")
for a in ax.ravel():
    a.set_xticks([]); a.set_yticks([])
fig.suptitle(f"{LENS} F814W — drizzle vs Bolton-2008 bilinear reduction "
             f"(central 20\", 0.05\"/px)", fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.98))
fig.savefig(f"{OUT_DIR}/bolton_{LENS[:5]}_compare.png", dpi=130)
print(f"wrote bolton_{LENS[:5]}_compare.png")
