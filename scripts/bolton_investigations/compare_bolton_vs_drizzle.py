#!/usr/bin/env python
"""Quantify how the Bolton-style bilinear science image differs from our drizzled one.

STANDALONE comparison -- NOT part of the pipeline. Both products are 0.05"/px North-up on
J1023+4230 F814W, but centred slightly differently, on different background pedestals, and
(possibly) different flux normalisation, so a fair comparison must first register + match
them. We solve for shift + flux-scale + background offset by least-squares on the galaxy,
then report the residual and the PSF sharpness (the expected real difference: bilinear
resampling is softer than drizzle).

Outputs (output/bolton_investigations/):
  bolton_vs_drizzle.png   - images, difference, radial profiles, star-PSF cuts
  (prints the metrics table)
"""
import os
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.modeling import models, fitting
from scipy.ndimage import shift as nd_shift
from scipy.optimize import minimize
from photutils.detection import DAOStarFinder
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.visualization import AsinhStretch, ImageNormalize, PercentileInterval

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = f"{REPO}/output/bolton_investigations"
os.makedirs(OUT, exist_ok=True)

drz = fits.getdata(f"{REPO}/data/cutouts/slacs_gold/J1023+4230/f814W/cutout_cr_sci.fits").astype(float)
bol_full = fits.getdata(f"{OUT}/bolton_J1023_sci.fits").astype(float)
H = drz.shape[0]
c0 = (bol_full.shape[0] - H) // 2
bol = bol_full[c0:c0 + H, c0:c0 + H]                       # central 20" == the cutout FOV
yy, xx = np.mgrid[0:H, 0:H]
r = np.hypot(xx - H / 2, yy - H / 2)

# ── background: sigma-clipped sky from the outer region, subtract from each ──────────
sky_drz = sigma_clipped_stats(drz[r > 150], sigma=3)[1]
sky_bol = sigma_clipped_stats(bol[r > 150], sigma=3)[1]
drz -= sky_drz
bol -= sky_bol
bg_rms_drz = sigma_clipped_stats(drz[r > 150], sigma=3)[2]
bg_rms_bol = sigma_clipped_stats(bol[r > 150], sigma=3)[2]

# ── register + flux-scale + residual background on the galaxy (r<120) ───────────────
fitreg = r < 120


def cost(p):
    dx, dy, s, b = p
    m = nd_shift(bol, [dy, dx], order=3, mode="nearest")
    return np.nansum(((drz - (s * m + b))[fitreg]) ** 2)


res = minimize(cost, [0, 0, 1.0, 0.0], method="Powell")
dx, dy, scale, bkg = res.x
bol_reg = scale * nd_shift(bol, [dy, dx], order=3, mode="nearest") + bkg
diff = drz - bol_reg
print(f"registration  dx={dx:+.2f} dy={dy:+.2f} px   flux scale (bol->drz) = {scale:.4f}")

# ── PSF sharpness: fit a Gaussian to the brightest isolated star in each ─────────────
finder = DAOStarFinder(fwhm=3.0, threshold=8 * bg_rms_drz)
srcs = finder(drz)
srcs = srcs[(np.hypot(srcs["xcentroid"] - H / 2, srcs["ycentroid"] - H / 2) > 60)
            & (srcs["xcentroid"] > 20) & (srcs["xcentroid"] < H - 20)
            & (srcs["ycentroid"] > 20) & (srcs["ycentroid"] < H - 20)]
srcs.sort("flux")
sx, sy = float(srcs[-1]["xcentroid"]), float(srcs[-1]["ycentroid"])


def fwhm_at(img, x, y, half=7):
    xi, yi = int(round(x)), int(round(y))
    stamp = img[yi - half:yi + half + 1, xi - half:xi + half + 1]
    gy, gx = np.mgrid[0:stamp.shape[0], 0:stamp.shape[1]]
    g0 = models.Gaussian2D(amplitude=stamp.max(), x_mean=half, y_mean=half,
                           x_stddev=2, y_stddev=2)
    gf = fitting.LevMarLSQFitter()(g0, gx, gy, stamp)
    s = 0.5 * (abs(gf.x_stddev.value) + abs(gf.y_stddev.value))
    return 2.3548 * s, gf, stamp


fw_drz, gd, st_drz = fwhm_at(drz, sx, sy)
fw_bol, gb, st_bol = fwhm_at(bol_reg, sx, sy)
print(f"stellar FWHM  drizzle {fw_drz:.3f} px ({fw_drz*0.05:.3f}\")   "
      f"bilinear {fw_bol:.3f} px ({fw_bol*0.05:.3f}\")   "
      f"bilinear is {100*(fw_bol/fw_drz-1):+.1f}% broader")

# ── residual level and photometry ───────────────────────────────────────────────────
core = r < 60
peak = drz[core].max()
resid_rms_core = np.sqrt(np.nanmean(diff[core] ** 2))
flux_drz = drz[core].sum()
flux_bol = bol_reg[core].sum()
print(f"core residual RMS = {resid_rms_core:.4g} e/s = {100*resid_rms_core/peak:.2f}% of peak")
print(f"aperture flux (r<60) drizzle {flux_drz:.1f}  bilinear {flux_bol:.1f}  "
      f"ratio {flux_bol/flux_drz:.4f}")
print(f"background RMS  drizzle {bg_rms_drz:.4g}   bilinear {bg_rms_bol:.4g} e/s  "
      f"(bilinear/drizzle {bg_rms_bol/bg_rms_drz:.2f})")

# ── azimuthally-averaged radial profile of the deflector ────────────────────────────
rb = np.arange(0, 120, 1.0)
prof_d = np.array([np.nanmean(drz[(r >= a) & (r < a + 1)]) for a in rb])
prof_b = np.array([np.nanmean(bol_reg[(r >= a) & (r < a + 1)]) for a in rb])

# ── figure ───────────────────────────────────────────────────────────────────────
def sn(a):
    return ImageNormalize(drz, interval=PercentileInterval(99.4), stretch=AsinhStretch())

fig = plt.figure(figsize=(16.5, 10))
gs = fig.add_gridspec(2, 3)
a00 = fig.add_subplot(gs[0, 0]); a00.imshow(drz, norm=sn(drz), cmap="inferno", origin="lower")
a00.set_title("Drizzle science"); a00.plot(sx, sy, "co", mfc="none", ms=12)
a01 = fig.add_subplot(gs[0, 1]); a01.imshow(bol_reg, norm=sn(drz), cmap="inferno", origin="lower")
a01.set_title("Bolton bilinear (registered, flux-matched)")
dl = np.nanpercentile(np.abs(diff[r < 120]), 99)
a02 = fig.add_subplot(gs[0, 2]); im = a02.imshow(diff, vmin=-dl, vmax=dl, cmap="RdBu_r", origin="lower")
a02.set_title("Difference (drizzle − bilinear)"); fig.colorbar(im, ax=a02, fraction=0.046)
for a in (a00, a01, a02):
    a.set_xticks([]); a.set_yticks([])

a10 = fig.add_subplot(gs[1, 0])
a10.plot(rb * 0.05, prof_d, label="drizzle", lw=2)
a10.plot(rb * 0.05, prof_b, label="bilinear", lw=2, ls="--")
a10.set_yscale("symlog", linthresh=1e-3); a10.set_xlabel("radius (arcsec)")
a10.set_ylabel("mean flux (e/s)"); a10.legend(); a10.set_title("Deflector radial profile")

a11 = fig.add_subplot(gs[1, 1])
mid = st_drz.shape[0] // 2
a11.plot(st_drz[mid], "o-", label=f"drizzle  FWHM {fw_drz:.2f}px")
a11.plot(st_bol[mid], "s--", label=f"bilinear FWHM {fw_bol:.2f}px")
a11.set_title("Star PSF — central row cut"); a11.legend(); a11.set_xlabel("px")

a12 = fig.add_subplot(gs[1, 2]); a12.axis("off")
txt = (f"J1023+4230 F814W  (0.05\"/px, 20\" stamp)\n\n"
       f"registration shift: dx={dx:+.2f} dy={dy:+.2f} px\n"
       f"flux scale (bilinear→drizzle): {scale:.3f}\n\n"
       f"stellar FWHM:\n  drizzle  {fw_drz:.2f} px ({fw_drz*0.05:.3f}\")\n"
       f"  bilinear {fw_bol:.2f} px ({fw_bol*0.05:.3f}\")\n"
       f"  bilinear {100*(fw_bol/fw_drz-1):+.1f}% broader\n\n"
       f"core residual RMS: {100*resid_rms_core/peak:.2f}% of peak\n"
       f"aperture flux ratio (r<3\"): {flux_bol/flux_drz:.3f}\n"
       f"background RMS ratio: {bg_rms_bol/bg_rms_drz:.2f}×")
a12.text(0.0, 0.95, txt, va="top", ha="left", family="monospace", fontsize=12)
fig.suptitle("Bolton bilinear vs drizzle science — after registration + flux match",
             fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig(f"{OUT}/bolton_vs_drizzle.png", dpi=130)
print("wrote bolton_vs_drizzle.png")
