#!/usr/bin/env python
"""Hybrid: drizzle image + bilinear-style healing of the dead-column stripe *only*
where it crosses the lens.

STANDALONE comparison script -- NOT part of the pipeline.

Motivation: in the AstroDrizzle product the dead-column artifact lives almost entirely
in the NOISE map (= 1/sqrt(WHT)). The science image is already stripe-free, because
inverse-variance weighting fills each striped output pixel from the good, dithered frames.
So to get "drizzle image quality + no stripe through the lens" we don't re-reduce anything:
we keep the drizzle science untouched and interpolate the noise map *across* the stripe in
the lens region only -- i.e. locally replace the weight-deficit pixels with the smooth,
source-structure-preserving local level (a 2-D interpolation across the thin diagonal).

Detection is from the noise map itself: the stripe is a thin diagonal of elevated noise on
top of a smoothly varying (source-dependent) noise field, so it stands out as a local
excess over a median-filtered version that bridges the ~1-3 px stripe but preserves the
core.

IMPORTANT caveat (printed at the end too): healing the noise map *down* asserts those
pixels are as good as their neighbours, when they really had one fewer contributing frame.
That is fine for a clean image/likelihood through the arc (it reproduces the legacy look on
a drizzle image) but it mildly OVER-weights genuinely-noisier pixels. The opposite, more
conservative modelling choice is to inflate/mask them (Etherington scalable-noise). This
script does the heal-down the question asks for; `--mask-up` flips it to the conservative
version for comparison.

Applied to J1023+4230 F814W (bad column through the deflector core).

Outputs (diagnostics/bolton_test_outputs/ — tracked on this branch while testing):
  hybrid_J1023_noise.fits    - drizzle noise map, stripe healed in the lens region
  hybrid_J1023_stripe.fits   - the detected stripe mask (1 = healed)
  hybrid_J1023_compare.png   - 3x3 (rows: drizzle / healed / difference;
                               cols: signal / noise / S/N; stripe detection overlaid on noise)
"""
import argparse
import glob
import os
import sys
import warnings
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import median_filter, binary_dilation, binary_opening, rotate
from astropy.visualization import AsinhStretch, ImageNormalize, PercentileInterval
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CUT = f"{REPO}/data/cutouts/slacs_gold/J1023+4230/f814W"
CAL = f"{REPO}/data/calibrated/slacs_gold/J1023+4230/f814W"
OUT = f"{REPO}/diagnostics/bolton_test_outputs"
os.makedirs(OUT, exist_ok=True)

ap = argparse.ArgumentParser()
ap.add_argument("--excess", type=float, default=0.05, help="min along-line excess over bg")
ap.add_argument("--linelen", type=int, default=41, help="along-stripe median length (px)")
ap.add_argument("--mfsize", type=int, default=9, help="across-stripe interpolation footprint (px)")
ap.add_argument("--mask-up", action="store_true", help="inflate noise instead of healing down")
ap.add_argument("--inflate", type=float, default=1e8, help="sentinel noise for --mask-up")
args = ap.parse_args()

sci = fits.getdata(f"{CUT}/cutout_cr_sci.fits").astype(np.float64)
noise = fits.getdata(f"{CUT}/cutout_cr_noise.fits").astype(np.float64)
cut_wcs = WCS(fits.getheader(f"{CUT}/cutout_cr_sci.fits"))

# ── seed angle: PA of the ACS detector-Y axis (bad columns run along Y), from the
#    FLC full-distortion WCS -> stamp pixel coords. Used only to seed the search. ─────
from stwcs.wcsutil import HSTWCS
flc = sorted(glob.glob(f"{CAL}/*flc.fits"))[0]
with fits.open(flc) as h:
    ext = next(e for e in range(1, len(h)) if h[e].header.get("EXTNAME") == "SCI"
               and h[e].header.get("CCDCHIP") == 1)          # lens is on chip 1
    wflc = HSTWCS(h, ext=ext)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        (r0, d0) = wflc.all_pix2world([[2048, 1024]], 0)[0]
        (r1, d1) = wflc.all_pix2world([[2048, 1124]], 0)[0]
x0s, y0s = cut_wcs.all_world2pix(r0, d0, 0)
x1s, y1s = cut_wcs.all_world2pix(r1, d1, 0)
theta_seed = np.degrees(np.arctan2(y1s - y0s, x1s - x0s)) % 180
print(f"detector-Y PA in stamp = {theta_seed:.1f} deg (search seed)")

# ── Radon-style detection on a HIGH-PASS noise map: subtract a large median filter to
#    kill the core wings + any large-scale gradient, leaving only thin features. Binning
#    the residual by perpendicular offset makes the stripe a clean narrow peak while
#    compact sources average out (they don't align along a whole offset line). ─────────
H, W = noise.shape
yy, xx = np.mgrid[0:H, 0:W].astype(float)
cy, cx = (H - 1) / 2, (W - 1) / 2
r = np.hypot(xx - cx, yy - cy)
hp = noise - median_filter(noise, size=21)                  # thin-feature residual
# clip out the BRIGHTEST thin features (satellite/CR trails also make weight-deficit lines,
# but far brighter than a bad column) so they can't dominate the offset profile; the faint
# bad-column stripe survives the clip.
hi_clip = np.nanpercentile(hp[r > 55], 95)
valid = (r > 55) & (hp < hi_clip)


def profile(theta_deg):
    t = np.radians(theta_deg)
    off = -(xx - cx) * np.sin(t) + (yy - cy) * np.cos(t)     # perpendicular offset (px)
    b = np.round(off - off.min()).astype(int)
    n = b.max() + 1
    s = np.bincount(b[valid], weights=hp[valid], minlength=n)
    c = np.bincount(b[valid], minlength=n).astype(float)
    prof = np.where(c > 30, s / np.maximum(c, 1), np.nan)
    sd = np.nanstd(prof)
    return off, b, prof, (np.nanmax(prof) / sd if sd > 0 else 0)


# The bad-column stripe angle is FIXED by the optics (detector-Y projected North-up), so
# use the exact WCS value -- no search. A search only lets a brighter, nearly-parallel
# field feature (the satellite trail here, ~39 deg) capture the fit.
theta = theta_seed
off, b, prof, score = profile(theta)
print(f"angle {theta:.2f} deg (fixed, geometric)  high-pass peak/sigma = {score:.1f}")

# stripe = offset bins whose high-pass residual sits >4 sigma above zero. Each is one
# full-length line along the stripe direction -> selecting the bins selects the stripe.
sd = np.nanstd(prof)
stripe_bins = np.where(prof > 4.0 * sd)[0]
in_bin = np.isin(b, stripe_bins)
smooth = median_filter(noise, size=args.mfsize)             # across-stripe interp value
stripe = binary_dilation(in_bin, iterations=1)             # thicken to fully cover the line
print(f"detected {stripe.sum()} stripe pixels in {len(stripe_bins)} offset bins "
      f"({100*stripe.mean():.2f}% of the stamp)")

# ── heal: interpolate the noise across the stripe (replace with the smooth local level) ──
healed = noise.copy()
if args.mask_up:
    healed[stripe] = args.inflate
    tag = "noise INFLATED on stripe (conservative / mask-up)"
else:
    healed[stripe] = smooth[stripe]
    tag = "noise HEALED across stripe (interpolated down)"

# ── write products ───────────────────────────────────────────────────────────────
hdr = fits.getheader(f"{CUT}/cutout_cr_noise.fits")
hdr["STRIPEHL"] = (not args.mask_up, "dead-column stripe healed in lens region")
fits.writeto(f"{OUT}/hybrid_J1023_noise.fits", healed.astype(np.float32), hdr, overwrite=True)
fits.writeto(f"{OUT}/hybrid_J1023_stripe.fits", stripe.astype(np.uint8), hdr, overwrite=True)

# ── figure: signal / noise / S/N for drizzle (before) vs healed (after) + differences ──
def sig_norm(a):
    return ImageNormalize(a, interval=PercentileInterval(99.3), stretch=AsinhStretch())


def noise_lim(a):
    lo, hi = np.nanpercentile(noise, [2, 98])          # shared scale from the drizzle noise
    return dict(vmin=lo, vmax=hi, cmap="inferno", origin="lower")


def snr_norm(a):
    return ImageNormalize(a, interval=PercentileInterval(99.0), stretch=AsinhStretch())


def snr_map(s, n):
    with np.errstate(divide="ignore", invalid="ignore"):
        r = s / n
    r[~np.isfinite(r)] = 0.0
    return r


def diff_lim(d):
    v = d[np.isfinite(d)]
    v = v[np.abs(v) > 0]                                # diffs here are sparse (stripe only)
    m = float(np.nanpercentile(np.abs(v), 99)) if v.size else 1e-12
    return dict(vmin=-(m or 1e-12), vmax=(m or 1e-12), cmap="RdBu_r", origin="lower")


snr_before = snr_map(sci, noise)
if args.mask_up:                                       # blank the down-weighted pixels
    noise_disp = np.where(stripe, np.nan, healed)
    snr_after = np.where(stripe, np.nan, snr_map(sci, healed))
    d_noise = np.where(stripe, np.nanpercentile(noise, 99) * 5, 0.0)   # sentinel-capped
else:
    noise_disp = healed
    snr_after = snr_map(sci, healed)
    d_noise = noise - healed                           # removed excess (positive on stripe)
d_snr = snr_after - snr_before                         # S/N change from the heal
sn, nl, rn = sig_norm(sci), noise_lim(noise), snr_norm(snr_before)

fig, ax = plt.subplots(3, 3, figsize=(13.5, 13.5))
# row 0 — drizzle (before), stripe detection overlaid on the noise panel
ax[0, 0].imshow(sci, cmap="inferno", origin="lower", norm=sn)
ax[0, 0].set_title("Drizzle: signal (kept, untouched)")
ax[0, 1].imshow(noise, **nl)
ov = np.ma.masked_where(~stripe, stripe)
ax[0, 1].imshow(ov, origin="lower", cmap="cool", vmin=0, vmax=1, alpha=0.9)
ax[0, 1].set_title("Drizzle: noise + detected stripe (magenta)")
ax[0, 2].imshow(snr_before, cmap="inferno", origin="lower", norm=rn)
ax[0, 2].set_title("Drizzle: S/N")
# row 1 — after the heal (science identical)
ax[1, 0].imshow(sci, cmap="inferno", origin="lower", norm=sn)
ax[1, 0].set_title("Healed: signal (identical)")
ax[1, 1].imshow(noise_disp, **nl)
ax[1, 1].set_title(tag)
ax[1, 2].imshow(snr_after, cmap="inferno", origin="lower", norm=rn)
ax[1, 2].set_title("Healed: S/N")
# row 2 — differences (before − after)
im = ax[2, 0].imshow(sci - sci, **diff_lim(sci - sci))
ax[2, 0].set_title("Δ signal = 0 (image untouched)")
im = ax[2, 1].imshow(d_noise, **diff_lim(d_noise))
ax[2, 1].set_title("Δ noise (drizzle − healed) = the stripe")
fig.colorbar(im, ax=ax[2, 1], fraction=0.046, pad=0.04)
im = ax[2, 2].imshow(d_snr, **diff_lim(d_snr))
ax[2, 2].set_title("Δ S/N (healed − drizzle)")
fig.colorbar(im, ax=ax[2, 2], fraction=0.046, pad=0.04)
for a in ax.ravel():
    a.set_xticks([]); a.set_yticks([])
fig.suptitle("J1023+4230 F814W  —  drizzle image kept, stripe interpolated in the lens region "
             "(rows: drizzle / healed / difference; cols: signal / noise / S/N)", fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig(f"{OUT}/hybrid_J1023_compare.png", dpi=130)
print("wrote hybrid_J1023_noise.fits / _stripe.fits / _compare.png")
print("NOTE: heal-down under-estimates noise on genuinely lower-coverage pixels; "
      "use --mask-up for the conservative alternative.")
