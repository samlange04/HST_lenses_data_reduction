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

Outputs (output/bolton_investigations/):
  hybrid_J1023_noise.fits    - drizzle noise map, stripe healed in the lens region
  hybrid_J1023_stripe.fits   - the detected stripe mask (1 = healed)
  hybrid_J1023_compare.png   - before / mask / after
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
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CUT = f"{REPO}/data/cutouts/slacs_gold/J1023+4230/f814W"
CAL = f"{REPO}/data/calibrated/slacs_gold/J1023+4230/f814W"
OUT = f"{REPO}/output/bolton_investigations"
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

# ── figure ───────────────────────────────────────────────────────────────────────
def nn(a):
    lo, hi = np.nanpercentile(noise, [2, 98])
    return dict(vmin=lo, vmax=hi, cmap="inferno", origin="lower")

fig, ax = plt.subplots(1, 3, figsize=(16.5, 5.6))
ax[0].imshow(noise, **nn(noise)); ax[0].set_title("Drizzle noise (stripe through core)")
ov = np.ma.masked_where(~stripe, stripe)
ax[1].imshow(noise, **nn(noise))
ax[1].imshow(ov, origin="lower", cmap="cool", vmin=0, vmax=1, alpha=0.9)
ax[1].set_title("Detected stripe (interpolate these)")
disp = healed if not args.mask_up else np.where(stripe, np.nan, healed)
ax[2].imshow(disp, **nn(disp)); ax[2].set_title(tag)
for a in ax:
    a.set_xticks([]); a.set_yticks([])
fig.suptitle("J1023+4230 F814W  —  drizzle image kept, stripe interpolated in the lens region",
             fontsize=13)
fig.tight_layout(rect=(0, 0, 1, 0.97))
fig.savefig(f"{OUT}/hybrid_J1023_compare.png", dpi=130)
print("wrote hybrid_J1023_noise.fits / _stripe.fits / _compare.png")
print("NOTE: heal-down under-estimates noise on genuinely lower-coverage pixels; "
      "use --mask-up for the conservative alternative.")
