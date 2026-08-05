#!/usr/bin/env python
"""Input-level bad-column fix: fill the DQ-flagged bad columns in the FLC frames, un-flag
them, then RE-DRIZZLE with the pipeline's exact settings so the weight map comes out
uniform by construction -- no post-hoc stripe detection at all.

STANDALONE comparison script -- NOT part of the pipeline. This is the principled sibling of
stripe_heal.py: instead of interpolating the stripe out of the finished noise map, it
removes the cause. The ACS dead columns are DQ bit 128 (bad column) + bit 4 (bad detector
pixel) -- neither is in this repo's ACS "good" set {16,64,256}, so AstroDrizzle drops them
and the affected output pixels get contributions from fewer of the dithered frames -> the
weight deficit that shows as a diagonal stripe in 1/sqrt(WHT).

Per frame, per chip: the flagged pixels (whole detector columns) are filled by linear
interpolation ACROSS the column (from the good pixels either side in the same row), the
ERR array is filled the same way, and DQ bits 4|128 are cleared so drizzle now weights
them like any other pixel. The frames are then re-drizzled with the *identical* no-CR
AstroDrizzle call the pipeline uses (final_bits='256,64,16', ERR weighting, 0.05"/px,
North-up at the lens, num_cores=1, mmap write workaround). A baseline drizzle of the
UNMODIFIED copies is run the same way, so the only difference between the two products is
the bad-column fill -> the comparison is clean and attributable.

Applied to J1023+4230 F814W (bad column through the deflector core).

Usage (the no-arg form orchestrates all stages, each drizzle in its own process):
  python redrizzle_bcfill.py                 # prep -> drizzle baseline -> drizzle filled -> compare
  python redrizzle_bcfill.py drizzle baseline|filled   # (internal) one drizzle pass
  python redrizzle_bcfill.py compare                    # rebuild the figure only

Outputs (output/bolton_investigations/):
  redrizzle_baseline_{sci,noise}.fits  - standard drizzle (stripe present)
  redrizzle_filled_{sci,noise}.fits    - bad-columns filled pre-drizzle (no stripe)
  redrizzle_bcfill_compare.png         - 2x3 comparison
"""
import glob
import os
import shutil
import subprocess
import sys
import warnings

import numpy as np
from astropy.io import fits

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CAL = f"{REPO}/data/calibrated/slacs_gold/J1023+4230/f814W"
OUT = f"{REPO}/output/bolton_investigations"
os.makedirs(OUT, exist_ok=True)
WORK = f"{OUT}/redrizzle_work"
BADBITS = 4 | 128                      # ACS bad detector pixel (4) + bad column (128)
REGION = 200                           # half-size of the comparison crop (px); 200 -> 20"

# lens position (same source the drizzle scripts use)
sys.path.insert(0, f"{REPO}/info")
from slacs_coords import slacs_coords
from astropy.coordinates import SkyCoord
import astropy.units as u
_lc = SkyCoord(*slacs_coords["J1023+4230"], unit=(u.hourangle, u.deg))
RA0, DEC0 = float(_lc.ra.deg), float(_lc.dec.deg)


def fill_bad_columns(dst):
    """In-place: interpolate SCI+ERR across DQ 4|128 columns and clear those DQ bits."""
    n_filled = 0
    with fits.open(dst, mode="update") as h:
        for sci_e, err_e, dq_e in [(1, 2, 3), (4, 5, 6)]:          # both WFC chips
            sci, err, dq = h[sci_e].data, h[err_e].data, h[dq_e].data
            bad = (dq & BADBITS) > 0
            for row in np.where(bad.any(axis=1))[0]:
                b = bad[row]
                xg = np.where(~b)[0]
                xb = np.where(b)[0]
                if len(xg) < 2:
                    continue
                sci[row, xb] = np.interp(xb, xg, sci[row, xg])     # bilinear across column
                err[row, xb] = np.interp(xb, xg, err[row, xg])
            h[dq_e].data = dq & ~BADBITS                            # un-flag -> full weight
            n_filled += int(bad.sum())
        h.flush()
    return n_filled


def prep():
    for kind in ("baseline", "filled"):
        d = f"{WORK}/{kind}"
        if os.path.isdir(d):
            shutil.rmtree(d)
        os.makedirs(d)
        for src in sorted(glob.glob(f"{CAL}/*flc.fits")):
            dst = f"{d}/{os.path.basename(src)}"
            shutil.copy(src, dst)                                   # copy: drizzle edits DQ
            if kind == "filled":
                nf = fill_bad_columns(dst)
                print(f"  filled {os.path.basename(src)}: {nf} bad-column px "
                      f"interpolated + un-flagged")
    print("prep done")


def run_drizzle(kind):
    """One no-CR AstroDrizzle pass, matching drizzle_acs_wfc.py lines 196-213."""
    sys.path.insert(0, f"{REPO}/scripts")
    import mmap_fits_write                                          # macOS write-hang fix
    mmap_fits_write.install()
    ref = f"{REPO}/data/reference_files"
    os.environ["CRDS_SERVER_URL"] = "https://hst-crds.stsci.edu"
    os.environ["CRDS_PATH"] = ref
    os.environ["jref"] = os.path.join(ref, "references", "hst", "acs") + os.sep
    from drizzlepac import astrodrizzle

    os.chdir(f"{WORK}/{kind}")
    flc = sorted(glob.glob("*flc.fits"))
    astrodrizzle.AstroDrizzle(
        flc, output=kind,
        preserve=False, build=False, context=False,
        skysub=True, skymethod="localmin",
        driz_sep_wcs=True, driz_sep_scale=0.05,
        driz_sep_bits="256,64,16", driz_sep_fillval=-1,
        median=False, blot=False, driz_cr=False,
        resetbits=4096,
        final_fillval=None, final_bits="256,64,16",
        final_wcs=True, final_scale=0.05, final_wht_type="ERR",
        final_rot=0.0, final_ra=RA0, final_dec=DEC0,
        num_cores=1)


def _crop(kind):
    sci = fits.getdata(f"{WORK}/{kind}/{kind}_drc_sci.fits")
    with fits.open(f"{WORK}/{kind}/{kind}_drc_wht.fits") as h:
        wht = h[0].data
        cx = int(round(h[0].header["CRPIX1"])) - 1                 # lens pixel = CRPIX
        cy = int(round(h[0].header["CRPIX2"])) - 1
    sl = (slice(cy - REGION, cy + REGION), slice(cx - REGION, cx + REGION))
    s = sci[sl]
    w = wht[sl]
    with np.errstate(divide="ignore"):
        noise = np.where(w > 0, 1.0 / np.sqrt(w), np.nan)          # ACS K=1
    return s, noise


def compare():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from astropy.visualization import (AsinhStretch, ImageNormalize, PercentileInterval)

    b_sci, b_noise = _crop("baseline")
    f_sci, f_noise = _crop("filled")
    for kind, (s, n) in [("baseline", (b_sci, b_noise)), ("filled", (f_sci, f_noise))]:
        hdr = fits.Header({"BUNIT": "ELECTRONS/S"})
        fits.writeto(f"{OUT}/redrizzle_{kind}_sci.fits", s.astype(np.float32), hdr, overwrite=True)
        fits.writeto(f"{OUT}/redrizzle_{kind}_noise.fits", n.astype(np.float32), hdr, overwrite=True)

    def snorm(a):
        return ImageNormalize(a, interval=PercentileInterval(99.3), stretch=AsinhStretch())

    def nlim(a):
        lo, hi = np.nanpercentile(b_noise, [2, 98])
        return dict(vmin=lo, vmax=hi, cmap="inferno", origin="lower")

    dnoise = b_noise - f_noise
    dsci = b_sci - f_sci
    dl = np.nanpercentile(np.abs(dnoise), 99)
    fig, ax = plt.subplots(2, 3, figsize=(16.5, 11))
    ax[0, 0].imshow(b_sci, norm=snorm(b_sci), cmap="inferno", origin="lower")
    ax[0, 0].set_title("Standard drizzle: science")
    ax[0, 1].imshow(b_noise, **nlim(b_noise))
    ax[0, 1].set_title("Standard drizzle: NOISE — stripe through core")
    ax[0, 2].imshow(dnoise, vmin=-dl, vmax=dl, cmap="RdBu_r", origin="lower")
    ax[0, 2].set_title("noise difference (standard − filled) = the removed stripe")
    ax[1, 0].imshow(f_sci, norm=snorm(f_sci), cmap="inferno", origin="lower")
    ax[1, 0].set_title("Bad-cols filled pre-drizzle: science")
    ax[1, 1].imshow(f_noise, **nlim(f_noise))
    ax[1, 1].set_title("Filled + re-drizzled: NOISE — uniform weight, no stripe")
    dsl = np.nanpercentile(np.abs(dsci), 99) or 1e-6
    ax[1, 2].imshow(dsci, vmin=-dsl, vmax=dsl, cmap="RdBu_r", origin="lower")
    ax[1, 2].set_title("science difference (standard − filled) ≈ 0")
    for a in ax.ravel():
        a.set_xticks([]); a.set_yticks([])
    fig.suptitle("J1023+4230 F814W — input-level bad-column fill + re-drizzle "
                 "(identical AstroDrizzle settings; only the fill differs)", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(f"{OUT}/redrizzle_bcfill_compare.png", dpi=130)
    # report + tidy the bulky FLC copies (keep the drc products)
    m = np.isfinite(b_noise) & np.isfinite(f_noise)
    print(f"median noise  standard {np.nanmedian(b_noise[m]):.5f}  ->  filled "
          f"{np.nanmedian(f_noise[m]):.5f} e/s")
    print(f"science max |diff| = {np.nanmax(np.abs(dsci)):.3g} e/s  "
          f"(median |diff| {np.nanmedian(np.abs(dsci)):.3g})")
    for kind in ("baseline", "filled"):
        for f in glob.glob(f"{WORK}/{kind}/*flc.fits"):
            os.remove(f)
    print("wrote redrizzle_{baseline,filled}_{sci,noise}.fits + redrizzle_bcfill_compare.png")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    if mode == "drizzle":
        run_drizzle(sys.argv[2])
    elif mode == "prep":
        prep()
    elif mode == "compare":
        compare()
    else:                                                          # orchestrate all stages
        prep()
        for kind in ("baseline", "filled"):
            print(f"\n=== drizzling {kind} (subprocess) ===")
            subprocess.run([sys.executable, __file__, "drizzle", kind], check=True)
        compare()
