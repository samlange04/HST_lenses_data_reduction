#!/usr/bin/env python
"""WFPC2/WF3 sibling of redrizzle_bcfill.py: fill the DQ-flagged bad COLUMNS in the
extracted WF3 frames, un-flag them, rebuild the per-frame IVM from the filled SCI, then
RE-DRIZZLE with the pipeline's exact WFPC2 settings so the weight map comes out uniform by
construction -- no post-hoc stripe detection.

STANDALONE comparison script -- NOT part of the pipeline. Lens is set by the BCFILL_LENS
env var (default J0822+2652, the primary WFPC2 split-visit); any single-visit WFPC2 F606W
lens works too. WF3 carries genuine full-height dead columns (flagged by WFPC2 c1m bits
2 = calibration/mask defect + 256; verified in the detector-frame DQ maps -- ~57 columns on
J0822). final_rot=0 rotates those detector-vertical columns by the exposure roll, so they
show as THICK diagonal stripes in 1/sqrt(WHT).

Two things make this different from the ACS version (redrizzle_bcfill.py), and both matter:
  1. The fill bits are WFPC2's (2|256), NOT ACS's (4|128) -- WFPC2 uses the c1m DQ
     convention. The pipeline "good" set is {8,1024}; a bad-column bit is dropped, so those
     output pixels get fewer of the dithered frames -> the weight deficit / stripe.
  2. WFPC2's noise map is a CUSTOM IVM model (var = SCI/gain + floor^2; DrizzlePac ignores
     the WFPC2 ERR array entirely -- see drizzle_wfpc2_wf3.py). So after filling SCI the IVM
     is REBUILT from the filled SCI, or the noise map would be inconsistent with the fill.
The fill is INTERIOR-ONLY (a flagged pixel is filled only if it has good pixels on both
sides in its row), which naturally leaves the vignetted WF3 border -- also bit-2-flagged --
untouched, since a border pixel has no good pixel further out to interpolate from.

Inputs are the pipeline's already-prepared work-dir files (WF3 extracted + updatewcs
distortion + LACosmic DQ 4096 + per-frame IVM), reused so the baseline re-drizzle matches
the science product and the ONLY baseline-vs-filled difference is the column fill:
    data/drizzle_files/slacs_gold/J0822+2652/f606W/{wf3_*_flt.fits, ivm_wf3_*_flt.fits}
Requires those present -- re-run the F606W drizzle first if that dir was cleared.

Usage (no-arg form orchestrates all stages, each drizzle in its own process):
  python redrizzle_bcfill_wfpc2.py                            # default lens (J0822+2652)
  BCFILL_LENS=J0252+0039 python redrizzle_bcfill_wfpc2.py     # any single-visit F606W lens
  python redrizzle_bcfill_wfpc2.py drizzle baseline|filled    # (internal) one drizzle pass
  python redrizzle_bcfill_wfpc2.py compare                     # rebuild the figure only

Outputs (bolton_test_outputs/, tagged by lens so runs don't clobber):
  redrizzle_wfpc2_<lens>_baseline_{sci,noise}.fits  - standard WF3 drizzle (stripe present)
  redrizzle_wfpc2_<lens>_filled_{sci,noise}.fits    - bad-columns filled pre-drizzle (no stripe)
  redrizzle_wfpc2_<lens>_bcfill_compare.png         - 3x3 (rows: standard / filled / difference;
                                                      cols: signal / noise / S/N)
"""
import glob
import os
import shutil
import subprocess
import sys

import numpy as np
from astropy.io import fits

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# Lens is set via the BCFILL_LENS env var (default the split-visit demonstrator), which the
# per-drizzle subprocesses inherit automatically. Outputs + scratch are tagged by lens so
# runs on different lenses don't clobber each other. Any single-visit WFPC2 F606W lens works
# as-is; a split-visit lens (J0728+3835 / J0822+2652) is per-visit -- point SRC at the visit
# whose drizzle_files you want (the primary f606W dir here).
LENS = os.environ.get('BCFILL_LENS', 'J0822+2652')
SRC = f"{REPO}/data/drizzle_files/slacs_gold/{LENS}/f606W"    # pipeline-prepared WF3 inputs
OUT = f"{REPO}/bolton_test_outputs"                           # tracked (final products)
os.makedirs(OUT, exist_ok=True)
WORK = f"{REPO}/output/redrizzle_wfpc2_work/{LENS}"           # git-ignored scratch, per lens

BADCOL = 2 | 256                       # WFPC2 c1m: calibration/mask defect (2) + bit 256
REGION = 180                           # half-size of the comparison crop (px); 180 -> 18"
WF3_NATIVE_SCALE, OUT_SCALE, OUT_PIXFRAC = 0.0996, 0.05, 1.0
WF3_READNOISE = 5.2                    # electrons (WFPC2 IHB Table 4.2), for the IVM floor fallback

# lens position (same source + North-up common-WCS the drizzle scripts use)
sys.path.insert(0, f"{REPO}/info")
from slacs_coords import slacs_coords
from astropy.coordinates import SkyCoord
import astropy.units as u
_lc = SkyCoord(*slacs_coords[LENS], unit=(u.hourangle, u.deg))
RA0, DEC0 = float(_lc.ra.deg), float(_lc.dec.deg)


def fill_bad_columns(wf3_file, ivm_file):
    """In-place: interpolate SCI across DQ 2|256 columns (INTERIOR ONLY), clear those bits
    on the pixels filled, and update the IVM at exactly those pixels (reusing the pipeline's
    stored per-frame gain+floor, IVMGAIN/IVMFLOOR). Updating only the filled pixels -- rather
    than rebuilding the whole IVM, whose from-scratch floor estimate shifts slightly once the
    LACosmic DQ 4096 flags are present -- keeps every other pixel byte-identical to the
    pipeline IVM, so the baseline-vs-filled difference is EXACTLY the dead columns and nothing
    else (matching the ACS redrizzle_bcfill.py 'only the fill differs' design). Returns the
    fill count."""
    filled = []
    with fits.open(wf3_file, mode="update") as h:
        sci, dq = h['SCI', 1].data, h['DQ', 1].data
        bad = (dq & BADCOL) > 0
        for row in np.where(bad.any(axis=1))[0]:
            b = bad[row]
            xg = np.where(~b)[0]
            xb = np.where(b)[0]
            if len(xg) < 2:
                continue
            interior = (xb > xg.min()) & (xb < xg.max())     # skip the vignetted border
            xbi = xb[interior]
            if xbi.size == 0:
                continue
            sci[row, xbi] = np.interp(xbi, xg, sci[row, xg])
            dq[row, xbi] &= ~BADCOL                            # un-flag -> full weight
            filled.append((row, xbi))
        h['SCI', 1].data = sci
        h['DQ', 1].data = dq
        h.flush()

    # Update IVM only at the filled pixels, using the pipeline's own gain/floor.
    with fits.open(ivm_file, mode="update") as h:
        ivm = h['IVM', 1].data
        gain = float(h['IVM', 1].header['IVMGAIN'])
        floor2 = float(h['IVM', 1].header['IVMFLOOR'])
        with fits.open(wf3_file) as hs:
            sci = hs['SCI', 1].data
        for row, xbi in filled:
            var = np.clip(sci[row, xbi].astype(np.float64), 0.0, None) / gain + floor2
            ivm[row, xbi] = np.where(var > 0, 1.0 / var, 0.0).astype(ivm.dtype)
        h['IVM', 1].data = ivm
        h.flush()
    return int(sum(x.size for _, x in filled))


def prep():
    for kind in ("baseline", "filled"):
        d = f"{WORK}/{kind}"
        if os.path.isdir(d):
            shutil.rmtree(d)
        os.makedirs(d)
        wf3 = sorted(glob.glob(f"{SRC}/wf3_u*_flt.fits"))
        if not wf3:
            sys.exit(f"no WF3 inputs in {SRC} -- re-run the F606W drizzle first "
                     "(drizzle_wfpc2_wf3.py --lens J0822+2652 --filt f606W --pa <primary>)")
        assoc = f"{d}/wf3_ivm_association.lst"
        with open(assoc, "w") as af:
            for src in wf3:
                base = os.path.basename(src)
                dst = f"{d}/{base}"
                ivm = f"{d}/ivm_{base}"
                shutil.copy(src, dst)                          # copy: drizzle edits DQ/sky
                shutil.copy(f"{SRC}/ivm_{base}", ivm)          # start from the pipeline IVM
                if kind == "filled":
                    nf = fill_bad_columns(dst, ivm)            # updates SCI/DQ + IVM at fills
                    print(f"  filled {base}: {nf} interior bad-column px "
                          f"(SCI + IVM updated at those pixels only)")
                af.write(f"{base} ivm_{base}\n")
    print("prep done")


def run_drizzle(kind):
    """One WFPC2 CR-pass AstroDrizzle, matching drizzle_wfpc2_wf3.py's LACosmic pass
    (the wf3 files already carry the LACosmic DQ 4096; resetbits=0 keeps it, plain
    weighted mean)."""
    sys.path.insert(0, f"{REPO}/scripts")
    import mmap_fits_write                                     # macOS write-hang fix
    mmap_fits_write.install()
    ref = f"{REPO}/data/reference_files"
    os.environ["CRDS_SERVER_URL"] = "https://hst-crds.stsci.edu"
    os.environ["CRDS_PATH"] = ref
    os.environ["uref"] = os.path.join(ref, "references", "hst", "wfpc2") + os.sep
    from drizzlepac import astrodrizzle

    os.chdir(f"{WORK}/{kind}")
    astrodrizzle.AstroDrizzle(
        "@wf3_ivm_association.lst", output=kind,
        preserve=False, build=False, context=False,
        skysub=True, skymethod="localmin",
        driz_sep_wcs=True, driz_sep_scale=WF3_NATIVE_SCALE,
        driz_sep_bits="8,1024", driz_sep_fillval=-1,
        median=False, blot=False, driz_cr=False,
        resetbits=0,                                           # keep the LACosmic 4096 flags
        final_fillval=None, final_bits="8,1024",
        final_wcs=True, final_scale=OUT_SCALE, final_pixfrac=OUT_PIXFRAC,
        final_wht_type="IVM",
        final_rot=0.0, final_ra=RA0, final_dec=DEC0,
        num_cores=1)


def _crop(kind):
    sci_f = glob.glob(f"{WORK}/{kind}/{kind}_dr?_sci.fits")[0]
    wht_f = glob.glob(f"{WORK}/{kind}/{kind}_dr?_wht.fits")[0]
    sci = fits.getdata(sci_f)
    with fits.open(wht_f) as h:
        wht = h[0].data
        cx = int(round(h[0].header["CRPIX1"])) - 1             # lens pixel = CRPIX (final_ra/dec)
        cy = int(round(h[0].header["CRPIX2"])) - 1
    sl = (slice(max(0, cy - REGION), cy + REGION), slice(max(0, cx - REGION), cx + REGION))
    s = sci[sl]
    w = wht[sl]
    with np.errstate(divide="ignore"):
        noise = np.where(w > 0, 1.0 / np.sqrt(w), np.nan)      # WFPC2 IVM: K=1
    return s, noise


def compare():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from astropy.visualization import AsinhStretch, ImageNormalize, PercentileInterval

    b_sci, b_noise = _crop("baseline")
    f_sci, f_noise = _crop("filled")
    for kind, (s, n) in [("baseline", (b_sci, b_noise)), ("filled", (f_sci, f_noise))]:
        hdr = fits.Header({"BUNIT": "COUNTS/S"})
        fits.writeto(f"{OUT}/redrizzle_wfpc2_{LENS}_{kind}_sci.fits", s.astype(np.float32), hdr, overwrite=True)
        fits.writeto(f"{OUT}/redrizzle_wfpc2_{LENS}_{kind}_noise.fits", n.astype(np.float32), hdr, overwrite=True)

    def snorm(a):
        return ImageNormalize(a, interval=PercentileInterval(99.3), stretch=AsinhStretch())

    def nlim(a):
        lo, hi = np.nanpercentile(b_noise, [2, 98])
        return dict(vmin=lo, vmax=hi, cmap="inferno", origin="lower")

    def rnorm(a):
        return ImageNormalize(a, interval=PercentileInterval(99.0), stretch=AsinhStretch())

    def snr_map(s, n):
        with np.errstate(divide="ignore", invalid="ignore"):
            r = s / n
        r[~np.isfinite(r)] = 0.0
        return r

    def dlim(d):
        v = d[np.isfinite(d)]
        m = float(np.nanpercentile(np.abs(v), 99)) if v.size else 1e-12
        return dict(vmin=-(m or 1e-12), vmax=(m or 1e-12), cmap="RdBu_r", origin="lower")

    b_snr, f_snr = snr_map(b_sci, b_noise), snr_map(f_sci, f_noise)
    sn, rn = snorm(b_sci), rnorm(b_snr)
    d_sci, d_noise, d_snr = b_sci - f_sci, b_noise - f_noise, b_snr - f_snr

    fig, ax = plt.subplots(3, 3, figsize=(13.5, 13.5))
    ax[0, 0].imshow(b_sci, norm=sn, cmap="inferno", origin="lower")
    ax[0, 0].set_title("Standard WF3 drizzle: signal")
    ax[0, 1].imshow(b_noise, **nlim(b_noise))
    ax[0, 1].set_title("Standard: noise — dead-column stripes")
    ax[0, 2].imshow(b_snr, norm=rn, cmap="inferno", origin="lower")
    ax[0, 2].set_title("Standard: S/N")
    ax[1, 0].imshow(f_sci, norm=sn, cmap="inferno", origin="lower")
    ax[1, 0].set_title("Bad-cols filled: signal")
    ax[1, 1].imshow(f_noise, **nlim(f_noise))
    ax[1, 1].set_title("Filled: noise — uniform weight, no stripe")
    ax[1, 2].imshow(f_snr, norm=rn, cmap="inferno", origin="lower")
    ax[1, 2].set_title("Filled: S/N")
    for col, (d, lab) in enumerate([(d_sci, "signal"), (d_noise, "noise"), (d_snr, "S/N")]):
        im = ax[2, col].imshow(d, **dlim(d))
        ax[2, col].set_title(f"Δ {lab} (standard − filled)"
                             + (" = removed stripes" if lab == "noise" else
                                " ≈ 0" if lab == "signal" else ""))
        fig.colorbar(im, ax=ax[2, col], fraction=0.046, pad=0.04)
    for a in ax.ravel():
        a.set_xticks([]); a.set_yticks([])
    fig.suptitle(f"{LENS} F606W (WFPC2/WF3) — input-level bad-column fill + re-drizzle "
                 "(identical AstroDrizzle+IVM settings; only the fill differs; "
                 "rows: standard / filled / difference, cols: signal / noise / S/N)",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(f"{OUT}/redrizzle_wfpc2_{LENS}_bcfill_compare.png", dpi=130)
    m = np.isfinite(b_noise) & np.isfinite(f_noise)
    print(f"median noise  standard {np.nanmedian(b_noise[m]):.5f}  ->  filled "
          f"{np.nanmedian(f_noise[m]):.5f}")
    print(f"science max |diff| = {np.nanmax(np.abs(d_sci)):.3g}  "
          f"(median |diff| {np.nanmedian(np.abs(d_sci)):.3g})")
    print(f"wrote redrizzle_wfpc2_{LENS}_{{baseline,filled}}_{{sci,noise}}.fits + "
          f"redrizzle_wfpc2_{LENS}_bcfill_compare.png")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    if mode == "drizzle":
        run_drizzle(sys.argv[2])
    elif mode == "prep":
        prep()
    elif mode == "compare":
        compare()
    else:
        prep()
        for kind in ("baseline", "filled"):
            print(f"\n=== drizzling {kind} (subprocess) ===")
            subprocess.run([sys.executable, __file__, "drizzle", kind], check=True)
        compare()
