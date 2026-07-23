# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Environment

All scripts must run inside the `stenv` conda environment (STScI pipeline stack):

```bash
conda run -n stenv python scripts/<script>.py --lens <LENS> --filt <FILTER>
```

To activate interactively: `conda activate stenv`

## macOS write-hang workaround (required)

On this Mac, AstroDrizzle's large FITS output writes hit a kernel lost-wakeup in the
buffered-write path (`tofile → write() → cluster_write → copyin → lck_rw_sleep`) that
wedges the process into an unkillable U-state (~17–26 s CPU in; only a reboot clears it).

`scripts/mmap_fits_write.py` fixes this by monkeypatching `astropy.io.fits.file._array_to_file`
to write via `mmap`+`memcpy` (the `vm_fault`-on-mapped-file path), which dodges the hang.
It is a no-op off macOS and byte-identical to stock astropy. Every drizzle script
imports it and calls `install()` before AstroDrizzle — keep this wiring.
Set `DRIZZLE_MMAP_DEBUG=1` to log each mmap write. `num_cores=1` in all scripts is a related,
separate requirement (parallel `fork` triggers the same U-state).

## What This Repo Does

This is an HST image reduction pipeline for gravitational lens samples (SLACS and BELLS). For each lens+filter combination it:
1. Downloads calibrated exposures from MAST
2. Downloads CRDS reference files
3. Runs `updatewcs` → `TweakReg` (alignment) → `AstroDrizzle` (combination)
4. Produces a no-CR-rejection drizzled mosaic, and a CR-rejected one too for WFPC2 or when `--cr` is given
5. Updates three JSON tracking files in `info/`

## Running a Single Lens

```bash
conda run -n stenv python scripts/drizzle_wfpc2_wf3.py --lens J0008-0004 --filt f606W --sample slacs
conda run -n stenv python scripts/drizzle_acs_wfc.py  --lens J0008-0004 --filt f814W --sample slacs
conda run -n stenv python scripts/drizzle_wfc3_ir.py  --lens J0008-0004 --filt f160W --sample slacs
conda run -n stenv python scripts/drizzle_nic2.py     --lens J0008-0004 --filt f160W --sample slacs  # deprioritised, see below
conda run -n stenv python scripts/drizzle_acs_wfc.py  --lens J0216-0813 --filt f555W --sample slacs
```

All scripts are idempotent: they skip MAST download if calibrated files are already present, and skip the entire drizzle if the final output already exists in `data/drizzled/`. To force a re-run, delete the lens's directory under `data/drizzled/` (and `data/drizzle_files/`).

The ACS, WFC3/IR, and NICMOS scripts accept `--cr` to enable the CR-rejection drizzle pass (disabled by default). Without `--cr`, only the no-CR-rejection pass runs. WFPC2 always runs both passes.

## Running All Lenses (WFPC2 / SLACS)

```bash
bash scripts/run_wfpc2_wf3.sh   # WF3 F606W drizzle + cutout, all 22 lenses
bash scripts/run_all_lenses.sh  # WF3 F606W drizzle only, with a retry pass
```

Both drive `drizzle_wfpc2_wf3.py`. `run_wfpc2_wf3.sh` also runs `make_cutouts.py`
after each drizzle; `run_all_lenses.sh` does not, but retries failures once.
Neither carries an exclusion list — the drizzle script measures each lens's dither
coverage itself and skips any lens that cannot reach 0.05″/px.

`scripts/drizzle_wfpc2_pc.py` is **superseded** and refuses to run (it extracts the
wrong chip *and* deletes the WF3 products on its way). Override with
`ALLOW_SUPERSEDED_WFPC2_PC=1` only if you know why you want it.

## Data Flow and Directory Layout

```
data/
  calibrated/<sample>/<lens>/<filter>/   ← downloaded FLT/FLC/CAL files land here
  drizzle_files/<sample>/<lens>/<filter>/ ← working directory; AstroDrizzle runs here
      run.log                              ← stdout+stderr from the full run
      shift_*.txt                          ← TweakReg shift files
      *_single_sci.fits / *_single_wht.fits
      *.png                                ← diagnostic plots
  drizzled/<sample>/<lens>/<filter>/      ← final products copied here
      <prefix>_cr_d*_sci.fits / _d*_wht.fits
      <prefix>_nocrrej_d*_sci.fits / _d*_wht.fits
  cutouts/<sample>/<lens>/<filter>/       ← recentred stamps from make_cutouts.py
      cutout_sci.fits / cutout_noise.fits / cutout.png
  run_logs/                               ← per-lens stdout logs from the batch runners
  reference_files/                        ← CRDS reference files (auto-downloaded once)
```

## Instrument-Specific Script Details

| Script | Input files | MAST product | Ref env var | Pixel scale | Output suffix |
|---|---|---|---|---|---|
| `drizzle_wfpc2_wf3.py` | `u*flt.fits` | FLT / CALWFPC2 | `uref` | 0.0996″ native → 0.05″ out | `_drw_` |
| `drizzle_acs_wfc.py`  | `*flc.fits`  | FLC / CALACS   | `jref` | 0.05″   | `_drc_` |
| `drizzle_wfc3_ir.py`  | `*flt.fits`  | FLT / CALWF3   | `iref` | 0.1283″ native → 0.06″ out | `_drz_` |
| `drizzle_nic2.py`     | `*cal.fits`  | CAL / CALNIC   | `nref` | 0.0756″ | `_drz_` |

The WFPC2 script extracts only the WF3 chip (SCI/ERR/DQ extension 3) into `wf3_`-prefixed files before drizzling. The other instruments are multi-extension MEF files that DrizzlePac handles natively — no chip extraction needed.

### WFPC2: the lens is on WF3, not the PC

For all 22 SLACS lenses with WFPC2 F606W data the lens galaxy falls on **WF3**
(extension 3) at ~pixel (435, 424) — never on the PC. `DETECTOR = PC` in the
primary header, and therefore the `"WFPC2/PC"` label on MAST and in
`lens_instrument.json`, describes the *aperture*, not the chip the target lands
on; the standard WFPC2 full-field aperture centres the target on WF3. Only the
per-extension `DETECTOR` values identify chips (1=PC, 2=WF2, 3=WF3, 4=WF4).

The superseded `drizzle_wfpc2_pc.py` extracted extension 1 and so produced 22
mosaics of blank sky ~79″ from the lens, which is why `make_cutouts.py` failed
on f606W with "catalogue position falls outside the mosaic".

Two consequences for `drizzle_wfpc2_wf3.py`:

- **Chip renumbering.** DrizzlePac indexes chips positionally as `(SCI, 1..N)`,
  so the extracted WF3 extension must be rewritten to `EXTVER=1` or
  `WFPC2InputImage` raises `KeyError: Extension ('SCI', 1) not found`. This does
  not mislabel the chip: `detnum` comes from the `DETECTOR` keyword (stwcs
  `instruments.py` `set_chip`), which stays 3, so the WF3 gain/readnoise row is
  still the one used. DQ bits `8,1024` are WFPC2-wide and carry over unchanged.
- **Output scale.** WF3 is 0.0996″/px. The `WFPC2-BOX` 4-point pattern (spacing
  0.559″) dithers by half a pixel in both axes — note `POSTARG1/2` are all zero,
  the offsets live in the WCS — which supports a ~2× finer grid. Output is
  0.05″/px, `pixfrac=0.8`, matching the ACS mosaics so the optical bands land on
  a common grid. `dither_phase_counts()` measures the phase coverage from the
  exposure WCSs at runtime and the script **exits without writing anything** if
  either axis has fewer than 2 distinct phases — a native-scale mosaic is not
  wanted. There is no hardcoded exclusion list: which exposures a lens has depends
  on the MAST query, so a fixed table goes stale (J0728+3835 looked unusable with
  2 exposures and is fine with 6).

### WFPC2: COPY visits, C0M-only exposures, exposure-time filtering

Three traps in the WFPC2 F606W archive, all of which silently cost exposures:

- **`-COPY` targets are genuine repeat visits, not duplicates.** J0728+3835 has
  2×1100s on 2007-09-14 plus 4×1100s "-COPY" on 2007-11-05, every frame with a
  distinct `t_min`. The old "prefer non-COPY" rule discarded two thirds of the
  data. `drizzle_wfpc2_wf3.py` keeps both and filters on exposure time instead
  (`MIN_EXPTIME = 10s`), which also covers the `EXPTIME=0` case the non-COPY rule
  originally existed for.
- **`WFPC2/WFC`-labelled observations usually ship no FLT** — only raw `C0M`
  (science) + `C1M` (DQ). Broadening `instrument_name` alone therefore changes
  nothing. The script downloads C0M/C1M and converts them via
  `drizzlepac.wfpc2Data.wfpc2_to_flt`. J1218+0830 is the one lens whose extra
  exposures are genuinely WFC-labelled 1100s science frames; on J0728+3835,
  J0822+2652 and J1142+1001 the WFC frames are 0.5s check shots, dropped by
  `MIN_EXPTIME`.
- **Combining visits needs a wider search radius.** Two visits means two guide-star
  solutions and a roll difference (J0822+2652: PA_V3 101.85 vs 87.92), needing
  `searchrad=3`, not 1.

AstroDrizzle output suffix is determined by input file type, not the output name. Use `_drc_` for FLC (ACS), `_drw_` for WFPC2 FLT, `_drz_` for everything else.

## WCS alignment: `--align`, and why it differs by instrument

**Do not unify the three scripts on this point.** Every drizzle script takes
`--align {mast,tweakreg}`, but the correct default is not the same for all of them.
Verified by stacked stellar FWHM (the metric that reflects the actual product):

| Instrument | Default | What it does | Metric |
|---|---|---|---|
| ACS/WFC | `mast` | no `updatewcs`, no TweakReg | FWHM **0.234″** vs 0.286″ (J0330-0020) |
| WFC3/IR | `mast` | no `updatewcs`, no TweakReg | FWHM **0.364″** vs 0.512″ (J0728+3835) |
| WFPC2/WF3 | **`mast`** (per-lens, see audit) | `updatewcs(use_db=True)`, no TweakReg | core-registration scatter, below |

**Why ACS and WFC3 must not re-solve.** MAST delivers those files already fitted to
GSC 2.4.2 or GAIA eDR3 (`WCSNAME = IDC_*-FIT_REL_GSC242` / `-FIT_REL_GAIAeDR3`), with
relative astrometry across a dither sequence good to ~0.05–0.8 px. TweakReg aligns
every frame onto the *first* frame, so for dithered exposures it measures the dither
itself as an error and removes it. On J0330-0020 (`ACS-WFC-DITHER-BOX`, POSTARG ±0.19″
≈ ±3.7 px) its reported `XSH/YSH` of (−5.04, −1.65) px are exactly the POSTARG offsets;
afterwards all four WCSs mapped a sky position to the same detector pixel to 0.01 px.
AstroDrizzle then stacked four dithered frames as if they shared a pointing. Per-frame
WCS error went from spread 0.82 px (a harmless common offset) to 3.61 px — **this is
what smeared point sources and split lensed arcs into offset copies.**
`updatewcs(use_db=False)` compounds it by stripping the `-FIT_REL_*` refinement.

**WFPC2 still needs `updatewcs`, but NOT TweakReg.** It is the only instrument that
cannot skip `updatewcs`: AstroDrizzle needs the NPOL (`NPOLFILE`) and detector-to-image
(`D2IMFILE`) distortion arrays, and the WF3 chip extraction does not carry them over —
without it AstroDrizzle stops and prompts for the missing DGEO correction. `--align
mast` runs `updatewcs(use_db=True)` (which restores the `GSC240` fit) and then stops —
no TweakReg. The `GSC240` solution is only ~0.5″ off in *absolute* astrometry, but its
*relative* frame-to-frame registration is good to ~0.02–0.03″, which is what the stack
needs; the absolute offset is fixed afterwards by `align_wfpc2_to_acs.py`.

**Reversal (2026-07-23): WFPC2 must NOT re-solve with TweakReg either** — it is the same
dither-erasing bug as ACS. Earlier guidance (and an FWHM comparison) had WFPC2 default
to TweakReg, but a **per-lens core-registration audit** overturned it: for the mostly
single-visit F606W lenses TweakReg scatters the frames by ~0.7″, splitting the deflector
core into ~4 offset knots (visible on J0252+0039), while MAST registers the same frames
to ~0.02″. The audit drizzles each frame separately onto the common grid and measures
the deflector-core position scatter (LACosmic-masked so cosmic rays don't inflate it);
the earlier FWHM metric was misleading because it centroided a single extended galaxy.
The per-lens choice is stored in **`info/wfpc2_alignment.json`** and read by the batch
runner. Result: **all 22 lenses use `mast`**. See `[[wfpc2_tweakreg_misregisters]]`.

**Multi-visit lenses are split, not TweakReg'd.** J0728+3835 and J0822+2652 each have
two visits at a ~14–16° roll difference (two guide-star solutions). Rather than let
TweakReg cross-register them (its only historical justification), they are drizzled as
**separate per-visit datasets** — `--pa <PA_V3> --out-suffix _v1/_v2`, each visit
single-guide-star so each uses `mast`. Each product lands in `f606W_v1`/`f606W_v2`.
Outcome: **J0822+2652** = `f606W_v1` (2×1100s) + `f606W_v2` (4×1100s); **J0728+3835** =
`f606W_v2` only (its 2-frame visit has just 1 x-dither-phase and cannot reach 0.05″, so
it is dropped). **J1142+1001 stays combined** — its two visits share a roll (PA 119.00 vs
118.87), so there is no offset to separate. `align_wfpc2_to_acs.py --f606-dir f606W_v1`
and `make_cutouts.py --filt f606W_v1` handle the suffixed dirs.

**Diagnosing this class of bug.** Compare the *spread* of the per-frame WCS error
(predicted vs actual position of a source), not its magnitude — a common offset is a
harmless absolute-astrometry shift, frame-to-frame scatter is what smears a stack. The
cleanest metric is per-frame core-registration scatter on the drizzled common grid (with
CRs masked); confirm against the visible product (one core, not knots). Stacked stellar
FWHM can mislead here — it centroids an extended galaxy and rewards TweakReg's internal
self-consistency even when the deflector is split.

### F606W absolute astrometry: `align_wfpc2_to_acs.py` (run after drizzle, before cutouts)

WFPC2 F606W carries only a GSC 2.4.0 (`GSC240`) solution, ~0.3–1″ off absolute; ACS
F814W and WFC3/IR F160W carry GAIA eDR3 / GSC242 (<0.02″, agree with each other to
~0.01″). So the F606W mosaic sits ~0.5–0.9″ from the other bands (0.66″ measured on
J0252+0039) — a whole-mosaic shift, not a smear, so it does not hurt the F606W stack
but breaks cross-band registration.

`scripts/align_wfpc2_to_acs.py` fixes it by tying the **deflector light-centroid** to
the GAIA-accurate ACS F814W: an iterative windowed `centroid_com` (`stable_centroid`,
robust to the roughly symmetric ring) in both bands, then shifts the F606W `CRVAL1/2`
so the centroids coincide and stamps `GSC240FX=True`. It is idempotent (re-measures the
residual, applies ~0), refuses any tie implying a shift > `MAX_SHIFT = 1.5″`, and where
F160W exists uses it as an **independent** check (not in the fit). Verified on
J0252+0039: 0.66″ → 0.009″ vs F814W and 0.009″ vs F160W.

    conda run -n stenv python scripts/align_wfpc2_to_acs.py --lens J0252+0039   # or --all

Run order: ACS + WFPC2 drizzles → `align_wfpc2_to_acs.py` → `make_cutouts.py`.

## Weight maps: `final_wht_type=ERR` (default), not `EXP`

All four drizzle scripts take `--wht-type {ERR,IVM,EXP}`, default **`ERR`**, and write
the chosen type into every AstroDrizzle pass and forward it to the no-CR subprocess.
`cutout_noise.fits` is `1/sqrt(WHT)`, so the weight type *is* the noise model:

- **`ERR`** — full inverse-variance: source Poisson + sky + read + dark. The correct
  per-pixel σ for modelling. On J1143 the core/sky noise ratio is ~4.0 (source shot
  noise present, as it must be).
- **`EXP`** — effective-exposure-time map, **uncalibrated** (no source shot noise): the
  core/sky ratio comes out a flat ~1.04, i.e. it claims the bright core is as noisy as
  blank sky. Do not use for a likelihood.
- **`IVM`** — background noise only, no object Poisson.

(DrizzlePac Handbook pp.103,139.) The blank-sky floor is **included** in ERR, not
something to add back. See "Drizzle correlated noise" for the covariance caveat that
ERR does *not* capture.

## Common output WCS across filters (orientation + centre)

Every band is pinned to the **same output grid geometry** so the filters co-register
pixel-for-pixel without a later reprojection: the drizzle scripts pass
`final_rot=0.0` (North-up) and `final_ra`/`final_dec` at the lens catalogue position
(the tangent point) into every AstroDrizzle call. The coordinates come from
`info/slacs_coords.py`; a lens missing from that table falls back to the native drizzle
WCS (`_common_wcs = {}`) with a printed warning. This is what lets `make_cutouts.py`
cut all bands on a shared centre — it aligns orientation and tangent point, not the
absolute astrometry, which is still the job of the delivered WCS (and, for F606W,
`align_wfpc2_to_acs.py`).

## TweakReg `threshold` is per-instrument

`threshold` is in **image data units**, so it does not transfer between detectors.
A single value of 200 was tuned on the WFPC2 PC chip and then copied into every
script, where it silently starved `minobj` and broke alignment:

| Script | threshold | searchrad | Why |
|---|---|---|---|
| `drizzle_wfpc2_wf3.py` | 100 | 3 | 200 left too few WF3 sources: 1 of 6 frames nan, another matched spuriously at (40, −65) px. 100 and 50 both give 6/6 and a coherent ~9.5 px visit offset; ≤20 lets false matches back in |
| `drizzle_wfc3_ir.py` | 20 | 1 | FLTs are ELECTRONS/S — sky ~0.7, 99.9th pct ~8. A 200 e/s cut found 4–6 objects/image, below `minobj=7`, so no shiftfile was written and the script died on `FileNotFoundError: shift_flt.txt`. 20 gives ~40–50 objects |
| `drizzle_acs_wfc.py` | (unchanged) | | works on F814W and F555W as-is |

When TweakReg fails with nan shifts, no shiftfile, or one wild outlier, check the
`FINAL number of objects` lines in the log **before** touching `searchrad` — too
few sources is the usual cause. Note the fit is not fragile: on WFC3/IR every
threshold from 50 down to 1.5 gave the identical 4.61 px solution.

These values only take effect under `--align tweakreg`. Since the per-lens audit put all
22 WFPC2 lenses on `mast` (and ACS/WFC3 default to `mast` too), TweakReg is off in the
standard pipeline for every band — this tuning is now comparison-run-only dead code,
kept in case a future lens's audit picks `tweakreg`.

## Cosmic-ray rejection: LACosmic, not `driz_cr` (ACS **and** WFPC2)

`drizzle_acs_wfc.py --cr` masks cosmic rays per frame with LACosmic
(`--cr-method lacosmic`, the default) and then drizzles a plain weighted mean
(`median=False, blot=False, driz_cr=False`). `--cr-method drizcr` restores the old
AstroDrizzle route. **`drizzle_wfpc2_wf3.py` now uses the same LACosmic route for its
CR pass** (`--cr-method lacosmic` default; `--lacosmic-sigclip 4.5`, `--lacosmic-objlim
5.0`; gain/readnoise/saturate taken from the WF3 header — `ATODGAIN`, `WF3_READNOISE=5.2`,
`SATURATE`). The mask is written to DQ bit 4096 with `resetbits=0`, exactly as for ACS.
This replaced `driz_cr` for WFPC2 because `driz_cr` ate the deflector core the same way
(WFPC2 F606W core retention: LACosmic peak 1.000, 1″ 0.979, 2″ 0.933 — matches ACS).

`driz_cr` detects CRs by comparing each frame to a blotted median of the stack. On a
steep PSF core that reference is systematically low, so the core's residual reads as a
cosmic ray: on J0330-0020 it flagged 113–206 real pixels inside the 1″ core of *every*
frame and destroyed 37% of the deflector flux. Loosening `driz_cr_snr`/`driz_cr_scale`
made it **worse** — the fault is the biased reference, not the threshold.

| method | core flux vs no-CR | detections in the 10″ stamp |
|---|---|---|
| no CR rejection | 1.000 | 112 |
| `driz_cr`, `combine_type=minmed` | 0.628 | 0 |
| `driz_cr`, `combine_type=median` | 0.846 | — |
| **LACosmic, `objlim=5`** | **0.993** | **7** |

**`resetbits=0` is mandatory on the LACosmic CR pass.** It defaults to 4096 and would
clear the DRIZ_CR bit the mask is written into, silently producing an un-masked
drizzle that still looks plausible — the first run of this reported a flawless
`core=1.000` that was really the no-CR image scored against itself. Conversely the
no-CR passes pin `resetbits=4096` explicitly, because they run *after* the CR pass on
the same FLCs and would otherwise inherit the mask and become silently CR-rejected.

## NICMOS is deprioritised

Do not generate or propose NICMOS (NIC2) products unless explicitly asked. The
field of view is far too small to be useful — 258x256 px at 0.0756″ is ~19″ across,
against 2318x2052 at 0.06″ (~139″) for WFC3/IR — and the pipeline may be unsound.
F160W coverage questions should be answered from WFC3/IR.

All NICMOS data has been **deleted** (2026-07-21): drizzled products, working dirs,
108 `*cal.fits` exposures, run logs, and the NICMOS CRDS reference cache — 472 MB.
The 24 affected lenses now carry `f160W: null` in all three tracking JSONs. This is
reversible: `scripts/drizzle_nic2.py` is intentionally kept, and re-running it
re-downloads from MAST and re-fetches the CRDS refs automatically. That also means
running it by accident will quietly bring all of it back.

## Output pixel scales

Every band is drizzled to a sub-native grid where the dither supports it. The scale
is chosen by measurement (weight-map uniformity + stellar FWHM), not convention:

| Band | Instrument | Native | Output | pixfrac |
|---|---|---|---|---|
| F606W | WFPC2/WF3 | 0.0996″ | 0.05″ | 1.0 |
| F814W, F555W | ACS/WFC | 0.05″ | 0.05″ (native) | default |
| F160W | WFC3/IR | 0.1283″ | 0.06″ | 1.0 |

F160W sits on 0.06″ while the optical bands are on 0.05″, and this is **kept** — do
not re-drizzle F160W to 0.05″ to grid-match. PyAutoLens ingests each band at its
native scale and pixel-matches at the modelling stage (it even fits sub-pixel
inter-band `grid_offset`/`grid_rotation_angle`), so a common drizzle grid buys nothing.
Measured: 0.05″ opens no *empty* pixels on any of the 13 lenses but worsens weight
non-uniformity (8–19% of pixels below half-median weight; J1430+4105, J1029+0420,
J0841+3824 worst at ~17–19%) with no resolution gain, and pixfrac does not fix it
(the under-coverage is set by the 4-point dither geometry, not the drop size).

**F160W uses pixfrac 1.0, not 0.8.** With the noise now carried on calibrated ERR
weight maps, the drizzle goal is a *uniform, low-correlation* noise map for the
likelihood, not the sharpest PSF. On J0252+0039 (0.06″, correctly registered),
pixfrac 1.0 drops the adjacent-pixel noise correlation from 7.3% to 2.4% (texture
5.4%→3.3%); the PSF is marginally softer but PyAutoLens fits the PSF explicitly. The
older pixfrac-0.8 tuning (chosen on stacked FWHM) predates ERR weighting and the shift
to prioritising noise-map covariance for the modelling.

**F606W also uses pixfrac 1.0, for the same reason.** It is 2× oversampled (0.0996″ →
0.05″) just like F160W, so pixfrac 1.0 likewise lowers the noise-map correlation for the
likelihood at a small PSF cost that PyAutoLens absorbs. ACS F814W/F555W keep the default
pixfrac — they are drizzled at native 0.05″ (no oversampling), so the correlation
penalty is already negligible and there is nothing to buy by lowering it.

**F555W needs no new script** — `drizzle_acs_wfc.py` is filter-agnostic
(`--filt f555W`). All 16 F555W lenses are already reduced; they are exactly the 16
lenses that have no WFPC2 F606W data (proposals 10494/10798).

## Drizzle correlated noise (matters for strong-lens modelling)

Drizzling onto a **finer-than-native grid** with `pixfrac < 1` correlates adjacent
output pixels. Each output pixel is a weighted sum of the input pixels its shrunken
"drop" overlaps; neighbouring output pixels draw on overlapping sets of input pixels,
so their noise is **covariant**. This is a property of drizzle resampling (Fruchter &
Hook 2002; Casertano et al. 2000), not a defect, and it is stronger the more the grid
is oversampled relative to native.

This shows up directly as visible **texture in the noise maps**, measured as the
fractional pixel-to-pixel RMS of a blank-sky region of the weight-derived noise map:

| Band | Native → Output | Oversample | pixfrac | Noise-map texture |
|---|---|---|---|---|
| F814W | 0.05″ → 0.05″ | 1.0× | default | ~6% |
| F160W | 0.1283″ → 0.06″ | 2.1× | 1.0 | ~3% (was ~5–9% at pixfrac 0.8) |
| F606W | 0.0996″ → 0.05″ | 2.0× | 1.0 | ~4% (see correction below) |

The texture tracks the oversampling, but **less than a first look suggested, and
pixfrac 1.0 suppresses most of it.** Two corrections to the earlier reading of this:

- **The F606W band was never the ~19% outlier the first table recorded.** That figure
  was dominated by **TweakReg misregistration**, not drizzle oversampling: TweakReg
  scattered the (mostly single-visit) F606W frames by ~0.7″, so the "texture" was
  really coverage patchwork from mis-stacked frames (this is the "large patches that
  ruin the S/N" seen on J0252). Once each lens is drizzled on its correctly-registered
  WCS (MAST for single-visit, TweakReg only where a multi-visit roll needs it — see
  the alignment audit), the blank-sky texture drops to ~4% at pixfrac 0.8 and lower at
  pixfrac 1.0 — comparable to F160W, **not** a 3× outlier. See `[[wfpc2_tweakreg_misregisters]]`.
- **F606W and F160W now ship at pixfrac 1.0**, chosen specifically to minimise this
  correlation for the likelihood (F160W adjacent-pixel correlation 7.3%→2.4% going
  0.8→1.0 on J0252; the PSF is fit explicitly so the marginal softening is free). The
  older pixfrac-0.8 numbers above are the pre-change state, kept for comparison.

`final_wht_type=ERR` makes the **per-pixel variance** correct (it carries source
Poisson + sky + read + dark), but it says nothing about the **off-diagonal covariance**
between neighbouring pixels — a per-pixel σ map understates the true correlated noise.
On J1143 blank sky (F606W, pixfrac 0.8) the empirical pixel RMS was ~**1.47×** what the
ERR map predicted; pixfrac 1.0 shrinks this factor but does not remove it — any
oversampled band still has some residual correlation the diagonal σ map misses.

**Modelling implication.** A per-pixel independent-Gaussian likelihood (diagonal
covariance) mis-estimates parameter *uncertainties* — it does not bias the best-fit
point, but it makes error bars/evidence wrong, and worse in the oversampled
F606W/F160W bands than in native-scale F814W. Options: build the pixel covariance into
the likelihood; inflate the per-pixel noise by the measured correlation factor as a
crude correction; or follow Bayer et al. (2023), who characterise the correlated noise
with azimuthally-averaged blank-sky power spectra and fold that into the noise model.
Prefer native-scale F814W where a clean per-pixel noise model matters most, and be
explicit about the correlation whenever F606W/F160W drive the constraint.

## Cutouts (`scripts/make_cutouts.py`)

```bash
conda run -n stenv python scripts/make_cutouts.py --lens J0029-0055 --filt f606W --sample slacs
```

Cuts a square stamp (default 10″) from `data/drizzled/` into `data/cutouts/`,
writing `cutout_sci.fits`, `cutout_noise.fits` (from the weight map) and a
3-panel PNG.

The stamp is recentred on the galaxy rather than the catalogue position, and the
peak search runs on the **CR-rejected** mosaic when one exists, even though the
science stamp is cut from the no-CR pass. The no-CR mosaics contain cosmic rays
by construction and a brightest-pixel search locks onto them. Suppressing CRs
purely by widening the median window needs ~21 pix ≈ 1″ at 0.05″/px, wide enough
to bias the peak of a compact galaxy, so centring uses data with no CRs instead.
The CR pass loses core flux and is unfit for science, but its core is still the
local maximum, which is all centring needs. Do not diagnose a bad recentre by
reaching for `--median-size` first — check that a `*_cr_*` mosaic is present.

Offsets around 1–2″ are not necessarily failures: several lenses (J0912+0029,
J0956+5100) have genuine multi-knot morphology, so the brightest pixel is a knot
rather than the catalogue centroid, and the stamp is still correctly placed.

**Shared centre across bands (`--center-band`, default `f814W`).** All bands are cut
about a single centre measured from the `--center-band` mosaic (F814W: highest S/N,
GAIA-accurate, cleanest), so the stamps co-register across filters instead of each
band recentring on its own peak (which would leave them offset by their individual
morphology). `--center-self` restores per-band recentring; a band whose center-band
products are missing falls back to its own peak with a warning. The peak search itself
still prefers the CR mosaic (`--pass {auto,cr,nocrrej}`, `auto` = CR if present).
Combined with the common output WCS and `align_wfpc2_to_acs.py`, the three optical/IR
bands land on a common pixel grid.

## Tracking JSONs in `info/`

Three files are updated automatically by every script run:

- **`lens_products.json`** — `{lens: {filter: [obsid, ...]}}` — HST rootnames downloaded from MAST
- **`lens_instrument.json`** — `{lens: {filter: "INSTRUME/DETECTOR"}}` — e.g. `"ACS/WFC"`, `"WFPC2/WF3"`, `"WFC3/IR"`
- **`lens_exptime.json`** — `{lens: {filter: exptime_seconds}}` — from the CR-rejected drizzle header

If a lens has no data for a filter, the value is stored as `null`. All three files
carry a key for every band (`f160W, f555W, f606W, f814W`) on every lens, so a
missing key means the file is out of sync, not that data is absent.

**`null` is overloaded for F160W.** The 24 F160W nulls do *not* mean "no data on
MAST" — those lenses have NICMOS F160W observations that were deliberately deleted
(see *NICMOS is deprioritised*). Only the 13 WFC3/IR lenses have F160W data that is
wanted. `drizzle_nic2.py` also writes `null` when a MAST query genuinely returns
nothing, so the JSONs alone cannot distinguish "dropped on purpose" from "never
existed" — this note is the record. Re-running `drizzle_nic2.py` would silently
repopulate those 24 entries and re-download the data.

Both levels are kept sorted: lenses across the file, and **filters within each lens
entry**. `_update_info_json` and the `lens_products` write in every drizzle script
re-sort on each update, so the ordering survives partial runs.

Note `lens_instrument.json` records `WFPC2/WF3` for F606W, not the `WFPC2/PC` string
MAST uses — it names the chip the data actually came from.

## Lens Samples

- **SLACS** (38 lenses): HST proposals 10886, 11202, 10494, 10798. Coverage as of the last MAST survey:

  | Band | Instrument | Lenses | Notes |
  |---|---|---|---|
  | F814W | ACS/WFC | 38 | all |
  | F606W | WFPC2/WF3 | 22 | the other 16 have no WFPC2 data in any filter |
  | F555W | ACS/WFC | 16 | exactly the 16 without F606W (props 10494/10798) |
  | F160W | WFC3/IR | 13 | all proposal 11202, ~2397 s each |
  | F160W | NICMOS/NIC2 | 24 | deprioritised — data deleted, entries are `null` |

  So every lens has F814W, and the sample splits cleanly into 22 with WFPC2 F606W
  and 16 with ACS F555W. F160W is WFC3/IR for 13; the remaining 24 have only
  NICMOS, which is not wanted.
- **BELLS** (16 lenses): WFC3/UVIS multi-band (F225W, F275W, F438W, F606W, F814W). Reduction scripts not yet written for this sample.
- **GALLERY**: HST proposals 14189, 16734.

Target names on MAST follow the pattern `SDSS<LENS>` (e.g. `SDSSJJ0008-0004`). The MAST query uses `target_name=f'SDSS{lens}%'` with a wildcard to handle minor naming variations.

**COPY handling differs by instrument, deliberately:**

- **ACS** (`drizzle_acs_wfc.py`) filters COPY observations out in favour of non-COPY when both exist — **except** for lenses in `mast_target_names.FORCE_COPY_LENSES`, whose non-COPY observations are unusable. Currently `{J1032+5322}` (F814W): the non-COPY frames are `EXPTIME=0`, so `force_copy(lens)` selects the COPY frames. J1032+5322 is the only ACS lens with COPY data at all — surveyed across F814W and F555W, nothing else has any.
- **WFPC2** (`drizzle_wfpc2_wf3.py`) keeps **both**, because there the COPY sets are genuine repeat visits carrying most of the usable exposure time (see above). It rejects junk on `MIN_EXPTIME` instead of on target name.

Do not "unify" these two policies without re-checking the archive — they encode different facts about different datasets.

### Non-standard MAST target names (`GAL-*`)

Some lenses are **not** on MAST under an `SDSS<LENS>` name — they use a `GAL-<plate>-<mjd>-<fiber>` designation instead, so the default `SDSS{lens}%` query returns no observations. The drizzle scripts resolve this automatically via `scripts/mast_target_names.py`: for a lens in the table below they query the `GAL-*` name first and fall back to `SDSS{lens}%`. **Output/directory names always stay in the J convention** (the left column) regardless of the MAST name used to fetch. Not all of these are on the current lens list yet; recorded here for future additions.

| Output name (J convention) | MAST target name |
|---|---|
| J0216-0813 | GAL-0668-52162-428 |
| J0737+3216 | GAL-0541-51959-145 |
| J0912+0029 | GAL-0472-51955-429 |
| J0956+5100 | GAL-0902-52409-068 |
| J0959+0410 | GAL-0572-52289-495 |
| J1205+4910 | GAL-0969-52442-134 |
| J1250+0523 | GAL-0847-52426-549 |
| J1402+6321 | GAL-0605-52353-503 |
| J1420+6019 | GAL-0788-52338-605 |
| J1627-0053 | GAL-0364-52000-084 |
| J1630+4520 | GAL-0626-52057-518 |
| J2238-0754 | GAL-0722-52224-442 |
| J2300+0022 | GAL-0677-52606-520 |
| J2303+1422 | GAL-0743-52262-304 |

## AstroDrizzle Key Parameters

Two drizzle passes are defined, over the same input files:
- **CR pass** (`median=True, blot=True, driz_cr=True`): cosmic-ray cleaned science image
- **No-CR pass** (`median=False, blot=False, driz_cr=False`): uncleaned, useful for comparison

Only WFPC2 runs both unconditionally. ACS, WFC3/IR and NICMOS run the no-CR pass
only unless `--cr` is passed. Current on-disk state reflects that: F606W, F814W and
F555W have both passes, but **the 13 WFC3/IR F160W mosaics have no CR pass** (the 24
CR products under `f160W` are all NICMOS).

That has a consequence for `make_cutouts.py`, which prefers the CR mosaic for
recentring: on those 13 it falls back to the science pass. Acceptable here because
WFC3/IR FLTs are already up-the-ramp CR-rejected, so the no-CR mosaic is not
CR-infested the way a WFPC2 or ACS one is — but re-run with `--cr` if a recentre
looks wrong.

DQ bits treated as good pixels:
- WFPC2: `8,1024` (full-well saturated, cosmic ray corrected)
- ACS/WFC: `256,64,16` (full-well saturated, warm pixels, stable hot pixels)
- WFC3/IR: `64,512` (warm pixels, blobs)
- NICMOS: `2,4,8` (uncertain linearity/dark/flat corrections — acceptable calibration imperfections; all defects, saturation, CRs excluded. Per NICMOS Data Handbook Table 2.3)
