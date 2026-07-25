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
3. Aligns and combines with `AstroDrizzle`. By default all bands trust the delivered
   MAST WCS (`--align mast`): ACS/WFC3 skip `updatewcs` and TweakReg entirely; WFPC2 runs
   `updatewcs` (for the distortion arrays) but **not** TweakReg. TweakReg is opt-in
   (`--align tweakreg`) and, per the alignment audit, not used by any lens — see the WCS
   alignment section for why re-solving smears the stack.
4. Produces a no-CR-rejection drizzled mosaic, and a CR-rejected one (LACosmic) too for
   WFPC2 or when `--cr` is given
5. Updates three JSON tracking files in `info/`

## Running a Single Lens

`--sample` defaults to **`slacs_gold`** everywhere and is shown below only for the
record. It sets the `<sample>` level of every `data/` path
(`data/drizzled/slacs_gold/<lens>/<filt>/`), so passing the wrong one silently writes a
correct product into the wrong tree.

```bash
conda run -n stenv python scripts/drizzle_wfpc2_wf3.py --lens J0008-0004 --filt f606W --sample slacs_gold
conda run -n stenv python scripts/drizzle_acs_wfc.py  --lens J0008-0004 --filt f814W --sample slacs_gold
conda run -n stenv python scripts/drizzle_wfc3_ir.py  --lens J0008-0004 --filt f160W --sample slacs_gold
conda run -n stenv python scripts/drizzle_acs_wfc.py  --lens J0216-0813 --filt f555W --sample slacs_gold
```

All scripts are idempotent: they skip MAST download if calibrated files are already present, and skip the entire drizzle if the final output already exists in `data/drizzled/`. To force a re-run, delete the lens's directory under `data/drizzled/` (and `data/drizzle_files/`).

A lens with no data for the requested instrument+filter prints `=== NO DATA: ...`,
records `null` in the tracking JSONs and **exits 0** — see *Lens Samples*.

The ACS, WFC3/IR, and NICMOS scripts accept `--cr` to enable the CR-rejection drizzle pass (disabled by default). Without `--cr`, only the no-CR-rejection pass runs. WFPC2 always runs both passes.

## Running All Lenses

Every runner takes an optional sample argument and defaults to `slacs_gold`:

```bash
bash scripts/run_acs_all.sh                  # ACS/WFC F814W + F555W
bash scripts/run_wfc3_all.sh                 # WFC3/IR F160W
bash scripts/run_wfpc2_wf3.sh                # WFPC2/WF3 F606W: drizzle -> align -> cutout
bash scripts/run_cutouts_all.sh              # stamps for whatever products exist
bash scripts/run_acs_all.sh slacs_other      # any runner, any sample
```

All four take their lens roster from `info/lens_samples.json` via
`scripts/mast_target_names.py` (`run_cutouts_all.sh` is the deliberate exception — it
globs `data/drizzled/`, because a stamp can only be cut from a mosaic that exists). They
report `ok` / `no data` / `FAILED` separately, so the 16 `slacs_gold` lenses with no
WFPC2 data are not mistaken for errors.

**`run_wfpc2_wf3.sh` is the single WFPC2 driver, and it runs the full three-stage order**
(`drizzle_wfpc2_wf3.py` → `align_wfpc2_to_acs.py` → `make_cutouts.py`) for every lens in
the sample. It reads the per-lens alignment mode from `info/wfpc2_alignment.json`, expands
the two split-visit lenses into their per-visit products, and retries failures once. It
carries no exclusion list — the drizzle script measures each lens's dither coverage and
skips any lens that cannot reach 0.05″/px, which the runner reports as
`SKIPPED (dither phase)` rather than counting as a failure. It stops before the align and
cutout stages on a `no data` lens, since both would fail on a product that does not exist.

**It was rebuilt on 2026-07-26 because it had drifted into producing quietly wrong
products** — worth knowing, because the same three traps apply to anything else that
drives this pipeline:

- it passed **no `--align`**, so it took the drizzle script's default, which was
  still `tweakreg` — the mode the per-lens audit rejected for all 22 lenses. **The
  script default is now `mast`**; CLAUDE.md had documented the reversal but the
  `argparse` default was never changed to match, so the audit result was live only
  for anyone who passed the flag by hand.
- it had **no split-visit handling**, so J0728+3835 and J0822+2652 would have been
  drizzled as single combined datasets across a ~15° roll difference *and* would have
  rewritten their per-visit tracking-JSON keys back into a bogus combined `f606W`.
- it **skipped `align_wfpc2_to_acs.py`**, which is not optional: a re-drizzle discards
  the tie (the tie is a `CRVAL1/2` edit on the drizzled product), so it has to be
  re-applied after *every* drizzle. Skipping it gives stamps that look perfect alone
  and are ~0.3–0.9″ off the other bands.

`scripts/stale_scripts/run_all_lenses.sh` is **retired** — it was a second, independently
maintained driver with all three faults above, and keeping two is what allowed the drift.
It now refuses to run and points at `run_wfpc2_wf3.sh`, which absorbed its retry pass.

`scripts/run_cutouts_all.sh` globs `<filt>*`, not `<filt>`, so the per-visit product
directories are included. With a bare `<filt>` the two split-visit lenses matched
nothing and their stamps were silently never regenerated — no error, just two lenses
left on stale cutouts.

### `scripts/stale_scripts/` (moved there 2026-07-26)

Four scripts live in `scripts/stale_scripts/` rather than alongside the live pipeline.
They are retained on purpose, not deleted — but **nothing in the pipeline invokes
them**, and the path change is itself part of the guard:

- `drizzle_wfpc2_pc.py` — **superseded**: extracts the wrong chip *and* `rmtree`s the
  good WF3 products on its way. Override `ALLOW_SUPERSEDED_WFPC2_PC=1`.
- `drizzle_nic2.py` — **deprioritised**: an accidental run re-downloads ~472 MB from
  MAST and repopulates the `f160W` NICMOS entries. Override `ALLOW_NICMOS=1`.
- `run_all_lenses.sh` — retired WFPC2 driver (see above); refuses to run.
- `rebuild_stenv_arm64.sh` — one-off plan to rebuild `stenv` as native arm64, written
  when the U-state write hang was thought to be a Rosetta artefact. Never used: the
  cause was the buffered-write path, fixed in-process by `scripts/mmap_fits_write.py`.

The first two also **raise `NotImplementedError` on import**. That raise is deliberate —
it fails with a traceback and a non-zero status a batch runner cannot mistake for a
clean skip, and it fires before any MAST or CRDS network call.

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
  0.05″/px (matching the ACS mosaics' scale) at `pixfrac=1.0` (raised from 0.8 to lower
  noise-map correlation; see the drizzle-correlated-noise and pixel-scale sections).
  `dither_phase_counts()` measures the phase coverage from the
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
- **Multi-visit lenses are now split, not TweakReg-combined.** Two visits means two
  guide-star solutions and a roll difference (J0822+2652: PA_V3 101.85 vs 87.92). The old
  approach cross-registered them with TweakReg (`searchrad=3`); the current pipeline
  instead drizzles each visit as a separate MAST dataset (`--pa`/`--out-suffix`) — see the
  split-visit note in the WCS alignment section. (The wider `searchrad=3` still applies if
  you ever force `--align tweakreg` on a combined multi-visit run.)

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

That claim was aspirational until 2026-07-26: the runner passed no `--align` at all
and the script's `argparse` default was still `tweakreg`, so **the audit result was
live only for someone who passed the flag by hand.** Both are fixed — the default is
`mast`, and `run_wfpc2_wf3.sh` genuinely reads the JSON per lens. A lens absent from
the file falls back to `mast`; `tweakreg` is never a safe fallback, because it erases
the dither it is asked to align.

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
    # takes --sample too (default slacs_gold); --all globs that sample's f606W products

Run order: ACS + WFPC2 drizzles → `align_wfpc2_to_acs.py` → `make_cutouts.py`.

## Weight maps: `final_wht_type=ERR` (default), not `EXP`

All four drizzle scripts take `--wht-type {ERR,IVM,EXP}` and write the chosen type
into every AstroDrizzle pass and forward it to the no-CR subprocess. Default is
**`ERR`** for ACS/WFC3/NICMOS, but **`IVM`** for WFPC2 — see the override below.
`cutout_noise.fits` is `1/sqrt(WHT)`, so the weight type *is* the noise model:

- **`ERR`** — full inverse-variance: source Poisson + sky + read + dark. The correct
  per-pixel σ for modelling. On J1143 the core/sky noise ratio is ~4.0 (source shot
  noise present, as it must be).
- **`EXP`** — effective-exposure-time map, **uncalibrated** (no source shot noise): the
  core/sky ratio comes out a flat ~1.04, i.e. it claims the bright core is as noisy as
  blank sky. Do not use for a likelihood.
- **`IVM`** — inverse-variance map. DrizzlePac's own auto-generated IVM (used when no
  file is supplied) is background noise only, no object Poisson — but see the WFPC2
  override below, which supplies a real per-pixel IVM instead.

(DrizzlePac Handbook pp.103,139.) The blank-sky floor is **included** in ERR, not
something to add back. See "Drizzle correlated noise" for the covariance caveat that
ERR does *not* capture.

**WFPC2 override: `drizzle_wfpc2_wf3.py` defaults to `IVM`, not `ERR`.**
`drizzlepac.wfpc2Data.WFPC2InputImage` hardcodes `self.errExt = None` unconditionally
(standard WFPC2 pipeline products never carry an ERR extension), so
`imageObject.buildERRmask()` always takes the "WFPC2 not supported" branch and
silently falls back to **exposure-time-only weighting** — confirmed firing in every
WFPC2 `astrodrizzle.log` ("No ERR weighting will be applied ... WFPC2 data is not
supported by this weighting type"). `--wht-type ERR` on this script is therefore a
trap, not a calibrated noise model: measured core/sky noise ratio was exactly
**1.000** (flat) on J0330-0020 and J1213+6708, vs ~3.5 for ACS F814W on the same
lenses.

Fix: `build_ivm_files()` in `drizzle_wfpc2_wf3.py` builds a per-frame IVM file and
feeds it to AstroDrizzle via a two-column `@`-association file (irafglob's
`atfile_ivm` convention), which DrizzlePac *does* honor for WFPC2.

### The WFPC2 ERR array is not an error array: `ERR == sqrt(SCI)` (2026-07-25)

**Do not build the IVM from the file's ERR extension** — an earlier version of
`build_ivm_files()` did (`IVM = 1/ERR^2`) on the belief that it was "the genuine
Poisson-correct MAST `ERR,3` array". It is not. Verified as an *identity*, not a fit:
`ERR / sqrt(SCI)` = 1.000000 for **100.00%** of good pixels in every frame. That is
Poisson statistics applied to **data numbers as though they were electrons** — the
gain conversion is missing entirely and there is no read-noise term at all.

With SCI in DN (`BUNIT = COUNTS`) and `ATODGAIN = 7.0 e/DN`, it overstates the true
noise by **2.11×** at sky level. Three independent references agree on what the truth
is, on J0330-0020 (sky 17.5 DN):

| | value |
|---|---|
| delivered ERR | 4.184 DN |
| measured, unclipped MAD of blank sky | 1.998 DN |
| measured, adjacent-pixel differences / √2 | 1.981 DN |
| physics, `sqrt(N/g + (RN/g)²)`, RN 5.2 e | 1.747 DN |

The two measurements are model-free and agree to 1%; the physics line is only a
sanity check on them. Anchor any such fix on the **measurement**.

**The model now used is `var_DN = SCI/gain + floor²`**, with the Poisson slope from
the header gain and the additive floor measured **per frame** from its own blank sky
(`measure_noise_floor()`). The floor comes out ~8–9 e, above the nominal 5.2 e read
noise, the excess being dark current plus scattered-light structure. Results on
J0330-0020, blank-sky block-sum test (1.0 = calibrated):

| product | 0.24″ | 0.48″ | 0.96″ |
|---|---|---|---|
| before (`IVM = 1/ERR²`) | 0.43 | 0.51 | 0.53 |
| **after, CR pass** | **0.90** | **1.07** | **1.09** |
| after, no-CR pass | 0.95 | 1.21 | 1.18 |

Peak SNR 20.5 → **53.7**; sky σ fell 2.13× against the 2.11 predicted. The residual
~1.1–1.2 is the same drizzle correlated noise ACS (1.24) and F160W (1.17) show — note
the 0.96″ row rests on only ~60 blocks, so treat it as 1.1 ± 0.1, not a precise
constant. **`K = 1` is correct for WFPC2** and is no longer a placeholder: DrizzlePac's
IVM branch leaves the supplied IVM unscaled but sets `wt_scl = exptime²/scale⁴` (vs
the ERR branch's `1/scale⁴`), so the `exptime²` the ERR mask carries internally is
supplied by `wt_scl` instead and both paths land on a σ in DN/s.

**Two traps this exposed, both worth generalising:**

- **Never select blank sky by pixel *value*.** Clipping at ±2 MAD and then measuring
  the MAD of the survivors biases the width low by ~4%; that is how the error was
  first overstated as 2.28× when it is 2.11×. Select by *position* (a blank box), or
  use a differencing estimator, which needs no selection at all.
- **A model that fits at sky level can still be wrong at the core.** On one lens's
  four frames (sky spanning only 1.2×) the additive-floor, multiplicative, and
  wrong-effective-gain models fit equally well yet diverge ~20% at core brightness.
  Resolved by widening the lever: across **all 92 WF3 frames** the sky spans 7.6–47 DN
  (6×), and there the multiplicative model dies outright — `sqrt(floor)/N` runs 0.116
  to 0.031 where it would have to be constant. A global fit gives
  `var = 0.176·N + 0.536 DN²`, whose intercept is **5.1 e against a nominal read noise
  of 5.2 e**.

  **Do not read that 0.176 as the within-frame slope.** It is a trend *across* fields
  with different scattered light; using it as a signal dependence inside one frame is a
  category error. `var = SCI/gain + floor²` is algebraically
  `(N − N_sky)/gain + var_sky_measured` — physics slope for photons above sky, anchored
  on measurement at sky. The per-frame floor spread (6.2–10.5 e, median 8.5) is why the
  floor is measured per frame and not assumed.

**Verified across the whole sample, not just J0330-0020:** `ERR == sqrt(SCI)` to 0.1%
on **100% of pixels in all 92 frames of all 22 lenses**, `ATODGAIN = 7.0` and
`BUNIT = COUNTS` uniformly. That includes J1218+0830, whose frames come via
`wfpc2_to_flt` from C0M/C1M rather than a MAST FLT — both provenances build the same
bogus ERR. The floor measurement never hit its nominal-read-noise fallback on any
frame, and the floor carries 24–44% of the sky variance (so it is not negligible).

Products carry `IVMMODEL = 'SCI/gain+floor^2'` in the primary header. Three
generations of WFPC2 weight map exist and are indistinguishable by inspection, so
`make_cutouts.py` keys its warning on that keyword; absent means pre-fix. A pre-fix
product **cannot be rescued by a scale factor** — the error is in the noise model,
not the units — it needs a re-drizzle.

### F606W `BUNIT` said `counts` while holding count rates (fixed 2026-07-26)

AstroDrizzle writes count **rates** (it records `D001OUUN = 'cps'`) but does not
rewrite `BUNIT`, which stays at whatever the input carried. For WFPC2 that is
`BUNIT = 'COUNTS'` (DN), so every F606W product claimed counts while holding DN/s —
an **EXPTIME-sized (4400x) error** for anything that reads `BUNIT` to set units. ACS
is unaffected (`ELECTRONS/S`), which is why the discrepancy only shows up when the
bands are compared.

`drizzle_wfpc2_wf3.py` now stamps `BUNIT = 'COUNTS/S'` on the SCI products in the
same loop that stamps `IVMMODEL` (the WHT map is correctly `UNITLESS` and is left
alone). All 92 existing F606W sci/noise files — drizzled and cutouts — were
rewritten in place, header-only, guarded on `D001OUUN == 'cps'`; pixel data verified
byte-identical against `HEAD`.

Note the two optical bands are still in **different unit systems**: F606W in DN/s
(instrumental) against F814W in e/s. That is self-consistent per band because
`PHOTFLAM` is defined per data unit, but any cross-band flux comparison has to go
through `PHOTFLAM` — it cannot compare raw pixel values.

### `1/sqrt(WHT)` is only a sigma map if the input ERR is in counts (2026-07-25)

`make_cutouts.py` turns the weight map into `cutout_noise.fits`. The classic recipe
`sigma = 1/sqrt(WHT)` assumes WHT is a true inverse-variance map, and **whether it is
depends on the units of the calibrated ERR array**, because DrizzlePac computes

    weight = (EXPTIME / ERR)**2          # imageObject.buildERRmask, line ~815

- **ACS FLC** carries SCI *and* ERR in **ELECTRONS**, so `EXPTIME/ERR` is exactly
  `1/sigma_rate` and 1/sqrt(WHT) is a calibrated sigma for the ELECTRONS/S output.
- **WFC3/IR FLT** carries SCI and ERR in **ELECTRONS/S already**, so the same
  expression is `EXPTIME/sigma_rate`: every weight inflated by EXPTIME², and the
  noise map came out a factor **EXPTIME (599.2 s)** too small — SNR ~60,000.

`weight_to_sigma_scale()` in `make_cutouts.py` applies `K = per-frame EXPTIME` for
detectors in `_ERR_IN_RATE_UNITS` (currently `{('WFC3','IR')}`) and `K = 1` otherwise,
and **refuses to guess** if the `D00nDEXP` per-frame exposure times are unequal — the
correction is a single constant only for equal-length frames. All 13 SLACS F160W
lenses are 4 × 599.232 s. `D001WTSC = 1/scale⁴` does **not** enter: it cancels against
the finer output grid, which is why ACS (WTSC 1.0) and WFC3/IR (WTSC 20.9) share one
formula.

**The diagnostic is the blank-sky block-sum test**, not a per-pixel comparison: sum
the background in N×N blocks, compare its scatter to `sqrt(sum(sigma²))` over the same
block, and increase N past the drizzle correlation length. A per-pixel MAD comparison
is *not* usable — drizzle correlation drives it to 0.73 (ACS) / 0.36 (F160W) even when
the map is correct. Measured on J0841+3824:

| band | 0.24″ | 0.96″ | 1.44″ | 1.92″ |
|---|---|---|---|---|
| F814W (ACS, K=1) | 1.04 | 1.24 | 1.24 | 1.28 |
| F160W **before** | 477 | 655 | 725 | 704 |
| F160W **after** (K=599.2) | 0.80 | 1.09 | 1.21 | 1.18 |

Peak SNR went 60,000 → 219, against 225 for F814W on the same lens.

**The residual ~1.2 is real and shared by every band** — it is the drizzle correlated
noise (see the correlated-noise section), i.e. a diagonal-covariance likelihood such
as PyAutoLens understates integrated-flux uncertainties by that much. It differs per
band (1.24 ACS, 1.17 F160W, ~1.1 F606W), so ignoring it also mis-weights bands
*relative to each other* in a joint fit. `make_cutouts.py --corr-factor` applies it;
**default 1.0 (off)**, leaving a pure per-pixel sigma. Both constants are stamped into
the noise FITS header as `NOISEK` and `NOISECOR`.

**WFPC2 F606W is now calibrated too**, by a different fix — its ERR array was never a
real error array. See the `ERR == sqrt(SCI)` section above. `K = 1` for WFPC2 is
correct, not a placeholder.

**Not regenerated**: only **J0841+3824** F160W and **J0330-0020** F606W have
calibrated noise maps. The other 12 WFC3/IR lenses, the other 21 F606W lenses and the
`data/mosaics/` QC grids do not.

### What the calibrated ERR array actually contains (2026-07-25)

Measured on J0841+3824, not taken from the handbooks, by fitting `ERR² = a·SCI + b` in
bins of SCI over the DQ==0 pixels of a single exposure:

| | slope `a` | floor `√b` (read+dark) | sky level | sky Poisson σ |
|---|---|---|---|---|
| ACS/WFC F814W (FLC) | 0.92 | 5.86 e | 59.8 e | 7.73 e |
| WFC3/IR F160W (FLT) | 1.35 | 19.6 e | 466 e | 21.6 e |

**Sky shot noise is included, and dominates the background σ** — it exceeds the
read/dark floor in both bands. ERR is built from total counts *before* sky subtraction;
AstroDrizzle later subtracts `MDRIZSKY` (59.71 e on the ACS frame) from SCI and leaves
ERR alone, which is right — subtracting a constant removes the mean, not the variance.
This is also why the block-sum test lands at ~1.2 and not ~1.6: had the sky term been
dropped, the noise map would be short by that much. Do not "add the sky back".

The ACS slope of 0.92 rather than 1.0 is the flat-field division correlating with the
SCI binning, not a missing term. The IR slope of 1.35 is the up-the-ramp fit's own
read-noise penalty; it falls to 1.17 at the bright end as Poisson takes over.

**Flat-field error is NOT worth adding.** The flat reference files do carry ERR
extensions — fractional error **0.269%** (`qb12257pj_pfl`, ACS F814W) and **0.129%**
(`4ac1921li_pfl`, WFC3/IR F160W) — and it is not resolvable from the data whether
CALACS/CALWF3 propagate them. It does not matter, because at the deflector core:

| | core counts/exposure | Poisson σ | flat stat σ | σ inflation if added |
|---|---|---|---|---|
| F814W | 10,446 e | 102.2 e (0.98%) | 28.1 e | +3.7% |
| F160W | 27,077 e | 164.6 e (0.61%) | 34.9 e | +2.2% |

and it reaches parity with Poisson only at ~138,000 e (ACS) / ~601,000 e (IR) per
exposure, 13× and 22× the core. Two reasons not to add it: it is several times smaller
than the drizzle correlated-noise factor (~1.24) already unaccounted for, and — more
fundamentally — flat error is a **fixed detector pattern, not a random per-pixel draw**.
Four dither positions only partially decorrelate it, so what survives is structured
residual that a diagonal σ map cannot represent; putting it in σ files a systematic into
a statistical slot and misstates the correlation structure to the likelihood.

**At the core the dominant unmodelled term is PSF error, not the flat.** Poisson σ is
0.98% of the signal at the F814W core, so a 1% PSF mismatch already equals the entire
photon-noise budget. If a lens fit shows χ² concentrated on the deflector centre, look
there — not at the noise map.

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
The 24 affected lenses carried `f160W: null` in all three tracking JSONs until the
2026-07-26 full reset emptied those files (see *Tracking JSONs*). This is
reversible: `scripts/stale_scripts/drizzle_nic2.py` is intentionally kept, and re-running it
re-downloads from MAST and re-fetches the CRDS refs automatically. That is exactly
why it can no longer be run by accident — since 2026-07-26 it raises
`NotImplementedError` on import unless `ALLOW_NICMOS=1` is set. Use that override
when the re-download is what you actually want.

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
  WCS (MAST for all 22 lenses per the alignment audit; the two multi-visit lenses are
  split per visit rather than TweakReg-bridged), the blank-sky texture drops to ~4% at
  pixfrac 0.8 and lower at pixfrac 1.0 — comparable to F160W, **not** a 3× outlier.
  See `[[wfpc2_tweakreg_misregisters]]`.
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
conda run -n stenv python scripts/make_cutouts.py --lens J0029-0055 --filt f606W --sample slacs_gold
```

Cuts a square stamp (default 20″) from `data/drizzled/` into `data/cutouts/`,
writing a sci FITS, a noise FITS (from the weight map) and a 3-panel PNG.

**Which pass, and the output prefix.** `--pass {auto,cr,nocrrej}` (default `auto`) picks
the CR pass when one exists, else no-CR. The prefix encodes it so the two never clobber
and can coexist: **`cutout_cr_*` for the CR pass, `cutout_*` for no-CR.** Consequence by
instrument: **WFPC2 F606W always has a CR pass (LACosmic), so `auto` cuts it → the
science stamp is `cutout_cr_*`.** This is fine — LACosmic preserves the core (unlike the
old `driz_cr` route; see `[[acs_cr_pass_eats_core]]`), so its CR pass *is* science-grade.
ACS/WFC3 by default run no CR pass, so `auto` → `cutout_*`. (F160W has no CR pass either
— its FLTs are already ramp-CR-rejected.)

The stamp is recentred on the galaxy rather than the catalogue position, and the peak
search prefers the **CR-rejected** mosaic when one exists — a brightest-pixel search on a
CR-laden no-CR mosaic locks onto cosmic rays. Suppressing CRs purely by widening the
median window needs ~21 pix ≈ 1″ at 0.05″/px, wide enough to bias the peak of a compact
galaxy, so centring uses CR-free data instead. Do not diagnose a bad recentre by reaching
for `--median-size` first — check that a `*_cr_*` mosaic is present.

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

> **Full reset, 2026-07-26.** All three files were emptied to `{}` and every product
> under `data/` was deleted (`calibrated/`, `drizzle_files/`, `drizzled/`, `cutouts/`,
> `mosaics/`, `run_logs/`) as a deliberate clean restart — the accumulated products
> predated too many pipeline fixes to be trusted. Kept: `data/reference_files/` (CRDS
> cache) and `data/pre_drizzled/` (46 MAST-delivered mosaics, not pipeline output).
> The sample was also renamed in the same pass: what these files called `slacs` is now
> **`slacs_gold`**, so any surviving `data/*/slacs/` path is pre-reset.
> **Every "current on-disk state" claim elsewhere in this file is therefore stale** —
> including the two *Not regenerated* notes and the F160W-`null` note below, which
> describe a tracking state that no longer exists. The *reasoning* in those sections
> is still valid and is why the reruns must be done with the current scripts; only the
> inventory is void.

Three files are updated automatically by every script run:

- **`lens_products.json`** — `{lens: {filter: [obsid, ...]}}` — HST rootnames downloaded from MAST
- **`lens_instrument.json`** — `{lens: {filter: "INSTRUME/DETECTOR"}}` — e.g. `"ACS/WFC"`, `"WFPC2/WF3"`, `"WFC3/IR"`
- **`lens_exptime.json`** — `{lens: {filter: exptime_seconds}}` — from the CR-rejected drizzle header

If a lens has no data for a filter, the value is stored as `null`.

**The key is the product directory, not the filter (2026-07-26).** For nearly every
lens those coincide (`f606W`), but a **split-visit lens is keyed per visit** —
`f606W_v1` / `f606W_v2`, matching `data/drizzled/<lens>/<key>/` — and carries **no
bare `f606W` key at all**. So the old invariant ("all three files carry a key for
every band on every lens; a missing key means the file is out of sync") no longer
holds for J0728+3835 and J0822+2652. Check keys against the product directories
instead.

That change fixed a live error. Both split lenses had been recording a combined
`f606W`: 6 obsids and **6600 s** for each, describing a product that has never
existed on disk. J0728+3835's real product is `f606W_v2` alone — 4400 s, 4 frames —
its 2-frame visit having been dropped for want of dither phase, and J0822+2652 has
two products (2200 s + 4400 s). **The failure was invisible because the key was
present and plausible**, which is the opposite of the "missing key = out of sync"
check above and the reason that check is not sufficient on its own. Root cause:
`drizzle_wfpc2_wf3.py` keyed its JSON writes on the bare filter while writing
products to a `--out-suffix` directory. It now keys on `product_key = filt +
out_suffix`.

**`lens_products.json` records the frames that reached the drizzle, not the whole
download.** All three scripts now write it after the input list is settled, so a
re-run on already-downloaded data refreshes it (it used to sit inside the download
block and never update). The two sets differ for two reasons:

- **WFPC2**: `--pa` selects one visit out of two, and a lens that exits for want of
  dither phase writes nothing at all — correct, since no product exists.
- **ACS/WFC3**: AstroDrizzle silently drops `EXPTIME=0` frames, so recording the
  download overstated four entries — J0008-0004 f814W (4 recorded vs 3 drizzled),
  J0912+0029 f555W (8 vs 5) and f814W (8 vs 7), J1213+6708 f814W (5 vs 4). Those
  frames are now excluded from the record but **not deleted**, unlike WFPC2 where
  `MIN_EXPTIME` removes them outright because they would otherwise reach the drizzle.

Exposure times were correct throughout — they come from the drizzle header.

When auditing obsid counts against `NDRIZIM`, note **ACS/WFC FLCs are 2-chip MEFs**,
so `NDRIZIM = 2 x exposures` there and 1× for WFPC2/WFC3. Comparing them directly
reports 54 false mismatches.

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

**`info/lens_samples.json` is the single source of truth** for which lenses are in which
sample, and for the per-lens MAST quirks (`mast_target`, `force_copy`). Read it only
through `scripts/mast_target_names.py` — never parse it elsewhere and never keep a second
lens list in a script or runner. `info/list_of_lenses.txt` was exactly such a second copy
and was deleted on 2026-07-26.

```bash
conda run -n stenv python scripts/mast_target_names.py --list        # samples + sizes
conda run -n stenv python scripts/mast_target_names.py slacs_gold    # lens names
```

| Sample | Lenses | What it is |
|---|---|---|
| **`slacs_gold`** | 38 | The lenses reduced so far — the working sample and **the default `--sample` of every script** (`mast_target_names.DEFAULT_SAMPLE`) |
| **`slacs_other`** | 93 | The rest of SLACS from Bolton et al. 2008 Table 4 (`info/slacs_coords.py`, 131 total). Not yet reduced |
| **`gallery`** | 16 | BELLS GALLERY, HST proposals 14189, 16734. WFC3/UVIS multi-band (F225W, F275W, F438W, F606W, F814W). Reduction scripts not yet written |

`slacs_gold` band coverage, as of the last MAST survey:

| Band | Instrument | Lenses | Notes |
|---|---|---|---|
| F814W | ACS/WFC | 38 | all |
| F606W | WFPC2/WF3 | 22 | the other 16 have no WFPC2 data in any filter |
| F555W | ACS/WFC | 16 | exactly the 16 without F606W (props 10494/10798) |
| F160W | WFC3/IR | 13 | all proposal 11202, ~2397 s each |
| F160W | NICMOS/NIC2 | 24 | deprioritised — data deleted, see *NICMOS is deprioritised* |

So every `slacs_gold` lens has F814W, and the sample splits cleanly into 22 with WFPC2
F606W and 16 with ACS F555W. HST proposals: 10886, 11202, 10494, 10798.

**Two caveats on the other two samples.** `slacs_other` has **no `GAL-*` names surveyed**,
so for those lenses "no observations" may only mean the lens is archived under a
designation nobody has looked up — weaker evidence than the same result on `slacs_gold`,
where all 14 GAL names are known. And `gallery` lenses are **not in `info/slacs_coords.py`**
(it is SLACS only), so they get no common output WCS: the drizzle scripts fall back to the
native WCS with a printed warning. Both need work before those samples are reduced.

**Every lens is tried on every run; only the ones with data are downloaded.** A lens with
no data for an instrument+filter is an ordinary outcome, not a failure — the drizzle
script records `null` in the tracking JSONs, prints a line beginning `=== NO DATA:` and
**exits 0**. The batch runners count those separately from failures. This is why the
runners iterate the sample roster rather than globbing `data/calibrated/`: a glob only
re-runs what is already on disk, so it does nothing after a wipe and never picks up a
newly added lens. A genuine download error still exits non-zero — the two are kept apart
by `mast_target_names.NoMastData`, a dedicated exception the download block's broad
`except Exception` cannot swallow. Getting that wrong would silently record a network
failure as "this lens has no data".

Target names on MAST follow the pattern `SDSS<LENS>` (e.g. `SDSSJ0008-0004`). The MAST
query uses `target_name=f'SDSS{lens}%'` with a wildcard to handle minor naming variations.

**COPY handling differs by instrument, deliberately:**

- **ACS** (`drizzle_acs_wfc.py`) filters COPY observations out in favour of non-COPY when both exist — **except** for lenses carrying `"force_copy": true` in `info/lens_samples.json`, whose non-COPY observations are unusable. Currently only `J1032+5322` (F814W): the non-COPY frames are `EXPTIME=0`, so `mast_target_names.force_copy(lens)` selects the COPY frames. J1032+5322 is the only ACS lens in `slacs_gold` with COPY data at all — surveyed across F814W and F555W, nothing else has any.
- **WFPC2** (`drizzle_wfpc2_wf3.py`) keeps **both**, because there the COPY sets are genuine repeat visits carrying most of the usable exposure time (see above). It rejects junk on `MIN_EXPTIME` instead of on target name.

Do not "unify" these two policies without re-checking the archive — they encode different facts about different datasets.

### Non-standard MAST target names (`GAL-*`)

Some lenses are **not** on MAST under an `SDSS<LENS>` name — they use a `GAL-<plate>-<mjd>-<fiber>` designation instead, so the default `SDSS{lens}%` query returns no observations. The drizzle scripts resolve this automatically via `scripts/mast_target_names.py`: for a lens carrying a `mast_target` in `info/lens_samples.json` they query the `GAL-*` name first and fall back to `SDSS{lens}%`. **Output/directory names always stay in the J convention** (the left column) regardless of the MAST name used to fetch.

The table below is a **copy for reading**; the live values are the `mast_target` entries
in `info/lens_samples.json`. Edit the JSON, not this table. All 14 are in `slacs_gold` —
no `GAL-*` names have been surveyed for `slacs_other` or `gallery`, which is why a
no-data result on those samples is not conclusive.

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
- **CR pass**: cosmic-ray-cleaned science image. Default method is **LACosmic** — mask CRs
  per frame, then a plain-mean drizzle (`median=False, blot=False, driz_cr=False`,
  `resetbits=0`). `--cr-method drizcr` restores the old AstroDrizzle route
  (`median=True, blot=True, driz_cr=True`), which eats the core — see the CR-rejection
  section.
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
- WFC3/IR: `512` only (IR flat-field "blobs") — see the quadrupled-defects note below.
  Write it as `'512'`, **never** `''`
- NICMOS: `2,4,8` (uncertain linearity/dark/flat corrections — acceptable calibration imperfections; all defects, saturation, CRs excluded. Per NICMOS Data Handbook Table 2.3)

### WFC3/IR quadrupled defects, and the `bits=''` trap (2026-07-25)

The user reported F160W noise maps where "each spot turned into 4". Two rounds of
diagnosis; the first was wrong, and its "fix" made things worse. Both are recorded
because the failure mode is the transferable lesson.

**`final_bits=''` / `driz_sep_bits=''` disables DQ masking entirely.**
`astropy.nddata.bitmask.interpret_bit_flags('')` returns `None`, and AstroDrizzle then
logs `bits : None` and keeps *every* flagged pixel as good — the exact opposite of the
intent. Measured: of 32,210 DQ-flagged pixels per frame only 3,328 ended up masked, and
those came from the static mask, not DQ. It also silently voids any CR flag written into
DQ 4096 — 97.7% of a LACosmic run's flags were discarded this way, which is why an
earlier "LACosmic doesn't help / `driz_cr` doesn't help" conclusion was meaningless.
**`0` is the spelling that means "reject everything flagged"; `''` and `None` mean
"reject nothing".** `drizzle_wfc3_ir.py` now asserts this at import.

**Why the defects quadruple.** The no-CR pass (the only pass F160W uses) runs
`median=False, blot=False, driz_cr=False` — no cross-frame rejection at all. A
detector-fixed defect kept as "good" is drizzled as-is at each dither position, and the
standard WFC3-IR-DITHER-BOX has 4, so one detector pixel becomes 4 sky replicas.
Confirmed on J0841+3824: 39 single-frame defects in the 20″ stamp mapped back to only
**15 distinct detector pixels**, each appearing 3–5 times, and 38 of 39 carried
**DQ = 48** (16 hot + 32 unstable).

**The right bit list is `'512'`,** measured rather than assumed:

| DQ class | n/frame | median \|residual\| | 90th pct | how it behaves |
|---|---|---|---|---|
| DQ==0 (clean) | 996k | 0.67 MAD | 1.7 | reference |
| 16 hot | 7,140 | 1.93 MAD | 16.1 | **bright — must mask** |
| 32 unstable | 14,434 | 2.26 MAD | 62.4 | **bright — must mask** |
| 512 blob | 12,480 | 0.70 MAD | 1.8 | indistinguishable from clean |
| 64 warm | **0** | — | — | never set in these FLTs |

So 512 is kept as good (free — masking it would raise zero-coverage sky pixels from
0.19% to 0.77%), 64 is dropped from the string because it never fires, and everything
else is rejected. Result on J0841+3824, fresh whole-field census outside r>3″:
single-frame defects **28% → 5%** of detected peaks (42% → 9% among the brighter half),
with all 67 real sources retained.

**Do not add CR rejection to F160W.** Both routes damage point sources, because at
0.1283″/px the IR PSF is ~1 px FWHM and looks like an outlier:

| route | defects left | field star | deflector core |
|---|---|---|---|
| no CR (current) | 5 of 92 peaks | intact | 1.000 |
| `driz_cr` | 1 of 87 | peak −10% | F(1″) 0.983 |
| LACosmic | 0 of 96 | **zeroed** (59 px at WHT=0) | F(1″) 1.000 |

`objlim` does not rescue LACosmic — at 5/15/20/30/50/100 it still clips the star, while
never catching more than 12% of the unflagged >5σ outliers it exists to remove.
`--cr-method {lacosmic,drizcr}` exists on the script for comparison runs only.

**The noise map still shows the quadruplets, and that is now correct.** Masking a defect
in 1 of 4 frames leaves 3 frames of coverage there, so σ rises by √(4/3) = 15.5%. The
dots are an honest statement that those pixels have 3 frames instead of 4, not injected
flux. Before the fix the same defects were contaminating the *science* image; now they
are confined to the weight map. Nothing short of more dither positions removes the
genuine ones, and inflating the weights to hide them would be a lie to the likelihood.

### `--dq-refine`: most of the inherited flags are not justified by the data

The speckle above was still much heavier than it needed to be, because WFC3/IR hot and
unstable flags are inherited from the **dark reference file**, which characterises a
pixel across a whole anneal cycle — a pixel that misbehaved once is flagged in every
exposure of the cycle. Measured on J0841+3824, fraction of pixels whose residual from a
5×5 local median exceeds 3 MAD *in the same exposure that flags them* (unflagged pixels
sit at 1.2% for scale):

| DQ | n/frame | actually deviant |
|---|---|---|
| 16 hot | 7,140 | 37% |
| 32 unstable | 14,434 | 43% |
| 8 deviant zero-read | 3,140 | 35% |
| 4 bad detector px | 3,815 | 68% — **keep**, a permanent defect |
| 512 blob | 12,480 | 1.9% — already kept as good |

So ~60% of what was being masked was indistinguishable from clean sky in the very frame
where it was discarded. `refine_dq_flags()` in `drizzle_wfc3_ir.py` clears bits **8, 16,
32** (`_SOFT_DQ`) on pixels that are not deviant at `--dq-refine` σ (default 3.0; 0
disables), per exposure, so a pixel stays masked in the frames where it actually
misbehaves. It errs conservative on real sources — a flagged pixel on a steep PSF has a
large local-median residual and keeps its flag. It runs on the **copies in
`data/drizzle_files/`**, never on `data/calibrated/`. It also clears **DQ 4096**, which
is not a calibration flag at all: driz_cr/LACosmic write it and it persists in the file,
so a rejected experiment's flags leak into later runs (5,730 px/frame were still present
from the abandoned LACosmic tests).

Measured effect: masked pixels **1.96% → 0.82%** per frame, and in the 20″ stamp the
blank-sky σ excess >5% went **3.88% → 1.39%** of blank-sky pixels, in **480 → 190
clumps**. This is the floor for the approach: the pixels still masked are genuinely bad
(median deviation **9.2σ**, 81% >3σ, 22% >30σ), and with 4 dither positions each one
necessarily leaves a 3-of-4 coverage mark.

**Two framings to avoid when judging this.** (1) The weight map is a **continuum, not
quantised** into 4/3/2/1 frames — binning `WHT/plateau` into coverage steps overstates
the damage badly (it read "26% of the mosaic degraded", which is not a meaningful
number). (2) Most of the σ structure in the stamp is **legitimate ERR weighting**: with
`final_wht_type=ERR` the weight correctly drops wherever there is source flux, so ~48%
of on-source pixels carry >5% σ excess and should. Judge this by the **blank-sky** clump
count only, with on-source pixels masked out — that is the number that isolates coverage
deficits from correct inverse-variance behaviour.

### Why F814W shows no "quadrupling" (it does — as stripes)

A natural question once the F160W dots are understood, and the answer is **not** that ACS
has cleaner data. Measured on J0841+3824, same frame count (4), same metric:

| | F814W (ACS/WFC) | F160W (WFC3/IR) |
|---|---|---|
| masked per frame | **1.58%** | 0.81% |
| footprint of one masked input px | **1.00 output px** | **4.57 output px** |
| stable hot px (16) | 41,150 — **kept as good** | masked |
| warm px (64) | 49,222 — **kept as good** | never set |
| dominant masked population | CRs, 72,162 (random per frame) | hot/unstable (detector-fixed) |
| detector-fixed masked population | bad columns, 52,422 | hot/unstable |

Four causes, only one of which is about the defects:

1. **No footprint amplification.** ACS drizzles native 0.05″ → 0.05″, so a masked input
   pixel marks exactly one output pixel. F160W's 0.1283″ → 0.06″ wipes 4.57. This is what
   turns an isolated bad pixel into a visible *blob* instead of an invisible single pixel,
   and it is the single biggest difference between the bands.
2. **ACS deliberately keeps the population that would quadruple.** `final_bits='256,64,16'`
   keeps stable-hot and warm pixels as good — exactly what dominates F160W's speckle. The
   opposite trade: let them into the science image rather than punch coverage holes. On a
   CCD that is defensible (stable, dark-corrected); WFC3/IR's "unstable" flag by definition
   means the pixel is not reproducible. **Do not "unify" the two bit lists.**
3. **What ACS masks is mostly cosmic rays** (72,162/frame, 27% of its masked total). CRs
   land at random positions per frame, so they never replicate — they scatter as
   single-pixel marks.
4. **The one detector-fixed thing ACS masks is bad columns**, and those *do* quadruple.
   With `final_rot=0.0` the detector columns rotate by the roll angle, so the 4 dither
   replicas appear as **4 parallel diagonal stripes** in the noise map — same physics as
   the F160W dots, rendered as lines because the defect is a line. See the ACS bad-column
   note; this is what those stripes are.

**After `--dq-refine`, F160W has the more uniform blank sky of the two.** σ excess on
blank sky: F814W 50th −4.4%, 90th **+13.1%**, 99th **+22.8%**; F160W 50th −3.7%, 90th
−0.5%, 99th +6.8%. F814W's larger blank-sky excess by *area* (22.65% vs 1.39%) is one
6,178-px connected region (the stripes) plus a haze of 2-px CR marks, against F160W
clumps that top out at 17 px. The F160W dots are more *legible*, not more numerous — so
do not rank the two bands' noise maps by eye.

### Scale/pixfrac scan: 0.06″ + pixfrac 1.0 is right, and coarsening does not help

Run on J0841+3824's DQ-refined frames, drizzling the same exposures onto six grids and
cutting the same 20″ field from each. **Speckle area is in arcsec², never pixels** —
pixel counts are not comparable across output scales, which is the trap in this whole
comparison:

| scale | pixfrac | speckle area | clumps | texture | adj corr | PSF FWHM |
|---|---|---|---|---|---|---|
| **0.06″** | **1.0** (current) | 2.69 □″ | 195 | 4.81% | 63.1% | **0.253″** |
| 0.06″ | 0.8 | 5.32 □″ | 931 | 6.94% | 53.5% | 0.235″ |
| 0.08″ | 1.0 | 2.20 □″ | 171 | 4.74% | 55.2% | 0.264″ |
| 0.08″ | 0.8 | 7.76 □″ | 1,023 | 6.22% | 44.2% | 0.247″ |
| 0.10″ | 1.0 | 1.81 □″ | 132 | 4.68% | 49.1% | 0.293″ |
| 0.1283″ (native) | 1.0 | 1.35 □″ | 75 | 4.63% | 49.8% | 0.310″ |

**Coarsening the grid does not remove defects.** The sky area a masked pixel costs is
set by the **input** pixel size (0.1283″); a coarser output grid just renders the same
patch in fewer, bigger pixels. Hence 0.06″ → 0.08″ cuts the affected area only **18%**,
not the ~43% a naive `(0.1283/scale)²` argument predicts. That ratio describes pixel
*count* — it is why F160W defects read as blobs and ACS's as single pixels — and it must
not be used to argue for a scale change.

**pixfrac 0.8 is decisively wrong**, and this is much stronger evidence than the original
stacked-FWHM argument: it doubles the speckle area at 0.06″ (2.69 → 5.32 □″) and triples
it at 0.08″ (2.20 → 7.76 □″), with clump counts 195 → 931 and 171 → 1,023, and it puts a
visible waffle texture into the noise map. It does lower adjacent-pixel correlation by
~10 points at fixed scale — the known trade — but nowhere near enough to justify that.

**The resolution cost is real**: FWHM 0.253″ → 0.264″ at 0.08″ (+4%) → 0.310″ at native
(+22%, and visibly pixelated — 1.2 px per FWHM is undersampled for PSF convolution).
So 0.08″ is defensible (4% resolution for 18% less speckle) but marginal, and native is
not. **Keep 0.06″ / pixfrac 1.0** — settled by the user on 2026-07-26 after reviewing
this scan. It is a decision, not a provisional default: do not reopen it or re-drizzle
F160W to another scale without being asked.

Caveat on the `adj corr` column: cross-scale values are not strictly comparable, because
the high-pass filter used to isolate the noise spans a different physical size at each
scale. The fixed-scale pixfrac comparison is clean.

**Not regenerated**: only **J0841+3824** carries post-fix products. The other 12 WFC3/IR
F160W lenses (and the `data/mosaics/` QC grids) still have pre-fix products.
