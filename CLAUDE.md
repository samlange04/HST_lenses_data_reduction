# CLAUDE.md

Guidance for Claude Code working in this repo. This file is the **operational**
reference: commands, layout, decisions, and the traps that silently produce wrong
products. The *evidence* behind each decision (measurement tables, derivations) lives in
the persistent memory store — pointers below read `→ memory: <slug>`. Don't re-derive
what a pointer already settles; don't restate the tables here.

> **CLAUDE.md documents intent, not always state.** This doc has repeatedly described
> decisions the code never implemented (argparse defaults, runners reading a JSON, a
> pipeline stage). Before running any batch campaign or reporting what's outstanding,
> verify against source — argparse defaults, what the runner actually passes, tracking
> JSONs vs what's on disk. → memory: claude-md-documents-intent-not-state

## Environment

All scripts run inside the `stenv` conda environment (STScI pipeline stack):

```bash
conda run -n stenv python scripts/<script>.py --lens <LENS> --filt <FILTER>
```

Interactive: `conda activate stenv`.

## macOS write-hang workaround (required — keep this wiring)

On this Mac, AstroDrizzle's large buffered FITS writes hit a kernel lost-wakeup
(`tofile → write() → cluster_write → copyin → lck_rw_sleep`) that wedges the process into
an unkillable U-state (only a reboot clears it). `scripts/mmap_fits_write.py` monkeypatches
`astropy.io.fits.file._array_to_file` to write via `mmap`+`memcpy` (the `vm_fault` path),
dodging the hang; it is a no-op off macOS and byte-identical to stock astropy. Every
drizzle script imports it and calls `install()` before AstroDrizzle. `DRIZZLE_MMAP_DEBUG=1`
logs each mmap write. **`num_cores=1` in all scripts is a related, separate requirement** —
parallel `fork` triggers the same U-state. → memory: env_rosetta_x86, feedback_num_cores

## What this repo does

HST image reduction for gravitational-lens samples (SLACS, BELLS). Per lens+filter:
1. Download calibrated exposures from MAST
2. Download CRDS reference files
3. Align + combine with `AstroDrizzle`. All bands trust the delivered MAST WCS
   (`--align mast`) by default: ACS/WFC3 skip `updatewcs` and TweakReg; WFPC2 runs
   `updatewcs` (for the distortion arrays) but **not** TweakReg. TweakReg is opt-in
   (`--align tweakreg`) and used by no lens — see *WCS alignment*.
4. Produce a CR-rejected mosaic (LACosmic) — now the default and the science product for
   ACS and WFPC2. The no-CR-rejection mosaic is opt-in (`--nocrrej`) for comparison only.
   WFC3/IR F160W has no CR pass (see *Cosmic-ray rejection*).
5. Update three JSON tracking files in `info/`.

## Running a single lens

`--sample` defaults to **`slacs_gold`** everywhere; it sets the `<sample>` level of every
`data/` path, so a wrong value silently writes a correct product into the wrong tree.

```bash
conda run -n stenv python scripts/drizzle_wfpc2_wf3.py --lens J0008-0004 --filt f606W
conda run -n stenv python scripts/drizzle_acs_wfc.py  --lens J0008-0004 --filt f814W
conda run -n stenv python scripts/drizzle_wfc3_ir.py  --lens J0008-0004 --filt f160W
conda run -n stenv python scripts/drizzle_acs_wfc.py  --lens J0216-0813 --filt f555W
```

Scripts are **idempotent**: they skip MAST download if calibrated files exist, and skip
the whole drizzle if the final output already exists. To force a re-run, delete the lens's
dir under `data/drizzled/` (and `data/drizzle_files/`).

A lens with no data for the requested instrument+filter prints `=== NO DATA: ...`, records
`null` in the tracking JSONs, and **exits 0** — see *Lens Samples*.

The **CR-rejection pass (LACosmic) is the default** for ACS and WFPC2, and is the product
downstream reads. The no-CR pass is opt-in for comparison via `--nocrrej` (`--no-nocrrej`
is the default). ACS also accepts `--no-cr` to skip CR (e.g. with `--nocrrej` for a no-CR
only run). WFC3/IR F160W has no CR pass at all — do not add one (see *Cosmic-ray
rejection*).

## Running all lenses

Each runner takes an optional sample arg, defaulting to `slacs_gold`:

```bash
bash scripts/run_acs_all.sh                  # ACS/WFC F814W + F555W
bash scripts/run_wfc3_all.sh                 # WFC3/IR F160W
bash scripts/run_wfpc2_wf3.sh                # WFPC2/WF3 F606W: drizzle -> align -> cutout
bash scripts/run_gallery_uvis_all.sh         # WFC3/UVIS F225W/F275W/F438W/F606W/F814W (gallery only)
bash scripts/run_cutouts_all.sh              # stamps for whatever products exist
bash scripts/run_psf_all.sh                  # PSF kernels for whatever products exist
bash scripts/run_acs_all.sh slacs_other      # any runner, any sample
```

All runners take the roster from `info/lens_samples.json` via `scripts/mast_target_names.py`
— **except `run_cutouts_all.sh` and `run_psf_all.sh`, which glob `data/drizzled/`** (a stamp
or PSF needs a mosaic that exists). They report `ok` / `no data` / `FAILED` separately, so
the 16 `slacs_gold` lenses with no WFPC2 data aren't mistaken for errors. `run_cutouts_all.sh`
globs `<filt>*` (not `<filt>`) so per-visit split-visit dirs are included, and covers every
sample's filters (SLACS `f606W f814W f555W f160W` plus gallery's UV/blue bands `f438W f275W
f225W`) in one runner. `run_gallery_uvis_all.sh` defaults to `gallery`, the only sample with
WFC3/UVIS data, but (like the other runners) accepts any sample as its first arg.

### `run_wfpc2_wf3.sh` — the single WFPC2 driver (three traps it exists to avoid)

It runs the full three-stage order (`drizzle_wfpc2_wf3.py` → `align_wfpc2_to_acs.py` →
`make_cutouts.py`) for every lens, reads per-lens alignment from
`info/wfpc2_alignment.json`, expands the two split-visit lenses into per-visit products,
and retries failures once. It carries **no exclusion list** — the drizzle script measures
each lens's dither coverage and skips any that can't reach 0.05″/px, reported as
`SKIPPED (dither phase)`, not a failure. It stops before align+cutout on a `no data` lens.

It was rebuilt (2026-07-26) because a prior version drifted into quietly-wrong products.
The three traps generalise to anything driving this pipeline:

- **`--align` default matters.** The old runner passed no `--align`, taking the script
  default — which was still `tweakreg`, the mode the audit rejected for all 22 lenses. The
  script default is now `mast`; the runner reads the JSON per lens; a lens absent from the
  file falls back to `mast`. `tweakreg` is never a safe fallback — it erases the dither it
  is asked to align.
- **Split-visit handling.** J0728+3835 and J0822+2652 must be drizzled per visit, not as a
  combined dataset across a ~15° roll — else the stack smears *and* the tracking JSON keys
  get rewritten to a bogus combined `f606W`.
- **The align step is not optional.** A re-drizzle discards the astrometric tie (it's a
  `CRVAL1/2` edit on the drizzled product), so `align_wfpc2_to_acs.py` must run after
  *every* drizzle. Skipping it gives stamps that look perfect alone and are ~0.3–0.9″ off
  the other bands.

→ memory: claude-md-documents-intent-not-state (all three were once documented but
unimplemented).

### `scripts/stale_scripts/` — retained, but nothing invokes them

- `drizzle_wfpc2_pc.py` — **superseded**: extracts the wrong chip *and* `rmtree`s the good
  WF3 products. Override `ALLOW_SUPERSEDED_WFPC2_PC=1`.
- `drizzle_nic2.py` — **deprioritised**: an accidental run re-downloads ~472 MB and
  repopulates NICMOS entries. Override `ALLOW_NICMOS=1`.
- `run_all_lenses.sh` — retired WFPC2 driver; refuses to run (its retry pass moved into
  `run_wfpc2_wf3.sh`).
- `rebuild_stenv_arm64.sh` — one-off, never used; the write hang was fixed in-process by
  `mmap_fits_write.py`, not by a native-arm64 rebuild.

The first two also **raise `NotImplementedError` on import** — a deliberate guard: it fails
loudly with a non-zero status a batch runner can't mistake for a clean skip, before any
network call.

### Total-exposure-time gate

All four drizzle scripts (`drizzle_acs_wfc.py`, `drizzle_wfc3_ir.py`,
`drizzle_wfpc2_wf3.py`, `drizzle_wfc3_uvis.py`) sum `EXPTIME` over the frames that would
actually reach the drizzle (post `EXPTIME=0`/`MIN_EXPTIME` filtering, post `--pa` visit
selection for WFPC2) and gate on the total before doing the expensive drizzle work — added
because `slacs_other` runs generally shorter total exposures than `slacs_gold`.

- **`BLOCK_EXPTIME = 500s`** — no product is written. Same outcome/shape as no MAST
  data: tracking JSONs get `null`, the script prints `=== BLOCKED (exptime): ... ===`
  and exits 0, so a batch runner counts it separately from a failure (`run_acs_all.sh` /
  `run_wfc3_all.sh` / `run_gallery_uvis_all.sh` track it in a `blocked` counter;
  `run_wfpc2_wf3.sh` reports `blocked (exptime)` and skips the align/cutout stages, same
  as `no data`).
- **`WARN_EXPTIME = 1200s`** — the drizzle proceeds; the script prints
  `  EXPTIME WARNING: ...` and the batch runners report `OK (low exptime)`.

No current `slacs_gold` or `gallery` product falls under either threshold. The gate has
fired for real in `slacs_other`: its F814W visits are mostly Bolton-era legacy exposures,
and **16 of the 27 `slacs_other` lenses are `BLOCK_EXPTIME`-gated at F814W** (3 succeed, 8
have no ACS data at all) — the first sample where this isn't a no-op. F606W (24/27),
F160W (6/27) and gallery's five UVIS bands cleared the gate everywhere they had data.

## Data flow and directory layout

```
data/
  calibrated/<sample>/<lens>/<filter>/    ← downloaded FLT/FLC/CAL files
  drizzle_files/<sample>/<lens>/<filter>/ ← working dir; AstroDrizzle runs here (run.log, shift_*.txt, *_single_*.fits, *.png)
  drizzled/<sample>/<lens>/<filter>/      ← final products (<prefix>_cr_*/_nocrrej_* sci+wht)
  cutouts/<sample>/<lens>/<filter>/       ← cutout_sci.fits / cutout_noise.fits / cutout.png / cutout_[cr_]psf.fits
  psf/<sample>/<lens>/<filter>/           ← archival PSF products (psf_kernel.fits / psf_epsf.fits / psf.png)
  mosaics/<sample>/                       ← QC mosaics tiling every lens's cutouts/PSFs (make_mosaics.py, make_psf_mosaics.py)
  pre_drizzled/                           ← 46 MAST-delivered mosaics, kept for reference; not pipeline output
  run_logs/                               ← per-lens batch-runner logs
  reference_files/                        ← CRDS reference files (auto-downloaded once)
```

## Instrument-specific scripts

| Script | Input | MAST product | Ref env | Pixel scale | Suffix |
|---|---|---|---|---|---|
| `drizzle_wfpc2_wf3.py` | `u*flt.fits` | FLT / CALWFPC2 | `uref` | 0.0996″ → 0.05″ | `_drw_` |
| `drizzle_acs_wfc.py`  | `*flc.fits`  | FLC / CALACS   | `jref` | 0.05″ | `_drc_` |
| `drizzle_wfc3_ir.py`  | `*flt.fits`  | FLT / CALWF3   | `iref` | 0.1283″ → 0.06″ | `_drz_` |
| `drizzle_wfc3_uvis.py` | `*flc.fits` | FLC / CALWF3   | `iref` | 0.0396″ (native) | `_drc_` |
| `drizzle_nic2.py`     | `*cal.fits`  | CAL / CALNIC   | `nref` | 0.0756″ | `_drz_` |

The output suffix is set by **input file type**, not output name: `_drc_` for FLC (ACS,
WFC3/UVIS), `_drw_` for WFPC2 FLT, `_drz_` for everything else. The WFPC2 script extracts
only the WF3 chip (SCI/ERR/DQ ext 3) into `wf3_`-prefixed files first; the others are MEF
files DrizzlePac handles natively. → memory: instrument_drizzle_ref, crds_bestrefs_always_run
(never skip `bestrefs` when the ref dir is non-empty).

`drizzle_wfc3_uvis.py` is the BELLS GALLERY driver — see *BELLS GALLERY: WFC3/UVIS
reduction* below for its defaults, alignment, and current coverage.

### WFPC2: the lens is on WF3, not the PC

For all 22 SLACS WFPC2 F606W lenses the lens galaxy falls on **WF3** (ext 3) at
~(435, 424). `DETECTOR = PC` in the primary header (hence the `WFPC2/PC` MAST label) names
the *aperture*, not the chip; the full-field aperture centres the target on WF3. Only the
per-extension `DETECTOR` identifies chips (1=PC…4=WF4). The superseded `drizzle_wfpc2_pc.py`
extracted ext 1 → 22 blank-sky mosaics ~79″ off the lens. → memory: wfpc2_target_on_wf3

Two consequences for `drizzle_wfpc2_wf3.py`:
- **Chip renumbering.** DrizzlePac indexes chips positionally `(SCI, 1..N)`, so the
  extracted WF3 ext must be rewritten to `EXTVER=1` or `WFPC2InputImage` raises
  `KeyError: Extension ('SCI', 1) not found`. `detnum` still comes from `DETECTOR` (stays
  3), so the WF3 gain/readnoise row is unaffected. DQ bits `8,1024` carry over unchanged.
- **Output scale.** WF3 is 0.0996″/px; the `WFPC2-BOX` pattern half-pixel-dithers both
  axes (offsets in the WCS — `POSTARG1/2` are zero), supporting 0.05″/px at pixfrac 1.0.
  `dither_phase_counts()` measures phase coverage at runtime and the script **exits without
  writing** if either axis has <2 distinct phases. No hardcoded exclusion list — coverage
  depends on the MAST query (J0728+3835 looked unusable at 2 exposures, is fine at 6).

### WFPC2 archive traps (each silently costs exposures)

→ memory: wfpc2_copy_visits_and_c0m

- **`-COPY` targets are genuine repeat visits**, not duplicates — they carry most of the
  usable exposure time. The script keeps both and filters on `MIN_EXPTIME = 10s` instead of
  a "prefer non-COPY" rule (which also covers the `EXPTIME=0` case).
- **`WFPC2/WFC`-labelled obs usually ship no FLT** — only raw `C0M`+`C1M`. The script
  downloads them and converts via `drizzlepac.wfpc2Data.wfpc2_to_flt`. J1218+0830 is the
  one lens whose extra WFC frames are real 1100s science; elsewhere they're 0.5s check
  shots dropped by `MIN_EXPTIME`.
- **Multi-visit lenses are split, not TweakReg-combined** — see *WCS alignment*.

## WCS alignment: `--align`, and why it differs by instrument

**Do not unify the scripts on this.** Each takes `--align {mast,tweakreg}`; the correct
default is not the same for all. Verified by stacked FWHM and (for WFPC2) core-registration
scatter. → memory: wfpc2_tweakreg_misregisters

| Instrument | Default | What `mast` does |
|---|---|---|
| ACS/WFC | `mast` | no `updatewcs`, no TweakReg |
| WFC3/IR | `mast` | no `updatewcs`, no TweakReg |
| WFC3/UVIS (gallery) | `mast` (no per-lens audit yet) | no `updatewcs`, no TweakReg |
| WFPC2/WF3 | `mast` (per-lens audit) | `updatewcs(use_db=True)`, no TweakReg |

**Why ACS/WFC3 must not re-solve.** MAST delivers them fitted to GSC 2.4.2 / GAIA eDR3
(`WCSNAME = *-FIT_REL_GSC242/-GAIAeDR3`), relative astrometry good to ~0.05–0.8 px.
TweakReg aligns every frame onto the *first*, so on dithered exposures it measures the
dither as an error and removes it — the reported `XSH/YSH` come out equal to the POSTARG
offsets, and AstroDrizzle then stacks dithered frames as if they shared a pointing. **This
is what smears point sources and splits lensed arcs into offset copies.**

**WFPC2 needs `updatewcs` but NOT TweakReg.** It's the only instrument that can't skip
`updatewcs`: AstroDrizzle needs the NPOL/D2IM distortion arrays, which the chip extraction
doesn't carry over (without them it stops on the missing DGEO correction). `--align mast`
runs `updatewcs(use_db=True)` (restoring the `GSC240` fit) then stops. `GSC240` is only
~0.5″ off in *absolute* astrometry but its frame-to-frame registration is ~0.02–0.03″,
which is what the stack needs; the absolute offset is fixed afterwards by
`align_wfpc2_to_acs.py`. A per-lens core-registration audit (LACosmic-masked) put **all 22
lenses on `mast`** — TweakReg scatters the mostly-single-visit frames ~0.7″, splitting the
core into ~4 knots. The choice is stored in `info/wfpc2_alignment.json`.

**Multi-visit lenses are split, not TweakReg'd.** J0728+3835 and J0822+2652 each have two
visits at a ~14–16° roll (two guide-star solutions). They're drizzled as separate per-visit
datasets (`--pa <PA_V3> --out-suffix _v1/_v2`, each single-guide-star → `mast`) into
`f606W_v1`/`f606W_v2`. Outcome: **J0822+2652** = `f606W_v1` (2×1100s) + `f606W_v2`
(4×1100s); **J0728+3835** = `f606W_v2` only (its 2-frame visit has 1 x-dither-phase, can't
reach 0.05″, dropped); **J1142+1001 stays combined** (visits share roll, PA 119.00 vs
118.87). `align_wfpc2_to_acs.py --f606-dir` and `make_cutouts.py --filt f606W_v1` handle
the suffixed dirs.

**Diagnosing this class of bug:** compare the *spread* of per-frame WCS error, not its
magnitude — a common offset is a harmless absolute-astrometry shift; frame-to-frame scatter
is what smears a stack. Cleanest metric: per-frame core-registration scatter on the
drizzled common grid, CRs masked; confirm against the visible product (one core, not
knots). Stacked stellar FWHM can mislead — it centroids an extended galaxy and rewards
TweakReg's self-consistency even when the deflector is split.

TweakReg `threshold` is in **image data units** and does not transfer between detectors
(WFPC2/WF3 100, WFC3/IR 20, ACS default). This only takes effect under `--align tweakreg`,
which no lens uses — dead code kept for comparison runs. → memory:
tweakreg_threshold_per_instrument

### `align_wfpc2_to_acs.py` — F606W absolute astrometry (after drizzle, before cutouts)

WFPC2 F606W carries only GSC 2.4.0, ~0.3–1″ off absolute; ACS F814W and WFC3/IR F160W carry
GAIA eDR3 / GSC242 (<0.02″, agree to ~0.01″). So F606W sits ~0.5–0.9″ off the other bands —
a whole-mosaic shift (harmless to the F606W stack, breaks cross-band registration). The
script ties the deflector light-centroid to ACS F814W via an iterative windowed
`centroid_com` (robust to the ring), shifts F606W `CRVAL1/2` so the centroids coincide, and
stamps `GSC240FX=True`. Idempotent; refuses any tie implying > `MAX_SHIFT = 1.5″`; uses
F160W as an *independent* check where present. Verified J0252+0039: 0.66″ → 0.009″.

```bash
conda run -n stenv python scripts/align_wfpc2_to_acs.py --lens J0252+0039   # or --all
```

**Run order: ACS + WFPC2 drizzles → `align_wfpc2_to_acs.py` → `make_cutouts.py`.**

## Weight maps and noise: `final_wht_type` and `cutout_noise.fits`

All scripts take `--wht-type {ERR,IVM,EXP}`. `cutout_noise.fits` is `1/sqrt(WHT)`, so the
weight type *is* the noise model. Default **`ERR`** for ACS/WFC3/NICMOS, **`IVM`** for
WFPC2 (override below).

- **`ERR`** — full inverse-variance (source Poisson + sky + read + dark): the correct
  per-pixel σ for modelling. The blank-sky floor is *included* — don't "add it back".
- **`EXP`** — uncalibrated exposure-time map, no source shot noise (core/sky ratio ~1.04).
  Not for a likelihood.
- **`IVM`** — inverse-variance map; DrizzlePac's auto one is background-only, but the WFPC2
  override supplies a real per-pixel IVM.

ERR captures the per-pixel variance but *not* the off-diagonal covariance from drizzle
resampling — see *Drizzle correlated noise*.

### WFPC2 noise: two fixes, both required (`IVM`, and a real noise model)

→ memory: wfpc2_err_weighting_not_supported

1. **DrizzlePac ignores ERR for WFPC2.** `WFPC2InputImage` hardcodes `errExt = None`, so
   `--wht-type ERR` silently falls back to exposure-time weighting (measured core/sky ratio
   exactly 1.000, vs ~3.5 for ACS). Fix: `build_ivm_files()` builds a per-frame IVM and
   feeds it via a two-column `@`-association file, which DrizzlePac *does* honor for WFPC2.
2. **The WFPC2 ERR array is bogus:** `ERR == sqrt(SCI)` as an *identity* on 100% of pixels
   in all 92 frames — Poisson on DN as if they were electrons, no gain, no read-noise term
   (overstates true noise 2.11× at sky). Do **not** build IVM from it. The model used is
   `var_DN = SCI/gain + floor²` (Poisson slope from header `ATODGAIN=7.0`; additive floor
   measured per frame from its own blank sky via `measure_noise_floor()`). `K=1` is correct
   for WFPC2. Products carry `IVMMODEL='SCI/gain+floor^2'`; `make_cutouts.py` warns when the
   keyword is absent (three indistinguishable generations of WFPC2 weight map exist). A
   pre-fix product **cannot be rescued by a scale factor** — it needs a re-drizzle.

Two transferable traps: never select blank sky by pixel *value* (biases the width low —
select by position or difference adjacent pixels); a model that fits at sky level can still
be ~20% wrong at the core, so widen the lever across frames before trusting it.

### `1/sqrt(WHT)` is a σ map only if the input ERR is in counts

→ memory: noise_map_err_units. DrizzlePac computes `weight = (EXPTIME/ERR)²`:
- **ACS FLC** ERR in ELECTRONS → `EXPTIME/ERR = 1/σ_rate` → K=1, correct.
- **WFC3/IR FLT** ERR in ELECTRONS/S → weights inflated by EXPTIME², noise map came out
  **EXPTIME (599.2 s) too small** (SNR ~60,000).

`weight_to_sigma_scale()` applies `K = per-frame EXPTIME` for detectors in
`_ERR_IN_RATE_UNITS = {('WFC3','IR')}`, K=1 otherwise, and refuses to guess if `D00nDEXP`
per-frame times are unequal. **Diagnose with the blank-sky block-sum test, never per-pixel**
(drizzle correlation drives per-pixel MAD to 0.73 ACS / 0.36 F160W even when correct).
Constants stamped as `NOISEK`/`NOISECOR`.

### F606W `BUNIT`, and cross-band units

AstroDrizzle writes count *rates* (records `D001OUUN='cps'`) but doesn't rewrite `BUNIT`.
For WFPC2 that left `BUNIT='COUNTS'` on products holding DN/s — an EXPTIME-sized error for
anything reading `BUNIT`. The script now stamps `BUNIT='COUNTS/S'` on WFPC2 SCI products
(the WHT map stays `UNITLESS`), guarded on `D001OUUN=='cps'`. The two optical bands remain
in **different unit systems** — F606W in DN/s (instrumental), F814W in e/s — self-consistent
per band via `PHOTFLAM`; any cross-band flux comparison must go through `PHOTFLAM`, not raw
pixel values.

## Common output WCS across filters

Every band is pinned to the same output grid geometry so filters co-register pixel-for-pixel
with no later reprojection: `final_rot=0.0` (North-up) and `final_ra`/`final_dec` at the
lens catalogue position (from `info/slacs_coords.py` for SLACS, `info/gallery_coords.py` for
`gallery`) go into every AstroDrizzle call. A lens missing from that table falls back to the
native drizzle WCS with a warning. This aligns orientation + tangent point, not absolute
astrometry (still the delivered WCS's job, and `align_wfpc2_to_acs.py`'s for F606W — gallery
has no equivalent cross-band tie script, and doesn't need one; see *BELLS GALLERY* below).

## Cosmic-ray rejection: LACosmic, not `driz_cr` (ACS **and** WFPC2)

→ memory: acs_cr_pass_eats_core. `--cr` masks CRs per frame with LACosmic (default
`--cr-method lacosmic`) then drizzles a plain weighted mean (`median=False, blot=False,
driz_cr=False`); `--cr-method drizcr` restores the old route. WFPC2 uses the same LACosmic
route (`--lacosmic-sigclip 4.5 --lacosmic-objlim 5.0`; gain/readnoise/saturate from the WF3
header). Mask written to DQ bit 4096 with `resetbits=0`.

`driz_cr` compares each frame to a blotted median; on a steep PSF core that reference is
systematically low, so real core pixels read as CRs — it destroyed ~37% of the deflector
core flux, and loosening thresholds made it *worse* (the reference is the fault, not the
cut). LACosmic preserves the core (peak 1.000, 1″ 0.979).

**`resetbits=0` is mandatory on the LACosmic pass** — it defaults to 4096 and would clear
the DRIZ_CR bit the mask lives in, silently producing an un-masked drizzle that still looks
plausible (the first run reported a flawless `core=1.000` that was the no-CR image scored
against itself). Conversely the no-CR passes pin `resetbits=4096` explicitly, else they
inherit the mask and become silently CR-rejected.

**Do not add CR rejection to F160W.** At 0.1283″/px the IR PSF is ~1 px FWHM and looks like
an outlier: LACosmic zeroes field stars, `driz_cr` costs ~10% of a star's peak, at every
`objlim`. WFC3/IR FLTs are already up-the-ramp CR-rejected. → memory: wfc3ir_quadrupled_defects

## Output pixel scales

Chosen by measurement (weight-map uniformity + FWHM + noise-map correlation), not
convention. → memory: drizzle_correlated_noise

| Band | Instrument | Native | Output | pixfrac |
|---|---|---|---|---|
| F606W | WFPC2/WF3 | 0.0996″ | 0.05″ | 1.0 |
| F814W, F555W | ACS/WFC | 0.05″ | 0.05″ (native) | default |
| F160W | WFC3/IR | 0.1283″ | 0.06″ | 1.0 |

- **F160W stays at 0.06″** — do not re-drizzle to 0.05″ to grid-match. PyAutoLens ingests
  each band at its native scale and fits sub-pixel inter-band offsets, so a common drizzle
  grid buys nothing; 0.05″ opens no empty pixels but worsens weight non-uniformity. **Kept
  at 0.06″ / pixfrac 1.0, settled by the user 2026-07-26 after a scale/pixfrac scan — do
  not reopen without being asked.** → memory: wfc3ir_quadrupled_defects
- **F160W and F606W use pixfrac 1.0, not 0.8.** With ERR-based noise the goal is a uniform,
  low-correlation noise map for the likelihood, not the sharpest PSF (which PyAutoLens fits
  explicitly). pixfrac 1.0 drops adjacent-pixel noise correlation markedly at a small PSF
  cost. ACS F814W/F555W keep the default pixfrac (native 0.05″, no oversampling).
- **F555W needs no new script** — `drizzle_acs_wfc.py` is filter-agnostic. The 16 F555W
  lenses are exactly the 16 without WFPC2 F606W (props 10494/10798).

## Drizzle correlated noise (matters for strong-lens modelling)

→ memory: drizzle_correlated_noise. Drizzling onto a finer-than-native grid correlates
adjacent output pixels (each output pixel is a weighted sum of overlapping input pixels;
neighbours share input pixels, so their noise is covariant) — intrinsic to drizzle
resampling, stronger with oversampling. `final_wht_type=ERR` makes the **per-pixel
variance** correct but says nothing about the **off-diagonal covariance**: on J1143 blank
sky the empirical pixel RMS was ~1.47× the ERR-map prediction (pixfrac 0.8; pixfrac 1.0
shrinks but doesn't remove it).

**Modelling implication:** a per-pixel independent-Gaussian likelihood (diagonal covariance)
does not bias the best-fit but mis-estimates *uncertainties/evidence*, worse in the
oversampled F606W/F160W than in native-scale F814W. The residual correlation factor is
band-dependent (~1.24 ACS, ~1.17 F160W, ~1.1 F606W, ~1.5–1.6 WFC3/UVIS gallery — higher
than native ACS despite also being native-scale, plausibly UVIS's larger geometric
distortion; see *BELLS GALLERY* below), so ignoring it also mis-weights bands relative to
each other in a joint fit. `make_cutouts.py --corr-factor` applies it (default 1.0, off).
Prefer native-scale F814W where a clean per-pixel noise model matters most.

## Cutouts (`scripts/make_cutouts.py`)

```bash
conda run -n stenv python scripts/make_cutouts.py --lens J0029-0055 --filt f606W
```

Cuts a square stamp (default 20″) from `data/drizzled/` into `data/cutouts/`: a sci FITS, a
noise FITS (from the weight map), and a 3-panel PNG.

- **Pass + prefix.** `--pass {auto,cr,nocrrej}` (default `auto`) picks the CR pass when one
  exists, else no-CR. The prefix encodes it so the two coexist: **`cutout_cr_*` for CR,
  `cutout_*` for no-CR.** ACS F814W/F555W and WFPC2 F606W now default to a LACosmic CR pass
  → `auto` cuts `cutout_cr_*`, which is science-grade (LACosmic preserves the core). WFC3/IR
  F160W has no CR pass, so `auto` falls back to `cutout_*` there. → memory:
  cutout_centring_on_cr_pass
- **Recentring.** The stamp recentres on the galaxy, and the peak search prefers the CR
  mosaic — a brightest-pixel search on a CR-laden no-CR mosaic locks onto cosmic rays. Don't
  reach for `--median-size` for a bad recentre; check a `*_cr_*` mosaic is present.
  Offsets of 1–2″ aren't always failures — some lenses (J0912+0029, J0956+5100) have genuine
  multi-knot morphology.
- **Shared centre (`--center-band`, default `f814W`).** All bands are cut about a single
  centre from the `--center-band` mosaic (highest S/N, GAIA-accurate) so stamps co-register
  across filters. `--center-self` restores per-band recentring; a band whose center-band
  products are missing falls back to its own peak with a warning. **Never use
  `--center-self` on gallery's UV bands** (F225W, F275W) — the early-type deflector is
  UV-dark and a self-centred peak search locks onto the bright lensed arc instead; see
  *BELLS GALLERY* below.

## PSF generation (`scripts/make_psf.py`, `scripts/psf_models.py`)

```bash
conda run -n stenv python scripts/make_psf.py --lens J0252+0039 --filt f814W
bash scripts/run_psf_all.sh            # every drizzled product; globs data/drizzled/ like run_cutouts_all.sh
bash scripts/run_psf_all.sh --models-only   # skip existing empirical builds; rebuild only the model tier
```

**One filter-agnostic script, per-instrument defaults** (keyed on `INSTRUME/DETECTOR`; WFPC2
is keyed on `INSTRUME` alone because the drizzled primary header says `DETECTOR=PC`, the
aperture — see *WFPC2: the lens is on WF3*). Builds an **empirical ePSF** from field stars in
the full drizzled mosaic (photutils `DAOStarFinder → extract_stars → EPSFBuilder`, a
production port of `old_notebooks/load_data.ipynb`), and **falls back to a model** when a
field is too star-poor — for WFPC2 F606W a **native ePSF from the MAST PSF database**, else
an STScI STDPSF model (`psf_models.py`). **Two products, two homes:**
- **`data/psf/<sample>/<lens>/<filt>/`** (archival characterisation): **`psf_kernel.fits`**
  — the **full** kernel (whole ePSF footprint binned to image scale, `star_size` px:
  35/41/51), **`psf_epsf.fits`** (oversampled ePSF), **`psf.png`** (4-panel QA: star montage
  · kernel linear · kernel **log** · radial profile).
- **`data/cutouts/<sample>/<lens>/<filt>/`** (modelling-ready): **`cutout_[cr_]psf.fits`** —
  the **trimmed** kernel, cut to the amplitude-`--trim-threshold` (default 1e-3 of peak)
  radius and pass-matched to `cutout_[cr_]sci.fits`. Written by `make_psf.py` itself.

Records `info/lens_psf.json` (`{sample:{lens:{filt:{method,n_stars,fwhm_pix,oversample,
kernel_size,cutout_kernel_size,trim_threshold}}}}`, `null` + exit 0 on no data). Run
**after** the drizzles (a PSF needs a mosaic). → memory: psf_generation, psf_kernel_sizing

**Kernel size is an amplitude cut, not enclosed-energy.** The trimmed modelling kernel is cut
where the azimuthally-averaged PSF drops below `trim_threshold`×peak — the extent over which
the PSF still spreads flux above ~that level, which is what a *convolution* kernel needs. A
95% enclosed-energy cut is **wrong** here: EE is area-weighted (a photometry criterion) and
integrates the noisy empirical wing, so it under-sizes — on J0252 F606W it truncates at ~18px
while the PSF is still ~1% of peak, vs ~31px for amplitude-1e-3. The cut is band-adaptive by
construction (sharp ACS F814W → 19–27px; broad F160W → 29–41px; F606W → 31px).
→ memory: psf_kernel_sizing

**Star selection is fully automatic; per-lens tweaks live in `info/psf_stars.json`**
(`{sample:{lens:{filt:{...}}}}`; absent sample/lens/filt ⇒ automatic). Overrides:
`include`/`exclude` coords or boxes, and any parameter
(`max_stars`, `threshold_scale`, `min_snr`, `oversample`, …). The same knobs exist as CLI
flags. Precedence: instrument default < JSON < CLI. This replaces the notebook's hand-typed
NaN rectangles and manual star deletion.

**The traps (all cost a silently-wrong PSF):**
- **A 5σ DAO detection is not a PSF star.** On star-poor fields the only round detections are
  ~5σ noise blobs (measured peak-S/N 4–6 on J0252 F606W) that build a *pure-noise* ePSF which
  passes a naive "peak is centred" check. Gate on an **absolute peak-S/N floor** (`min_snr=30`)
  *and* validate the ePSF core stands ≥15× above its own outskirt noise; else fall back to the
  model. Star/galaxy separation is a per-candidate Gaussian fit on a ~6·FWHM window (galaxies'
  extended light inflates the fitted FWHM past the cut; an 11px window only sees the nucleus).
- **EPSFBuilder diverges with few stars on a large stamp** (peak drifts to a corner). WFPC2
  uses `star_size=35` (converges at 3 stars where 51 diverged); ACS/WFC3 keep 51/41.
- **`oversample=2`, not 4.** SLACS fields yield only ~5–20 stars — a 4×-oversampled ePSF has
  too many pixels per star and comes out noisy. Raise per lens only where a field is star-rich.
- **Flux floor.** EPSFBuilder normalises each star by its flux, so a faint star amplifies its
  background noise into the ePSF wings; drop stars fainter than 5% of the brightest kept.

**ACS/WFC model = focus-diverse ePSF, not STDPSF.** When an ACS/WFC field is too star-poor
for an empirical build (or `--method model`), `psf_models.acs_focus_diverse_psf` retrieves
the **observation-matched, focus-corrected ePSF** (`acstools.focus_diverse_epsfs`; Bellini
et al. ACS ISR 2018-08 / 2023-06) for each contributing exposure (rootnames from
`lens_products.json`), interpolates each to the lens's detector position via the FLC WCS
(chip from `CCDCHIP`: 1→WFC1, 2→WFC2; grid centre if the FLC is absent), averages across
exposures (the drizzled PSF is the exposure-average), and resamples the 4×-supersampled
detector grid to the output scale. This is **native F555W and F814W** — no filter
substitution — and matched to the HST focus/breathing of the actual exposures, so it's a
strictly better ACS model than STDPSF. Records method `model_acs_fdpsf`; grids cached under
`data/reference_files/acs_fdpsf/`. Detector-frame (omits drizzle broadening) but **rotated
per-exposure to North-up before averaging** (see *Model PSFs are rotated to North-up*), so the
empirical ePSF is preferred when stars exist; **STDPSF stays the fallback-of-the-fallback**
if a retrieval fails. Verified J0252 F814W: 0.100″ (STDPSF F814W 0.096″, empirical 0.125″).

**WFPC2 F606W model = native ePSF from the MAST PSF database (`psf_models.py`,
`scripts/mast_api_psf.py`; Dauphin et al., ISR WFC3 2021-12).** STDPSF has no WFPC2 F606W
grid, so the STDPSF path substitutes WFPC2 F555W (right chip, wrong filter) — historically
the least-verified product in the pipeline. Instead `wfpc2_f606w_db_epsf()` queries the MAST
PSF database for good-quality (`qfit<0.05`), unsaturated (`n_sat_pixels=0`) **WF3 (chip 3)**
F606W star cutouts near the lens position (`x_cal/y_cal` within ±200px of ~435,424 — every
F606W lens puts its target at the same WF3 spot), downloads the `c0m[3]` cutouts, gates them
(inner-window crop + centroid to drop edge neighbours/warm pixels in the un-CR-cleaned c0m),
and builds **one shared native-F606W WF3 ePSF** (~147 stars) cached under
`data/reference_files/wfpc2_f606w_psfdb/`. Every F606W product (incl. split-visit `_v1/_v2`)
reuses it; records method `model_wfpc2_psfdb`. The cached ePSF is detector-frame; each lens
**rotates it to North-up at resample time via that lens's WF3 exposure CD** (so one shared
build serves every roll — see *Model PSFs are rotated to North-up*), and the pedestal
subtraction removes its ~1e-3 wing floor. Detector-frame still omits drizzle broadening, so
the lens's-own-field empirical build is preferred where stars exist — this is the model
**fallback, slotted ABOVE the STDPSF F555W proxy** (now fallback-of-the-fallback). **Traps:**
the MAST PSF DB `chip` column is the WFPC2 CCD/FITS ext (1=PC…4=WF4) — **query chip 3 for
WF3**, not chip 1 (PC has a different pixel scale *and* PSF); split-visit filter keys
(`f606W_v1`) must be normalised to the base filter before any STDPSF/pivot lookup
(`_base_filter`), or the model path KeyErrors. FWHM ~0.22″ (post-rotation). → memory:
wfpc2_f606w_mast_psf_db

**STDPSF fallback (`psf_models.py`, Anderson & King 2000; Dauphin et al. 2021; Anderson
2016):** 4×-supersampled 101×101 grids read by photutils
`GriddedPSFModel.read(..., format='stdpsf')`, cached under `data/reference_files/stdpsf/`.
Used for WFPC2/WFC3, and for ACS only when the focus-diverse retrieval fails.
- **Neither our WFPC2 F606W nor our ACS F555W has an exact STDPSF grid** — the library skips
  them. `_resolve_filter` substitutes the nearest published band by pivot wavelength (WFPC2
  F606W→**F555W**, ACS F555W→**F606W**) with a printed NOTE. WFC3/IR F160W and ACS F814W are
  exact. (The ACS F555W substitution only bites if the focus-diverse path above also fails.)
- **WFPC2 grids are per-chip** (3×3 fiducials × 4 chips) — select WF3 with `detector_id=3`;
  ACS/WFC has two chips (`detector_id=1`); WFC3/IR is single-detector.
- The STDPSF is defined on the **detector native scale**, so it is resampled to the drizzled
  output scale (WFC3/IR 0.1283″→0.06″); skipping this makes the model the wrong size.
- It is the **detector-frame** ePSF and omits AstroDrizzle broadening (Anderson 2016), so it
  runs slightly sharp — the **empirical ePSF is the true drizzled PSF and always preferred**
  when enough stars exist. It is **rotated to North-up** like the other model tiers (see
  below). WFPC2 fields are usually star-poor (A&K build ePSFs from globular clusters, which
  these are not), so F606W falls to the model — now the native MAST-PSF-DB ePSF above, with
  this F555W proxy only as its fallback.

### Model PSFs are rotated to North-up; empirical ones already are

→ memory: psf_model_northup_rotation. Every **model** ePSF (STDPSF, ACS focus-diverse, WFPC2
MAST-DB) is a *detector-frame* build: its axes are the exposure's detector axes, so its
diffraction spikes / asymmetric wings sit at the exposure roll (`ORIENTAT`, up to ~105° for
SLACS) relative to the North-up drizzled science image (`final_rot=0.0`). `psf_models`
resamples each into the output frame through the exposure **CD matrix** (`_northup_M` —
rotation + parity + scale in one map, not an Euler angle), the CD read from a contributing
exposure by `make_psf.representative_input_cd`: WFPC2 from the extracted WF3 file in
`data/drizzle_files` (single WF3 chip, correct **per-visit** roll for the split lenses),
WFC3/IR from the calibrated FLT SCI, ACS per-exposure from each FLC (rotated **before** the
focus-diverse exposure-average, so multi-roll is handled). The **empirical ePSF needs no
rotation** — it is cut from the North-up mosaic itself, which is exactly *why* its orientation
can be trusted. Verified against the WCS chain on all three detectors (ACS/WFPC2/WFC3-IR).
Without a CD the model is left unrotated with a printed NOTE (a degraded but not wrong-scale
fallback).

### Drizzle-broadened model PSF by injection (`make_psf_inject.py`, the Anderson 2016 route)

Every model tier above is a **detector-frame** ePSF that `psf_models` resamples/rotates to
North-up *analytically* — reproducing orientation and scale but **not** the extra blur
AstroDrizzle's resampling puts on a point source, so the analytic model kernel runs sharp
(Anderson 2016, WFC3/IR ISR). `make_psf_inject.py` implements the rigorous fix that ISR
names: inject the model PSF as artificial stars into the individual exposures and **re-drizzle
them exactly as the science frames were**, so the drizzled star carries the broadening,
North-up orientation, and exposure-average weighting for free — all produced by the real
drizzle, not emulated. **Model tier only** (empirical builds are already the drizzled PSF, cut
from the mosaic). It is a **parallel comparison product** — nothing downstream changes:

- **Reuses the persisted `data/drizzle_files/<sample>/<lens>/<filt>/` inputs** the science
  drizzle consumed (ACS `<root>_flc`, WFC3/IR `<root>_flt`, WFPC2 extracted `wf3_<root>_flt`
  + per-frame IVM + the two-column `@`-association), inheriting every prep step (WF3
  extraction+renumbering, distortion, updatewcs, IVM weighting). **Requires those inputs
  present** — re-run the band's drizzle first if `drizzle_files` was cleared.
- Per frame: **zeroes SCI, adds the detector-frame model PSF** at the lens position (that
  frame's WCS), keeping ERR/DQ/IVM untouched (weighting identical to science; the injected
  star is clean so no CR pass). Model source mirrors make_psf: ACS focus-diverse per exposure
  (STDPSF fallback), WFC3/IR exact-filter STDPSF, WFPC2 the shared MAST-DB ePSF (STDPSF F555W
  fallback) — recorded as `inject_acs_fdpsf` / `inject_stdpsf` / `inject_wfpc2_psfdb`.
- Re-drizzles onto the **same output grid** (`final_rot=0`, `final_ra/dec`=lens, per-band
  `final_scale`/`pixfrac`/`bits`/`wht_type` lifted from the drizzle scripts). The drizzled
  star is the broadened North-up PSF at the modelling scale; the kernel is cut/centred at
  image scale, then `subtract_pedestal`/`trim_kernel_to_amplitude` as usual.

Outputs alongside make_psf, non-clobbering: `data/psf/<...>/psf_kernel_injected.fits` +
`psf_injected.png`; `data/cutouts/<...>/cutout_[cr_]psf_injected.fits` (trimmed); records
`info/lens_psf_injected.json` (`null`+exit 0 on no data). Runner:
`bash scripts/run_psf_inject_all.sh` (model tier only; `--all` also runs empirical products
for validation). **Validated** on F160W where an empirical truth exists on the same lens:
analytic-STDPSF 3.15px → **injected-STDPSF 3.5–3.9px ≈ empirical 3.8px** (drizzle broadening
recovered; the residual is the real optical wing STDPSF underestimates, not a broadening
error). ACS broadens little (native 0.05″, small resampling — FD 1.65→1.76px); F160W and
WFPC2 broaden more (0.1283→0.06″, 0.0996→0.05″). → memory: psf_injection_drizzle_broadening

### Pedestal subtraction (`subtract_pedestal`) — every kernel

→ memory: psf_kernel_sizing. EPSFBuilder leaves a small flat DC floor in the ePSF wings. Left
in, a ~1e-3-of-peak pedestal across the kernel is several % of the (renormalised) flux as a
spurious *uniform background*, **and** it stops the amplitude trim from ever crossing 1e-3 —
the trimmed kernel then caps at the full `star_size` (this was the real cause of the old
"F160W full kernel too tight", not a genuinely larger IR PSF; the clean STDPSF F160W trims to
27–29px). The outer-annulus median is subtracted from every full kernel before the archival
write **and** the trim, recorded as `PSFPED` / `pedestal_frac`. Near no-op on the sharp ACS
bands (~1e-4); real on F606W (~1e-3) and F160W empirical (~1e-3). Applied to all methods.

### F160W hybrid quality gate

A *validated* empirical ePSF whose wing pedestal or scatter exceeds
`pedestal_bad`/`scatter_bad` (both 3e-3, in `_BASE`) is dropped to the model under `--method
auto` — star-poor oversampled F160W fields build noisy wings that still pass the core checks.
Thresholds chosen to pass every clean ACS/F555W build (worst scatter ~2e-3) and the good F160W
empirical builds (≤1.2e-3) while dropping the noisy ones (J0936 7.8e-3, J0946 6.6e-3). Keeps
the drizzle-broadened empirical PSF (FWHM ~3.8px, which the detector-frame model lacks at
3.15px) where the build is clean, and the reproducible model where it isn't.

### `run_psf_all.sh --models-only`

`run_psf_all.sh [SAMPLE] [--all|--models-only]`. Default `--all` rebuilds everything. Use
`--models-only` after changing only model-PSF code (e.g. the rotation): it skips any product
whose *existing* `lens_psf.json` method is `empirical` (that code never touches the empirical
builds — cut from the North-up mosaic; the pedestal is a ~1e-4 no-op on ACS) and rebuilds the
model tier plus any not-yet-built product. The skip set is "empirical AND already recorded",
so a new lens still builds.

**Current state and open limitations.** `run_psf_all.sh` **has been run across all of
`slacs_gold`** — `info/lens_psf.json` holds **90 products, 0 failures** (2026-07-28), each
carrying `pedestal_frac`; the model tier is rotated to North-up. Method breakdown: **F814W**
33 empirical + 5 focus-diverse; **F555W** 16 empirical; **F606W** 23/23 `model_wfpc2_psfdb`
(native MAST-DB, no more F555W proxy); **F160W** 4 empirical + 9 exact-filter STDPSF `model`
(the hybrid gate dropped J0936/J0946 from the previous 6 empirical). No capped trims remain
(all F160W 27–29px). The gate thresholds (`min_snr=30`, core ≥15× outskirt, flux floor 5%,
`fwhm_tol_hi=1.4`, `star_size=35` for WFPC2, `pedestal_bad`/`scatter_bad`=3e-3) generalised
fine — no retune needed. **Not yet run for `slacs_other`** (reduced 2026-07-29, same three
instruments so `run_psf_all.sh slacs_other` needs no code change) **or `gallery`** (needs
WFC3/UVIS support in `make_psf.py`/`psf_models.py` first — see *BELLS GALLERY* below). Open
items:
- **No PSF uncertainty is propagated.** The kernel is a point estimate; an ePSF built from
  as few as 3–5 stars has no accompanying error/covariance for downstream lens modelling.
  Genuinely new work — not addressed by having more lenses.
- **Model-PSF rotation uses cubic resampling** (`_resample_centered`, order 3, was order
  1/bilinear — bilinear softened the model kernel slightly, e.g. F606W DB FWHM 0.212″→0.221″).
  Code and the regenerated model-tier products (`model`, `model_acs_fdpsf`,
  `model_wfpc2_psfdb`) are **committed together** (commit 92b1003, 2026-07-28); `info/lens_psf.json`
  and the on-disk `cutout_[cr_]psf.fits` kernels reflect the cubic resample. (`data/psf/` archival
  kernels are gitignored, not tracked.) No longer an open item.
- F160W's 9 STDPSF models are **exact-filter** (F160W has a real grid), so acceptable. Using the
  MAST PSF DB to build a *native* WFC3/IR F160W ePSF for injection instead of STDPSF was
  **prototyped and measured (2026-07-28): no gain** — the native-DB kernel is indistinguishable
  from exact-filter STDPSF (FWHM 3.97 vs 3.93px on J0728+3835, wings agree to ≤1e-3) and neither
  closes the ~2–3% wing deficit vs the empirical drizzled truth, because a DB build averages
  breathing out just like STDPSF. Unlike the WFPC2 F606W DB win (which replaced a *wrong-filter*
  proxy), F160W already has the right-filter grid, so there's nothing to fix. Only focus-matched
  retrieval (à la ACS focus-diverse) could close the wing gap — not the DB. → memory:
  f160w_mast_db_injection_no_gain

## Tracking JSONs in `info/`

> **Full reset, 2026-07-26.** All three files were emptied to `{}` and every product under
> `data/` deleted (`calibrated/`, `drizzle_files/`, `drizzled/`, `cutouts/`, `mosaics/`,
> `run_logs/`) as a deliberate clean restart. Kept: `data/reference_files/` (CRDS cache) and
> `data/pre_drizzled/` (46 MAST-delivered mosaics, not pipeline output). The sample was
> renamed `slacs` → **`slacs_gold`** in the same pass. **Any surviving `data/*/slacs/` path,
> and any "current on-disk state" / "Not regenerated" claim in memory, is pre-reset and
> void** — the *reasoning* in those notes stands and is why reruns use the current scripts;
> only the inventory is stale.

> **Split by sample, 2026-07-29.** Every tracking JSON below (plus `lens_psf.json`,
> `lens_psf_injected.json`, `wfpc2_alignment.json`, `psf_stars.json`) was flat `{lens:
> {...}}`, mixing slacs_gold/slacs_other/gallery lenses in one namespace — harmless only
> because no lens name has ever collided across samples. All are now nested `{sample:
> {lens: {...}}}`, matching `lens_samples.json`'s own top-level-by-sample layout and the
> `data/<sample>/<lens>/...` directory convention. Every read/write site now goes through
> `scripts/info_json.py` (`load`/`update`), which consolidated what had been 4 near-
> identical copies of the same read-modify-write helper across the drizzle scripts. The
> migration was a pure reshape (verified by round-tripping every value back to its
> pre-migration flat form); no data changed.

Updated automatically by every run:
- **`lens_products.json`** — `{sample: {lens: {key: [rootname, ...]}}}` — frames that
  reached the drizzle (not the whole download).
- **`lens_instrument.json`** — `{sample: {lens: {key: "INSTRUME/DETECTOR"}}}` — records
  `WFPC2/WF3` for F606W (the chip), not MAST's `WFPC2/PC`.
- **`lens_exptime.json`** — `{sample: {lens: {key: seconds}}}` — from the CR-rejected
  drizzle header.

No data for a filter → value `null`.

- **The key is the product directory, not the filter.** Usually they coincide (`f606W`), but
  a split-visit lens is keyed per visit (`f606W_v1`/`f606W_v2`) with **no bare `f606W` key**.
  The old "missing key = out of sync" check no longer holds for J0728+3835/J0822+2652 —
  check keys against product directories instead. This fixed a live error: both split lenses
  had recorded a plausible-but-nonexistent combined `f606W` (6 obsids, 6600 s), invisible
  precisely because the key was present. Root cause: JSON writes keyed on the bare filter
  while writing to a `--out-suffix` dir; now keyed on `product_key = filt + out_suffix`.
- **Records drizzled frames, not the download.** WFPC2 `--pa` selects one visit; a lens that
  exits for want of dither phase writes nothing. ACS/WFC3 silently drop `EXPTIME=0` frames,
  so those are excluded from the record but not deleted (unlike WFPC2, where `MIN_EXPTIME`
  removes them outright). Exposure times were always correct (from the header).
- Auditing obsid counts vs `NDRIZIM`: **ACS/WFC FLCs are 2-chip MEFs**, so `NDRIZIM =
  2×exposures` there and 1× for WFPC2/WFC3 — comparing directly reports 54 false mismatches.
- All three levels stay sorted (sample, lens, and filter/key within each lens) across
  partial runs.

## Lens Samples

**`info/lens_samples.json` is the single source of truth** for sample membership and per-lens
MAST quirks (`mast_target`, `force_copy`). Read it only through
`scripts/mast_target_names.py` — never parse it elsewhere, never keep a second lens list
(`info/list_of_lenses.txt` was exactly that and was deleted 2026-07-26).

```bash
conda run -n stenv python scripts/mast_target_names.py --list        # samples + sizes
conda run -n stenv python scripts/mast_target_names.py slacs_gold    # lens names
```

| Sample | Lenses | What it is |
|---|---|---|
| **`slacs_gold`** | 38 | The working sample; **default `--sample` of every script** |
| **`slacs_other`** | 27 | Rest of SLACS restricted to Bolton et al. 2008 Table 4 class E-S-A/L-S-A (incl. the `*` variants). Reduced |
| **`gallery`** | 15 | BELLS GALLERY (props 14189, 16734), restricted to Shu et al. 2016 Table 1 class E-S-A; WFC3/UVIS multi-band. Reduced |

`slacs_gold` coverage: **F814W** (ACS/WFC, all 38), **F606W** (WFPC2/WF3, 22 — the other 16
have no WFPC2 data), **F555W** (ACS/WFC, 16 — exactly those 16), **F160W** (WFC3/IR, 13 —
prop 11202). NICMOS F160W (24 lenses) is deprioritised and its data deleted. HST props:
10886, 11202, 10494, 10798.

`slacs_other` coverage (27 lenses, reduced 2026-07-29): **F606W** (WFPC2/WF3, 24/27),
**F814W** (ACS/WFC, 3/27 — the rest are `BLOCK_EXPTIME`-gated or absent, see *Total-
exposure-time gate*), **F160W** (WFC3/IR, 6/27), **F555W** (ACS/WFC, 0/27 — none of these
lenses fall in props 10494/10798). Only 3 lenses (J0959+4416, J1153+4612, J1416+5136) have
both F606W and F814W, so `align_wfpc2_to_acs.py` has only tied those three to absolute
astrometry — the other 21 F606W products carry their delivered GSC240 absolute WCS
(~0.3–1″ off) untied. `info/wfpc2_alignment.json` has no per-lens `--align` audit for
`slacs_other` (it only covers the 22 `slacs_gold` WFPC2 lenses); every `slacs_other` WFPC2
lens falls back to the documented default, `mast`. No `lens_psf.json` entries yet — the PSF
pipeline (`run_psf_all.sh`) has not been run for this sample.

`gallery` coverage (15 lenses, reduced 2026-07-29): **F606W** on all 15 (the primary band),
**F814W**/**F438W** on 6, **F275W** on 5, **F225W** on 1 (J2342-0120). See *BELLS GALLERY:
WFC3/UVIS reduction* below for the pipeline and its caveats. Gallery lenses are **not in
`info/slacs_coords.py`** — they use `info/gallery_coords.py` instead, read by
`drizzle_wfc3_uvis.py` for the common output WCS; a lens absent from that table falls back
to native drizzle WCS with a warning. No PSF products: `make_psf.py`/`psf_models.py` have no
WFC3/UVIS support (only ACS/WFC, WFC3/IR, WFPC2/WF3 are keyed).

Caveat carried over from before both samples were reduced: `slacs_other`'s naming is
settled — all lenses resolve under plain `SDSS{lens}%`, no `GAL-*` overrides needed
(verified; see *Non-standard MAST target names*), so a "no observations" result there is a
real absence, not a naming gap.

**Every lens is tried every run; only ones with data download.** No-data is an ordinary
outcome (`null`, `=== NO DATA:` line, exit 0), counted separately from failures by the
runners — which is why they iterate the roster rather than globbing `data/calibrated/` (a
glob does nothing after a wipe and never picks up a new lens). A genuine download error
exits non-zero, kept apart by `mast_target_names.NoMastData` (a dedicated exception the
download block's broad `except` can't swallow).

MAST target names follow `SDSS<LENS>`; the query uses `target_name=f'SDSS{lens}%'` (wildcard
for naming variations).

**COPY handling differs by instrument, deliberately — don't unify without re-checking the
archive:**
- **ACS** filters COPY out in favour of non-COPY, **except** lenses with `"force_copy":
  true` (only `J1032+5322` F814W, whose non-COPY frames are `EXPTIME=0`). → memory:
  j1032_exptime_zero
- **WFPC2** keeps both (COPY sets are genuine repeat visits) and rejects junk on
  `MIN_EXPTIME`.

### Non-standard MAST target names (`GAL-*`)

Some lenses aren't on MAST under `SDSS<LENS>` — they use `GAL-<plate>-<mjd>-<fiber>`. Scripts
resolve this via `mast_target_names.py`: a lens with a `mast_target` in `lens_samples.json`
queries the `GAL-*` name first, falls back to `SDSS{lens}%`. **Output/directory names stay
in the J convention.** The live values are the `mast_target` entries in the JSON — edit
those, not any table. **Verified fact, not just current state:** all 14 `GAL-*` overrides are
in `slacs_gold`, and `slacs_other` (surveyed at its original 93-lens definition, before the
2026-07-29 restriction to the 27-lens Bolton E-S-A/L-S-A subset — a strict subset, so the
survey still covers it) has been confirmed to need none — no SLACS lens outside `slacs_gold`
should ever need a `mast_target` override for this reason. If one shows up in `slacs_other`,
treat it as a bug/regression to investigate, not a new legitimate case.

`gallery` (BELLS GALLERY) has a separate, also-verified naming mismatch: those targets sit on
MAST under their **full SDSS coordinate designation** (e.g. `SDSSJ002927.38+254401.7`), not a
short name at all — `GAL-*` doesn't apply here. `info/lens_samples.json`'s `gallery.lenses`
entries carry that full designation as `mast_target` for every one of the 15 lenses (queried
by `drizzle_wfc3_uvis.py` via `mast_target_names.py`, same mechanism as SLACS `GAL-*`), while
`info/`, `data/`, and every output path stay keyed on the short J-name (`J0029+2544`),
matching the SLACS convention.

## AstroDrizzle key parameters

- **CR pass** (default for ACS + WFPC2): **LACosmic** — mask CRs per frame, then a
  plain-mean drizzle (`median=False, blot=False, driz_cr=False`, `resetbits=0`). This is the
  product downstream reads. `--cr-method drizcr` restores the AstroDrizzle route (which eats
  the core).
- **No-CR pass** (`median=False, blot=False, driz_cr=False`): uncleaned, opt-in via
  `--nocrrej` for comparison. When both passes run they share one crop bbox (union of the two
  wht>0 boxes) so they register pixel-for-pixel.

ACS, WFPC2, and WFC3/UVIS default to **CR-only** (`--nocrrej` adds the no-CR pass; ACS/UVIS
`--no-cr` skips CR). WFC3/IR F160W has no CR pass at all — `make_cutouts.py` falls back to
the science pass for recentring (acceptable: FLTs are already up-the-ramp CR-rejected;
re-run with `--cr` if a recentre looks wrong).

**DQ bits treated as good** (do not unify these — they encode different detector facts):
- WFPC2: `8,1024`
- ACS/WFC: `256,64,16` (saturated, warm, stable hot)
- WFC3/UVIS: `256,64,16` (same meanings as ACS/WFC; STScI's own UVIS examples use `80`
  = 16+64 instead, dropping saturation from "good" — revisit if saturated cores prove a
  problem, not yet needed on gallery)
- WFC3/IR: `512` only — write it as `'512'`, **never `''`** (see the trap below)
- NICMOS: `2,4,8`

### WFC3/IR: quadrupled defects and the `bits=''` trap

→ memory: wfc3ir_quadrupled_defects, drizzle_bits_empty_string_trap

- **`final_bits=''` / `driz_sep_bits=''` disables DQ masking entirely.**
  `interpret_bit_flags('')` → `None`, and AstroDrizzle then keeps *every* flagged pixel as
  good (the opposite of intent) and silently voids any CR flag in DQ 4096. **`0` means
  "reject everything flagged"; `''`/`None` mean "reject nothing".** The script asserts this
  at import. Grep `bits :` in `astrodrizzle.log` to check.
- **Why defects quadruple.** The F160W no-CR pass does no cross-frame rejection, so a
  detector-fixed defect kept as "good" drizzles at each of the 4 dither positions → 4 sky
  replicas. Culprit is DQ 16+32 (hot+unstable), not 512 (blob, indistinguishable from
  clean) or 64 (never set). `_DQ_GOOD='512'` fixed it (single-frame defects 28%→5% of peaks
  outside r>3″, all real sources retained). Residual quadruplets in the *noise* map are
  correct (masking 1 of 4 frames raises σ by √(4/3)=15.5%); only more dither positions
  remove the genuine ones.
- **`--dq-refine`** (`refine_dq_flags()`, default 3σ, 0 disables) clears bits 8/16/32 on
  pixels not deviant from a 5×5 local median *in the same frame* — the dark-ref flags are
  ~60% unjustified per exposure — and clears stale DQ 4096. Runs on the copies in
  `data/drizzle_files/`, never `data/calibrated/`. Masked pixels 1.96%→0.82%/frame; judge
  speckle on **blank-sky clump count only** (ERR weighting legitimately raises σ on sources;
  the weight map is a continuum, so binning into coverage steps overstates the damage).

ACS shows the same physics as diagonal **stripes**, not dots: it drizzles native 0.05″
(one masked input px → one output px, vs 4.57 for F160W) and deliberately keeps
stable-hot/warm pixels as good; what it *does* mask that replicates is bad columns, which
`final_rot=0.0` rotates into 4 parallel diagonal stripes in the noise map. → memory:
acs_bad_column_stripes

## BELLS GALLERY: WFC3/UVIS reduction (`scripts/drizzle_wfc3_uvis.py`)

The `gallery` sample (15 lenses, Shu et al. 2016 Table 1 class E-S-A) is BELLS GALLERY
(props 14189, 16734), imaged in WFC3/UVIS across five filters: **F225W, F275W, F438W,
F606W, F814W**. One filter-agnostic script (`--filt`), modeled closely on
`drizzle_acs_wfc.py` — UVIS is a two-CCD optical detector like ACS/WFC, so it inherits the
same alignment/CR reasoning, not yet independently re-audited on gallery data (a per-lens
`--align` override table exists in the script but is empty so far; an explicit `--align`
always wins).

```bash
conda run -n stenv python scripts/drizzle_wfc3_uvis.py --lens J1110+3649 --filt f606W
bash scripts/run_gallery_uvis_all.sh                   # all 15 lenses, all 5 filters
```

`run_gallery_uvis_all.sh` iterates `f606W f814W f438W f275W f225W` in that order — F606W
first (the primary band and cutout `--center-band` proxy target), F814W second and before
the UV filters (see below), and does **not** `rm` the output dir first (unlike
`run_acs_all.sh`), relying on the idempotent skip so a multi-GB campaign is resumable.

- **Alignment default `mast`** — same reasoning as ACS/WFC3 (*WCS alignment* above):
  TweakReg erases the dither. Delivered WCS fit type is **not uniform across gallery
  observations** — bare IDCTAB (no absolute fit), GAIA eDR3, and GSC240 all occur, varying
  per lens *and* per filter — but relative (frame-to-frame) alignment is good regardless of
  which one a given exposure carries, verified by the registration-QC block passing on
  every product drizzled so far. No cross-band tie like `align_wfpc2_to_acs.py` is needed:
  a lens's bands share one field, so any absolute offset is common across them and they
  co-register with each other. → memory: gallery_uvis_idctab_only_wcs
- **Native pixel scale, `pixfrac=0.7`** (`FINAL_SCALE=0.0396`, chosen from a pixfrac scan
  on J1110+3649 F606W trading correlated noise against weight-map uniformity) — the
  opposite lever from the oversampled F606W (WFPC2)/F160W bands, which chose `pixfrac=1.0`
  to *reduce* correlation; UVIS is already at native scale, where a smaller drop shrinks
  the input footprint instead of opening coverage holes. Residual correlation is higher
  than native ACS (~1.5–1.6× integrated vs ~1.24) — use `make_cutouts.py --corr-factor
  ~1.6` for a diagonal-covariance likelihood. → memory: gallery_uvis_pixfrac
- **CR pass, ERR weighting**: same LACosmic-then-plain-mean route as ACS/WFPC2
  (`resetbits=0` on the CR pass, `4096` on no-CR), same `--wht-type ERR` with `K=1` (UVIS
  FLC ERR is in electrons, like ACS, not electrons/s like WFC3/IR).
- **UV bands show almost no deflector.** F225W (and to a lesser extent F275W) deflectors
  are early-type and UV-dark, while the lensed source stays bright (often star-forming) —
  a self-centred cutout recentre locks onto the arc or blank sky, not the galaxy. Always
  cut UV bands with the default `--center-band f814W` in `make_cutouts.py`, never
  `--center-self`, and make sure F814W is drizzled first (`run_gallery_uvis_all.sh` already
  orders filters this way). → memory: gallery_uvis_uv_deflector_faint
- **No PSF support yet.** `make_psf.py`/`psf_models.py` key off `ACS/WFC`, `WFC3/IR`, and
  `WFPC2` — WFC3/UVIS isn't wired in, so gallery has no `psf_kernel.fits` /
  `cutout_[cr_]psf.fits` products and no `info/lens_psf.json` entries.

Current state (reduced 2026-07-29): all 15 lenses have F606W; F814W/F438W on 6 each,
F275W on 5, F225W on 1 (J2342-0120) — matches the sparse per-lens filter coverage BELLS
GALLERY actually has on MAST, not a pipeline gap. `run_cutouts_all.sh` was extended to glob
the UV/blue bands (`f438W f275W f225W`) alongside the SLACS filters so it stays one runner
for every sample.

## QC mosaics (`scripts/make_mosaics.py`, `scripts/make_psf_mosaics.py`)

```bash
conda run -n stenv python scripts/make_mosaics.py --sample slacs_gold
conda run -n stenv python scripts/make_psf_mosaics.py --sample slacs_gold
```

Read-only QC: tile every lens's existing cutout (or trimmed PSF kernel) onto a 5-wide grid,
one mosaic per (filter group, panel type), written to `data/mosaics/<sample>/`. Nothing is
re-drizzled, re-cut, or rebuilt — pure visualization over what's already on disk in
`data/cutouts/`.

- **Filter groups** come from `scripts/mosaic_groups.py`, shared by both scripts so they
  stay in sync. For `slacs_gold`/`slacs_other`: `f814W`, `f606W_f555W` (WFPC2 F606W —
  including the split-visit `f606W_v1`/`f606W_v2` keys — merged with ACS F555W per lens,
  since no SLACS lens has both), `f160W`. For `gallery`: one group per UVIS filter (no
  cross-filter merging). A sample not listed there falls back to one group per filter
  subdirectory found on disk, so a new sample/filter still produces mosaics with no code
  change.
- `make_mosaics.py` writes `{group}_signal.png` / `{group}_noise.png` / `{group}_snr.png`
  (inferno + asinh stretch, the astropy convention for smoothly showing negative
  background-noise pixels alongside bright cores). Multi-instrument groups
  (`f606W_f555W`) split the colourbar per instrument — WFPC2 F606W is far noisier than ACS
  F555W (both in raw flux scale and true SNR), so a shared scale washes one out.
- `make_psf_mosaics.py` writes `{group}_psf.png`, peak-normalised per panel (kernels are
  unit-sum normalised by `trim_kernel_to_amplitude`, so raw peak reflects kernel *size* as
  much as sharpness), each panel tagged `emp` (empirical ePSF) or `mod` (STDPSF /
  focus-diverse / MAST PSF DB) from the `PSFMETH` keyword. Only `slacs_gold` has PSF
  products today; `slacs_other`/`gallery` mosaics for the `_psf` panel type will appear
  once `run_psf_all.sh` is run for those samples (no code change needed for `slacs_other`,
  same three instruments; `gallery` needs WFC3/UVIS PSF support first).

## NICMOS is deprioritised

Do not generate or propose NICMOS (NIC2) products unless explicitly asked — the FOV is far
too small (~19″ vs ~139″ for WFC3/IR) and the pipeline may be unsound. Answer F160W coverage
from WFC3/IR. All NICMOS data was deleted 2026-07-21 (472 MB); re-runnable via
`scripts/stale_scripts/drizzle_nic2.py` (raises `NotImplementedError` on import unless
`ALLOW_NICMOS=1`), which re-downloads from MAST. → memory: feedback_ignore_nicmos,
nicmos_tweakreg_tuning
