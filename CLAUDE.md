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
4. Produce a no-CR-rejection mosaic, plus a CR-rejected one (LACosmic) for WFPC2 or when
   `--cr` is given.
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

ACS/WFC3/NICMOS take `--cr` to enable the CR-rejection pass (off by default). WFPC2 always
runs both passes.

## Running all lenses

Each runner takes an optional sample arg, defaulting to `slacs_gold`:

```bash
bash scripts/run_acs_all.sh                  # ACS/WFC F814W + F555W
bash scripts/run_wfc3_all.sh                 # WFC3/IR F160W
bash scripts/run_wfpc2_wf3.sh                # WFPC2/WF3 F606W: drizzle -> align -> cutout
bash scripts/run_cutouts_all.sh              # stamps for whatever products exist
bash scripts/run_acs_all.sh slacs_other      # any runner, any sample
```

All runners take the roster from `info/lens_samples.json` via `scripts/mast_target_names.py`
— **except `run_cutouts_all.sh`, which globs `data/drizzled/`** (a stamp needs a mosaic
that exists). They report `ok` / `no data` / `FAILED` separately, so the 16 `slacs_gold`
lenses with no WFPC2 data aren't mistaken for errors. `run_cutouts_all.sh` globs `<filt>*`
(not `<filt>`) so per-visit split-visit dirs are included.

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

## Data flow and directory layout

```
data/
  calibrated/<sample>/<lens>/<filter>/    ← downloaded FLT/FLC/CAL files
  drizzle_files/<sample>/<lens>/<filter>/ ← working dir; AstroDrizzle runs here (run.log, shift_*.txt, *_single_*.fits, *.png)
  drizzled/<sample>/<lens>/<filter>/      ← final products (<prefix>_cr_*/_nocrrej_* sci+wht)
  cutouts/<sample>/<lens>/<filter>/       ← cutout_sci.fits / cutout_noise.fits / cutout.png
  run_logs/                               ← per-lens batch-runner logs
  reference_files/                        ← CRDS reference files (auto-downloaded once)
```

## Instrument-specific scripts

| Script | Input | MAST product | Ref env | Pixel scale | Suffix |
|---|---|---|---|---|---|
| `drizzle_wfpc2_wf3.py` | `u*flt.fits` | FLT / CALWFPC2 | `uref` | 0.0996″ → 0.05″ | `_drw_` |
| `drizzle_acs_wfc.py`  | `*flc.fits`  | FLC / CALACS   | `jref` | 0.05″ | `_drc_` |
| `drizzle_wfc3_ir.py`  | `*flt.fits`  | FLT / CALWF3   | `iref` | 0.1283″ → 0.06″ | `_drz_` |
| `drizzle_nic2.py`     | `*cal.fits`  | CAL / CALNIC   | `nref` | 0.0756″ | `_drz_` |

The output suffix is set by **input file type**, not output name: `_drc_` for FLC (ACS),
`_drw_` for WFPC2 FLT, `_drz_` for everything else. The WFPC2 script extracts only the WF3
chip (SCI/ERR/DQ ext 3) into `wf3_`-prefixed files first; the others are MEF files
DrizzlePac handles natively. → memory: instrument_drizzle_ref, crds_bestrefs_always_run
(never skip `bestrefs` when the ref dir is non-empty).

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
lens catalogue position (from `info/slacs_coords.py`) go into every AstroDrizzle call. A
lens missing from that table falls back to the native drizzle WCS with a warning. This
aligns orientation + tangent point, not absolute astrometry (still the delivered WCS's job,
and `align_wfpc2_to_acs.py`'s for F606W).

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
band-dependent (~1.24 ACS, ~1.17 F160W, ~1.1 F606W), so ignoring it also mis-weights bands
relative to each other in a joint fit. `make_cutouts.py --corr-factor` applies it (default
1.0, off). Prefer native-scale F814W where a clean per-pixel noise model matters most.

## Cutouts (`scripts/make_cutouts.py`)

```bash
conda run -n stenv python scripts/make_cutouts.py --lens J0029-0055 --filt f606W
```

Cuts a square stamp (default 20″) from `data/drizzled/` into `data/cutouts/`: a sci FITS, a
noise FITS (from the weight map), and a 3-panel PNG.

- **Pass + prefix.** `--pass {auto,cr,nocrrej}` (default `auto`) picks the CR pass when one
  exists, else no-CR. The prefix encodes it so the two coexist: **`cutout_cr_*` for CR,
  `cutout_*` for no-CR.** So WFPC2 F606W (always has a LACosmic CR pass) → `auto` cuts
  `cutout_cr_*`, which is science-grade (LACosmic preserves the core). ACS/WFC3 default to
  no CR pass → `cutout_*`. → memory: cutout_centring_on_cr_pass
- **Recentring.** The stamp recentres on the galaxy, and the peak search prefers the CR
  mosaic — a brightest-pixel search on a CR-laden no-CR mosaic locks onto cosmic rays. Don't
  reach for `--median-size` for a bad recentre; check a `*_cr_*` mosaic is present.
  Offsets of 1–2″ aren't always failures — some lenses (J0912+0029, J0956+5100) have genuine
  multi-knot morphology.
- **Shared centre (`--center-band`, default `f814W`).** All bands are cut about a single
  centre from the `--center-band` mosaic (highest S/N, GAIA-accurate) so stamps co-register
  across filters. `--center-self` restores per-band recentring; a band whose center-band
  products are missing falls back to its own peak with a warning.

## Tracking JSONs in `info/`

> **Full reset, 2026-07-26.** All three files were emptied to `{}` and every product under
> `data/` deleted (`calibrated/`, `drizzle_files/`, `drizzled/`, `cutouts/`, `mosaics/`,
> `run_logs/`) as a deliberate clean restart. Kept: `data/reference_files/` (CRDS cache) and
> `data/pre_drizzled/` (46 MAST-delivered mosaics, not pipeline output). The sample was
> renamed `slacs` → **`slacs_gold`** in the same pass. **Any surviving `data/*/slacs/` path,
> and any "current on-disk state" / "Not regenerated" claim in memory, is pre-reset and
> void** — the *reasoning* in those notes stands and is why reruns use the current scripts;
> only the inventory is stale.

Updated automatically by every run:
- **`lens_products.json`** — `{lens: {key: [rootname, ...]}}` — frames that reached the
  drizzle (not the whole download).
- **`lens_instrument.json`** — `{lens: {key: "INSTRUME/DETECTOR"}}` — records `WFPC2/WF3`
  for F606W (the chip), not MAST's `WFPC2/PC`.
- **`lens_exptime.json`** — `{lens: {key: seconds}}` — from the CR-rejected drizzle header.

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
- Both levels stay sorted (lenses, and filters within each lens) across partial runs.

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
| **`slacs_other`** | 93 | The rest of SLACS (Bolton et al. 2008 Table 4). Not yet reduced |
| **`gallery`** | 16 | BELLS GALLERY (props 14189, 16734); WFC3/UVIS multi-band. Scripts not written |

`slacs_gold` coverage: **F814W** (ACS/WFC, all 38), **F606W** (WFPC2/WF3, 22 — the other 16
have no WFPC2 data), **F555W** (ACS/WFC, 16 — exactly those 16), **F160W** (WFC3/IR, 13 —
prop 11202). NICMOS F160W (24 lenses) is deprioritised and its data deleted. HST props:
10886, 11202, 10494, 10798.

Caveats on the other samples: `slacs_other` has **no `GAL-*` names surveyed**, so "no
observations" there may just mean the lens is archived under an unlooked-up designation;
`gallery` lenses are **not in `info/slacs_coords.py`** so they get no common output WCS
(native-WCS fallback with a warning). Both need work before reduction.

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
those, not any table. All 14 are in `slacs_gold`; none surveyed for the other samples (which
is why a no-data result there isn't conclusive).

## AstroDrizzle key parameters

Two passes over the same inputs:
- **CR pass**: default **LACosmic** — mask CRs per frame, then a plain-mean drizzle
  (`median=False, blot=False, driz_cr=False`, `resetbits=0`). `--cr-method drizcr` restores
  the AstroDrizzle route (which eats the core).
- **No-CR pass** (`median=False, blot=False, driz_cr=False`): uncleaned, for comparison.

Only WFPC2 runs both unconditionally; ACS/WFC3/NICMOS run no-CR only unless `--cr`. WFC3/IR
F160W therefore has no CR pass — `make_cutouts.py` falls back to the science pass for
recentring (acceptable: FLTs are already up-the-ramp CR-rejected; re-run with `--cr` if a
recentre looks wrong).

**DQ bits treated as good** (do not unify these — they encode different detector facts):
- WFPC2: `8,1024`
- ACS/WFC: `256,64,16` (saturated, warm, stable hot)
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

## NICMOS is deprioritised

Do not generate or propose NICMOS (NIC2) products unless explicitly asked — the FOV is far
too small (~19″ vs ~139″ for WFC3/IR) and the pipeline may be unsound. Answer F160W coverage
from WFC3/IR. All NICMOS data was deleted 2026-07-21 (472 MB); re-runnable via
`scripts/stale_scripts/drizzle_nic2.py` (raises `NotImplementedError` on import unless
`ALLOW_NICMOS=1`), which re-downloads from MAST. → memory: feedback_ignore_nicmos,
nicmos_tweakreg_tuning
