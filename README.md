# HST lens image reduction

Downloads calibrated HST exposures from MAST and drizzles them into science-ready
mosaics + cutouts for gravitational-lens samples (SLACS, BELLS). One product per
**lens + filter**.

> This README is the operator's quick-start: what to run, in what order, and how the
> instruments differ. `CLAUDE.md` is the deep reference — every non-obvious decision and
> the trap it avoids. When in doubt about *why*, read that; when the doc and the code
> disagree, trust the code.

## Setup

All scripts run inside the `stenv` conda environment (the STScI pipeline stack):

```bash
conda activate stenv
# or, per-command:
conda run -n stenv python scripts/<script>.py --lens <LENS> --filt <FILTER>
```

CRDS reference files download themselves on first use into `data/reference_files/`.

## The pipeline, end to end

For each lens+filter, a drizzle script does:

1. **Download** calibrated exposures from MAST (skipped if already present).
2. **Download** CRDS reference files, run `bestrefs`.
3. **Drizzle** — align + combine with AstroDrizzle onto a common North-up grid.
4. Produce a **no-CR mosaic**, plus a **CR-rejected** one (WFPC2 always; others with `--cr`).
5. Update three tracking JSONs in `info/`.

Then, across bands:

6. `align_wfpc2_to_acs.py` — tie F606W absolute astrometry to ACS (**F606W only**).
7. `make_cutouts.py` — cut co-registered stamps + noise maps + preview PNGs.

**Run order matters: ACS + WFPC2 drizzles → `align_wfpc2_to_acs.py` → `make_cutouts.py`.**
A re-drizzle discards the F606W astrometric tie, so the align step must run again after it.

Scripts are **idempotent**: they skip the MAST download if calibrated files exist, and skip
the whole drizzle if the final product exists. To force a re-run, delete the lens's dir
under `data/drizzled/` (and `data/drizzle_files/`).

## Running things

### One lens

```bash
conda run -n stenv python scripts/drizzle_acs_wfc.py   --lens J0008-0004 --filt f814W
conda run -n stenv python scripts/drizzle_acs_wfc.py   --lens J0216-0813 --filt f555W
conda run -n stenv python scripts/drizzle_wfpc2_wf3.py --lens J0008-0004 --filt f606W
conda run -n stenv python scripts/drizzle_wfc3_ir.py   --lens J0008-0004 --filt f160W
conda run -n stenv python scripts/make_cutouts.py      --lens J0029-0055 --filt f606W
```

`--sample` defaults to **`slacs_gold`** everywhere. It sets the `<sample>` level of every
`data/` path, so passing the wrong value silently writes a good product into the wrong tree.

A lens with no data for the requested instrument+filter prints `=== NO DATA: ...`, records
`null` in the tracking JSONs, and **exits 0** — no-data is a normal outcome, not a failure.

### All lenses

```bash
bash scripts/run_acs_all.sh                  # ACS/WFC  F814W + F555W
bash scripts/run_wfc3_all.sh                 # WFC3/IR  F160W
bash scripts/run_wfpc2_wf3.sh                # WFPC2/WF3 F606W: drizzle -> align -> cutout
bash scripts/run_cutouts_all.sh              # stamps for whatever products exist
bash scripts/run_acs_all.sh slacs_other      # any runner takes an optional sample arg
```

Runners take the roster from `info/lens_samples.json` (via `mast_target_names.py`) and try
**every** lens — only ones with data download. They report `ok` / `no data` / `FAILED`
separately. Exception: `run_cutouts_all.sh` globs `data/drizzled/` instead (a stamp needs a
mosaic that already exists).

**`run_wfpc2_wf3.sh` is the one WFPC2 driver** — it runs the full three-stage order per lens,
reads per-lens alignment from `info/wfpc2_alignment.json`, expands the two split-visit lenses
into per-visit products, and retries failures once. Don't drive the WFPC2 stages by hand
unless you know why; the runner exists to avoid three silent-wrong-product traps (align
default, split visits, the mandatory align step).

## Directory layout

```
data/
  calibrated/<sample>/<lens>/<filter>/    downloaded FLT/FLC/CAL files
  drizzle_files/<sample>/<lens>/<filter>/ AstroDrizzle working dir (logs, shift files, PNGs)
  drizzled/<sample>/<lens>/<filter>/      final mosaics (_cr_ / _nocrrej_, sci + wht)
  cutouts/<sample>/<lens>/<filter>/       cutout_sci.fits / cutout_noise.fits / cutout.png
  reference_files/                        CRDS cache (auto-downloaded once)
  run_logs/                               per-lens batch-runner logs
info/
  lens_samples.json     single source of truth for sample membership + MAST quirks
  wfpc2_alignment.json  per-lens WFPC2 align mode (mast) and split-visit handling
  slacs_coords.py       catalogue positions -> common output WCS
  lens_products.json / lens_instrument.json / lens_exptime.json   auto-updated tracking
```

## The four instruments — how they differ

| Script | Input | Detector | Native → output scale | Suffix | Default CR |
|---|---|---|---|---|---|
| `drizzle_acs_wfc.py`  | `*flc.fits` | ACS/WFC   | 0.05″ → 0.05″     | `_drc_` | off (`--cr` to enable) |
| `drizzle_wfc3_ir.py`  | `*flt.fits` | WFC3/IR   | 0.1283″ → 0.06″   | `_drz_` | **never add CR** |
| `drizzle_wfpc2_wf3.py`| `u*flt.fits`| WFPC2/WF3 | 0.0996″ → 0.05″   | `_drw_` | always (LACosmic) |

The suffix is set by **input file type**, not output name (`_drc_` FLC, `_drw_` WFPC2 FLT,
`_drz_` everything else). Key per-instrument differences an operator should know:

### ACS/WFC — F814W, F555W (the workhorse)
- Native 0.05″, no oversampling; cleanest per-pixel noise model of the three.
- `--align mast` (default): trusts the delivered GAIA/GSC WCS, no re-solve. **Do not run
  TweakReg** — on dithered frames it removes the dither and smears the stack.
- Noise: `--wht-type ERR` (default) → correct per-pixel σ.
- One filter-agnostic script covers both F814W and F555W.

### WFC3/IR — F160W
- Drizzled to **0.06″, pixfrac 1.0** — settled by measurement; do not re-drizzle to 0.05″
  to grid-match other bands (PyAutoLens fits inter-band offsets itself).
- **No CR pass, ever** — the ~1px IR PSF looks like an outlier, so LACosmic/driz_cr eat point
  sources. FLTs are already up-the-ramp CR-rejected.
- DQ subtlety: treat **only** bit `512` as good, written as `'512'` — **never `''`** (empty
  string disables masking entirely and quadruples detector defects).
- ERR is in ELECTRONS/**S**, so the noise-map scaling needs a per-frame EXPTIME factor
  (handled by `make_cutouts.py`); the raw `1/sqrt(WHT)` would be ~600× too small.

### WFPC2/WF3 — F606W (the fiddly one)
- **The lens sits on WF3 (ext 3), not the PC** — `DETECTOR=PC` in the header names the
  aperture, not the chip. The script extracts the WF3 chip first.
- **Two noise fixes, both required**: DrizzlePac ignores ERR for WFPC2 (hand-built IVM fed
  via an `@`-file), and the WFPC2 ERR array is bogus (`ERR==sqrt(SCI)`), so the model is
  `var = SCI/gain + floor²`. Products stamp `IVMMODEL`; a pre-fix product needs a re-drizzle.
- **`--align mast` runs `updatewcs` but NOT TweakReg** — it's the only instrument that can't
  skip `updatewcs` (needs the distortion arrays), but TweakReg still scatters its
  mostly-single-visit frames into a multi-knot core.
- **Absolute astrometry is off** (~0.5–0.9″); `align_wfpc2_to_acs.py` ties it to ACS F814W
  afterwards. This is why the align step is mandatory and non-optional.
- **Two lenses are split per visit** (J0728+3835, J0822+2652) → `f606W_v1`/`f606W_v2` product
  dirs, not a combined `f606W`. The runner handles this.
- Units: SCI is DN/s (`BUNIT='COUNTS/S'`), F814W is e/s — cross-band flux comparison must go
  through `PHOTFLAM`, never raw pixel values.

### NICMOS — deprioritised
Do not generate or propose NICMOS products unless explicitly asked. FOV is too small and the
pipeline may be unsound; answer F160W coverage from WFC3/IR. `scripts/stale_scripts/` holds
retired scripts that refuse to run without an override env var — leave them alone.

## Samples & coverage

`info/lens_samples.json` is the single source of truth. Query it via:

```bash
conda run -n stenv python scripts/mast_target_names.py --list      # samples + sizes
conda run -n stenv python scripts/mast_target_names.py slacs_gold  # lens names
```

| Sample | Lenses | Status |
|---|---|---|
| `slacs_gold`  | 38 | Working sample; default `--sample`. F814W (38), F606W (22), F555W (16), F160W (13) |
| `slacs_other` | 93 | Rest of SLACS; not yet reduced |
| `gallery`     | 16 | BELLS GALLERY; scripts not yet written |

## Tracking JSONs

Every run updates three files in `info/`, keyed `{lens: {product_key: value}}`:
`lens_products.json` (drizzled frame rootnames), `lens_instrument.json` (`INSTRUME/DETECTOR`),
`lens_exptime.json` (seconds). No data → `null`. The key is the **product directory**, so
split-visit lenses are keyed `f606W_v1`/`f606W_v2` with no bare `f606W`.
