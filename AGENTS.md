# AGENTS.md

Guidance for any coding agent (and for humans) working in this repo — agent-agnostic;
`CLAUDE.md` is a pointer to this file. This file is the **operational**
reference: commands, layout, decisions, and the traps that silently produce wrong
products. The *evidence* behind each decision (measurement tables, derivations) lives in
the agent's persistent memory store — pointers below read `→ memory: <slug>`. Don't
re-derive what a pointer already settles; don't restate the tables here. An agent with no
such store should treat a `→ memory:` pointer as "this was measured and settled, don't
re-litigate it"; the slug names the note, not a file in this repo.

> **AGENTS.md documents intent, not always state.** This doc has repeatedly described
> decisions the code never implemented (argparse defaults, runners reading a JSON, a
> pipeline stage). Before running any batch campaign or reporting what's outstanding,
> verify against source — argparse defaults, what the runner actually passes, tracking
> JSONs vs what's on disk.

> **bcfill reverted; the standard drizzle is the modelling input (2026-09-22).** The
> bad-column fill only removes the *honest* coverage stripe from the noise map (see the
> decision under *Bad-column fill*), so `main` now carries the **standard** stamps for all
> 105 ACS/WFPC2 bands that had been bcfill, and the bcfill products live on the **`bcfill`
> branch** for collaborators. To make that split safe, **the stamp filename now carries the
> reduction**: `cutout_[cr_]{sci,noise}.fits` is standard, `cutout_cr_bcfill_sci.fits` /
> `cutout_cr_bcfill_crfill_sci.fits` are variants (`cutout_paths.stamp_name`). Masks, arc
> masks, positions and PSFs keep the bare `cutout_[cr_]` prefix — they describe the grid,
> which is identical across reductions — so a mask commit on `main` cherry-picks onto
> `bcfill` without dragging a stamp along, and a merge that lands two reductions in one dir
> makes every reader raise `AmbiguousStampError` instead of silently picking one. Readers
> resolve stamps through `cutout_paths.find_stamp` / `find_prefix`, never by spelling the
> name; `make_cutouts.py` refuses to put a second reduction beside an existing one unless
> `--force` (which then removes the other's sci/noise/png). The `info/lens_cutout_qc.json`
> record gained `stamp` (the filename) and `variant`. The 2026-09-14 "one science tree"
> merge otherwise stands: one dir per band, no `cutouts_bcfill/` tree.
>
> **Bolton-interpolation investigation → bcfill productionized.**
> `scripts/bolton_investigations/` (+ `README.md`) holds the standalone exploratory scripts for
> the ACS dead-column noise-stripe question (Bolton-2008 bilinear reduction, post-hoc stripe
> heal/mask, input-level bad-column fill + re-drizzle). The investigation settled on **option 3
> (input-level bad-column fill), now productionized as `--bcfill`** — see *Bad-column fill*
> below. The other scripts stay standalone validation-only tools writing to
> `diagnostics/bolton_test_outputs/`, *not* `data/`. **That folder was deleted 2026-09-13**
> (24 files, 31 MB), closing the "keep it long-term?" question: the investigation is settled
> and bcfill is in production, so the demonstrator figures no longer earn their place in the
> repo. Every script recreates the folder on run (`makedirs(..., exist_ok=True)`), and the
> deleted figures are recoverable from git history if a claim ever needs re-checking.

## Environment

All scripts run inside a **uv-managed** virtual environment (`pyproject.toml` + `uv.lock`
+ `.python-version`, all tracked in git), which replaced the `stenv` conda environment
(2026-07-30):

```bash
uv sync                                               # one-time, or after a .venv wipe;
                                                       # downloads the pinned CPython itself
uv run python scripts/<script>.py --lens <LENS> --filt <FILTER>
```

Interactive: `source .venv/bin/activate`.

**The repo root is derived from `__file__`, never hard-coded.** Every script sets
`ws_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))` (three `dirname`s
in `scripts/stale_scripts/` and `scripts/bolton_investigations/`), and the shell runners
take `WS="$(dirname "$SD")"` off `BASH_SOURCE`, so the checkout can be moved or cloned
anywhere with no edit — done 2026-09-08, when an absolute `/Users/.../HST_lenses_data_reduction`
literal in 18 scripts broke on a move. Don't reintroduce one. A move *does* invalidate
`.venv` (its `activate` and console-script shebangs carry the old absolute path, and
`uv sync` alone won't rewrite them): `rm -rf .venv && uv sync`. Historical logs under
`data/run_logs/` and `data/drizzle_files/*/run.log` still quote the old path — correct as
records of where they ran, deliberately left alone.

**Fully uv-managed, cross-platform reproducible — deliberately, 2026-07-30.** `.python-version`
pins the exact stenv interpreter version (`3.12.13`); `pyproject.toml` sets
`tool.uv.python-preference = "only-managed"` so `uv sync` always downloads that exact
CPython build itself (via `python-build-standalone`) rather than reusing whatever system
Python happens to be on `PATH` — no external installer, no conda, works identically on
macOS (arm64 or x86_64) and Linux with the one command above. The only 12 packages
actually imported by the pipeline (astropy, astroquery, drizzlepac, stwcs, photutils,
matplotlib, numpy, scipy, acstools, astroscrappy, requests, crds) are pinned in
`pyproject.toml` to their exact stenv versions; `uv.lock` resolves the rest of the tree
fresh (not a byte-for-byte freeze of all ~250 stenv packages — the unused
jupyter/dask/ginga/easyocr/torch bulk of the conda env was dropped as dead weight). All
have Linux + macOS (arm64/x86_64) wheels on PyPI, verified at conversion time.

**Superseded x86_64 pin (2026-07-30):** an earlier version pinned x86_64 (all numeric
results were originally measured under stenv/Rosetta), needing a manual python.org
bootstrap. Abandoned for cross-platform reproducibility once confirmed the codebase has no
architecture-specific logic — the one platform-sensitive piece (`mmap_fits_write.py`) gates
on `sys.platform == 'darwin'`, not architecture. Native arm64/Linux is an accepted, not
fully re-verified, change — spot-check a known product against a prior measurement if
something looks numerically off.

## Long batch runs: `find -newermt` is NOT a safe completion test (2026-09-16)

Chaining batch stages on a "has stage 1 finished?" poll loop cost a **10-hour idle overnight**
here, and the trap is entirely portable. On this Mac the interactive shell's `find` is a
**function shimming `bfs`**, which parses `-newermt '2026-09-15T21:10:00'` as **local time**; a
script gets `/usr/bin/find` (BSD), which parses the same string as **UTC**. The same command
counted **90** files by hand and **44** inside the script, so a condition verified
interactively could never be true where it actually ran — and the loop slept on, indefinitely.
Two rules:
- **Key a wait on the runner's own exit**, not on a reconstructed timestamp condition
  (`cmd_a && cmd_b`, or wait on the PID). If a count is unavoidable, count files without a
  time filter and compare against the expected total.
- **A poll loop that can never satisfy its condition is indistinguishable from one still
  working** — same silence, same live process. Any watcher on a long job must emit on failure
  and timeout, not only on success, or "still running" will be the last thing you hear.

## macOS write-hang workaround (required — keep this wiring)

On this Mac, AstroDrizzle's large buffered FITS writes hit a kernel lost-wakeup
(`tofile → write() → cluster_write → copyin → lck_rw_sleep`) that wedges the process into
an unkillable U-state (only a reboot clears it). `scripts/mmap_fits_write.py` monkeypatches
`astropy.io.fits.file._array_to_file` to write via `mmap`+`memcpy` (the `vm_fault` path),
dodging the hang; it is a no-op off macOS and byte-identical to stock astropy. Every
drizzle script imports it and calls `install()` before AstroDrizzle. `DRIZZLE_MMAP_DEBUG=1`
logs each mmap write. **`num_cores=1` in all scripts is a related, separate requirement** —
parallel `fork` triggers the same U-state.

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
uv run python scripts/drizzle_wfpc2_wf3.py --lens J0008-0004 --filt f606W
uv run python scripts/drizzle_acs_wfc.py  --lens J0008-0004 --filt f814W
uv run python scripts/drizzle_wfc3_ir.py  --lens J0008-0004 --filt f160W
uv run python scripts/drizzle_acs_wfc.py  --lens J0216-0813 --filt f555W
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
bash scripts/run_cutouts_all.sh              # stamps for whatever products exist (12", default)
bash scripts/run_cutouts_all.sh slacs_gold 20  # same, at a 20" stamp size (parallel tree)
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

(All three were once documented in AGENTS.md but unimplemented — a concrete instance of
the "documents intent, not state" warning at the top of this file.)

### `scripts/stale_scripts/` — retained, but nothing invokes them

- `drizzle_wfpc2_pc.py` — **superseded**: extracts the wrong chip *and* `rmtree`s the good
  WF3 products. Override `ALLOW_SUPERSEDED_WFPC2_PC=1`.
- `drizzle_nic2.py` — **deprioritised**: an accidental run re-downloads ~472 MB and
  repopulates NICMOS entries. Override `ALLOW_NICMOS=1`.
- `run_all_lenses.sh` — retired WFPC2 driver; refuses to run (its retry pass moved into
  `run_wfpc2_wf3.sh`).

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
  drizzled[_bcfill]/<sample>/<lens>/<filter>/ ← final products (<prefix>_cr_*/_nocrrej_* sci+wht); _bcfill is the bad-column-filled re-drizzle
  cutouts/<sample>/<lens>/<filter>/       ← THE science tree: one stamp per band — cutout_[cr_]{sci,noise,psf,psf_err,mask}.fits + cutout_[cr_][_dataset].png
  cutouts_<S>arcsec/<sample>/...          ← the same stamps cut at a non-default --size (see Cutouts)
  psf/<sample>/<lens>/<filter>/           ← archival PSF products (psf_kernel.fits / psf.png; model-tier also carries psf_kernel_analytic.fits / psf_analytic.png)
  mosaics/<sample>/                       ← QC mosaics tiling every lens's cutouts/PSFs (make_mosaics.py, make_psf_mosaics.py)
  mosaics_<S>arcsec/<sample>/             ← QC mosaics of the <S>" stamps
  run_logs/                               ← per-lens batch-runner logs
  reference_files/                        ← CRDS reference files (auto-downloaded once)

diagnostics/                            ← QC + investigation figures, GITIGNORED IN FULL
  arc_detection/<sample>/               ← the only subfolder a script writes (detect_arcs.py)
  masks/                                ← is the contaminant-mask boundary right?
  feature_vetting/                      ← is this feature lensed?
  cr_defects/                           ← instrumental: CRs, LACosmic scans, dead columns
  gui/                                  ← screenshots of the interactive tools
```

Nothing under `diagnostics/` survives a clone, so every figure there must be regenerable and
the **conclusions belong in this file, not in the pixels** — never cite a figure without the
AGENTS.md note that goes with it. It was a flat listing until **2026-09-21**; the per-file
convention `<lens>_<filt>_<topic>.png` is unchanged, and `diagnostics/README.md` restates the
split for anyone browsing the folder directly.

**What goes in which `diagnostics/` subfolder.** `arc_detection/` is the only one a script
owns — `detect_arcs.py` writes it (`--out-dir` defaults to `arc_detection/<sample>/`), so leave
its layout alone. The other four are for hand analysis and split **by the question asked**, not
by lens: `masks/` = *is the contaminant-mask boundary in the right place* (mask checks,
proposals, pull-backs); `feature_vetting/` = *is this feature lensed* (colour against the
deflector, counter-image searches, shear/elongation, diffuse-arc and galaxy-model residuals);
`cr_defects/` = *instrumental, not astrophysical* (cosmic-ray tracks and their weight/noise
residue, LACosmic scans, dead columns, per-frame drizzle comparisons); `gui/` = screenshots of
the interactive tools. A figure that answers two questions goes with the question it was made
to answer, and any AGENTS.md note that cites it names the full path.

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
files DrizzlePac handles natively. Never skip the CRDS `bestrefs` sync when the ref dir is
non-empty — references can update independently of what's already cached locally.

`drizzle_wfc3_uvis.py` is the BELLS GALLERY driver — see *BELLS GALLERY: WFC3/UVIS
reduction* below for its defaults, alignment, and current coverage.

### WFPC2: the lens is on WF3, not the PC

For all 22 SLACS WFPC2 F606W lenses the lens galaxy falls on **WF3** (ext 3) at
~(435, 424). `DETECTOR = PC` in the primary header (hence the `WFPC2/PC` MAST label) names
the *aperture*, not the chip; the full-field aperture centres the target on WF3. Only the
per-extension `DETECTOR` identifies chips (1=PC…4=WF4). The superseded `drizzle_wfpc2_pc.py`
extracted ext 1 → 22 blank-sky mosaics ~79″ off the lens.

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
scatter.

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
datasets (`--pa <PA_V3> --out-suffix ''/_v2`, each single-guide-star → `mast`); the
**longer-exptime visit is the primary product, keyed bare `f606W`** (no suffix), any shorter
visit `_v2` (`_v3`, ... if a lens ever has more). Outcome: **J0822+2652** = `f606W` (4×1100s)
+ `f606W_v2` (2×1100s); **J0728+3835** = `f606W` only (its other, 2-frame visit has 1
x-dither-phase, can't reach 0.05″, dropped — never gets a product or suffix at all);
**J1142+1001 stays combined** (visits share roll, PA 119.00 vs 118.87). Renamed 2026-08-03
from the original `f606W_v1`/`f606W_v2` (arbitrary visit order, no bare `f606W` key) to this
exptime-ranked scheme — `align_wfpc2_to_acs.py --f606-dir` and `make_cutouts.py --filt
f606W_v2` handle the suffixed dir; the primary needs no flag, same as any other filter.

**J1142+1001 is the only combined-visit lens** — verified (2026-08-03) across every WFPC2
F606W product in `slacs_gold` (38) and `slacs_other` (27): grouping `info/lens_products.json`
frame lists by 6-char rootname prefix, only J1142+1001's `f606W` mixes two visits (`ua1l38` +
`ua1lc8`, 6 frames); the split keys J0728+3835/J0822+2652 are each single-visit. The
shared-roll combine criterion is a registration-consistency proxy (core-registration scatter),
not a PSF/noise-homogeneity check — those are handled separately (per-frame IVM downweights a
noisier visit; the PSF is measured *after* combination from the actual stack). A one-off
spot-check put J1142's injected-kernel FWHM/pedestal/scatter within the single-visit spread —
no smearing penalty.

**Diagnosing this class of bug:** compare the *spread* of per-frame WCS error, not its
magnitude — a common offset is a harmless absolute-astrometry shift; frame-to-frame scatter
is what smears a stack. Cleanest metric: per-frame core-registration scatter on the
drizzled common grid, CRs masked; confirm against the visible product (one core, not
knots). Stacked stellar FWHM can mislead — it centroids an extended galaxy and rewards
TweakReg's self-consistency even when the deflector is split.

TweakReg `threshold` is in **image data units** and does not transfer between detectors
(WFPC2/WF3 100, WFC3/IR 20, ACS default). This only takes effect under `--align tweakreg`,
which no lens uses — dead code kept for comparison runs.

### `align_wfpc2_to_acs.py` — F606W absolute astrometry (after drizzle, before cutouts)

WFPC2 F606W carries only GSC 2.4.0, ~0.3–1″ off absolute; ACS F814W and WFC3/IR F160W carry
GAIA eDR3 / GSC242 (<0.02″, agree to ~0.01″). So F606W sits ~0.5–0.9″ off the other bands —
a whole-mosaic shift (harmless to the F606W stack, breaks cross-band registration). The
script ties the deflector light-centroid to ACS F814W via an iterative windowed
`centroid_com` (robust to the ring), shifts F606W `CRVAL1/2` so the centroids coincide, and
stamps `GSC240FX=True`. Idempotent; refuses any tie implying > `MAX_SHIFT = 1.5″`; uses
F160W as an *independent* check where present. Verified J0252+0039: 0.66″ → 0.009″.

```bash
uv run python scripts/align_wfpc2_to_acs.py --lens J0252+0039   # or --all
```

**Run order: ACS + WFPC2 drizzles → `align_wfpc2_to_acs.py` → `make_cutouts.py`.**

**`--target`/`--ref`: the mis-fit band is not always F606W.** A GAIA-grade solution is a
property of the *delivered WCS*, not the instrument, and can fail on ACS/WFC3-IR too.
J1016+3859 (`slacs_other`) came down with a bare `-GSC240` fit on its ACS F814W (the
`force_copy` COPY visit), 0.60″ off its own F160W; since `make_cutouts.py` centres every band
on `--center-band f814W`, that dragged all three stamps ~0.6–0.9″ off the deflector (the
reference-band trap). So `align_lens` takes `--target` (band whose CRVAL moves, default F606W)
and `--ref` (the frame, default F814W); products record `ASTROREF` alongside `GSC240FX`. Fixed
2026-08-04 by tying F814W→F160W then F606W→F814W (→ memory: acs_wcs_not_always_gaia). **Check
`WCSNAME` before assuming which band is truth** — a `-GSC240` suffix on ACS/WFC3-IR (vs
`-FIT_REL_GSC242`/`-GAIAeDR3`) means that band is the ~0.5″ one, whatever the instrument.
Defaults are unchanged and stay correct for every other lens.

**It varies per FILTER within one visit, and `WCSNAME` predicts it exactly — 4 `slacs_gold`
F555W products were tied 2026-09-05.** J1402+6321 (1.327″), J1205+4910 (0.605″), J1420+6019
(0.546″) and J0959+0410 (0.464″) each sat that far off their *own lens's* F814W despite being
the same day, same visit, same proposal (10494): MAST fitted an absolute solution to the F814W
exposures and not the F555W ones. **The rule is exact over all 16 F555W products: `WCSNAME`
carrying `FIT_REL_*`/`FIT_IMG_*` is fine (all 12 within 0.048″); a bare IDCTAB or a `-GSC240`
suffix is offset (all 4).** J1402+6321's was bare `IDC_4bb1536oj` — no absolute fit at all —
and its 1.33″ sits close to `MAX_SHIFT = 1.5″`. Fixed with the tool as-is (`--target f555W
--ref f814W`, run in **both** trees, `--bcfill` second), then `make_cutouts.py --lens L --filt
f555W` in both — the crop window moves, so the stamps must be re-cut. F606W (0/22, already
tied) and F160W (0/13, max 0.078″) were clean. **Symptom to recognise:** a mask reprojected
from F814W lands with the *right shape in the wrong place*, because `reproject_mask_bool`
correctly puts it at the right *sky* position while the band's own WCS points that sky
position at the wrong pixels — and for the same reason the stamp itself is cut off-centre
(J1402's deflector sat at pixel 145 rather than 120).
- **Verify a tie with the median VECTOR residual over field sources, never a scalar
  separation.** Post-tie cross-matched separations read 0.06–0.18″ median and look like a
  residual; they are cross-filter centroid scatter. The median (dRA·cos δ, dDec) came to
  0.006–0.028″ on all four, matching never-offset J0216-0813 (0.025″). Field sources are the
  right check because the tie uses only the deflector centroid, so they are independent of it.
- A whole-sample sweep of centroid offset vs `WCSNAME` costs one read-only survey; worth
  re-running after any re-drizzle, and before trusting a cross-band mask broadcast.

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

DrizzlePac computes `weight = (EXPTIME/ERR)²`:
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

`--cr` masks CRs per frame with LACosmic (default
`--cr-method lacosmic`) then drizzles a plain weighted mean (`median=False, blot=False,
driz_cr=False`); `--cr-method drizcr` restores the old route. WFPC2 uses the same LACosmic
route (`--lacosmic-sigclip 4.5 --lacosmic-objlim 5.0`; gain/readnoise/saturate from the WF3
header). Mask written to DQ bit 4096 with `resetbits=0`.

**Per-lens LACosmic params: `info/lens_cr_params.json`.** The `--lacosmic-sigclip`/
`--lacosmic-objlim` defaults are `None` in all four drizzle scripts, resolved after parse to
this JSON (`{sample:{lens:{filt:{"lacosmic_sigclip":…,"lacosmic_objlim":…}}}}`) else the
hardcoded `4.5`/`5.0`; an explicit CLI flag always wins (precedence: default < JSON < CLI,
same as `info/psf_stars.json`). This **persists** a lens that needed non-default CR tuning so
a blind re-drizzle keeps it instead of reverting to the default. Currently one entry:
**J1420+6019 f814W = `sigclip 10 / objlim 12`**, because the default eroded a real all-frames
arc knot into a weight-0 noise hole (→ memory: todo_lacosmic_erodes_real_features; a
2026-08-18 sample-wide scan confirmed J1420 is the *only* affected product, so the pipeline
default is unchanged). A run whose lens/filt has an entry prints `=== LACosmic params for …
(info/lens_cr_params.json) ===`. Read-only; nothing writes it — edit by hand.

`driz_cr` compares each frame to a blotted median; on a steep, undersampled PSF core that
reference can read low, flagging real core pixels as CRs. **But this erosion is conditional on
dither quality, not intrinsic to driz_cr** (re-measured 2026-08-10, → memory:
drizcr_erosion_is_alignment_conditional): on well-dithered J0330-0020 under MAST alignment
driz_cr keeps the 1″ core to 0.998 of the no-CR pass, flat across `driz_cr_snr`/`scale`/combine
settings; the erosion (core→0.563, the old "~37% loss") only appears when the dither is forced
to zero (the pre-`--align mast` TweakReg default). **LACosmic stays the default** because it is
robust *regardless* of dither quality (per-frame, object-protected, no stacked reference), so
poorly-dithered/single-visit lenses get a clean core too.

**`resetbits=0` is mandatory on the LACosmic pass** — it defaults to 4096 and would clear
the DRIZ_CR bit the mask lives in, silently producing an un-masked drizzle that still looks
plausible (the first run reported a flawless `core=1.000` that was the no-CR image scored
against itself). Conversely the no-CR passes pin `resetbits=4096` explicitly, else they
inherit the mask and become silently CR-rejected.

**Do not add CR rejection to F160W.** At 0.1283″/px the IR PSF is ~1 px FWHM and looks like
an outlier: LACosmic zeroes field stars, `driz_cr` costs ~10% of a star's peak, at every
`objlim`. WFC3/IR FLTs are already up-the-ramp CR-rejected.

## Output pixel scales

Chosen by measurement (weight-map uniformity + FWHM + noise-map correlation), not
convention (see *Drizzle correlated noise* below).

| Band | Instrument | Native | Output | pixfrac |
|---|---|---|---|---|
| F606W | WFPC2/WF3 | 0.0996″ | 0.05″ | 1.0 |
| F814W, F555W | ACS/WFC | 0.05″ | 0.05″ (native) | default |
| F160W | WFC3/IR | 0.1283″ | 0.06″ | 1.0 |

- **F160W stays at 0.06″** — do not re-drizzle to 0.05″ to grid-match. PyAutoLens ingests
  each band at its native scale and fits sub-pixel inter-band offsets, so a common drizzle
  grid buys nothing; 0.05″ opens no empty pixels but worsens weight non-uniformity. **Kept
  at 0.06″ / pixfrac 1.0, settled by the user 2026-07-26 after a scale/pixfrac scan — do
  not reopen without being asked.**
- **F160W and F606W use pixfrac 1.0, not 0.8.** With ERR-based noise the goal is a uniform,
  low-correlation noise map for the likelihood, not the sharpest PSF (which PyAutoLens fits
  explicitly). pixfrac 1.0 drops adjacent-pixel noise correlation markedly at a small PSF
  cost. ACS F814W/F555W keep the default pixfrac (native 0.05″, no oversampling).
- **F555W needs no new script** — `drizzle_acs_wfc.py` is filter-agnostic. The 16 F555W
  lenses are exactly the 16 without WFPC2 F606W (props 10494/10798).

## Drizzle correlated noise (matters for strong-lens modelling)

Drizzling onto a finer-than-native grid correlates
adjacent output pixels (each output pixel is a weighted sum of overlapping input pixels;
neighbours share input pixels, so their noise is covariant) — intrinsic to drizzle
resampling, stronger with oversampling. `final_wht_type=ERR` makes the **per-pixel
variance** correct but says nothing about the **off-diagonal covariance**: on J1143 blank
sky the empirical pixel RMS was ~1.47× the ERR-map prediction (pixfrac 0.8; pixfrac 1.0
shrinks but doesn't remove it).

**Measured 2026-09-21, and the ratio above needs checking.** On blank sky (r > 4.5″, masked
pixels excluded) the empirical per-pixel rms comes out **BELOW** the noise map, by the amount
the drizzle correlation predicts — which is the opposite direction to the "1.47× the ERR-map
prediction" quoted above. Two independent estimators agree exactly (a 25 px local-median
subtraction, and a lag-8 px pair-difference that is immune to the correlation), and the deficit
tracks the **`CASR`** factor `make_cutouts.py` already stamps in each noise header:

| band / detector | measured 1/(emp÷map) | `CASR` in the header | adjacent-pixel corr |
|---|---|---|---|
| F606W WFPC2/PC | **2.45** | 2.39 | 0.58 |
| F555W ACS/WFC | 1.31 | 1.50 | 0.24 |
| F814W ACS/WFC | 1.35 | 1.50 | 0.24 |
| F606W WFC3/UVIS | 1.20 | 1.30 | 0.17 |
| F814W WFC3/UVIS | 1.12 | 1.30 | 0.16 |

So **the noise maps are behaving as designed** — `1/sqrt(WHT)` is the uncorrelated-equivalent
σ, and a blank-sky pixel rms reads low by ~R because drizzle smooths. Practical consequences:
(a) **do not "correct" a noise map because blank sky looks quiet** — that is the expected
signature, not an error; (b) the analytic `CASR` is ~10–13% high against the measurement for
ACS and UVIS and spot-on for WFPC2; (c) **WFPC2 F606W is by far the most correlated band**
(R≈2.4, neighbouring pixels 58% correlated), which is worth remembering before quoting a
per-pixel S/N on a WFPC2-only lens like J1403+0006. **Open item:** reconcile the 1.47×
sentence above with these numbers — it may have been measured against a different predictor,
but as written it points the correction the wrong way.

**Modelling implication:** a per-pixel independent-Gaussian likelihood (diagonal covariance)
does not bias the best-fit but mis-estimates *uncertainties/evidence*, worse in the
oversampled F606W/F160W than in native-scale F814W. The residual correlation factor is
band-dependent (~1.24 ACS, ~1.17 F160W, ~1.1 F606W, ~1.5–1.6 WFC3/UVIS gallery — higher
than native ACS despite also being native-scale, plausibly UVIS's larger geometric
distortion; see *BELLS GALLERY* below), so ignoring it also mis-weights bands relative to
each other in a joint fit. `make_cutouts.py --corr-factor` applies it (default 1.0, off).
Prefer native-scale F814W where a clean per-pixel noise model matters most.

## Bad-column fill (`--bcfill`) — the dead-column noise stripe

ACS and WFPC2 dead columns are DQ-flagged and dropped by AstroDrizzle, so on dithered frames
the affected output pixels get contributions from fewer exposures — a weight deficit that
shows as a **diagonal noise stripe** in `1/sqrt(WHT)` (diagonal because `final_rot=0` rotates
the detector-vertical columns by the exposure roll). The stripe is *correct* noise (those
pixels really do have fewer independent frames; measured amplitude ~1.1–1.4×, the √(N/N-eff)
coverage penalty — → memory: `legacy-slacs-bolton-bilinear-no-stripes`). `--bcfill` is an
**opt-in** reduction that removes it at the source: per input frame, linearly interpolate the
flagged columns across in SCI (and ERR/IVM) and clear the DQ bits *before* drizzling, so the
weight map comes out uniform by construction. Productionized 2026-08-13 from the
`scripts/bolton_investigations/redrizzle_bcfill*` prototypes (the investigation's chosen
option 3).

- **`drizzle_acs_wfc.py --bcfill`**: fills ACS bits **4|128** (bad detector pixel + bad
  column) in SCI+ERR. **`drizzle_wfpc2_wf3.py --bcfill`**: fills WF3 bits **2|256**
  (WFPC2 c1m calibration/mask defect + 256) in SCI, **interior-only** (spares the vignetted
  chip border), *before* `build_ivm_files` so the IVM noise model is rebuilt consistently on
  the filled columns. Per-visit like every WFPC2 product. Both stamp `BCFILL=True` in the
  product header.
- **A parallel tree at the DRIZZLE layer only, keyed via `cutout_paths.py`'s `variant` axis**
  (orthogonal to `--size`): `data/drizzled_bcfill/`, work dir `data/drizzle_files_bcfill/`.
  Both gitignored (big mosaics / scratch); the standard drizzle is never touched.
  **`make_cutouts.py --bcfill` cuts from that tree into the SAME `data/cutouts/` band dir
  under a TAGGED name** (`cutout_cr_bcfill_{sci,noise}.fits`, `cutout_cr_bcfill.png`;
  2026-09-22, replacing the 2026-09-14 header-only scheme). The header still carries
  **`BCFILL=T`** and `info/lens_cutout_qc.json` still records `bcfill`, but the name is what
  keeps the two reductions apart across branches. One stamp per band is enforced at the cut:
  a `--bcfill` cut into a dir holding a standard stamp (or vice versa) is **refused unless
  `--force`**, which then removes the other reduction's sci/noise/png. Every reader goes
  through `cutout_paths.find_stamp`, so `make_mosaics.py`, `make_masks.py` etc. need no
  `--bcfill` flag and work unchanged on either branch.
- **NOT for WFC3/IR F160W** — an IR array has **no bad columns** (verified on J0822: 0
  columns >50% dropped vs 57 for WF3); its noise-map dots are quadrupled hot-pixel replicas,
  a different, already-correct artifact AGENTS.md warns against altering. `--bcfill` exists
  only on the ACS and WFPC2 drizzle scripts.
- **Caveat (document it downstream):** the filled pixels carry no independent information, so
  the noise is **optimistic by ~√(3/4)≈13%** on those columns — a cosmetic/uniformity choice,
  not the science default. To keep an honest noise map with a clean *image*, inflate those
  pixels by the per-pixel √(WHT_filled/WHT_baseline) correction (≈ recovers the standard
  drizzle stripe) rather than masking to a huge value.
- Validated: ACS J1023+4230 (6.3% of the cutout is stripe, ratio ≤1.23) and WFPC2
  J0252+0039/J0822+2652 (4.7% of cutout, ratio median 1.12 / max 1.41; fill count
  detector-fixed at ~2020 px/frame, single- and split-visit identical). The standalone
  comparison figures (`redrizzle[_wfpc2]_*bcfill_compare.png`) were deleted with
  `diagnostics/bolton_test_outputs/` on 2026-09-13; rerun the scripts, or read them out of git
  history, to see them again.

**DECIDED 2026-09-22: `--bcfill` is NOT the modelling input. `main` carries the standard
stamps; the bcfill products live on the `bcfill` branch.** The reasoning, checked against
source, not memory:

- The fit's likelihood (`PyAutoArray/autoarray/fit/fit_util.py`) is
  `log L = -0.5 * [ sum((data-model)^2/sigma^2) + sum(log(2*pi*sigma^2)) ]`. Each pixel's
  residual is weighted by 1/sigma^2; the normalization term depends on sigma only, so it
  shifts absolute log L between reductions but leaves delta-log-L between ladder rungs on
  one dataset unchanged. Nothing in the fit *needs* a uniform noise map.
- The standard (non-bcfill) drizzle already has a stripe-free **science** image — the
  dithered good frames fill every affected pixel — and an **honest** noise stripe (1.1-1.4x,
  the sqrt(N/N_eff) coverage penalty). bcfill barely changes the science (median |delta|
  ~1.8e-6 e/s on J1023) and only removes the stripe from the noise map, making those pixels
  ~13% optimistic, i.e. over-weighted by ~1.3x in chi^2, on ~5-6% of a cutout, sometimes
  through the ring/core. That is a bias in the wrong direction for a fit; the stripe "looking
  worse" is the map being correct.
- Conclusion: **the standard drizzle stamp is the more correct modelling input; bcfill is
  cosmetic and not needed downstream.** Before the revert `info/lens_cutout_qc.json`
  recorded bcfill=True on f814W 42/48, f606W 47/62 (+1 v2), f555W 16/16; F160W and UV bands
  were always standard.
- **What was done.** (1) Stamp names now carry the reduction (see the top-of-file note and
  `cutout_paths.py`), so the two branches cannot silently swap stamps. (2) The `bcfill`
  branch was cut from `main` with the 105 bcfill stamps renamed to their tagged names and
  pushed; it is the collaborators' copy. (3) On `main` the same 105 bands were re-cut from
  `data/drizzled/` (standard, `--force`), mosaics and dataset subplots rebuilt. Geometry is
  identical (same drizzle call, same grid) so every mask and positions file carried over —
  verified by comparing `CRPIX`/`CRVAL` of each new stamp against the bcfill one, not
  assumed. (4) The two `--crfill` bands (J1213+6708, J1250+0523 f814W) were `bcfill_crfill`;
  see *`--crfill`* below for what replaced them.
- **Moving mask work between branches:** `git cherry-pick` the mask commit, never `git
  merge` — a merge would also carry the other branch's stamps across (they would then sit
  beside this branch's under a different name, and every reader stops with
  `AmbiguousStampError` until one is removed). A cherry-picked commit that also rebuilt the
  QC PNGs renders the other reduction's stamp underneath; regenerate them on the branch.
- Alternative that keeps bcfill's image with an honest noise map, if ever wanted: inflate
  the noise by the per-pixel sqrt(WHT_filled/WHT_baseline) correction above, which
  reproduces the standard stripe anyway. The modelling repo's data-prep docs should say the
  fits read the standard (`main`) stamps.
- Not changed by this: the cosmic-ray reasoning below (same physics, same direction).

**DO NOT `--bcfill` A COSMIC RAY — the two defects are opposites (asked and measured
2026-09-21, J1213+6708 and J1250+0523 f814W).** A dead column is flagged at the *same detector
position in every frame*, so no dither recovers it and interpolation is the only way to make
the weight uniform; a cosmic ray hits *one frame at a random position*, and the other frames
cover that sky perfectly well. Filling a CR track would therefore add **no information at
all** — it would only make the weight map claim 4 frames where 3 exist, leaving the noise
optimistic by √(4/3) ≈ 15% along the track. Both tracks here run through the Einstein ring,
which is exactly where an optimistic noise map biases the likelihood, so the trade is worse
than the cosmetic one bcfill already documents. What LACosmic left behind was measured and is
small and *honest*:
- **J1213+6708 f814W**: one frame of four (`j9op28baq`) carries a ~12″ track at image PA 74.5°
  passing **0.68″ from the deflector**. Residual noise ridge ~3 px wide: **+11% inside
  r < 1.8″** (noise/local 1.152 on-track vs 1.039 off), ~+1% median further out.
- **J1250+0523 f814W**: a fatter, shorter bar, **2.1″ × 0.6″ at r = 0.87–1.71″** (i.e. across
  the ring, screen up-right), **+9%** noise (1.112 vs 1.019), weight loss 19% median / 45%
  peak.
- **Both were cleanly REMOVED from the science image** — there is nothing to mask. On-track
  radial-subtracted flux matches off-track (J1250 +0.69σ vs +0.62σ control; J1213 has *fewer*
  >3σ pixels on-track than off, 4.3% vs 5.4%), and **neither feature exists in the lens's
  other band** (J1213 f606W and J1250 f555W read on-track = off-track in every statistic),
  which is the cheap band test that proves a per-exposure event rather than sky.
- **The discriminator to reuse**: with `final_wht_type=ERR` a CR pixel already carries a huge
  ERR and so contributes almost no weight, which is why rejecting it costs only ~10% and why
  `wht_cr/wht_nocr` under-reads the defect. Compare **noise / local-median** on the two bands
  of the same lens — source structure appears in both, a CR only in one. Figures:
  `diagnostics/cr_defects/J1213+6708_f814W_cr_streak.png`, `J1250+0523_defect_zoom.png`,
  `*_single_frames.png` (the per-frame drizzles, where the offending exposure is obvious).
- Caveat on the per-frame `*_single_sci.fits` in `data/drizzle_files/`: they are **not
  mutually registered** in this repo's work dirs (the galaxy lands tens of px apart between
  frames), so use them to identify *which exposure* carries a CR, never to measure *where*.

### `--crfill` — the CR-fill, wired up 2026-09-21 so the trade can be measured

`drizzle_acs_wfc.py --crfill` and `make_cutouts.py --crfill`. **Same mechanic as `--bcfill`,
different DQ bit and a different trade.** After `run_lacosmic` writes DQ 4096 and before
anything drizzles, `fill_cosmic_rays()` interpolates SCI *and* ERR across each flagged pixel
and clears the bit, so the drizzle weights it like any other and the track leaves no
weight/noise residue. Products go to `data/drizzled_crfill/`, or
`data/drizzled_bcfill_crfill/` with both (the `_variant` tag composes), and carry
**`CRFILL=True`** in the header exactly as bcfill carries `BCFILL`; the cut stamp is tagged
the same way (`cutout_cr_crfill_sci.fits`, `cutout_cr_bcfill_crfill_sci.fits`).

**STREAK-REPAIRED PRODUCTS ARE NOW THE STANDARD ON TWO BANDS (2026-09-22, user's call), and
`info/lens_crfill.json` is the list.** Read that file before using either stamp; it records the
exact command, the tracks kept (frame, chip, area, closest approach), pixels filled, the ridge
before and after, and what moved downstream.

| sample / lens / band | tracks | px filled | ridge | astrometry moved? |
|---|---|---|---|---|
| `slacs_gold` **J1213+6708 f814W** | 1 of 220,590 | **210** of 941,414 | **+11.0% → +1.6%** | **YES — 180 mas** |
| `slacs_gold` **J1250+0523 f814W** | 4 of 258,884 | **256** of 1,113,026 | **+9.2% → +1.9%** | no (1 mas) |

- **What you are getting.** The cosmic-ray weight/noise ridge across the science region is
  gone, and **the filled pixels' noise is optimistic by √(4/3) = 1.155** — 210 px on J1213,
  256 px on J1250, and nothing else. That confinement is measured, not assumed: a raw
  subtraction against a zero-fill null on an identical grid changed **392 px of 57,600**, all
  on the track. Nothing else in either stamp moved, and neither arc is touched (net flux
  +0.1%).
- **Settings differ per lens, on evidence.** J1213 needed one track (a single track already
  matched a frame-wide fill there, +2.4% vs +2.3%); J1250 needed four, because three smaller
  tracks also sit on its ring and one track only reached +2.8% against +1.4% frame-wide.
  **Do not copy one lens's flags to the other** — count the streaks in the noise map first.
- **J1213+6708 also changed ASTROMETRY, and that is the bigger deal.** It is the `--align
  tweakreg` lens, so the re-drizzle lands **180 mas (3.6 px)** from the archived product; the
  measured content shift was +1.47 rows, +0.01 cols. Its `cutout_cr_mask.fits` and
  `cutout_cr_mask_arcs.fits` were **shifted by the integer (+1, 0) and NOT redrawn**, leaving
  a **0.47 px residual** (`MASKSHFT` / `MASKSHRS` header cards). That lens has no positions
  file, so nothing else needed re-deriving — but any *future* f814W-proposed mask on it starts
  from a mask that is half a pixel off.
- **J1250+0523 cost nothing**: `CRPIX` is identical to the archived stamp, so its contaminant
  mask, arc mask and `cutout_cr_positions.json` all stay valid unchanged.
- Rebuilt on both: `cutout_cr_dataset.png`, `cutout_cr_mask_arcs.png`, and
  `data/mosaics/slacs_gold/`. Figure: `diagnostics/cr_defects/as_is_vs_streak_repaired_snr.png`.
- **The superseded products are the plain bcfill stamps from the 2026-08-13 campaign.** They
  are recoverable by re-running the same drizzle without `--crfill` — except on J1213+6708,
  where a re-run does **not** reproduce the archived astrometry (see the re-drizzle trap).

**`--crfill` is otherwise wired up, not adopted.** There is no campaign behind it, nothing downstream cuts
`--crfill` by default, and the argument against it is directly above: a CR hits one frame at
a random position, so filling it adds no information and only makes the weight map claim N
frames where N−1 exist. **Use it to compare against the other two options, which are
"drop the affected exposure" (honest noise, real exposure lost) and "keep as-is" (a +9–11%
noise ridge that is correct).**

- **Guarded**: `--crfill` without `--cr --cr-method lacosmic` is a hard error, because
  `driz_cr` writes 4096 *during* the drizzle (too late to fill) and `--no-cr` never writes it
  — otherwise you would get a `crfill` tree that is just the standard reduction renamed.
- **The interpolation is NOT bcfill's.** A dead column is narrow in every *row* it crosses, so
  bcfill interpolates row-wise, always. A cosmic ray is an arbitrary blob at an arbitrary
  angle: a near-horizontal track is narrow in every *column* and enormous in every row, and
  filling it row-wise would smear a real gradient across tens of pixels. So **each pixel is
  filled along whichever axis its own contiguous gap is shorter**, and the run prints the
  **widest gap actually crossed** — if that is large the fill is guessing, and the frame is a
  candidate for dropping instead. Unit-tested on synthetic vertical/horizontal/diagonal tracks:
  exact on a linear ramp, correct axis chosen in each case.
- **A pixel with no usable anchors on either axis is LEFT FLAGGED**, not fabricated — it keeps
  the honest weight deficit — and the count is printed.
- **`make_cutouts.py` guards both cards in both directions now.** A `crfill` stamp differs
  from the standard one *only in the noise map, along one track*, so overwriting either way is
  invisible in the science image. Cutting without `--crfill` over a `CRFILL=T` stamp (or
  without `--bcfill` over `BCFILL=T`) is refused unless `--force`. `info/lens_cutout_qc.json`
  gains a `crfill` key beside `bcfill`.

**MEASURED on J1213+6708 f814W (2026-09-21), the three options side by side.** Ran
`--bcfill --crfill`, cut a stamp from it and compared against the tracked stamp in each one's
own geometry (the re-drizzle lands on a slightly different grid — see the trap below), on the
same physical line and the same lens centre:

**All four options built and measured, both lenses (2026-09-21).** Figure:
`diagnostics/cr_defects/cr_options_snr_fourway.png` — S/N side by side, common stretch per row.

| | NDRIZIM | EXPTIME | defect ridge | **median S/N inside 2″** | stamp median noise | noise map true? |
|---|---|---|---|---|---|---|
| **J1213+6708 f814W** | | | | | | |
| as-is | 8 | 2240 s | **+11.0%** | **32.95** | 0.009360 | **yes** |
| frame drop (`j9op28baq`) | 6 | 1680 s | +1.3% | **28.83 (−12.5%)** | 0.010828 (+15.7%) | yes |
| crfill frame-wide | 8 | 2240 s | +2.3% | 33.32 | 0.009286 (−0.8%) | no, **everywhere** |
| crfill streak-targeted | 8 | 2240 s | +2.4% | 33.20 | **0.009351 (unchanged)** | no, on 210 px |
| **J1250+0523 f814W** | | | | | | |
| as-is | 8 | 2232 s | **+9.2%** | **21.70** | 0.008934 | **yes** |
| frame drop (`j9c713d0q`) | 6 | 1674 s | +1.2% | **19.18 (−11.6%)** | 0.010226 (+14.5%) | yes |
| crfill frame-wide | 8 | 2232 s | +1.4% | 22.33 | 0.008858 (−0.9%) | no, **everywhere** |
| crfill streak-targeted | 8 | 2232 s | +2.8% | 22.26 | **0.008934 (unchanged)** | no, on 194 px |

**Does `--crfill` change anything outside the track? NO — proved by direct subtraction
(2026-09-21).** The honest test is to subtract one product from the other and check the
difference is confined to the fill. Done on J1213+6708 f814W between two **fresh** runs of the
crfill path that differ *only* in what was filled — a zero-fill null (`--crfill-radius 0.001`,
0 px) and the streak-targeted run (210 px). They share a grid exactly (CRPIX identical), so
this is a raw subtraction with **no resampling**:

| | on the CR track | off the track |
|---|---|---|
| pixels whose SCI changed at all | 212 | 180 |
| max SCI change | 3.19σ | 3.30σ |
| max noise change | 14.1% | 14.8% |

**392 px of 57,600 changed, and they form ONE connected region sitting on the track** (331 px
at r = 0.88″, >0.5% noise change; the "off track" 180 px are the fill spilling just outside a
3 px-wide definition of the track, not a second region). The deflector moves **1 mas**. The arc
is untouched: **net flux +0.1%**, peak S/N identical.

**The earlier "the arc looks dimmer under crfill" was MY comparison artefact, now retracted.**
It came from comparing fresh products against the *archived* stamp, which on this lens sits
180 mas away (see the re-drizzle trap above) — through nearest-neighbour resampling, which
smears a compact feature. Against a fresh reference the four options read:

| variant | arc net flux | arc peak S/N | deflector moves |
|---|---|---|---|
| **fresh reference (0 px filled)** | +0.0% | 32.3 | — |
| crfill streak-targeted | **+0.1%** | 32.3 | 1 mas |
| crfill frame-wide | **+0.2%** | 32.3 | 2 mas |
| frame drop | **−4.8%** | 21.7 | 44 mas |
| *archived stamp* | *+43.2%* | *106.8* | *180 mas* |

Only the frame drop touches the arc, and the archived row is the outlier for the astrometric
reason above, not a reduction reason. Figure:
`diagnostics/cr_defects/cr_options_snr_fourway.png` (bottom-right panels are the raw
subtraction: blank everywhere except the track).

**Does any of this touch the ARC? No — checked on J1213+6708, and the arc is real.** Asked
because the arc looks very faint against the deflector in every variant. It does, and that is
astrophysics: on the **elliptical-isophote** residual (b/a 0.91, PA 111°) the arc is **366 px at
r = 1.85″ = 1.30 θ_E, peak 13.0σ**, tangential (84° from radial, b/a 0.49), while the deflector
peaks in the hundreds — a linear S/N stretch set by the galaxy cannot show both. Three things
make it a detection rather than a residual, and the third is the one that matters:
1. it clears the model-error floor by **4.6×** in F814W (floor 2.82× the photon noise) and
   **3.3×** in F606W (floor 1.84×);
2. it is tangential in **both** bands (84°, 86°) at the same radius (1.30, 1.27 θ_E) and the
   same sky position to 0.15″ — and F606W is WFPC2, sharing **no exposures** with ACS;
3. **it is NOT mirror-symmetric** — 13.0σ against **+1.0σ** at its mirror. A galaxy-model
   residual (dipole or quadrupole) *is* symmetric under 180° rotation, so cross-band agreement
   alone would prove nothing (both bands share the same galaxy and the same model error);
   the asymmetry is what rules the artefact out.
   **There is also a counter-image candidate**: 24 px at **0.68 θ_E**, 8.9σ, tangential 80°,
   on the opposite side (image PA 305° against the arc's 147°) — a double, marginally above
   the floor, worth confirming with a real galaxy model.
**Do not use the circular radial median on this lens**: it manufactures a quadrupole whose
lobes sit near both features. Figure:
`diagnostics/cr_defects/J1213+6708_f814W_arc_vs_cr_options.png`.
**Across the four options the arc survives every one; only the frame drop dims it** — arc
integrated S/N **106σ → 94σ (−12%)**, exactly the exposure it threw away. Caveat on comparing
against the as-is stamp here: this is the `--align tweakreg` lens, and a fresh drizzle recentres
differently (0.354″ vs 0.479″), which is worth ~±10% on an integrated number. Compare the fresh
products with each other.

**What the two right-hand columns settle.** The **frame drop is the only option that costs real
signal** — ~12% of the median S/N inside 2″, over the whole stamp, to remove a ridge on 0.3–0.4%
of it. Both CR-fills keep the S/N. The difference between them is the *stamp median noise*: the
frame-wide fill pulls it down ~0.9% across the whole stamp (that number **is** the collateral —
a global lie in the noise map), while the **streak-targeted fill leaves it bit-for-bit at the
as-is value** and confines the optimism to the ~200 px it actually filled. Streak-targeting is
therefore the only version of `--crfill` worth running.

(J1250+0523's defect is the 2.1″×0.6″ bar at r = 0.87–1.71″, 247 px = 0.43% of the stamp;
J1213's is a ~12″ track passing 0.68″ from the deflector, 171 px inside 1.8″.)

Read the table this way. **`--crfill` works** — it removes 75–85% of the ridge on both lenses.
But it converts a *correct* +10% ridge on 0.3–0.4% of the stamp into an *incorrect* −13% one:
after filling, the true noise on those pixels is still √(4/3) = 1.155× what the map says. And
dropping the exposure costs **+15.5% over the whole stamp** to remove ~+10% on 0.4% of it — the
worst of the three on both lenses. **So "keep as-is" remains the recommendation**, now on
measured grounds rather than argued ones; `--crfill` exists for the case where a track is so
bad it dominates a fit, and that case has not appeared yet.
Figure: `diagnostics/cr_defects/crfill_before_after.png`.

**`--crfill-radius ARCSEC` makes the fill TARGETED (2026-09-21, at the user's suggestion), and
it is what you want.** Plain `--crfill` fills *every* LACosmic flag in every frame — ~1.4% of
all pixels — so the noise map goes optimistic wherever any cosmic ray landed, not just on the
track crossing the science region. With a radius, only the tracks that reach within that
distance of the lens are filled. **Measured on J1213+6708 f814W, `--crfill-radius 4`:**

| | frame-wide `--crfill` | `--crfill-radius 4` |
|---|---|---|
| pixels filled | **941,414** (100%) | **989** (0.1%) |
| CR tracks filled | 220,590 | 193 |
| track ridge, r<1.8″ | +10.4% → **+2.7%** | +10.4% → **+2.7%** |
| ridge at 1.8–6″ | +1.1% → +0.2% | +1.1% → +0.2% |
| CR speckle left **outside** 4″ (% px with noise/local > 1.05) | 7.1% → **4.1%** (erased) | 7.1% → **8.4%** (untouched) |

**Identical repair, 1/950th of the collateral.** Figure:
`diagnostics/cr_defects/crfill_targeted_vs_framewide.png`. The per-frame breakdown also names
the culprit for free — `j9op28baq` contributed 411 px against 180/190/208 for the other three.
**Selection is by connected component, never clipped at the circle**: a track that reaches into
the region is filled along its whole length, because a track half filled and half flagged would
put a step in the weight map partway along it — a worse artifact than the ridge being removed.
`--crfill-radius` without `--crfill` is an error, and it needs the lens in
`info/slacs_coords.py`.

**`--crfill-min-area PX` targets the STREAK instead of an aperture (2026-09-21).** A radius
alone still fills every incidental CR hit that happens to share the aperture — 193 tracks on
J1213+6708, 283 on J1250+0523, almost all a few pixels each — when the thing worth repairing is
the one long track you can see in the noise map. Pixel count separates them cleanly: ordinary
ACS hits are <20 px, a long track is hundreds. **`--crfill-radius 5 --crfill-min-area 50` on
J1250+0523 f814W selected exactly ONE track of 258,884** — 194 px in `j9c713d0q` chip 2,
closest approach 1.02″ — i.e. it found the bar, named the exposure that carries it, and filled
**194 px instead of 1,113,026**. The run prints every track it keeps (area, chip, closest
approach) so the selection is auditable rather than trusted.

| targeting | J1213+6708 px filled | J1250+0523 px filled | tracks kept |
|---|---|---|---|
| none (frame-wide) | 941,414 | 1,113,026 | all (220,590 / 258,884) |
| `--crfill-radius 4` | 989 | 1,377 | 193 / 283 |
| `--crfill-radius 5 --crfill-min-area 50` | **210** | **194** | **1 each** |

On both lenses that single track is in a single exposure — `j9op28baq` (J1213, closest approach
1.17″) and `j9c713d0q` (J1250, 1.02″) — so the selection also **identifies the exposure to drop**
if you would rather take that option. `--crfill-min-area` alone (no radius) is legitimate too:
it fills every long track in the frame and leaves the pinprick hits alone.

**`--crfill-n-tracks N` — "fill the N biggest streaks" (2026-09-21), and J1250+0523 needs it.**
A size threshold is awkward to pick per lens; a count is not. `N` is ranked **globally across
every frame and chip**, not N-per-frame, because the streaks you can count in a stamp generally
sit in *different exposures*. It needs a first pass over the frames to rank them, which is why
it is its own flag.

**It found what the eye found on J1250+0523: FOUR cosmic-ray tracks on the ring, not one.**
`--crfill-radius 2 --crfill-n-tracks 4` ranked 56 candidates and kept four, from **two**
exposures — `j9c713d0q` chip 2 (the 194 px bar, 1.02″) and `j9c713cuq` chip 2 (31, 17 and 14 px
at 0.72, 0.76 and 0.89″). **256 px filled of 1,113,026.** The three small ones are why the
single-track fill only got the ridge to +2.8%:

| J1250+0523 f814W | px filled | tracks | bar ridge |
|---|---|---|---|
| as-is | 0 | — | **+9.2%** |
| `--crfill-min-area 50` | 194 | 1 | +2.8% |
| `--crfill-radius 2 --crfill-n-tracks 4` | **256** | **4** | **+1.9%** |
| frame-wide | 1,113,026 | all | +1.4% |

Four tracks at 256 px get within 0.5% of what a million-pixel blanket fill achieves. Figure:
`diagnostics/cr_defects/J1250+0523_f814W_crfill_n_tracks.png`.

**Did the ranking pick the ones a human would? On this lens, yes — checked, don't assume.**
The user named four from `diagnostics/cr_defects/J1250+0523_f814W_cr_inventory.png`: the bar,
the two inside θ_E, and the extended feature by #5. Verified against the *stamp's* noise map
(the fraction by which each pixel's noise dropped), which is what the eye actually judges:

| inventory # | stamp position | r | noise drop | repaired? |
|---|---|---|---|---|
| #1 (the bar) | col 135, row 100 | 1.25″ | 13.8% | **yes** |
| #7 | col 126, row 105 | 0.79″ | 12.2% | **yes** |
| #6 | col 137, row 113 | 0.93″ | 14.9% | **yes** |
| #5 | col 119, row 142 | 1.13″ | 6.4% | **yes** |

All four. But **size-rank is a proxy for "what you can see", not the same thing**, and the two
can part company: the ranking sorts *native-frame* track pixel counts in individual exposures,
while prominence in the stamp depends on where the track lands (a small track on the bright
ring outranks a bigger one in blank sky) and on how the four exposures combine. Here #6 was
not even a separate repaired region — it merged with the bar's footprint in the stamp. **So
treat the count as a convenience, not an oracle**: the run prints every track it keeps (area,
chip, closest approach) precisely so you can check the list against the inventory figure before
accepting the product. If a lens ever needs a specific set that ranking will not produce, the
right primitive is naming positions, not a bigger N — that has not been built, because it has
not been needed.

**`--exclude-frames ROOTNAME[,...]` builds the frame-drop product (2026-09-21)**, so the third
option is measurable rather than hypothetical. The named exposures never reach the work
directory, so bestrefs, alignment, LACosmic, both drizzle passes and the provenance JSONs all
see a genuinely short visit — which is what a real frame drop is. Products land in a `_drop`
variant tree (`data/drizzled_bcfill_drop/…`) and carry **`DROPFRMS`** in the header;
`make_cutouts.py --drop` cuts from it. J1213+6708 f814W went **NDRIZIM 8 → 6, EXPTIME 2240 →
1680 s**, which is the cost in one line.

**A region-targeted FRAME DROP, though, is not a third option — it is what we already have.**
"Drop that exposure only where the cosmic ray is" *is* flagging that frame's CR pixels and
letting the drizzle weight them out, which is exactly what LACosmic + `final_bits` already do;
the +10.4% ridge **is** its cost, already paid and already honest. The only knob a "targeted
drop" could add is dropping a *larger* region of that frame, which is strictly worse: excluding
one of four frames everywhere inside 4″ would put **+15.5% on every pixel inside 4″** against
today's +10.4% on the 171 px of the track. So the real menu is two items, not three —
keep the honest ridge, or fill it with `--crfill --crfill-radius` and accept a noise map that
is √(N/(N−1)) optimistic on the track.

**TRAP found doing this: `make_cutouts.py` writes `info/lens_cutout_qc.json` even with
`--output` pointing somewhere else.** A scratch cut therefore overwrites the tracked stamp's
provenance (here `weight_uniformity` 0.315083 → 0.314273, `offset_arcsec` 0.4794 → 0.3537, and
a spurious `crfill: true`) while the tracked stamp itself is untouched — so the JSON silently
starts describing a file that is not in the tree. Restored by hand. **Either pass `--output`
only when you are willing to fix the JSON afterwards, or fix the script to skip the
`info_json.update` when `--output` is set.**

**Second trap, and it is worse than "not bit-identical": A RE-DRIZZLE OF J1213+6708 f814W
LANDS 180 mas (3.6 px) AWAY FROM THE ARCHIVED PRODUCT.** Measured 2026-09-21 by centroiding the
deflector in each stamp and comparing *sky* positions, not pixels. This is the one lens with
`--align tweakreg` (`ALIGN_OVERRIDES`), so TweakReg re-solving the alignment is the obvious
suspect; the drizzle grid shifts sub-pixel with it, so the two stamps are not even related by
an integer shift. Consequences, all of which bit me before I caught it:
- **Any comparison of a fresh product against an archived stamp on a tweakreg lens is
  confounded at the 3–4 px level**, and nearest-neighbour resampling onto the archived grid
  then makes compact features (an arc, say) look systematically fainter. That is a comparison
  artefact, not a reduction difference.
- **A fresh reference is mandatory.** Re-run the pipeline with the *same* flags plus a no-op
  (`--crfill --crfill-radius 0.001` fills 0 px and is a perfect null), cut that, and compare
  everything against it. Two runs of the crfill path then share a grid **exactly** (CRPIX
  identical to 0.000000) and can be differenced with no resampling at all.
- **`--sample`-wide re-drizzles will not reproduce archived stamps on this lens**, so do not
  treat a diff against it as a regression.
- Beware the skip: `drizzle_acs_wfc.py` prints *"drizzled products already exist, skipping"*
  and exits 0. A "control re-run" that silently skipped will look perfectly reproducible
  because it is the same file. **Check the log says it actually drizzled.**

**`run_acs_all.sh` and `run_wfpc2_wf3.sh` are `--bcfill`-aware** — pass `--bcfill` as a
second arg (after the optional sample) and it threads through the *drizzle* stages into the
parallel `data/drizzled_bcfill/` tree, with `_bcfill`-tagged logs. `run_acs_all.sh --bcfill`
*drizzles* f814W+f555W (no cutout step — it never cuts); `run_wfpc2_wf3.sh --bcfill` runs the
full drizzle → tie → cutout for F606W (the tie is against the bcfill tree's OWN f814W, so the
ACS bcfill drizzle must run first). A whole-sample campaign is thus: `run_acs_all.sh <sample>
--bcfill` → `run_wfpc2_wf3.sh <sample> --bcfill` → `make_cutouts.py --bcfill` per ACS product
(the cutout step the ACS runner skips) → `make_mosaics.py --sample <sample>`. **The cutout and
mosaic steps land in the ordinary `data/cutouts/`/`data/mosaics/`** — bcfill supersedes the
standard stamp in place — so `make_mosaics.py` takes no `--bcfill` and the bcfill `--size`
tag composes with nothing. The other runners
(`run_wfc3_all.sh`, `run_gallery_uvis_all.sh`, `run_cutouts_all.sh`, `run_psf_all.sh`) are
**not** `--bcfill`-aware — bcfill is ACS+WFPC2 only, so there's nothing for them to pass.

**Campaign state (2026-08-13, 0 failures):** ran the full recipe across `slacs_gold` (f814W
38, f555W 16, f606W 22 incl. split J0822+2652) and `slacs_other` (f814W 4 — the only ones past
`BLOCK_EXPTIME`; f606W 24), into the `*_bcfill` trees (`BCFILL=True`; F606W tied,
`GSC240FX=True`). **`gallery` has no bcfill** — entirely WFC3/UVIS, and `drizzle_wfc3_uvis.py`
has no `--bcfill` (porting the fill to UVIS would be new dev + validation, not a campaign).

**Cross-tree astrometry trap the bcfill F606W tie can't self-heal (J1016+3859, slacs_other).**
`run_wfpc2_wf3.sh --bcfill` ties F606W to the bcfill tree's OWN F814W (a raw re-drizzle with
the delivered WCS), so if that ACS WCS is itself off (`acs_wcs_not_always_gaia`: J1016+3859's
F814W came down `-GSC240`, ~0.6–0.9″ off), bcfill F606W ties to a wrong frame and recentres
~0.7″/14px off — while the standard tree's F814W→F160W fix never propagated (there's no bcfill
F160W; IR has no bad columns). Fix (the script's `--ref` can't cross variants): tie **bcfill
F814W → standard-tree F160W** manually (centroid shift via `align_wfpc2_to_acs`'s
`stable_centroid`/`find_product`, apply dRA/dDec to bcfill F814W CRVAL, stamp
`ASTROREF='f160W'`), then re-run `align_wfpc2_to_acs.py --bcfill` and re-cut both bands. Done
2026-08-13 (0.194″/3.9px). Watch for this on any bcfill lens whose standard F814W has a
non-default `ASTROREF`.

## Cutouts (`scripts/make_cutouts.py`)

```bash
uv run python scripts/make_cutouts.py --lens J0029-0055 --filt f606W
```

Cuts a square stamp (default 12″) from `data/drizzled/` into `data/cutouts/`: a sci FITS, a
noise FITS (from the weight map), and a 3-panel PNG.

- **Stamp size is a parallel product tree, not an overwrite (`--size`, `scripts/cutout_paths.py`).**
  The cutout filenames carry no size, so a re-cut at a different size would be an invisible
  clobber. `--size` therefore *derives* the output tree: 12″ keeps `data/cutouts/` +
  `data/mosaics/` + `info/lens_cutout_qc.json`; any other size S writes
  `data/cutouts_<S>arcsec/`, `data/mosaics_<S>arcsec/` and `info/lens_cutout_qc_<S>arcsec.json`.
  Both scripts take it (`make_cutouts.py --size`, `make_mosaics.py --size`) and
  `run_cutouts_all.sh [SAMPLE] [SIZE]` passes it through; every stamp carries `CUTSIZE` in its
  header. **A 20″ set exists for all three samples alongside the 12″ one** (2026-08-04; 90 +
  34 + 33 products, 0 failures). PSF kernels are *not* duplicated per size — the kernel is
  trimmed by amplitude, so it's a property of the band, and a size-variant stamp pairs with
  the same `cutout_[cr_]psf.fits` from the default tree — which is the whole of what
  `cutout_paths.psf_cutout_dir()` now asserts (before 2026-09-14 it also had to pick between
  the bcfill and standard trees). See *PSF generation*.
  **`data/cutouts/` (12″) is the tree
  tracked in git**, made the pipeline default 2026-08-11 (swapped in place from the old 20″
  default, incl. the 422 `cutout_[cr_]psf*.fits`/mosaic-panel files `make_psf.py` and
  `make_psf_mosaics.py` then hardcoded into `data/cutouts/`/`data/mosaics/` rather than going
  through `cutout_paths.py` — the PSF half of that now goes through `psf_cutout_dir()`)
  because it's what `scripts/make_masks.py` targets and what
  downstream modelling reads. Other size variants (now including `data/cutouts_20arcsec/`)
  and their QC JSONs are gitignored (regenerable in one runner call, and the tree is ~550 MB,
  mostly PNGs whose size follows the figure, not the stamp) — **unless a tree carries
  non-regenerable content**, which is why `data/cutouts/` is the exception: `make_masks.py`
  writes its hand-drawn GUI masks there, and unlike every other file in the tree those aren't
  regenerable by any script, so losing them means redoing manual work — see *Masks* below.
  `data/mosaics_20arcsec/` stays gitignored; nothing non-regenerable lives there.

- **Pass + prefix.** `--pass {auto,cr,nocrrej}` (default `auto`) picks the CR pass when one
  exists, else no-CR. The prefix encodes it so the two coexist: **`cutout_cr_*` for CR,
  `cutout_*` for no-CR.** ACS F814W/F555W and WFPC2 F606W now default to a LACosmic CR pass
  → `auto` cuts `cutout_cr_*`, which is science-grade (LACosmic preserves the core). WFC3/IR
  F160W has no CR pass, so `auto` falls back to `cutout_*` there.
- **Recentring.** The stamp recentres on the galaxy, and the peak search prefers the CR
  mosaic — a brightest-pixel search on a CR-laden no-CR mosaic locks onto cosmic rays. Don't
  reach for `--median-size` for a bad recentre; check a `*_cr_*` mosaic is present.
  Offsets of 1–2″ aren't always failures — some lenses (J0912+0029, J0956+5100) have genuine
  multi-knot morphology.
- **Shared centre (`--center-band`, default `f814W`).** All bands are cut about a single
  centre from the `--center-band` mosaic (highest S/N, GAIA-accurate) so stamps co-register
  across filters. `--center-self` restores per-band recentring; a band whose center-band
  products are missing falls back to its own peak with a warning. **Never use
  `--center-self` on gallery's UV bands** (F225W, F275W) — a self-centred peak search locks
  onto noise or an unrelated field source, not the lens. These bands are in fact unusable
  for lens science sample-wide (confirmed, not just hard to centre) — see *BELLS GALLERY*
  below.
- **QC diagnostics + provenance (`info/lens_cutout_qc.json`, 2026-07-31).** Every run
  computes and prints two diagnostics (ported from a comparison against the sibling
  `PyAutoReduce` framework), purely informational — neither changes any output:
  - **Weight-map uniformity** (`weight_uniformity()`): the STScI RMS/median rule-of-thumb
    on the **cutout's** own weight-map region (not the full mosaic, which mixes coverage
    tiers across the whole union footprint). `<=0.2` is the accepted threshold; stamped as
    `WHTUNIF`/`WHTUNIFL` in `cutout_[cr_]noise.fits`'s header.
  - **Analytic Casertano R** (`casertano_r(pixfrac, scale_ratio)`): a closed-form
    correlated-noise factor from `D001PIXF`/`D001ISCL`/`D001SCAL` alone (verified against
    the textbook values: R=1.5 at p=1,s=1; 1.364 at p=0.8,s=1; 1.25 at p=0.6,s=1), reported
    as a **cross-check** next to `--corr-factor`, not a source of truth for it — it doesn't
    capture geometric distortion or real dither-pattern effects, so it can (and does, e.g.
    ACS: analytic 1.5 vs the empirically-measured 1.24 in *Drizzle correlated noise* above)
    disagree with the measured factor. Stamped as `CASR`.
  - Both diagnostics, plus the recentring offset, drizzle pass, centring source, and a
    **`bcfill` flag naming the reduction the stamp was cut from** (mirroring its `BCFILL`
    header card), are recorded per (sample, lens, filt) in `info/lens_cutout_qc.json` —
    structured provenance for what shaped a given cutout, not just the sci/noise arrays
    themselves. `info/lens_masks.json`, `info/lens_arc_masks.json` and
    `info/lens_dataset_subplots.json` carry the same `bcfill` key, for the same reason.

## Masks (`scripts/make_masks.py`)

```bash
uv run python scripts/make_masks.py --sample slacs_gold                # best band per lens
uv run python scripts/make_masks.py --lens J0008-0004 --filt f606W     # force the draw band
uv run python scripts/make_masks.py --lens J0008-0004 --filt f160W --force  # redo one band
uv run python scripts/make_masks.py --sample slacs_gold --filt f555W   # review f814W's mask
uv run python scripts/make_masks.py --lens J0008-0004 --propose-from none  # blank canvas
```

Interactive, manual tool — not part of the automated per-lens pipeline above. For each lens
in a sample it opens PyAutoLens's `Scribbler` GUI (a Tk window; scribble the areas to
**REMOVE** from the fit — painted = excluded, unpainted = kept, the *opposite* polarity to
`autolens_workspace:scripts/imaging/data_preparation/gui/mask.py`; `Esc` to finish) over the
cutout's science image, writes `cutout_[cr_]mask.fits`, and records provenance per (sample,
lens, filt) in `info/lens_masks.json`. A lens that already has a mask is skipped (`--force`
to redraw), so a long GUI session across a sample is resumable.

**GITIGNORED BANDS ARE SKIPPED BY EVERY SWEEP (`--include-ignored` opts back in, 2026-09-22).**
`make_masks.discover_targets` — the single discovery used by `make_masks.py`,
`make_arc_masks.py`, `make_positions.py` and `detect_arcs.py`, with a matching guard in
`make_dataset_subplots.py` — now drops any cutout directory that `.gitignore` excludes, and
says so rather than shortening the list silently. **A gitignored band is by construction not a
science product**, so hand-drawn work written there is non-regenerable *and* invisible to every
clone. Today that means exactly the three **420 s SLACS SNAP f814W diagnostics**
(`slacs_other` J1134+6027, J1403+0006, J1538+5817) and nothing else — measured across all three
samples: `slacs_gold` 90→90 bands, `gallery` 33→33, `slacs_other` 37→34.

**The guard exists because it already went wrong.** J1538+5817's image positions were **marked
on its SNAP f814W band** and broadcast from there to f606W across a **688 mas (13.8 px)** band
misregistration — eight times the worst tie documented anywhere in this file (J0029-0055 f160W,
84 mas). Both copies were deleted 2026-09-22 and the lens dropped from
`info/lens_positions.json`; it is now unmarked and should be re-marked on f606W, the 4400 s
band with the visible ring. Nothing else was affected: no mask, arc mask or positions file
lives on any of the three SNAP bands (checked), and J1538+5817's f606W arc mask is `source:
drawn` with `proposal_from: None`, so no detector proposal was built off the SNAP frame either.

- **Implementation.** `cutout_paths.gitignored(paths)` — one batched `git check-ignore --stdin`
  call, not one per directory (a few hundred cutout dirs makes per-dir subprocesses the
  dominant cost). It lives in `cutout_paths` rather than `make_masks` because
  `make_dataset_subplots` needs it too and `make_masks` already imports *that* module — putting
  it in either would be an import cycle. **Any git failure returns the empty set**, so the
  guard can only narrow a run when it works and can never block one when it cannot; note
  `git check-ignore` exits **1** for "nothing matched", which is a normal answer, so only
  exit >1 counts as failure.
- **Consequence for `detect_arcs.py`:** those three lenses lose their red band, so the detector
  now reports `SKIPPED -- need a red and a blue band, have ['f606W']` instead of building a
  colour from a 420 s frame. That is the behaviour this file already asked for in prose
  (*"Do not quote F814W colours for this lens"* — J1538+5817, where the SNAP's own elliptical
  model puts a ring image at S/N −17.7). Verified to exit 0, not crash.

**Closing the GUI without painting writes NOTHING** — no FITS, no JSON entry — so the lens
stays pending and the next run re-offers it. An all-False mask would be a legitimate "exclude
nothing" mask, and that is exactly the ambiguity to avoid: an absent file means "not drawn
yet", where an empty one would silently mark the lens finished. (Under `--force` on a lens
that already has a mask, drawing nothing likewise writes nothing, so the existing mask
survives untouched.) `make_arc_masks.py` refuses an empty draw for a different reason — there
an empty arc region would mask out the whole stamp.

**All hand-drawn masks were deleted 2026-09-03 at the user's request** (54 files: 39
`slacs_gold` ACS — f814W plus one f555W — and 15 `gallery` f606W; at the time those ACS ones
lived in the separate `data/cutouts_bcfill/` tree), and `info/lens_masks.json` reset to `{}`. The old ones are
recoverable from git history if ever wanted; masking restarted under the current per-band
defaults.

**Current state (2026-09-13): 94 masks — `slacs_gold` is COMPLETE on every band**: f814W
38/38, f606W 22/22, f555W 16/16, f160W 13/13 — all in `data/cutouts/` beside the sci they
were drawn on, all recorded in `info/lens_masks.json` (whose per-band `bcfill` key says which
reduction that was). Most bands after the first came through the reviewed-proposal route
rather than fresh draws (f606W 21 `edited_from_f814W` + 1 `drawn`; f160W 9 + 4). **The last
gap, `J0822+2652 f606W_v2`, is now drawn too** — the split-visit second visit is its own
product directory and so needed its own draw, which is why it trailed the f606W sweep.
**`slacs_other` f814W is now drawn as well (4/4)** — that sample's ACS band, and its only
band with a bcfill reduction behind it. Its f606W (0/24) and f160W (0/6) are still open. **`gallery` f606W is now
COMPLETE too (15/15, drawn 2026-09-14/15)** — its F814W/F438W (0/12) are still open. Per-lens
notes follow, `slacs_gold`/`slacs_other` first, then `gallery`:
- **J1451-0239 was pulled back from its brighter lensed image** (2026-09-04): masked pixels
  at `r < 2.0″` of the deflector and `r < 0.5″` of image A were cleared (8930 → 8531 px),
  because the drawn edge sat 0.16″ from image A's centroid and covered 18–25% of the pixels
  within 0.3–0.5″ of it — clipping that image's PSF wing asymmetrically. **Worth checking on
  any lens whose mask reaches inward past ~2 θ_E**: measure the distance from the mask to
  each marked/known lensed image, not just how the stamp looks. The diffuse SW "fluff" the
  mask covers is correctly excluded — vetted as not lensed (→ memory:
  contaminant_vetting_j0728_j0841).
- **J1451-0239 f555W is UNCERTAIN and should be treated as provisional** (flagged by the
  user on drawing it, 2026-09-13). Where the mask boundary belongs on this lens is a genuine
  judgement call, not a settled one, so **any result that is sensitive to it needs a
  refit with the mask varied** before it is quoted. Provenance is consistent with that:
  it is `edited_from_f814W`, and the f814W proposal was heavily reworked rather than accepted
  — 8456 px proposed → 7222 px kept, with 2633 px erased and only 259 added, the largest
  erase of any f555W mask in the sample. It is also the same lens as the note above, whose
  f814W mask had to be pulled back off lensed image A; the two cautions compound, so check
  the f555W boundary against image A's centroid too.
- **J1020+1122: the mask SIZE is correct, the BOUNDARY is uncertain — and most uncertain in
  f160W.** Two things that are easy to conflate, so keep them apart. *The size is settled*:
  **two large nearby galaxies contaminate the right side of the frame**, so a mask covering
  ~47% of the f814W stamp (`n_excluded_px` 27173) is that big because the contamination is —
  confirmed by the user, and it retires the older "needs redrawing" flag, which was raised on
  stamp fraction alone before anyone looked at what was in the frame. **`n_excluded_px` is not
  on its own evidence of a bad mask; check the image before flagging a large one.** *Where the
  edge falls is not settled*: the user flags this lens's masks as uncertain across bands and
  **f160W most of all** (2026-09-13), so treat them as provisional and refit with the mask
  varied before quoting anything sensitive to it. The provenance shows the disagreement — all
  three bands were **drawn fresh, each rejecting the f814W proposal outright** (`source:
  drawn`, `n_erased_px` 0), and they do not agree on sky area: **67.93 arcsec² excluded in
  f814W, 53.39 in f606W, 66.10 in f160W**, i.e. f606W excludes ~21% less sky than f814W for
  the same contaminants. Per-band `note` keys in `info/lens_masks.json` carry this.

- **J1016+3859 (`slacs_other`) f814W: the obvious "arc" is not one — check before you spare
  it.** A compact red object 3.68″ from the deflector reads as a lensed image and was left
  unmasked for that reason; it is now masked. It really is **tangentially elongated** (b/a
  0.70, major axis 89° from the lens-radial direction, stable across annuli and surviving
  removal of the deflector's radial envelope in both f814W and f160W) — but that is **weak
  shear, not strong lensing**. Three things settle it: **no counter-image** anywhere at
  2.5–5.5″; it is **~3× REDDER than the deflector** in f160W/f814W (and 3× fainter relatively
  in f606W) where a lensed SLACS source is a blue star-forming galaxy; and 3.68″ is a wild
  θ_E for a single early-type. The consistent reading is a background galaxy just outside the
  caustic, singly imaged and sheared ~20%. That is self-checking: the observed ellipticity
  implies **θ_E ≈ 1.3–1.4″** (SIS, γ = θ_E/2r), and the faint f606W knots at **r = 1.4–1.8″**
  — the *real* lensed features, invisible in f814W — sit at that same radius. (Three-band and
  mask-check figures were written to the now-untracked `diagnostics/feature_vetting/`; regenerate rather than
  expect them in a clone.)
  **Generalise the method, not the verdict**: tangential elongation alone does not make an
  arc, and colour-relative-to-the-deflector plus a counter-image search is the cheap test.
  Compare colours as a *ratio against the deflector's own* — raw cross-band flux is
  meaningless here (F606W is WFPC2 DN/s, F814W ACS e/s; see *F606W `BUNIT`* above).
  - **A SECOND sheared galaxy, object B, at r=3.94″ (sky NW; dRA −3.67, dDec +1.43),
    masked in f606W 2026-09-13.** Same verdict, reached the same way, and it sharpens the
    method. Two objects at nearly the same radius on the same side invite a **fold-pair**
    reading — merging images straddling the critical curve, which is what a knotty arc looks
    like — and **colour kills it outright: A is redder than the deflector and B is bluer**
    (f606W/f814W relative to the lens 0.30 vs 1.42). *Two images of one source must share a
    colour*, so A and B cannot be images of anything. B also has no counter-image (f606W S/N
    −0.2 at its mirror position vs 3.7 at B), and its elongation (b/a 0.75 at 82° from radial
    in f814W, 0.56 at 70° in f606W) matches the **0.71 that weak shear predicts at 3.94″ for
    θ_E 1.35″**. B is faint (S/N 5–6) so its shape is noisy — the counter-image null and the
    radius carry it. Bonus: two independent galaxies at r≈3.7–3.9″ both sheared ~20% is a
    small weak-lensing measurement, and it agrees with the θ_E the f606W knots give.
    **B is detected in all three bands but masked only in f606W** — add it to f814W when that
    mask is next touched (it is outside the default 3.5″ aperture, so nothing is wrong today).
  - **J1538+5817 (`slacs_other`): two faint features near the ring, C1 and C2, are
  UNRESOLVED — left unmasked deliberately.** The lens itself is a **ring at θ_E ≈ 1.0″**
  (three images, PA ≈ 340/145/260, f606W S/N 37/33/17) and the drawn mask correctly clears
  it (0% masked inside 1.2″, nearest masked px 1.20″). Two marginal features sit outside it:
  **C1 (col 136, row 143; r=1.44″)** and **C2 (col 117, row 161; r=2.08″)**, both raw peak
  S/N ≈ 4, smoothed 2.7. C2 is tangential (80° from radial) but at twice the ring radius, so
  it cannot be another image of the ring source; C1 at 47° is neither radial nor tangential.
  **Spiral arm of the deflector vs lensed emission is NOT settled.** The F814W SNAP frame was
  recovered to try the colour test and **it failed**: at 420 s it is ~10× shallower than the
  4400 s F606W, the ring is at or below its limit, and its own elliptical model over-subtracts
  — one ring image reads **S/N −17.7**, which is model error, not flux. Do not quote colours
  for this lens from that product. Both features are left unmasked: C1 is only 0.44″ outside
  the ring, so masking a marginal detection there risks eating real lensed flux.
- **J1416+5136 (`slacs_other`): θ_E is 1.37″ and the two bright objects are NOT at it — but
  whether they are lensed is still OPEN.** Measured θ_E = 6.08 kpc = **1.37″** (Auger+2009;
  the system has a successful lens model, so an arc definitely exists). The two objects sit at
  **r = 2.2″ and 2.6″, i.e. 1.6–1.9 × θ_E**, neither has a counter-image (mirror of the bright
  one reads S/N −4.4 against 264 at the object), and both are resolved galaxies (2.08× and
  1.55× the PSF, b/a ≈ 0.92). A **weak bridge joins them** and survives MGE subtraction: +2σ
  at the midpoint against +0.1σ 1.2″ perpendicular, with the ends dominated by each object's
  own halo. **Open hypothesis, to be tested by modelling:** that this diffuse material is
  *extended source emission offset from the main source centre* — an extended source does
  image over a range of radii, so 1.6–1.9 × θ_E is large but not impossible. The plan is to
  fit both with the object-2 knots included and excluded and compare. **Do not treat the
  "companion galaxies" reading as settled**, and keep any mask outside **~1.6″** so it cannot
  touch the 1.37″ arc.

- **`slacs_other` f606W, three masks re-checked against the Auger+2009 θ_E (2026-09-21):
  J1134+6027, J1251-0208, J1403+0006. All three are CLEAR OF THE RING — 0% masked inside
  1.2 θ_E on every one — and only J1403 needs anything.** Figures:
  `diagnostics/masks/<lens>_f606W_mask_check_2026-09-21.png` (mask boundary + θ_E and 1.5 θ_E
  circles over as-observed / radial-subtracted / 0.3″-smoothed panels).
  - **J1134+6027 (θ_E 1.10″) — clean, nothing to do.** 9252 px (16% of stamp), nearest masked
    pixel **1.60″ = 1.45 θ_E**. **No unmasked S/N>4, ≥6 px detection anywhere out to 6″**: the
    one bright thing near the lens, a galaxy at **2.14″ = 1.94 θ_E, peak S/N 89**, is masked
    and has **no counter-image** (mirror S/N −0.6). No arc is detected at θ_E, and what the
    0.3″-smoothed panel shows there is a **quadrupole** (positive at PA 20–80° *and* 220–300°,
    negative at 100–180° *and* 320–20°) — i.e. the ellipticity residual the radial subtraction
    always leaves, not lensed emission. Do not mistake it for an arc.
  - **J1251-0208 (θ_E 0.84″) — clean, and the mask is nowhere near the ring.** 10419 px (18%),
    nearest masked pixel **2.85″ = 3.39 θ_E**; 0% masked inside 2 θ_E. The **arc is intact and
    unmasked**: col 118 row 144, **r = 1.23″ = 1.46 θ_E, 384 px, peak S/N 21, tangential
    (72° from radial), b/a 0.53, mirror S/N −0.5**. There IS unmasked S/N>4 material at
    r = 3.3–4.3″ (the biggest 472 px at 3.61″), which the `gallery` audit rule would have
    masked — but **the deflector here is a SPIRAL**, and the 0.3″-smoothed panel resolves that
    material into a two-armed pattern centred on the deflector, not field objects. Leaving it
    is right; the consequence is that **this lens's light model has to cope with spiral
    structure**, which an elliptical Sérsic will not.
  - **J1403+0006 (θ_E 0.83″) — FIXED 2026-09-21; the mask now stops at 1.00″ = 1.20 θ_E.**
    Re-analysed on an **elliptical-isophote** model of the deflector (b/a 0.89, PA 151° image, fitted with the mask *and* the neighbour excluded),
    not the circular radial median — which matters, see the retraction below. Figure:
    `diagnostics/masks/J1403+0006_f606W_neighbour_tradeoff.png`.
    - **This lens has ONE science band.** f606W is WFPC2/WF3, 4400 s; the f814W beside it is
      the untracked 420 s SNAP diagnostic (→ memory: `j1403-f814w-recentre-box`). So there is
      **no colour test available here**, and any argument that needs one cannot be made.
    - **The neighbour**: centroid **r = 1.82″ = 2.19 θ_E**, σ_major 0.66″ (**5.2× the PSF**,
      which is σ 0.128″/FWHM 0.301″), ~1100 px, peak S/N 15. A resolved galaxy, not a knot.
    - **The 1.20″ inner edge cut straight through it**, leaving **427 px of its halo inside
      the fit region** — 3.45% of the deflector's flux inside 1″, brightest leftover pixel
      **11.5σ**. This is the *mirror image* of the J1451/J2342/J1116 problem: there a mask
      clipped a lensed image, here a mask stopped in the middle of a **contaminant**.
      **The fix applied**: added 308 px — the neighbour's own **1σ isophote between 1.00″ and
      2.0″** — so 3123 → 3431 px. Nothing was cleared and the hand-drawn **outer** boundary is
      untouched; the edit is scripted, recorded in `info/lens_masks.json` and in the FITS header
      (`MASKEDIT`, `MASKEDDT`, `MASKEDPX`, `MASKEDRI`), and reversible with `git checkout`.
      `cutout_cr_dataset.png` was rebuilt (`make_dataset_subplots --force`) and the mosaics
      re-run, per the chain above.

      | inner floor | added | left inside 2″ | brightest leftover | % of deflector <1″ | cost 0.8–1.3 θ_E | cost 0.8–1.5 θ_E |
      |---|---|---|---|---|---|---|
      | 1.20″ (was) | — | 427 px | 11.5σ | 3.45% | 0.0% | 2.8% |
      | 1.10″ | +237 | 190 px | 10.7σ | 1.83% | 0.0% | 8.6% |
      | **1.00″ (applied)** | **+308** | **119 px** | **9.1σ** | **0.88%** | **6.4%** | **13.7%** |
      | 0.90″ | +358 | 69 px | 4.5σ | 0.38% | 11.8% | 17.3% |

      **Why 1.00″ and not 1.10″ or 0.90″.** The stopping criterion is the *systematics floor*,
      not the pixel count: this lens is galaxy-model-limited at 3.2× the photon noise (below),
      so nothing under ~10σ is believable anyway. 1.00″ is the loosest floor that pushes the
      brightest leftover under that floor (9.1σ); 1.10″ leaves it at 10.7σ, still above it, and
      0.90″ buys little more while reaching 1.10 θ_E — closer to the ring than any other mask
      in this repo (J1134 stops at 1.45 θ_E). The 6.4% it costs is **entirely in image
      PA 150–300°** (≈23% of the three sectors 180–270°, 0% in every other sector) — the only
      direction the neighbour reaches in. Figure:
      `diagnostics/masks/J1403+0006_f606W_mask_yours_vs_mine.png`. To move the edge, re-run the
      same script with a different floor; the table above is the price list.
    - **RETRACTED: the "blue tangential knot at 0.99 θ_E" reported on 2026-09-20 is not a
      feature.** It is at r = 0.77″, image PA 181°, **0.22″ from the neighbour's 1σ isophote**,
      and on the elliptical model it reads **+2.9σ against a 3.4σ model-error rms at that
      radius** — i.e. below the systematics. It was an artefact of **circular** radial-median
      subtraction on a lens with a bright overlapping neighbour, and its "4.8× bluer than the
      deflector" rested on a 2.2σ f814W measurement in the 420 s SNAP. *Generalise this*: on a
      lens with a close bright neighbour, `radial_median_subtract` is a finding aid that
      manufactures features at the neighbour's azimuth — fit an elliptical model before
      believing anything it shows inside ~2 θ_E.
    - **"No arc is detected here" STANDS, and now for a measured reason.** Inside 1.2″ the
      residual scatter is **3.2× the photon noise** — the image is *galaxy-model-limited, not
      noise-limited* (rms/noise by radius: 3.9 at 0.3–0.6″, 3.1 at 0.6–0.9″, 3.6 at 0.9–1.2″,
      falling to 0.9 by 3–4″). A feature therefore needs **≳10σ** to be a 3σ excess over model
      error, and the brightest thing in the 0.6–1.3 θ_E annulus outside the neighbour is
      **8.7σ**. Not a detection — and not a non-detection either: **this data cannot settle it**.
      What would: a second band (there is none), or a real Sérsic/MGE deflector fit rather
      than an isophote model.

- **`gallery`: every mask was audited against the PUBLISHED lens model (2026-09-15), and the
  audit is cheap to repeat.** All 15 F606W masks were checked with three tests — (1) is the
  object inside the image radii the Shu+2016 SIE model predicts, (2) is it tangentially
  elongated, (3) does it share the lensed source's colour — see *BELLS GALLERY → published
  lens models* below for where the model parameters come from. Result: **13 of 15 clean**; no
  masked lensed image, and **no unmasked object out to 6″** anywhere in the sample (every
  S/N>4, ≥6-px detection beyond 1.3× the outermost predicted image radius is masked on all 15
  lenses). Two lenses needed the mask pulled back off an image; both were fixed the same day
  — see the next two bullets.
  *Method note that made the audit trustworthy:* centre on the **catalogue position**, not the
  brightest smoothed peak — on J0201+3228, J0742+3341 and J0755+3445 a brighter neighbour
  hijacks the peak and every radius comes out nonsense.
- **`gallery` J2342-0120 f606W: the mask CLIPS the NE ring knot — pull it back.** The object
  at the screen bottom (col 140, row 198; r = 1.88″ = 1.7 θ_E) is correctly masked: it is
  **red** (F606W/F814W = 0.75 where all four ring knots read 1.8–2.8), radially rather than
  tangentially elongated (11° from the radius vector against 75–86° for the knots), and
  Shu+2016's own residual map leaves it unfitted. But its mask reaches to **0.15″ (1.1 PSF
  FWHM) of the ring knot at col 143, row 176**, masking 21% of the pixels within 0.3″ of that
  knot and **29% of its flux inside 0.5″** — worse than the documented J1451-0239 case.
  **Fixed 2026-09-15** the same way (scripted, not redrawn): masked pixels within r < 0.5″ of
  the knot centroid were cleared, **2818 → 2683 px**, nearest masked pixel **0.15″ → 0.50″**,
  flux of the knot inside 0.5″ masked **27% → 0%**, at a cost of 3.7% of the contaminant's own
  flux — it stays masked. Provenance is in `info/lens_masks.json` under `edited`.
- **`gallery` J1116+0915 f606W: right call on the contaminant, but the mask sits 0.27″ from
  the demagnified counter-image.** Of the two bright objects straddling the deflector, only
  the one at **col 155, row 111 (r = 1.58″, screen TOP = sky South)** is lensed; the one at
  **col 158, row 184 (r = 1.35″)** is a contaminant and is now masked — correctly: the model
  predicts images at 0.47″ and 1.65–1.70″ only, a source imaged at 1.35″ would need a partner
  at ~0.71″ only ~2× fainter (nothing is there), and Shu+2016's residual map leaves it
  standing. **The third feature, col 150, row 163 (r = 0.50″, peak S/N 6.8), is the
  demagnified counter-image — keep it.** The mask's nearest pixel was 0.27″ (2.0 PSF FWHM)
  from it, masking 18% of the pixels within 0.5″ and 6.5% of its flux. **Fixed 2026-09-15**:
  masked pixels within r < 0.5″ of col 150, row 163 cleared, **7599 → 7506 px**, nearest
  masked pixel **0.27″ → 0.51″**, flux masked **7% → 0%**, costing 0.6% of the contaminant's.
  It mattered more than the raw numbers suggest: that image's flux is the model's
  magnification-ratio constraint.
- **`gallery` J0918+5104 f606W: two real features near the quad — neither is a cosmic ray,
  and neither should be masked (2026-09-15).** (i) The object at **col 185, row 106
  (r = 2.25″)** is the **fourth image** of this cusp quad, not a CR: it is detected in three
  independent filter stacks at the same sky position to 0.4 px (F606W 26.0, F814W 9.4, F438W
  8.2 peak per-pixel S/N), is *broader* than the PSF (σ 0.076″ vs 0.055″ — CRs are sharper,
  never broader), is tangentially elongated (b/a 0.41, 89° from the radius vector), and sits
  at the radius the model predicts for image 4 (2.16″). (ii) The **diffuse patch ~0.9″ from
  it, around col 163–170, row 92–100 (r ≈ 2.3″), is real but its nature is OPEN.** Real:
  positive in **all four** input `flc` frames (145/111/236/182 e, 2.3–5.0σ each, no CR flags)
  where a blank-sky control scatters about zero, and +7.9σ over the local systematic floor in
  F606W. Open: its colour matches the arc (F606W/F814W ≈ 2.2 vs the 4th image's 2.9), but it
  is elongated only **33° from radial** where a lensed feature at 1.5 θ_E should be
  tangential, and de-lensing it puts its source **0.4″** from the real source — a separate
  clump, whose required counter-image (col 163, row 173, μ = 0.9) would be 3× fainter than
  the patch and so is undetectable either way. **Left unmasked deliberately**: fit it, look at
  the residual, and mask only if the reconstruction cannot produce it.
  *Reusable lesson:* **the per-exposure test settles cosmic-ray questions outright** — a CR
  is in exactly one input frame, so measure the same aperture in each `flc` against a
  blank-sky control in the same frames, instead of arguing from the drizzled stack.
  - **Two of my own claims here were wrong and are retracted**: an azimuthal-median model
    showed a "diffuse arc connecting the objects" that is mostly model residual (MGE removes
    most of it), and an apparent **"Einstein ring at r=0.83″, radius constant to 0.018″"** was
    a **centring dipole**, not a ring — the annulus reads +8σ at PA 0–100 and **−6σ** at PA
    140–220. *A real arc is positive on one side and ~0 opposite; it is never negative.* A
    radius that constant should itself have raised suspicion.
  - The galaxy model is the limiting factor, not the data: even with a 36-component MGE with a
    fitted centre (residuals within a few percent), a ~10σ dipole remains, because the MGE uses
    a single PA and real ellipticals twist. At 1.37″ the residual is +5σ at PA 240–280 and −6σ
    at PA 40–60 — the arc is in there and cannot yet be cleanly separated from model error.
    A Sérsic-plus-twist fit, or a real MGE code (`mgefit` is not installed), is the way on.
- **Orientation, which bit twice while working on this lens: on screen the stamps are NOT
    North-up-East-left.** The WCS is North-up in the usual sense (+column = **West**,
    +row = **North**), but the display puts row 0 at the top, so **North appears DOWN and East
    appears LEFT** — a *vertical flip* of the usual convention, **not** a 180° rotation.
    (Re-verified 2026-09-15 on J0237-0641 f606W: stepping +1″ East moves −25.2 columns, so
    East is the −column direction, i.e. screen LEFT. An earlier version of this bullet said
    "East appears RIGHT / a 180° rotation", which contradicted its own "+column = West" and
    was wrong.) Quote sky offsets (ΔRA/ΔDec), or **(column, row) pixels**, rather than compass
    words, and derive them from the WCS rather than from how the panel looks.
  - **Checked across the whole archive (2026-09-15): the orientation is uniform — there is no
    lens or band that breaks the rule.** All **317** cutouts in `data/cutouts/` *and*
    `data/cutouts_20arcsec/` (slacs_gold 180, slacs_other 71, gallery 66; ACS/WFC, WFPC2/PC,
    WFC3/IR and WFC3/UVIS) have, from the CD matrix: **PA(+row) = 0.000000° (North),
    PA(+col) = 270.000000° (West), axes orthogonal to 1e-6°, det(CD) < 0 everywhere**, and
    **no SIP or distortion-lookup keywords at all**. Pixel scales are uniform per instrument
    (ACS/WFC and WFPC2/PC 0.05″, WFC3/IR 0.06″, WFC3/UVIS 0.0396″). So `col = c0 − ΔRA/scale`,
    `row = r0 + ΔDec/scale` is exact on every product, and East-left/North-down holds
    everywhere on screen.
  - **TRAP that faked a violation on first pass: `RADESYS` is MIXED across the archive.**
    **61 products on 18 lenses are `FK5`** (J0216-0813, J0737+3216, J0912+0029, J0956+5100,
    J0959+0410, J1134+6027, J1143-0144, J1205+4910, J1250+0523, J1402+6321, J1403+0006,
    J1420+6019, J1538+5817, J1627-0053, J1630+4520, J2238-0754, J2300+0022, J2303+1422); the
    other 256 are `ICRS`. Building an offset position as a bare `SkyCoord(ra, dec)` — which
    defaults to **ICRS** — and feeding it to `world_to_pixel` on an FK5 header applies the
    ICRS↔FK5 frame bias, a **constant ~0.03″ (0.5 px) shift**. Differenced against an
    unshifted origin that reads as a **~1.5° rotation and a 2.6% scale error** that are not
    there. It cost a full wrong answer here before the CD matrices settled it. **Take the
    orientation from the CD matrix, and if you must build offset coordinates, build them in
    the header's own frame** (`SkyCoord(..., frame=c.frame)`), or difference two positions
    that both went through the same transform. The same bias matters for any cross-band or
    cross-catalogue position comparison at the tens-of-mas level.

**Where the arc should be: use the MEASURED Einstein radius, not an estimate.** For SLACS,
Auger+2009 (SLACS IX) modelled these systems and VizieR carries the result (for `gallery` it is
Shu+2016 Table 2, hand-entered — see *Tracking JSONs*):
`J/ApJ/705/1099/lenses`, keyed on `SDSS` (e.g. `J1416+5136`), column **`RE` — the Einstein
radius in kpc**, alongside `zlens`, `zsrc`, `sigma`, `MType`. Convert with
θ_E["] = RE[kpc] / D_A(z_lens)[kpc] × 206265. **A lens appearing in that table has a
successful lens model, which is itself the answer to "does this system really have an arc?"**
Fetch it directly (`curl` the `asu-tsv` endpoint); the Bolton+2008 table
`J/ApJ/682/964/table4` has the classification and σ but **no Einstein radius**.

Measured values for the lenses worked on here, against what image analysis found:

| lens | θ_E measured | found in the images |
|---|---|---|
| J1538+5817 | **1.00″** | ring measured at 1.00″ — exact |
| J1100+5329 | **1.52″** | images at 1.50/1.63″ — 2% |
| J1016+3859 | **1.09″** | 1.35″ inferred from shear — 24% high |
| J1251-0208 | **0.84″** | blue arc knots at 1.22″ — 45% high |
| J1403+0006 | **0.83″** | no arc found; masked companion at 1.77″ = **2.1×θ_E** |
| J1416+5136 | **1.37″** | no arc found; two objects at 2.2/2.6″ = **1.6–1.9×θ_E** |
| J1134+6027 | **1.10″** | no arc found; the N/S pair at 2.00/3.79″ |

**This is the cheapest way to settle "is that bright thing the lensed source?"** — on the last
three it confirmed, independently of colour or counter-image searches, that the eye-catching
objects sit well outside θ_E. Two cautions: an *extended* source arcs outside θ_E (J1251 at
1.45×), so treat θ_E as a floor not an exact locus; and a σ-based SIS estimate is a poor
substitute where the table has no row — **measured over all 89 rows, θ_E/θ_SIS has median
1.11 with a 16–84% range of 0.84–1.40** (worst: J1110+3649 0.25, J1100+5329 2.63), so the SIS
runs ~10% low *on average* and is useless per-lens. **Any note claiming a fixed correction —
e.g. "the SIS runs 30–40% low" — is wrong and is being corrected where it appears.**

**AUDIT of `info/lens_einstein_radii.json` (2026-09-21) — the file itself is NOT stale.**
Three checks, all clean, so do not re-derive it:
1. **Fresh against the source.** Re-fetched `J/ApJ/705/1099/lenses` from VizieR and diffed:
   **0 differences** across all 74 Auger rows and every field (`theta_e_arcsec`, `RE_kpc`,
   `zlens`, `zsrc`, `sigma_kms`); the 15 Shu rows were correctly left alone by the merge.
2. **Internally consistent.** All 74 Auger rows reproduce θ_E from `RE_kpc` on flat
   H0=70/**Ωm=0.3** to <0.6 mas, and all 15 Shu rows reproduce `RE_kpc` from θ_E on flat
   H0=70/**Ωm=0.274** to <15 pc. The two cosmologies are deliberate and opposite in direction
   (SLACS kpc native, gallery arcsec native) — do not "fix" one to match the other.
3. **Complete** for every lens with cutouts, bar the three already documented
   (`J1259+6134`, `J2141-0001`, `J2302-0840`).
   Also checked: all **54** `detector_theta_e_arcsec` values recorded in
   `info/lens_arc_masks.json` still match the catalogue exactly, so no arc proposal was built
   on a stale radius.

**What WAS stale was θ_E quoted OUTSIDE the catalogue, in two hand-written mask notes** — both
now carry a dated `CORRECTION` in `info/lens_masks.json`:
- **J1403+0006 f606W**: reasoned from SIS 0.75″ + a 30–40% correction to "expect the arc near
  ~1.0″". Measured θ_E is **0.83″** (the SIS was 11% low, not 30–40%), so the mask's 1.20″
  edge is **1.45 θ_E**. See the per-lens note above.
- **J1016+3859 f814W**: "the implied shear gives θ_E ~1.3–1.4″, matching the faint f606W knots
  at r = 1.4–1.8″". Measured θ_E is **1.09″**, so the shear inference ran **24% high** and
  those knots are at 1.3–1.65 θ_E, not at θ_E. The verdict on object A is unaffected (it rests
  on the counter-image null and the colour), but the 1.3–1.4″ is not this lens's θ_E.

**The lesson both share: a θ_E *inferred* from shear or from σ is not a measurement.** Look the
lens up in `info/lens_einstein_radii.json` first, and only fall back to an estimate if it has
no row — the catalogue has been there since 2026-09-17 and both notes predate it.

**Two ways of looking that have each produced a WRONG mask verdict here. Both are easy to
repeat, so check against them before recommending anything.**

- **An integrated-aperture S/N invents sources that are not there.** Summing
  `flux/sqrt(sum(noise^2))` over a ~100-pixel aperture turns a faint positive background
  residual — an imperfectly subtracted envelope gradient, say — into a confident-looking
  detection. It produced a phantom "counter-image" at J1403+0006 PA300 (integrated S/N 35.9,
  actual **peak S/N 1.9**, zero pixels above S/N 4) and another at J1251-0208 PA240
  (integrated 11.5–12.8, **peak 3.0**, one pixel above S/N 4), and it inflated the J1538+5817
  colour work. **Judge a detection on per-pixel S/N and the count of connected pixels above
  threshold**, never on an aperture sum alone; quote the peak alongside any integrated number.
- **Compact-source detection is structurally blind to a diffuse arc.** A `S/N > 4` per-pixel
  cut with a 3x3 opening finds knots and misses low-surface-brightness lensed emission
  entirely. At J1416+5136 that produced "no ring anywhere" and a recommendation to mask two
  objects — when smoothing the deflector-subtracted image to 0.3" reveals a coherent ridge
  over **~220 degrees of azimuth**, fitting a circle centred **0.58" from the deflector**
  (compare J1100+5329's neighbour galaxy at **4.13"**, which really was unlensed), with the
  diffuse light measuring **9% bluer than the deflector** once the bright knots are excluded.
  **Before concluding a lens has no arc, smooth and look for extended structure**, and fit the
  ridge curvature — an arc curves about the deflector, a tidal bridge or neighbour does not.
  Beware the converse: heavy smoothing of an imperfectly subtracted galaxy makes coherent
  residuals of its own, so confirm with colour and with a better galaxy model.

**One draw per BAND; the next band starts from a REVIEWED proposal, not a blank canvas
(`--propose-from`, default `auto`, added 2026-09-04).** Each run draws on one band — the
highest-S/N available (priority `f814W>f606W>f555W>f160W>…`, shared with `make_positions`
via `make_masks.pick_display_filt`) unless `--filt` forces it — and writes the mask for
**that band only**; run again with `--filt <band>` for the next, the skip check being
per-band so a sample sweep resumes cleanly. **Per-band drawing reverses the 2026-09-03
broadcast-by-default model** (changed the same day, at the user's call): a mask marks a *sky*
region, so sharing it is well-defined, but what a mask should *exclude* is not in practice
band-independent — a contaminant can be bright in one filter and absent in another, and each
band's depth and PSF move where the sensible boundary falls.

That is why a mask already drawn for another band is **proposed and reviewed** rather than
written out blind. `find_proposal_mask` reprojects it onto the band about to be drawn and
`draw_mask_gui` **outlines** it (1-px boundary, `mask_boundary`) over the usual two panels:
*radial-subtracted* (where the arcs are visible at all, so you can see what the inherited mask
may be clipping) | *as-observed* in this band (the contaminant's real extent and the galaxy
envelope). **A third "proposal APPLIED" panel — the same view with the mask interior blanked —
was built and then dropped the same day at the user's call**: it differed from the as-observed
panel only by the fill, which showed nothing the boundary does not already carry while hiding
the very pixels being judged and costing the other panels a third of the window. Where the
edge falls against this band's structure is the whole question, so the outline is the review.
- **Two brushes, from `al.Scribbler`'s two built-in scribble segments** (nothing new to
  maintain): `'1'` GREEN **adds** to the mask, `'2'` RED **erases** from it, on whichever
  panel you like. `show_mask()` returns only segment 1, so the code reads `get_scribble_masks()`
  after the blocking constructor returns. The mask written is
  `(proposal | painted) & ~erased`. `make_arc_masks.py` inherits the eraser (its arc region
  is `painted & ~erased`), which is the only change that tool needed.
- **The decision comes AFTER the GUI closes**, so it is made having seen the mask over this
  band's own image: `[a]` apply proposal+edits (default), `[d]` keep only what you drew
  (rejecting the proposal outright), `[s]` skip and write nothing. A non-interactive stdin or
  Ctrl-C gets `skip` — never a silently-written mask nobody approved.
- **Source (`auto`)**: the band's **own** existing mask when `--force`-redrawing it (so
  refining a mask is not redrawing it — this is how to fix J1020+1122 below), else the
  highest-priority other band that has one. `--propose-from <band>` forces it, `none` disables.
- `--subtract-radial` is now **tri-state**: on when reviewing a proposal, off when drawing
  from scratch (a contaminant mask is judged against the real sky), explicit flag always wins.
- Provenance per band records `source` = `drawn` / `accepted_from_<band>` /
  `edited_from_<band>` / `reprojected_from_<band>`, plus `proposal_from`, `proposal_px`,
  `n_added_px`, `n_erased_px`, `n_excluded_px`. **`accepted_from_x` with zero edits is a real,
  deliberate outcome** and is deliberately distinguishable from a hand-drawn mask that
  happens to resemble one.
- **`--broadcast` is the old unreviewed route and is now mutually exclusive with
  `--propose-from`** (argparse errors; pass `--propose-from none` alongside it) — combining
  them would write a mask reviewed for f555W back over f814W's own.

Both routes share the same transfer, a rigorous **per-pixel WCS reprojection**
(`reproject_mask_bool`, nearest-neighbour via `scipy.ndimage.map_coordinates`, source
sci-header WCS → target sci-header WCS), **not** an array-index copy — because only the 0.05″
optical bands (f814W/f606W/f555W) share a grid (to <1px); **f160W is a genuinely different grid**
(0.06″/px, 200px vs 240px for a 12″ stamp — a common sky point is ~20px away by index), and only
reprojection places the mask correctly there (verified: a test box regrids 0.05″→0.06″
area-preserving, centroid within 0.048″ ≈ sub-pixel).
- Per band, the mask is written into that band's one cutout dir, beside the sci it was drawn
  on, and a lens is skipped if the *draw band* already has a mask there.

**Side-by-side display (`--side-by-side`, default on, active only with
`--subtract-radial`).** Shows both views at once — radial-subtracted left, as-observed right,
separated by a blank gutter — because each answers a different question: the subtracted panel
is where the arcs are visible at all, the as-observed panel is where a contaminant's real
extent and the galaxy envelope are. **You may scribble on either panel**; the two halves are
read back and UNIONed onto the single-band mask, so the same stroke lands at the same sky
position from either side (verified: identical masks from a left-panel and a right-panel
stroke). Each panel is stretched independently — their dynamic ranges differ by orders of
magnitude, so a shared scale would flatten one. Panel read-back goes through `fold_panels`
(stride `panel_n_x + _PANEL_GAP`), which is written for N panels though only 1 or 2 are ever
built — a reviewed proposal adds an outline to both panels, not a panel of its own.
`make_arc_masks.py` has the same flag, also
on by default. A panel-labelling title is applied on the first in-axes mouse move, because
`al.Scribbler` builds its figure and then blocks inside `__init__` — there is no
post-construction hook.

**Defaults to `cutout_paths.DEFAULT_SIZE` (12″) — the pipeline's standard
tree**, and the reason `data/cutouts/` is the one size-variant tree kept tracked in git (see
*Cutouts* above and `.gitignore`): everything else in it is regenerable from the drizzled
mosaics, but a hand-drawn mask is not, so it needs the same durability as a tracked product.
Pass `--size 20` to mask the untracked, regenerable `data/cutouts_20arcsec/` tree instead.

**There is nothing to pick — one cutout dir per band** (`--variant` was removed 2026-09-14
along with the parallel trees). For ACS f814W/f555W + WFPC2 f606W you are scribbling on the
bcfill image (the same geometry with its dead-column stripes filled — a cleaner image to draw
on); for f160W, the gallery and un-bcfilled lenses, on the standard one. Either way the mask
lands beside that sci and there is no second copy to keep in step.

### Dataset QC subplots (`scripts/make_dataset_subplots.py`) — auto-rebuilt with every mask

```bash
uv run python scripts/make_dataset_subplots.py --sample slacs_gold           # all lenses
uv run python scripts/make_dataset_subplots.py --lens J1451-0239 --filt f814W --force
```

Assembles the four FITS products a lens model actually consumes — sci, noise, PSF, mask — into
an `al.Imaging`, applies the mask, and writes PyAutoLens's own 3x3
`aplt.subplot_imaging_dataset` as `{prefix}_dataset.png` beside the sci it read. **It is the
only QC PNG in the pipeline that shows the MASK laid over the data/noise/S-N**, i.e. what the
fit will really see. Provenance per (sample, lens, filt) in `info/lens_dataset_subplots.json`;
all four components are required, so an unmasked band is reported and skipped, not an error.

- **All four components come from the band's one cutout dir** — there is no per-component
  tree search left to get wrong, and no way for the mask to describe a different reduction
  than the sci beside it. **The psf is the one exception, and only along `--size`**: it is
  resolved through `cutout_paths.psf_cutout_dir()`, which always points at the default-size
  tree because the kernel is trimmed by amplitude and is not size-keyed, so a `--size 20` run
  pairs its stamps with that same kernel rather than reporting it missing.
- Panels use the pipeline's own look (inferno + `PercentileInterval(99)` + `AsinhStretch(0.1)`)
  via an `Axes.imshow` patch, since autolens's vendored plotter only knows linear/log norms; the
  dedicated log10 panels stay log10.
- **`make_masks.py` regenerates this PNG automatically for every band it writes a mask to**
  (`--no-dataset-subplot` opts out), so the subplot always describes the mask on disk rather
  than a previous draw. The hook is deliberately best-effort — the hand-drawn mask is the
  non-regenerable product, so a plotting failure prints the retry command and is swallowed
  rather than costing a mask that was just drawn. `make_arc_masks.py` does **not** trigger it:
  its product is `cutout_[cr_]mask_arcs.fits`, which the subplot doesn't read.

## Image positions (`scripts/make_positions.py`, `scripts/run_positions_all.sh`)

```bash
uv run python scripts/make_positions.py --sample slacs_gold            # best band per lens
uv run python scripts/make_positions.py --lens J0008-0004 --filt f606W # force the band
bash scripts/run_positions_all.sh                                      # every sample
bash scripts/run_positions_all.sh slacs_gold --force                   # one sample, re-mark
```

Interactive, manual tool (the positions analogue of `make_masks.py`). For each lens it opens
PyAutoLens's `al.Clicker` GUI over the lens's best cutout; you **double-click** each lensed
image of the source (2 for a double, 4 for a quad), each click snapping to the brightest pixel
within `--search-box-size` (default 5). The clicks are saved as an `al.Grid2DIrregular` in
`cutout_[cr_]positions.json` (via `al.output_to_json`) — the file a modelling script loads to
build a positions likelihood penalty (`al.PositionsLH`) that rejects mass models mapping the
images too far apart in the source plane — plus a `cutout_[cr_]positions.png` QC overlay per
band. Provenance per (sample, lens, marked-filt) in `info/lens_positions.json`.

**The broadcast goes THROUGH SKY, not by copying arcsec (fixed 2026-09-17).** The old
transfer wrote the same `(y, x)` into every band on the reasoning below — which is *nearly*
true and not exactly true. Measured on the 12 f160W pairs, the stamps' **tangent points differ
by 21–56 mas**, so the same arcsec offset is a different sky position in each band, up to a
whole f160W pixel. `reproject_positions` now does pixel → world → pixel per band, which
removes that term exactly for any pixel scale, stamp size or tangent point (verified: identity
within a band, exact round-trip, and the correction equals the independently measured
tangent-point offset). Applied to the 34 already-marked lenses with
`make_positions.py --rebroadcast`, which re-derives every non-marked band from the marked
band's own JSON without reopening the GUI — the clicks were fine, the transfer was not.
Corrections landed at **12–71 mas (0.25–1.41 px)**.

**F160W IS NOW TIED to each lens's reference band (done 2026-09-17/18).**
`align_wfpc2_to_acs.py` ties WFPC2 F606W to ACS; F160W was left trusting the delivered MAST
WCS because that script's docstring asserted ACS and WFC3/IR "agree with each other to
~0.01"". They do not, unless the f160W itself carries a GAIA fit, which most SLACS F160W does
not — so **10 of the 19 f160W bands were shifted onto their lens's reference band**:

| | |
|---|---|
| gold, ref f814W | J0029-0055 72 mas, J1430+4105 83, J1029+0420 51, J1023+4230 49, J1020+1122 35, J1032+5322 34, J0903+4116 21 |
| other, ref f606W | J1636+4707 98 mas; J1251-0208 and J1134+6027 measured 0 (already tied) |

Cutout-level tie across **all 19 f160W bands is now median 0.01 px (0.7 mas)**, from a
pre-alignment median of 8.6 mas. The one remaining outlier is **J1134+6027 at 46 mas
(0.77 px), and only against its own f814W — which is the untracked 420s SNAP diagnostic**,
not a science product; its f160W↔f606W tie, the one every tracked product on that lens uses,
is 0.

**MEASURE THE TIE WITH `align_wfpc2_to_acs.stable_centroid`, NOT A FIXED-BOX CENTROID.** A
plain flux centroid in a 0.6″ box is pulled by the lensed ring and by the broad IR PSF, and
over-reported these ties by tens of mas — it called J1251-0208 65 mas and J1134+6027 45 mas
when the ring-robust iterative centroid puts both at 0, and it under-called J1032+5322 (0.42
vs 0.56 px). The iterative windowed centroid re-centres each pass and converges independently
of the starting guess; it is the estimator the alignment itself uses, so it is the one that
says whether alignment is needed.

**Re-cutting after an alignment is an EXACT INTEGER PIXEL TRANSLATION of the stamp — verified
0.00e+00 max difference over ~39,800 overlapping pixels on all 8 re-cuts.** The WCS is
corrected continuously but the cut snaps to whole mosaic pixels, so a 0.35–1.63 px astrometric
correction produces a 0–2 px stamp shift, and **the two do not match**. Consequences:
- **Anything already drawn on that stamp must move by the ARRAY shift, not by a WCS
  reprojection.** `reproject_mask_bool` is the wrong tool here — it assumes only the grid
  changed, whereas re-cutting moved the astrometry *and* the cut window. On J1636+4707 it
  would have shifted the mask 1.63 px when the pixels had not moved at all, and on
  J1430+4105 it disagreed in sign. The 8 f160W contaminant masks were shifted by the
  cross-correlation integer instead (`MASKSHFT` header card records it).
- Positions were re-derived with `make_positions.py --rebroadcast` (f160W moved 14–52 mas).
- **`make_cutouts.py` does NOT write `cutout_dataset.png` — `make_masks.py` does**, so a
  re-cut (or any mask edited outside the GUI) leaves that QC subplot showing the OLD stamp
  and the OLD mask position, silently. Rebuild it explicitly:
  `make_dataset_subplots.py --lens <L> --filt f160W --force`.
- **`data/mosaics/` is tracked and is built from the cutouts**, so the per-sample QC grids go
  stale too: `make_mosaics.py --sample <sample>`.

  **The full chain after an alignment is: align → `make_cutouts` → mask shift →
  `make_dataset_subplots --force` → `make_mosaics` → `make_positions --rebroadcast`.**
  Skipping either of the middle two leaves a tracked PNG that disagrees with the FITS beside
  it, which is exactly the kind of quiet inconsistency this file exists to prevent.
- `GSC240FX=True` in the drizzled header marks a corrected band. `WCSNAME` does **not** —
  the script shifts CRVAL without renaming, so a corrected band still advertises its original
  solution.

Re-run the alignment safely at any time: it is idempotent (re-measures and applies ~0).

**What the broadcast CANNOT fix: a misregistered band.** Separately from the tangent point,
the same physical deflector sits at different sky coords in different bands' WCS — 9–84 mas
across the f160W pairs. That is an astrometry problem upstream, not a broadcast problem, so
`process_lens` measures the deflector tie per band and **warns** (`deflector_tie_mas` /
`deflector_tie_px` in the provenance) rather than absorbing it. Worst cases, all f160W:
**J0029-0055 84 mas (1.40 px)**, J1023+4230 51, J1029+0420 51 — J0029-0055 confirmed against
an independent **field source** (2.06 f160W px), which is the check `detect_arcs` already tells
you to run before believing a band tie. Do not re-mark these; the clicks are right and the
pixels are shifted. The `f555W`/`f606W` ties are tight (median 0.13–0.14 px), so this is an
f160W problem.

**One mark per lens, then broadcast to every band.**

**Four `slacs_gold` marks were DELETED at the user's request (2026-09-21): J0912+0029,
J1142+1001, J1143-0144, J2300+0022.** Both the per-band `cutout_cr_positions.json`/`.png` (8
JSONs, 8 PNGs across f814W plus f555W/f606W) and the `info/lens_positions.json` entries are
gone; `slacs_gold` is 33 → 29 marked lenses. **They were never committed**, so unlike the
2026-09-03 mask deletion there is nothing to recover from git history — re-marking is a
fresh `make_positions.py` GUI session. Each was a 2-image (double) mark, and three of the
four (all but J0912+0029) still carried `search_box_size: 5`, i.e. they predate the
2026-09-17 drop to 2 that exists precisely because a 5 px box lets a click move ~0.25″ onto
the deflector envelope or a neighbouring knot. Re-mark under the current default.

**`--search-box-size` default is 2, not 5 (2026-09-17).** 5 let a click move up to ~0.25″ on
ACS — far enough to land on the deflector envelope or a neighbouring knot instead of the image
you aimed at. Each click prints the distance the snap moved it, and `snap_moved_px_max` /
`snap_moved_px_mean` are recorded per lens. Note `al.Clicker`'s box is `range(p-n, p+n)`, i.e.
2n wide and half a pixel off-centre — upstream behaviour, unchanged here. Image positions are band-**independent**
sky coordinates: every band is cut about the same shared centre (`--center-band f814W`) and
pinned to the same output WCS (`final_rot=0`, tangent point at the lens), so a position in
**arcsec** relative to the stamp centre is identical in every band — even where pixel scales
differ (f160W 0.06″), because arcsec, not pixels, is the shared frame PyAutoLens works in. So
you mark once on the best band and the same `Grid2DIrregular` is written to every band, one
file per band into that band's cutout dir (exactly as masks route), under that band's pass
prefix. This is *simpler* than the
mask broadcast — positions need no reprojection, the arcsec values transfer verbatim; only
masks (per-pixel booleans on a specific grid) need the WCS regrid. **Positions still broadcast
unconditionally** (there is no `--broadcast`/`--no-broadcast` there, unlike the two mask tools):
a marked image is a sky coordinate that is literally the same number in every band, so there is
no per-band judgement to preserve. Trees/skip/`--size` behave as in `make_masks.py`; `run_positions_all.sh [SAMPLE] [flags…]` sweeps the
samples (all three by default). These hand-marked positions are non-regenerable, so like the
masks they live in the git-tracked `data/cutouts/` tree.

**Of the two files per band, only the JSON is a product: the PNG is gitignored
(2026-09-17).** `cutout_[cr_]positions.json` is the hand-marked work and is tracked like a
mask; `cutout_[cr_]positions.png` is a QC *view* of that JSON, one per band (~160), and is
rebuilt from the JSON plus the sci stamp with **no GUI and no re-clicking**:

```bash
uv run python scripts/make_positions.py --sample slacs_gold --overlays-only
```

That flag exists so the ignore rule is honest — without it "regenerable" would have meant
re-marking by hand, since `--force` reopens the Clicker. Verified byte-identical to what the
GUI wrote. Note this is the **one** ignored PNG under `data/cutouts/`: `cutout_cr.png`,
`_dataset.png` and `_mask_arcs.png` are all still tracked, so do not generalise the rule.

**Two panels by default, and the band's contaminants blanked (2026-09-17).** The GUI shows
the **radial-subtracted view LEFT and the as-observed view RIGHT** (`--no-side-by-side` for the
single subtracted panel, and the pair is only meaningful with `--subtract-radial`), because
each answers a different question: the subtracted panel is where an image buried in the galaxy
envelope shows up at all, the as-observed panel is where you judge whether that blob is really
there or is a subtraction artefact. **Double-click either — they are the same sky.**
- **Unlike `make_masks.py`, the panels are real matplotlib AXES, not an hstacked array**, and
  that is the whole reason no folding is needed here. The Scribbler reads one array back, so
  there the panels must be composited and the scribble folded up afterwards
  (`make_masks.fold_panels`). `al.Clicker` reads only `event.xdata/ydata`, so giving each panel
  its own axes over the *same* arcsec extent means a click already carries panel-local sky
  coordinates and is forwarded untouched. Verified: the same click on the left panel, the right
  panel, and the old single-panel GUI all snap to identical positions.
- **`--hide-masked` (default ON) blanks what `cutout_[cr_]mask.fits` already calls a
  contaminant**, outlines it, and drops it from *both* the stretch percentiles and the radial
  profile. Three reasons, not one: a neighbour bright enough to be masked owns the colour
  scale, invites a mis-click as a lensed image, and — the non-obvious one — drags up the
  per-radius median so the subtraction stamps a **dark ring clean across the arcs at its own
  radius**. Silently inactive on a band with no mask.
- **This is the one place a display lever can reach the saved product, so it is checked rather
  than trusted.** The snap runs on true flux, so a click at the edge of a blanked region can
  still land inside it; `process_lens` tests every saved position against the mask and prints a
  per-position WARNING naming the offender, with `n_positions_on_contaminant` in the
  provenance. Note this is the OPPOSITE call to `make_arc_masks.py`, which *outlines* the
  contaminant mask without filling it — there you are judging those very pixels, here they are
  the one region whose verdict is already in.
- `make_masks.radial_median_subtract` gained an optional `exclude=` for this; the default
  reproduces the unexcluded profile **exactly** (verified byte-identical), so `make_masks.py`
  and `make_arc_masks.py` are untouched.

**Three display-only levers make the arcs visible under the deflector light** (added
2026-09-03). The deflector sets the colour scale on its own, so on most lenses the images
you are trying to click are the least visible thing on screen. All three transform *only*
the imshow'd array — `al.Clicker` still snaps on the untouched flux array, so the saved
positions are bit-identical either way — and each is recorded in `info/lens_positions.json`
so a stamp's provenance says how it was marked:
- **`--subtract-radial` (DEFAULT ON since 2026-09-03; `--no-subtract-radial` opts out)** —
  subtract the deflector's azimuthally-averaged (median) radial
  profile. **The strongest of the three**: an elliptical's
  light is nearly a function of radius alone while the arcs are not, so the galaxy vanishes
  and the images stand out, *including the ones buried inside the envelope* that a central
  blank would hide along with the galaxy. Verified on J0330-0020 (all four images obvious)
  and J0008-0004. A finding aid, never photometry: real ellipticity leaves a quadrupole
  residual, and an arc biases the median at its own radius (self-subtraction).
- **`--mask-center ARCSEC`** — blank a disc at the stamp centre (which *is* the deflector,
  since `make_cutouts.py` recentres there). The blanked pixels are also **dropped from the
  stretch percentiles** — that second half is most of the gain, since the core otherwise
  pins `vmax` far above the arcs. The hidden disc is outlined so it is obvious what it
  covers. Note the snap is flux-based, so a click at the blank's edge can still land on a
  core pixel within `--search-box-size`.
- **`--vmax-value V`** — saturate above an absolute value in the displayed base's units
  (S/N by default, so `--vmax-value 10` is "everything above S/N 10 is red"); overrides
  `--vmax-percent`. This one is also applied to the QC overlay PNG, while `--mask-center`
  deliberately is not — the overlay should still show the deflector so you can judge the
  positions against the whole system.

The shared machinery lives in `make_masks.py` next to the other display helpers
(`stretched_display` gained `exclude=`/`vmax_value=`, plus new `radial_median_subtract`
and `central_disc`); its default call path is byte-identical to the pre-2026-09-03
implementation. `make_masks.py` takes `--subtract-radial` too but **OFF by default** — a
contaminant mask is normally judged against the real sky — while `make_arc_masks.py` (below)
has it **on**.

## Arc masks (`scripts/make_arc_masks.py`)

```bash
uv run python scripts/make_arc_masks.py --sample slacs_gold            # best band per lens
uv run python scripts/make_arc_masks.py --lens J0330-0020 --force      # one lens, redraw
uv run python scripts/make_arc_masks.py --lens J0330-0020 --detect-snr 3.0   # looser proposal
uv run python scripts/make_arc_masks.py --lens J0330-0020 --propose-from none  # blank canvas
```

The second hand-drawn mask tool, and the **opposite polarity to `make_masks.py`**: you paint
*only* the arcs / multiple images, and everything unpainted is masked out — the deflector
included. The product isolates the lensed source, for a source-only fit, an arc S/N
measurement, or a source-plane analysis.

| tool | you paint | product |
|---|---|---|
| `make_masks.py` | what to REMOVE (contaminants) | `cutout_[cr_]mask.fits` |
| `make_arc_masks.py` | the ARCS to KEEP | `cutout_[cr_]mask_arcs.fits` |

- **The saved array is in `al.Mask2D` convention (`True` = EXCLUDED) in both files**, so they
  load through the same reader and can't be swapped by polarity — but that means the arc
  mask is the **inverse of what you painted** (`saved = ~painted`; `~mask` is the arc
  region). Header records it: `MASKTYPE='ARCS'`, `MASKPOL`, `NARCPX`. A `make_masks.py`
  scribble, by contrast, is saved verbatim.
- **`--subtract-radial` is ON by default here** (the arcs are the subject and are usually
  invisible under the galaxy envelope), and images already marked by `make_positions.py` are
  **marked with a dark CROSS in the display** as a drawing guide (`--no-show-positions` off,
  `--cross-size` for the arm length; silently inactive where no positions file exists). Both
  are display-only.
  - **Crosses, not rings (2026-09-21, user's call).** Every other annotation on this GUI is a
    *closed boundary* — the arc proposal is a solid white outline and the contaminant mask a
    dashed dark one — so a dark ring read as a third small closed region, i.e. as one more
    thing masked OUT, which is the opposite of what a marked image means. Four open ticks
    cannot be read as an enclosed area at all, so the three annotations no longer have to be
    told apart by line style alone. The arms stop 2 px short of the centre for the reason the
    contaminant outline is not filled: the marked pixel is the one you are judging, so the
    marker points at it rather than covering it (verified: centre pixel untouched on all four
    J1250+0523 positions). Kept ~3 px thick so it is about as many display pixels as the 6 px
    ring it replaced (~50 vs ~38) — a cross of hairlines reads fainter than a ring of the same
    radius, and it has to survive being drawn over a bright arc.
- Refuses to write an empty draw — an empty arc region would mask out the whole stamp.
- **Where `--broadcast` is used, it reprojects the ARC REGION, not the saved array.**
  Outside-footprint pixels
  then default to "not arc" (masked out), the safe side of the stamp edge; reprojecting the
  inverted array would instead leave edge slivers unmasked. Verified across the f160W regrid
  (0.06″/200px vs 0.05″/240px): same sky area to 0.06%, same annulus radii to half a pixel.
- Everything else — one draw per **band** (`--broadcast` to share it, same default and same
  reasoning as `make_masks.py`), band priority, `--size` routing, skip/`--force`
  — is `make_masks.py`'s, imported rather than re-implemented. It gets the **red erase brush**
  (`'2'`) for free from the shared GUI helper — the arc region is `painted & ~erased`, so an
  over-painted stroke is trimmed rather than undone whole. Provenance in
  `info/lens_arc_masks.json`; a per-band QC PNG (`cutout_[cr_]mask_arcs.png`) outlines the
  region over the radial-subtracted image.
- **`--propose-from` IS now wired up, default `auto` (2026-09-15; union added 2026-09-16)** — the arc region you
  review comes from the automatic detector (`scripts/detect_arcs.py`, below), so arc masking
  is *correcting an outline*, not drawing one. The review loop and its code are
  `make_masks.py`'s: outline on every panel, green adds / red erases, and `confirm_proposal`
  after the GUI closes — `[a]` apply proposal+edits, `[d]` keep only what you drew, `[s]` skip
  and write nothing. Provenance gains `source` = `accepted_from_detect` / `edited_from_detect`
  / `drawn`, `proposal_px`, `n_added_px`, `n_erased_px`, and the detector's own numbers
  (`detector_bands`, `detector_theta_e_arcsec`, `detector_band_offset_px`, `detector_core_dipole`, per-component
  radii and peak S/N). `--propose-from <band>` proposes another band's arc mask instead —
  still the judgement call it always was (an arc detected in one filter is often simply absent
  in another), which is why it is reviewed and never broadcast — and `--propose-from none` is
  the old blank canvas. `--broadcast` is mutually exclusive with it, as in `make_masks.py`.
  `--detect-*` pass parameters through; `--detect-snr` (default 3.5) is the one to reach for.
- **Sources COMBINE BY UNION, which is what makes the second band cheap (2026-09-16).**
  `auto` = the detector **plus** the highest-priority other band that already has an arc mask,
  so the first band of a lens starts from the detector alone and every band after it starts
  from the detector *and your own reviewed region*. Both belong on screen: the detector sees
  what is blue in this lens now, your mask carries a judgement it cannot make, and neither is
  authoritative — an arc detected in one filter is often absent in another, which is precisely
  why the union is outlined for review and trimmed with the red brush rather than broadcast.
  An explicit comma-separated list works too (`--propose-from detect,f814W`), `none` is the
  blank canvas, and under `--force` the draw band's own mask joins the union so refining a mask
  is not redrawing it. Provenance records `proposal_from` (`detect+f814W`),
  `proposal_px_by_source` and `proposal_overlap_px`, and the breakdown is printed per lens.
  Verified with a hand-edited mask, not just a round-trip: a 900px box added by hand on f814W
  reappeared in the f606W proposal, and union = detect + inherited − overlap held exactly.
- **The CONTAMINANT mask is outlined while you paint arcs (`--show-contaminants`, default on,
  2026-09-16).** The thing most likely to be mistaken for a lensed image is a neighbour or
  field source that has *already* been judged a contaminant on that very band, so
  `cutout_[cr_]mask.fits` is drawn as a **dashed dark** boundary — which keeps three
  annotations apart at a glance: arc proposal = solid WHITE outline, marked positions = DARK
  CROSSES (rings until 2026-09-21, see above), contaminants = DASHED DARK. The interior is deliberately not filled, the same
  call recorded for `make_masks.py`: blanking hides the pixels you are judging. Display only
  and verified so (722 display pixels changed, written arc region untouched), silently
  inactive on a band with no contaminant mask. The detector already excludes those pixels, so
  a proposal never lands on one — checked on J1020+1122, whose mask covers ~47% of the stamp.
- **Two orientation traps, both found by test and both silent:** autoarray's native grid puts
  row 0 at **+y** (`row = cy - y/scale`), so a `cy + y/scale` position marker lands on the
  *mirror* of each image — plausibly near the lens, and on nothing; and `plt.contour` with an
  `extent` defaults to row 0 at the **bottom** while `imshow` puts it at the top, so the QC
  outline needs `origin='upper'` or it is drawn vertically mirrored.
- **`slacs_gold` f606W arc masks are PARKED OUT OF GIT (2026-09-21, user's call).** Its f814W
  arc masks (38/38) are committed; the 10 f606W ones drawn so far are held out by a
  `.gitignore` rule on `data/cutouts/slacs_gold/*/f606W/cutout_cr_mask_arcs.{fits,png}` —
  the only `_mask_arcs` files under `data/cutouts/` that are ignored (gallery's and
  `slacs_other`'s f606W, and `slacs_gold`'s own f814W, are all tracked, so do not generalise
  the rule). **Nothing was deleted**: the files stay on disk beside the sci they were drawn
  on and their provenance stays in `info/lens_arc_masks.json`, so `make_arc_masks.py` still
  counts those 10 bands done and skips them (`--force` to redraw) — the rule governs git, not
  the campaign. Delete the two lines to commit them.

## Automatic arc detection (`scripts/detect_arcs.py`)

```bash
uv run python scripts/detect_arcs.py --sample slacs_gold      # QC figures + summary, no products
uv run python scripts/detect_arcs.py --lens J0330-0020 --snr 3.0   # tune one lens
uv run python scripts/detect_arcs.py --fetch-theta-e          # refresh the Einstein-radius cache
```

Finds the arcs and hands the region to `make_arc_masks.py --propose-from detect` as a
**proposal**. It writes no mask and no product of its own — QC figures and a summary JSON
under `diagnostics/arc_detection/<sample>/` (untracked), plus a **contact sheet** of the whole
sample on one page, which is how to judge a sweep before opening any GUI.

**The discriminant is COLOUR, not brightness.** The deflector is a red early-type, the source a
blue star-forming galaxy, so the detector builds a **blue-excess image** `be = B' - k·R'`:
- `B'`, `R'` are the blue (f606W/f555W) and red (f814W) stamps **PSF-matched** by convolving
  each with the *other* band's kernel, so both carry `PSF_R ⊛ PSF_B` and the difference is not
  a sharp-minus-blurry dipole at every steep gradient.
- **`k` is keyed on the red band's own BRIGHTNESS, not on radius** — the deflector's colour on
  each isophote, measured as the median blue/red ratio in log-spaced brightness bins. Two
  earlier versions failed and say why the third works: a single scalar `k` leaves the galaxy's
  colour *gradient* behind as a broad disc exactly where the arcs are; per-**radius** binning
  needs an isophote shape, and real ellipticals twist and grow discs, so fixed elliptical
  annuli left signed lobes along the major axis (J0841+3824) and a noisy per-annulus ratio
  printed **bullseye rings** round the core. Brightness bins assume no shape at all. They must
  be **log-spaced**, not equal-population: quantile bins collapse the whole core into one bin
  and leave a broad negative bowl over the middle of every lens, which is where the inner arcs
  are. Because `k` is a measured ratio, the bands' different units (WFPC2 DN/s vs ACS e/s)
  cancel and are never converted.
- Candidates are kept only in an annulus around the **measured** θ_E, cached in
  `info/lens_einstein_radii.json` (see *Tracking JSONs* for the two catalogues and their
  differing cosmologies). **All 38 `slacs_gold` and all 15 `gallery` lenses have a row**;
  `slacs_other` is missing three. The same prior kills both classic false positives: the
  central residual inside, field galaxies outside. **No row → `propose()` raises**, so a
  missing θ_E takes the whole sample's `--propose-from detect` down rather than degrading.
- Components are judged on **peak per-pixel S/N and connected area**, never an aperture sum,
  and detection runs on a 0.07″-smoothed map so a diffuse arc is not missed — the two traps
  recorded above. A component must also be visible in the blue band on its own, so an
  over-subtracted red patch cannot pose as a source.
- **The mirror test, straight from this file's own rule that a real arc is never negative
  opposite.** A centring or PSF mismatch is antisymmetric about the deflector; lensing is not —
  opposite an arc lies a counter-image or empty sky, never a deficit. So a component whose
  180° mirror reads below `MIRROR_DEFICIT_SNR` (−2σ) is dropped as the bright half of a dipole,
  and `core_dipole` reports the antisymmetric fraction for the lens as a whole. It earns its
  keep: it removes a component on 24 of 38 lenses, and what it removes is overwhelmingly
  **inner (0.4–0.85 θ_E), RADIALLY elongated, with a −2 to −8σ hole opposite** — the dipole
  signature, not an arc's.

**The band astrometry is FINE, and a claim here that it was not is retracted (2026-09-15).**
An earlier version of this section said WFPC2 f606W sits up to 2.2 px = 0.11″ off ACS f814W.
**That was an artefact of the measurement, not a property of the data**, and the retraction is
worth more than the claim was: the number came from a flux-weighted centroid in a small
aperture, which in the blue band walks onto the arc, run on PSF-MATCHED frames, which are not
the sky. Three independent measurements agree the tie from `align_wfpc2_to_acs.py` holds:
- **cross-matched field sources** (43–175 per lens — the only check the deflector cannot bias,
  and the one this repo already prescribes): median vector offset **0.002–0.026″ = 0.04–0.5 px**;
- **deflector centroid by windowed centre-of-mass** (the estimator `align_wfpc2_to_acs.py`
  itself uses): **0.05–0.42 px**, and *identical in the mosaic and in the cutout*, so
  `make_cutouts.py` carries the WCS through correctly;
- the detector's own per-lens diagnostic: median **0.31 px**, max 0.69 px.

So **there is nothing to re-register**, and `detect_arcs.py` deliberately does not: it measures
the offset, records it (`band_offset_px`), and flags anything over 0.8 px for
`align_wfpc2_to_acs.py` — the tool that owns the WCS — rather than shifting pixels itself.
Fitting a shift by minimising the colour residual was tried too and is worse than useless: with
an imperfect PSF match a shift can always buy a smaller residual, so it ran to its limit on
half the test lenses and did not reduce the dipole it was meant to remove. **Rule: a
measurement that says the astrometry is broken is nearly always the estimator.** Check it
against field sources before believing it, and never against a PSF-matched frame.

**What WAS broken is the PSF kernels — and the cause was an x/y swap, fixed 2026-09-15.**
`photutils.centroid_com` returns **(x, y), column first**. Four sites unpacked it as `(y, x)`:
`make_psf.oversampled_to_kernel`, `psf_models._resample_centered`, and both centroid calls in
`make_psf_inject.py` (`_render_from_array`, `extract_kernel`). Each then recentred by the
**transposed** offset — shifting by `(dx, dy)` where it meant `(dy, dx)` — so every kernel kept
an antisymmetric residual, `dy ≈ −dx`. Before the fix, over all 149 cutout PSFs repo-wide:
median offset **0.33 px, max 1.73 px, 37 off by >0.5 px**, and `dy` vs `dx` anticorrelated at
**−0.76** — that anticorrelation *is* the transposition's signature, and it is what identified
the bug. `align_wfpc2_to_acs.py` had the same call **correct** (`cx, cy = centroid_com(d)`),
which is why the astrometric ties were never affected.

**Why it matters beyond the kernels: convolving with an off-centre kernel TRANSLATES the
image by that offset.** It is what the fictitious 2.2″ band offset above was really made of (in
a PSF kernel swap each band moves by the *other* band's kernel offset), and, far more
importantly, **any fit that convolves a model with these kernels was shifted the same way** —
1.7 px is 0.085″, which a lens model absorbs by moving the deflector, differently in each band.
**All 154 PSF products were regenerated after the fix (2026-09-16)**, measured three ways
(7×7 COM about the peak, whole-kernel flux centroid, quadratic core fit):

| sample | n | before: median / max / >1px | after: median / max / >1px |
|---|---|---|---|
| `slacs_gold` | 90 | 0.308 / 1.616 / **4** | 0.285 / 0.682 / **0** |
| `slacs_other` | 37 | 0.768 / 1.725 / **12** | **0.099** / 0.554 / **0** |
| `gallery` | 27 | 0.287 / 0.984 / 0 | 0.239 / 0.394 / **0** |
| **all** | **154** | 0.326 / 1.725 / **16** | 0.258 / 0.682 / **0** |

`slacs_other` was the worst affected and gained the most (median 0.768 → 0.099 px). **The
~0.26 px that remains is the measurement floor, not a residual error**: only **5 of 154**
kernels exceed 0.3 px on all three estimators (worst 0.68 px), i.e. at that level the
estimators disagree with each other about where an undersampled PSF's centre is.

Knock-on effects of the rebuild, all checked: **70 injected-tier `psf_err` maps** rebuilt with
`make_psf_err_injected.py` (`make_psf.py` only rewrites the empirical tier's, so skipping this
leaves error maps describing the old off-centre kernel — no stale ones remain); **2 products
changed tier**, both *to* `empirical` (`gallery/J0918+5104 f606W`, `slacs_other/J2302-0840
f606W`), where the centred ePSF now passes the core-vs-outskirt quality gate; 3
`slacs_other` f814W products still have **no** `psf_err` map, which is pre-existing (`no
calibration for slacs_other:f814W`), not something the rebuild removed. Dataset QC subplots
were regenerated too, since they embed the PSF. **Hand-drawn masks were NOT affected and did
not need redrawing** — `make_masks.py`/`make_arc_masks.py` read only `_sci`/`_noise` and other
masks, and a mask is read back as brush pixel positions. **`--bcfill` is likewise unaffected**:
`make_psf.py` has no variant key and `make_psf_inject.py` stages from `data/drizzle_files/`, so
every kernel comes from the standard tree, which is right — bcfill is an input-level
dead-column repair drizzled onto the same output grid, not an optical or resampling change.

**When a recentring step looks right but its product is not centred, check the coordinate
order before anything else** — and note a symmetric test PSF cannot catch this, because a
Gaussian centred on its own peak has equal x and y centroids inside the window: the test input
must be asymmetric, or real.

**One more data fact, silent and general:**
- **A few noise maps flag zero-coverage pixels with `1e8`** — 12 products repo-wide, a handful
  of pixels each (J0252+0039 f606W, J0728+3835 both bands, J0737+3216 f814W, J1451-0239 f814W,
  and 5 in `slacs_other`). Squared and pushed through an FFT that one value **destroys the
  whole frame**. Any convolution or smoothing of a noise map here must detect and heal them
  first (`SENTINEL_FACTOR`).

**What it achieves on `slacs_gold` (38/38 run, ~15 s for the sample), measured after the PSF
kernel fix:** every lens gets a proposal; **27/38 have their largest component at 0.6–1.5 θ_E**,
and the median lens has **100%** of its proposed area in that range — the proposal is the ring
and little else. Median proposal 501 px = 1.25 arcsec². The PSF fix moved these: before it,
25/38 on-ring and one lens (J0841+3824) proposed nothing at all.

**Where it does NOT work, and why that is visible rather than hidden.** The method's one
failure mode is a deflector subtraction leaving an antisymmetric residual, and `core_dipole`
measures exactly that: median **1.43**, above 2 on **11 of 38** lenses. On those the inner
region is model error, the mirror test correctly refuses to propose there, and **the inner
images have to be drawn by hand**. Centring the PSF kernels cut this materially (median was
1.80 and 16/38 exceeded 2 beforehand) — i.e. a good part of what looked like an irreducible
subtraction artefact was the off-centre kernels. What remains is the last few tenths of a
pixel of band mismatch plus genuine WFPC2 PSF error; a better f606W kernel would buy more here
than any detector tuning. **This is a proposal generator, not a detector to be believed** —
nothing is written without the GUI review.

**The band-offset diagnostic must compare RAW against RAW.** It reads the regridded blue stamp
against the red stamp on its own grid, *neither* PSF-matched. An earlier version compared the
PSF-*matched* red against the raw blue and therefore measured the kernels as much as the sky —
median 0.31 px, max 1.15 px. Like-for-like it reads **median 0.10 px, max 0.87 px**, and it
tracks the independent field-source measurement (J1430+4105 is the worst lens by both: 0.87 px
here, 0.68 px from 111 cross-matched field sources — still well under the f606W PSF, and
flagged rather than corrected).

**Other samples degrade gracefully, not silently:** a lens with no blue band or no measured
θ_E is reported and skipped (`gallery` is f606W-only and is not in Auger+2009 — Shu+2016
Table 2 has θ_E for all 15 if it is ever wanted; 4 of 24 `slacs_other` lenses run today).

## PSF generation (`scripts/make_psf.py`, `scripts/psf_models.py`)

```bash
uv run python scripts/make_psf.py --lens J0252+0039 --filt f814W
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
  35/41/51; block-reduced from the in-memory oversampled ePSF, which is not itself written to
  disk — nothing reads it back), **`psf.png`** (4-panel QA: star montage · kernel linear ·
  kernel **log** · radial profile).
- **the band's cutout dir** (modelling-ready): **`cutout_[cr_]psf.fits`** —
  the **trimmed** kernel, cut to the amplitude-`--trim-threshold` (default 1e-3 of peak)
  radius and pass-matched to `cutout_[cr_]sci.fits`. Written by `make_psf.py` itself.

  **Which cutout dir is `cutout_paths.psf_cutout_dir()`'s call, not the caller's**: the one
  cutout dir for the band, `data/cutouts/<sample>/<lens>/<filt>/`, so the four products a fit
  consumes — sci, noise, **psf**, mask — sit together by construction. *(Between 2026-09-04
  and 2026-09-14 this function also had to choose between the parallel bcfill and standard
  trees — 308 kernel/err/analytic files were moved into the bcfill one on 2026-09-04, 114
  stayed. The trees are now merged, so all it still asserts is the size-independence below.)*
  - There is still exactly **one kernel per band**, built from `data/drizzled/`, the standard
    mosaic — a kernel sitting beside a bcfill stamp is not a separately-measured bcfill PSF.
    **Building from `data/drizzled_bcfill/` instead was measured (2026-09-04) and is not
    worth it**: half the products come out bit-identical (no PSF star's stamp touches a
    filled column), the rest differ by ≈1× the kernel's own bootstrap error in the direction
    of a *bias* (interpolating across a PSF core, a sharp ridge, is exactly where the fill's
    local-linearity assumption fails) — and the bcfill tree has **no no-CR pass**, which the
    empirical build depends on (+74% stars). The injected tier is structurally identical
    (noiseless injected frames ⇒ dropping frames doesn't bias the drizzle mean; 1 of 60
    products has any zero-weight pixel within 1.5″ of the lens). → memory:
    psf_from_bcfill_no_gain
  - It is **orthogonal to `--size`**, which is now its whole job: `psf_cutout_dir()`
    deliberately takes no size and always resolves in the default-size tree, so a `--size 20`
    stamp still pairs with the one kernel instead of looking for one that by design is not
    there.
  - Every writer and reader goes through it: `make_psf.py`, `make_psf_inject.py` (incl.
    `_promote`'s move-aside), `make_psf_err_injected.py`, `make_cutouts.py --psf-err`,
    `make_psf_mosaics.py`, and `make_dataset_subplots.py` (which resolves the psf component
    through it, so the psf does not obey its `--size`).

For **empirical** builds, `make_psf.py` also writes a per-pixel **PSF error map** alongside
each kernel — `psf_kernel_err.fits` (full) and `cutout_[cr_]psf_err.fits` (trimmed, same grid
as the kernel) — from a bootstrap/jackknife over the star sample. See the *PSF uncertainty*
open-item bullet below for the method, JSON scalars, and the `make_cutouts.py --psf-err`
effective-noise folding that propagates it into a fit. Model-tier builds do not (yet) carry
one.

For a **model-tier** build, `make_psf.py` immediately auto-chains
`make_psf_inject.run_injection(..., promote=True)` (see *Drizzle-broadened model PSF by
injection* below) — so `psf_kernel.fits` / `cutout_[cr_]psf.fits` end up holding the
**drizzle-broadened injected kernel**, the more correct product, not the sharper analytic
model this section just described; that analytic build is kept alongside as
`psf_kernel_analytic.fits` / `cutout_[cr_]psf_analytic.fits` for comparison. An **empirical**
build is already the true drizzled PSF (cut from the mosaic) and is never touched by this.

Records `info/lens_psf.json` (`{sample:{lens:{filt:{method,n_stars,fwhm_pix,oversample,
kernel_size,cutout_kernel_size,trim_threshold}}}}`, `null` + exit 0 on no data). Run
**after** the drizzles (a PSF needs a mosaic). → memory: psf_kernel_sizing

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
- **`photutils.centroid_com` returns (x, y) — COLUMN FIRST — and unpacking it `(y, x)` put
  every kernel off its own centre (found and fixed 2026-09-15).** Four sites here did that and
  so recentred by the transposed offset: `oversampled_to_kernel`, `psf_models._resample_centered`
  and both centroid calls in `make_psf_inject.py`. The residual is antisymmetric (`dy ≈ −dx`),
  it reached **1.73 px**, and an off-centre kernel **translates whatever it is convolved with**
  — i.e. it shifts the model in every fit that uses it. Full account under *Automatic arc
  detection*; the short version is that a centred-looking recentring step is worth testing with
  an **asymmetric** input, since a symmetric one cannot fail. All PSF products were regenerated
  after the fix.
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
- **Crowding: DAOStarFinder's own `min_separation` is not enough** (`crowd_sep_frac`,
  `make_psf.reject_crowded`, 2026-07-31). `min_separation` (~25px for ACS/WFPC2) only stops
  DAOStarFinder double-detecting one blended peak; it's smaller than the `extract_stars()`
  stamp width (`star_size`, 51px for ACS), so two accepted detections can still be close
  enough for their stamps to overlap and mix a neighbour's flux into the ePSF wings. A
  greedy-by-flux cut (brightest of any pair within `crowd_sep_frac x star_size`, default
  `crowd_sep_frac=1.0` — stamps never overlap) runs after the flux floor. QA'd as
  `n_crowded_rejected` (console + `info/lens_psf.json`); intended to reduce, not replace,
  `info/psf_stars.json` manual exclusions on genuinely crowded fields.

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
`data/reference_files/wfpc2_f606w_psfdb/`. Every F606W product (incl. split-visit `f606W_v2`)
reuses it; records method `model_wfpc2_psfdb`. The cached ePSF is detector-frame; each lens
**rotates it to North-up at resample time via that lens's WF3 exposure CD** (so one shared
build serves every roll — see *Model PSFs are rotated to North-up*), and the pedestal
subtraction removes its ~1e-3 wing floor. Detector-frame still omits drizzle broadening, so
the lens's-own-field empirical build is preferred where stars exist — this is the model
**fallback, slotted ABOVE the STDPSF F555W proxy** (now fallback-of-the-fallback). **Traps:**
the MAST PSF DB `chip` column is the WFPC2 CCD/FITS ext (1=PC…4=WF4) — **query chip 3 for
WF3**, not chip 1 (PC has a different pixel scale *and* PSF); split-visit filter keys
(`f606W_v2`) must be normalised to the base filter before any STDPSF/pivot lookup
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
from the mosaic).

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

**The injected kernel is CANONICAL for the model tier, not a side product** (decided
2026-07-30, superseding the original "parallel comparison, nothing downstream changes"
design). Validated on F160W where an empirical truth exists on the same lens:
analytic-STDPSF 3.15px → **injected-STDPSF 3.5–3.9px ≈ empirical 3.8px** (drizzle broadening
recovered; the residual is the real optical wing STDPSF underestimates, not a broadening
error). ACS broadens little (native 0.05″, small resampling — FD 1.65→1.76px); F160W and
WFPC2 broaden more (0.1283→0.06″, 0.0996→0.05″). Since the injected build is strictly closer
to the true drizzled PSF, downstream modelling should read it, not the sharp analytic model.
→ memory: psf_injection_drizzle_broadening

`make_psf_inject.run_injection(lens, filt, sample, promote=None)` decides promotion by
reading the CURRENT `info/lens_psf.json` method for that product: `promote=True` for a
model-tier primary (`model...` not yet promoted, or `inject...` already promoted — either
way this is the model tier and injection is what it needs), `promote=False` for an
`empirical` primary. `make_psf.py` calls it with `promote=True` explicitly right after
building a model-tier product, so:

- **`data/psf/<...>/psf_kernel.fits` / `<psf cutout dir>/cutout_[cr_]psf.fits`** (the cutout
  dir being `psf_cutout_dir()`'s — the band's, always default-size; see *PSF generation*) — the
  **canonical** files, now the drizzle-broadened injected kernel (`PSFINJ=True` in the
  header) for every model-tier lens.
- **`data/psf/<...>/psf_kernel_analytic.fits` / `cutout_[cr_]psf_analytic.fits`** — the
  pre-broadening analytic model (STDPSF / focus-diverse / MAST-DB) this promotion moved
  aside, kept for comparison. `psf.png`/`psf_analytic.png` mirror the split.
- `info/lens_psf.json`'s `method` for these products is now `inject_acs_fdpsf` /
  `inject_stdpsf` / `inject_wfpc2_psfdb` — it describes what's actually IN
  `cutout_[cr_]psf.fits`, same rule as every other product in this file.
- `info/lens_psf_injected.json` still records every injection run's own metadata
  (`fwhm_pix`, `wing_scatter`, `n_frames`, `null`+exit 0 on no data) independent of
  promotion — an audit trail of what injection produced, kept even though its content now
  usually matches `lens_psf.json` for that product.

The **only** case that still gets the old parallel `*_injected`-suffixed names
(`psf_kernel_injected.fits`, `psf_injected.png`, `cutout_[cr_]psf_injected.fits`, no
promotion) is running injection on an **empirical** primary — `run_psf_inject_all.sh --all`,
purely for the validation comparison above; nothing canonical changes there because the
empirical build is already correct.

Runners: `bash scripts/run_psf_all.sh` already leaves every model-tier lens promoted (the
auto-chain runs inside `make_psf.py`); `bash scripts/run_psf_inject_all.sh` (model tier only
by default, matching `model...` or `inject...` methods; `--all` also runs empirical products
for validation) is for re-promoting after an injection-only code change, or building the rare
lens whose injected product doesn't exist yet, without re-running the analytic model build.
A failed injection (e.g. `data/drizzle_files/` was cleared) degrades inside `make_psf.py` to
the **cheap analytic drizzle-broadening fallback** (`make_psf_inject.analytic_broadened_fallback`
/ `psf_models.analytic_drop_broaden`, 2026-07-31): rather than silently leaving the sharp,
un-broadened analytic model canonical, it convolves that already-North-up analytic kernel
with a box the drizzle drop projects onto the output grid — `pixfrac * native_scale /
out_scale`, further scaled by `1/sqrt(n_frames)` since dithered exposures at different
sub-pixel phases average down a single-frame box's blur (the same "well-dithered" limit
`casertano_r` assumes; `n_frames` from the drizzle header's `NDRIZIM`, halved for 2-chip
ACS/WFC3-UVIS MEFs) — and promotes *that* as canonical instead, same promotion logic as a
real injected build (analytic model moved aside to `*_analytic`). Calibrated against
J0008-0004 F606W, where a real injected build exists to check against: the naive n=1 box
overshoots (4.35→4.67px vs the true re-drizzled 4.42px), the dither-corrected box (n=4)
lands at 4.43px. Distinguished from a real injected build by `PSFINJ=False`/`PSFBROAD=True`
and method `broadened_<...>` — never silently mistaken for the rigorous re-drizzled kernel,
and only used at all when injection itself fails; if the fallback also fails, the analytic
model stays canonical as before.

Promoted 2026-07-30 for all 34 model-tier products then on disk. `psf_epsf.fits` was dropped
from every lens (~150 files) in the same pass — nothing ever read it back from disk; the
archival characterisation now keeps only `psf_kernel.fits`.

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

**The same 3e-3 gate was validated on WFC3/UVIS (gallery, 2026-07-30) and needs no
UVIS-specific retune.** Across 25 gallery products, passing empirical builds top out at 2.26e-3,
with a clean gap to the three F606W builds dropped to the model (J1110+2808 3.1e-3, J0237-0641
3.4e-3, J0918+5104 4.9e-3) — the threshold sits in that gap. Despite UVIS's higher correlated
noise (~1.5–1.6× native ACS), its clean wings are no noisier than ACS/F555W at the gate. Two
marginal drops (J0237-0641, J1110+2808) were contaminant/faint-star problems, not a too-strict
gate: both were rescued to clean empirical builds via `info/psf_stars.json` (leave-one-out
method → memory: uvis_scatter_gate_validated). Gallery is now 19 empirical + 6 injected-model
(J0918+5104 stays model: 4.9e-3, genuinely noisy).

### `run_psf_all.sh --models-only`

`run_psf_all.sh [SAMPLE] [--all|--models-only]`. Default `--all` rebuilds everything. Use
`--models-only` after changing only model-PSF code (e.g. the rotation): it skips any product
whose *existing* `lens_psf.json` method is `empirical` (that code never touches the empirical
builds — cut from the North-up mosaic; the pedestal is a ~1e-4 no-op on ACS) and rebuilds the
model tier plus any not-yet-built product. The skip set is "empirical AND already recorded",
so a new lens still builds.

### Empirical PSFs are built from the no-CR pass, not the CR pass

→ memory: empirical_psf_from_nocr_pass. The LACosmic CR pass flags sharp **field-star
cores** as cosmic rays and masks them in most frames — the drizzled star gets a hole in
its core and many stars fall below the `min_snr`/shape gates and are dropped. The extended
deflector is preserved (why LACosmic beats driz_cr), so the CR *science* image keeps the
true PSF while the CR-*built* star PSF is the corrupted one (J0008 f814W: 3 stars + a core
hole vs 11 clean stars from no-CR). So `make_psf.py` builds the empirical ePSF from the
**least-CR-rejected mosaic on disk**: `star_pass` defaults to a `*_nocrrej_*` pass when
present, else the output pass. It's decoupled from the cutout name — the kernel is still
`cutout_cr_psf.fits` (`PSFPASS=cr, PSFSTARP=nocrrej`) because the point-source PSF is
**pass-independent** (verified: CR/no-CR sky identical, deflector core byte-identical, only
CR-hit pixels differ — so the no-CR PSF is the *correct*, consistent PSF for the CR
science). A per-lens `"psf_star_pass"` override forces a pass. This default is the **standard for
every LACosmic dataset** (the rule is instrument-generic). It helps where the field has
real stars whose cores the CR pass ate: **ACS/WFC** (big win, below). It was **tested on
WFPC2/WF3 F606W (2026-07-29) and does not help** — 0 stars pass the gates on *both* passes
(fields genuinely star-poor, why the MAST-DB tier exists; and with 2–6 frames the no-CR
mosaic is CR-infested), so **F606W stays on the MAST-DB model and gets no no-CR pass**.
(WFPC2 also has no `--no-cr` flag and its skip fires on the cr product, so generating no-CR
there rebuilds *both* passes and loses the `align_wfpc2_to_acs` tie.) **WFC3/UVIS (gallery)**
is a LACosmic dataset that should benefit like ACS, and will the moment `make_psf` gains
UVIS support and no-CR UVIS passes exist — no further code change needed.

A **no-CR pass exists for every ACS product** (generated 2026-07-29 by the safe move-aside
method — `mv` the cr files out, drizzle `--no-cr --nocrrej`, `mv` back; **never** run
`--no-cr` on a lens whose cr you want kept in place, the drizzle script `rmtree`s its output
dir on every non-skipped run). Grids are pixel-identical between passes (crop = dither
footprint), so `psf_stars.json` exclude boxes transfer unchanged. `info/psf_stars.json`
(was empty) now carries per-lens star exclusions/overrides — box order is
`[xmin,xmax,ymin,ymax]`; → memory: psf_stars_exclusion_traps.

**Current state and open limitations.** `run_psf_all.sh`/`make_psf.py` has run across all
three samples; `info/lens_psf.json` holds all products, 0 failures. ACS empirical PSFs were
rebuilt from the no-CR pass (2026-07-29): ACS stars 344→599 (+74%), 4 model→empirical
conversions (J0157, J1023, J1525, J2341 f814W). `slacs_gold` method breakdown (model-tier
methods read `inject_*` post-promotion): **F814W** 37 empirical + 1 `inject_acs_fdpsf`
(J1213+6708, star-poor at 2 stars); **F555W** 16 empirical; **F606W** `inject_wfpc2_psfdb`
(native MAST-DB); **F160W** 3 empirical + 10 exact-filter `inject_stdpsf` (F160W has no CR
pass; the hybrid gate dropped the noisy empirical builds). Each product carries
`pedestal_frac`; the model tier is rotated to North-up. Vet after any no-CR rebuild: more
stars can surface a close double / galaxy (fixed J1451-0239 f814W); high-count builds dilute a
single bad star, low-count (≤~6, e.g. J0029 at 3) ones don't. Gate thresholds (`min_snr=30`,
core ≥15× outskirt, flux floor 5%, `fwhm_tol_hi=1.4`, `star_size=35` for WFPC2,
`pedestal_bad`/`scatter_bad`=3e-3) generalised fine — no retune needed. `slacs_other`/`gallery`
also run (2026-08-01, incl. crowding cut / error-map / cutout-QC; gallery excludes J1110+2808
F814W/F438W): no empirical/model method flips. Open items:
- **PSF uncertainty — coverage is complete; interpretation is the remaining caveat. All 149
  PSF products across all three samples carry an error map** (`psf_kernel_err.fits` +
  `cutout_[cr_]psf_err.fits`, no `null` `err_method`; verified 2026-08-03) — **but the maps do
  not all mean the same thing** (ensemble scatter vs measured model-vs-truth; see the
  injected-kernel bullet below and always check `err_lower_bound` before pooling them in a
  fit). Breakdown: **79 empirical**
  (62 `bootstrap` + 17 `jackknife`), **49 `ensemble_broadened`** (47 WFPC2 MAST-DB + 2 ACS
  focus-diverse — a lower bound), **21 `calibrated_vs_empirical`** (STDPSF, the measured
  model-vs-truth error). The empirical ePSF ships a per-pixel **error
  map** (`make_psf.py`, 2026-07-31): the star sample is resampled and the ePSF rebuilt —
  bootstrap-with-replacement (`--n-boot`, default 100), or leave-one-out jackknife when
  `< JACKKNIFE_MAX_STARS`=6 stars (bootstrap draws degenerate at tiny N) — and the per-pixel
  std of the (unit-sum, co-registered) ensemble is the error map. Written as parallel
  single-HDU files, *not* an ERR extension (keeps every `[0]`/`al.Kernel2D` reader intact):
  `data/psf/.../psf_kernel_err.fits` (full) and `data/cutouts/.../cutout_[cr_]psf_err.fits`
  (trimmed, cropped to the primary kernel's own trim window so it matches `cutout_[cr_]psf.fits`
  pixel-for-pixel). `info/lens_psf.json` empirical entries gain `err_method`, `err_source`
  (`stars`/`exposures`/`db_stars`), `n_boot_valid`, `psf_err_frac` (integrated `sqrt(Σ err²)`)
  and `fwhm_pix_err`. **The scalar tracks star count as intended** — J0008 f814W (11 stars)
  `psf_err_frac`≈0.013, J0029 f814W (3 stars)≈0.103, ~8× larger, so a star-poor build is now
  quantifiably down-weightable. `--no-psf-err` (or `psf_stars.json` `"psf_err": false`) skips
  it for any tier; the point-estimate kernel is untouched either way; `< 2` valid members
  degrades gracefully (no map, null scalars).
  - **Model tier (2026-07-31): ACS focus-diverse and WFPC2 MAST-DB now have error maps too,
    each from its own natural ensemble, not a contrived one.** `psf_models.acs_focus_diverse_psf`
    / `wfpc2_f606w_db_epsf` take `return_ensemble=False` (default, unchanged return — safe
    for any caller that doesn't ask for it) and, when `True`, additionally return
    `(ensemble, method)` in the *exact* convention `make_psf.psf_error_map`/`_fwhm_spread`
    already expect (they're generic over any oversampled-kernel ensemble, not
    empirical-specific — no duplicate statistics code). **ACS**: the per-exposure North-up
    kernels already computed before averaging (one per contributing exposure) are reduced via
    `psf_models._reduce_ensemble` — leave-one-exposure-out jackknife (almost always, since
    SLACS ACS visits run 2–8 exposures, under the same `_JACKKNIFE_MAX_STARS`=6 threshold) or
    bootstrap if ever ≥6. **WFPC2**: bootstrap-with-replacement (n=100) over the ~147 shared
    archival DB stars — a genuinely new, one-time-built ensemble (`psf_models.
    _wfpc2_f606w_db_ensemble`, cached as a 3D FITS cube, `data/reference_files/
    wfpc2_f606w_psfdb/wf3_f606w_epsf_ensemble.fits`, ~12 min to build via ~100 EPSFBuilder
    reruns at ~7s each), shared and re-resampled/rotated per lens exactly like the point
    estimate. Verified: `J1213+6708 f814W` (ACS, 4 exposures, jackknife) `psf_err_frac`
    3.83e-04; all 23 WFPC2 F606W lenses cluster tightly at 0.0555–0.0593 (expected — same
    shared ensemble, only the per-lens North-up rotation differs).
  - **These write to `psf_kernel_analytic_err.fits` / `cutout_[cr_]psf_analytic_err.fits`,
    not the canonical `_err` names** — because a model-tier build immediately auto-chains the
    injection promotion (see *Drizzle-broadened model PSF by injection* below), which replaces
    the canonical kernel with the drizzle-broadened injected one. `make_psf_inject._promote`
    now moves the analytic error files aside alongside the analytic kernel/PNG/cutout it
    already moved, so a stale `_err` file never sits under the canonical name describing a
    kernel it doesn't match. Rolled out sample-wide via `run_psf_all.sh
    slacs_gold --models-only` (2026-07-31): 24/34 model-tier products got an error map (1
    `inject_acs_fdpsf` + 23 `inject_wfpc2_psfdb`); the other 10 are STDPSF, correctly `null`.
  - **The injected (canonical) kernel now has its own error map too** —
    `scripts/make_psf_err_injected.py`, 2026-08-02, closing what this file previously called
    a deliberate gap. All 70 injected products carry `psf_kernel_err.fits` +
    `cutout_[cr_]psf_err.fits` describing the *canonical* kernel (same single-HDU,
    non-renormalised convention as the empirical tier). **Two sources, never
    interchangeable — the header `PSFERR`/`PSFEBND` and the JSON `err_method`/
    `err_lower_bound` always say which:**
    - **`ensemble_broadened` (49: 47 WFPC2 MAST-DB + 2 ACS focus-diverse) — a LOWER BOUND,
      `PSFEBND=True`.** The analytic ensemble map convolved with the same drop box the
      kernel got (`psf_models.drop_convolve_box`; the box preserves the sum, so no
      renormalisation and the map stays in kernel amplitude units). Assumes the ensemble
      perturbations are correlated on scales ≥ the box, so σ transforms like the signal
      rather than in quadrature — true for smooth focus/star-sampling perturbations. It is
      a lower bound because a bootstrap over N estimates the standard error of the *mean*
      ePSF, which shrinks as 1/√N, while the model's error against a specific observation
      is set by focus/breathing mismatch and does not shrink at all. **This is the only
      thing measurable for WFPC2 F606W** — its fields are too star-poor for any empirical
      build to exist to compare against. Measured 0.053–0.059 for the 47 WFPC2 products;
      the 2 ACS focus-diverse ones sit ~2 orders of magnitude lower (J1213+6708 3.83e-04,
      J1016+3859 5.71e-04) because a jackknife over 4–8 exposures measures how stable the
      exposure-average is, not how far the focus-diverse model sits from truth — a
      **very** loose lower bound there, and one whose tier exists precisely because those
      fields had too few stars for an empirical build to check against.
    - **`calibrated_vs_empirical` (21 STDPSF: 14 F160W + 7 gallery UVIS) — the measured
      model-vs-truth error**, the quantity the ensemble route cannot see. Built from
      `scripts/psf_model_error*.py` (*PSF model-error calibration*, below): radial *shape*
      from the measured residual
      profile, *magnitude* from the group's median quadrature-corrected residual. A flat map
      would be badly wrong — the residual is ~95% core-dominated while the wings carry a
      6–48% flux deficit, so `psf_err_frac_core`/`psf_err_frac_wing` are recorded alongside
      the total. Necessarily a **group-level transfer by filter**, since a model-tier lens
      has no empirical build of its own (that is why it fell back to the model): F160W
      0.0358 (n=5, real injected comparisons), gallery F606W 0.0370 (n=11), F814W 0.0242
      (n=4), F438W 0.0188 (n=3) — the gallery three from the analytic stand-in, so group
      medians only. Empirical products are never touched (only `method` starting `inject`
      is considered); verified 79 empirical maps unchanged. → memory:
      stdpsf_model_error_measured
  - **STDPSF has no *ensemble* uncertainty and never will** — a single static detector-frame
    grid with no per-lens ensemble to resample (no per-exposure retrieval, no per-star
    archive); it would need a focus-perturbation grid STScI doesn't publish per-filter. That
    is why its error map is the calibrated-vs-empirical one above rather than a contrived
    ensemble. → memory: psf_uncertainty_empirical, stdpsf_model_error_measured
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

### PSF model-error calibration (`psf_model_error.py`, `psf_model_error_injected.py`)

```bash
uv run python scripts/psf_model_error.py                       # analytic stand-in, all tiers
uv run python scripts/psf_model_error_injected.py --filt f160W # real injected builds
uv run python scripts/make_psf_err_injected.py --dry-run       # consume it; writes nothing
```

Read-only QC that measures the **model tier against empirical drizzled truth** — the
model-vs-accuracy number an ensemble bootstrap structurally cannot produce (it estimates the
standard error of the *mean* ePSF, which shrinks as 1/√N; the error against a specific
observation is set by focus/breathing mismatch and doesn't shrink at all). Targets every
empirical build whose model tier would be STDPSF: F160W and all of `gallery`. Metric of
record is `sqrt(Σ(mod−emp)²)` on aligned unit-sum kernels — the same construction as
`psf_err_frac`, so the two are directly comparable. Measured 2026-08-02: **median 2.0× the
bootstrap scatter on well-constrained truth (n_stars≥6), up to 7×**; F160W 0.036 vs bootstrap
0.0086 (~4×). This is what `make_psf_err_injected.py` turns into the STDPSF tier's error map.

**Three traps, all of which silently corrupt the measurement:**
- **`make_psf.measure_fwhm` is unreliable on sharp, near-critically-sampled kernels** — its
  2D-Gaussian fit returns a *larger* FWHM for a *higher*-peak kernel on the unbroadened
  analytic model. Use `psf_model_error.fwhm_radial` (half-max crossing on the radial profile)
  there. **Not a live pipeline bug**: the pipeline only ever fits the broader
  drizzled/injected kernels, and every recorded `fwhm_pix` was checked and is sane.
- **Residuals must be sub-pixel aligned first** — raw model-vs-empirical centroid offsets
  reach ~1.9px and swamp the metric, turning a shape measurement into a centring one.
- **The L2 residual is 94–97% core-dominated**, so it barely sees the wing deficit (model
  enclosed wing flux runs 6–48% low). Report/propagate core and wing separately; one scalar
  hides the part that matters most for extended arc flux.

`psf_model_error.py` stands in for the injection promotion with
`psf_models.analytic_drop_broaden`, which is **right in the median but unreliable per-lens**
(individual F160W lenses moved up to 2.7× in both directions when checked against real
injection, while the group median moved only 0.032→0.037) — so quote its **group medians
only**. `psf_model_error_injected.py` does the honest per-lens comparison wherever a real
`psf_kernel_injected.fits` exists. Building those is non-destructive on an empirical primary:
`run_injection(promote=None)` auto-detects `promote=False` and writes only the parallel
`*_injected` names (verified by mtime). → memory: stdpsf_model_error_measured

## Tracking JSONs in `info/`

> **Full reset, 2026-07-26.** All three files were emptied to `{}` and every product under
> `data/` deleted (`calibrated/`, `drizzle_files/`, `drizzled/`, `cutouts/`, `mosaics/`,
> `run_logs/`) as a deliberate clean restart. Kept: `data/reference_files/` (CRDS cache) and
> `data/pre_drizzled/` (46 MAST-delivered mosaics, not pipeline output; **that folder was
> deleted 2026-09-14** — it had been empty for some time, so any claim elsewhere that the
> repo carries MAST-delivered mosaics is void). The sample was
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

- **`lens_crfill.json`** — `{sample: {lens: {band: {...}}}}` — **the list of bands whose
  science stamp is a STREAK-REPAIRED (`--crfill`) product**, written by hand when one is
  promoted, not by a run. Two entries so far (J1213+6708 and J1250+0523 f814W, 2026-09-22).
  Each records the drizzle and cutout commands verbatim, the selection flags, every track
  kept (frame, chip, area, closest approach), pixels filled against pixels flagged, the
  noise ridge before and after, the √(N/(N−1)) noise caveat, whether the astrometry moved
  and what was shifted downstream. **Check it before trusting a per-pixel noise value on
  those bands** — the stamp's own `CRFILL=T` card says *that* it is filled, this file says
  *where* and *how much*. Keep it in step by hand: nothing validates it against the headers.

Not written by a pipeline run, and the odd one out in `info/`:
- **`lens_einstein_radii.json`** — `{lens: {theta_e_arcsec, RE_kpc, zlens, zsrc, sigma_kms,
  source}}` — **flat by lens, not nested by sample**, because it is a *catalogue* and not a
  record of this repo's products: a lens lands here whatever sample it sits in. **89 rows from
  two papers, and the `source` field is the only thing that says which** — read it before
  pooling or comparing them:
  - **74 SLACS** from Auger+2009 (SLACS IX) via VizieR, `RE` kpc → arcsec on flat
    H0=70/**Ωm=0.3**. Refreshed on demand (`detect_arcs.py --fetch-theta-e`).
  - **15 gallery** from **Shu+2016 (BELLS GALLERY IV, ApJ 833, 264) Table 2** `bSIE`
    (Sérsic foreground subtraction), z/σ from its Table 1, **entered by hand 2026-09-17** —
    that paper is **not on VizieR** (only Shu+2016a, the parent candidate list, is), so
    `--fetch-theta-e` can never supply them. `RE_kpc` is *derived* here, on Shu's own flat
    H0=70/**Ωm=0.274**, the reverse of the SLACS rows where kpc is native and arcsec derived.
    The join was made on the full SDSS names already in `info/gallery_coords.py`, not by hand
    — Table 2 also lists **J0918+4518**, one digit from gallery's **J0918+5104**, and both are
    grade-A. Values validated end-to-end: every detected gallery component lands at
    0.64–1.37 θ_E, and the dashed θ_E circle traces the arcs on all 6 contact-sheet panels.
  - **`--fetch-theta-e` MERGES, and must keep doing so.** It used to write the fetched table
    wholesale, which would now silently delete all 15 gallery rows and take
    `--propose-from detect` down on that whole sample — the failure would surface only as
    "no measured Einstein radius" much later. Rows the catalogue does not carry are kept; a
    fetched row wins over a hand-entered one of the same name.

  Tracked in git so the detector runs offline. A lens with no row is one neither paper
  modelled: **`J1259+6134`, `J2141-0001`, `J2302-0840`** in `slacs_other` (of which only the
  last two have cutouts at all).

No data for a filter → value `null`.

- **The key is the product directory, not the filter.** Usually they coincide (`f606W`), but
  a split-visit lens is keyed per visit: `f606W` for its primary (longer-exptime) visit,
  `f606W_v2` for a shorter one — so a bare `f606W` key does not by itself mean "the whole
  filter, unsplit" for J0728+3835/J0822+2652 the way it does for every other WFPC2 lens.
  Check keys against product directories, not against the filter name, when auditing these
  two. This once caught a live error (pre-2026-08-03 naming): both split lenses had recorded
  a plausible-but-nonexistent combined `f606W` (6 obsids, 6600 s) under the visit-order
  `_v1`/`_v2` scheme, invisible precisely because a `f606W` key was present but meant
  something else. Root cause was JSON writes keyed on the bare filter while writing to a
  `--out-suffix` dir; keyed on `product_key = filt + out_suffix` since.
- **Records drizzled frames, not the download.** WFPC2 `--pa` selects one visit; a lens that
  exits for want of dither phase writes nothing. ACS/WFC3 silently drop sub-`MIN_EXPTIME`
  frames (same 10s floor as WFPC2 since 2026-08-03 — see *Lens Samples* below), so those are
  excluded from the record and from the drizzle work directory. Exposure times were always
  correct (from the header).
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
uv run python scripts/mast_target_names.py --list        # samples + sizes
uv run python scripts/mast_target_names.py slacs_gold    # lens names
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
**F814W** (ACS/WFC, 4/27 — the rest are `BLOCK_EXPTIME`-gated or absent, see *Total-
exposure-time gate*; `J1016+3859` recovered 2026-08-03 via `force_copy`, see above),
**F160W** (WFC3/IR, 6/27), **F555W** (ACS/WFC, 0/27 — none of these lenses fall in props
10494/10798). 4 lenses (J0959+4416, J1016+3859, J1153+4612, J1416+5136) have both F606W and
F814W; `align_wfpc2_to_acs.py` has tied all 4 (`J1016+3859` on 2026-08-04, and it also needed
its *F814W* tied to F160W first — see the `--target`/`--ref` note in *align_wfpc2_to_acs.py*)
— the other 21 F606W products carry their delivered GSC240 absolute WCS (~0.3–1″ off) untied. `info/wfpc2_alignment.json` has no
per-lens `--align` audit for `slacs_other` (it only covers the 22 `slacs_gold` WFPC2 lenses);
every `slacs_other` WFPC2 lens falls back to the documented default, `mast`.
`run_psf_all.sh` has run for this sample (2026-08-01): 34/34 PSF products,
`info/lens_psf.json` populated (24 `inject_wfpc2_psfdb` at F606W; empirical +
`inject_stdpsf`-family across F814W/F160W, incl. J1016+3859's star-poor `inject_acs_fdpsf`) —
see *PSF generation* above.

`gallery` coverage (15 lenses, reduced 2026-07-29): **F606W** on all 15 (the primary band),
**F814W**/**F438W** on 6, **F275W** on 5, **F225W** on 1 (J2342-0120). See *BELLS GALLERY:
WFC3/UVIS reduction* below for the pipeline and its caveats. Gallery lenses are **not in
`info/slacs_coords.py`** — they use `info/gallery_coords.py` instead, read by
`drizzle_wfc3_uvis.py` for the common output WCS; a lens absent from that table falls back
to native drizzle WCS with a warning. PSF products exist (WFC3/UVIS is keyed in
`make_psf.py`/`psf_models.py`; `run_psf_all.sh gallery` run 2026-08-01, 25/25 — see *PSF
generation*). **F225W/F275W are confirmed
unusable for lens science across the whole sample** (arc undetected, not just the deflector
— see *BELLS GALLERY* below); **J1110+2808's F275W is pure noise, but its F814W/F438W DO show the lensed
images** — the older "F606W only" claim was wrong and is corrected below. No further reduction
effort on F225W/F275W anywhere — reduced correctly, just not lensing-useful.

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
- **ACS/WFC3-IR/WFC3-UVIS** filter COPY out in favour of non-COPY, **except** lenses with
  `"force_copy": true` in `lens_samples.json`, needed whenever the preferred non-COPY visit
  turns out to be entirely dead frames: `J1032+5322` F814W (non-COPY frames `EXPTIME=0`) and
  `J1016+3859` F814W (non-COPY frames `EXPTIME=0`, recovered 2026-08-03 — see the
  per-exposure floor note below). `drizzle_wfc3_ir.py` gained `force_copy` support in the
  same pass, so a WFC3/IR lens hitting this can now use it too.
- **WFPC2** keeps both (COPY sets are genuine repeat visits) and rejects junk on
  `MIN_EXPTIME`.
- **Per-exposure exptime floor (2026-08-03).** The ACS/WFC3-IR/WFC3-UVIS per-frame filter
  was `EXPTIME > 0`, which only catches an exact-zero frame — a dead/aborted exposure can
  carry a small nonzero `EXPTIME` (0.5s seen on real archive data) that used to slip through
  into the drizzle work directory at a small nonzero weight. Raised to `MIN_EXPTIME = 10s` in
  all three scripts, matching WFPC2's existing constant; the copy-to-work-dir step now
  excludes sub-floor frames directly, not just the provenance record.

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
`final_rot=0.0` rotates into 4 parallel diagonal stripes in the noise map.

## BELLS GALLERY: WFC3/UVIS reduction (`scripts/drizzle_wfc3_uvis.py`)

The `gallery` sample (15 lenses, Shu et al. 2016 Table 1 class E-S-A) is BELLS GALLERY
(props 14189, 16734), imaged in WFC3/UVIS across five filters: **F225W, F275W, F438W,
F606W, F814W**. One filter-agnostic script (`--filt`), modeled closely on
`drizzle_acs_wfc.py` — UVIS is a two-CCD optical detector like ACS/WFC, so it inherits the
same alignment/CR reasoning. **That inheritance was independently audited on gallery data
2026-09-13 and both halves hold** (see *Gallery audit* below); the per-lens `--align`
override table in the script is still empty, which the audit now justifies rather than
merely reflects (an explicit `--align` always wins).

```bash
uv run python scripts/drizzle_wfc3_uvis.py --lens J1110+3649 --filt f606W
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
- **F225W/F275W are unusable for lens science, sample-wide — confirmed, not just a faint
  deflector.** A 2026-07-29 same-stretch S/N check found the *arc* undetected too, not only the
  UV-dark deflector: J1110+2808 F275W (15768s, the deepest UV exposure) is pure noise at the
  lens; J0742+3341's deflector+ring sits at S/N~2.5 (vs ~14 in F606W); J2342-0120 F225W shows
  nothing at the lens centroid. **Do not spend further effort on F225W/F275W** — no PSF work,
  no recentring/alignment tuning. Products stay on disk as a correctly-reduced record, not
  science-ready cutouts. → memory: gallery_uv_bands_unusable
  - The centring mechanism that motivated the original caveat is still real if these are
    ever re-cut: **never use `--center-self`** on F225W/F275W — a self-centred peak search
    locks onto noise or a field source, not the lens. Cut with `--center-band f814W`
    (default) so the stamp geometry is at least correct even though nothing but noise
    should be expected there.
- **PSF support: WFC3/UVIS is wired in** (`make_psf.py`/`psf_models.py` key off `ACS/WFC`,
  `WFC3/IR`, `WFPC2`, and `WFC3/UVIS`), and `run_psf_all.sh gallery` has been run (2026-08-01,
  excluding J1110+2808's F814W/F438W — **that exclusion rests on a claim since retracted; those
  two PSFs need building if the bands are used**, see the per-lens bullet): 25/25 products ok, 18 empirical +
  7 `inject_stdpsf` (no exact-filter WFC3/UVIS STDPSF grid substitution issue — UVIS uses the
  same ACS/WFC3 STDPSF machinery). See *PSF generation* above (F160W hybrid quality gate /
  `uvis_scatter_gate_validated`) for the empirical/model split rationale.
- **Per-lens caveat: J1110+2808 — F275W is noise, but F814W and F438W DO detect the lensed
  images. The 2026-07-29 "F606W only" verdict was wrong and is RETRACTED (re-measured
  2026-09-14).** Its F814W/F438W run ~2× the exposure of every other gallery lens (and F275W
  15768s, the deepest). The system is a double (published theta_E 0.98", images predicted at
  0.73"/1.28"); the observed pair sits at 0.73" (col 142, row 136) and 1.23" (col 172, row 175).
  Measured on the **elliptical-profile-subtracted** stamp, per-pixel S/N peaks and the count of
  pixels above S/N 4 are: F606W 11.3 (23 px) and 13.8 (32 px); **F814W 5.5 (4 px) and 8.6
  (12 px)**; **F438W 6.5 (6 px)** for the inner image, 2.9 (0 px) for the outer; F275W 2.1/2.0
  (nothing). Against control apertures placed at the same radius with all real structure
  excluded, the F814W detections are +6.3 and +7.4 sigma above the local systematic floor and
  the F438W inner image +11.3. **Why the original check missed them:** it was a same-stretch
  *visual* S/N comparison, and both images sit on the deflector envelope, where they are
  invisible until the galaxy is subtracted. **Consequence:** F814W is usable for this lens
  (positions, colours); its F814W/F438W PSFs were deliberately skipped and must be built before
  the band is modelled.

### Published lens models — every gallery lens has one (Shu et al. 2016, Table 2)

**Before deciding whether a blob is a lensed image or a contaminant, solve the published
model.** All 15 gallery lenses are grade-A in Shu et al. 2016, ApJ 833, 264 (*BELLS IV:
Smooth Lens Models*), which fits an SIE (+shear where needed) to this same F606W imaging.
**It is not on VizieR** — only Shu+2016a (`J/ApJ/824/86/table2`, the parent candidate list)
is there. Get it from the paper:

```bash
curl -sL https://arxiv.org/pdf/1608.08707 -o /tmp/shu.pdf && pdftotext -layout /tmp/shu.pdf
```

Table 2 = lens parameters (`bSIE` = θ_E, q, PA east of north, centroid offset, shear γ/φ_γ,
mean magnification µ); the second table labelled "TABLE 2" (source parameters, referenced in
the text as Table 4) = per-component source offsets, q, n, R_eff, m_AB; Table 3 = per-lens
notes (perturbers, quads, which systems need shear). θ_E in arcsec:

| J0029+2544 | J0201+3228 | J0237-0641 | J0742+3341 | J0755+3445 | J0856+2010 | J0918+5104 | J1110+2808 |
|---|---|---|---|---|---|---|---|
| 1.34 | 1.70 | 0.65 | 1.22 | 2.05 (γ=0.24) | 0.98 | 1.60 (γ=0.18) | 0.98 |

| J1110+3649 | J1116+0915 | J1141+2216 | J1201+4743 | J1226+5457 | J2228+1205 | J2342-0120 |
|---|---|---|---|---|---|---|
| 1.16 | 1.03 | 1.27 | 1.18 | 1.37 (γ=0.15) | 1.28 | 1.11 |

**Use the predicted image RADII and the image COUNT; do not trust the predicted PAs.** Feed
each source component through the SIE+shear lens equation and keep the radii — those are
convention-free and they are what settles "is this thing too far out to be an image". The
paper's ΔR.A. sign convention did **not** reproduce observed configurations consistently
across lenses here, so a predicted position angle can come out mirrored. Validated against
the data on three lenses: J0237-0641 predicts two images at 0.61/0.79″ and the stamp shows a
pair at 0.56/0.74″ (166° apart, mean 0.65″ = θ_E exactly); J1110+2808 predicts 0.73/1.28″ and
shows 0.73/1.23″; J0918+5104 predicts a quad spanning 1.16–2.16″ and its outermost observed
image sits at 2.25″. **The paper's own figures are a second, independent check**: pages 8–12
show data / lens-light-subtracted / model / **residual** per lens, and an object the authors
did not treat as lensed is left standing in the residual panel — that is how J1116+0915's
1.35″ object and J2342-0120's 1.88″ object were confirmed as contaminants. Render one with
`pdftoppm -png -r 400 -f <page> -l <page> /tmp/shu.pdf out`.

### Gallery audit (2026-09-13) — the inherited ACS reasoning, checked

Read-only; nothing was re-drizzled or rewritten.

**Alignment: no cross-band tie is needed — verified, not assumed.** Field sources
cross-matched between F606W and F814W on all 6 lenses that have both (8–140 matches each)
give median VECTOR offsets of **≤0.04″** (largest component 0.039″). So a lens's bands do
co-register and `align_wfpc2_to_acs.py` has nothing to fix here, exactly as claimed above.
WCSNAME does vary per lens and per filter as documented, but unlike ACS F555W it does **not**
predict a real offset.

- **TRAP: do not point `align_wfpc2_to_acs.py` at gallery.** Its deflector-centroid test —
  the measurement that correctly drove the F555W tie in 882b813 — reports **0.16–0.63″** on
  these same well-aligned lenses. The centroid genuinely moves between bands because
  gallery's lensed arc is blue and bright in F606W while the deflector is red, so the light
  centroid is not a fixed point across filters. Acting on that number would shift CRVAL by up
  to 0.6″ and **break** astrometry that is currently good. The SLACS F555W case was safe from
  this only because its arcs are faint in both ACS bands. **Cross-band astrometry on gallery
  must use field sources, not the deflector.**

**CR rejection: the LACosmic default stands; no per-lens tuning is warranted.** Core flux
(r<1″) is preserved to a median **0.985** of the no-CR pass across F606W/F814W. The arc
annulus (1–2″) sits at 0.932 vs 0.989 for a slacs_gold F814W control, which looked like
erosion and is **fully explained by the arcs being fainter** (see next bullet) — a
`scripts/lacosmic_erosion_scan.py` sweep over
sigclip 4.5–10 / objlim 5–12 on the two worst lenses (J2228+1205, J0237-0641) shows their
CR-flags-on-real-signal fraction barely responds to the thresholds (−13%, −10%), where the
documented ACS erosion case J1420+6019 halves (−54%). **The best-behaved gallery lens scores
the *highest* absolute fraction**, so the quantity does not track arc loss at all. Read that
script's docstring before repeating this: the intuitive "flagged in every frame" metric
returns zero even on the known-eroded product, and the absolute fraction is not comparable
across instruments — only the *response* to the thresholds is.

- **Why the arc annulus reads 0.932 on gallery and 0.989 on SLACS — same CR flux, fainter
  arcs, nothing eroded.** The decisive control is a **blank-sky annulus (8–10″)** in the same
  images, where there is no source to shave. Per lens, *(CR flux removed at the arc) / (CR
  flux removed in blank sky)* comes out at **gallery median 1.01** (range 0.70–2.04, n=15)
  and **slacs median 1.11** (0.64–1.76, n=6): the source region loses exactly what empty sky
  loses, in both samples, and gallery is if anything the cleaner of the two. So the numerator
  is cosmic rays, uniformly distributed. What differs is the **denominator** — annulus surface
  brightness is **5.8 e/s/arcsec² on gallery vs 38.5 on slacs_gold, 6.7× fainter** — while the
  CR flux removed is essentially the same (0.38 vs 0.48 e/s/arcsec² in blank sky). Same
  subtraction, much smaller thing to subtract it from. Independent confirmation from frame
  count: a CR is diluted 1/N by the drizzle average, and the two NDRIZIM=24 gallery lenses
  lose **0.154** e/s/arcsec² against **0.381** for the NDRIZIM=8 majority — a factor 2.5
  against the 3.0 predicted. **The lesson for any future check: a CR/no-CR flux RATIO is not
  a measure of erosion, because it divides by source brightness.** Compare the absolute flux
  removed against blank sky in the same image instead.
- Two measurement traps this audit walked into, both worth avoiding: **selecting pixels on
  the no-CR pass and then measuring CR/no-CR biases the ratio low** (F438W read 0.29–0.76
  that way and 0.96–1.01 selecting the other way — the truth is ~1); and **ratios of
  background-dominated sums are meaningless**, which is what made the first F438W and
  arc-annulus numbers look alarming.

**F438W detects the ARC but not the deflector** — sharpen "science-ready" accordingly. Peak
S/N within 0.5″ of the lens centre is **2.2–3.5** on all six lenses, indistinguishable from
the unusable UV bands; but over the inner 3″ it reaches **S/N 5.0–19.7**. The band is doing
exactly what a blue filter should on a red lens with a blue source: the deflector is absent
and the lensed arc is not. So F438W is usable for source/arc work and useless for deflector
light or for anything keyed on the deflector centroid — including the tie check above.

Current state (reduced 2026-07-29): all 15 lenses have F606W; F814W/F438W on 6 each,
F275W on 5, F225W on 1 (J2342-0120) — matches the sparse per-lens filter coverage BELLS
GALLERY actually has on MAST, not a pipeline gap. `run_cutouts_all.sh` was extended to glob
the UV/blue bands (`f438W f275W f225W`) alongside the SLACS filters so it stays one runner
for every sample. **F225W/F275W across every lens are not lensing-useful** (see bullets
above) — treat F606W (all 15 lenses) and F814W/F438W (**all 6** of the lenses that have them,
J1110+2808 included — see its corrected per-lens bullet above) as science-ready, with the
F438W qualification from the audit above: **arc yes, deflector no.**

## QC mosaics (`scripts/make_mosaics.py`, `scripts/make_psf_mosaics.py`)

```bash
uv run python scripts/make_mosaics.py --sample slacs_gold
uv run python scripts/make_mosaics.py --sample slacs_gold --size 20   # the 20" stamp tree
uv run python scripts/make_psf_mosaics.py --sample slacs_gold
```

Read-only QC: tile every lens's existing cutout (or trimmed PSF kernel) onto a 5-wide grid,
one mosaic per (filter group, panel type), written to `data/mosaics/<sample>/`. Nothing is
re-drizzled, re-cut, or rebuilt — pure visualization over what's already on disk in
`data/cutouts/`.

- **Filter groups** come from `scripts/mosaic_groups.py`, shared by both scripts so they
  stay in sync. For `slacs_gold`/`slacs_other`: `f814W`, `f606W_f555W` (WFPC2 F606W —
  including the split-visit `f606W_v2` key, tried only if the lens has no bare `f606W` —
  merged with ACS F555W per lens, since no SLACS lens has both), `f160W`. For `gallery`: one
  group per UVIS filter (no
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
  much as sharpness), on a **log stretch over 1e-4..1** (`pooled_log_norm`, matching the
  'log' wing panel in each lens's `psf.png` — the wings are what QC here is about; passed to
  `make_mosaics.plot_mosaic` via its `norm_fn` hook, which still defaults to pooled asinh for
  the signal/noise mosaics). Each panel tagged `emp` (empirical ePSF) or `mod` (STDPSF /
  focus-diverse / MAST PSF DB) from the `PSFMETH` keyword. All three samples have PSF
  products and `_psf`-panel mosaics now; regenerate both mosaic scripts after any PSF
  campaign (fast, read-only — a few seconds per sample) so the QC PNGs match the kernels
  currently on disk.

## NICMOS is deprioritised

Do not generate or propose NICMOS (NIC2) products unless explicitly asked — the FOV is far
too small (~19″ vs ~139″ for WFC3/IR) and the pipeline may be unsound. Answer F160W coverage
from WFC3/IR. All NICMOS data was deleted 2026-07-21 (472 MB); re-runnable via
`scripts/stale_scripts/drizzle_nic2.py` (raises `NotImplementedError` on import unless
`ALLOW_NICMOS=1`), which re-downloads from MAST.
