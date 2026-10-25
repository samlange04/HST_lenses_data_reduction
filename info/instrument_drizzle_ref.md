# HST Drizzle & PSF Reference

Condensed from the STScI handbooks and PSF papers behind the pipeline's DQ-bit, pixel-scale,
read-noise and PSF choices, with what each driver actually does today. Page citations are
written as (DHB p###) and refer to the **PDF page** of the edition listed here. The PDFs are
not part of the repo; download them from the links if you want to check a citation. The
handbooks are living documents, and a newer edition may renumber pages.

| cited as | document | download |
|---|---|---|
| `acs_dhb.pdf` | ACS Data Handbook v14.0 (Cycle 34, 2026), 240 pp. | https://hst-docs.stsci.edu/acsdhb |
| `wfc3dhb2024_final.pdf` | WFC3 Data Handbook v6.0 (June 2024), 296 pp. | https://hst-docs.stsci.edu/wfc3dhb |
| `The_DrizzlePac_Handbook_Version3.pdf` | The DrizzlePac Handbook v3.0, 154 pp. | https://hst-docs.stsci.edu/drizzpac |
| `wfpc2_dhb.pdf` | WFPC2 Data Handbook v5.0 (July 2010), 176 pp. | https://www.stsci.edu/files/live/sites/www/files/home/hst/instrumentation/legacy/wfpc2/_documents/wfpc2_dhb.pdf |
| `wfpc2_ihb.pdf` | WFPC2 Instrument Handbook v10.0 (August 2008), 338 pp.; Table 4.2 read noise and gain (PDF p91) | https://www.stsci.edu/files/live/sites/www/files/home/hst/documentation/_documents/wfpc2/wfpc2_ihb_cycle17.pdf |
| ACS ISR 2017-02 | Hoffmann & Avila 2017 (Nov 2018 update, Table A), saturation bit 256 in the MDRIZTAB | https://www.stsci.edu/hst/instrumentation/acs/documentation/instrument-science-reports-isrs |
| ACS ISR 2017-05 | Borncamp et al. 2017, unstable-pixel bit 32 | same ISR index |
| WFC3 ISR 2016-12 | Anderson 2016, "Empirical Models for the WFC3/IR PSF" | https://www.stsci.edu/hst/instrumentation/wfc3/documentation/instrument-science-reports-isrs |
| WFC3 ISR 2021-12 | Dauphin et al. 2021, "The WFPC2 and WFC3 PSF Database" | same ISR index |
| Anderson & King 2000 | "Toward High-Precision Astrometry with WFPC2. I. Deriving an Accurate PSF", PASP 112, 1360 | https://arxiv.org/abs/astro-ph/0006325 |

The NICMOS Data Handbook was consulted for the retired `stale_scripts/drizzle_nic2.py` and is
no longer cited by anything in use.

---

## Quick-reference table: what the four drivers do

| Driver | Instrument | Input | Pipeline / ref env | Native → output scale | pixfrac | DQ bits kept as good | Weights | CR rejection | Suffix |
|---|---|---|---|---|---|---|---|---|---|
| `drizzle_acs_wfc.py` | ACS/WFC | `_flc` | CALACS / `jref` | 0.05 → 0.05″ | 1.0 | `256,64,16` | ERR | LACosmic per frame, then plain-mean drizzle | `_drc_` |
| `drizzle_wfc3_uvis.py` | WFC3/UVIS | `_flc` | CALWF3 / `iref` | 0.0396 → 0.0396″ | 1.0 | `256,64,16` | ERR | LACosmic, as ACS | `_drc_` |
| `drizzle_wfc3_ir.py` | WFC3/IR | `_flt` | CALWF3 / `iref` | 0.1283 → 0.06″ | 1.0 | `512` | ERR (K = EXPTIME) | **none** (up-the-ramp in calwf3) | `_drz_` |
| `drizzle_wfpc2_wf3.py` | WFPC2/WF3 | `_flt` (WF3 chip only) | CALWFPC2 / `uref` | 0.0996 → 0.05″ | 1.0 | `8,1024` | IVM, built per frame | LACosmic, as ACS | `_drw_` |

Common to all four: `--align mast` (trust the delivered WCS; TweakReg is opt-in, used only
by J1213+6708 f814W), `final_rot=0` and `final_ra/dec` at the catalogue position, `num_cores=1`, the
exposure-time gate (block below 500 s, warn below 1200 s) and a 10 s per-frame floor. The
reasoning behind each choice is in `AGENTS.md`; this file records the handbook facts it rests
on. The output suffix follows the input file type (FLC → `_drc_`, WFPC2 FLT → `_drw_`, else
`_drz_`).

AstroDrizzle's own cosmic-ray route (`driz_cr` against a blotted median) is kept as
`--cr-method drizcr` for comparison runs only. Its tuning values, where the drivers still
carry them, are the handbook defaults: `driz_cr_snr='3.5 3.0'`, `driz_cr_scale='1.2 0.7'`
for ACS, UVIS and IR; `driz_cr_snr='15.0 10.0'`, `driz_cr_scale='1.5 1.0'` for WFPC2.

---

## ACS/WFC (`drizzle_acs_wfc.py`)
**Source:** ACS Data Handbook (`acs_dhb.pdf`)

| Parameter | Value | Notes |
|---|---|---|
| File type | `_flc.fits` | CTE-corrected; FLT also exists but FLC is used |
| MAST project filter | `CALACS` | `HAP-MVM`, `HAP-SVM` products also appear and are filtered out |
| Pixel scale | 0.05″/px | Native ~0.049 (ACS DHB p12); IDCTAB and pipeline use 0.05 (drizzle headers, p40–42) |
| `driz_sep_bits` / `final_bits` | `'256,64,16'` = 336 | STScI MDRIZTAB default since 2017 (DrizzlePac DHB p132; ACS ISR 2017-02, 2018 update Table A; ACS ISR 2017-05 for bit 32) |
| Full well | 84 700 e⁻ | Used by the PSF star saturation gate (`make_psf._SATLEVEL`) |
| TweakReg `conv_width` / `ylimit` | 3.5 / 0.2 | Only under `--align tweakreg` |

**Why bit 256 (saturation) is included:** Hoffmann & Avila 2017 (ACS ISR 2017-02) added 256
to `driz_sep_bits`/`final_bits` because AstroDrizzle was masking full-well-saturated pixels
and losing flux: "charge is not lost but moves along the columns for saturated pixels". The
same ISR raised `driz_cr_scale` to 1.5 1.2 for 2- and 4-image associations. Borncamp et al.
2017 (ACS ISR 2017-05) then redefined bit 32 as *unstable* pixel and swapped 32 → 16 in the
MDRIZTAB, giving 16 + 64 + 256 = 336 (ACS ISR 2017-02, November 2018 update, Table A).
Stable hot pixels (16) are dark-subtracted and kept; unstable ones (32) stay excluded.

**DQ flag table (Table 3.4, ACS DHB p87):**

| Flag | Meaning | Kept? |
|---|---|---|
| 4 | Bad/vignetted pixel | — (filled by `--bcfill` on the `bcfill` branch) |
| 8 | Masked by aperture feature | — |
| 16 | Hot pixel (dark > 0.14 e⁻/s), **stable** | ✓ |
| 32 | Unstable dark / RTS noise | — |
| 64 | Warm pixel (dark 0.06–0.14 e⁻/s) | ✓ |
| 128 | Bad column | — (filled by `--bcfill` on the `bcfill` branch) |
| 256 | Full-well saturation | ✓ (MDRIZTAB 2017) |
| 512 | Bad pixel in reference file | — |
| 1024 | Sink pixel / charge trap | — |
| 2048 | A-to-D saturation | — |
| 4096 | Cosmic ray (written by LACosmic here; AstroDrizzle's bit) | — |
| 8192 | acsrej CR | — |

---

## WFC3/UVIS (`drizzle_wfc3_uvis.py`)
**Source:** WFC3 Data Handbook 2024 (`wfc3dhb2024_final.pdf`)

| Parameter | Value | Notes |
|---|---|---|
| File type | `_flc.fits` | Pixel-based CTE correction applied; like ACS, not like WFC3/IR |
| MAST project filter | `CALWF3` | |
| Pixel scale | 0.0396″/px, kept native | pixfrac 1.0 (a 0.7 scan bought nothing; `AGENTS.md` *BELLS GALLERY*) |
| `driz_sep_bits` / `final_bits` | `'256,64,16'` | Same meanings as ACS/WFC (WFC3 DHB Table 3.3, UVIS column). STScI's UVIS examples use 80 = 16 + 64, dropping saturation; 256 is kept here as for ACS |
| Full well | ~70 000 e⁻ | WFC3 DHB 2024 §1.2; used by the PSF star saturation gate |
| ERR units | electrons | So `1/sqrt(WHT)` is a σ map with K = 1, as for ACS |

Bands reduced: F606W, F814W, F438W (BELLS GALLERY) and F475X, F600LP (SL2S). F225W and F275W
were reduced, found to detect neither deflector nor arc, and deleted (2026-09-29).

---

## WFC3/IR (`drizzle_wfc3_ir.py`)
**Source:** WFC3 Data Handbook 2024 (`wfc3dhb2024_final.pdf`)

| Parameter | Value | Notes |
|---|---|---|
| File type | `_flt.fits` | IR detector; no CTE; no FLC |
| MAST project filter | `CALWF3` | `HAP-MVM`, `HAP-SVM`, `HLA` also appear |
| Pixel scale | 0.1283 → 0.06″/px | Handbook says 0.13 (WFC3 DHB p16); 0.1283 is the precise value. 0.06″ / pixfrac 1.0 settled by a scale/pixfrac scan (`AGENTS.md` *Output pixel scales*) |
| `driz_sep_bits` / `final_bits` | `'512'` | Handbook default is `'64,512'` (WFC3 DHB p210); bit 64 is unused on IR, so dropping it changes nothing. Written as `'512'`, never `''` (see below) |
| ERR units | electrons/s | So `1/sqrt(WHT)` must be rescaled by EXPTIME (`make_cutouts.weight_to_sigma_scale`), or by nothing for unequal-length stacks drizzled with an explicit IVM |
| CR rejection | none | Up-the-ramp fitting in calwf3 already rejects CRs; a CR pass eats the ~1 px PSF of field stars |
| TweakReg `conv_width` / `ylimit` | 2.5 / 0.2 | Only under `--align tweakreg` |

**Handbook quote (WFC3 DHB p210):** "The current default for IR channel data is for
AstroDrizzle to ignore (treat as bad) all pixels with any flag **except 512 and 64**."

**On IR CR rejection (DrizzlePac DHB p98):** WFC3/IR cosmic rays are already flagged during
up-the-ramp fitting. AstroDrizzle's CR steps can be turned off for IR; the handbook notes they
catch additional artefacts, but on these undersampled fields they also flag real star cores,
which is why the pipeline runs none and relies on correct DQ masking instead.

**`bits=''` trap:** `interpret_bit_flags('')` returns `None`, which disables DQ masking
entirely. `0` means "reject everything flagged"; `''`/`None` mean "reject nothing". The driver
asserts against it at import.

**DQ flag table (Table 3.3, WFC3 DHB PDF p89), IR channel:**

| Flag | Meaning | Kept? |
|---|---|---|
| 16 | Stable hot pixel | — (cleared per frame by `--dq-refine` where not locally deviant) |
| 32 | Unstable pixel | — (same) |
| 64 | WARMPIX, "unused"; §7.9.1: "currently not attached to any type of bad pixel" | — (harmless either way) |
| 256 | Full-well saturation | — |
| 512 | Bad/uncertain flat value, incl. IR blobs | ✓ (flat-corrected) |
| 1024 | SPIKE: CR spike detected during ramp fitting | — |
| 2048 | Signal in zero read | — |
| 4096 | AstroDrizzle CR | — (stale flags cleared by `--dq-refine`) |
| 8192 | calwf3 ramp-fit CR (ima only; not in flt) | — |

---

## WFPC2/WF3 (`drizzle_wfpc2_wf3.py`)
**Sources:** WFPC2 Data Handbook (`wfpc2_dhb.pdf`); WFPC2 Instrument Handbook (`wfpc2_ihb.pdf`)

| Parameter | Value | Notes |
|---|---|---|
| File type | `_flt.fits`, or `C0M`+`C1M` converted with `wfpc2_to_flt` | The WF3 chip (ext 3) is extracted and renumbered to `EXTVER=1` before drizzling |
| MAST project filter | `CALWFPC2` | `HAP-SVM`, `HLA` also appear |
| Which chip | **WF3**, not the PC | `DETECTOR=PC` in the primary header names the aperture; the full-field aperture (and SL2S's `WF3-FIX`) put the target on WF3 at ~(435, 424) |
| Pixel scale | 0.0996 → 0.05″/px | WF chips 0″.0996/px, PC 0″.0455/px (Table 1.1, WFPC2 DHB p17). The half-pixel `WFPC2-BOX` dither supports 0.05″ at pixfrac 1.0; the driver measures the dither phases and refuses to write otherwise |
| `driz_sep_bits` / `final_bits` | `'8,1024'` | bit 8 = A/D saturation; bit 1024 = repaired warm pixel (Table 3.3, WFPC2 DHB p54) |
| Weights | per-frame IVM via an `@`-association | DrizzlePac ignores ERR for WFPC2, and the ERR array is `sqrt(SCI)` with no gain or read noise. The IVM is `var_DN = SCI/gain + floor²`, gain from `ATODGAIN` (7.0), floor measured on each frame's blank sky |
| Gain / read noise | WF3 gain 6.90 ± 0.32 e⁻/DN at the "7" setting (header `ATODGAIN=7.0` is used); read noise 5.22 ± 0.28 e⁻ | WFPC2 IHB Table 4.2 (PDF p91) |
| A/D full well | 4095 DN | Used by the PSF star saturation gate |
| TweakReg `conv_width` / `ylimit` | 3.0 / 1 | Only under `--align tweakreg` |

`updatewcs(use_db=True)` runs even under `--align mast`, because the chip extraction drops
the NPOL/D2IM distortion arrays AstroDrizzle needs. The delivered WCS is GSC 2.4.0 only
(~0.3–1″ off absolute), fixed after drizzling by `align_wfpc2_to_acs.py`.

**On bit 8:** WFPC2 bit 8 is A/D converter saturation, not the aperture mask it means on
ACS. A/D saturation is rare on these lensing targets.

**DQ flag table (Table 3.3, WFPC2 DHB p54):**

| Flag | Meaning | Kept? |
|---|---|---|
| 2 | Calibration file defect (charge transfer traps) | — (filled by `--bcfill` on the `bcfill` branch) |
| 8 | A/D converter saturation (unrecoverable) | ✓ |
| 32 | Bad pixel | — |
| 256 | Questionable pixel (above charge trap) | — (filled by `--bcfill`) |
| 512 | Unrepaired warm pixel | — |
| 1024 | Repaired warm pixel | ✓ |
| 2048 | Uncorrected bias level < 100 DN | — |

The PC-chip extraction (`stale_scripts/drizzle_wfpc2_pc.py`) is superseded: it drizzled the
wrong chip, ~79″ from the lens.

---

## AstroDrizzle key parameters (DrizzlePac Handbook)
**Source:** The DrizzlePac Handbook Version 3 (`The_DrizzlePac_Handbook_Version3.pdf`)

### pixfrac and final_scale
- `pixfrac=1.0` is the pipeline default. MAST drizzled products use the native plate scale
  with pixfrac 1.0 (DrizzlePac DHB p135).
- Smaller pixfrac sharpens the PSF but lowers low-surface-brightness sensitivity (p131);
  recommended 0.7–1.0 for a few images (p148); keep ≥ 0.7 with 2–3 dithers (p130).
- The output PSF is the optics PSF convolved with the input pixel footprint scaled by pixfrac
  and with the drizzle kernel (p36, p131).
- Suggested fine output scales (p131): 0.0333″/px for ACS/WFC and WFC3/UVIS, 0.0667″/px for
  WFC3/IR.
- Weight-map check: rms/median of the weight map should stay below ~0.2–0.3 (p131, p148).
  `make_cutouts.py` records it per stamp as `WHTUNIF`.

**What this repo does with that.** Every band runs at pixfrac 1.0: ACS and UVIS at native
scale, WFC3/IR at 0.06″ and WFPC2 at 0.05″. Scans on real stamps showed a smaller pixfrac
moves noise between a pixel and its neighbours without changing what a fit sees, while
worsening weight uniformity (`AGENTS.md` *Output pixel scales*, *Drizzle correlated noise*).
Casertano et al. 2000's R = 1/r is stamped on every noise map as `CASR`; the separate
`NOISECOR` factor corrects the ERR arrays' understated per-pixel noise, which is not a
drizzle effect.

### Cosmic-ray rejection
- AstroDrizzle's test (p102): `|data − blotted| > scale × deriv + SNR × noise`, with
  `driz_cr_snr` and `driz_cr_scale` each a pair (primary detection, adjacent pixels).
- `driz_sep_pixfrac=1.0` for the separate drizzles the median is built from (p129).
- This repo instead masks CRs per frame with LACosmic (astroscrappy) into DQ bit 4096 and
  drizzles a plain weighted mean with `median=False, blot=False, driz_cr=False`. The reason is
  that the blotted-median reference reads low on a steep undersampled core when the dither is
  poor, flagging real core pixels. `resetbits=0` on that pass is mandatory, or AstroDrizzle
  clears bit 4096 and the mask is lost; the opt-in no-CR pass pins `resetbits=4096`.
- Per-lens LACosmic settings live in `lens_cr_params.json`.

### Sky subtraction
- Default `skymethod='localmin'` suits most fields (p128). For fields with large extended
  sources, `skystat='mode'` or `skysub=False` (p128). Incorrect sky subtraction biases CR
  identification and leaves tile artefacts (p128).

### DQ bit philosophy (p132, p148)
- `final_bits=0` rejects every flagged pixel and is too aggressive; which bits to keep
  depends on frame count and dither pattern.
- Holes in the weight map mean more bits need keeping or a larger pixfrac.
- `resetbits=4096` clears old AstroDrizzle CR flags on reprocessing.

### Drizzle kernels (p130, p138)
- `square` (default) is general purpose and is what every driver uses.
- `lanczos3` preserves the PSF on undithered data but rings around bad pixels and CRs.
- The point kernel (`pixfrac=0`) suits only very large, well-dithered sets.

---

## PSF notes

### Sampling by instrument (DrizzlePac DHB p19)

| Instrument | Sampling | Notes |
|---|---|---|
| ACS/WFC | ~Nyquist at λ > 6000 Å | Pixel ≈ PSF FWHM at optical wavelengths |
| WFC3/UVIS | ~Nyquist | Similar to ACS |
| WFC3/IR | Under-sampled | 0.13″/px against a ~0.15″ FWHM; why no CR pass is run |
| WFPC2/WF3 | Severely under-sampled | 0.0996″/px; the PC at 0.0455″ is also under-sampled |

### Sub-pixel dithering (DrizzlePac DHB p23–27)
- A 2-point half-pixel diagonal dither samples like a 45°-rotated array at √2 finer scale.
- A 4-point half-pixel dither on both axes recovers nearly all sub-pixel information; this is
  the `WFPC2-BOX` pattern that lets WF3 reach 0.05″.
- An 8-point dither is the choice for very accurate PSFs across orbits.

### Effective PSF (ePSF)
- Standard centroiding on under-sampled detectors carries pixel-phase bias up to ±0.1 px
  (DrizzlePac DHB p126); the ePSF method (Anderson & King 2000, 2006; Bellini et al. 2018)
  removes it.
- `photutils.psf.EPSFBuilder` implements it; `make_psf.py` builds an empirical ePSF from
  field stars in each mosaic and falls back to a model tier (ACS focus-diverse ePSF, STDPSF,
  or the MAST PSF database for WFPC2 F606W; Dauphin et al. 2021, WFC3 ISR 2021-12).
- A model ePSF is detector-frame and omits drizzle broadening (Anderson 2016, WFC3 ISR
  2016-12). `make_psf_inject.py` injects it into the exposures and re-drizzles, so the
  canonical model-tier kernel carries the broadening.

### Drizzle and the PSF
- `pixfrac > 0` broadens the output PSF by the input pixel footprint (p36); `pixfrac=1` is
  shift-and-add, `pixfrac=0` pure interlacing.
- Reducing `final_scale` below the native pixel creates no information but reduces aliasing.
