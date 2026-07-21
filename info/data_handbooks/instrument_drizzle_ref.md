# HST Drizzle & PSF Reference
Condensed from STScI data handbooks in this directory. Sources cited as (DHB p###).

---

## Quick-reference table

| Instrument | File | Scale (″/px) | `driz_sep_bits` / `final_bits` | Output suffix |
|------------|------|-------------|--------------------------------|---------------|
| ACS/WFC | FLC | 0.05 | `256,64,16` (=336) | `_drc_` |
| WFC3/IR | FLT | 0.1283 | `64,512` | `_drz_` |
| NICMOS/NIC2 | CAL | 0.0756 | `2,4,8` | `_drz_` |
| WFPC2/PC | FLT | 0.0455 | `8,1024` | `_drw_` |

---

## ACS/WFC (`drizzle_acs_wfc.py`)
**Source:** ACS Data Handbook (`acs_dhb.pdf`)

| Parameter | Value | Notes |
|-----------|-------|-------|
| File type | `_flc.fits` | CTE-corrected; FLT also produced but use FLC |
| Pipeline | CALACS | `jref` env var → `references/hst/acs` |
| MAST project filter | `CALACS` | Also `HAP-MVM`, `HAP-SVM` in products — filter correctly |
| Output suffix | `_drc_` | FLC → drc (distortion + CTE corrected) |
| Pixel scale | 0.05 ″/px | Native ~0.049 (ACS DHB p12); pipeline/IDCTAB uses 0.05 (confirmed from drizzle headers p40–42) |
| `driz_sep_bits` / `final_bits` | `'256,64,16'` = 336 | **Official STScI MDRIZTAB default since June 2017** (DrizzlePac DHB p132, ACS ISR 2017-05) |
| `driz_cr_snr` | `'3.5 3.0'` | |
| `driz_cr_scale` | `'1.2 0.7'` | |
| TweakReg `conv_width` | 3.5 | |
| TweakReg `ylimit` | 0.2 | |

**Why bit 256 (saturation) is included:** The 2017 MDRIZTAB update (ACS ISR 2017-05) explicitly sets `driz_sep_bits=336=16+64+256`. Saturated pixels are included so AstroDrizzle's median can span all input pixels (including bright saturated galaxy centres), enabling proper comparison and CR rejection. Unstable hot pixels (bit 32) remain excluded.

**DQ Flag Table (Table 3.4, ACS DHB p87):**
| Flag | Meaning | In bits? |
|------|---------|----------|
| 4 | Bad/vignetted pixel | — |
| 8 | Masked by aperture feature | — |
| 16 | Hot pixel (dark > 0.14 e⁻/s), **stable** | ✓ accepted |
| 32 | Unstable dark / RTS noise | — excluded |
| 64 | Warm pixel (dark 0.06–0.14 e⁻/s) | ✓ accepted |
| 256 | Full-well saturation | ✓ accepted (per MDRIZTAB 2017) |
| 512 | Bad pixel in reference file | — |
| 1024 | Sink pixel / charge trap | — |
| 2048 | A-to-D saturation | — |
| 4096 | AstroDrizzle CR | — |
| 8192 | acsrej CR | — |

---

## WFC3/IR (`drizzle_wfc3_ir.py`)
**Source:** WFC3 Data Handbook 2024 (`wfc3dhb2024_final.pdf`)

| Parameter | Value | Notes |
|-----------|-------|-------|
| File type | `_flt.fits` | IR detector; no CTE; no FLC |
| Pipeline | CALWF3 | `iref` env var → `references/hst/wfc3` |
| MAST project filter | `CALWF3` | Also `HAP-MVM`, `HAP-SVM`, `HLA` in products |
| Output suffix | `_drz_` | FLT → drz |
| Pixel scale | 0.1283 ″/px | Handbook says 0.13 (WFC3 DHB p16); 0.1283 is more precise |
| `driz_sep_bits` / `final_bits` | `'64,512'` | **Handbook-explicit default** (WFC3 DHB p210) |
| `driz_cr_snr` | `'3.5 3.0'` | |
| `driz_cr_scale` | `'1.2 0.7'` | |
| TweakReg `conv_width` | 2.5 | |
| TweakReg `ylimit` | 0.2 | |

**Handbook quote (WFC3 DHB p210):** "The current default for IR channel data is for AstroDrizzle to ignore (treat as bad) all pixels with any flag **except 512 and 64**."

**Note on IR CR rejection (DrizzlePac DHB p98):** WFC3/IR cosmic rays are already flagged during up-the-ramp fitting (calwf3). AstroDrizzle's CR rejection steps can be turned off for IR, but running them catches additional detector artifacts not present in the DQ arrays.

**DQ Flag Table (Table 3.3, WFC3 DHB p39) — IR channel:**
| Flag | Meaning | In bits? |
|------|---------|----------|
| 16 | Stable hot pixel | — |
| 32 | Unstable pixel | — |
| 64 | (Obsolete: Warm pixel — not used for IR) | ✓ accepted (safe; unused) |
| 256 | Full-well saturation | — |
| 512 | Bad/uncertain flat value, incl. IR blobs | ✓ accepted (flat-corrected) |
| 1024 | (Unused for IR) | — |
| 2048 | Signal in zero read | — |
| 4096 | AstroDrizzle CR | — |
| 8192 | calwf3 ramp-fit CR (ima only; not in flt) | — |

---

## NICMOS/NIC2 (`drizzle_nic2.py`)
**Source:** NICMOS Data Handbook (`nic_dhb.pdf`)

| Parameter | Value | Notes |
|-----------|-------|-------|
| File type | `_cal.fits` | IR detector; pipeline produces CAL not FLT |
| Pipeline | CALNICA | `nref` env var → `references/hst/nicmos` |
| MAST project filter | `CALNIC` | **Not** `CALNICA` — confirmed from MAST product tables |
| Output suffix | `_drz_` | CAL → drz |
| Pixel scale | 0.0756 ″/px | Geometric mean of X=0.075948, Y=0.075355 (NIC DHB Table 5.2) |
| `driz_sep_bits` / `final_bits` | `'2,4,8'` | Uncertain calibration imperfections — acceptable |
| `driz_cr_snr` | `'4.0 3.5'` | Higher threshold due to lower S/N vs. optical |
| `driz_cr_scale` | `'1.2 0.7'` | |
| TweakReg `conv_width` | 2.5 | |
| TweakReg `ylimit` | 0.2 | |
| TweakReg `minobj` | 5 | Fewer stars in IR fields |
| `CAMERA` header keyword | integer | Read as `str(header.get('CAMERA', 2))` — not a string |

**Note:** NICMOS images often contain pedestal effects (differing bias levels between quadrants) requiring extra processing. Per DrizzlePac DHB p122: "NICMOS data may require special attention: images often contain additional signal in the sky, persistence or pedestal effects."

**DQ Flag Table (NIC DHB Table 2.3):**
| Flag | Meaning | In bits? |
|------|---------|----------|
| 1 | Telemetry error | — |
| 2 | Uncertain linearity correction | ✓ accepted |
| 4 | Uncertain dark correction | ✓ accepted |
| 8 | Uncertain flat correction | ✓ accepted |
| 16 | Grot (large dust particle) | — |
| 32 | Defective pixel | — |
| 64 | Saturated pixel | — |
| 128 | Missing data | — |
| 256 | Bad calibration pixel | — |
| 512 | Cosmic ray | — |
| 1024 | Source | — |
| 2048 | 0th-read signal | — |
| 4096 | CR from MultiDrizzle | — |

---

## WFPC2/PC (`drizzle_wfpc2_pc.py`)
**Source:** WFPC2 Data Handbook (`wfpc2_dhb.pdf`)

| Parameter | Value | Notes |
|-----------|-------|-------|
| File type | `_flt.fits` | Extract PC chip (group 1) before drizzling |
| Pipeline | CALWFPC2 | `uref` env var → `references/hst/wfpc2` |
| MAST project filter | `CALWFPC2` | Also `HAP-SVM`, `HLA` in products |
| Output suffix | `_drw_` | WFPC2 FLT → drw convention |
| Pixel scale (PC) | 0.0455 ″/px | Confirmed Table 1.1 (WFPC2 DHB p17): "PC 36″×36″, 0″.0455/pixel" |
| Pixel scale (WF) | 0.0996 ″/px | For reference; script uses PC chip only |
| `driz_sep_bits` / `final_bits` | `'8,1024'` | bit 8 = A/D saturation; bit 1024 = repaired warm pixel |
| `driz_cr_snr` | `'15.0 10.0'` | High threshold; WFPC2 noise characteristics |
| `driz_cr_scale` | `'1.5 1.0'` | |
| TweakReg `conv_width` | 3.0 | |
| TweakReg `ylimit` | 1 | |

**Note on bit 8 (WFPC2):** WFPC2 DQ bit 8 = A/D converter saturation (Table 3.3, WFPC2 DHB p54) — different from ACS where bit 8 = aperture mask. For SLACS/BELLS lensing targets, A/D saturation is rare so this is unlikely to affect results in practice.

**DQ Flag Table (Table 3.3, WFPC2 DHB p54):**
| Flag | Meaning | In bits? |
|------|---------|----------|
| 2 | Calibration file defect (charge transfer traps) | — |
| 8 | A/D converter saturation (unrecoverable) | ✓ accepted |
| 32 | Bad pixel | — |
| 256 | Questionable pixel (above charge trap) | — |
| 512 | Unrepaired warm pixel | — |
| 1024 | Repaired warm pixel | ✓ accepted |
| 2048 | Uncorrected bias level < 100 DN | — |

---

## AstroDrizzle Key Parameters (DrizzlePac Handbook)
**Source:** The DrizzlePac Handbook Version 3 (`The_DrizzlePac_Handbook_Version3.pdf`)

### pixfrac and final_scale
- **`pixfrac=1.0`** is the pipeline default (conservative). Drizzled products from MAST use native plate scale with pixfrac=1.0 (DrizzlePac DHB p135).
- For dithered data, **smaller pixfrac improves resolution but reduces low-surface-brightness sensitivity** (p131).
- Recommended range: **0.7–1.0** for routine observations with a few images (p148).
- PSF is convolved **in quadrature** with both `final_scale` and `pixfrac` (p131).
- With few images (2–3 dithers): keep `pixfrac >= 0.7` to ensure good output coverage (p130).
- With many dithers (>6): can reduce pixfrac toward 0.5 and `final_scale` toward half-native for better resolution.
- Suggested output scales (p131): 0.03333″/px for ACS/WFC3-UVIS; 0.06666″/px for WFC3/IR.
- **Check weight map:** std/mean of weight map should stay **above ~0.2–0.3** to avoid S/N loss (p131, p148).

### CR rejection
- CR formula (p102): `|data - blotted| > scale × deriv + SNR × noise`
- `driz_cr_snr`: two values — first for primary CR detection, second (lower) for adjacent pixels.
- `driz_sep_pixfrac=1.0` recommended for the separate drizzle step used in CR detection (p129).
- For IR: up-the-ramp fitting (calwf3) already flags CRs, but running AstroDrizzle CR steps still useful for catching additional artifacts (p98).

### Sky subtraction
- Default `skymethod='localmin'` works for most fields (p128).
- For fields with large extended sources (like lens galaxies): consider `skystat='mode'` instead of default `median`, or turn off sky subtraction (`skysub=False`) (p128).
- Incorrect sky subtraction can bias CR identification and create tile artefacts (p128).

### DQ bit mask philosophy (p132, p148)
- Default `driz_sep_bits=0` / `final_bits=0` in AstroDrizzle is over-aggressive — rejects all flagged pixels.
- Choice of bits depends on number of input frames and dither pattern.
- "Holes" in weight map → may need to include more DQ bits or increase pixfrac.
- During pipeline reprocessing, use `resetbits=4096` to clear old AstroDrizzle CR flags before re-running.

### Drizzle kernels (p130, p138)
- **`square`** (default): general purpose.
- **`lanczos3`**: best for single or non-dithered images; preserves PSF well but creates ringing artefacts around bad pixels/CRs.
- Point kernel (`pixfrac=0`): only for very large, well-dithered datasets (HUDF-type).

---

## PSF Notes (DrizzlePac Handbook)

### HST PSF sampling by instrument
| Instrument | Sampling | Notes |
|------------|----------|-------|
| ACS/WFC | ~Nyquist at λ > 6000 Å | Pixel ≈ PSF FWHM at optical |
| WFC3/UVIS | ~Nyquist | Similar to ACS |
| WFC3/IR | Under-sampled | 0.13″/px vs PSF FWHM ~0.15″ |
| WFPC2/PC | Under-sampled | 0.046″/px, WF chip severely under-sampled |
| NICMOS/NIC2 | Under-sampled | 0.0756″/px |

- WFC3/UVIS, WFC3/IR, ACS/WFC: pixel widths comparable to PSF FWHM — not fully Nyquist (p19).
- WFPC2 WF chips severely under-sampled; WFPC2 PC and NIC2 also under-sampled (p19).

### Sub-pixel dithering and PSF recovery
- **2-point dither** (1/2-pixel offset on diagonal): provides sampling equivalent to array rotated 45° at √2 finer scale (p23).
- **4-point dither** (1/2-pixel offsets on both axes): recovers nearly all sub-pixel information (p25). Recommended default.
- **8-point dither** (4-point × 2-point secondary): best for very accurate PSFs across orbits (p27).
- Smaller `pixfrac` in drizzle gives sharper effective PSF in the output image (p36, p131).

### ePSF (effective PSF) method
- For under-sampled detectors (WFPC2, NIC3, WFC3/IR), standard centroiding in TweakReg produces systematic residuals up to ±0.1 pixel due to pixel-phase bias (p126).
- The **ePSF method** (Anderson & King 2006; Bellini et al. 2018) accounts for under-sampling and gives much better alignment (p126).
- WFC3 PSF reference: WFC3 IHB Sections 6.6.1 and 7.6 (DrizzlePac DHB p137).
- ACS PSF reference: ACS IHB Section 5.6 (DrizzlePac DHB p137).
- Tools: `photutils.psf.EPSFBuilder` (Python) implements the Anderson & King ePSF method.

### Drizzle effect on PSF
- Setting `pixfrac > 0` broadens the output PSF by convolving with input pixel footprint (p36).
- `pixfrac=0` = pure interlacing (no broadening), requires very uniform dither coverage.
- `pixfrac=1` = shift-and-add (maximum broadening).
- The PSF in the drizzled image is the **convolution** of: optics PSF ⊗ pixel response(pixfrac) ⊗ drizzle kernel.
- Reducing `final_scale` below the native pixel does not create information but does reduce aliasing.
