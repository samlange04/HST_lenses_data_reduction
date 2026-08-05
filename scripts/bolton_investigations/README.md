# Bolton-style interpolation investigations

Standalone, exploratory scripts (**not** part of the automated pipeline) that probe how
the *original* SLACS imaging pipeline (Bolton et al. 2008, SLACS V, arXiv:0805.1931)
handled — or rather, never had to handle — the **dead-column diagonal stripes** that this
repo's AstroDrizzle reduction correctly shows in its ACS noise maps, and what it would take
to reproduce that clean look while keeping drizzle's advantages.

All scripts here read the repo's real products under `data/` and write their outputs to
**`output/bolton_investigations/`** (git-ignored, regenerable). They are self-contained and
hard-wired to the demonstrator lens **J1023+4230 F814W**, whose bad detector column runs
straight through the deflector core — the worst case, and the clearest for a side-by-side.

Run any of them with `uv run python scripts/bolton_investigations/<script>.py`.

---

## Background: why our noise maps stripe and the legacy ones don't

The stripe is real and **correctly represented** — it is not a bug and `final_bits` is
already right (`'256,64,16'`; see CLAUDE.md → *AstroDrizzle key parameters* /
*Drizzle correlated noise*, and the `legacy-slacs-bolton-bilinear-no-stripes` memory).

- ACS dead columns are DQ **bit 128** (bad column, ~19 near-full columns on chip 1) plus
  **bit 4** (bad detector pixel, 2 columns). Neither is in this repo's ACS "good" set
  `{16,64,256}`, so AstroDrizzle drops them.
- Because the frames are **dithered**, a given striped output pixel is contributed to by
  fewer of the exposures (typically 3 of 4). Fewer frames → lower weight → higher noise in
  `1/sqrt(WHT)` → a diagonal stripe (diagonal because `final_rot=0` rotates the
  detector-vertical columns by the exposure roll).
- The stripe lives almost **entirely in the noise/weight map**. The *science* image is
  already stripe-free, because inverse-variance weighting fills each affected pixel from the
  good, dithered frames.

The "clean" legacy comparison maps (Nightingale/Etherington SLACS datasets) trace to
Bolton 2008's pipeline, which **masked only CRs + cold pixels** (no bad-column DQ handling),
rectified frames onto the 0.05″ grid by **bilinear interpolation** (not drizzle), and
**never built a per-pixel weight/noise map** — so a masked column has nothing to imprint on.
Confirmed from the papers themselves (Etherington 2022 arXiv:2202.09201 and Nightingale 2022
arXiv:2209.10566 both cite Bolton 2008 for their SLACS reduction). `data/pre_drizzled/` is a
*different* thing — a genuine AstroDrizzle/HLA re-drizzle — and it **does** stripe, like ours.

---

## The four options we built

| # | Script | What it does | Stripe? | Cost |
|---|---|---|---|---|
| 1 | `bolton_reduce.py` | Full Bolton bilinear reduction from the FLCs | none (by method) | seconds |
| 2 | `stripe_heal.py` | Keep drizzle image, interpolate the noise map across the stripe | removed post-hoc | instant |
| 3 | `redrizzle_bcfill.py` | Fill bad columns in the FLCs, un-flag, re-drizzle | none (by construction) | one re-drizzle (~40 s) |
| 4 | `compare_bolton_vs_drizzle.py` | Quantify Bolton-bilinear vs drizzle **science** | — (analysis) | seconds |

### 1. Full Bolton bilinear reduction — `bolton_reduce.py`
Reimplements the SLACS V recipe: mask **only** CRs (LACosmic/astroscrappy) + cold pixels
(no DQ bad-column masking); bilinear-rectify each FLC onto a North-up 0.05″/px grid
(`stwcs` full-distortion WCS → `scipy.ndimage.map_coordinates` order 1); combine with a
nan-aware sigma-clipped mean; build the noise map **from the combined counts** + a measured
background RMS (there is no weight map). Result: clean, no stripe — but a fundamentally
*different kind* of product (count-derived noise that traces flux; ~6 % softer PSF). This is
what the legacy datasets actually are.
Outputs: `bolton_J1023_{sci,noise}.fits`, `bolton_J1023_compare.png`.

### 2. Post-hoc stripe heal on the cutout — `stripe_heal.py`
Keeps the drizzle science untouched and interpolates the **noise map** across the stripe
where it crosses the lens. Detection uses the *geometric* detector-column angle (from the
FLC WCS, ~33.5° here — fixed, not searched, so a brighter nearly-parallel field feature such
as the satellite/CR trail can't hijack it) plus a Radon-style perpendicular-offset binning
that isolates the thin coherent line from compact sources. Two modes:
- default (**heal-down**): replace the stripe noise with the local across-stripe level →
  clean, legacy-style look.
- `--mask-up`: **inflate** the stripe noise (Etherington scalable-noise style) → the
  conservative choice; those pixels contribute ≈ nothing to a fit.

Outputs: `hybrid_J1023_noise.fits`, `hybrid_J1023_stripe.fits` (the mask),
`hybrid_J1023_compare.png`.

### 3. Input-level bad-column fill + re-drizzle — `redrizzle_bcfill.py`
The principled version — removes the *cause*, no detection at all. Per frame/chip: linearly
interpolate SCI **and** ERR *across* each DQ 4|128 column (from the good pixels either side),
clear those DQ bits, then re-drizzle with the pipeline's **exact** no-CR AstroDrizzle call
(`final_bits='256,64,16'`, ERR weighting, 0.05″/px, North-up at the lens, `num_cores=1`, mmap
write workaround). A baseline drizzle of the **unmodified** copies runs identically, so the
only difference between the two products is the fill → the weight map comes out uniform
*by construction*. The noise-difference panel shows the removed effect is *exactly* the two
stripes and nothing else; the science is ≈ unchanged (median |Δ| ~1.8e-6 e/s).
Each drizzle runs in its own subprocess (the repo's pattern, to keep a clean memory slate on
macOS). Outputs: `redrizzle_{baseline,filled}_{sci,noise}.fits`, `redrizzle_bcfill_compare.png`.

### 4. Bolton-bilinear vs drizzle **science** — `compare_bolton_vs_drizzle.py`
Registers (shift), background-matches and flux-matches the two images, then reports PSF
sharpness, residual, photometry and radial profiles. Output: `bolton_vs_drizzle.png`.

---

## Answers to the follow-up questions

### Q: "Is there a way to apply the interpolation to columns that go through the lens but have the image primarily produced by drizzle?"

**Yes** — and because the stripe is essentially *only* in the noise map, the drizzle image
is kept 100 % either way. Two routes, both built:

- **Post-hoc** (`stripe_heal.py`, option 2): detect the weight-deficit lines and interpolate
  the noise across them on the finished cutout. Instant, but relies on stripe detection.
- **Input-level** (`redrizzle_bcfill.py`, option 3): fill the flagged columns in the FLCs and
  re-drizzle → uniform weight by construction, no detection. Costs a re-drizzle; more robust.

Both share the same **honesty caveat**: because each striped sky pixel already had ~3 of 4
real frames, and this adds back an *interpolated* 4th (or heals the noise down to match), the
resulting noise map is optimistic by ~√(3/4) ≈ 13 % on those columns — the interpolated pixel
carries no independent information. That's why the science barely moves. For lens modelling,
if you don't trust those pixels, **inflating/masking** them (`stripe_heal.py --mask-up`) is
the conservative alternative. So the two orthogonal choices are:
- **where** to intervene — post-hoc cutout (2) vs input re-drizzle (3) vs full alt reduction (1);
- **what** to do to the flagged pixels — interpolate the noise *down* (make it look clean) vs
  inflate it *up* (down-weight it in the fit).

### Q: "How different are the Bolton-style science images to our drizzled science images?"

Very close. On J1023+4230 F814W, after registering + flux-matching (`compare_bolton_vs_drizzle.py`):

| Metric | Result |
|---|---|
| Stellar FWHM | drizzle 3.39 px (0.169″) vs bilinear 3.61 px (0.180″) → **bilinear ~6 % broader** |
| Aperture photometry (r<3″) | ratio 0.989 → agree to **~1 %** |
| Background noise | 1.02× → essentially identical |
| Core residual after alignment | **0.09 % of peak** |
| Deflector radial profile | overlap almost perfectly |

So photometrically and structurally the two are almost indistinguishable; the one real,
expected difference is that **bilinear resampling gives a ~6 % softer PSF**. That is the
whole reason the legacy images look fine for modelling but you'd still prefer drizzle
(sharper PSF *and* a real weight/noise map).

**Two caveats on generality:**
- J1023 F814W is **native-scale** ACS (0.05″→0.05″), where drizzle resamples very little, so
  bilinear is close to best-case here. The **oversampled** bands would soften more under
  bilinear: F160W (0.128″→0.06″) and WFPC2 F606W (0.0996″→0.05″) — exactly the bands where
  drizzle also wins most on correlated noise.
- The ~1 % flux offset is at the level of this quick reduction's sky/combine choices, not
  necessarily a real bias.

---

## If we productionize one of these

The two natural candidates to wire into the pipeline (the point of the `bolton_interpolation`
branch):
- **Option 3** as an opt-in flag on `drizzle_acs_wfc.py` (fill bad columns pre-drizzle) — the
  more correct, detection-free route.
- **Option 2** as an option in `make_cutouts.py` (post-hoc heal or `--mask-up`) — cheaper, no
  re-drizzle, and the `--mask-up` mode is the conservative modelling choice.

Both should carry the ~13 % optimism caveat above in their docs, and any "heal-down" default
should be justified against just masking those pixels instead.
