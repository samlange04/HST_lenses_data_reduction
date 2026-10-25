# HST lens image reduction — `bcfill` branch (alternative reduction)

> **If you use this repository, please credit the repository
> ([github.com/samlange04/HST_lenses_data_reduction](https://github.com/samlange04/HST_lenses_data_reduction))
> and Samuel C. Lange (The Chinese University of Hong Kong).**

> **This is not the main branch, and development on it stopped on 2026-09-22.** On that date
> the `main` branch went back to the standard reduction, and this branch was kept only as an
> alternative for comparison. Afterwards it received some of `main`'s changes, the last on
> 2026-10-21, but not the BELLS and SL2S samples, not the noise-map correction, and nothing
> since. **For science, use [`main`](https://github.com/samlange04/HST_lenses_data_reduction/tree/main).**

This repository turns Hubble Space Telescope (HST) archive images of **strong gravitational
lenses** into clean, ready-to-model image cutouts.

A strong lens is a massive foreground galaxy (the *deflector*) whose gravity bends the light
of a more distant galaxy (the *source*) behind it, so the source appears as arcs, a ring, or
two to four separate images around the deflector. Modelling those arcs measures how mass is
distributed in the deflector. A model needs more than a picture, though: it needs to know how
uncertain every pixel is, how blurred the telescope makes a point of light, which pixels to
ignore, and roughly where the multiple images are. This repository produces all of that.

**Recommended analysis tool: [PyAutoLens](https://github.com/PyAutoLabs/PyAutoLens)**, from
the [PyAutoLabs](https://github.com/PyAutoLabs) software family. The products are written in
the formats it reads directly (its masks, image positions and PSF kernels), and its own GUI
tools drew the masks and marked the positions.

## What this branch does differently

Two of HST's cameras (ACS and WFPC2) have **dead detector columns**. The standard reduction
on `main` simply leaves them out when it combines the exposures, so the pixels they would
have covered rest on fewer exposures. That shows up as faint diagonal stripes of slightly
higher noise in the noise map — which is the honest answer, since those pixels really are
measured less well.

This branch instead **fills the dead columns by interpolating across them** before combining
(*bcfill* = bad-column fill). The noise map comes out smooth, but the filled pixels are about
13% more confident than the data justify, so a lens model over-weights them. That is why
`main` uses the standard reduction and this branch is kept only for comparison.

Also unlike `main`, **the noise maps here carry no measured correction**: HST's own error
estimates, which come out a few to ~25% too small on these images, are used as they are
(`NOISECOR` = 1.0 in every header).

## What's in the box

| Sample folder | Lenses | Survey | HST filters (number of lenses with each) |
|---|---|---|---|
| `slacs_gold` | 38 | SLACS (Bolton et al. 2008; Auger et al. 2009), the Etherington et al. (2022) "Gold" lenses — see below | F814W (38), F606W (22), F555W (16), F160W (13) |
| `slacs_other` | 24 | SLACS, the remaining lenses with a clean early-type deflector | F606W (24), F160W (6), F814W (4) |
| `gallery` | 15 | BELLS GALLERY (Shu et al. 2016) | F606W (15), F814W (6), F438W (6) |

That is 151 lens-and-filter images. The column fill applies to the 105 ACS and WFPC2 images
(F814W, F555W, F606W in the two SLACS samples); the WFC3 images — F160W and all of
`gallery` — have no dead columns and are identical to the standard reduction.

**Where the name `slacs_gold` comes from.** Etherington et al. (2022, *Automated
galaxy-galaxy strong lens modelling: no lens left behind*, MNRAS 517, 3275) modelled 43 SLACS
lenses automatically with PyAutoLens and graded each fit Gold, Silver or Bronze. `slacs_gold`
is exactly the 38 SLACS lenses graded **Gold** there. Four of the paper's other five SLACS
lenses (J0959+4416, J1016+3859, J1153+4612, J1416+5136) are in `slacs_other` instead.

A filter name gives its approximate central wavelength in nanometres: F438W is blue, F555W
and F606W are visible light, F814W is near-infrared, and F160W is infrared. Every cutout
covers the same 12″ × 12″ patch of sky, centred on the deflector, at 0.05″ per pixel (ACS and
WFPC2), 0.0396″ (WFC3/UVIS, `gallery`) or 0.06″ (WFC3/IR, F160W).

## Using the products

Each lens and filter has one folder, `data/cutouts/<sample>/<lens>/<filter>/`. The image and
noise map are named after the reduction they come from:

| File | What it is |
|---|---|
| `cutout_cr_bcfill_sci.fits` | **The image**, column-filled (ACS and WFPC2 bands). Brightness per pixel in electrons per second (counts per second for WFPC2 F606W — check `BUNIT`). |
| `cutout_cr_bcfill_noise.fits` | **The noise map**: the 1σ uncertainty of every pixel, same units as the image. |
| `cutout_cr_sci.fits`, `cutout_cr_noise.fits` | The same, standard reduction, in `gallery` bands (nothing to fill). F160W bands use `cutout_sci.fits` / `cutout_noise.fits`. |
| `cutout_cr_bcfill_crfill_sci.fits` | Only J1213+6708 and J1250+0523 in F814W: a cosmic-ray streak crossing the lensed arc was also filled in (details in `info/lens_crfill.json`). |
| `cutout_cr_psf.fits` | **The PSF** (point-spread function): how the telescope blurs a single point of light, on the image's pixel grid, normalised to sum to 1. |
| `cutout_cr_psf_err.fits` | The per-pixel uncertainty of that PSF. |
| `cutout_cr_mask.fits` | **The mask**: 1 marks a pixel to *exclude* from a fit (nearby galaxies, stars, defects), 0 a pixel to use. Hand-drawn; identical to `main`'s. |
| `cutout_cr_positions.json` | **Image positions**: where the multiple images of the lensed source are, hand-marked. |
| `cutout_cr_dataset.png` | A quick-look figure of the image, noise, PSF and mask together. |

The `cr` in the names means cosmic-ray hits were removed. Masks, positions and PSFs describe
the pixel grid, which is identical in both reductions, so they are shared with `main`.

**Orientation.** North is towards increasing row number and West towards increasing column
number, so with Python's usual display (row 0 at the top) **North appears at the bottom and
East on the left** — a vertical flip of the usual astronomical view, not a rotation. The FITS
headers carry full sky coordinates (WCS).

**Image positions** are in PyAutoLens's format: `(y, x)` offsets in arcseconds from the
cutout centre, with +y pointing South and +x pointing West. To get pixels:
`row = row_centre − y / pixel_size`, `col = col_centre + x / pixel_size` (the centre of a
240-pixel cutout is 119.5). **An empty positions file is deliberate**: on 16 lenses the
images were too uncertain to mark.

```python
from astropy.io import fits

band = "data/cutouts/slacs_gold/J0008-0004/f814W/"
image = fits.getdata(band + "cutout_cr_bcfill_sci.fits")
noise = fits.getdata(band + "cutout_cr_bcfill_noise.fits")
psf   = fits.getdata(band + "cutout_cr_psf.fits")
mask  = fits.getdata(band + "cutout_cr_mask.fits").astype(bool)   # True = exclude
```

Overview figures tiling every lens of a sample are in `data/mosaics/<sample>/`. Per-lens notes
on the masks (which nearby objects are contaminants, which masks are provisional) are in
[`AGENTS.md`](AGENTS.md) and in the `note` fields of `info/lens_masks.json`.

## How the reduction works

HST takes several short exposures of each target, shifted slightly between exposures, each
spoiled by cosmic-ray hits and detector defects. For every lens and filter the pipeline:

1. **Downloads** the calibrated exposures from the HST archive (MAST) and the matching
   calibration files.
2. **Fills dead detector columns** by interpolation (this branch only; ACS and WFPC2).
3. **Removes cosmic rays** exposure by exposure (F160W's standard calibration already does
   this).
4. **Combines the exposures** with STScI's *drizzle* algorithm onto a common North-up grid
   centred on the lens, producing the image and a weight map, from which the noise map is made.
5. **Aligns the filters** with each other, so the same pixel is the same point on the sky in
   every colour.
6. **Cuts out** the 12″ stamp around the deflector and measures the **PSF** from stars in the
   image, or from STScI's PSF models where there are too few.
7. **Masks and marks** by hand in small interactive tools.

Every choice behind these steps is documented in [`AGENTS.md`](AGENTS.md), a technical
reference that is also the guide for AI coding assistants working in this repository. On this
branch it describes the reduction as of the stop date.

## Re-running the reduction

Install [uv](https://docs.astral.sh/uv/), then run `uv sync` in the repository folder; it
installs the exact Python and package versions used. The column fill is the `--bcfill` option
of the ACS and WFPC2 scripts, e.g. for lens J0008-0004 in F814W:

```bash
uv run python scripts/drizzle_acs_wfc.py --sample slacs_gold --lens J0008-0004 --filt f814W --bcfill
uv run python scripts/make_cutouts.py    --sample slacs_gold --lens J0008-0004 --filt f814W --bcfill
```

The batch runners take it as a second argument (`bash scripts/run_acs_all.sh slacs_gold
--bcfill`, likewise `run_wfpc2_wf3.sh`; ACS first, because WFPC2 is aligned to it). A full
re-reduction downloads tens of gigabytes. The run order is in `AGENTS.md`.

## Credit and acknowledgements

If you use this repository, please credit the repository
([github.com/samlange04/HST_lenses_data_reduction](https://github.com/samlange04/HST_lenses_data_reduction))
and **Samuel C. Lange (The Chinese University of Hong Kong)**.

Please also cite the survey papers for the lenses you use (SLACS: Bolton et al. 2008, Auger
et al. 2009, and Etherington et al. 2022 for the `slacs_gold` selection; BELLS GALLERY: Shu et
al. 2016), and acknowledge the HST data, e.g. *"Based on observations made with the NASA/ESA
Hubble Space Telescope, obtained from the Mikulski Archive for Space Telescopes (MAST) at the
Space Telescope Science Institute, which is operated by the Association of Universities for
Research in Astronomy, Inc., under NASA contract NAS 5-26555."*

The reduction is built on STScI's [DrizzlePac](https://drizzlepac.readthedocs.io),
[Astropy](https://www.astropy.org), [photutils](https://photutils.readthedocs.io),
[astroscrappy](https://github.com/astropy/astroscrappy) (cosmic-ray removal) and
[PyAutoLens](https://github.com/PyAutoLabs/PyAutoLens) (mask and position tools).
