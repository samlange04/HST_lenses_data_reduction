# HST lens image reduction

> **If you use this repository, please credit the repository
> ([github.com/samlange04/HST_lenses_data_reduction](https://github.com/samlange04/HST_lenses_data_reduction))
> and Samuel C. Lange (The Chinese University of Hong Kong).**

This repository turns Hubble Space Telescope (HST) archive images of **strong gravitational
lenses** into clean, ready-to-model image cutouts.

A strong lens is a massive foreground galaxy (the *deflector*) whose gravity bends the light
of a more distant galaxy (the *source*) behind it, so the source appears as arcs, a ring, or
two to four separate images around the deflector. Modelling those arcs measures how mass is
distributed in the deflector. A model needs more than a picture, though: it needs to know how
uncertain every pixel is, how blurred the telescope makes a point of light, which pixels to
ignore, and roughly where the multiple images are. This repository produces all of that, for
133 lenses from five published lens surveys.

**Most people only need the finished products in `data/cutouts/`** — they are in git, so a
plain clone has them. Re-running the reduction is only needed to change how it was done.

**Recommended analysis tool: [PyAutoLens](https://github.com/PyAutoLabs/PyAutoLens)**, from
the [PyAutoLabs](https://github.com/PyAutoLabs) software family. The products are written in
the formats it reads directly (its masks, image positions and PSF kernels), and its own GUI
tools drew the masks and marked the positions.

## What's in the box

| Sample folder | Lenses | Survey | HST filters (number of lenses with each) |
|---|---|---|---|
| `slacs_gold` | 38 | SLACS (Bolton et al. 2008; Auger et al. 2009), the Etherington et al. (2022) "Gold" lenses | F814W (38), F606W (22), F555W (16), F160W (13) |
| `slacs_other` | 24 | SLACS, the remaining lenses with a clean early-type deflector | F606W (24), F160W (6), F814W (4) |
| `gallery` | 15 | BELLS GALLERY (Shu et al. 2016) | F606W (15), F814W (6), F438W (6) |
| `bells` | 21 | BELLS (Brownstein et al. 2012) | F814W (21) |
| `sl2s_dinos_a` | 9 | SL2S (Sonnenfeld et al. 2013) + Project Dinos (Sheu et al. 2025) | F475X (9), F606W (8), F814W (3), F600LP (1), F160W (1) |
| `sl2s_dinos_x` | 11 | as above, F475X-only lenses | F475X (11), F600LP (3) |
| `sl2s_wf` | 10 | SL2S with older HST WFPC2 imaging | F606W (10), F475X (2), F600LP (2) |
| `sl2s_snap` | 5 | SL2S with short "snapshot" imaging | F475X (5), F600LP (5) |

That is 232 lens-and-filter images in total. A filter name gives its approximate central
wavelength in nanometres: F438W is blue, F606W and F475X/F600LP are visible light, F814W is
near-infrared, and F160W is infrared. Seeing the same lens in several colours helps tell the
(usually blue) lensed source apart from the (red) deflector.

The data come from four HST cameras. They differ in pixel size, so the cutouts do too:

| Camera | Filters here | Pixel size | Cutout size in pixels |
|---|---|---|---|
| ACS/WFC | F814W, F555W, F606W (one lens) | 0.05″ | 240 × 240 |
| WFPC2 (its WF3 chip) | F606W | 0.05″ (resampled from 0.1″) | 240 × 240 |
| WFC3/UVIS | F438W, F475X, F600LP, F606W, F814W | 0.0396″ | 303 × 303 |
| WFC3/IR | F160W | 0.06″ (resampled from 0.13″) | 200 × 200 |

Every cutout covers the same 12″ × 12″ patch of sky, centred on the deflector.

## Using the finished products

Each lens and filter has one folder, `data/cutouts/<sample>/<lens>/<filter>/`, for example
`data/cutouts/slacs_gold/J0008-0004/f814W/`. The files in it are:

| File | What it is |
|---|---|
| `cutout_cr_sci.fits` | **The image.** Brightness per pixel, in electrons per second (counts per second for WFPC2 F606W — check the `BUNIT` header keyword before comparing filters). |
| `cutout_cr_noise.fits` | **The noise map**: the 1σ uncertainty of every pixel, in the same units as the image. |
| `cutout_cr_psf.fits` | **The PSF** (point-spread function): how the telescope blurs a single point of light, on the same pixel grid and orientation as the image, normalised to sum to 1. A lens model convolves with this. |
| `cutout_cr_psf_err.fits` | The per-pixel uncertainty of that PSF. |
| `cutout_cr_mask.fits` | **The mask**: 1 marks a pixel to *exclude* from a fit (nearby galaxies, stars, detector defects), 0 a pixel to use. The lens itself is never masked. Hand-drawn. |
| `cutout_cr_positions.json` | **Image positions**: where the multiple images of the lensed source are, hand-marked (see below). |
| `cutout_cr_dataset.png` | A quick-look figure of the image, noise, PSF and mask together — the fastest way to check a lens. |

The `cr` in the names means cosmic-ray hits were removed (below). F160W images have no
cosmic-ray step, so their files are named `cutout_sci.fits`, `cutout_noise.fits` and so on.
A few lenses carry extra files beside these (for example `cutout_psf_injected.fits`, an
alternative PSF kept for comparison); the files in the table are the ones to use.

**Orientation.** North is towards increasing row number and West towards increasing column
number, so with Python's usual display (row 0 at the top) **North appears at the bottom and
East on the left**. That is a vertical flip of the way astronomical images are usually shown,
not a rotation. The FITS headers carry full sky coordinates (WCS) if you need them.

**Image positions** are stored in the format of the [PyAutoLens](https://github.com/PyAutoLabs/PyAutoLens)
lens-modelling software, which reads them directly. Each entry is a `(y, x)` offset in
arcseconds from the cutout centre. Following PyAutoLens's convention, +y points South and +x
points West. To get pixels: `row = row_centre − y / pixel_size` and
`col = col_centre + x / pixel_size`, where the centre of a 240-pixel cutout is 119.5.
**An empty positions file is deliberate**: on 24 lenses the images were too uncertain to
mark, or the counter-images are too faint to see. It does not mean the job was skipped.

A minimal example in Python:

```python
from astropy.io import fits

band = "data/cutouts/slacs_gold/J0008-0004/f814W/"
image = fits.getdata(band + "cutout_cr_sci.fits")
noise = fits.getdata(band + "cutout_cr_noise.fits")
psf   = fits.getdata(band + "cutout_cr_psf.fits")
mask  = fits.getdata(band + "cutout_cr_mask.fits").astype(bool)   # True = exclude

signal_to_noise = image / noise
```

**Overview figures** tiling every lens of a sample, one sheet per filter, are in
`data/mosaics/<sample>/` (`*_signal.png`, `*_snr.png`, `*_noise.png`, `*_psf.png`).

### Things worth knowing before you model

- **The noise maps already include a measured correction.** HST's own error estimates come
  out too small on most of these images: by 4–26% depending on camera and filter, and up to
  ~40% on a few low-background F475X fields (the infrared F160W needs none). Each noise map
  has been scaled up by the measured factor, recorded in the `NOISECOR` header keyword.
  Neighbouring pixels also share some noise (a side effect of combining several exposures).
  The noise maps are correct for anything summed over a few pixels or more, which is what
  matters for a lensed arc.
- **Some masks are judgement calls.** A short verdict for every lens where it mattered —
  which nearby objects are contaminants and which are lensed images, which masks are
  provisional, and which lenses need something extra in the model (a second source, a
  nearby perturbing galaxy) — is in [`AGENTS.md`](AGENTS.md) under *Masks → Per-lens
  verdicts*, and in the `note` fields of `info/lens_masks.json`.
- **Einstein radii, redshifts and velocity dispersions** from the survey papers are collected
  in `info/lens_einstein_radii.json`, each row naming its source.
- **Arc masks** (masks isolating just the lensed light) are optional and not included.
- **Some lenses were removed** from the published samples, mostly because a nearby galaxy is
  massive enough to need its own lens model, or because no second image of the source can be
  seen. [`AGENTS.md`](AGENTS.md) lists each one and why.

## How the reduction works

HST does not deliver one finished picture. It takes several short exposures of each target,
shifted slightly between exposures (*dithered*), each spoiled by cosmic-ray hits and detector
defects. For every lens and filter the pipeline:

1. **Downloads** the calibrated exposures from the HST archive (MAST) and the matching
   calibration files.
2. **Removes cosmic rays** — the sharp streaks and dots left when charged particles hit the
   detector — exposure by exposure (the infrared camera's standard calibration already
   removes them, so F160W skips this).
3. **Combines the exposures** into one sharper, deeper image with STScI's *drizzle*
   algorithm, onto a common North-up grid centred on the lens. The same step produces a
   weight map, from which the noise map is made.
4. **Aligns the filters** with each other, so the same pixel is the same point on the sky in
   every colour.
5. **Cuts out** the 12″ stamp around the deflector and writes the noise map.
6. **Measures the PSF** from stars in the same image where there are enough, or from STScI's
   PSF models otherwise.
7. **Masks and marks** — the masks and image positions are drawn and clicked by hand in
   small interactive tools, then checked against the published lens models and Einstein
   radii.

Every choice behind these steps — and the mistakes each one avoids — is documented in
[`AGENTS.md`](AGENTS.md). It is written as a technical reference (it is also the guide for
AI coding assistants working in this repository), so start here and go there for the detail.

## Re-running the reduction

You only need this to change how the data were processed, or to rebuild the large
intermediate files that are not stored in git.

**Setup.** Install [uv](https://docs.astral.sh/uv/), then from the repository folder:

```bash
uv sync
```

This installs the exact Python version and package versions the pipeline was run with, on
macOS or Linux, with nothing else to install. Run any script as `uv run python scripts/<name>.py`.

**Rebuild everything git does not carry.** The full reduction downloads tens of gigabytes
and takes hours per sample. One command rebuilds the missing layers in the right order, and
can be re-run to resume after an interruption:

```bash
uv run python scripts/regenerate_files.py --sample slacs_gold --dry-run   # show the plan, run nothing
uv run python scripts/regenerate_files.py --sample slacs_gold --light     # quick-look PNGs only: minutes, no download
uv run python scripts/regenerate_files.py --sample slacs_gold             # download + full re-reduction
```

**Reduce one lens in one filter**, e.g. lens J0008-0004 in F814W:

```bash
uv run python scripts/drizzle_acs_wfc.py --sample slacs_gold --lens J0008-0004 --filt f814W
uv run python scripts/make_cutouts.py    --sample slacs_gold --lens J0008-0004 --filt f814W
```

There is one combining script per camera (`drizzle_acs_wfc.py`, `drizzle_wfpc2_wf3.py`,
`drizzle_wfc3_uvis.py`, `drizzle_wfc3_ir.py`) and a batch runner per camera
(`scripts/run_*_all.sh <sample>`). The order they must run in, and why, is in `AGENTS.md`. A
lens with no data in the requested filter prints `NO DATA` and is skipped; that is normal.

## Repository layout

```
data/cutouts/<sample>/<lens>/<filter>/   the finished products (in git)
data/mosaics/<sample>/                   overview figures (in git)
data/...                                 downloads and intermediate files (not in git, rebuildable)
info/                                    sample lists, lens coordinates, catalogue values, and a
                                         record of how every product was made (JSON)
scripts/                                 the pipeline
AGENTS.md                                full technical reference
```

`info/lens_samples.json` is the master list of which lenses are in which sample.

**Branches.** `main` holds the standard reduction described here and is the one to use. The
`bcfill` branch holds an alternative in which dead detector columns are filled in by
interpolation before combining. It looks cleaner but makes the noise map slightly
optimistic, so it is kept only for comparison.

## Credit and acknowledgements

If you use this repository, please credit the repository
([github.com/samlange04/HST_lenses_data_reduction](https://github.com/samlange04/HST_lenses_data_reduction))
and **Samuel C. Lange (The Chinese University of Hong Kong)**. GitHub's **"Cite this
repository"** button (right-hand sidebar, from `CITATION.cff`) gives a ready-made citation.

**Licences.** The code is under the [MIT License](LICENSE). The data products — everything
in `data/`, and the JSON/CSV/Markdown records in `info/` — are under
[CC BY 4.0](LICENSE-DATA): you may share and adapt them for any purpose, as long as you give
credit as above.

Please also cite the survey papers for the lenses you use (SLACS: Bolton et al. 2008, Auger
et al. 2009, and Etherington et al. 2022 for the `slacs_gold` selection; BELLS: Brownstein et
al. 2012; BELLS GALLERY: Shu et al. 2016; SL2S:
Sonnenfeld et al. 2013; Project Dinos: Sheu et al. 2025), and acknowledge the HST data, e.g.
*"Based on observations made with the NASA/ESA Hubble Space Telescope, obtained from the
Mikulski Archive for Space Telescopes (MAST) at the Space Telescope Science Institute, which
is operated by the Association of Universities for Research in Astronomy, Inc., under NASA
contract NAS 5-26555."*

The reduction is built on STScI's [DrizzlePac](https://drizzlepac.readthedocs.io),
[Astropy](https://www.astropy.org), [photutils](https://photutils.readthedocs.io),
[astroscrappy](https://github.com/astropy/astroscrappy) (cosmic-ray removal) and
[PyAutoLens](https://github.com/PyAutoLabs/PyAutoLens) (mask and position tools).
