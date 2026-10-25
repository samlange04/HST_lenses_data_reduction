# info/

The repo's catalogue and provenance layer: which lenses are in which sample, where each lens
is on the sky, what every pipeline stage did to every band, and the per-lens overrides the
scripts consult. Everything here is tracked in git and is small; it is what lets
`scripts/regenerate_files.py` rebuild the gitignored FITS layers on a fresh clone.

**Layout rule.** Every tracking JSON is nested `{sample: {lens: {band: record}}}`, matching
`info/lens_samples.json` and the `data/<layer>/<sample>/<lens>/<band>/` tree. Scripts read and
write them only through `scripts/info_json.py` (`load` / `update`), which keeps all three
levels sorted. The one exception is `lens_einstein_radii.json`, a catalogue keyed by lens
alone. A band key is the product directory name, so a split-visit WFPC2 lens has `f606W` and
`f606W_v2`; `null` for a band means "tried, no data or gated out", and an absent key means
"never run".

## Where to look

| question | file |
|---|---|
| Which lenses are in a sample, and how is each found on MAST? | `lens_samples.json` |
| Where is a lens on the sky? | `lens_coords.py` (merges the four `*_coords.py`) |
| What is its Einstein radius and redshift? | `lens_einstein_radii.json` |
| Which exposures went into a band, from which instrument, for how long? | `lens_products.json`, `lens_instrument.json`, `lens_exptime.json` |
| How was the stamp cut, and how good is its noise map? | `lens_cutout_qc.json` |
| What noise inflation does a band carry? | `noise_corr_factors.json` |
| Who drew the mask, from what, and what did they note? | `lens_masks.json` |
| Which images were marked, and on which band? | `lens_positions.json` |
| What PSF does a band have and how uncertain is it? | `lens_psf.json`, `lens_psf_injected.json` |
| Were any bands re-registered after drizzling? | `lens_field_ties.json`, `wfpc2_alignment.json` |
| Does a lens have non-default CR or PSF-star settings? | `lens_cr_params.json`, `psf_stars.json` |
| Why these DQ bits, scales and kernels? | `instrument_drizzle_ref.md` |

## Sample definition and positions (hand-maintained)

- **`lens_samples.json`** is the single source of truth for sample membership. Per sample:
  `description`, `coords` (which positions table), and for the SL2S samples `uvis_filters`
  (which WFC3/UVIS bands the UVIS runner loops) and `center_bands` (the order of preference
  for the band every stamp is centred on). Per lens: `mast_target` where the MAST name is not
  `SDSS<J-name>` (the `GAL-*` names, full SDSS designations, `SL2SJ…`), `force_copy` (use the
  `-COPY` visit because the primary has dead frames), `proposals` (restrict a band to listed
  programs), `archive_proposals` (take a band from another program's data), `cutout_box`
  (narrower peak search where arcs outshine the deflector), and a `note` saying why the lens
  is in. Read it through `scripts/mast_target_names.py`, never directly.
- **`slacs_coords.py`, `gallery_coords.py`, `bells_coords.py`, `sl2s_coords.py`** hold each
  sample's catalogue positions, transcribed from the source paper named in the file header
  (Bolton+2008 Table 4; Shu+2016 Table 1; Brownstein+2012 Table 2; Gavazzi+2014 Tables 1 and
  3 for SL2S, since the SL2S designations are truncated). **`lens_coords.py`** merges them into
  one dict keyed on the short J-name and raises on a collision. Every drizzle pins its output
  tangent point to this position.
- **`lens_einstein_radii.json`** (flat, keyed by lens): `theta_e_arcsec`, `theta_e_err_arcsec`
  where the paper gives one, `RE_kpc`, `zlens`, `zsrc`, `sigma_kms`, and a `source` naming the
  paper and table. Four series with different provenance: Auger+2009 for SLACS (fetched from
  VizieR by `detect_arcs.py --fetch-theta-e`, which merges and never overwrites a non-null
  θ_E with a null), Shu+2016 Table 2 for GALLERY, Brownstein+2012 Table 5 for BELLS, and
  Sheu+2025 Table 1 or Sonnenfeld+2013 Table 6 for SL2S (hand-entered). Rows with no lens
  model carry `theta_e_arcsec: null` but real redshifts; every reader must allow that.

## Written by the drizzle drivers

- **`lens_products.json`**: the exposure rootnames that reached the drizzle for each band
  (not the whole download; sub-10 s frames and dropped visits are excluded).
- **`lens_instrument.json`**: `INSTRUME/DETECTOR` per band, e.g. `ACS/WFC`, `WFC3/UVIS`,
  `WFPC2/WF3` (the chip, not MAST's `WFPC2/PC` label). Also the guard against two instruments
  claiming one band key: a driver that finds another instrument recorded here exits with no
  data.
- **`lens_exptime.json`**: total exposure time per band from the drizzled header.

## Written by the cutout, mask, positions and QC stages

- **`lens_cutout_qc.json`**: how each 12″ stamp was cut. `stamp` (the reduction name, e.g.
  `cutout_cr`), `drizzle_pass` (`cr` or `nocrrej`), `center_source` / `center_from` (which
  band or stamp set the centre), `offset_arcsec` (how far the stamp centre moved from the
  catalogue position), `noise_k` (the EXPTIME factor for WFC3/IR), `corr_factor` and
  `corr_factor_source` (the applied noise inflation and where it came from), and two noise-map
  diagnostics: `weight_uniformity` against `weight_uniformity_limit` 0.2, and
  `casertano_r_analytic`. The `bcfill` / `crfill` / `dropped_frames` flags name variant
  reductions that live only on the `bcfill` branch. A 20″ re-cut writes its own
  `lens_cutout_qc_20arcsec.json`, which is gitignored.
- **`lens_masks.json`**: provenance of each band's contaminant mask. `source` is `drawn`,
  `accepted_from_<band>`, `edited_from_<band>` or `reprojected_from_<band>`; with
  `proposal_from`, `proposal_px`, `n_added_px`, `n_erased_px` and the final `n_excluded_px`.
  Hand-added keys carry the review: `note` (the per-lens verdict), `edited` (a scripted
  pull-back or addition), `field_tie_shift` and `recentred` (integer pixel rolls applied when
  the stamp under the mask moved). **A `--force` redraw replaces the record, so re-add the
  hand keys afterwards.**
- **`lens_arc_masks.json`**: the same for the optional arc masks (`n_arc_px`, the detector's
  proposal under `detector_*`, union bookkeeping under `proposal_px_by_source`). The arc-mask
  FITS files themselves are untracked, so a clone has these records but not the masks.
- **`lens_positions.json`**: one record per lens under the band it was marked on.
  `n_positions`, `positions_arcsec`, `broadcast_to` (the other bands that received the same
  sky positions), the snap statistics (`snap_moved_px_*`), the display levers in use when the
  marks were made, and `n_positions_on_contaminant`. A lens judged unmarkable has
  `n_positions: 0` and a `deliberately_empty` reason; downstream code must skip the positions
  likelihood for it rather than build one from zero points.
- **`lens_dataset_subplots.json`**: provenance of the `_dataset.png` QC panel per band
  (`frac_masked`, stretch, output path).

## Written by the PSF stages

- **`lens_psf.json`**: per band, `method` (`empirical`, `empirical_donor`, `empirical_median`,
  `inject_stdpsf`, `inject_acs_fdpsf`, `inject_wfpc2_psfdb`, or a `broadened_*` fallback),
  `n_stars`, `n_crowded_rejected`, `fwhm_pix`, kernel sizes, `pedestal_frac`, and the error
  map's `err_method` (`bootstrap`, `jackknife`, `ensemble_broadened`,
  `calibrated_vs_empirical`) with `psf_err_frac` (and its `_core` / `_wing` split where
  measured) and `err_lower_bound`. Donor-built kernels list their `donors`; those lists can
  name lenses since removed from the roster, whose PSF stars remain valid.
- **`lens_psf_injected.json`**: an audit trail of every injection run (the drizzle-broadened
  model-tier kernel), independent of whether it was promoted to canonical.

## Alignment records

- **`wfpc2_alignment.json`**: `{sample: {lens: "mast" | "tweakreg"}}`, the per-lens choice
  `run_wfpc2_wf3.sh` passes as `--align`. Only `slacs_gold` was audited lens by lens; every
  other WFPC2 lens falls back to `mast`, which is also every audited lens's value.
- **`lens_field_ties.json`**: the field-source ties `align_field_sources.py --apply` made
  between an SL2S lens's secondary bands and its F475X: `ref`, and a `history` of each pass
  with the measured offset (`dE_mas`, `dN_mas`, `err_mas`, `n_sources`) and the residual
  after.

## Hand-maintained overrides the scripts consult

- **`lens_cr_params.json`**: per-lens LACosmic settings (`lacosmic_sigclip`,
  `lacosmic_objlim`, and `lacosmic_relax_radius`, the radius about the lens inside which the
  relaxed values apply). Precedence: script default < this file < CLI flag.
- **`psf_stars.json`**: per-lens PSF-star selection overrides (`include` / `exclude`
  positions or boxes `[xmin, xmax, ymin, ymax]`, and any `make_psf` parameter such as
  `max_stars`, `min_snr`, `satlevel`, `psf_star_pass`, `psf_err`). Same precedence.
- **`noise_corr_factors.json`**: the noise-map inflation `make_cutouts.py` applies, keyed by
  instrument/detector then filter, plus `_lens_overrides` per sample and lens. The `_method`
  string records how the values were measured and when. These correct the ERR arrays'
  understated per-pixel noise, not drizzle correlation; see `AGENTS.md` *Drizzle correlated
  noise*.
- **`lens_crfill.json`**: the two cosmic-ray-filled stamps, with the exact commands, the
  tracks filled and the noise caveat. Both live only on the `bcfill` branch (`on_main:
  false`).

## Reference material

- **`instrument_drizzle_ref.md`**: the handbook facts behind each driver's DQ bits, pixel
  scales, weights and kernels, with page citations to the editions listed at its top and
  download links for the PDFs (which are not part of the repo).
