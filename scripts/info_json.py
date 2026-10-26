#!/usr/bin/env python
"""
Sample-nested tracking-JSON helpers shared by the drizzle and PSF scripts.

Every tracking JSON in info/ (lens_products.json, lens_exptime.json,
lens_instrument.json, lens_psf.json, lens_psf_injected.json, wfpc2_alignment.json,
psf_stars.json) is nested {sample: {lens: {key: value}}}, mirroring both
info/samples/lens_samples.json's own top-level-by-sample layout and the
data/<sample>/<lens>/<filt>/ directory convention every script already uses. Before
this module, the tracking JSONs were flat {lens: {...}} with no sample namespace --
harmless only because no lens name has ever collided across slacs_gold/slacs_other/
gallery; nesting by sample removes that latent trap and lets a JSON diff show which
sample changed.

**Where each file lives (2026-10-26).** info/ is split by who writes a file, and `path()`
below is the one place that knows which subfolder:

    samples/      lens_samples.json                         the sample definition
    coords/       lens_coords.py + one *_coords.py per sample
    catalogues/   lens_einstein_radii.json                  literature values per lens
    overrides/    hand-edited per-lens settings the scripts consult
    provenance/   records the pipeline writes as it runs
"""
import json
import os

INFO_LAYOUT = {
    'samples': ('lens_samples.json',),
    'catalogues': ('lens_einstein_radii.json',),
    'overrides': ('lens_cr_params.json', 'psf_stars.json', 'noise_corr_factors.json',
                  'wfpc2_alignment.json', 'lens_crfill.json'),
    'provenance': ('lens_products.json', 'lens_instrument.json', 'lens_exptime.json',
                   'lens_stamp_qc.json', 'lens_masks.json', 'lens_arc_masks.json',
                   'lens_positions.json', 'lens_psf.json', 'lens_psf_injected.json',
                   'lens_field_ties.json', 'lens_dataset_subplots.json'),
}
_SUBFOLDER = {name: sub for sub, names in INFO_LAYOUT.items() for name in names}


def path(ws_path, name):
    """Absolute path of info file `name` (e.g. 'lens_psf.json') in its subfolder.
    Size variants of the stamp-QC record ('lens_stamp_qc_20arcsec.json') go with it."""
    sub = _SUBFOLDER.get(name) or ('provenance' if name.startswith('lens_stamp_qc') else None)
    if sub is None:
        raise KeyError(f'{name} is not a known info/ file (see info_json.INFO_LAYOUT)')
    return os.path.join(ws_path, 'info', sub, name)


def load(path):
    """Whole tracking JSON, or {} if it doesn't exist yet."""
    try:
        with open(path) as f:
            return json.load(f)
    except FileNotFoundError:
        return {}


def update(path, sample, lens, key, value):
    """Merge {sample: {lens: {key: value}}} into `path`, sorted at every level."""
    data = load(path)
    sample_data = data.setdefault(sample, {})
    lens_entry = sample_data.setdefault(lens, {})
    lens_entry[key] = value
    sample_data[lens] = dict(sorted(lens_entry.items()))
    data[sample] = dict(sorted(sample_data.items()))
    with open(path, 'w') as f:
        json.dump(dict(sorted(data.items())), f, indent=4)


def band_owner_conflict(path, sample, lens, key, instrument_prefix):
    """The instrument recorded for this lens's band in `path` (lens_instrument.json) if it is
    a DIFFERENT instrument from `instrument_prefix`, else None.

    The pipeline keys products by filter alone (data/<tree>/<sample>/<lens>/<filt>/, and
    every tracking JSON), so one filter can hold one instrument's product per lens. On
    2026-10-04 the UVIS runner, looping F814W over SL2S, found J0959+0206's ACS F814W FLCs in
    that band's calibrated dir, took them as its own download and replaced the ACS product,
    and its F606W "no data" exit nulled the WFPC2 F606W records. The drizzle scripts call
    this first and leave a band owned by another instrument untouched."""
    owner = entry(path, sample, lens, key)
    if owner and not str(owner).startswith(instrument_prefix):
        return owner
    return None


def entry(path, sample, lens, key, default=None):
    """Read {sample: {lens: {key: ...}}} from `path`, or `default` if any level
    is absent. The read-only counterpart to update(); used for per-lens override
    JSONs (e.g. info/overrides/lens_cr_params.json) that a script consults but never writes."""
    return load(path).get(sample, {}).get(lens, {}).get(key, default)
