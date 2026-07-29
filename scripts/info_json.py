#!/usr/bin/env python
"""
Sample-nested tracking-JSON helpers shared by the drizzle and PSF scripts.

Every tracking JSON in info/ (lens_products.json, lens_exptime.json,
lens_instrument.json, lens_psf.json, lens_psf_injected.json, wfpc2_alignment.json,
psf_stars.json) is nested {sample: {lens: {key: value}}}, mirroring both
info/lens_samples.json's own top-level-by-sample layout and the
data/<sample>/<lens>/<filt>/ directory convention every script already uses. Before
this module, the tracking JSONs were flat {lens: {...}} with no sample namespace --
harmless only because no lens name has ever collided across slacs_gold/slacs_other/
gallery; nesting by sample removes that latent trap and lets a JSON diff show which
sample changed.
"""
import json


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
