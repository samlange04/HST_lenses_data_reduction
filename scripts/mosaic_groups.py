"""Shared per-sample (group_name -> filter-precedence-list) definitions for
make_mosaics.py and make_psf_mosaics.py, so a lens contributes exactly one panel per
group. SLACS' WFPC2 F606W and ACS F555W never coexist on one lens, so they merge into
one group with the filter recorded per-panel (and the split-visit f606W_v1/v2 keys
preferred behind the combined f606W); gallery's WFC3/UVIS filters never merge.
"""
import glob
import os

_SLACS_GROUPS = {
    'f814W':       ['f814W'],
    'f606W_f555W': ['f606W', 'f606W_v2', 'f606W_v1', 'f555W'],
    'f160W':       ['f160W'],
}
_GALLERY_FILTERS = ['f225W', 'f275W', 'f438W', 'f606W', 'f814W']

SAMPLE_GROUPS = {
    'slacs_gold':  _SLACS_GROUPS,
    'slacs_other': _SLACS_GROUPS,
    'gallery':     {f: [f] for f in _GALLERY_FILTERS},
}


def groups_for_sample(sample, cutouts_dir):
    """(group_name -> precedence list) for `sample`. Falls back to one singleton group
    per filter subdirectory actually found on disk, so a sample not yet in
    SAMPLE_GROUPS -- or a new filter added later -- still produces mosaics instead of
    silently being skipped."""
    if sample in SAMPLE_GROUPS:
        return SAMPLE_GROUPS[sample]
    filts = set()
    for lens_dir in glob.glob(os.path.join(cutouts_dir, '*')):
        for filt_dir in glob.glob(os.path.join(lens_dir, '*')):
            if os.path.isdir(filt_dir):
                filts.add(os.path.basename(filt_dir))
    return {f: [f] for f in sorted(filts)}
