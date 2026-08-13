"""Stamp-size- and variant-keyed output paths, shared by the drizzle/cutout/mosaic scripts.

The science stamps are 12" square (make_cutouts.py --size default) -- switched from 20"
2026-08-11, since the 12" tree is the one carrying scripts/make_masks.py's hand-drawn,
non-regenerable masks and is what downstream modelling reads. A second set cut at a
different size must not overwrite it, so every path that depends on the stamp size is
derived here rather than in each script:

    12" (the default)   data/cutouts/<sample>/...        data/mosaics/<sample>/
                        info/lens_cutout_qc.json
    any other size S    data/cutouts_<S>arcsec/<sample>/ data/mosaics_<S>arcsec/<sample>/
                        info/lens_cutout_qc_<S>arcsec.json

Keying the tree on --size itself, rather than on an independent --output flag the caller
has to remember to set, is deliberate: `make_cutouts.py --size 20` on its own then cannot
silently clobber the 12" product set. That is exactly the class of quietly-wrong-product
failure CLAUDE.md warns about, and the cutout FITS names (cutout_[cr_]{sci,noise}.fits)
carry no size in them, so a clobbered stamp is indistinguishable from a correct one on
inspection. An explicit --output still wins, for one-off work.

`variant` is a second, orthogonal keying axis for a whole alternate *reduction* (not just a
different crop of the same mosaic). The one variant so far is 'bcfill' -- the ACS
bad-column-filled re-drizzle (drizzle_acs_wfc.py --bcfill; see
scripts/bolton_investigations/redrizzle_bcfill.py and CLAUDE.md). It suffixes the drizzled
tree too, since the mosaics themselves differ, and composes before the size tag:

    variant 'bcfill', 12"   data/drizzled_bcfill/  data/cutouts_bcfill/  data/mosaics_bcfill/
                            info/lens_cutout_qc_bcfill.json
    variant 'bcfill', 20"   data/cutouts_bcfill_20arcsec/  (etc.)

The default-size bcfill QC JSON (lens_cutout_qc_bcfill.json) does NOT match the gitignore's
`*arcsec.json`, so it is tracked -- deliberately, since data/cutouts_bcfill/ is a tracked
tree (a small alternate science product worth versioning, like data/cutouts/). A bcfill
tree at a non-default size stays untracked, same as any other size variant.

The PSF products in data/cutouts/ (cutout_[cr_]psf*.fits) are NOT size-keyed and are not
duplicated into a size tree: the kernel is trimmed by amplitude (CLAUDE.md, *PSF
generation*), so it is a property of the band, not of the stamp it will be convolved
with. A size-variant stamp pairs with the same kernel from the default tree. (make_psf.py
and friends hardcode data/cutouts/data/mosaics rather than going through this module, so
they always target whichever tree is currently unsuffixed.)
"""
import os

# The pipeline's standard stamp size, in arcsec. Products at this size keep the
# unsuffixed paths every other script and every downstream reader already expects.
DEFAULT_SIZE = 12.0


def size_tag(size):
    """'' for the default size, else '_<size>arcsec' (e.g. '_12arcsec')."""
    return '' if float(size) == DEFAULT_SIZE else f'_{float(size):g}arcsec'


def variant_tag(variant):
    """'' for the standard reduction, else '_<variant>' (e.g. '_bcfill')."""
    return f'_{variant}' if variant else ''


def drizzled_root(ws_path, variant=''):
    """data/drizzled[_<variant>] -- the tree holding the drizzled mosaics.

    Unlike cutouts/mosaics this has no size axis (a mosaic is size-independent); it is
    variant-keyed only, so the bad-column-filled re-drizzle lands in data/drizzled_bcfill/
    alongside the standard data/drizzled/ instead of overwriting it.
    """
    return os.path.join(ws_path, 'data', f'drizzled{variant_tag(variant)}')


def cutouts_root(ws_path, size=DEFAULT_SIZE, variant=''):
    """data/cutouts[_<variant>][_<S>arcsec] -- stamps for `variant` at `size`."""
    return os.path.join(ws_path, 'data', f'cutouts{variant_tag(variant)}{size_tag(size)}')


def mosaics_root(ws_path, size=DEFAULT_SIZE, variant=''):
    """data/mosaics[_<variant>][_<S>arcsec] -- the QC mosaics tiling those stamps."""
    return os.path.join(ws_path, 'data', f'mosaics{variant_tag(variant)}{size_tag(size)}')


def qc_json_path(ws_path, size=DEFAULT_SIZE, variant=''):
    """info/lens_cutout_qc[_<variant>][_<S>arcsec].json.

    A separate file per (variant, size), not a deeper nesting inside one: info_json.update
    is {sample: {lens: {key: value}}} throughout info/, and the per-cutout diagnostics it
    records (weight_uniformity in particular) are measured over the cutout region, so
    they genuinely differ between reductions/sizes and must not overwrite each other.
    """
    return os.path.join(ws_path, 'info',
                        f'lens_cutout_qc{variant_tag(variant)}{size_tag(size)}.json')
