"""Stamp-size-keyed output paths, shared by make_cutouts.py and make_mosaics.py.

The science stamps are 20" square (make_cutouts.py --size default). A second set cut at
a different size must not overwrite them, so every path that depends on the stamp size
is derived here rather than in each script:

    20" (the default)   data/cutouts/<sample>/...        data/mosaics/<sample>/
                        info/lens_cutout_qc.json
    any other size S    data/cutouts_<S>arcsec/<sample>/ data/mosaics_<S>arcsec/<sample>/
                        info/lens_cutout_qc_<S>arcsec.json

Keying the tree on --size itself, rather than on an independent --output flag the caller
has to remember to set, is deliberate: `make_cutouts.py --size 12` on its own then cannot
silently clobber the 20" product set. That is exactly the class of quietly-wrong-product
failure CLAUDE.md warns about, and the cutout FITS names (cutout_[cr_]{sci,noise}.fits)
carry no size in them, so a clobbered stamp is indistinguishable from a correct one on
inspection. An explicit --output still wins, for one-off work.

The PSF products in data/cutouts/ (cutout_[cr_]psf*.fits) are NOT size-keyed and are not
duplicated into a size tree: the kernel is trimmed by amplitude (CLAUDE.md, *PSF
generation*), so it is a property of the band, not of the stamp it will be convolved
with. A size-variant stamp pairs with the same kernel from the default tree.
"""
import os

# The pipeline's standard stamp size, in arcsec. Products at this size keep the
# unsuffixed paths every other script and every downstream reader already expects.
DEFAULT_SIZE = 20.0


def size_tag(size):
    """'' for the default size, else '_<size>arcsec' (e.g. '_12arcsec')."""
    return '' if float(size) == DEFAULT_SIZE else f'_{float(size):g}arcsec'


def cutouts_root(ws_path, size=DEFAULT_SIZE):
    """data/cutouts[_<S>arcsec] -- the tree holding every sample's stamps at `size`."""
    return os.path.join(ws_path, 'data', f'cutouts{size_tag(size)}')


def mosaics_root(ws_path, size=DEFAULT_SIZE):
    """data/mosaics[_<S>arcsec] -- the QC mosaics tiling those stamps."""
    return os.path.join(ws_path, 'data', f'mosaics{size_tag(size)}')


def qc_json_path(ws_path, size=DEFAULT_SIZE):
    """info/lens_cutout_qc[_<S>arcsec].json.

    A separate file per size, not a fourth nesting level inside one: info_json.update is
    {sample: {lens: {key: value}}} throughout info/, and the per-cutout diagnostics it
    records (weight_uniformity in particular) are measured over the cutout region, so
    they genuinely differ between sizes and must not overwrite each other.
    """
    return os.path.join(ws_path, 'info', f'lens_cutout_qc{size_tag(size)}.json')
