"""Stamp-size-keyed output paths, shared by the drizzle/cutout/mosaic scripts.

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
failure AGENTS.md warns about, and the cutout FITS names (cutout_[cr_]{sci,noise}.fits)
carry no size in them, so a clobbered stamp is indistinguishable from a correct one on
inspection. An explicit --output still wins, for one-off work.

**One stamp per (sample, lens, filt) -- the bcfill/standard split is a drizzle-layer axis
only** (2026-09-14). `variant` ('bcfill') keys a whole alternate *reduction*: the ACS/WFPC2
bad-column-filled re-drizzle (drizzle_acs_wfc.py --bcfill; see AGENTS.md, *Bad-column
fill*). The mosaics it produces genuinely differ, so `drizzled_root` is variant-keyed:

    standard   data/drizzled/          work dir data/drizzle_files/
    bcfill     data/drizzled_bcfill/   work dir data/drizzle_files_bcfill/

But the *cutouts* are not. bcfill supersedes standard for the bands it covers -- its stamps
are what the masks were drawn on and what the fits read -- so cutting a band with --bcfill
replaces that band's stamp in the one tree rather than starting a parallel one. There is
therefore exactly one science stamp per band, no priority-tree search to resolve, and the
four products a fit consumes (sci, noise, psf, mask) always sit together by construction.
Which reduction a stamp came from is read from its **BCFILL header card** (and mirrored in
the `bcfill` key of info/lens_cutout_qc.json), not from its path.

The guard that replaces the old parallel tree lives in make_cutouts.py: a standard cut
refuses to overwrite a stamp whose header says BCFILL=True unless --force, so re-cutting a
band without remembering --bcfill cannot silently downgrade it.

The PSF products (cutout_[cr_]psf*.fits) are NOT size-keyed and are not duplicated into a
size tree: the kernel is trimmed by amplitude (AGENTS.md, *PSF generation*), so it is a
property of the band, not of the stamp it will be convolved with. A size-variant stamp
pairs with the same kernel from the default tree -- which is what `psf_cutout_dir` below
exists to express.
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

    Unlike cutouts/mosaics this has no size axis (a mosaic is size-independent), and unlike
    them it IS variant-keyed: the bad-column-filled re-drizzle is a different mosaic, so it
    lands in data/drizzled_bcfill/ alongside data/drizzled/ instead of overwriting it. Both
    feed the one cutout tree; see the module docstring.
    """
    return os.path.join(ws_path, 'data', f'drizzled{variant_tag(variant)}')


def cutouts_root(ws_path, size=DEFAULT_SIZE):
    """data/cutouts[_<S>arcsec] -- the single science-stamp tree at `size`."""
    return os.path.join(ws_path, 'data', f'cutouts{size_tag(size)}')


def mosaics_root(ws_path, size=DEFAULT_SIZE):
    """data/mosaics[_<S>arcsec] -- the QC mosaics tiling those stamps."""
    return os.path.join(ws_path, 'data', f'mosaics{size_tag(size)}')


def qc_json_path(ws_path, size=DEFAULT_SIZE):
    """info/lens_cutout_qc[_<S>arcsec].json.

    A separate file per size, not a deeper nesting inside one: info_json.update is
    {sample: {lens: {key: value}}} throughout info/, and the per-cutout diagnostics it
    records (weight_uniformity in particular) are measured over the cutout region, so they
    genuinely differ between sizes and must not overwrite each other. Each band record
    carries a `bcfill` flag naming the reduction its stamp came from.
    """
    return os.path.join(ws_path, 'info', f'lens_cutout_qc{size_tag(size)}.json')


def cutout_dir(ws_path, sample, lens, filt, size=DEFAULT_SIZE):
    """The one cutout dir for this band: data/cutouts[_<S>arcsec]/<sample>/<lens>/<filt>."""
    return os.path.join(cutouts_root(ws_path, size), sample, lens, filt)


def psf_cutout_dir(ws_path, sample, lens, filt):
    """The cutout dir holding this band's modelling PSF products (cutout_[cr_]psf*.fits).

    Deliberately takes no `size`: the kernel is not size-keyed (see the module docstring),
    so it is always resolved in the DEFAULT_SIZE tree and a stamp cut at another size pairs
    with that same kernel. Taking a size here would let a --size 20 run look for a kernel
    that by design does not exist there. That size-independence is the whole of what this
    wrapper now asserts -- before 2026-09-14 it also had to pick between the bcfill and
    standard trees, which no longer exist separately.

    Used by every writer and reader of those files: make_psf.py / make_psf_inject.py /
    make_psf_err_injected.py write here, make_cutouts.py --psf-err and make_psf_mosaics.py
    read here, and make_dataset_subplots.py resolves the psf component through it.
    """
    return cutout_dir(ws_path, sample, lens, filt)
