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
failure AGENTS.md warns about, and the cutout FITS names (cutout_[cr_][<variant>_]{sci,noise}.fits)
carry no size in them, so a clobbered stamp is indistinguishable from a correct one on
inspection. An explicit --output still wins, for one-off work.

**One stamp per (sample, lens, filt), and the stamp's filename says which reduction it is**
(2026-09-22). `variant` ('bcfill', 'crfill', 'bcfill_crfill', 'drop', ...) keys a whole
alternate *reduction* at the drizzle layer (drizzle_acs_wfc.py --bcfill / --crfill; see
AGENTS.md, *Bad-column fill*), so `drizzled_root` is variant-keyed:

    standard   data/drizzled/          work dir data/drizzle_files/
    bcfill     data/drizzled_bcfill/   work dir data/drizzle_files_bcfill/

The cutouts stay in ONE tree, but a stamp cut from a variant drizzle carries that variant
in its NAME, not just in a header card:

    standard      cutout_cr_sci.fits           cutout_cr_noise.fits           cutout_cr.png
    bcfill        cutout_cr_bcfill_sci.fits    cutout_cr_bcfill_noise.fits    cutout_cr_bcfill.png
    bcfill+crfill cutout_cr_bcfill_crfill_sci.fits  ...

Only the sci/noise stamps and their 3-panel PNG are tagged. The mask, arc-mask, positions
and PSF files keep the bare `cutout_[cr_]` prefix: they describe the band's *grid*, which
is identical across reductions (same drizzle call, same output WCS), so they are shared by
every reduction of that band and can be moved between branches on their own.

Why the tag is in the name and not only in the header (it used to be header-only,
2026-09-14..22): the repo's `main` branch carries the standard stamps and the `bcfill`
branch carries the bcfill ones. With one shared filename, merging or cherry-picking
between those branches silently replaced one reduction's bytes with the other's -- no
conflict, no error, and nothing downstream could tell. With distinct names a merge can at
worst put BOTH stamps in one dir, and `find_stamp` below then refuses to pick, which is
the loud failure we want.

The invariant is still ONE stamp per band: make_cutouts.py refuses to write a second
variant beside an existing one unless --force, and with --force it removes the other.
Every reader resolves the stamp through `find_stamp(cutout_dir, prefix, kind)` rather than
spelling the filename, so a reader never has to know which reduction a band carries.

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


# The reduction axes a drizzle/cut can be run with, in the order their tags compose:
# --bcfill --crfill gives 'bcfill_crfill' (and data/drizzled_bcfill_crfill/), never
# 'crfill_bcfill'. make_cutouts.py builds its variant from this; find_stamp parses with it.
VARIANT_COMPONENTS = ('bcfill', 'crfill', 'drop')
STAMP_KINDS = ('sci', 'noise')
PREFIXES = ('cutout_cr', 'cutout')     # CR-rejected pass first: the preferred one


def variant_from_flags(**flags):
    """'bcfill_crfill' from bcfill=True, crfill=True, drop=False -- the canonical join."""
    unknown = set(flags) - set(VARIANT_COMPONENTS)
    if unknown:
        raise ValueError(f'unknown variant component(s) {sorted(unknown)}')
    return '_'.join(c for c in VARIANT_COMPONENTS if flags.get(c))


def _is_variant(tag):
    """True if `tag` is a canonical composition of VARIANT_COMPONENTS ('' included)."""
    if tag == '':
        return True
    parts = tag.split('_')
    order = [VARIANT_COMPONENTS.index(x) for x in parts if x in VARIANT_COMPONENTS]
    return len(order) == len(parts) and order == sorted(order) and len(set(parts)) == len(parts)


def stamp_name(prefix, kind, variant=''):
    """'cutout_cr_bcfill_sci.fits' -- the tagged stamp filename. `kind` is 'sci' or
    'noise'; the 3-panel QC PNG is stamp_png_name."""
    return f'{prefix}{variant_tag(variant)}_{kind}.fits'


def stamp_png_name(prefix, variant=''):
    return f'{prefix}{variant_tag(variant)}.png'


class AmbiguousStampError(RuntimeError):
    """More than one reduction's stamp sits in a cutout dir -- see the module docstring."""


def list_stamps(cutout_dir, prefix, kind='sci'):
    """Every `{prefix}[_<variant>]_<kind>.fits` in cutout_dir, as {variant: path}.

    Matches on the canonical variant grammar, so 'cutout' does not swallow 'cutout_cr_*'
    and a stray file with an unknown middle is not mistaken for a stamp.
    """
    import glob
    found = {}
    for path in glob.glob(os.path.join(cutout_dir, f'{prefix}*_{kind}.fits')):
        base = os.path.basename(path)
        middle = base[len(prefix):-len(f'_{kind}.fits')]
        if middle == '':
            found[''] = path
        elif middle.startswith('_') and _is_variant(middle[1:]):
            found[middle[1:]] = path
    return found


def find_stamp(cutout_dir, prefix, kind='sci', variant=None):
    """The one `{prefix}[_<variant>]_<kind>.fits` in cutout_dir, or None if there is none.

    With `variant=None` (the normal case) the band must carry exactly ONE reduction's
    stamp; two or more raise AmbiguousStampError naming them, because choosing silently is
    precisely the failure the tagged names exist to prevent. Pass `variant` to ask for a
    specific one (returns None if absent).
    """
    stamps = list_stamps(cutout_dir, prefix, kind)
    if variant is not None:
        return stamps.get(variant)
    if len(stamps) > 1:
        names = ', '.join(os.path.basename(p) for _, p in sorted(stamps.items()))
        raise AmbiguousStampError(
            f'{cutout_dir}: {len(stamps)} {prefix}_{kind} stamps from different reductions '
            f'({names}). One stamp per band: remove the one that does not belong '
            f'(or re-cut with --force), then re-run.')
    return next(iter(stamps.values()), None)


def stamp_variant(path):
    """'bcfill' from '.../cutout_cr_bcfill_sci.fits', '' for an untagged stamp."""
    base = os.path.basename(path)
    for prefix in PREFIXES:
        for kind in STAMP_KINDS:
            if base.startswith(prefix) and base.endswith(f'_{kind}.fits'):
                middle = base[len(prefix):-len(f'_{kind}.fits')]
                if middle == '':
                    return ''
                if middle.startswith('_') and _is_variant(middle[1:]):
                    return middle[1:]
    raise ValueError(f'{base} is not a stamp filename')


def find_prefix(cutout_dir, drizzle_pass='auto'):
    """Pick 'cutout_cr' or 'cutout' for this band, mirroring make_cutouts.py --pass auto:
    prefer the CR-rejected pass, fall back to no-CR (F160W has no CR pass). Returns None
    if the requested pass has no sci stamp here, whichever reduction it is from. Raises
    AmbiguousStampError if a pass has stamps from two reductions (find_stamp).
    """
    has_cr = find_stamp(cutout_dir, 'cutout_cr', 'sci') is not None
    has_nocr = find_stamp(cutout_dir, 'cutout', 'sci') is not None
    if drizzle_pass == 'cr':
        return 'cutout_cr' if has_cr else None
    if drizzle_pass == 'nocrrej':
        return 'cutout' if has_nocr else None
    return 'cutout_cr' if has_cr else ('cutout' if has_nocr else None)


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


def gitignored(paths, ws_path=None):
    """The subset of `paths` that .gitignore excludes, as a set of absolute paths.

    Lives here because BOTH make_masks and make_dataset_subplots need it and make_masks
    already imports make_dataset_subplots -- putting it in either would be a cycle.

    One batched `git check-ignore --stdin` call rather than one per directory: there are a
    few hundred cutout dirs and a subprocess each would dominate the run. Any failure --
    git missing, not a repository, an unexpected exit code -- returns the EMPTY set, i.e.
    "nothing is ignored", so the guard can only ever make a run narrower when it works and
    never blocks one when it cannot. (`git check-ignore` exits 1 for "nothing matched",
    which is a normal answer, not an error; only >1 is a real failure.)
    """
    import subprocess
    paths = [os.path.abspath(p) for p in paths]
    if not paths:
        return set()
    root = ws_path or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    try:
        proc = subprocess.run(['git', '-C', root, 'check-ignore', '--stdin'],
                              input='\n'.join(paths), capture_output=True, text=True)
    except (OSError, ValueError):
        return set()
    if proc.returncode > 1:
        return set()
    return {os.path.abspath(line) for line in proc.stdout.splitlines() if line.strip()}
