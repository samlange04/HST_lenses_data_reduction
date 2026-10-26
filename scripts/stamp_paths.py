"""Stamp paths and filenames, shared by the drizzle, stamp, mask, PSF and mosaic scripts.

**One tree, the size in the filename.** Every science stamp lives in

    data/stamps/<sample>/<lens>/<band>/

and its filenames start `stamp<S>_`, where S is the stamp size in arcsec (`stamp12_` for the
pipeline's standard 12" stamps, `stamp20_` for a 20" re-cut). Stamps of different sizes sit
side by side in the same band directory and cannot overwrite each other, because the size
is part of every name. Only the 12" set is tracked in git; other sizes are gitignored by
pattern.

**Each lens has one MODELLING size**: its `stamp_size` in info/samples/lens_samples.json, else
DEFAULT_SIZE. Every tool's `--size` defaults to it (`size=None` below), so a lens given smaller
stamps needs no flag anywhere; .gitignore must then track its `stamp<S>_*` too. Before
2026-10-26 the size lived in the directory (`data/cutouts/` for 12",
`data/cutouts_20arcsec/` otherwise) and the files were named `cutout_cr_*` / `cutout_*`.

    stamp12_sci.fits            science image
    stamp12_noise.fits          1-sigma noise map
    stamp12.png                 3-panel QC render (gitignored, rebuilt by make_stamps.py)
    stamp12_mask.fits           contaminant mask (hand-drawn, make_masks.py)
    stamp12_mask_arcs.fits      optional arc mask (make_arc_masks.py; gitignored)
    stamp12_positions.json      lensed-image positions (hand-marked, make_positions.py)
    stamp12_dataset.png         sci/noise/PSF/mask QC panel (make_dataset_subplots.py)
    psf.fits                    modelling PSF kernel (trimmed)
    psf_err.fits                its per-pixel error map
    psf_analytic[_err].fits     model tier only: the pre-broadening analytic kernel
    psf_injected.fits           empirical tier, validation only: the injected comparison

**The PSF carries no size and no pass.** The kernel is trimmed by amplitude (AGENTS.md,
*PSF generation*), so it is a property of the band rather than of the stamp it is
convolved with, and the point-source PSF is the same in the CR-rejected and no-CR passes.
Every stamp size in the band directory pairs with the one `psf.fits`.

**The drizzle pass is not in the name.** The science stamp of a band is cut from the
CR-rejected pass where one exists (ACS, WFPC2, WFC3/UVIS) and from the no-CR pass where none
does (WFC3/IR F160W); either way it is `stamp<S>_*`. Which pass it came from is recorded in
the stamp header (`DRZPASS`, written since 2026-10-26; `stamp_pass()` infers it for older
stamps) and in info/provenance/lens_stamp_qc.json. Only a no-CR comparison stamp cut beside a
CR science stamp (`make_stamps.py --pass nocrrej` on a band that has a CR pass) is tagged in
the name, `stamp<S>_nocrrej_*`, so it can never be mistaken for the science stamp.

**A variant reduction is tagged in the name too.** `variant` ('bcfill', 'crfill',
'bcfill_crfill', 'drop', ...) keys a whole alternate drizzle (drizzle_acs_wfc.py --bcfill
etc.; AGENTS.md, *Bad-column fill*), so `drizzled_root` is variant-keyed:

    standard   data/drizzled/          work dir data/drizzle_files/
    bcfill     data/drizzled_bcfill/   work dir data/drizzle_files_bcfill/

and a stamp cut from it carries the variant after the size: `stamp12_bcfill_sci.fits`,
`stamp12_bcfill_crfill_sci.fits`. Only the sci/noise stamps and their PNG are tagged; the
mask, positions and PSF describe the band's grid, which every reduction shares. The
invariant is ONE reduction's stamp per band and size: make_stamps.py refuses to write a
second beside an existing one unless --force, and every reader resolves the stamp through
`find_stamp` rather than spelling the filename, which raises AmbiguousStampError if two
reductions ever share a directory. The tag lives in the name, not only in a header card,
because a merge or cherry-pick between branches carrying different reductions would
otherwise replace one's bytes with the other's silently.
"""
import os

# The pipeline's standard stamp size, in arcsec: the tracked, hand-masked set.
DEFAULT_SIZE = 12.0

# The reduction axes a drizzle/cut can be run with, in the order their tags compose:
# --bcfill --crfill gives 'bcfill_crfill' (and data/drizzled_bcfill_crfill/), never
# 'crfill_bcfill'. make_stamps.py builds its variant from this; find_stamp parses with it.
VARIANT_COMPONENTS = ('bcfill', 'crfill', 'drop')
STAMP_KINDS = ('sci', 'noise')
NOCR_TAG = 'nocrrej'   # tags a no-CR comparison stamp cut beside a CR science stamp

# Filenames of the size- and pass-independent PSF products in a band directory.
PSF_FILES = {
    'psf': 'psf.fits',
    'psf_err': 'psf_err.fits',
    'psf_analytic': 'psf_analytic.fits',
    'psf_analytic_err': 'psf_analytic_err.fits',
    'psf_injected': 'psf_injected.fits',
}


def size_label(size=DEFAULT_SIZE):
    """'12' for 12.0, '20' for 20.0, '7.5' for 7.5 -- the size as it appears in names."""
    return f'{float(size):g}'


def modelling_size(lens, sample=None):
    """The lens's modelling stamp side in arcsec: its lens_samples.json stamp_size, else
    DEFAULT_SIZE. A lens the roster does not know (a scratch dir) gets DEFAULT_SIZE."""
    import mast_target_names
    try:
        return float(mast_target_names.stamp_size(lens, sample) or DEFAULT_SIZE)
    except (KeyError, ValueError, SystemExit):
        return DEFAULT_SIZE


def resolve_size(stamp_dir, size=None):
    """`size`, or the modelling size of the lens whose band dir this is (data/stamps/<sample>/
    <lens>/<band>) when size is None."""
    if size is not None:
        return float(size)
    lens_dir = os.path.dirname(os.path.normpath(stamp_dir))
    return modelling_size(os.path.basename(lens_dir), os.path.basename(os.path.dirname(lens_dir)))


def size_tag(size):
    """'' for each lens's modelling size (pass None), else '_<size>arcsec'. Keys the
    size-dependent things that are not stamps: the QC mosaics tree and the stamp-QC JSON."""
    return '' if size is None else f'_{size_label(size)}arcsec'


def variant_tag(variant):
    """'' for the standard reduction, else '_<variant>' (e.g. '_bcfill')."""
    return f'_{variant}' if variant else ''


def stamp_prefix(size=DEFAULT_SIZE, comparison=False):
    """'stamp12' -- the prefix of a band's science stamp at `size` and of everything drawn
    on its grid (mask, arc mask, positions, dataset PNG). `comparison=True` gives the
    tagged no-CR comparison prefix, 'stamp12_nocrrej'."""
    return f'stamp{size_label(size)}' + (f'_{NOCR_TAG}' if comparison else '')


def prefix_size(prefix):
    """12.0 from 'stamp12' or 'stamp12_nocrrej'."""
    return float(prefix[len('stamp'):].split('_')[0])


def prefixes(size=DEFAULT_SIZE):
    """(science, comparison) prefixes at `size`: ('stamp12', 'stamp12_nocrrej')."""
    return stamp_prefix(size), stamp_prefix(size, comparison=True)


PREFIXES = prefixes()   # at the default size, science first


def psf_name(kind='psf'):
    """'psf.fits', 'psf_err.fits', ... -- see PSF_FILES."""
    return PSF_FILES[kind]


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
    """'stamp12_bcfill_sci.fits' -- the tagged stamp filename. `kind` is 'sci' or 'noise';
    the 3-panel QC PNG is stamp_png_name."""
    return f'{prefix}{variant_tag(variant)}_{kind}.fits'


def stamp_png_name(prefix, variant=''):
    return f'{prefix}{variant_tag(variant)}.png'


class AmbiguousStampError(RuntimeError):
    """More than one reduction's stamp sits in a band dir -- see the module docstring."""


def list_stamps(stamp_dir, prefix, kind='sci'):
    """Every `{prefix}[_<variant>]_<kind>.fits` in stamp_dir, as {variant: path}.

    Matches on the canonical variant grammar, so 'stamp12' does not swallow
    'stamp12_nocrrej_*' (or 'stamp120_*'), and a stray file with an unknown middle is not
    mistaken for a stamp.
    """
    import glob
    found = {}
    for path in glob.glob(os.path.join(stamp_dir, f'{prefix}*_{kind}.fits')):
        base = os.path.basename(path)
        middle = base[len(prefix):-len(f'_{kind}.fits')]
        if middle == '':
            found[''] = path
        elif middle.startswith('_') and _is_variant(middle[1:]):
            found[middle[1:]] = path
    return found


def find_stamp(stamp_dir, prefix, kind='sci', variant=None):
    """The one `{prefix}[_<variant>]_<kind>.fits` in stamp_dir, or None if there is none.

    With `variant=None` (the normal case) the band must carry exactly ONE reduction's
    stamp; two or more raise AmbiguousStampError naming them, because choosing silently is
    precisely the failure the tagged names exist to prevent. Pass `variant` to ask for a
    specific one (returns None if absent).
    """
    stamps = list_stamps(stamp_dir, prefix, kind)
    if variant is not None:
        return stamps.get(variant)
    if len(stamps) > 1:
        names = ', '.join(os.path.basename(p) for _, p in sorted(stamps.items()))
        raise AmbiguousStampError(
            f'{stamp_dir}: {len(stamps)} {prefix}_{kind} stamps from different reductions '
            f'({names}). One stamp per band: remove the one that does not belong '
            f'(or re-cut with --force), then re-run.')
    return next(iter(stamps.values()), None)


def stamp_variant(path):
    """'bcfill' from '.../stamp12_bcfill_sci.fits', '' for an untagged stamp."""
    base = os.path.basename(path)
    if not base.startswith('stamp'):
        raise ValueError(f'{base} is not a stamp filename')
    for kind in STAMP_KINDS:
        if base.endswith(f'_{kind}.fits'):
            parts = base[:-len(f'_{kind}.fits')].split('_')[1:]   # drop 'stamp<S>'
            if parts[:1] == [NOCR_TAG]:
                parts = parts[1:]
            tag = '_'.join(parts)
            if _is_variant(tag):
                return tag
    raise ValueError(f'{base} is not a stamp filename')


def stamp_pass(sci_path):
    """'cr' or 'nocrrej': the drizzle pass a science stamp was cut from.

    Read from the DRZPASS header card where make_stamps.py wrote one (2026-10-26 on).
    Older stamps are inferred from the name and the instrument: a `_nocrrej` tag means
    no-CR, and otherwise the science stamp is the CR pass on every instrument except
    WFC3/IR, which has no CR pass (AGENTS.md, *Cosmic-ray rejection*).
    """
    from astropy.io import fits
    hdr = fits.getheader(sci_path)
    if hdr.get('DRZPASS'):
        return str(hdr['DRZPASS']).strip()
    if f'_{NOCR_TAG}_' in os.path.basename(sci_path):
        return 'nocrrej'
    inst = (str(hdr.get('INSTRUME', '')).strip().upper(), str(hdr.get('DETECTOR', '')).strip().upper())
    return 'nocrrej' if inst == ('WFC3', 'IR') else 'cr'


def find_prefix(stamp_dir, drizzle_pass='auto', size=None):
    """The stamp prefix to use for this band at `size` (None: the lens's modelling size), or
    None if there is no such stamp.

    'auto' (the normal case) returns the science stamp's prefix, 'stamp<S>', falling back to
    a lone comparison stamp 'stamp<S>_nocrrej'. 'cr' and 'nocrrej' ask for a specific pass:
    the science stamp if it was cut from that pass (`stamp_pass`), else the tagged
    comparison stamp for 'nocrrej'. Raises AmbiguousStampError if a prefix has stamps from
    two reductions (find_stamp).
    """
    size = resolve_size(stamp_dir, size)
    science, comparison = stamp_prefix(size), stamp_prefix(size, comparison=True)
    sci = find_stamp(stamp_dir, science, 'sci')
    has_cmp = find_stamp(stamp_dir, comparison, 'sci') is not None
    if drizzle_pass == 'auto':
        return science if sci else (comparison if has_cmp else None)
    if sci and stamp_pass(sci) == drizzle_pass:
        return science
    if drizzle_pass == NOCR_TAG and has_cmp:
        return comparison
    return None


def drizzled_root(ws_path, variant=''):
    """data/drizzled[_<variant>] -- the tree holding the drizzled mosaics.

    A mosaic is size-independent, but it IS variant-keyed: the bad-column-filled re-drizzle
    is a different mosaic, so it lands in data/drizzled_bcfill/ alongside data/drizzled/
    instead of overwriting it. Both feed the one stamp tree; see the module docstring.
    """
    return os.path.join(ws_path, 'data', f'drizzled{variant_tag(variant)}')


def stamps_root(ws_path):
    """data/stamps -- the one science-stamp tree, every size."""
    return os.path.join(ws_path, 'data', 'stamps')


def mosaics_root(ws_path, size=None):
    """data/mosaics[_<S>arcsec] -- the QC mosaics tiling the stamps; None: each lens's
    modelling size."""
    return os.path.join(ws_path, 'data', f'mosaics{size_tag(size)}')


def qc_json_path(ws_path, size=None):
    """info/provenance/lens_stamp_qc[_<S>arcsec].json; None: the modelling-size record.

    A separate file per size, not a deeper nesting inside one: info_json.update is
    {sample: {lens: {key: value}}} throughout info/, and the per-stamp diagnostics it
    records (weight_uniformity in particular) are measured over the stamp region, so they
    genuinely differ between sizes and must not overwrite each other. Only the modelling-size
    file is tracked.
    """
    return os.path.join(ws_path, 'info', 'provenance', f'lens_stamp_qc{size_tag(size)}.json')


def stamp_dir(ws_path, sample, lens, filt):
    """The one band dir: data/stamps/<sample>/<lens>/<filt>. Every stamp size, the PSF and
    everything drawn on the stamps live here."""
    return os.path.join(stamps_root(ws_path), sample, lens, filt)


# The PSF lives in the same band dir as every stamp size; kept as a name so the call sites
# that write or read PSF products say what they are after (distinct from data/psf/, the
# archival full-kernel tree those scripts call `psf_dir`).
psf_stamp_dir = stamp_dir


def gitignored(paths, ws_path=None):
    """The subset of `paths` that .gitignore excludes, as a set of absolute paths.

    Lives here because BOTH make_masks and make_dataset_subplots need it and make_masks
    already imports make_dataset_subplots -- putting it in either would be a cycle.

    One batched `git check-ignore --stdin` call rather than one per directory: there are a
    few hundred band dirs and a subprocess each would dominate the run. Any failure --
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
