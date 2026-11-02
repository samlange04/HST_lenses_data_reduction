"""Lens samples, their directory layout, and MAST target-name resolution.

``info/samples/lens_samples.json`` is the single source of truth for *which lenses exist in which
sample*, *where each lens's data live*, and the per-lens MAST quirks. This module is the only
reader of it; nothing else should parse that file or keep a second copy of a lens list. Six
samples are defined (catalogue positions for all of them are merged in
``info/coords/lens_coords.py``), named after their surveys since 2026-11-02:

  ``slacs``    SLACS (Bolton et al. 2008; Auger et al. 2009 grade A) -- the default ``--sample``.
  ``swells``   SWELLS disk-deflector lenses (Treu et al. 2011; Brewer et al. 2012).
  ``gallery``  BELLS GALLERY (Shu et al. 2016), WFC3/UVIS.
  ``bells``    BELLS (Brownstein et al. 2012), ACS F814W.
  ``dinos``    SL2S lenses with the deep GO-17130 F475X (Project Dinos II, Sheu et al. 2025).
  ``sl2s``     the other SL2S lenses (Sonnenfeld et al. 2013): WFPC2 F606W or the snapshots.

**Every lens sits in ``<sample>/<band_mode>/<category>/<lens>``** under every data layer
(``data/calibrated``, ``data/drizzle_files``, ``data/drizzled``, ``data/psf``, ``data/stamps``;
``stamp_paths.lens_dir``). Two per-lens keys, both REQUIRED, set it:

  ``band_mode``   ``single_band`` or ``multi_band``: the number of distinct FILTERS reduced
                  (a split-visit ``f606W_visit2`` is not a second band). Stored, not derived,
                  so reducing a new band never moves a lens silently: ``check_layout.py``
                  reports a lens whose bands on disk no longer match.
  ``category``    ``prime``, ``complex`` (harder morphology, blending, low S/N, ...),
                  ``limited`` (positions deliberately empty, or no lens or source redshift)
                  or ``multi_galaxy`` (more than one deflector galaxy). Precedence
                  multi_galaxy > limited > complex > prime.

Per-lens entries carry only what is *unusual* about a lens; an empty ``{}`` means "no
quirks, query MAST as ``SDSS<LENS>``" (beyond the two layout keys). Eight keys are understood:

  ``mast_target``  the lens is archived under a ``GAL-<plate>-<mjd>-<fiber>`` designation
                   rather than ``SDSS<LENS>``, so the default query returns nothing.
                   ``target_patterns()`` tries this name first, then the SDSS fallback.
  ``force_copy``   the lens's *non*-COPY observations are unusable (e.g. all EXPTIME=0),
                   so the ``-COPY`` duplicates must be used instead. Normally non-COPY is
                   preferred. See ``force_copy()``.
  ``proposals``    ``{filt: [proposal_id, ...]}``: for that filter, keep only observations
                   from these HST programs. SL2S uses it to take F475X from the deep
                   GO-17130 (Project Dinos II) visit alone where an older GO-11588 snapshot
                   also exists. See ``allowed_proposals()``.
  ``archive_proposals`` ``{filt: [proposal_id, ...]}``: for that filter, ignore the
                   target-name query and take the lens's data from these programs, keeping
                   the visits whose footprint covers the lens (``mast_archive.py``). For
                   survey data filed under other target names (CANDELS 'ANY', COSMOS ...).
  ``stamp_box``   side (px) of make_stamps.py's brightest-pixel search box around the
                   catalogue position, overriding the default 100. For a lens whose arc
                   outshines the deflector in its centre band (SL2S J0904-0059 in F606W).
  ``stamp_size``   the lens's MODELLING stamp side in arcsec, overriding the default 12".
                   Every tool's --size defaults to it (stamp_paths.modelling_size), and its
                   stamp<S>_* files need a !data/stamps/**/stamp<S>_* line in .gitignore.
                   A sample-level ``stamp_size`` sets it for every lens of a sample; a
                   per-lens value wins. The multi_galaxy lenses carry 30" (2026-10-30).
  ``stamp_center`` ``["HH:MM:SS.sss", "+DD:MM:SS.ss"]``: the sky position every band's
                   stamp is centred on, skipping make_stamps.py's peak search. For a
                   multi-galaxy deflector, the light-weighted centroid of its member galaxies
                   (the multi_galaxy lenses, 2026-10-29), so the whole group sits in the stamp.
  ``note``         free text recording why the lens is in the sample; not read by code.

Output and directory names always stay in the J convention (the lens keys) regardless of
the MAST name used to fetch.

**Every lens in a sample is meant to be tried on every run.** A lens with no data for a
given instrument/filter is not an error: the drizzle scripts record ``null`` in the
tracking JSONs and exit 0, so a batch runner sweeps the whole sample and only downloads
where MAST actually has something. For SLACS the naming has been surveyed: the 14 GAL-*
overrides are all on Etherington Gold lenses, and every other SLACS lens resolves under
``SDSS<LENS>``, so "no observations" is trustworthy there. ``gallery`` and SL2S carry their
full designations as ``mast_target``.
"""

import json
import os

_INFO_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'info')
SAMPLES_JSON = os.path.join(_INFO_DIR, 'samples', 'lens_samples.json')

DEFAULT_SAMPLE = 'slacs'
BAND_MODES = ('single_band', 'multi_band')
CATEGORIES = ('prime', 'complex', 'limited', 'multi_galaxy')


class NoMastData(Exception):
    """Raised inside a drizzle script's download block when the MAST query returned zero
    observations for this lens+instrument+filter.

    This is an *ordinary* outcome, not a failure -- every lens in a sample is tried on
    every run and most have no data in most bands. It is a distinct exception type purely
    so the download block's broad ``except Exception`` cannot swallow it and report an
    empty archive as a broken download; the caller records ``null`` in the tracking JSONs
    and exits 0.
    """


_cache = None


def _load():
    global _cache
    if _cache is None:
        with open(SAMPLES_JSON) as f:
            _cache = json.load(f)
    return _cache


def samples():
    """Names of all defined samples, in file order."""
    return list(_load().keys())


def _sample(sample):
    data = _load()
    try:
        return data[sample]
    except KeyError:
        raise KeyError(
            f'Unknown sample {sample!r}. Defined samples: {", ".join(data)}. '
            f'Add it to {SAMPLES_JSON} rather than hardcoding a lens list.'
        ) from None


def lenses(sample=DEFAULT_SAMPLE):
    """Sorted lens names in ``sample`` (J convention). Raises KeyError if unknown."""
    return sorted(_sample(sample)['lenses'])


def sample_of(lens):
    """The sample a lens belongs to, or None. Ambiguous only if a lens were listed twice,
    which the generator forbids."""
    for name, block in _load().items():
        if lens in block['lenses']:
            return name
    return None


def _entry(lens, sample=None):
    """Per-lens quirk dict. Searches ``sample`` if given, else every sample."""
    if sample is not None:
        return _sample(sample)['lenses'].get(lens, {})
    for block in _load().values():
        if lens in block['lenses']:
            return block['lenses'][lens]
    return {}


def band_mode(lens, sample=None):
    """``single_band`` or ``multi_band`` (lens_samples.json). Raises KeyError if unset."""
    mode = _entry(lens, sample).get('band_mode')
    if mode not in BAND_MODES:
        raise KeyError(f'{lens}: band_mode {mode!r} not one of {BAND_MODES} in {SAMPLES_JSON}')
    return mode


def category(lens, sample=None):
    """``prime``, ``complex``, ``limited`` or ``multi_galaxy`` (lens_samples.json)."""
    cat = _entry(lens, sample).get('category')
    if cat not in CATEGORIES:
        raise KeyError(f'{lens}: category {cat!r} not one of {CATEGORIES} in {SAMPLES_JSON}')
    return cat


def lens_rel(lens, sample=None):
    """``<sample>/<band_mode>/<category>/<lens>``: the lens's path below any data layer."""
    owner = sample or sample_of(lens)
    if owner is None:
        raise KeyError(f'{lens} is in no sample of {SAMPLES_JSON}')
    return os.path.join(owner, band_mode(lens, owner), category(lens, owner), lens)


LAYOUT_TAGS = {'single_band': 'SB', 'multi_band': 'MB',
               'prime': 'P', 'complex': 'C', 'limited': 'L', 'multi_galaxy': 'MG'}


def layout_tag(lens, sample=None):
    """Short sub-type tag for QC panels, e.g. 'MB-P' (multi-band, prime) or 'SB-L'."""
    return f'{LAYOUT_TAGS[band_mode(lens, sample)]}-{LAYOUT_TAGS[category(lens, sample)]}'


def force_copy(lens, sample=None):
    """True if this lens must be fetched from its ``-COPY`` MAST observations (its
    non-COPY frames are unusable, e.g. EXPTIME=0)."""
    return bool(_entry(lens, sample).get('force_copy', False))


def uvis_filters(sample):
    """WFC3/UVIS filters run_gallery_uvis_all.sh drizzles for this sample. A sample-level
    ``uvis_filters`` list overrides the gallery default. SL2S sets [f475X, f600LP] so the UVIS
    runner never touches the F606W/F814W bands that its lenses take from WFPC2/ACS (one
    instrument per filter key, see info_json.band_owner_conflict)."""
    return list(_sample(sample).get('uvis_filters', ['f606W', 'f814W', 'f438W']))


def center_bands(sample):
    """Ordered bands whose mosaic may define a lens's shared stamp centre: make_stamps.py
    uses the first one drizzled for that lens. A sample-level ``center_bands`` list in
    lens_samples.json overrides the default ``["f814W"]`` (which keeps the old behaviour: a
    lens without F814W centres each band on itself). The SL2S samples set
    f814W > f606W > f600LP > f475X, because most have no F814W and F475X shows little
    deflector light."""
    return list(_sample(sample).get('center_bands', ['f814W']))


def allowed_proposals(lens, filt, sample=None):
    """HST program IDs (strings) this lens's ``filt`` observations are restricted to, or
    None for no restriction. ``filt`` is the product key as passed to ``--filt``."""
    props = _entry(lens, sample).get('proposals', {}).get(filt)
    return [str(p) for p in props] if props else None


def archive_proposals(lens, filt, sample=None):
    """HST program IDs (strings) to take ``filt`` from by footprint, or None."""
    props = _entry(lens, sample).get('archive_proposals', {}).get(filt)
    return [str(p) for p in props] if props else None


def stamp_box(lens, sample=None):
    """Per-lens stamp search-box override (px), or None."""
    box = _entry(lens, sample).get('stamp_box')
    return int(box) if box else None


def stamp_size(lens, sample=None):
    """Modelling stamp side (arcsec): the lens's own ``stamp_size``, else its sample's,
    else None for the default."""
    size = _entry(lens, sample).get('stamp_size')
    if not size:
        owner = sample or sample_of(lens)
        size = _sample(owner).get('stamp_size') if owner else None
    return float(size) if size else None


def stamp_center(lens, sample=None):
    """Per-lens fixed stamp centre as (RA, Dec) sexagesimal strings, or None."""
    c = _entry(lens, sample).get('stamp_center')
    return tuple(c) if c else None


def target_patterns(lens, sample=None):
    """Ordered list of MAST ``target_name`` wildcard patterns to query for ``lens``.

    For a lens with a known ``GAL-*`` designation the GAL name is tried first, then the
    default ``SDSS<LENS>`` pattern; otherwise only the ``SDSS<LENS>`` pattern is returned.
    An unknown lens is not an error -- it falls through to the SDSS pattern, so a lens can
    be queried before it is added to a sample.
    """
    patterns = []
    gal = _entry(lens, sample).get('mast_target')
    if gal:
        patterns.append(f'{gal}%')
    patterns.append(f'SDSS{lens}%')
    return patterns


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description='Print the lenses in a sample, one per line.')
    p.add_argument('sample', nargs='?', default=DEFAULT_SAMPLE,
                   help=f'sample name (default {DEFAULT_SAMPLE}); '
                        f'"--list" prints the defined samples')
    p.add_argument('--list', action='store_true', help='list sample names and sizes')
    p.add_argument('--uvis-filters', action='store_true',
                   help='print the sample\'s WFC3/UVIS filters (space-separated) for the UVIS runner')
    p.add_argument('--rel', metavar='LENS',
                   help='print LENS\'s <sample>/<band_mode>/<category>/<lens> path, for the shell runners')
    p.add_argument('--print-sample', action='store_true',
                   help='print only the resolved sample name, so a shell runner can echo '
                        'the default without hardcoding a second copy of it')
    a = p.parse_args()
    if a.list:
        for s in samples():
            print(f'{s:10s} {len(lenses(s)):3d} lenses')
    elif a.rel:
        print(lens_rel(a.rel, a.sample))
    elif a.uvis_filters:
        print(' '.join(uvis_filters(a.sample)))
    elif a.print_sample:
        _sample(a.sample)          # validate before printing
        print(a.sample)
    else:
        for lens in lenses(a.sample):
            print(lens)
