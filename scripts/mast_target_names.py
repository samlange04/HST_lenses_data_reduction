"""Lens samples and MAST target-name resolution.

``info/lens_samples.json`` is the single source of truth for *which lenses exist in which
sample* and for the per-lens MAST quirks. This module is the only reader of it; nothing
else should parse that file or keep a second copy of a lens list. Three samples are
defined:

  ``slacs_gold``   the 38 SLACS lenses reduced so far -- the working sample, and the
                   default ``--sample`` of every script.
  ``slacs_other``  the other 93 SLACS lenses from Bolton et al. 2008 Table 4
                   (``info/slacs_coords.py``). Not yet reduced.
  ``gallery``      the 16 BELLS GALLERY lenses (HST proposals 14189, 16734).

Per-lens entries carry only what is *unusual* about a lens; an empty ``{}`` means "no
quirks, query MAST as ``SDSS<LENS>``". Two keys are understood:

  ``mast_target``  the lens is archived under a ``GAL-<plate>-<mjd>-<fiber>`` designation
                   rather than ``SDSS<LENS>``, so the default query returns nothing.
                   ``target_patterns()`` tries this name first, then the SDSS fallback.
  ``force_copy``   the lens's *non*-COPY observations are unusable (e.g. all EXPTIME=0),
                   so the ``-COPY`` duplicates must be used instead. Normally non-COPY is
                   preferred. See ``force_copy()``.

Output and directory names always stay in the J convention (the lens keys) regardless of
the MAST name used to fetch.

**Every lens in a sample is meant to be tried on every run.** A lens with no data for a
given instrument/filter is not an error: the drizzle scripts record ``null`` in the
tracking JSONs and exit 0, so a batch runner sweeps the whole sample and only downloads
where MAST actually has something. Note the asymmetry this leaves: for ``slacs_gold`` the
GAL-* names are known, so "no observations" is trustworthy; for ``slacs_other`` and
``gallery`` no GAL-* names have been surveyed, so a no-data result may just mean the lens
is archived under a name we have not looked up yet.
"""

import json
import os

_INFO_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'info')
SAMPLES_JSON = os.path.join(_INFO_DIR, 'lens_samples.json')

DEFAULT_SAMPLE = 'slacs_gold'


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


def force_copy(lens, sample=None):
    """True if this lens must be fetched from its ``-COPY`` MAST observations (its
    non-COPY frames are unusable, e.g. EXPTIME=0)."""
    return bool(_entry(lens, sample).get('force_copy', False))


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
    p.add_argument('--print-sample', action='store_true',
                   help='print only the resolved sample name, so a shell runner can echo '
                        'the default without hardcoding a second copy of it')
    a = p.parse_args()
    if a.list:
        for s in samples():
            print(f'{s:14s} {len(lenses(s)):3d} lenses')
    elif a.print_sample:
        _sample(a.sample)          # validate before printing
        print(a.sample)
    else:
        for lens in lenses(a.sample):
            print(lens)
