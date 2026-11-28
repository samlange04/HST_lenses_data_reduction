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

import argparse
import difflib
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


def sample_arg(name):
    """argparse ``type`` for --sample: ``name`` if it is a defined sample, else an error naming
    the closest one. Without it a mistyped sample reads as a real sample with no stamps, and a
    GUI reports "no stamps found" instead of the typo. (Shared with the JWST repo's
    ``lens_samples`` -- the profile seam of the merged GUI tools.)"""
    known = samples()
    if name in known:
        return name
    close = difflib.get_close_matches(name.replace('-', '_'), known, n=1, cutoff=0.5)
    raise argparse.ArgumentTypeError(f'unknown sample {name!r}'
                                     + (f" (did you mean '{close[0]}'?)" if close else '')
                                     + f'; defined samples: {", ".join(known)}')


def resolve_sample(parser, args, required=False, default=None):
    """The sample a script works on, from its parsed ``args`` (``args.sample``, ``args.lens``:
    one name or a list), so --sample is never needed alongside --lens.

      * --lens given: the sample is the one the lens is in (lens_samples.json; a lens is in one
        sample only). An unknown lens is an error naming the closest one; --sample, if also
        given, must agree; lenses from different samples are an error.
      * no --lens: --sample, else ``default`` (the script's old default sample; None for "every
        sample"), or an error when ``required``.

    Errors go through ``parser.error``, so they print the usage line and exit 2. The result is
    also written back to ``args.sample``. (The JWST repo's signature and behaviour, adopted
    2026-11-28 as the shared GUI seam; it replaces the old ``resolve_sample(sample, lens)``.)"""
    names = args.lens if isinstance(args.lens, (list, tuple)) else ([args.lens] if args.lens else [])
    sample = args.sample
    if names:
        found = {}
        for lens in names:
            s = sample_of(lens)
            if s is None:
                every = [l for b in _load().values() for l in b['lenses']]
                close = difflib.get_close_matches(lens, every, n=1, cutoff=0.6)
                parser.error(f'unknown lens {lens!r}' + (f" (did you mean '{close[0]}'?)" if close else '')
                             + f'; it is in no sample of {SAMPLES_JSON}')
            found[lens] = s
        if sample is not None:
            wrong = [f'{l} is in {s!r}' for l, s in found.items() if s != sample]
            if wrong:
                parser.error(f'{"; ".join(wrong)}, not --sample {sample!r} (leave --sample out: it '
                             f'is found from --lens)')
        elif len(set(found.values())) > 1:
            parser.error('the lenses are in different samples ('
                         + ', '.join(f'{l}: {s}' for l, s in found.items()) + '); run one sample at a time')
        else:
            sample = next(iter(found.values()))
    elif sample is None:
        if required:
            parser.error('give --lens, or --sample for every lens of a sample')
        sample = default
    args.sample = sample
    return sample


def _entry(lens, sample=None):
    """Per-lens quirk dict. Searches ``sample`` if given, else every sample."""
    if sample is not None:
        return _sample(sample)['lenses'].get(lens, {})
    for block in _load().values():
        if lens in block['lenses']:
            return block['lenses'][lens]
    return {}


# Public name of _entry, for the merged GUI tools (the JWST repo's profile calls it ``entry``).
# Behaviour differs from JWST's on an UNKNOWN lens only: this one returns {}, theirs raises.
entry = _entry


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


def deflector_radec(lens, sample=None):
    """(ra, dec, exact) of the deflector when it is NOT the stamp centre, else None.

    Only the group-centred multi_galaxy stamps (``stamp_center``) qualify: every other stamp
    is cut on the deflector itself, and the isophote fit refines from the stamp centre instead.
    The position is the catalogue one (``info/coords/lens_coords.py``), parsed as plain units
    -- the caller builds the SkyCoord in the stamp WCS's own frame (the FK5/ICRS trap,
    AGENTS.md *Orientation*) -- so ``exact`` is False: the fit should centroid-refine from it
    (a catalogue position can sit ~1" off the core on a group stamp). Profile seam of the
    merged GUI tools (isophote_model.deflector_pixel).
    """
    if not stamp_center(lens, sample):
        return None
    import sys
    coords_dir = os.path.join(_INFO_DIR, 'coords')
    if coords_dir not in sys.path:
        sys.path.insert(0, coords_dir)
    from lens_coords import lens_coords
    if lens not in lens_coords:
        return None
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    c = SkyCoord(*lens_coords[lens], unit=(u.hourangle, u.deg), frame='icrs')
    return (float(c.ra.deg), float(c.dec.deg), False)


RADII_JSON = os.path.join(_INFO_DIR, 'catalogues', 'lens_einstein_radii.json')
PRIMARY_BAND_JSON = os.path.join(_INFO_DIR, 'provenance', 'lens_primary_band.json')

_radii_cache = None


def _radii():
    """The Einstein-radius catalogue (flat by lens; see AGENTS.md *Tracking JSONs*)."""
    global _radii_cache
    if _radii_cache is None:
        with open(RADII_JSON) as f:
            _radii_cache = json.load(f)
    return _radii_cache


def theta_e(lens, sample=None):
    """The lens's measured Einstein radius in arcsec, or None (no row, or a null theta_E --
    Auger+2009 has no lens model for 13 SLACS lenses, and J1037+3517 has no Brewer+2012 model)."""
    v = _radii().get(lens, {}).get('theta_e_arcsec')
    return None if v is None else float(v)


def redshifts(lens, sample=None):
    """(zlens, zsrc, zlens_kind, zsrc_kind) from the catalogue; None where unknown.
    After the 2026-11-06 audit every non-null value is spectroscopic ('spec')."""
    row = _radii().get(lens, {})
    return (row.get('zlens'), row.get('zsrc'), row.get('zlens_kind'), row.get('zsrc_kind'))


def redshift_pair(lens, sample=None):
    """(z_lens, z_source), None where unknown -- the JWST repo's ``redshifts`` arity, for the
    merged GUI tools. (This repo's ``redshifts`` returns the kinds too; the JWST catalogue
    keeps them in the sample file.)"""
    zl, zs, _, _ = redshifts(lens, sample)
    return zl, zs


def redshift_kind_of(kind_text):
    """'spectroscopic' or 'photometric' from a ``zlens_kind``/``zsrc_kind`` string ('spec',
    'photo-z (...)'), else None. Shared with the JWST repo's ``lens_samples``."""
    k = (kind_text or '').strip().lower()
    if k.startswith('spec'):
        return 'spectroscopic'
    if k.startswith('photo'):
        return 'photometric'
    return None


def redshift_kinds(lens, sample=None):
    """(kind_lens, kind_source): 'spectroscopic', 'photometric', or None (no redshift).

    Only a SPECTROSCOPIC redshift is a confirmed redshift (the 2026-11-06 audit; Sam's rule).
    A photometric one may be shown, but always labelled as photometric. A redshift without a
    recognisable kind is an error, not silently 'confirmed'."""
    zl, zs, kl, ks = redshifts(lens, sample)
    out = []
    for z, k in ((zl, kl), (zs, ks)):
        if z is None:
            out.append(None)
            continue
        kind = redshift_kind_of(k)
        if kind is None:
            raise ValueError(f'{lens}: redshift {z} has no recognisable kind ({k!r}); '
                             f'it must start with "spec" or "photo-z"')
        out.append(kind)
    return tuple(out)


def redshift_text(z, kind, fmt='{:g}'):
    """A redshift for display, as precise as recorded: '0.587' (spectroscopic), '0.577 (photo-z)'
    (photometric), '?' (none). Shared with the JWST repo's ``lens_samples``."""
    if z is None:
        return '?'
    return fmt.format(z) + ('' if kind == 'spectroscopic' else ' (photo-z)')


def primary_band(lens, sample=None):
    """(band, source): the lens's PRIMARY band -- the one the mask and arc-mask GUIs draw first
    (the other bands start from its mask), the one positions are clicked on, and the one to
    model first. Ported from JWST_lenses_data_reduction 2026-11-10. `source` says where it
    came from:
      'visual'   ``best_band_visual`` in lens_samples.json, chosen by eye with
                 choose_best_band.py (always wins);
      'arc_snr'  ``best_band_arc_snr``, choose_primary_band.py's measurement in
                 info/provenance/lens_primary_band.json (arc S/N in the ring annulus,
                 discounted where the PSF under-resolves theta_E);
      (None, None) when neither exists, and callers fall back to their fixed band order."""
    e = _entry(lens, sample)
    if e.get('best_band_visual'):
        return e['best_band_visual'], 'visual'
    owner = sample or sample_of(lens)
    try:
        with open(PRIMARY_BAND_JSON) as f:
            rec = json.load(f).get(owner, {}).get(lens, {})
    except FileNotFoundError:
        rec = {}
    if rec.get('best_band_arc_snr'):
        return rec['best_band_arc_snr'], 'arc_snr'
    return None, None


def set_best_band_visual(lens, band, note, sample=None):
    """Record the band chosen BY EYE for `lens` as ``best_band_visual`` in lens_samples.json
    (choose_best_band.py), with a note saying how. Visual inspection only: not necessarily the
    band with the highest arc S/N (that is ``best_band_arc_snr`` in lens_primary_band.json).
    Dumps with ensure_ascii=False so the file keeps its literal unicode."""
    global _cache
    owner = sample or sample_of(lens)
    data = _load()
    e = data[owner]['lenses'][lens]
    e['best_band_visual'] = band
    e['best_band_visual_note'] = note
    with open(SAMPLES_JSON, 'w') as f:
        f.write(json.dumps(data, indent=4, ensure_ascii=False) + '\n')
    _cache = None


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
