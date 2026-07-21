"""MAST target-name resolution for the lens samples.

Most lenses are on MAST under an ``SDSS<LENS>`` target name (e.g. ``SDSSJ0008-0004``).
A subset are archived instead under a ``GAL-<plate>-<mjd>-<fiber>`` designation, so a
``SDSS{lens}%`` query returns no observations for them. ``target_patterns()`` returns the
ordered list of ``target_name`` wildcard patterns to try: the ``GAL-*`` name first (when
known) then the ``SDSS<LENS>`` fallback.

Output/directory names always stay in the J convention (the dict keys) regardless of the
MAST name used to fetch. See the "Non-standard MAST target names" section of CLAUDE.md.
"""

# J-convention output name -> GAL-<plate>-<mjd>-<fiber> MAST target name
GAL_TARGET_NAMES = {
    'J0216-0813': 'GAL-0668-52162-428',
    'J0737+3216': 'GAL-0541-51959-145',
    'J0912+0029': 'GAL-0472-51955-429',
    'J0956+5100': 'GAL-0902-52409-068',
    'J0959+0410': 'GAL-0572-52289-495',
    'J1205+4910': 'GAL-0969-52442-134',
    'J1250+0523': 'GAL-0847-52426-549',
    'J1402+6321': 'GAL-0605-52353-503',
    'J1420+6019': 'GAL-0788-52338-605',
    'J1627-0053': 'GAL-0364-52000-084',
    'J1630+4520': 'GAL-0626-52057-518',
    'J2238-0754': 'GAL-0722-52224-442',
    'J2300+0022': 'GAL-0677-52606-520',
    'J2303+1422': 'GAL-0743-52262-304',
}


# Lenses whose *non*-COPY MAST observations are unusable, so the "-COPY" duplicate
# observations must be used instead. Normally non-COPY is preferred; these are exceptions.
FORCE_COPY_LENSES = {
    'J1032+5322',  # non-COPY F814W ACS FLCs are all EXPTIME=0; the -COPY frames have real exposures
}


def force_copy(lens):
    """True if this lens must be fetched from its ``-COPY`` MAST observations (its
    non-COPY frames are unusable, e.g. EXPTIME=0). See FORCE_COPY_LENSES."""
    return lens in FORCE_COPY_LENSES


def target_patterns(lens):
    """Ordered list of MAST ``target_name`` wildcard patterns to query for ``lens``.

    For a lens with a known ``GAL-*`` designation the GAL name is tried first, then the
    default ``SDSS<LENS>`` pattern; otherwise only the ``SDSS<LENS>`` pattern is returned.
    """
    patterns = []
    if lens in GAL_TARGET_NAMES:
        patterns.append(f'{GAL_TARGET_NAMES[lens]}%')
    patterns.append(f'SDSS{lens}%')
    return patterns
