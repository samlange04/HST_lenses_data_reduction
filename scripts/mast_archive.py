"""Fetch a lens's observations from named HST programs by footprint, not by target name.

Archive data taken by other surveys (CANDELS, COSMOS, ...) is filed under their own target
names -- 'ANY', 'COSMOS33-11', 'UDS-V0R' -- so the drizzle scripts' SDSS/SL2S name query
never sees it. A lens's ``archive_proposals`` quirk in info/lens_samples.json names the
programs to take for a filter; this queries those programs and keeps the observations whose
footprint (``s_region``) contains the lens's catalogue position (info/lens_coords.py).

A position search alone is not used because MAST matches it on footprint centres: on the
SL2S inventory (2026-10-04) a 10" cone missed most calibrated visits that do cover the lens.
"""
import os
import sys

import numpy as np
from matplotlib.path import Path
from astropy.coordinates import SkyCoord
import astropy.units as u

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'info'))
from lens_coords import lens_coords
import mast_target_names


def _polygons(s_region):
    """[(ra, dec) vertex arrays] for every POLYGON in an STC-S s_region string."""
    polys, cur = [], None
    for tok in str(s_region).replace('(', ' ').replace(')', ' ').split():
        if tok.upper() == 'POLYGON':
            cur = []
            polys.append(cur)
        elif cur is not None:
            try:
                cur.append(float(tok))
            except ValueError:   # frame names (ICRS, J2000, ...) between POLYGON and vertices
                pass
    return [np.array(p).reshape(-1, 2) for p in polys if len(p) >= 6]


def footprint_contains(s_region, coord):
    """True if any polygon of ``s_region`` contains ``coord`` (tangent-plane test)."""
    ra0, dec0 = coord.ra.deg, coord.dec.deg
    for xy in _polygons(s_region):
        dx = (xy[:, 0] - ra0 + 180.0) % 360.0 - 180.0
        x, y = dx * np.cos(np.radians(dec0)), xy[:, 1] - dec0
        if Path(np.c_[x, y]).contains_point((0.0, 0.0)):
            return True
    return False


def archive_observations(Observations, lens, filt, instrument, sample=None):
    """The MAST observations for ``lens``/``filt`` from its ``archive_proposals`` programs
    whose footprint covers the lens, or None when the lens has no such quirk for ``filt``
    (the caller then falls back to its target-name query). May return an empty table."""
    props = mast_target_names.archive_proposals(lens, filt, sample)
    if not props:
        return None
    coord = SkyCoord(*lens_coords[lens], unit=(u.hourangle, u.deg))
    obs = Observations.query_criteria(proposal_id=props, obs_collection='HST',
                                      instrument_name=instrument, filters=[filt.upper()])
    keep = np.array([footprint_contains(r['s_region'], coord) for r in obs], dtype=bool)
    print(f'  archive programs {", ".join(props)}: {keep.sum()}/{len(obs)} {instrument} '
          f'{filt.upper()} observations cover {lens}')
    return obs[keep]
