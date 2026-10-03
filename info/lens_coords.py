"""Catalogue lens positions for every sample, merged into one lookup.

Each sample keeps its own table, transcribed from its source paper: slacs_coords.py
(Bolton et al. 2008), gallery_coords.py (Shu et al. 2016) and bells_coords.py (Brownstein
et al. 2012). All are keyed on the short J-name and the names don't collide, so scripts
resolve any lens through this one dict instead of importing a per-sample table. A new
sample's table is added to _TABLES; a name collision is an error, not a silent overwrite.
"""

from slacs_coords import slacs_coords
from gallery_coords import gallery_coords
from bells_coords import bells_coords

_TABLES = {
    'slacs_coords.py': slacs_coords,
    'gallery_coords.py': gallery_coords,
    'bells_coords.py': bells_coords,
}

lens_coords = {}
for _source, _table in _TABLES.items():
    for _lens, _coord in _table.items():
        if _lens in lens_coords:
            raise ValueError(f'{_lens} appears in more than one info/*_coords.py table '
                             f'(second time in {_source})')
        lens_coords[_lens] = _coord
