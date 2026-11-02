"""
Read and write `{prefix}_positions.json` files that hold one OR MORE position sets.

A position set is the lensed images of ONE source (or one source plane): the multiple images
a positions likelihood (`al.PositionsLH`) asks a mass model to map back to one point. Most
lenses have one. A double source plane, or two distinct sources behind the same deflector,
needs one set each, because their images trace back to DIFFERENT source-plane points
(2026-11-01). make_positions.py marks them with the number keys: 1 = source 1 (the default),
2 = source 2, and so on.

THE FORMAT IS BACKWARD COMPATIBLE BY CONSTRUCTION:
  * one set is written EXACTLY as before: a single `Grid2DIrregular` via `al.output_to_json`,
    so every existing file is a valid one-set file and every existing reader (`al.from_json`)
    still reads a one-set file the same way;
  * two or more sets are written as autoconf's native list type (`{"type": "list", "values":
    [Grid2DIrregular, ...]}`), which `al.from_json` turns back into a Python list.
`read_sets` always returns a list of sets, so a reader never has to care which one it got.
A deliberately empty file (no images marked, e.g. counter-images below the noise) stays a
single empty Grid2DIrregular: one set with no positions.

The same helpers are proposed upstream next to `al.Clicker` (PyAutoGalaxy
`feature/clicker-position-sets`); once released, these can import from there.
"""
import autolens as al
import numpy as np


def read_sets(path):
    """[[(y, x), ...], ...] arcsec, one inner list per position set, in file order."""
    obj = al.from_json(file_path=path)
    grids = obj if isinstance(obj, list) else [obj]
    return [[(float(y), float(x)) for (y, x) in np.asarray(g).reshape(-1, 2)] for g in grids]


def write_sets(sets, path):
    """Write position sets; one set is written in the original single-grid format."""
    # A (n, 2) float array, (0, 2) when empty: an empty Python list would serialise as
    # autoconf's "list" type rather than the "ndarray" every existing empty file holds.
    sets = [np.asarray(list(s), dtype=float).reshape(-1, 2) for s in sets] or [np.zeros((0, 2))]
    if len(sets) == 1:
        al.output_to_json(obj=al.Grid2DIrregular(values=sets[0]), file_path=path)
    else:
        al.output_to_json(obj=[al.Grid2DIrregular(values=s) for s in sets], file_path=path)


def flatten(sets):
    """Every position of every set, in order (for displays that do not care which set)."""
    return [p for s in sets for p in s]
