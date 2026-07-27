"""
WFPC2 / WFC3 PSF MAST-database query + cutout-download helpers.

Vendored (trimmed, tqdm-free) from the STScI hst_notebooks module of the same
name by Fred Dauphin (July 2024):
  spacetelescope/hst_notebooks notebooks/WFC3/mast_api_psf/mast_api_psf.py
We keep only the query + single-file download path used by make_psf.py's WFPC2
F606W model tier, and drop the multiprocessing/tqdm progress-bar downloaders so
the module has no dependency beyond astroquery + requests (both in stenv).

The database is described in Dauphin et al., "The WFPC2 and WFC3 PSF Database"
(ISR WFC3 2021-12). WFPC2's filter field is `filter_1` (not `filter`), and the
`chip` column is the WFPC2 CCD / FITS extension: 1=PC, 2=WF2, 3=WF3, 4=WF4 --
our lens galaxies fall on WF3 (chip 3). See CLAUDE.md 'WFPC2: the lens is on WF3'.
"""

import os
import requests
from astroquery.mast import Mast

REQUEST_URL_PREFIX = 'https://mast.stsci.edu/api/v0.1/Download'

_DETECTOR_SERVICE = {
    'UVIS': 'Wfc3Psf.Uvis',
    'IR': 'Wfc3Psf.Ir',
    'WFPC2': 'Wfpc2Psf.Uvis',
}


def set_filters(parameters):
    """Dict of {column: values} -> MAST filter list. `values` is a list of exact
    matches (strings) or a [{'min':a,'max':b}] range."""
    return [{'paramName': p, 'values': v} for p, v in parameters.items()]


def mast_query_psf_database(detector, filts, columns=('*',)):
    """Query a WFPC2/WFC3 PSF database. Returns an astropy Table.

    For WFPC2 the 'filter' column is auto-mapped to 'filter_1'.
    """
    detector = detector.upper()
    columns = list(columns)
    try:
        database = _DETECTOR_SERVICE[detector]
    except KeyError:
        raise ValueError(f'{detector} not in {list(_DETECTOR_SERVICE)}')
    service = f'Mast.Catalogs.Filtered.{database}'

    if detector == 'WFPC2':
        if 'filter' in columns:
            columns[columns.index('filter')] = 'filter_1'
        for param in filts:
            if param.get('paramName') == 'filter':
                param['paramName'] = 'filter_1'

    cols = '*' if '*' in columns else ','.join(columns)
    return Mast.service_request(service, {'columns': cols, 'filters': filts})


def make_dataURIs(obs, detector, file_suffix=('c0m',), unsat_size=51, sat_size=101):
    """Build fitscut dataURIs for each source. For WFPC2 the cutout is taken from
    the `chip` extension at (x_cal, y_cal). Returns a list of dataURI strings."""
    valid = {'raw', 'd0m', 'flt', 'c0m', 'flc'}
    file_suffix = list(file_suffix)
    for s in file_suffix:
        if s not in valid:
            raise ValueError(f'{s} not in {sorted(valid)}')
    detector = detector.upper()
    instrument = 'WFC3' if detector in ('UVIS', 'IR') else 'WFPC2'
    base = f'mast:{instrument}PSF/url/cgi-bin/fitscut.cgi'

    uris = []
    for row in obs:
        iden = row['id']
        root = row['rootname']
        filt = row['filter_1'] if detector == 'WFPC2' else row['filter']
        chip = row['chip']
        size = unsat_size if row['qfit'] > 0 else sat_size
        if detector == 'UVIS':
            fits_ext = 4 if (str(chip) == '1' and row['subarray'] == 0) else 1
        else:
            fits_ext = chip
        for suffix in file_suffix:
            coord = 'raw' if suffix in ('raw', 'd0m') else 'cal'
            x, y = row[f'x_{coord}'], row[f'y_{coord}']
            read = f'red={root}_{suffix}[{fits_ext}]'
            cut = f'size={size}&x={x}&y={y}&format=fits'
            save = f'{root}_{iden}_{filt}_{suffix}_cutout.fits'
            uris.append(f'{base}?{read}&{cut}/{save}')
    return uris


def download_request_file(dataURI, filename):
    """GET a single cutout dataURI to `filename`. Returns filename."""
    resp = requests.get(f'{REQUEST_URL_PREFIX}/file', params={'uri': dataURI})
    resp.raise_for_status()
    with open(filename, 'wb') as fh:
        fh.write(resp.content)
    return filename
