#!/usr/bin/env python
"""
Write info/manifest.csv: one row per lens and band, the index users filter to find data.

The directory layout answers one question (data/stamps/<sample>/<band_mode>/<category>/<lens>/
<band>/); the manifest answers the rest -- "every SLACS F606W band whatever its category",
"every lens with a source redshift", "every band with an empirical PSF" -- and is what
make_release_zips.py builds its per-filter archives from. Regenerate it after any change to
the roster, the stamps, the positions or the PSFs; it is tracked, so a stale one shows in the
diff.

Columns: sample, band_mode, category, category_reason, lens, band (the product directory, so
a split visit is its own row, f606W_visit2), filter (the band without a visit suffix),
instrument, exptime_s, stamp_size_arcsec, path (the band directory, repo-relative),
n_positions, n_position_sets, psf_method, zlens, zsrc, theta_e_arcsec, etherington_gold.

Usage:
    uv run python scripts/make_manifest.py
"""
import csv
import glob
import os
import sys

ws_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ws_path, 'scripts'))
import info_json  # noqa: E402
import mast_target_names  # noqa: E402
import positions_io  # noqa: E402
import stamp_paths  # noqa: E402

OUT = os.path.join(ws_path, 'info', 'manifest.csv')
COLUMNS = ['sample', 'band_mode', 'category', 'category_reason', 'lens', 'band', 'filter',
           'instrument', 'exptime_s', 'stamp_size_arcsec', 'path', 'n_positions',
           'n_position_sets', 'psf_method', 'zlens', 'zsrc', 'theta_e_arcsec',
           'etherington_gold']


def main():
    prov = os.path.join(ws_path, 'info', 'provenance')
    instrument = info_json.load(os.path.join(prov, 'lens_instrument.json'))
    exptime = info_json.load(os.path.join(prov, 'lens_exptime.json'))
    psf = info_json.load(os.path.join(prov, 'lens_psf.json'))
    theta = info_json.load(os.path.join(ws_path, 'info', 'catalogues', 'lens_einstein_radii.json'))
    rows = []
    for sample in mast_target_names.samples():
        for lens in mast_target_names.lenses(sample):
            entry = mast_target_names._entry(lens, sample)
            ldir = stamp_paths.lens_dir(ws_path, 'stamps', sample, lens)
            size = stamp_paths.modelling_size(lens, sample)
            cat = theta.get(lens, {})
            for band_dir in sorted(glob.glob(os.path.join(ldir, '*'))):
                band = os.path.basename(band_dir)
                prefix = stamp_paths.find_prefix(band_dir, 'auto', size)
                if prefix is None:
                    continue
                pos = os.path.join(band_dir, f'{prefix}_positions.json')
                sets = positions_io.read_sets(pos) if os.path.exists(pos) else []
                rows.append({
                    'sample': sample, 'band_mode': entry['band_mode'],
                    'category': entry['category'],
                    'category_reason': entry.get('category_reason', ''),
                    'lens': lens, 'band': band, 'filter': band.split('_')[0],
                    'instrument': (instrument.get(sample, {}).get(lens) or {}).get(band) or '',
                    'exptime_s': (exptime.get(sample, {}).get(lens) or {}).get(band) or '',
                    'stamp_size_arcsec': f'{size:g}',
                    'path': os.path.relpath(band_dir, ws_path),
                    'n_positions': len(positions_io.flatten(sets)) if os.path.exists(pos) else '',
                    'n_position_sets': len([s for s in sets if s]) if os.path.exists(pos) else '',
                    'psf_method': ((psf.get(sample, {}).get(lens) or {}).get(band) or {}).get('method', ''),
                    'zlens': cat.get('zlens', ''), 'zsrc': cat.get('zsrc', ''),
                    'theta_e_arcsec': cat.get('theta_e_arcsec', ''),
                    'etherington_gold': 'yes' if entry.get('etherington_gold') else '',
                })
    with open(OUT, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS, lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    print(f'{len(rows)} bands of {len({r["lens"] for r in rows})} lenses -> '
          f'{os.path.relpath(OUT, ws_path)}')


if __name__ == '__main__':
    main()
