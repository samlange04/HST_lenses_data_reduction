#!/usr/bin/env python
"""
Build the release archives: one zip per sample and one per (sample, filter), from the manifest.

The repo layout is data/stamps/<sample>/<band_mode>/<category>/<lens>/<band>/, so "all of a
sample" is one directory, but "all SLACS F606W, whatever its category" is spread across up to
eight. The per-filter zips gather exactly that (split visits included: f606W_visit2 counts as
F606W), so a user downloads one file per (sample, filter) they want. Every zip keeps the
repo-relative paths, so unzipping several archives (or one into a clone) fills one consistent
tree, and each carries the matching rows of info/manifest.csv as manifest.csv.

Only TRACKED files go in (git ls-files): the science stamps, noise maps, masks, positions,
PSFs and dataset QC PNGs. Gitignored renders are regenerable (scripts/regenerate_files.py) and
stay out. Regenerate info/manifest.csv first (make_manifest.py); this script refuses a band the
manifest lists that has no tracked files.

Output: release/<sample>.zip and release/<sample>_<filter>.zip (gitignored). The Zenodo upload
itself is manual (see the lensprep release plan).

Usage:
    uv run python scripts/make_release_zips.py                       # every sample and filter
    uv run python scripts/make_release_zips.py --sample slacs --filter f606W
    uv run python scripts/make_release_zips.py --dry-run
"""
import argparse
import csv
import io
import os
import subprocess
import sys
import zipfile

ws_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANIFEST = os.path.join(ws_path, 'info', 'manifest.csv')
OUT_DIR = os.path.join(ws_path, 'release')


def tracked_files():
    out = subprocess.run(['git', 'ls-files', 'data/stamps'], cwd=ws_path, check=True,
                         capture_output=True, text=True).stdout
    by_dir = {}
    for f in out.splitlines():
        by_dir.setdefault(os.path.dirname(f), []).append(f)
    return by_dir


def write_zip(path, rows, files_by_dir, header, dry_run):
    files = []
    for r in rows:
        band_files = files_by_dir.get(r['path'], [])
        if not band_files:
            raise SystemExit(f"{r['path']}: in the manifest but no tracked files "
                             f"(regenerate info/manifest.csv?)")
        files += band_files
    size = sum(os.path.getsize(os.path.join(ws_path, f)) for f in files)
    print(f'  {os.path.relpath(path, ws_path):34s} {len(rows):3d} bands  {len(files):5d} files  '
          f'{size / 1e6:7.1f} MB')
    if dry_run:
        return
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=header, lineterminator='\n')
    w.writeheader()
    w.writerows(rows)
    with zipfile.ZipFile(path, 'w', compression=zipfile.ZIP_DEFLATED) as z:
        z.writestr('manifest.csv', buf.getvalue())
        for f in files:
            z.write(os.path.join(ws_path, f), arcname=f)


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--sample', help='one sample (default: all)')
    p.add_argument('--filter', help='one filter, e.g. f606W (default: all, plus the whole-sample zip)')
    p.add_argument('--dry-run', action='store_true', help='list the archives and sizes only')
    a = p.parse_args()

    with open(MANIFEST) as f:
        reader = csv.DictReader(f)
        header, rows = reader.fieldnames, list(reader)
    files_by_dir = tracked_files()
    os.makedirs(OUT_DIR, exist_ok=True)
    for sample in sorted({r['sample'] for r in rows}):
        if a.sample and sample != a.sample:
            continue
        srows = [r for r in rows if r['sample'] == sample]
        if not a.filter:
            write_zip(os.path.join(OUT_DIR, f'{sample}.zip'), srows, files_by_dir, header,
                      a.dry_run)
        filts = sorted({r['filter'] for r in srows})
        for filt in filts:
            if a.filter and filt != a.filter:
                continue
            if len(filts) == 1 and not a.filter:
                continue          # one-filter sample (bells): the sample zip already is it
            write_zip(os.path.join(OUT_DIR, f'{sample}_{filt}.zip'),
                      [r for r in srows if r['filter'] == filt], files_by_dir, header, a.dry_run)


if __name__ == '__main__':
    main()
