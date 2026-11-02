#!/usr/bin/env python
"""
Check that every lens's data sit where its lens_samples.json keys say they should.

Each lens lives in data/<layer>/<sample>/<band_mode>/<category>/<lens>/ in every data layer
(calibrated, drizzle_files, drizzled, psf, stamps, and the variant trees such as
drizzled_bcfill). `band_mode` and `category` are stored in info/samples/lens_samples.json rather
than derived, so reducing a new band or re-classifying a lens never moves anything silently;
this script is what catches the drift. Reported, per lens:

  * a lens dir anywhere other than its expected path (a stale copy after a reclassification);
  * band_mode not matching the stamps on disk (distinct filters; `f606W_visit2` is the same
    band as `f606W`), e.g. after a new band is reduced;
  * a lens dir under no lens of the roster.

Read-only; exit status 1 if anything is wrong, so it can gate a commit.

Usage:
    uv run python scripts/check_layout.py
"""
import glob
import os
import sys

ws_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ws_path, 'scripts'))
import mast_target_names  # noqa: E402

LAYER_PREFIXES = ('calibrated', 'drizzle_files', 'drizzled', 'psf', 'stamps')


def main():
    data = os.path.join(ws_path, 'data')
    layers = sorted(d for d in os.listdir(data)
                    if d.startswith(LAYER_PREFIXES) and os.path.isdir(os.path.join(data, d)))
    expected = {}
    for sample in mast_target_names.samples():
        for lens in mast_target_names.lenses(sample):
            expected[lens] = mast_target_names.lens_rel(lens, sample)

    problems = []
    for layer in layers:
        root = os.path.join(data, layer)
        for d in glob.glob(os.path.join(root, '*', '*', '*', '*')):
            if not os.path.isdir(d):
                continue
            rel = os.path.relpath(d, root)
            lens = os.path.basename(d)
            if lens not in expected:
                problems.append(f'{layer}/{rel}: not a lens in the roster')
            elif rel != expected[lens]:
                problems.append(f'{layer}/{rel}: expected {layer}/{expected[lens]}')
        for d in glob.glob(os.path.join(root, '*', '*')):
            # anything left at the pre-2026-11-02 depth (<sample>/<lens>)
            if os.path.basename(d) in expected:
                problems.append(f'{layer}/{os.path.relpath(d, root)}: old layout')

    for lens, rel in sorted(expected.items()):
        bands = glob.glob(os.path.join(data, 'stamps', rel, '*'))
        filts = {os.path.basename(b).split('_')[0] for b in bands if os.path.isdir(b)}
        if not filts:
            continue
        mode = 'multi_band' if len(filts) > 1 else 'single_band'
        if mode != rel.split(os.sep)[1]:
            problems.append(f'{lens}: band_mode {rel.split(os.sep)[1]} but stamps hold '
                            f'{", ".join(sorted(filts))} -> {mode}')

    for p in problems:
        print('  ' + p)
    print(f'{len(expected)} lenses, {len(layers)} layers checked: '
          f'{len(problems) or "no"} problem(s)')
    return 1 if problems else 0


if __name__ == '__main__':
    sys.exit(main())
