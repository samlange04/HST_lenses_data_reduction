#!/usr/bin/env python
"""
Choose each lens's best (primary) band BY EYE: a GUI over the lens's band overview.

Merged 2026-11-28 across the HST and JWST repos: the file is identical in both (the repo
profile is the import seam below).

    uv run python scripts/choose_best_band.py --sample swells                      # every lens
    uv run python scripts/choose_best_band.py --sample swells --unreviewed
    uv run python scripts/choose_best_band.py --lens J0841+3824

Ported from JWST_lenses_data_reduction 2026-11-10. The primary band is the one where the lensed
images are best defined: the mask GUIs draw it first, positions are clicked on it, and it is the
band to model first. choose_primary_band.py measures a suggestion (arc S/N discounted for coarse
resolution), but the measurement is only a guide -- it is biased toward the red bands for bright,
compact lenses, where lens light leaks into the ring annulus. This makes the choice the
reviewer's.

For each lens it opens the band overview (make_band_overview.build_figure: every band as S/N on
ONE stretch shared by all bands, so the bands can be compared by eye; full stamp, zoom,
isophote-subtracted zoom), with each band's measured arc S/N and score under its name and the
current choice outlined.
  * click a band's column, or press 1-N (bands in wavelength order), to select it;
  * Enter accepts the selection: if it IS the measured suggestion, the record stays
    `arc_snr`-sourced and is marked `reviewed` (choose_primary_band.set_reviewed); if it is a
    DIFFERENT band, it is written to info/samples/lens_samples.json as the lens's
    `best_band_visual` (source 'visual', reviewed by construction), a visual inspection and
    not necessarily the highest arc S/N, with a note naming this script and the best arc-S/N
    band. Either way the RECORDS are refreshed at once (choose_primary_band.record), but the
    band overview PNGs are re-rendered in ONE pass at the end of the run: a render re-fits
    every band's isophote model, the slow step, and batching it keeps the GUI moving between
    lenses;
  * n (or closing the window) skips the lens, leaving it as it was; q quits. (Skip is NOT 's':
    matplotlib binds 's' to 'save figure' in every window.)
--unreviewed shows only lenses without a reviewed choice (no best_band_visual and no reviewed
arc_snr record). Single-band lenses are never shown: there is nothing to choose.
The GUI shows the RECORDED arc-S/N measurement (it reads lens_primary_band.json, like the
overview's column titles) rather than re-measuring -- re-run choose_primary_band.py first if
the stamps or masks changed.
"""
import argparse
import datetime
import glob
import os
import sys
import warnings

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import choose_primary_band
import make_band_overview
import stamp_paths
try:
    import lens_samples as mast_target_names   # the JWST repo's profile module
except ImportError:
    import mast_target_names                   # the HST repo's

ws_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PICK = '#d62728'


def choose(sample, lens):
    """Show the lens and return ('pick', band), ('skip', None) or ('quit', None)."""
    import matplotlib.pyplot as plt
    # The recorded measurement, not a fresh measure(): the overview's column titles come from
    # the same record, and measure() re-reads every band's FITS for no new information here.
    pb_rec = choose_primary_band.info_json.load(
        choose_primary_band.info_json.path(ws_path, 'lens_primary_band.json')).get(sample, {}).get(lens, {})
    suggested = pb_rec.get('best_band_arc_snr')
    current, source = mast_target_names.primary_band(lens, sample)
    fig = plt.figure()
    fig, info = make_band_overview.build_figure(sample, lens, fig=fig)
    if fig is None:
        plt.close('all')
        return 'skip', None
    filts, axes = info['filts'], info['axes']
    state = dict(sel=current if current in filts else (suggested or filts[-1]),
                 action='skip')

    for j, f in enumerate(filts):
        # build_figure's title already carries the measured arc S/N and score; the GUI adds
        # only the key number and the [best arc S/N] tag.
        tag = '  [best arc S/N]' if f == suggested else ''
        axes[0][j].set_title(f'{j + 1}: {axes[0][j].get_title()}{tag}', fontsize=9)

    def paint():
        for j, f in enumerate(filts):
            on = f == state['sel']
            for i in range(3):
                for sp in axes[i][j].spines.values():
                    sp.set_edgecolor(PICK if on else 'black')
                    sp.set_linewidth(3.5 if on else 0.8)
        rev = '' if source in (None, 'visual') else (', reviewed' if pb_rec.get('reviewed') else ', UNREVIEWED')
        fig.texts[-1].set_text(f'selected: {state["sel"]}    (current: {current or "none"}'
                               f'{f", {source}{rev}" if current else ""}; best arc S/N: {suggested or "none"})\n'
                               f'click a column or press 1-{len(filts)}; Enter accepts, n skips, q quits')
        fig.canvas.draw_idle()

    fig.tight_layout(rect=(0, 0.055, 1, 0.90))         # room for the two instruction lines
    fig.text(0.5, 0.008, '', ha='center', va='bottom', fontsize=10, color=PICK, weight='bold')

    def onclick(event):
        for j in range(len(filts)):
            if any(event.inaxes is axes[i][j] for i in range(3)):
                state['sel'] = filts[j]
                paint()

    def onkey(event):
        k = event.key
        if k and k.isdigit() and 1 <= int(k) <= len(filts):
            state['sel'] = filts[int(k) - 1]
            paint()
        elif k == 'enter':
            state['action'] = 'pick'
            plt.close(fig)
        elif k == 'n':
            state['action'] = 'skip'
            plt.close(fig)
        elif k == 'q':
            state['action'] = 'quit'
            plt.close(fig)

    fig.canvas.mpl_connect('button_press_event', onclick)
    fig.canvas.mpl_connect('key_press_event', onkey)
    paint()
    plt.show()
    return (state['action'], state['sel'] if state['action'] == 'pick' else None)


def n_bands(sample, lens):
    """How many bands have a stamp for this lens (single-band lenses have nothing to choose)."""
    lens_dir = stamp_paths.lens_dir(ws_path, 'stamps', sample, lens)
    prefix = stamp_paths.stamp_prefix(stamp_paths.modelling_size(lens, sample))
    n = 0
    for d in glob.glob(os.path.join(lens_dir, '*')):
        if (os.path.isdir(d) and os.path.basename(d).startswith(('f', 'F'))
                and stamp_paths.find_stamp(d, prefix, 'sci')):
            n += 1
    return n


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--sample', type=mast_target_names.sample_arg, default=None,
                   help='every lens of this sample (not needed with --lens: its sample is found)')
    p.add_argument('--lens', default=None)
    p.add_argument('--unreviewed', action='store_true', help='only lenses without a best_band_visual')
    a = p.parse_args()
    warnings.simplefilter('ignore')
    mast_target_names.resolve_sample(p, a, required=True)
    sample = a.sample
    lenses = [a.lens] if a.lens else mast_target_names.lenses(sample)
    pb = choose_primary_band.info_json.load(
        choose_primary_band.info_json.path(ws_path, 'lens_primary_band.json')).get(sample, {})
    n_shown = 0
    changed = []
    for lens in lenses:
        if n_bands(sample, lens) < 2:
            print(f'{lens}: single band -- nothing to choose, skipped')
            continue
        current, source = mast_target_names.primary_band(lens, sample)
        reviewed = source == 'visual' or bool(pb.get(lens, {}).get('reviewed'))
        if a.unreviewed and reviewed:
            continue
        n_shown += 1
        print(f'{lens}: current {current} ({source}) -- choose in the window')
        action, band = choose(sample, lens)
        if action == 'quit':
            print('quit')
            break
        if action == 'skip':
            print('  skipped, unchanged')
            continue
        # rerender=False: a render re-fits every band's isophote model (the slow step), so the
        # overviews are re-rendered in one pass at the end of the run instead of per lens.
        suggested = pb.get(lens, {}).get('best_band_arc_snr')
        if band is not None and band == suggested:
            # Accepted the measured suggestion: the record stays arc_snr-sourced, marked reviewed.
            choose_primary_band.record(sample, lens, rerender=False)
            choose_primary_band.set_reviewed(
                sample, lens, True,
                note='Confirmed by eye with choose_best_band.py (the measured band accepted).')
            print(f'  confirmed the measured band {band} (arc_snr, reviewed)')
        else:
            note = (f'Chosen by eye with choose_best_band.py ({datetime.date.today().isoformat()}); '
                    f'the best arc-S/N band (choose_primary_band.py) was '
                    f'{suggested or "none (arcs too faint to measure)"}.')
            mast_target_names.set_best_band_visual(lens, band, note, sample)
            choose_primary_band.record(sample, lens, rerender=False)
            print(f'  best_band_visual {band}')
        changed.append(lens)
    if changed:
        print(f're-rendering the band overviews of the {len(changed)} changed lens(es)...')
        for lens in changed:
            make_band_overview.rerender(sample, lens)


if __name__ == '__main__':
    sys.exit(main())
