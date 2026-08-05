#!/usr/bin/env python
"""
Tile every lens's cutout into per-filter-group QC mosaics.

Reads the postage stamps already written by make_cutouts.py
(data/cutouts/<sample>/<lens>/<filt>/cutout[_cr]_{sci,noise}.fits) and lays every
available lens out on a 5-wide grid, one mosaic per (filter group, data type). Nothing
is re-drizzled or re-cut - this only reads existing cutout FITS files.

Filter groups come from mosaic_groups.py, shared with make_psf_mosaics.py, so both
scripts stay in sync. For slacs_gold/slacs_other (SLACS: WFPC2 F606W and ACS F555W
never share a lens):
    f814W          - ACS/WFC
    f606W_f555W    - WFPC2 F606W (+ split-visit f606W_v2) merged with ACS F555W
    f160W          - WFC3/IR only (NICMOS F160W products were deleted, see CLAUDE.md)
gallery (WFC3/UVIS, no cross-filter merging) gets one group per filter: f225W, f275W,
f438W, f606W, f814W. A sample not listed in mosaic_groups.py falls back to one group
per filter subdirectory found on disk.

Each mosaic panel shows the full cutout as cut by make_cutouts.py (20" square by
default). Colour/stretch is inferno + an asinh stretch (astropy.visualization), the
standard astronomy image convention: it handles the negative background-noise pixels
smoothly (no NaN-masking artifacts) while still showing faint outskirts and bright
cores together.

--size mosaics a size-variant cutout tree (data/cutouts_<size>arcsec/) into its own
data/mosaics_<size>arcsec/ instead, mirroring make_cutouts.py so the two sizes' QC PNGs
never overwrite each other -- see cutout_paths.py.

Usage:
    uv run python scripts/make_mosaics.py --sample slacs_gold
    uv run python scripts/make_mosaics.py --sample slacs_gold --size 12
"""

import argparse
import glob
import math
import os

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from astropy.io import fits
from astropy.visualization import AsinhStretch, ImageNormalize, PercentileInterval

ws_path = '/Users/samlange/Code/HST_lenses_data_reduction'

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mast_target_names
import mosaic_groups
import cutout_paths

NCOLS = 5
CMAP = 'inferno'
NOISE_SENTINEL = 1.0e8
TITLE_FONTSIZE = 14
CBAR_FONTSIZE = 20
# Vertical room per colourbar, in INCHES, when the grid has no empty cell to put one
# in. Matches one panel cell (2.2 in), the minimum room the empty-cell path gets, so
# both routes render the bar + rotated ticks + axis label at the same proportions.
CBAR_BAND_INCHES = 2.2
# Horizontal room one tick label wants, in inches, at CBAR_FONTSIZE - sets how many
# ticks a bar of a given width carries comfortably (see spaced_ticks).
CBAR_INCHES_PER_TICK = 1.0
# Every colourbar carries at least this many ticks: two (a bare end-to-end range) reads
# as a legend rather than a scale. A bar too narrow to space that many at full size gets
# proportionally smaller tick labels instead - see spaced_ticks/style_colorbar.
MIN_CBAR_TICKS = 3
CBAR_TICK_MIN_FONTSIZE = 8


def find_cutout_pair(filt_dir):
    """Return (sci_path, noise_path) for a cutout dir, preferring the CR pass."""
    for prefix in ('cutout_cr', 'cutout'):
        sci = os.path.join(filt_dir, f'{prefix}_sci.fits')
        noise = os.path.join(filt_dir, f'{prefix}_noise.fits')
        if os.path.exists(sci) and os.path.exists(noise):
            return sci, noise
    return None, None


def short_filt(filt):
    """f606W[_v2] -> f6[_v2], f555W -> f5, for compact panel titles."""
    return filt.replace('f606W', 'f6').replace('f555W', 'f5')


def build_group(cutouts_dir, precedence):
    """Entries for one mosaic group. `precedence` is a list of filter-dir names tried
    in order per lens, so every lens contributes exactly one panel - needed for groups
    that merge filters from mutually-exclusive instruments (SLACS' WFPC2 F606W / ACS
    F555W, and the split-visit f606W/f606W_v2 keys). `group` (the per-panel colourbar-
    split key) is only set for multi-filter groups: the two sit on very different
    scales in every panel type (signal/noise: ~100x native flux-scale gap; SNR: WFPC2
    F606W is itself far noisier than ACS F555W), so those mosaics split the colourbar
    per instrument rather than share one."""
    entries = []
    for lens_dir in sorted(glob.glob(os.path.join(cutouts_dir, '*'))):
        lens = os.path.basename(lens_dir)
        for filt in precedence:
            sci, noise = find_cutout_pair(os.path.join(lens_dir, filt))
            if sci is None:
                continue
            if len(precedence) > 1:
                short = short_filt(filt)
                label, group = f'{lens} [{short}]', short.split('_')[0]
            else:
                label, group = lens, None
            entries.append((lens, label, sci, noise, group))
            break
    return entries


def load_signal_noise_snr(sci_path, noise_path):
    """Load a cutout pair and derive signal-to-noise.

    Zero-weight pixels carry the 1e8 sentinel from noise_map_via_weight_map_from
    (make_cutouts.py); they're masked to NaN here, same convention as
    make_cutouts.plot_cutouts. Wherever noise is NaN, snr is forced to 0 rather than
    NaN, so a masked pixel reads as "no signal" instead of "no data" on the SNR panel.
    """
    sci = fits.getdata(sci_path).astype(np.float64)
    noise = fits.getdata(noise_path).astype(np.float64)
    noise = np.where(noise >= NOISE_SENTINEL, np.nan, noise)

    with np.errstate(divide='ignore', invalid='ignore'):
        snr = sci / noise
    snr = np.where(np.isnan(noise), 0.0, snr)

    return sci, noise, snr


def pooled_asinh_norm(arrays, pct=99.0):
    """Shared asinh-stretch norm from the pooled finite pixels of every panel.

    Unlike a log norm, asinh is well-defined through zero and for negative values, so
    the background-noise pixels (real, negative-flux sky fluctuations) display directly
    instead of needing to be masked to NaN first.
    """
    pooled = np.concatenate([a[np.isfinite(a)].ravel() for a in arrays])
    vmin, vmax = PercentileInterval(pct).get_limits(pooled)
    return ImageNormalize(vmin=vmin, vmax=vmax, stretch=AsinhStretch(0.1))


INSTRUMENT_NAMES = {'f6': 'WFPC2 F606W', 'f5': 'ACS F555W'}


def tick_candidates(norm, scale, nbins):
    """Round tick values inside the norm's range, at most `nbins` intervals across it.

    `steps` is restricted to 1/2/5/10 and anything not a whole multiple of 0.1 scaled
    units is dropped, because style_colorbar's formatter prints one decimal of v/scale
    - a 0.25-scaled-unit tick would be labelled '0.2' and simply read as the wrong
    number. That matters here only because the denser candidate sets below exist to be
    thinned; the survivors still have to be honest values."""
    ticks = matplotlib.ticker.MaxNLocator(nbins=nbins, steps=[1, 2, 5, 10]).tick_values(
        norm.vmin, norm.vmax)
    return [t for t in ticks
            if norm.vmin <= t <= norm.vmax
            and abs(t / scale - round(t / scale, 1)) <= 1e-6]


def thin_to_gap(ticks, norm, min_gap):
    """Greedily drop ticks that land within `min_gap` (a fraction of the bar's length)
    of the last one kept."""
    kept, last = [], -np.inf
    for t in ticks:
        pos = float(norm(t))
        if pos - last >= min_gap:
            kept.append(t)
            last = pos
    return kept


def spaced_ticks(norm, scale, bar_w_in, min_ticks=MIN_CBAR_TICKS):
    """Ticks for a bar `bar_w_in` inches wide, at least `min_ticks` of them.

    Returns (ticks, squeeze), where squeeze <= 1 is how much tighter than comfortable
    the spacing had to get to reach min_ticks - the caller scales the tick font by it
    so the labels still don't touch.

    MaxNLocator picks values evenly in DATA space, but the bar is drawn in DISPLAY
    space, so under a nonlinear stretch (the log norm make_psf_mosaics.py passes) its
    values bunch up at the compressed end. Candidates are therefore mapped through the
    norm and thinned on their on-bar separation. Reaching min_ticks is tried first by
    offering DENSER candidates at the full comfortable spacing (a log bar has room for
    more ticks, just not at the round values a linear locator proposes), and only then
    by tightening the spacing - so a bar gives up label size only when its width really
    can't fit min_ticks any other way. Within one spacing the densest candidate set
    wins: the tightening loop is there to reach a floor on tick count, not impose one."""
    target = CBAR_INCHES_PER_TICK / max(bar_w_in, 1e-6)
    base_nbins = int(np.clip(bar_w_in / CBAR_INCHES_PER_TICK, 2, 4))

    best, best_squeeze = [], 1.0
    for squeeze in (1.0, 1 / 1.5, 0.5, 1 / 3, 0.25):
        densest = max(
            (thin_to_gap(tick_candidates(norm, scale, nbins), norm, target * squeeze)
             for nbins in (base_nbins, 6, 10, 16, 25)), key=len)
        if len(densest) >= min_ticks:
            return densest, squeeze
        if len(densest) > len(best):
            best, best_squeeze = densest, squeeze
    return best, best_squeeze


def fit_label_width(fig, text, max_w_in, min_fontsize=7):
    """Shrink `text` until it is at most `max_w_in` inches wide.

    The axis label is centred on its (possibly narrow) colourbar, so a long one -
    e.g. 'noise (counts/s)  (WFPC2 F606W)  (x10^-2)' over a bar squeezed into a
    single grid cell - overruns the edge of the figure and gets clipped. A double
    space in the label is a deliberate break point: wrapping there is tried first,
    so the text shrinks only as far as it has to."""
    renderer = fig.canvas.get_renderer()

    def width_in():
        return text.get_window_extent(renderer=renderer).width / fig.dpi

    if width_in() > max_w_in and '  ' in text.get_text():
        text.set_text(text.get_text().replace('  ', '\n', 1))

    for _ in range(12):
        w = width_in()
        size = text.get_fontsize()
        if w <= max_w_in or size <= min_fontsize:
            break
        text.set_fontsize(max(min_fontsize, size * max_w_in / w * 0.98))


def style_colorbar(cbar, norm, base_label, fontsize, bar_w_in, max_label_w_in,
                   rotate_ticks=False):
    """Tick labels at <=1 decimal place, with a shared '(x10^n)' multiplier in the
    label when the values need it - avoids the long decimal tick labels (e.g.
    0.0075/0.0100/0.0125) that overlap each other at this figure size.

    rotate_ticks angles the tick labels 60 degrees - used only for single-colourbar
    plots (one instrument per mosaic), where there's a full band of vertical room
    for the slanted labels; the stacked per-instrument bars don't have that room.

    bar_w_in / max_label_w_in are the bar's width and the room its centred label has
    before it runs off the figure, both in inches - everything that can overflow is
    sized against real space rather than a figure fraction, which is what let ticks
    and labels get clipped on the narrower mosaics."""
    magnitude = max(abs(norm.vmin), abs(norm.vmax))
    exponent = 0 if not np.isfinite(magnitude) or magnitude == 0 \
        else int(np.floor(np.log10(magnitude)))

    if -1 <= exponent <= 1:
        exponent = 0
    scale = 10.0 ** exponent

    ticks, squeeze = spaced_ticks(norm, scale, bar_w_in)
    cbar.set_ticks(ticks)
    cbar.ax.xaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda v, pos, s=scale: f'{v / s:.1f}'))

    label = base_label if exponent == 0 else f'{base_label}  (×10$^{{{exponent}}}$)'
    cbar.ax.tick_params(
        labelsize=max(CBAR_TICK_MIN_FONTSIZE, fontsize * squeeze))
    if rotate_ticks:
        plt.setp(cbar.ax.get_xticklabels(), rotation=60, ha='right', rotation_mode='anchor')
    cbar.set_label(label, fontsize=fontsize)
    fit_label_width(cbar.ax.figure, cbar.ax.xaxis.label, max_label_w_in)


def plot_mosaic(entries, arrays, title_key, out_path, panel_label, split_by=None,
                norm_fn=pooled_asinh_norm):
    """Render one mosaic. If split_by (a per-panel group key, e.g. instrument) is
    given, each group gets its own norm and its own horizontal colourbar, stacked in
    the empty grid cells (or in a band added below the grid when it divides evenly) -
    needed when the groups sit on very different native flux scales and a single
    shared norm would saturate one of them.

    norm_fn(arrays) -> ImageNormalize builds each group's norm; defaults to the
    pooled asinh stretch. make_psf_mosaics.py passes a log stretch instead, to match
    the 'log' wing panel in each lens's psf.png."""
    n = len(arrays)
    nrows = math.ceil(n / NCOLS)
    ncells = nrows * NCOLS

    if split_by is None:
        norm_of = {None: norm_fn(arrays)}
        panel_group = [None] * n
    else:
        groups = list(dict.fromkeys(split_by))
        norm_of = {g: norm_fn([a for a, gg in zip(arrays, split_by) if gg == g])
                   for g in groups}
        panel_group = split_by

    fig, axes = plt.subplots(nrows, NCOLS, figsize=(2.2 * NCOLS, 2.2 * nrows),
                              squeeze=False)

    im_of = {}
    for i, ax in enumerate(axes.flat):
        if i >= n:
            continue

        data = arrays[i]
        group = panel_group[i]
        im = ax.imshow(data, norm=norm_of[group], origin='upper', cmap=CMAP, aspect='equal')
        im_of[group] = im

        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_edgecolor('black')
            spine.set_linewidth(1.2)

        label = entries[i][title_key]
        ax.text(0.5, 0.97, label, transform=ax.transAxes, ha='center', va='top',
                color='white', fontsize=TITLE_FONTSIZE,
                path_effects=[pe.withStroke(linewidth=3, foreground='black')])

    fig.subplots_adjust(left=0.005, right=0.995, top=0.995, bottom=0.005,
                        wspace=0.02, hspace=0.02)

    # Put the colourbar(s) in the empty panel cells (grid cells with no lens) rather
    # than stealing a slice of figure margin. All colourbars are horizontal and span
    # the FULL combined width of every empty cell (not just one cell each). If the
    # grid divides evenly and there is no empty cell at all, the figure is instead
    # grown by a fixed band of inches below the grid (see below).
    empty_idx = list(range(n, ncells))
    for idx in empty_idx:
        axes.flat[idx].axis('off')

    groups_here = [None] if split_by is None else list(norm_of.keys())
    n_groups = len(groups_here)

    if empty_idx:
        boxes = [axes.flat[idx].get_position() for idx in empty_idx]
        x0 = min(b.x0 for b in boxes)
        x1 = max(b.x1 for b in boxes)
        y0 = min(b.y0 for b in boxes)
        y1 = max(b.y1 for b in boxes)
    else:
        # The grid divides evenly, so there is no spare cell to put the bar(s) in:
        # grow the figure downward by a fixed band of real INCHES per colourbar and
        # place them there. This must be in inches, not a figure fraction - a
        # fraction-sized band (the old bottom=0.12 / add_axes([.., 0.03, .., 0.03]))
        # shrinks with the mosaic, and on the short grids that actually hit this
        # branch (gallery's 1- and 3-row groups) it clipped the rotated tick labels,
        # the axis label, and part of the bar itself off the bottom of the figure.
        fig_w, fig_h = fig.get_size_inches()
        band_in = CBAR_BAND_INCHES * n_groups
        fig.set_size_inches(fig_w, fig_h + band_in)
        # subplots_adjust is in figure fractions, so redo it against the new height:
        # keep the panel grid at its original absolute size at the top, leaving the
        # new band clear underneath.
        shrink = fig_h / (fig_h + band_in)
        band_frac = band_in / (fig_h + band_in)
        fig.subplots_adjust(left=0.005, right=0.995,
                            top=1 - 0.005 * shrink, bottom=band_frac + 0.005 * shrink,
                            wspace=0.02, hspace=0.02)
        x0, x1 = 0.20, 0.80
        y0, y1 = 0.0, band_frac

    width, height = x1 - x0, y1 - y0
    pad = 0.06 * width
    band_h = height / n_groups
    fontsize = CBAR_FONTSIZE if n_groups == 1 else CBAR_FONTSIZE * 0.6
    # Tick count and label size follow the bar's real width: a bar confined to a
    # single empty grid cell (~1.9 in) fits ~2 intervals, and forcing the default 4
    # there overlapped the rotated labels into each other. The label is centred on
    # the bar, so the room it has before running off the figure is set by whichever
    # side is tighter.
    fig_w_in = fig.get_size_inches()[0]
    bar_w_in = (width - 2 * pad) * fig_w_in
    centre = x0 + width / 2
    max_label_w_in = 2 * min(centre, 1 - centre) * fig_w_in * 0.96
    for gi, group in enumerate(groups_here):
        # Bar sits near the top of its band; the rest of the band (~70%) is left
        # free for its tick labels + axis label so stacked bands don't collide.
        band_top = y1 - gi * band_h
        bar_h = 0.20 * band_h
        bar_top = band_top - 0.06 * band_h
        bar_bottom = bar_top - bar_h
        cax = fig.add_axes([x0 + pad, bar_bottom, width - 2 * pad, bar_h])
        base_label = panel_label if group is None else \
            f'{panel_label}  ({INSTRUMENT_NAMES.get(group, group)})'
        cbar = fig.colorbar(im_of[group], cax=cax, orientation='horizontal')
        style_colorbar(cbar, norm_of[group], base_label, fontsize,
                       bar_w_in, max_label_w_in, rotate_ticks=(n_groups == 1))

    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  wrote {out_path}  ({n} lenses, {nrows}x{NCOLS} grid)")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--sample', default=mast_target_names.DEFAULT_SAMPLE,
                   help='sample subdirectory under data/cutouts/ to mosaic. Defined in '
                        f'info/lens_samples.json (default {mast_target_names.DEFAULT_SAMPLE})')
    p.add_argument('--size', type=float, default=cutout_paths.DEFAULT_SIZE,
                   help='cutout size in arcsec to mosaic, matching make_cutouts.py --size '
                        f'(default {cutout_paths.DEFAULT_SIZE:g}). A non-default size reads '
                        'data/cutouts_<size>arcsec/ and writes data/mosaics_<size>arcsec/')
    a = p.parse_args()

    cutouts_dir = os.path.join(cutout_paths.cutouts_root(ws_path, a.size), a.sample)
    out_dir = os.path.join(cutout_paths.mosaics_root(ws_path, a.size), a.sample)
    os.makedirs(out_dir, exist_ok=True)
    print(f"cutouts: {cutouts_dir}\nmosaics: {out_dir}")

    groups = mosaic_groups.groups_for_sample(a.sample, cutouts_dir)
    for group_name, precedence in groups.items():
        entries_raw = build_group(cutouts_dir, precedence)
        if not entries_raw:
            print(f"{group_name}: no cutouts found under {cutouts_dir}, skipping")
            continue

        entries = [{'lens': lens, 'label': label} for lens, label, _, _, _ in entries_raw]
        groups = [group for _, _, _, _, group in entries_raw]
        # Only split the colourbar when the group actually mixes >1 instrument -
        # f814W/f160W are single-instrument, so groups is all-None there.
        split_by = groups if any(g is not None for g in groups) else None
        print(f"{group_name}: {len(entries)} lenses")

        signals, noises, snrs = [], [], []
        for lens, label, sci_path, noise_path, group in entries_raw:
            sig, noi, snr = load_signal_noise_snr(sci_path, noise_path)
            signals.append(sig)
            noises.append(noi)
            snrs.append(snr)

        plot_mosaic(entries, signals, 'label',
                    os.path.join(out_dir, f'{group_name}_signal.png'), 'signal (counts/s)',
                    split_by=split_by)
        plot_mosaic(entries, noises, 'label',
                    os.path.join(out_dir, f'{group_name}_noise.png'), 'noise (counts/s)',
                    split_by=split_by)
        # SNR is a ratio, so the flux-scale gap between instruments cancels out - but
        # WFPC2 F606W's true SNR is itself far lower than ACS F555W's (much noisier
        # detector/shorter exposures), so it still needs its own colourbar or its
        # panels wash out under ACS's much higher peak SNR.
        plot_mosaic(entries, snrs, 'label',
                    os.path.join(out_dir, f'{group_name}_snr.png'), 'signal / noise',
                    split_by=split_by)


if __name__ == '__main__':
    main()
