#!/usr/bin/env python
"""
Tile every lens's cutout into per-filter-group QC mosaics.

Reads the postage stamps already written by make_cutouts.py
(data/cutouts/<sample>/<lens>/<filt>/cutout[_cr]_{sci,noise}.fits) and lays every
available lens out on a 5-wide grid, one mosaic per (filter group, data type). Nothing
is re-drizzled or re-cut - this only reads existing cutout FITS files.

Three filter groups (SLACS-specific: WFPC2 F606W and ACS F555W never share a lens, so
they're combined into one grid with the filter named per-panel):
    f814W          - ACS/WFC, 38 lenses
    f606W_f555W    - WFPC2 F606W (+ split-visit f606W_v1/v2) merged with ACS F555W
    f160W          - WFC3/IR only (NICMOS F160W products were deleted, see CLAUDE.md)

Each mosaic panel shows the full cutout as cut by make_cutouts.py (10" square by
default). Colour/stretch is inferno + an asinh stretch (astropy.visualization), the
standard astronomy image convention: it handles the negative background-noise pixels
smoothly (no NaN-masking artifacts) while still showing faint outskirts and bright
cores together.

Usage:
    conda run -n stenv python scripts/make_mosaics.py --sample slacs
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

ws_path = '/Users/samlange/Code/data_reduction'

NCOLS = 5
CMAP = 'inferno'
NOISE_SENTINEL = 1.0e8
TITLE_FONTSIZE = 14
CBAR_FONTSIZE = 20


def find_cutout_pair(filt_dir):
    """Return (sci_path, noise_path) for a cutout dir, preferring the CR pass."""
    for prefix in ('cutout_cr', 'cutout'):
        sci = os.path.join(filt_dir, f'{prefix}_sci.fits')
        noise = os.path.join(filt_dir, f'{prefix}_noise.fits')
        if os.path.exists(sci) and os.path.exists(noise):
            return sci, noise
    return None, None


def build_group_f814w(cutouts_dir):
    entries = []
    for lens_dir in sorted(glob.glob(os.path.join(cutouts_dir, '*'))):
        lens = os.path.basename(lens_dir)
        filt_dir = os.path.join(lens_dir, 'f814W')
        sci, noise = find_cutout_pair(filt_dir)
        if sci is None:
            continue
        entries.append((lens, lens, sci, noise, None))
    return entries


def short_filt(filt):
    """f606W[_v1/_v2] -> f6[_v1/_v2], f555W -> f5, for compact panel titles."""
    return filt.replace('f606W', 'f6').replace('f555W', 'f5')


def build_group_f606w_f555w(cutouts_dir):
    # Precedence per lens: combined f606W, else the deeper split-visit f606W_v2,
    # else f606W_v1, else f555W - so every lens appears exactly once. `group` is the
    # instrument (f6* is WFPC2, f5 is ACS) - the two sit on very different scales in
    # every panel type (signal/noise: ~100x native flux-scale gap; SNR: WFPC2 F606W
    # is itself far noisier than ACS F555W), so all three mosaics split the colourbar
    # per instrument rather than share one.
    precedence = ['f606W', 'f606W_v2', 'f606W_v1', 'f555W']
    entries = []
    for lens_dir in sorted(glob.glob(os.path.join(cutouts_dir, '*'))):
        lens = os.path.basename(lens_dir)
        for filt in precedence:
            sci, noise = find_cutout_pair(os.path.join(lens_dir, filt))
            if sci is not None:
                short = short_filt(filt)
                group = short.split('_')[0]
                entries.append((lens, f'{lens} [{short}]', sci, noise, group))
                break
    return entries


def build_group_f160w(cutouts_dir):
    entries = []
    for lens_dir in sorted(glob.glob(os.path.join(cutouts_dir, '*'))):
        lens = os.path.basename(lens_dir)
        filt_dir = os.path.join(lens_dir, 'f160W')
        sci, noise = find_cutout_pair(filt_dir)
        if sci is None:
            continue
        entries.append((lens, lens, sci, noise, None))
    return entries


GROUPS = {
    'f814W':       build_group_f814w,
    'f606W_f555W': build_group_f606w_f555w,
    'f160W':       build_group_f160w,
}


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


def style_colorbar(cbar, norm, base_label, fontsize, rotate_ticks=False):
    """Tick labels at <=1 decimal place, with a shared '(x10^n)' multiplier in the
    label when the values need it - avoids the long decimal tick labels (e.g.
    0.0075/0.0100/0.0125) that overlap each other at this figure size.

    rotate_ticks angles the tick labels 60 degrees - used only for single-colourbar
    plots (one instrument per mosaic), where there's a full band of vertical room
    for the slanted labels; the stacked per-instrument bars don't have that room."""
    magnitude = max(abs(norm.vmin), abs(norm.vmax))
    exponent = 0 if not np.isfinite(magnitude) or magnitude == 0 \
        else int(np.floor(np.log10(magnitude)))

    if -1 <= exponent <= 1:
        exponent = 0
    scale = 10.0 ** exponent

    cbar.ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=4))
    cbar.ax.xaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda v, pos, s=scale: f'{v / s:.1f}'))

    label = base_label if exponent == 0 else f'{base_label}  (×10$^{{{exponent}}}$)'
    cbar.ax.tick_params(labelsize=fontsize)
    if rotate_ticks:
        plt.setp(cbar.ax.get_xticklabels(), rotation=60, ha='right', rotation_mode='anchor')
    cbar.set_label(label, fontsize=fontsize)


def plot_mosaic(entries, arrays, title_key, out_path, panel_label, split_by=None):
    """Render one mosaic. If split_by (a per-panel group key, e.g. instrument) is
    given, each group gets its own norm and its own horizontal colourbar, spread
    across the empty grid cells - needed when the groups sit on very different
    native flux scales and a single shared norm would saturate one of them."""
    n = len(arrays)
    nrows = math.ceil(n / NCOLS)
    ncells = nrows * NCOLS

    if split_by is None:
        norm_of = {None: pooled_asinh_norm(arrays)}
        panel_group = [None] * n
    else:
        groups = list(dict.fromkeys(split_by))
        norm_of = {g: pooled_asinh_norm([a for a, gg in zip(arrays, split_by) if gg == g])
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
    # the FULL combined width of every empty cell (not just one cell each) - falls
    # back to an external axis only if the grid divides evenly and no empty cell
    # exists at all.
    empty_idx = list(range(n, ncells))
    for idx in empty_idx:
        axes.flat[idx].axis('off')

    if not empty_idx:
        fig.subplots_adjust(bottom=0.12)
        cax = fig.add_axes([0.25, 0.03, 0.5, 0.03])
        group0 = panel_group[0]
        cbar = fig.colorbar(im_of[group0], cax=cax, orientation='horizontal', label=panel_label)
        style_colorbar(cbar, norm_of[group0], panel_label, CBAR_FONTSIZE,
                       rotate_ticks=(split_by is None))
    else:
        boxes = [axes.flat[idx].get_position() for idx in empty_idx]
        x0 = min(b.x0 for b in boxes)
        x1 = max(b.x1 for b in boxes)
        y0 = min(b.y0 for b in boxes)
        y1 = max(b.y1 for b in boxes)
        width, height = x1 - x0, y1 - y0
        pad = 0.06 * width

        groups_here = [None] if split_by is None else list(norm_of.keys())
        n_groups = len(groups_here)
        band_h = height / n_groups
        fontsize = CBAR_FONTSIZE if n_groups == 1 else CBAR_FONTSIZE * 0.6
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
                           rotate_ticks=(n_groups == 1))

    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"  wrote {out_path}  ({n} lenses, {nrows}x{NCOLS} grid)")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--sample', default='slacs',
                   help='sample subdirectory under data/cutouts/ to mosaic (default slacs; '
                        'other samples, e.g. bells, will work once their cutouts exist)')
    a = p.parse_args()

    cutouts_dir = os.path.join(ws_path, 'data', 'cutouts', a.sample)
    out_dir = os.path.join(ws_path, 'data', 'mosaics', a.sample)
    os.makedirs(out_dir, exist_ok=True)

    for group_name, builder in GROUPS.items():
        entries_raw = builder(cutouts_dir)
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
