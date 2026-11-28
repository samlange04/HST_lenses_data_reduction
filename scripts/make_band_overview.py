#!/usr/bin/env python
"""
One figure per lens with every band's stamp side by side:
data/stamps/<sample>/.../<lens>/stamp<S>_bands.png.

    uv run python scripts/make_band_overview.py                          # every lens of every sample
    uv run python scripts/make_band_overview.py --sample swells [--lens LENS]

Merged 2026-11-28 across the HST and JWST repos: the file is identical in both, the repo
profile is the import seam below, and the stamp layout is each repo's own `stamp_paths`.

The per-band QC (stamp<S>_dataset.png) shows one band's data, noise, PSF and mask; this shows
all of a lens's bands at once, so the lens and its arcs can be compared across wavelength.
Columns are the bands in wavelength order. Rows:
  1. the full stamp as S/N (data / noise), with the band's contaminant mask outlined (black on
     a white halo) where drawn, and the extra-galaxy centres (cyan x) where a multi_galaxy lens
     has them (make_extra_galaxies.py);
  2. a zoom on the lens (+-max(2.5 theta_E, 1.5"); at least a quarter of the stamp for the
     group-centred multi_galaxy lenses, whose groups a 2.5 theta_E zoom cuts apart), as S/N;
  3. the same zoom with the lens light removed by an ISOPHOTE model (isophote_model.py:
     elliptical isophotes free in ellipticity and PA; the positions GUI's panel), divided by
     the noise, smoothed and on a red-blue map (red = light the model lacks, masked pixels
     blank), with the Einstein radius (dashed, centred on the deflector when the profile knows
     its position). Image positions are NOT drawn: they have their own per-band QC overlays
     (make_positions.py, *_positions.png).
Every band shares ONE stretch on the S/N rows: their limits come from all bands' pixels combined
(PercentileInterval(99) + AsinhStretch, display_look.shared_norm). A per-panel stretch made a
faint band look as bright as a deep one, which hid exactly the arc-S/N and deblending
differences the best-band choice (choose_best_band.py) is made from. The residual row keeps a
per-band symmetric sigma limit (a shared one let the deepest band wash out the rest), and each
column's title carries the measured arc S/N and score (choose_primary_band.py's record, as in
the choose_best_band.py GUI). Every S/N panel has its OWN colourbar: the bar spans the one
shared stretch but is ticked at that panel's min-med-max (max = the panel's own 99th
percentile), so the deepest band's ticks span the bar and a fainter band's max tick lands
partway down it at its own value (2026-11-14: the shared stretch alone left a faint band's
panel unreadable); each residual panel has its own too.
The zoom windows are centred on the stamp centre (the shared deflector/group centre the stamps
are cut on). The isophote fit starts at the deflector when the profile knows its position
(isophote_model.deflector_pixel: the JWST measured centres, or the HST group stamps' catalogue
position), without the centroid refinement, which would creep toward arc light; otherwise it
refines from the stamp centre.
On screen North is DOWN and East is LEFT (drawn origin='upper', display_look.ORIGIN) -- the one
convention of both repos since 2026-11-27, shared with the positions/mask GUIs, the overlay
PNGs and the dataset subplots. (Until then the JWST side of this figure was drawn North-up,
which read flipped against the positions overlays mid-review.)
The primary band (primary_band(): reviewed, or choose_primary_band.py's suggestion) is starred,
and the header gives the lens and source redshifts (a photometric one labelled 'photo-z': only
a spectroscopic redshift is confirmed), theta_E and the primary band's source.

Tracked in git with the stamps (it is a stamp<S>_* file), so it costs ~0.3-0.6 MB of history
per render. make_masks.py and choose_primary_band.py re-render it after they write; a re-cut of
the stamps (make_stamps.py) does not, so run this after one.
"""
import argparse
import glob
import os
import sys
import warnings

import numpy as np
import matplotlib.cm
import matplotlib.patheffects
from matplotlib.ticker import FuncFormatter
from matplotlib.figure import Figure     # no pyplot: GUIs import this module (see build_figure)
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import display_look
import info_json
import isophote_model
import positions_io
import stamp_paths
try:
    import lens_samples as mast_target_names   # the JWST repo's profile module
except ImportError:
    import mast_target_names                   # the HST repo's

CORE_ARCSEC = 0.2      # the core left out of of the isophote panel's colour range (it saturates)
EXTRA_COLOUR = 'cyan'  # extra-galaxy centres (the JWST repo's make_extra_galaxies.py)
SMOOTH_PX = 1.5        # the isophote residual's smoothing (Gaussian sigma, px), as the positions GUI
HALO = [matplotlib.patheffects.withStroke(linewidth=3, foreground='white')]

ws_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Pivot wavelength (um) for the column order, both repos' filters. A split visit (f606W_visit2)
# sorts at its filter's wavelength via _base.
FILTER_UM = {'f435W': 0.435, 'f450W': 0.450, 'f475X': 0.475, 'f555W': 0.555, 'f600LP': 0.600,
             'f606W': 0.606, 'f814W': 0.814, 'f160W': 1.60,
             'F070W': 0.70, 'F090W': 0.90, 'F115W': 1.15, 'F150W': 1.50, 'F200W': 2.00,
             'F277W': 2.77, 'F356W': 3.56, 'F410M': 4.10, 'F444W': 4.44}


def _base(filt):
    """f606W_visit2 -> f606W: a split visit sorts at its filter's wavelength."""
    return filt.split('_')[0]


def load_band(stamp_dir, prefix):
    """The band's sci and noise (no-data pixels NaN), pixel scale, mask (or None), and
    extra-galaxy centres (or [])."""
    sci_f = stamp_paths.find_stamp(stamp_dir, prefix, 'sci')
    noise_f = stamp_paths.find_stamp(stamp_dir, prefix, 'noise')
    if sci_f is None:
        return None
    sci = fits.getdata(sci_f).astype(float)
    hdr = fits.getheader(sci_f)
    noise = fits.getdata(noise_f).astype(float) if noise_f is not None else None
    if noise is not None:
        sci = np.where(noise >= 1e7, np.nan, sci)          # the no-coverage sentinel is 1e8
    scale = float(proj_plane_pixel_scales(WCS(hdr).celestial)[0]) * 3600.0
    mask_f = os.path.join(stamp_dir, f'{prefix}_mask.fits')
    mask = fits.getdata(mask_f).astype(bool) if os.path.exists(mask_f) else None
    eg_f = os.path.join(stamp_dir, f'{prefix}_extra_galaxies_centres.json')
    extra = positions_io.flatten(positions_io.read_sets(eg_f)) if os.path.exists(eg_f) else []
    return dict(sci=sci, noise=noise, scale=scale, mask=mask, extra=extra, hdr=hdr,
                sci_path=sci_f)


def build_figure(sample, lens, fig=None):
    """Draw the overview onto `fig` (a new plain Figure by default) and return (fig, info), info
    holding the bands in column order and the axes grid, or (None, None) without stamps.

    Plain `Figure`, never pyplot: make_masks.py and choose_best_band.py import this module
    mid-GUI, and a pyplot import (or matplotlib.use) here would switch their backend, so every
    later GUI window would flash and return unanswered."""
    from autogalaxy.gui.display_util import mask_boundary
    prefix = stamp_paths.stamp_prefix(stamp_paths.modelling_size(lens, sample))
    lens_dir = stamp_paths.lens_dir(ws_path, 'stamps', sample, lens)
    bands = {}
    for d in sorted(glob.glob(os.path.join(lens_dir, '*'))):
        if not os.path.isdir(d) or not os.path.basename(d).startswith(('f', 'F')):
            continue
        b = load_band(d, prefix)
        if b is not None:
            bands[os.path.basename(d)] = b
    if not bands:
        return None, None
    filts = sorted(bands, key=lambda f: (FILTER_UM.get(_base(f), 99), f))
    te = mast_target_names.theta_e(lens, sample)
    zl, zs = mast_target_names.redshift_pair(lens, sample)
    kl, ks = mast_target_names.redshift_kinds(lens, sample)
    zl_s, zs_s = mast_target_names.redshift_text(zl, kl), mast_target_names.redshift_text(zs, ks)
    primary, psource = mast_target_names.primary_band(lens, sample)
    cat = mast_target_names.category(lens, sample)
    half_zoom = max(2.5 * (te or 0.6), 1.5)
    if cat == 'multi_galaxy':
        # group-centred stamps: a 2.5 theta_E zoom frames one member, not the group
        half_zoom = max(half_zoom, stamp_paths.modelling_size(lens, sample) / 4.0)

    ncol = len(filts)
    size = (max(3.3 * ncol, 8.5), 9.8)
    if fig is None:
        fig = Figure(figsize=size)
    else:
        fig.clf()
        fig.set_size_inches(*size)
    # panel columns alternate with slim colourbar columns, so a colourbar never shrinks a panel;
    # every row uses its band's slim column (the S/N rows' bars are per-panel -- see below)
    gs = fig.add_gridspec(3, 2 * ncol,
                          width_ratios=[w for _ in range(ncol) for w in (1.0, 0.05)])
    axes = np.empty((3, ncol), dtype=object)
    cax3 = []
    for j in range(ncol):
        for i in range(3):
            axes[i, j] = fig.add_subplot(gs[i, 2 * j])
        cax3.append(fig.add_subplot(gs[2, 2 * j + 1]))
    # The panels show S/N (data / noise), and every band is drawn on ONE shared stretch, so a
    # band where the arc is faint looks faint. A per-panel stretch made all bands look equally
    # bright, hiding exactly the arc-S/N and deblending differences the best-band choice
    # (choose_best_band.py) is made from.
    for f in filts:
        b = bands[f]
        sci, s, mask = b['sci'], b['scale'], b['mask']
        noise = b['noise'] if b['noise'] is not None else np.ones_like(sci)
        n = sci.shape[0]
        c = (n - 1) / 2.0
        # The deflector position, where the profile knows it (the JWST measured centres, or the
        # HST group stamps' catalogue position), as this band's pixel and as (West, South) arcsec
        # from the stamp centre -- the axes' frame (see the panel loop). None: the stamp centre.
        # `refine` says whether the isophote fit centroid-refines from it (a measured position
        # is exact, so refining would only creep toward arc light; a catalogue one needs it).
        got = isophote_model.deflector_pixel(b['sci_path'])
        dp = got[:2] if got is not None else None
        if dp is not None:
            b['dp_ws'] = ((dp[0] - c) * s, (c - dp[1]) * s)
            b['dp_exact'] = not got[2]
        bad = ~np.isfinite(sci) if mask is None else (mask | ~np.isfinite(sci))
        # the memoised helper: a GUI's display build and the end-of-run re-render then share
        # each band's isophote fit (the slow step) instead of fitting twice.
        sub = isophote_model.isophote_display(
            sci, s, exclude=bad, smooth_px=SMOOTH_PX, centre=dp,
            refine_centre=got[2] if got is not None else None)[0]
        with np.errstate(divide='ignore', invalid='ignore'):
            b['snr'] = np.where(np.isfinite(sci), sci / noise, np.nan)
            b['sub_snr'] = sub / noise
        zp = int(np.ceil(half_zoom / s)) + 1
        b['zoom'] = (slice(max(0, int(round(c)) - zp), min(n, int(round(c)) + zp + 1)),
                     slice(max(0, int(round(c)) - zp), min(n, int(round(c)) + zp + 1)))
        # this band's symmetric limit for the residual row, from the zoom outside the core (the
        # inner isophotes cannot follow the central cusp, and that residual would wash out the arcs)
        yy, xx = np.mgrid[b['zoom']]
        outer = np.hypot(xx - c, yy - c) * s > CORE_ARCSEC
        a = np.abs(b['sub_snr'][b['zoom']][outer])
        a = a[np.isfinite(a)]
        b['vmax3'] = float(np.percentile(a, 99.5)) if a.size else 1.0
    norm_snr = display_look.shared_norm([bands[f]['snr'] for f in filts])
    # the measured arc S/N and score per band (choose_primary_band.py's record; the
    # choose_best_band.py GUI shows the same numbers), for the column titles
    pb = info_json.load(info_json.path(ws_path, 'lens_primary_band.json')).get(sample, {}).get(lens, {})
    scores = pb.get('scores', {})
    for j, f in enumerate(filts):
        b = bands[f]
        s, mask = b['scale'], b['mask']
        n = b['sci'].shape[0]
        c = (n - 1) / 2.0
        for i, (img, zoom) in enumerate(((b['snr'], False), (b['snr'], True), (b['sub_snr'], True))):
            ax = axes[i, j]
            rs, cs = b['zoom'] if zoom else (slice(0, n), slice(0, n))
            a = img[rs, cs]
            diverging = i == 2
            # North DOWN, East LEFT (display_look.ORIGIN): pixel (col, row) sits at
            # x = (col - c) s (+col = West, to the right), y = (c - row) s (+row = North, so
            # North is at the bottom). This is autoarray's own arcsec frame, so overlay points
            # stored as (y, x) arcsec (extra-galaxy centres, positions) plot at (x, y) directly.
            ext = [(cs.start - 0.5 - c) * s, (cs.stop - 0.5 - c) * s,
                   (c - (rs.stop - 0.5)) * s, (c - (rs.start - 0.5)) * s]
            if diverging:
                d = np.array(a, dtype=float)
                if mask is not None:
                    d[mask[rs, cs]] = np.nan                     # blanked, as in the positions GUI
                b['im3'] = ax.imshow(d, origin=display_look.ORIGIN, cmap='RdBu_r',
                                     vmin=-b['vmax3'], vmax=b['vmax3'],
                                     extent=ext, interpolation='nearest')
            else:
                ax.imshow(display_look.normalised_with(a, norm_snr), origin=display_look.ORIGIN,
                          cmap=display_look.CMAP, vmin=0.0, vmax=1.0, extent=ext,
                          interpolation='nearest')
            if mask is not None and mask[rs, cs].any():
                edge = mask_boundary(mask[rs, cs])
                # black line on a thin white halo: red blended into inferno's bright end, and a
                # bare black line would vanish on its dark sky (the halo is the positions
                # crosses' convention)
                from scipy.ndimage import binary_dilation
                halo = binary_dilation(edge, iterations=1) & ~edge
                rgba = np.zeros(edge.shape + (4,))
                rgba[halo] = (1.0, 1.0, 1.0, 0.7)
                rgba[edge] = (0.0, 0.0, 0.0, 0.95)
                ax.imshow(rgba, origin=display_look.ORIGIN, extent=ext, interpolation='nearest')
            if zoom:
                xl, yl = (-half_zoom, half_zoom), (-half_zoom, half_zoom)
            else:
                h = (n / 2.0) * s
                xl, yl = (-h, h), (-h, h)
            if i == 2 and te:
                cx0, cy0 = b.get('dp_ws', (0.0, 0.0))
                t = np.linspace(0, 2 * np.pi, 200)
                ax.plot(cx0 + te * np.cos(t), cy0 + te * np.sin(t), ls='--', lw=0.8, color='0.3',
                        alpha=0.9)
            # extra-galaxy centres, every row (autoarray (y, x) arcsec about the stamp centre,
            # which plots at (x, y) in these axes -- see the extent comment)
            for y, x in b['extra']:
                ax.plot(x, y, 'x', color=EXTRA_COLOUR, ms=8, mew=1.8, path_effects=HALO, zorder=6)
            ax.set_xlim(*xl)
            ax.set_ylim(*yl)
            ax.set_xticks([]); ax.set_yticks([])
            if i == 0:
                star = ' ★' if f == primary else ''
                sc = scores.get(f) or {}
                meas = (f'arc S/N {sc.get("arc_snr")}  score {sc.get("score")}'
                        if sc.get('arc_snr') is not None else 'arc S/N not measured')
                ax.set_title(f'{f}{star}\n{meas}\n{s:.3f}"/px', fontsize=10,
                             color='darkred' if f == primary else 'black')
            if j == 0:
                ax.set_ylabel(['S/N, full stamp', 'S/N, zoom on the lens',
                               'residual / noise (zoom, smoothed)'][i], fontsize=9)
            # 1" scale bar, bottom right (x limits are normal now: +x = West to the right)
            (xa, xb), (ya, yb) = xl, yl
            wx, wy = xb - xa, yb - ya
            bx, by = xa + 0.93 * wx, ya + 0.08 * wy
            ink = 'black' if i == 2 else 'white'
            ax.plot([bx, bx - 1.0], [by, by], color=ink, lw=2)
            ax.text(bx - 0.5, by + 0.04 * wy, '1"', color=ink, ha='center', fontsize=8)
    # colourbars: every S/N PANEL gets its own bar, spanning the one shared stretch (so colours
    # stay comparable across bands) but ticked at that panel's own min-med-max, where a panel's
    # "max" is its own 99th percentile -- the stretch's statistic, so the deepest band's ticks
    # span the bar and every other band's max tick lands partway down it at its own value.
    # Residual panels keep their own per-band symmetric sigma bars.
    sm = matplotlib.cm.ScalarMappable(cmap=display_look.CMAP, norm=norm_snr)
    fmt = FuncFormatter(lambda v, _pos: f'{v:.3g}')
    lo, hi = norm_snr.vmin, norm_snr.vmax
    for j, f in enumerate(filts):
        b = bands[f]
        for i in (0, 1):
            cb = fig.colorbar(sm, cax=fig.add_subplot(gs[i, 2 * j + 1]))
            a = b['snr'] if i == 0 else b['snr'][b['zoom']]
            p99 = float(np.nanpercentile(a, 99))
            if np.isfinite(p99) and p99 > lo:
                top = min(p99, hi)
                # the mid tick: the full stamp is sky-dominated (median ~0), so it takes the
                # arithmetic midpoint of the panel's own range; in the ZOOM every band's core
                # can exceed the shared stretch top (a bright galaxy: all bands' zoom p99 > hi
                # pinned every bar identically), so a saturated panel takes its own MEDIAN --
                # where the bulk of the zoom sits is what discriminates the bands there
                if i == 1 and p99 > hi:
                    mid = float(np.nanpercentile(a, 50))
                    if not lo + 0.03 * (hi - lo) < mid < top - 0.03 * (hi - lo):
                        mid = None
                else:
                    mid = (lo + top) / 2.0
                ticks = [lo] + ([mid] if mid is not None else []) + [top]
                labels = [f'{v:.3g}' for v in ticks]
                if p99 > hi:
                    labels[-1] = f'>{p99:.3g}'     # the panel's own p99, beyond the bar top
                cb.set_ticks(ticks)
                cb.ax.set_yticklabels(labels)
            cb.ax.tick_params(labelsize=6)
    for j, f in enumerate(filts):
        cb = fig.colorbar(bands[f]['im3'], cax=cax3[j])
        cb.ax.tick_params(labelsize=6)
    # compass in the figure's top-left corner, hugging the edge so the (centred, long) suptitle
    # lines stay clear of it: N DOWN, E left (the repos' convention, like every other display)
    from matplotlib.patches import FancyArrowPatch
    W, H = fig.get_size_inches()
    ox, oy, L = 0.02, 0.93, 0.02                       # figure fractions; L along x
    for dx, dy, lab, ha, va in ((0.0, -L * W / H, 'N', 'center', 'top'), (-L, 0.0, 'E', 'right', 'center')):
        fig.add_artist(FancyArrowPatch((ox, oy), (ox + dx, oy + dy), transform=fig.transFigure,
                                       arrowstyle='->', mutation_scale=12, color='black', lw=1.3))
        fig.text(ox + 1.15 * dx, oy + 1.15 * dy, lab, ha=ha, va=va, fontsize=9)
    te_s = f'{te:.2f}"' if te else '?'
    reviewed = pb.get('reviewed')
    rev_s = '' if psource == 'visual' else (', reviewed' if reviewed else ', UNREVIEWED')
    prim = f'{primary} ({psource}{rev_s})' if primary else 'not set (no band measurable)'
    dp0 = bands[filts[0]].get('dp_ws')
    centre_s = ('zooms on the stamp centre (the shared deflector/group centre)' if dp0 is None
                else 'zooms on the stamp centre; the deflector sits at '
                     f'{1e3 * -dp0[0]:+.0f} mas E, {1e3 * -dp0[1]:+.0f} mas N of it '
                     f'({"measured" if bands[filts[0]].get("dp_exact") else "catalogue"})')
    fig.suptitle(f'{lens}  ({sample}{", " + cat if cat else ""})   z_l {zl_s}   z_s {zs_s}   θ_E {te_s}\n'
                 f'primary band: {prim}  (★)\n'
                 f'panels are S/N (data / noise) on ONE stretch shared by every band (per-panel colourbars); residual row per-band ±sigma\n'
                 f'{centre_s}\n'
                 f'black outline: contaminant mask   dashed: θ_E'
                 + ('   cyan x: extra-galaxy centres' if any(bands[f]['extra'] for f in filts) else ''),
                 fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    return fig, dict(filts=filts, axes=axes, lens_dir=lens_dir, prefix=prefix, primary=primary)


def render(sample, lens, out=None):
    """Write the lens's stamp<S>_bands.png. Returns its path, or None without stamps."""
    fig, info = build_figure(sample, lens)
    if fig is None:
        return None
    out = out or os.path.join(info['lens_dir'], f"{info['prefix']}_bands.png")
    fig.savefig(out, dpi=90)
    return out


def rerender(sample, lens):
    """Best-effort re-render for the GUIs' hooks: never fatal, a hand-drawn product must not be
    lost to a plotting error."""
    try:
        out = render(sample, lens)
        if out:
            print(f'  band overview: {os.path.relpath(out, ws_path)}')
    except Exception as exc:                              # deliberately broad -- see docstring
        print(f'  NOTE: band overview not re-rendered for {lens} ({exc!r}); run '
              f'uv run python scripts/make_band_overview.py --sample {sample} --lens {lens}')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--sample', type=mast_target_names.sample_arg, default=None,
                   help='default: every sample')
    p.add_argument('--lens', nargs='*', default=None, help='one or more lenses (default: all)')
    a = p.parse_args()
    for lens in a.lens or []:                 # each lens's sample is found; they may differ here
        mast_target_names.resolve_sample(p, argparse.Namespace(lens=lens, sample=a.sample))
    warnings.simplefilter('ignore')
    for sample in ([a.sample] if a.sample else mast_target_names.samples()):
        for lens in mast_target_names.lenses(sample):
            if a.lens and lens not in a.lens:
                continue
            out = render(sample, lens)
            print(f'{lens}: ' + (os.path.relpath(out, ws_path) if out else 'no stamps'))


if __name__ == '__main__':
    sys.exit(main())
