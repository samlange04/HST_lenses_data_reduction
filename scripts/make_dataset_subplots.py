#!/usr/bin/env python
"""
Render PyAutoLens `subplot_dataset` QC images from the cutout products.

For each (lens, filt) cutout, assembles a PyAutoLens `al.Imaging` dataset from the four
FITS products that a lens model actually consumes -- signal, noise, PSF, mask -- applies
the mask, and writes a 3x3 `aplt.subplot_imaging_dataset` PNG (data, data-log10, noise,
PSF, PSF-log10, S/N, and the two over-sample-size panels the mask induces). This is the
same dataset object and the same subplot a modelling run would build, so the image is a
direct QC of what will be fed to the fit -- crucially it shows the *mask* laid over the
data/noise/S/N, which none of the pipeline's other QC PNGs do.

All four inputs must exist for a (lens, filt) or it is skipped with a note listing what is
missing -- a subplot cannot be built from a partial dataset. In particular the MASK is
hand-drawn by scripts/make_masks.py and is not regenerable, so a lens that has not been
masked yet is simply reported and skipped, not an error.

Variant / component mixing (--variant, default 'auto'). The standard (data/cutouts/) and
bad-column-filled (data/cutouts_bcfill/) reductions share crop geometry EXACTLY (identical
NAXIS/CRPIX/CRVAL -- see scripts/cutout_paths.py and make_masks.py), so their arrays are
pixel-aligned and a component may be taken from whichever tree has it:
  - the hand-drawn MASK lives in the bcfill tree for ACS/WFPC2 lenses (make_masks.py's
    default), while
  - the PSF kernel is stored ONCE per band, not per tree: cutout_paths.psf_cutout_dir puts
    it in the bcfill cutout dir wherever that reduction exists (ACS f814W/f555W, WFPC2
    f606W) and in data/cutouts/ otherwise (f160W, gallery).
`--variant auto` (default) therefore resolves EACH of sci/noise/mask independently from an
ordered tree search (bcfill first, then standard); `--variant bcfill`/`standard` restrict
that search to one tree. The psf is resolved separately, through
cutout_paths.psf_cutout_dir, and so obeys neither --variant nor --size: a single stored
kernel is the same file whichever tree holds it, and restricting would report it missing
rather than reading the one that exists.

Output: `{prefix}_dataset.png` written into the cutout dir the SIGNAL was taken from (so it
sits beside that reduction's sci/noise), where `{prefix}` is `cutout_cr` or `cutout` exactly
as make_cutouts.py named the products. Provenance (which tree each of the four components
came from, pixel scale, masked-pixel fraction) is recorded per (sample, lens, filt) in
info/lens_dataset_subplots.json.

Usage:
    uv run python scripts/make_dataset_subplots.py --lens J0008-0004 --filt f814W
    uv run python scripts/make_dataset_subplots.py --sample slacs_gold
    uv run python scripts/make_dataset_subplots.py --sample slacs_gold --filt f814W
    uv run python scripts/make_dataset_subplots.py --sample slacs_gold --force
"""

import argparse
import contextlib
import glob
import os
import sys
import warnings

import numpy as np
from astropy.io import fits
from astropy.visualization import AsinhStretch, ImageNormalize, PercentileInterval
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from matplotlib.colors import LogNorm

# autonerves prints a workspace-version-mismatch UserWarning on import in this repo (it's
# not an autolens_workspace checkout) -- harmless, silence it so it doesn't repeat per lens.
warnings.filterwarnings('ignore', category=UserWarning, module='autonerves')

import autolens as al
import autolens.plot as aplt

ws_path = '/Users/samlange/Code/HST_lenses_data_reduction'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mast_target_names
import info_json
import cutout_paths

SUBPLOTS_JSON = os.path.join(ws_path, 'info', 'lens_dataset_subplots.json')

# The four FITS products a lens model consumes; all required to build the dataset.
COMPONENTS = ('sci', 'noise', 'psf', 'mask')

# Default appearance: match the pipeline's own QC plots (make_cutouts.py / make_mosaics.py),
# which use inferno + PercentileInterval(99.0) + AsinhStretch(0.1). PyAutoLens's vendored
# plot_array only knows linear / LogNorm norms, so the asinh stretch is injected via a
# context-managed Axes.imshow patch (asinh_norm_patch) rather than a public argument.
DEFAULT_CMAP = 'inferno'
ASINH_PCT = 99.0


@contextlib.contextmanager
def asinh_norm_patch(pct=ASINH_PCT):
    """Temporarily make matplotlib's ``Axes.imshow`` apply the same asinh stretch the
    pipeline's own plots use (``PercentileInterval(pct)`` limits + ``AsinhStretch(0.1)``),
    computed per panel from the array being shown.

    Only the linear panels are affected -- ``aplt.subplot_imaging_dataset``'s Data / Noise /
    S/N / over-sample panels, which autolens draws with matplotlib's default linear norm.
    The dedicated log10 panels (drawn with a ``LogNorm``) are left untouched, and so are RGB
    images and the translucent ``array_overlay`` imshow, so those keep their intended look.
    """
    from matplotlib.axes import Axes
    orig_imshow = Axes.imshow

    def patched(self, X, *args, **kwargs):
        norm = kwargs.get('norm', None)
        arr = np.asarray(X)
        skip = (
            isinstance(norm, LogNorm)             # log10 panels: keep the log stretch
            or arr.ndim != 2                      # RGB image
            or kwargs.get('alpha') is not None    # translucent array_overlay
            or kwargs.get('cmap') == 'Greys'      # array_overlay
        )
        if not skip:
            finite = arr[np.isfinite(arr)]
            if finite.size:
                vmin, vmax = PercentileInterval(pct).get_limits(finite)
                kwargs['norm'] = ImageNormalize(vmin=vmin, vmax=vmax,
                                                stretch=AsinhStretch(0.1))
                kwargs.pop('vmin', None)
                kwargs.pop('vmax', None)
        return orig_imshow(self, X, *args, **kwargs)

    Axes.imshow = patched
    try:
        yield
    finally:
        Axes.imshow = orig_imshow


def pixel_scale_from_header(hdr):
    """arcsec/pixel from a cutout FITS header's WCS (a plain float, as Array2D wants)."""
    return float(proj_plane_pixel_scales(WCS(hdr).celestial)[0] * 3600.0)


def find_prefix(cutout_dir, drizzle_pass='auto'):
    """Pick 'cutout_cr' or 'cutout', mirroring make_cutouts.py / make_masks.py --pass logic:
    prefer the CR-rejected pass, fall back to no-CR (e.g. F160W, no CR pass). Returns None
    if the requested pass has no sci file here.
    """
    has_cr = os.path.exists(os.path.join(cutout_dir, 'cutout_cr_sci.fits'))
    has_nocr = os.path.exists(os.path.join(cutout_dir, 'cutout_sci.fits'))
    if drizzle_pass == 'cr':
        return 'cutout_cr' if has_cr else None
    if drizzle_pass == 'nocrrej':
        return 'cutout' if has_nocr else None
    return 'cutout_cr' if has_cr else ('cutout' if has_nocr else None)


def discover_targets(trees, sample, lens=None, filt=None):
    """Yield (lens, filt, {variant: cutout_dir}) for every (lens, filt) present in any tree,
    sorted for a reproducible run order. `trees` is an ordered list of (variant, root) in
    search-preference order (bcfill before standard for --variant auto).
    """
    seen = {}
    for variant, root in trees:
        pattern = os.path.join(root, sample, lens or '*', filt or '*')
        for cutout_dir in sorted(glob.glob(pattern)):
            if not os.path.isdir(cutout_dir):
                continue
            if find_prefix(cutout_dir, 'auto') is None:
                continue
            this_filt = os.path.basename(cutout_dir)
            this_lens = os.path.basename(os.path.dirname(cutout_dir))
            # first tree in preference order wins for a given variant key
            seen.setdefault((this_lens, this_filt), {}).setdefault(variant, cutout_dir)
    for key in sorted(seen):
        yield key[0], key[1], seen[key]


def resolve_components(dirs_by_variant, tree_order, drizzle_pass, psf_dir=None):
    """Locate each of the four components across the trees in preference order.

    `dirs_by_variant` maps variant -> cutout_dir for this (lens, filt). `tree_order` is the
    ordered list of variant keys to search. The drizzle prefix (cr/nocrrej) is fixed once,
    from the first tree that has a sci for `drizzle_pass`, so all four components come from
    the same pass. Returns (prefix, {component: (path, variant)}, missing_list).

    `psf_dir` (cutout_paths.psf_cutout_dir) overrides the tree search for the psf component
    only. The kernel is stored once per band -- variant-placed, and NOT size-keyed -- so it
    is read from wherever that rule put it rather than from this run's --variant/--size
    trees, which would otherwise report it missing whenever the two disagree.
    """
    # Fix the prefix from the first available tree (both trees share the same pass set).
    prefix = None
    for variant in tree_order:
        cutout_dir = dirs_by_variant.get(variant)
        if cutout_dir is None:
            continue
        prefix = find_prefix(cutout_dir, drizzle_pass)
        if prefix is not None:
            break
    if prefix is None:
        return None, {}, list(COMPONENTS)

    resolved, missing = {}, []
    for comp in COMPONENTS:
        found = None
        if comp == 'psf' and psf_dir is not None:
            path = os.path.join(psf_dir, f'{prefix}_psf.fits')
            if os.path.exists(path):
                found = (path, cutout_paths.psf_cutout_variant(psf_dir))
        else:
            for variant in tree_order:
                cutout_dir = dirs_by_variant.get(variant)
                if cutout_dir is None:
                    continue
                path = os.path.join(cutout_dir, f'{prefix}_{comp}.fits')
                if os.path.exists(path):
                    found = (path, variant or 'standard')
                    break
        if found is None:
            missing.append(comp)
        else:
            resolved[comp] = found
    return prefix, resolved, missing


def build_and_plot(lens, filt, sample, prefix, resolved, output_dir, check_noise_map,
                   asinh=True, cmap=DEFAULT_CMAP):
    """Assemble the al.Imaging dataset, apply the mask, and write the subplot PNG. Returns
    a provenance dict (also recorded in the tracking JSON by the caller)."""
    sci_path = resolved['sci'][0]
    pixel_scales = pixel_scale_from_header(fits.getheader(sci_path))

    dataset = al.Imaging.from_fits(
        pixel_scales=pixel_scales,
        data_path=sci_path,
        noise_map_path=resolved['noise'][0],
        psf_path=resolved['psf'][0],
        check_noise_map=check_noise_map,
    )
    mask = al.Mask2D.from_fits(file_path=resolved['mask'][0], pixel_scales=pixel_scales)
    frac_masked = float(np.asarray(mask).astype(bool).mean())
    dataset = dataset.apply_mask(mask=mask)

    os.makedirs(output_dir, exist_ok=True)
    # asinh + inferno (the pipeline's own plot look) via the imshow patch; a no-op
    # nullcontext restores stock autolens linear/default when --stretch default is chosen.
    with (asinh_norm_patch() if asinh else contextlib.nullcontext()):
        aplt.subplot_imaging_dataset(
            dataset=dataset,
            output_path=output_dir,
            output_filename=f'{prefix}_dataset',
            output_format='png',
            colormap=cmap,
        )
    out_png = os.path.join(output_dir, f'{prefix}_dataset.png')
    print(f"  wrote {out_png}")
    return {
        'prefix': prefix,
        'drizzle_pass': 'cr' if prefix == 'cutout_cr' else 'nocrrej',
        'stretch': 'asinh' if asinh else 'default',
        'cmap': cmap,
        'pixel_scale_arcsec': round(pixel_scales, 6),
        'frac_masked': round(frac_masked, 6),
        'sci_variant': resolved['sci'][1],
        'noise_variant': resolved['noise'][1],
        'psf_variant': resolved['psf'][1],
        'mask_variant': resolved['mask'][1],
        'output': os.path.relpath(out_png, ws_path),
    }


def process(lens, filt, sample, dirs_by_variant, tree_order, drizzle_pass, force,
            check_noise_map, asinh=True, cmap=DEFAULT_CMAP):
    """Resolve components, build the dataset, and write the subplot. Returns True on write,
    False on skip (missing component or already present without --force)."""
    prefix, resolved, missing = resolve_components(
        dirs_by_variant, tree_order, drizzle_pass,
        psf_dir=cutout_paths.psf_cutout_dir(ws_path, sample, lens, filt))
    if prefix is None:
        print(f"{lens} {filt}: no cutout sci for --pass {drizzle_pass}, skipping")
        return False
    if missing:
        have = {c: resolved[c][1] for c in resolved}
        print(f"{lens} {filt} [{prefix}]: missing {missing} "
              f"(have {have or 'nothing'}) -- cannot build dataset, skipping")
        return False

    # Output lands beside the reduction the SIGNAL came from.
    output_dir = os.path.dirname(resolved['sci'][0])
    out_png = os.path.join(output_dir, f'{prefix}_dataset.png')
    if os.path.exists(out_png) and not force:
        print(f"{lens} {filt} [{prefix}]: {os.path.basename(out_png)} exists, "
              f"skipping (--force to redraw)")
        return False

    srcs = ', '.join(f"{c}:{resolved[c][1]}" for c in COMPONENTS)
    print(f"\n{lens} {filt} [{prefix}]  {srcs}")
    provenance = build_and_plot(lens, filt, sample, prefix, resolved, output_dir,
                                check_noise_map, asinh=asinh, cmap=cmap)
    info_json.update(SUBPLOTS_JSON, sample, lens, filt, provenance)
    return True


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--sample', default=mast_target_names.DEFAULT_SAMPLE,
                   help=f'sample subdirectory of the cutout trees (default '
                        f'{mast_target_names.DEFAULT_SAMPLE})')
    p.add_argument('--lens', default=None,
                   help='restrict to one lens; default every lens with cutouts in --sample')
    p.add_argument('--filt', default=None,
                   help='restrict to one filter; default every filter the lens has')
    p.add_argument('--pass', dest='drizzle_pass', choices=['auto', 'cr', 'nocrrej'],
                   default='auto',
                   help="which cutout pass to use, matching make_cutouts.py's --pass: "
                        "'auto' (default) prefers cutout_cr_*, falling back to cutout_* "
                        "where no CR pass exists (F160W)")
    p.add_argument('--size', type=float, default=cutout_paths.DEFAULT_SIZE,
                   help=f'cutout tree to read (default {cutout_paths.DEFAULT_SIZE:g}", i.e. '
                        'data/cutouts/ -- the tracked, default-size tree carrying the '
                        'hand-drawn masks; pass --size 20 for the 20" tree)')
    p.add_argument('--variant', choices=['auto', 'bcfill', 'standard'], default='auto',
                   help="which reduction(s) to draw components from: 'auto' (default) "
                        "resolves each of sci/noise/psf/mask from bcfill where present else "
                        "standard (the two share crop geometry exactly, so mixing is valid); "
                        "'bcfill' uses only data/cutouts_bcfill/; 'standard' uses only "
                        "data/cutouts/")
    p.add_argument('--stretch', choices=['asinh', 'default'], default='asinh',
                   help="colour stretch for the linear panels (Data/Noise/S-N/over-sample): "
                        "'asinh' (default) matches the pipeline's own plots "
                        "(make_cutouts.py / make_mosaics.py: PercentileInterval(99) + "
                        "AsinhStretch(0.1)); 'default' restores stock autolens linear scaling. "
                        "The dedicated log10 panels stay log10 either way")
    p.add_argument('--cmap', default=DEFAULT_CMAP,
                   help=f'matplotlib colormap for every panel (default {DEFAULT_CMAP}, matching '
                        f'the pipeline plots)')
    p.add_argument('--force', action='store_true', default=False,
                   help='redraw a subplot PNG that already exists (default: skip it)')
    p.add_argument('--check-noise-map', action='store_true', default=False,
                   help="run al.Imaging's noise-map sanity check (off by default: this "
                        "pipeline's ERR/IVM noise maps are legitimate and can trip the "
                        "heuristic at low-coverage edges)")
    a = p.parse_args()

    # Search order: bcfill first (cleaner image, where masks live), standard second.
    if a.variant == 'bcfill':
        variants = ['bcfill']
    elif a.variant == 'standard':
        variants = ['']
    else:  # auto
        variants = ['bcfill', '']
    trees = [(v, cutout_paths.cutouts_root(ws_path, a.size, variant=v)) for v in variants]
    tree_order = [v for v, _ in trees]

    targets = list(discover_targets(trees, a.sample, a.lens, a.filt))
    if not targets:
        roots = ', '.join(r for _, r in trees)
        raise SystemExit(f"no cutouts found under [{roots}] for sample {a.sample} matching "
                         f"lens={a.lens!r} filt={a.filt!r}")

    print(f"{len(targets)} lens/filter cutout(s) to process")
    made = skipped = 0
    for lens, filt, dirs_by_variant in targets:
        if process(lens, filt, a.sample, dirs_by_variant, tree_order, a.drizzle_pass,
                   a.force, a.check_noise_map, asinh=(a.stretch == 'asinh'), cmap=a.cmap):
            made += 1
        else:
            skipped += 1

    print(f"\nDone: {made} subplot(s) written, {skipped} skipped")


if __name__ == '__main__':
    main()
