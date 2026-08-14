#!/usr/bin/env python
"""
Interactive GUI mask-making, cycling through a sample's cutouts.

For each (lens, filt) cutout on disk, launches PyAutoLens's `Scribbler` GUI over the
cutout's science image -- the same tool as
autolens_workspace:scripts/imaging/data_preparation/gui/mask.py -- so you can scribble the
region to keep, then writes the result as cutout_[cr_]mask.fits alongside that cutout's
sci/noise/psf products.

Defaults to the pipeline's standard 12" cutout tree (data/cutouts/, cutout_paths.DEFAULT_SIZE)
-- the tree these hand-drawn masks are meant for and the only size-variant tree tracked in
git (see .gitignore) precisely because it carries them, which no script can regenerate.
Pass --size 20 to mask the (untracked, regenerable) 20" tree instead.

bcfill variant (--variant, default 'auto'). The bad-column-filled re-drizzle
(data/cutouts_bcfill/, also tracked) shares crop geometry with the standard tree EXACTLY
-- identical NAXIS/CRPIX/CRVAL -- so a mask drawn on one is pixel-valid on the other, and
bcfill differs only by having its dead-column noise stripes filled (a cleaner image to
scribble on). It covers ACS (f814W/f555W) + WFPC2 (f606W) only -- there is no bcfill for
f160W (WFC3/IR has no bad columns) or the gallery sample. So the default is:
    --variant auto  (default): per (lens, filt), a single mask is drawn on and written to
        the bcfill cutout if it exists, else the standard cutout. Because the geometry is
        shared, that one mask serves whichever reduction is modelled -- it is NOT duplicated
        across trees. A (lens, filt) is SKIPPED if a mask already exists in EITHER tree
        (bcfill or standard); with --force it is redrawn, again into the priority tree.
    --variant bcfill: bcfill tree only (skip lenses/filters with no bcfill cutout).
    --variant standard: standard tree only (the historic behaviour, data/cutouts/).

This is a manual, one-image-at-a-time tool (not a batch driver): each cutout blocks on its
own Tk window until you press Esc. Already-masked cutouts are skipped so a run can be
resumed across lenses; --force redraws.

Usage:
    uv run python scripts/make_masks.py --sample slacs_gold
    uv run python scripts/make_masks.py --lens J0008-0004 --filt f814W
    uv run python scripts/make_masks.py --sample slacs_gold --filt f814W   # one band, every lens
    uv run python scripts/make_masks.py --sample slacs_gold --force
"""

import argparse
import glob
import os
import sys
import warnings

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

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

MASKS_JSON = os.path.join(ws_path, 'info', 'lens_masks.json')


def pixel_scale_from_header(hdr):
    """arcsec/pixel from a cutout FITS header's WCS -- Array2D.from_fits wants a plain
    float, not the WCS itself."""
    return float(proj_plane_pixel_scales(WCS(hdr).celestial)[0] * 3600.0)


def find_prefix(cutout_dir, drizzle_pass='auto'):
    """Pick 'cutout_cr' or 'cutout', mirroring make_cutouts.py's --pass auto logic: prefer
    the CR-rejected pass, fall back to no-CR (e.g. F160W, which has no CR pass). Returns
    None if the requested pass has no sci file in this cutout_dir.
    """
    has_cr = os.path.exists(os.path.join(cutout_dir, 'cutout_cr_sci.fits'))
    has_nocr = os.path.exists(os.path.join(cutout_dir, 'cutout_sci.fits'))
    if drizzle_pass == 'cr':
        return 'cutout_cr' if has_cr else None
    if drizzle_pass == 'nocrrej':
        return 'cutout' if has_nocr else None
    return 'cutout_cr' if has_cr else ('cutout' if has_nocr else None)


def discover_targets(trees, sample, lens=None, filt=None):
    """Yield (lens, filt, display_dir, write_dirs) per lens/filt that has a cutout in any
    of `trees`, sorted for a reproducible run order.

    `trees` is an ordered list of (variant, cutouts_root_dir) in display-preference order
    (bcfill before standard for --variant auto). For each (lens, filt) present in any tree:
      - write_dirs   = every tree's cutout_dir that has a matching sci file (both trees
                       share crop geometry exactly, so one mask is valid for all of them),
                       tagged with its variant;
      - display_dir  = the first such tree in preference order (the image scribbled on).
    """
    seen = {}
    for variant, root in trees:
        pattern = os.path.join(root, sample, lens or '*', filt or '*')
        for cutout_dir in sorted(glob.glob(pattern)):
            if find_prefix(cutout_dir, 'auto') is None:
                continue
            this_filt = os.path.basename(cutout_dir)
            this_lens = os.path.basename(os.path.dirname(cutout_dir))
            seen.setdefault((this_lens, this_filt), []).append((variant, cutout_dir))
    for (this_lens, this_filt) in sorted(seen):
        write_dirs = seen[(this_lens, this_filt)]        # already in preference order
        display_variant, display_dir = write_dirs[0]
        yield this_lens, this_filt, display_dir, write_dirs


def make_mask_for(display_dir, write_dirs, lens, filt, sample, drizzle_pass, force,
                  brush_width):
    """Run the Scribbler GUI once and write ONE mask FITS, into the priority tree only
    (`write_dirs[0]` -- bcfill where its cutout exists, else standard).

    The trees share crop geometry exactly, so a single mask serves whichever reduction is
    modelled; it is not duplicated across trees. Skips (draws nothing) if a mask already
    exists in ANY of `write_dirs`, unless --force. Returns True if a mask was written.
    """
    display_prefix = find_prefix(display_dir, drizzle_pass)
    if display_prefix is None:
        print(f"{lens} {filt}: no cutout sci file for --pass {drizzle_pass}, skipping")
        return False

    # Existing masks anywhere across the trees -> skip (one mask suffices for both).
    existing = []
    for variant, cutout_dir in write_dirs:
        prefix = find_prefix(cutout_dir, drizzle_pass)
        if prefix is None:
            continue
        if os.path.exists(os.path.join(cutout_dir, f'{prefix}_mask.fits')):
            existing.append(variant or 'standard')
    if existing and not force:
        print(f"{lens} {filt}: mask already exists in [{', '.join(existing)}], "
              f"skipping (--force to redraw)")
        return False

    # Draw on, and write to, the priority tree only.
    write_variant = write_dirs[0][0] or 'standard'
    mask_path = os.path.join(display_dir, f'{display_prefix}_mask.fits')
    sci_path = os.path.join(display_dir, f'{display_prefix}_sci.fits')
    with fits.open(sci_path) as hdul:
        sci_hdr = hdul[0].header
    pixel_scales = pixel_scale_from_header(sci_hdr)

    print(f"\n{lens} {filt}  [{display_prefix} pass, {write_variant}]  "
          f"pixel_scale={pixel_scales:.4f}\"/pix")
    print(f"  {sci_path}")
    print("  Scribbler GUI: scribble the region to KEEP in the fit, press Esc when done.")

    data = al.Array2D.from_fits(file_path=sci_path, pixel_scales=pixel_scales)
    scribbler = al.Scribbler(image=data.native, brush_width=brush_width)
    scribbled = scribbler.show_mask()
    # Scribbler marks the scribbled (kept) region True; al.Mask2D's convention is the
    # opposite -- True means excluded from the fit -- so invert, exactly as
    # autolens_workspace:scripts/imaging/data_preparation/gui/mask.py does.
    mask = al.Mask2D(mask=np.invert(scribbled), pixel_scales=pixel_scales)

    aplt.fits_array(array=mask, file_path=mask_path, overwrite=True)
    print(f"  wrote {mask_path}")

    info_json.update(MASKS_JSON, sample, lens, filt, {
        'prefix': display_prefix,
        'drizzle_pass': 'cr' if display_prefix == 'cutout_cr' else 'nocrrej',
        'pixel_scale_arcsec': round(pixel_scales, 6),
        'brush_width': brush_width,
        'variant': write_variant,
    })
    return True


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--sample', default=mast_target_names.DEFAULT_SAMPLE,
                   help=f'sample subdirectory of data/cutouts/ (default '
                        f'{mast_target_names.DEFAULT_SAMPLE})')
    p.add_argument('--lens', default=None,
                   help='restrict to one lens; default every lens with cutouts in --sample')
    p.add_argument('--filt', default=None,
                   help='restrict to one filter; default every filter cutout the lens has')
    p.add_argument('--pass', dest='drizzle_pass', choices=['auto', 'cr', 'nocrrej'],
                   default='auto',
                   help="which cutout pass to mask, matching make_cutouts.py's --pass: "
                        "'auto' (default) prefers cutout_cr_*, falling back to cutout_* "
                        "where no CR pass exists (F160W)")
    p.add_argument('--size', type=float, default=cutout_paths.DEFAULT_SIZE,
                   help=f'cutout tree to draw masks for (default {cutout_paths.DEFAULT_SIZE:g}", '
                        'i.e. data/cutouts/ -- the tree this tool is meant for and '
                        'the only size-variant tree tracked in git, see .gitignore; pass '
                        '--size 20 for the untracked, regenerable 20" tree)')
    p.add_argument('--variant', choices=['auto', 'bcfill', 'standard'], default='auto',
                   help="which reduction's cutouts to mask: 'auto' (default) draws one mask "
                        "on the bcfill sci where it exists (cleaner image, shared geometry) "
                        "else the standard sci, skipping if a mask exists in EITHER tree; "
                        "'bcfill' uses only data/cutouts_bcfill/; 'standard' uses only "
                        "data/cutouts/")
    p.add_argument('--force', action='store_true', default=False,
                   help='redraw a mask that already exists (default: skip it)')
    p.add_argument('--brush-width', type=float, default=0.05,
                   help='Scribbler brush width, passed straight to al.Scribbler (default 0.05)')
    a = p.parse_args()

    # Ordered by display preference: bcfill first (cleaner image), standard second.
    if a.variant == 'bcfill':
        variants = ['bcfill']
    elif a.variant == 'standard':
        variants = ['']
    else:  # auto
        variants = ['bcfill', '']
    trees = [(v, cutout_paths.cutouts_root(ws_path, a.size, variant=v)) for v in variants]

    targets = list(discover_targets(trees, a.sample, a.lens, a.filt))
    if not targets:
        roots = ', '.join(r for _, r in trees)
        raise SystemExit(f"no cutouts found under [{roots}] for sample {a.sample} matching "
                         f"lens={a.lens!r} filt={a.filt!r}")

    print(f"{len(targets)} lens/filter cutout(s) to process")
    made = skipped = 0
    for lens, filt, display_dir, write_dirs in targets:
        if make_mask_for(display_dir, write_dirs, lens, filt, a.sample, a.drizzle_pass,
                         a.force, a.brush_width):
            made += 1
        else:
            skipped += 1

    print(f"\nDone: {made} mask(s) drawn, {skipped} skipped")


if __name__ == '__main__':
    main()
