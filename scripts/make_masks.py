#!/usr/bin/env python
"""
Interactive GUI mask-making, cycling through a sample's cutouts.

For each (lens, filt) cutout on disk, launches PyAutoLens's `Scribbler` GUI over the
cutout's science image -- the same tool as
autolens_workspace:scripts/imaging/data_preparation/gui/mask.py -- so you can scribble the
region to keep, then writes the result as cutout_[cr_]mask.fits alongside that cutout's
sci/noise/psf products.

Defaults to the 12" cutout tree (data/cutouts_12arcsec/, --size 12) rather than the
pipeline's usual 20" default: it's the only size-variant tree tracked in git (see
.gitignore) precisely because it now carries these hand-drawn masks, which no script can
regenerate. Pass --size 20 to mask the 20" tree instead.

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


def discover_targets(cutouts_root_dir, sample, lens=None, filt=None):
    """Yield (lens, filt, cutout_dir) for every lens/filt under the sample that has at
    least one cutout sci file on disk, sorted for a reproducible run order.
    """
    pattern = os.path.join(cutouts_root_dir, sample, lens or '*', filt or '*')
    for cutout_dir in sorted(glob.glob(pattern)):
        if find_prefix(cutout_dir, 'auto') is None:
            continue
        this_filt = os.path.basename(cutout_dir)
        this_lens = os.path.basename(os.path.dirname(cutout_dir))
        yield this_lens, this_filt, cutout_dir


def make_mask_for(cutout_dir, lens, filt, sample, drizzle_pass, force, brush_width):
    """Run the Scribbler GUI on one cutout's science image and write its mask FITS.

    Returns True if a mask was (re)written, False if skipped.
    """
    prefix = find_prefix(cutout_dir, drizzle_pass)
    if prefix is None:
        print(f"{lens} {filt}: no cutout sci file for --pass {drizzle_pass}, skipping")
        return False

    sci_path = os.path.join(cutout_dir, f'{prefix}_sci.fits')
    mask_path = os.path.join(cutout_dir, f'{prefix}_mask.fits')
    if os.path.exists(mask_path) and not force:
        print(f"{lens} {filt}: {os.path.basename(mask_path)} already exists, "
              f"skipping (--force to redraw)")
        return False

    with fits.open(sci_path) as hdul:
        sci_hdr = hdul[0].header
    pixel_scales = pixel_scale_from_header(sci_hdr)

    print(f"\n{lens} {filt}  [{prefix} pass]  pixel_scale={pixel_scales:.4f}\"/pix")
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
        'prefix': prefix,
        'drizzle_pass': 'cr' if prefix == 'cutout_cr' else 'nocrrej',
        'pixel_scale_arcsec': round(pixel_scales, 6),
        'brush_width': brush_width,
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
    p.add_argument('--size', type=float, default=12.0,
                   help='cutout tree to draw masks for (default 12", i.e. '
                        'data/cutouts_12arcsec/ -- the tree this tool is meant for and '
                        'the only size-variant tree tracked in git, see .gitignore; pass '
                        f'--size {cutout_paths.DEFAULT_SIZE:g} for the default 20" tree)')
    p.add_argument('--force', action='store_true', default=False,
                   help='redraw a mask that already exists (default: skip it)')
    p.add_argument('--brush-width', type=float, default=0.05,
                   help='Scribbler brush width, passed straight to al.Scribbler (default 0.05)')
    a = p.parse_args()

    cutouts_root_dir = cutout_paths.cutouts_root(ws_path, a.size)
    targets = list(discover_targets(cutouts_root_dir, a.sample, a.lens, a.filt))
    if not targets:
        raise SystemExit(f"no cutouts found under {cutouts_root_dir}/{a.sample} matching "
                         f"lens={a.lens!r} filt={a.filt!r}")

    print(f"{len(targets)} lens/filter cutout(s) to process")
    made = skipped = 0
    for lens, filt, cutout_dir in targets:
        if make_mask_for(cutout_dir, lens, filt, a.sample, a.drizzle_pass, a.force,
                         a.brush_width):
            made += 1
        else:
            skipped += 1

    print(f"\nDone: {made} mask(s) drawn, {skipped} skipped")


if __name__ == '__main__':
    main()
