#!/usr/bin/env python
"""
Interactive GUI ARC mask-making: paint the arcs / multiple images, keep only those.

The second of the two hand-drawn mask tools, and the OPPOSITE polarity to `make_masks.py`:

    make_masks.py       paint what to REMOVE  -> cutout_[cr_]mask.fits       (contaminants)
    make_arc_masks.py   paint the ARCS to KEEP -> cutout_[cr_]mask_arcs.fits (this script)

Here you scribble ONLY the lensed source -- the arcs and the multiple images -- and
everything you do not paint is masked out, deflector light included. The product is a mask
that isolates the source: the region a source-only fit, an arc S/N measurement, or a
source-plane analysis should look at.

WHAT IS ON DISK, AND ITS POLARITY. The saved array is in `al.Mask2D`'s own convention --
True = EXCLUDED from the fit -- exactly like `cutout_[cr_]mask.fits`, so both files load
through the same reader and can never be mixed up by polarity:

    mask = al.Mask2D.from_fits(file_path=".../cutout_cr_mask_arcs.fits", pixel_scales=0.05)
    # mask == True  everywhere EXCEPT the arcs   -> fitting under it fits only the arcs
    # ~mask         is the arc region you painted (what MASKARC/NARCPX count)

So what you PAINT is inverted before writing (`arc_region = painted`, `saved = ~painted`).
The header records this explicitly: MASKTYPE='ARCS', MASKPOL, NARCPX. `make_masks.py`'s
scribble, by contrast, is saved verbatim -- there the painted region IS the excluded region.

THE DISPLAY SUBTRACTS THE DEFLECTOR BY DEFAULT (`--subtract-radial`, on here, off in
make_masks.py). The arcs are the subject of this mask and on most lenses they are invisible
under the galaxy envelope until its azimuthally-averaged radial profile is removed. It is a
DISPLAY transform only -- the mask is read back as brush pixel positions, so it cannot
change what a stroke masks -- and a finding aid, not photometry (an elliptical leaves a
quadrupole residual; an arc partly self-subtracts at its own radius). Pass
--no-subtract-radial to draw against the raw image.

ALREADY-MARKED POSITIONS ARE SHOWN. Where `make_positions.py` has recorded
cutout_[cr_]positions.json for the band, each marked image is marked with a dark CROSS in
the display so the arc mask can be drawn around the same multiple images. A cross, not a
ring: the other two annotations here are closed boundaries (solid white arc proposal, dashed
dark contaminant mask), and a third closed shape read as one more region masked OUT.
--no-show-positions turns it off.

YOU CORRECT A PROPOSAL, YOU DO NOT DRAW FROM SCRATCH (--propose-from, default 'auto').
`scripts/detect_arcs.py` finds the arcs automatically -- a PSF-matched colour difference
against the deflector (the source is blue, the deflector is red), thresholded inside the
measured Auger+2009 Einstein radius -- and its region is OUTLINED over both panels for review.
You then add with the green brush, erase with the red one, and on closing the GUI choose
[a]pply proposal+edits, keep only what you [d]rew, or [s]kip. It is the same review loop
make_masks.py uses for cross-band mask proposals, and the same `confirm_proposal` code: the
detector proposes, nothing is ever written unapproved, and a band left unapproved stays
pending.

SOURCES COMBINE BY UNION, so the second band of a lens starts from BOTH. 'auto' (the default)
= the detector PLUS the highest-priority other band that already has an arc mask: the first
band you draw gets the detector alone, and every band after it gets the detector plus the
region you yourself reviewed, with the per-source breakdown printed (`detect 1333px + f814W
2233px (overlap 1333px)`). They answer different questions and both belong on screen -- the
detector sees what is blue in this lens now, your own mask carries a judgement it cannot make,
and neither is authoritative, since an arc detected in one filter is often ABSENT in another
and its usable extent changes with each band's depth and PSF. That is exactly why the union is
a proposal, outlined for review and trimmed with the red brush, and never a broadcast. Takes
an explicit comma-separated list too ('detect,f814W', 'f814W'), or 'none' for a blank canvas.
Under --force the draw band's own mask joins the union, so refining a mask is not redrawing it.
--detect-* pass detector parameters through (--detect-snr is the one to reach for: lower it to
loosen a lens that proposes too little).

Everything else -- ONE DRAW PER BAND by default (--broadcast opts into sharing the draw
across the lens's bands by WCS reprojection), the band-priority order, --size tree routing,
the skip/--force behaviour -- is `make_masks.py`'s, imported from it rather than
re-implemented, so the two tools stay in lock-step. Where --broadcast IS used, note the
reprojection is applied to the ARC REGION (True = arc), not to the saved array:
outside-footprint pixels then default to 'not arc', which is the safe side -- reprojecting
the inverted array would instead leave stamp-edge slivers unmasked.

Usage:
    uv run python scripts/make_arc_masks.py --sample slacs_gold             # best band per lens
    uv run python scripts/make_arc_masks.py --lens J0330-0020               # one lens
    uv run python scripts/make_arc_masks.py --lens J0330-0020 --filt f606W  # force the band
    uv run python scripts/make_arc_masks.py --sample slacs_gold --filt f606W  # 2nd band sweep
    uv run python scripts/make_arc_masks.py --lens J0330-0020 --propose-from detect,f814W
    uv run python scripts/make_arc_masks.py --lens J0330-0020 --detect-snr 3.0  # looser detect
    uv run python scripts/make_arc_masks.py --lens J0330-0020 --propose-from none  # blank
    uv run python scripts/make_arc_masks.py --lens J0330-0020 --no-subtract-radial --force
    uv run python scripts/make_arc_masks.py --lens J0330-0020 --broadcast --propose-from none
"""

import argparse
import os
import sys
import warnings

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

warnings.filterwarnings('ignore', category=UserWarning, module='autonerves')

import autolens as al
import autolens.plot as aplt
from matplotlib import pyplot as plt

ws_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mast_target_names
import info_json
import cutout_paths
# Every display/tree/discovery/broadcast helper is make_masks'. This script contributes the
# opposite polarity, the arc-specific display defaults, and its own product name.
import make_masks
import detect_arcs

ARC_MASKS_JSON = os.path.join(ws_path, 'info', 'lens_arc_masks.json')

ARC_SUFFIX = 'mask_arcs'          # cutout_[cr_]mask_arcs.fits, beside cutout_[cr_]mask.fits

# Geometry of the position crosses (see cross_overlay): 3 px-thick arms starting 2 px out, so
# the marked pixel itself is never covered and the marker is about as many display pixels as
# the 6 px ring it replaced (~50 vs ~38) -- a cross made of hairlines reads as fainter than a
# ring of the same radius, and it has to survive being drawn over a bright arc.
_CROSS_HALF_WIDTH_PX = 1.0
_CROSS_GAP_PX = 2.0


def _one_arc_source(filt_dirs, prefixes, display_filt, dst_wcs, dst_shape, lens, source,
                    detect_kw=None):
    """Resolve ONE proposal source to (region, extra) on the draw band's grid, or (None, {}).

    'detect'  -- run the automatic detector (`detect_arcs.py`): a colour difference against the
                 deflector, thresholded inside the measured Einstein radius. Built on the
                 lens's RED band and reprojected here, so it lands on the same sky whichever
                 band is being drawn.
    '<band>'  -- an arc mask already drawn for another band of this lens, reprojected. The
                 stored array is inverted (True = excluded), so the ARC REGION is `~mask`, and
                 it is the region, never the saved array, that gets reprojected --
                 outside-footprint pixels then default to 'not arc', the safe side of the
                 stamp edge.
    """
    if source == 'detect':
        region, components, maps = detect_arcs.propose(filt_dirs, prefixes, lens,
                                                       **(detect_kw or {}))
        extra = {'detector_components': components,
                 'detector_bands': f"{maps['blue_filt']}-{maps['red_filt']}",
                 'detector_theta_e_arcsec': round(maps['theta_e'], 3),
                 'detector_band_offset_px': list(maps['band_offset_px']),
                 'detector_core_dipole': maps['core_dipole']}
        if maps['red_filt'] != display_filt:
            region = make_masks.reproject_mask_bool(
                region, WCS(maps['hdr']).celestial, dst_wcs, dst_shape)
        return region, extra

    cutout_dir = filt_dirs.get(source)
    prefix = prefixes.get(source)
    if cutout_dir is None or prefix is None:
        return None, {}
    mask_path = os.path.join(cutout_dir, f'{prefix}_{ARC_SUFFIX}.fits')
    if not os.path.exists(mask_path):
        return None, {}
    arc_region = ~np.asarray(fits.getdata(mask_path), dtype=bool)
    if source == display_filt:
        return arc_region, {}
    src_hdr = fits.getheader(os.path.join(cutout_dir, f'{prefix}_sci.fits'))
    return make_masks.reproject_mask_bool(arc_region, WCS(src_hdr).celestial,
                                          dst_wcs, dst_shape), {}


def resolve_arc_sources(filt_dirs, display_filt, propose_from, force=False):
    """Expand `--propose-from` into an ordered, de-duplicated list of concrete sources.

    'none'   -> []            (blank canvas)
    'auto'   -> ['detect'] + the highest-priority OTHER band that already has an arc mask,
                so the first band of a lens starts from the detector alone and every band
                after it starts from the detector PLUS your own reviewed work. Under --force
                the draw band's own existing arc mask comes first, so refining a mask is not
                redrawing it -- the same rule make_masks.py uses.
    explicit -> a comma-separated list, e.g. 'detect,f814W' or 'f814W' or 'f814W,detect'.
    """
    if propose_from in (None, 'none'):
        return []
    if propose_from != 'auto':
        return [t.strip() for t in str(propose_from).split(',') if t.strip()]
    others = sorted((f for f in filt_dirs if f != display_filt),
                    key=lambda f: (make_masks._band_rank(f), f))
    return (['detect'] + ([display_filt] if force else []) + others)


def find_arc_proposal(filt_dirs, prefixes, display_filt, display_hdr, lens, propose_from,
                      detect_kw=None, force=False):
    """Find an arc region to PROPOSE for the band about to be drawn, on that band's grid.

    Returns (source_label, bool array, extra) or (None, None, {}).

    SOURCES COMBINE BY UNION, because they answer different questions and you want both on
    screen at once. The detector sees what is blue in THIS lens right now; a mask you already
    drew for another band carries a judgement it cannot make. Neither is authoritative: an arc
    detected in one filter is often simply ABSENT in another, and its usable extent changes
    with each band's depth and PSF -- which is exactly why the union is a proposal, outlined
    for review and trimmed with the red brush, and never a broadcast.

    With 'auto' (the default) the first band of a lens gets the detector alone, and every band
    after it gets the detector plus your own reviewed region, with the per-source breakdown
    printed so you can see which part came from where.
    """
    sources = resolve_arc_sources(filt_dirs, display_filt, propose_from, force=force)
    if not sources:
        return None, None, {}
    dst_wcs = WCS(display_hdr).celestial
    dst_shape = (display_hdr['NAXIS2'], display_hdr['NAXIS1'])

    region = None
    used, per_source, extra = [], {}, {}
    for source in sources:
        try:
            part, part_extra = _one_arc_source(filt_dirs, prefixes, display_filt, dst_wcs,
                                               dst_shape, lens, source, detect_kw)
        except Exception as exc:        # one bad source must not cost the others
            print(f"  NOTE: proposal source {source!r} unavailable for {lens} ({exc})")
            continue
        if part is None or not part.any():
            continue
        used.append(source)
        per_source[source] = int(part.sum())
        extra.update(part_extra)
        region = part if region is None else (region | part)
        # 'auto' takes the detector plus ONE inherited band -- the best one that actually has
        # a mask. Stop as soon as an inherited source has contributed.
        if propose_from == 'auto' and source != 'detect':
            break
    if region is None:
        return None, None, {}
    if len(used) > 1:
        extra['proposal_px_by_source'] = per_source
        extra['proposal_overlap_px'] = int(sum(per_source.values()) - region.sum())
    return '+'.join(used), region, extra


def load_positions(cutout_dir, prefix):
    """Return the band's marked image positions as a list of (y, x) arcsec, or [] if
    make_positions.py has not run for it. Read-only -- this tool never writes positions."""
    path = os.path.join(cutout_dir, f'{prefix}_positions.json')
    if not os.path.exists(path):
        return []
    try:
        return [(float(p[0]), float(p[1])) for p in al.from_json(file_path=path)]
    except Exception as exc:                       # a hand-edited/partial file must not block
        print(f"  NOTE: could not read {os.path.basename(path)} ({exc}); no position crosses")
        return []


def cross_overlay(positions, cross_size_px=6):
    """Return an `overlay(disp_array, pixel_scales)` callable for make_masks.draw_mask_gui
    that marks each `make_positions.py` image with a CROSS in the DISPLAYED array.

    A CROSS, NOT A RING (2026-09-21). Every other annotation on this GUI is a closed
    boundary -- the arc proposal is a solid white outline and the contaminant mask is a
    dashed dark one -- so a dark ring read as a third small closed region, i.e. as something
    masked out, which is the opposite of what a marked image means. Four open ticks cannot be
    read as an enclosed area at all, so the annotations no longer have to be told apart by
    line style alone.

    Burned in at the display MINIMUM (dark), as the ring was: dark reads clearly against both
    the mid-tone background and the bright arcs of a radial-subtracted image, where a
    maximum-valued marker would be indistinguishable from arc flux -- the one thing it must
    not hide. The arms stop short of the centre (`_CROSS_GAP_PX`) for the same reason the
    contaminant outline is not filled: the marked pixel is the one you are judging, so the
    marker points at it rather than covering it.

    Display only: the Scribbler returns brush positions, so this cannot enter the mask.
    """
    def _overlay(disp_array, pixel_scales):
        a = np.array(disp_array.native, dtype=float)
        n_y, n_x = a.shape
        cy, cx = (n_y - 1) / 2.0, (n_x - 1) / 2.0
        yy, xx = np.mgrid[0:n_y, 0:n_x]
        dark = float(a.min())
        for (y_as, x_as) in positions:
            # arcsec -> pixel. +x arcsec is +column, but +y arcsec is -ROW: autoarray's
            # native grid puts row 0 at the TOP (+y), verified against
            # al.Grid2D.uniform(...).native[0, 0] == (+5.975, -5.975) for a 240px/0.05"
            # stamp. Getting this sign wrong silently marks the mirror image of each
            # position -- it lands plausibly near the lens and on nothing.
            py, px = cy - y_as / pixel_scales, cx + x_as / pixel_scales
            dy, dx = np.abs(yy - py), np.abs(xx - px)
            arm = float(cross_size_px)
            vertical = (dx <= _CROSS_HALF_WIDTH_PX) & (dy >= _CROSS_GAP_PX) & (dy <= arm)
            horizontal = (dy <= _CROSS_HALF_WIDTH_PX) & (dx >= _CROSS_GAP_PX) & (dx <= arm)
            a[vertical | horizontal] = dark
        return al.Array2D.no_mask(values=a, pixel_scales=pixel_scales).native
    return _overlay


def contaminant_overlay(mask_bool):
    """Return an `overlay(disp_array, pixel_scales)` that outlines the CONTAMINANT mask
    (`cutout_[cr_]mask.fits`, what make_masks.py excludes) in the displayed array.

    Why it is here: painting arcs and painting contaminants are opposite jobs, and the thing
    most likely to be mistaken for a lensed image is a neighbour or field source that has
    ALREADY been judged a contaminant on this very band. Showing that verdict while you paint
    means you are not re-deciding it from memory.

    DASHED and burned at the display MINIMUM, which keeps three annotations apart at a glance:
    the arc proposal is a solid WHITE outline (display max), marked image positions are DARK
    CROSSES, and this is a DASHED DARK boundary. The interior is deliberately NOT filled --
    the same call `make_masks.py` made and recorded: blanking a region hides the very pixels
    you are judging it against. Display only; the Scribbler returns brush positions, so
    nothing here can enter the saved mask.
    """
    def _overlay(disp_array, pixel_scales):
        a = np.array(disp_array.native, dtype=float)
        edge = make_masks.mask_boundary(np.asarray(mask_bool, dtype=bool))
        if edge.shape == a.shape:
            yy, xx = np.mgrid[0:a.shape[0], 0:a.shape[1]]
            a[edge & (((yy + xx) % 3) != 0)] = float(a.min())    # dashed, not solid
        return al.Array2D.no_mask(values=a, pixel_scales=pixel_scales).native
    return _overlay


def compose_overlays(*fns):
    """Chain overlay callables in order; None entries are ignored. Returns None if none remain,
    so draw_mask_gui's `overlay is not None` fast path is preserved."""
    fns = [f for f in fns if f is not None]
    if not fns:
        return None
    if len(fns) == 1:
        return fns[0]

    def _overlay(disp_array, pixel_scales):
        for fn in fns:
            disp_array = fn(disp_array, pixel_scales)
        return disp_array
    return _overlay


def save_arc_qc_png(cutout_dir, prefix, arc_bool, pixel_scales, out_path, lens, filt,
                    stretch='asinh', vmin_percent=5.0, vmax_percent=99.5, asinh_a=0.1):
    """QC overlay: the arc region outlined over the band's own (radial-subtracted) image.

    Written per band, so a broadcast can be checked where it actually matters -- the f160W
    regrid (0.06"/px), where an array-index copy would land the region ~20px off.
    """
    base, _ = make_masks.load_display_base(cutout_dir, prefix, 'snr', pixel_scales)
    disp = make_masks.stretched_display(make_masks.radial_median_subtract(base), pixel_scales,
                                        stretch=stretch, vmin_percent=vmin_percent,
                                        vmax_percent=vmax_percent, asinh_a=asinh_a)
    a = np.asarray(disp.native)
    hw = int(a.shape[1] / 2) * pixel_scales
    ext = [-hw, hw, -hw, hw]
    fig = plt.figure(figsize=(7, 7))
    plt.imshow(a, cmap='jet', extent=ext)
    # origin='upper' is REQUIRED: with an extent, imshow draws row 0 at the top but
    # contour defaults to row 0 at the BOTTOM, so without it the outline is drawn
    # vertically mirrored -- a QC image that looks fine and checks nothing (verified).
    plt.contour(np.asarray(arc_bool, dtype=float), levels=[0.5], colors='white',
                linewidths=1.2, extent=ext, origin='upper')
    plt.title(f'{lens}  {filt}  --  arc mask region (white outline)')
    plt.xlabel('arcsec'); plt.ylabel('arcsec')
    plt.tight_layout()
    fig.savefig(out_path, dpi=110)
    plt.close(fig)


def write_arc_mask(arc_bool, pixel_scales, path):
    """Write the arc mask in al.Mask2D convention (True = EXCLUDED, so only the painted arc
    region survives) and stamp the polarity into the header, since this file's meaning is
    the inverse of what was painted."""
    excluded = ~np.asarray(arc_bool, dtype=bool)
    aplt.fits_array(array=al.Mask2D(mask=excluded, pixel_scales=pixel_scales),
                    file_path=path, overwrite=True)
    with fits.open(path, mode='update') as hdul:
        hdr = hdul[0].header
        hdr['MASKTYPE'] = ('ARCS', 'hand-drawn lensed-source (arc/multiple-image) mask')
        hdr['MASKPOL'] = ('True=excluded', 'al.Mask2D convention; ~mask is the arc region')
        hdr['NARCPX'] = (int(np.asarray(arc_bool).sum()), 'pixels inside the arc region')


def process_lens_arc_mask(lens, filt_dirs, sample, requested_filt, drizzle_pass, force,
                          broadcast, brush_radius, brush_width, display, stretch,
                          vmin_percent, vmax_percent, asinh_a, subtract_radial,
                          side_by_side, show_positions, cross_size_px,
                          propose_from='auto', detect_kw=None, show_contaminants=True):
    """Draw the arc mask once for one lens, on its best/forced band only unless
    `broadcast` is set (then also WCS-reprojected to the lens's other bands). Mirrors
    make_masks.process_lens_mask; differs only in polarity, product name, and the arc-specific
    display (radial subtraction + position crosses). Returns True if a mask was written.

    `propose_from` ('detect' by default) starts the draw from a REVIEWED PROPOSAL rather than
    a blank canvas: the proposed arc region is outlined on every panel, you correct it with
    the green (add) / red (erase) brushes, and on closing the GUI you apply, reject or skip.
    'none' restores the from-scratch draw.
    """
    filts = sorted(filt_dirs)
    display_filt = make_masks.pick_display_filt(filts, requested_filt)
    if display_filt is None:
        print(f"{lens}: requested --filt {requested_filt} not present "
              f"(have {', '.join(filts)}), skipping")
        return False

    display_dir = filt_dirs[display_filt]
    display_prefix = make_masks.find_prefix(display_dir, drizzle_pass)
    if display_prefix is None:
        print(f"{lens} {display_filt}: no cutout sci for --pass {drizzle_pass}, skipping")
        return False

    # Skip if the draw band already has an arc mask, unless --force.
    if (os.path.exists(os.path.join(display_dir, f'{display_prefix}_{ARC_SUFFIX}.fits'))
            and not force):
        print(f"{lens}: arc mask already exists ({display_filt}), skipping (--force to redraw)")
        return False

    positions = load_positions(display_dir, display_prefix) if show_positions else []
    overlay = cross_overlay(positions, cross_size_px) if positions else None
    if positions:
        print(f"  {len(positions)} marked position(s) crossed in the display "
              f"(from {display_prefix}_positions.json)")
    contaminants = None
    if show_contaminants:
        cont_path = os.path.join(display_dir, f'{display_prefix}_mask.fits')
        if os.path.exists(cont_path):
            contaminants = np.asarray(fits.getdata(cont_path), dtype=bool)
            if contaminants.any():
                print(f"  {int(contaminants.sum())}px contaminant-masked on this band "
                      f"(dashed dark outline; make_masks.py's verdict, and the detector "
                      f"already excluded them)")
                overlay = compose_overlays(overlay, contaminant_overlay(contaminants))
            else:
                contaminants = None
        else:
            print(f"  NOTE: no {display_prefix}_mask.fits on this band -- no contaminant "
                  f"outline to show (draw one with make_masks.py first if you want it)")

    display_hdr = fits.getheader(os.path.join(display_dir, f'{display_prefix}_sci.fits'))
    prefixes = {f: make_masks.find_prefix(d, drizzle_pass) for f, d in filt_dirs.items()}
    proposal_from = proposal = None
    proposal_extra = {}
    try:
        proposal_from, proposal, proposal_extra = find_arc_proposal(
            filt_dirs, prefixes, display_filt, display_hdr, lens, propose_from, detect_kw,
            force=force)
    except Exception as exc:                  # a proposal is a convenience, never a blocker
        print(f"  NOTE: no proposal for {lens} ({exc}) -- drawing from scratch")
    if proposal is not None and not proposal.any():
        print(f"  the detector proposed nothing for {lens} {display_filt} -- drawing from "
              f"scratch (lower --detect-snr to loosen it)")
        proposal_from, proposal = None, None
    if proposal is not None:
        where = '; '.join(f"r={c['r_over_theta_e']:.2f}th_E pk={c['peak_snr']:.0f}"
                          for c in proposal_extra.get('detector_components', [])[:4])
        by_src = proposal_extra.get('proposal_px_by_source')
        breakdown = ('  ' + ' + '.join(f'{k} {v}px' for k, v in by_src.items())
                     + f"  (overlap {proposal_extra['proposal_overlap_px']}px)"
                     if by_src else '')
        print(f"  reviewing the {proposal_from} arc proposal ({int(proposal.sum())}px"
              + (f"; {where}" if where else '') + ") -- edit it, then apply/reject")
        if breakdown:
            print(breakdown)

    prompt = [
        "  Scribbler GUI: paint ONLY the arcs / multiple images of the source.",
        "    Everything you do NOT paint is masked OUT (the deflector included) -- this is",
        "    the OPPOSITE polarity to make_masks.py, which paints what to remove.",
        "    keys '1' = GREEN brush, paint arc   |   '2' = RED brush, un-paint (erase)",
    ]
    (painted, erased, src_ps, src_hdr, brush_width, start_radius,
     source_label) = make_masks.draw_mask_gui(
        display_dir, display_prefix, lens, display_filt, brush_radius, brush_width,
        display, stretch, vmin_percent, vmax_percent, asinh_a,
        subtract_radial=subtract_radial, side_by_side=side_by_side,
        prompt=prompt, overlay=overlay, proposal=proposal,
        proposal_label=(f'{proposal_from} arc region' if proposal_from else None))

    painted = np.asarray(painted, dtype=bool)
    erased = np.asarray(erased, dtype=bool)
    # painted = the arcs to KEEP; the red brush trims an over-painted stroke back off it.
    drawn_only = painted & ~erased
    n_add, n_erase = int(painted.sum()), int(erased.sum())
    if proposal is None:
        arc_region, source_tag = drawn_only, 'drawn'
    else:
        reviewed = (proposal | painted) & ~erased
        choice = make_masks.confirm_proposal(lens, display_filt, proposal_from,
                                             int(proposal.sum()), n_add, n_erase,
                                             int(reviewed.sum()), int(drawn_only.sum()))
        if choice == 'skip':
            print(f"  skipped -- no arc mask written for {lens} {display_filt} (still pending)")
            return False
        if choice == 'drawn':
            arc_region, source_tag = drawn_only, 'drawn'
            print(f"  {proposal_from} proposal REJECTED -- keeping only what you painted")
        else:
            arc_region = reviewed
            source_tag = (f'edited_from_{proposal_from}' if (n_add or n_erase)
                          else f'accepted_from_{proposal_from}')
    if not np.asarray(arc_region).any():
        print(f"  nothing painted -- no arc mask written for {lens} (an empty arc region "
              f"would mask out the entire stamp)")
        return False
    src_wcs = WCS(src_hdr).celestial

    targets = (filts if broadcast else [display_filt])
    for filt in targets:
        cutout_dir = filt_dirs[filt]
        prefix = make_masks.find_prefix(cutout_dir, drizzle_pass)
        if prefix is None:
            continue
        band_hdr = fits.getheader(os.path.join(cutout_dir, f'{prefix}_sci.fits'))
        band_ps = make_masks.pixel_scale_from_header(band_hdr)
        is_draw = (cutout_dir == display_dir)
        # Reproject the ARC REGION, never the inverted array: outside-footprint pixels then
        # default to 'not arc' (masked out), the safe side of the stamp edge.
        band_arc = arc_region if is_draw else make_masks.reproject_mask_bool(
            arc_region, src_wcs, WCS(band_hdr).celestial,
            (band_hdr['NAXIS2'], band_hdr['NAXIS1']))
        mask_path = os.path.join(cutout_dir, f'{prefix}_{ARC_SUFFIX}.fits')
        write_arc_mask(band_arc, band_ps, mask_path)
        save_arc_qc_png(cutout_dir, prefix, band_arc, band_ps,
                        os.path.join(cutout_dir, f'{prefix}_{ARC_SUFFIX}.png'), lens, filt,
                        stretch=stretch, vmin_percent=vmin_percent,
                        vmax_percent=vmax_percent, asinh_a=asinh_a)
        n_arc = int(np.asarray(band_arc).sum())
        print(f"  wrote {mask_path}  ({source_tag if is_draw else f'reprojected<-{display_filt}'}, "
              f"{n_arc}px in the arc region, everything else masked out)")
        entry = {
            'prefix': prefix,
            'drizzle_pass': 'cr' if prefix == 'cutout_cr' else 'nocrrej',
            'pixel_scale_arcsec': round(band_ps, 6),
            'display': source_label,
            # Which reduction the stamp was drawn on, from its own header (cutout_paths.py).
            'bcfill': bool(band_hdr.get('BCFILL', False)),
            'source': source_tag if is_draw else f'reprojected_from_{display_filt}',
            'n_arc_px': n_arc,
            'polarity': 'saved True=excluded (al.Mask2D); ~mask = arc region',
        }
        if is_draw:
            entry['brush_width'] = round(brush_width, 6)
            entry['brush_radius_px'] = start_radius
            entry['n_positions_shown'] = len(positions)
            entry['subtract_radial'] = bool(subtract_radial)
            if proposal is not None:
                # What the arc region inherited from the detector vs what the hand changed --
                # 'accepted_from_detect' with no edits is a real, deliberate outcome, and the
                # detector's own numbers are kept so a proposal can be reproduced or blamed.
                entry['proposal_from'] = proposal_from
                entry['proposal_px'] = int(proposal.sum())
                entry['n_added_px'] = n_add
                entry['n_erased_px'] = n_erase
                entry.update(proposal_extra)
        info_json.update(ARC_MASKS_JSON, sample, lens, filt, entry)
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
                   help='force which band you DRAW on (e.g. f606W); default the best available '
                        'band per lens (f814W>f606W>f555W>f160W>...). Only that band gets an '
                        'arc mask unless --broadcast is passed')
    p.add_argument('--broadcast', action=argparse.BooleanOptionalAction, default=False,
                   help="ALSO write the arc mask to the lens's other bands, WCS-reprojected "
                        'onto each grid. Default OFF: run again with --filt for the next band. '
                        'An arc is often detected in one filter and not another, and its usable '
                        'extent differs with each band\'s depth and PSF, so the drawn region is '
                        'not in practice band-independent')
    p.add_argument('--pass', dest='drizzle_pass', choices=['auto', 'cr', 'nocrrej'],
                   default='auto',
                   help="which cutout pass to mask, matching make_masks.py/make_cutouts.py: "
                        "'auto' (default) prefers cutout_cr_*, falling back to cutout_* where "
                        "no CR pass exists (F160W)")
    p.add_argument('--size', type=float, default=cutout_paths.DEFAULT_SIZE,
                   help=f'cutout tree to draw arc masks for (default '
                        f'{cutout_paths.DEFAULT_SIZE:g}", data/cutouts/ -- the git-tracked '
                        'tree; like every hand-drawn product these are non-regenerable)')
    p.add_argument('--force', action='store_true', default=False,
                   help='redraw an arc mask that already exists (default: skip it)')
    p.add_argument('--brush-radius', type=int, default=4,
                   help='starting brush radius in PIXELS (default 4 -- finer than make_masks.py, '
                        "since arcs are thin). Resize live in the GUI with '='/'-'")
    p.add_argument('--brush-width', type=float, default=None,
                   help='explicit brush width as a FRACTION of image height; overrides '
                        '--brush-radius when given')
    p.add_argument('--display', choices=['snr', 'sci'], default='snr',
                   help="image to scribble on: 'snr' (default) = signal/noise, which tracks "
                        "statistical significance given the non-uniform noise map; 'sci' = raw "
                        'signal. Choice does not affect the saved mask')
    p.add_argument('--stretch', choices=list(make_masks._STRETCHES), default='asinh',
                   help='intensity stretch for the DISPLAYED image only (default asinh)')
    p.add_argument('--asinh-a', type=float, default=0.1,
                   help='asinh softening (smaller = stronger stretch; --stretch asinh only)')
    p.add_argument('--vmin-percent', type=float, default=5.0,
                   help='lower clip percentile for the display stretch (default 5.0)')
    p.add_argument('--vmax-percent', type=float, default=99.5,
                   help='upper clip percentile for the display stretch (default 99.5)')
    p.add_argument('--subtract-radial', action=argparse.BooleanOptionalAction, default=True,
                   help="subtract the deflector's azimuthally-averaged radial profile from the "
                        'DISPLAYED image (DEFAULT ON here -- the arcs are this mask\'s subject '
                        'and are usually invisible under the galaxy envelope without it). '
                        'Display only; --no-subtract-radial draws against the raw image')
    p.add_argument('--side-by-side', action=argparse.BooleanOptionalAction, default=True,
                   help='show BOTH views at once -- radial-subtracted left, as-observed right '
                        '-- and accept scribbles on either panel, combined onto the one mask. '
                        'The subtracted panel is where the arcs are visible; the as-observed '
                        'one is where their real extent against the galaxy envelope is '
                        '(default on)')
    p.add_argument('--show-positions', action=argparse.BooleanOptionalAction, default=True,
                   help='mark the images already marked by make_positions.py with a dark '
                        'cross in the display, as a guide for where the multiple images are '
                        '(default on; silently inactive for a band with no positions file)')
    p.add_argument('--cross-size', type=int, default=6,
                   help='arm length in pixels of those position crosses (default 6)')
    p.add_argument('--show-contaminants', action=argparse.BooleanOptionalAction, default=True,
                   help="outline this band's CONTAMINANT mask (cutout_[cr_]mask.fits, what "
                        'make_masks.py excludes) as a dashed dark boundary, so a neighbour '
                        'already judged a contaminant is not mistaken for a lensed image '
                        '(default on; display only, silently inactive where no mask exists)')
    p.add_argument('--propose-from', default='auto',
                   help="where the arc region you review comes from. 'auto' (default) = the "
                        "automatic colour detector (scripts/detect_arcs.py) PLUS the "
                        "highest-priority other band of this lens that already has an arc "
                        "mask, UNIONed -- so the first band you draw starts from the detector "
                        "alone and every band after it starts from the detector plus your own "
                        "reviewed work. Also takes an explicit comma-separated list "
                        "('detect,f814W', 'f814W'), a single band, or 'none' for a blank "
                        'canvas. You then add (green) / erase (red) and choose to apply, '
                        'reject or skip on closing the GUI')
    for name, val in detect_arcs.DEFAULTS.items():
        p.add_argument(f'--detect-{name.replace("_", "-")}', dest=f'detect_{name}',
                       type=type(val), default=val,
                       help=f'detector parameter passed through to detect_arcs (default {val})')
    p.add_argument('--include-ignored', action='store_true', default=False,
                   help='also offer bands whose cutout directory is GITIGNORED. Off by '
                        'default: such a band is not a science product (the three 420 s '
                        'SLACS SNAP f814W diagnostics are why the rule exists), and hand-drawn '
                        'work written there is non-regenerable and invisible to every clone. '
                        'Use only to work on a diagnostic deliberately.')
    a = p.parse_args()

    # Same rule as make_masks.py: --broadcast writes one draw to every band unseen, which is
    # the opposite of reviewing a proposal per band. Combining them would push a region
    # reviewed for one band over another band's own.
    if a.broadcast and a.propose_from != 'none':
        p.error('--broadcast is the unreviewed route and is mutually exclusive with '
                '--propose-from; pass --propose-from none alongside it')
    detect_kw = {k: getattr(a, f'detect_{k}') for k in detect_arcs.DEFAULTS}

    root = cutout_paths.cutouts_root(ws_path, a.size)

    lens_filts = {}
    for lens, filt, cutout_dir in make_masks.discover_targets(
            root, a.sample, a.lens, None, include_ignored=a.include_ignored):
        lens_filts.setdefault(lens, {})[filt] = cutout_dir

    if not lens_filts:
        raise SystemExit(f"no cutouts found under {root} for sample {a.sample} matching "
                         f"lens={a.lens!r}")

    print(f"{len(lens_filts)} lens(es) to process")
    made = skipped = 0
    for lens in sorted(lens_filts):
        if process_lens_arc_mask(lens, lens_filts[lens], a.sample, a.filt, a.drizzle_pass,
                                 a.force, a.broadcast, a.brush_radius, a.brush_width,
                                 a.display, a.stretch, a.vmin_percent, a.vmax_percent,
                                 a.asinh_a, a.subtract_radial, a.side_by_side,
                                 a.show_positions, a.cross_size,
                                 propose_from=a.propose_from, detect_kw=detect_kw,
                                 show_contaminants=a.show_contaminants):
            made += 1
        else:
            skipped += 1

    print(f"\nDone: {made} lens(es) arc-masked, {skipped} skipped")


if __name__ == '__main__':
    main()
