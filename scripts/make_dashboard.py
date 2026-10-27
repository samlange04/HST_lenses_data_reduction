#!/usr/bin/env python
"""Build the status dashboard: one self-contained HTML page of every sample, lens and band,
with thumbnails, completeness tables and the per-band provenance from the tracking JSONs.

    uv run python scripts/make_dashboard.py                 # all samples -> diagnostics/dashboard/
    uv run python scripts/make_dashboard.py --sample bells  # one sample only
    uv run python scripts/make_dashboard.py --no-thumbs     # rewrite index.html, keep images

Reads ONLY tracked inputs (info/*.json and the data/stamps stamps, resolved through
stamp_paths so no stamp name is ever spelled here) plus data/mosaics PNGs. Writes
diagnostics/dashboard/{index.html, thumbs/<sample>/<lens>/<band>.jpg,
mosaics/<sample>/<group>_<kind>.jpg}. The page opens from disk (relative image paths; the
lightbox then loads the full `_dataset.png` out of data/stamps) and is what gets published
as the shareable artifact (where the lightbox falls back to the thumbnail).

It is a QC/overview product, not a science product: nothing here feeds a fit.
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

ws_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ws_path, 'scripts'))
import stamp_paths  # noqa: E402
import mast_target_names  # noqa: E402
import info_json  # noqa: E402

THUMB_PX = 260
MOSAIC_W = 1600

# Display order for bands inside a lens card; anything unknown sorts after, alphabetically.
BAND_ORDER = ['f814W', 'f606W', 'f606W_v2', 'f555W', 'f475X', 'f600LP', 'f438W', 'f160W']
# lens_instrument.json label -> noise_corr_factors.json key
NOISE_INSTR = {'ACS/WFC': 'ACS/WFC', 'WFPC2/WF3': 'WFPC2', 'WFPC2/PC': 'WFPC2',
               'WFC3/IR': 'WFC3/IR', 'WFC3/UVIS': 'WFC3/UVIS'}


def load_json(name):
    path = info_json.path(ws_path, name)
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return json.load(f)


def band_key(b):
    return (BAND_ORDER.index(b) if b in BAND_ORDER else len(BAND_ORDER), b)


def base_filt(b):
    """'f606W_v2' -> 'f606W' (the split-visit suffix is not a filter)."""
    return b.split('_')[0]


# ───────────────────────────── data collection ──────────────────────────────

def collect(samples):
    inst_j = load_json('lens_instrument.json')
    prod_j = load_json('lens_products.json')
    masks_j = load_json('lens_masks.json')
    pos_j = load_json('lens_positions.json')
    psf_j = load_json('lens_psf.json')
    ein_j = load_json('lens_einstein_radii.json')
    noise_j = load_json('noise_corr_factors.json')
    cr_j = load_json('lens_cr_params.json')
    exp_j = load_json('lens_exptime.json')
    samples_j = load_json('lens_samples.json')

    overrides = noise_j.get('_lens_overrides', {})
    lenses_out = []
    for sample in samples:
        croot = stamp_paths.stamps_root(ws_path)
        for lens in mast_target_names.lenses(sample):
            ldir = os.path.join(croot, sample, lens)
            if not os.path.isdir(ldir):
                continue
            bands = []
            for band in sorted(os.listdir(ldir), key=band_key):
                bdir = os.path.join(ldir, band)
                if not os.path.isdir(bdir):
                    continue
                prefix = stamp_paths.find_prefix(bdir)
                if prefix is None:
                    continue
                sci = stamp_paths.find_stamp(bdir, prefix, 'sci')
                noise = stamp_paths.find_stamp(bdir, prefix, 'noise')
                variant = stamp_paths.stamp_variant(sci)
                have = {
                    'noise': noise is not None,
                    'psf': os.path.exists(os.path.join(bdir, f'psf.fits')),
                    'psf_err': os.path.exists(os.path.join(bdir, f'psf_err.fits')),
                    'mask': os.path.exists(os.path.join(bdir, f'{prefix}_mask.fits')),
                    'positions_file': os.path.exists(os.path.join(bdir, f'{prefix}_positions.json')),
                    'dataset_png': os.path.exists(os.path.join(bdir, f'{prefix}_dataset.png')),
                }
                instr = (inst_j.get(sample, {}).get(lens, {}) or {}).get(base_filt(band)) \
                    or (inst_j.get(sample, {}).get(lens, {}) or {}).get(band)
                frames = (prod_j.get(sample, {}).get(lens, {}) or {}).get(band) \
                    or (prod_j.get(sample, {}).get(lens, {}) or {}).get(base_filt(band)) or []
                exptime = (exp_j.get(sample, {}).get(lens, {}) or {}).get(band) \
                    or (exp_j.get(sample, {}).get(lens, {}) or {}).get(base_filt(band))
                psf = (psf_j.get(sample, {}).get(lens, {}) or {}).get(band) or {}
                mask = (masks_j.get(sample, {}).get(lens, {}) or {}).get(band) or {}
                nf_key = NOISE_INSTR.get(instr or '', None)
                nf = (noise_j.get(nf_key, {}) or {}).get(base_filt(band)) if nf_key else None
                nf_over = (overrides.get(sample, {}).get(lens, {}) or {}).get(base_filt(band))
                cr = (cr_j.get(sample, {}).get(lens, {}) or {}).get(base_filt(band))
                bands.append({
                    'band': band, 'prefix': prefix, 'variant': variant,
                    'cr_pass': stamp_paths.stamp_pass(sci) == 'cr',
                    'instrument': instr, 'n_frames': len(frames), 'exptime': exptime,
                    'have': have,
                    'psf_method': psf.get('method'), 'psf_n_stars': psf.get('n_stars'),
                    'psf_fwhm_pix': psf.get('fwhm_pix'), 'psf_err_frac': psf.get('psf_err_frac'),
                    'mask_source': mask.get('source'), 'mask_px': mask.get('n_excluded_px'),
                    'mask_note': mask.get('note'),
                    'noise_factor': nf_over if nf_over is not None else nf,
                    'noise_override': nf_over is not None,
                    'cr_params': cr,
                    'dataset_rel': os.path.relpath(os.path.join(bdir, f'{prefix}_dataset.png'), ws_path)
                    if have['dataset_png'] else None,
                    'thumb': f'thumbs/{sample}/{lens}/{band}.jpg',
                    'sci_path': sci,
                    'mask_path': os.path.join(bdir, f'{prefix}_mask.fits') if have['mask'] else None,
                })
            # positions are a lens-level decision (marked in one band, broadcast to the rest)
            pos_entry = pos_j.get(sample, {}).get(lens, {}) or {}
            pos_status, pos_n, pos_band = 'none', 0, None
            for b, p in pos_entry.items():
                if not isinstance(p, dict):
                    continue
                if p.get('n_positions', 0) > 0:
                    pos_status, pos_n, pos_band = 'marked', p['n_positions'], p.get('marked_filt', b)
                    break
                if p.get('deliberately_empty'):
                    pos_status, pos_band = 'empty', p.get('marked_filt', b)
            ein = ein_j.get(lens, {}) or {}
            lenses_out.append({
                'sample': sample, 'lens': lens, 'bands': bands,
                'theta_e': ein.get('theta_e_arcsec'), 'zlens': ein.get('zlens'),
                'zsrc': ein.get('zsrc'), 'sigma': ein.get('sigma_kms'), 'ein_source': ein.get('source'),
                'positions': pos_status, 'n_positions': pos_n, 'positions_band': pos_band,
            })
    meta = {
        'samples': {s: {'description': samples_j[s]['description'],
                        'n_roster': len(samples_j[s]['lenses'])} for s in samples},
        'noise_factors': {k: v for k, v in noise_j.items() if not k.startswith('_')},
        'noise_overrides': overrides,
        'noise_method': noise_j.get('_method', ''),
        'cr_params': cr_j,
    }
    return lenses_out, meta


# ───────────────────────────── images ──────────────────────────────

def render_thumb(sci_path, mask_path, out_path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from astropy.io import fits
    from astropy.visualization import AsinhStretch, ImageNormalize

    data = fits.getdata(sci_path).astype(float)
    finite = data[np.isfinite(data)]
    vmin, vmax = np.percentile(finite, 2.0), np.percentile(finite, 99.8)
    norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=AsinhStretch(0.08))
    fig = plt.figure(figsize=(THUMB_PX / 100, THUMB_PX / 100), dpi=100)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(data, norm=norm, origin='lower', cmap='inferno', interpolation='nearest')
    if mask_path and os.path.exists(mask_path):
        m = fits.getdata(mask_path)
        if np.any(m > 0):
            # black-stroked bright outline: plain cyan vanishes on inferno
            ax.contour(m > 0, levels=[0.5], colors='black', linewidths=2.2)
            ax.contour(m > 0, levels=[0.5], colors='#5ff2ff', linewidths=0.9)
    ax.set_axis_off()
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, dpi=100, pil_kwargs={'quality': 82})
    plt.close(fig)


def downsize_mosaics(samples, out_dir, bands_on_disk):
    """Mosaics whose every band still has stamps in the sample (gallery's F225W/F275W mosaics
    outlived the deleted UV products) -> JPEG at MOSAIC_W."""
    from PIL import Image
    found = {}
    for sample in samples:
        mdir = os.path.join(stamp_paths.mosaics_root(ws_path), sample)
        for png in sorted(glob.glob(os.path.join(mdir, '*.png'))):
            base = os.path.splitext(os.path.basename(png))[0]
            if not base.endswith(('_signal', '_snr', '_psf', '_noise')):
                continue
            group = base.rsplit('_', 1)[0]
            if not all(b in bands_on_disk.get(sample, set()) for b in group.split('_')):
                continue
            rel = f'mosaics/{sample}/{base}.jpg'
            dst = os.path.join(out_dir, rel)
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            if not os.path.exists(dst) or os.path.getmtime(dst) < os.path.getmtime(png):
                im = Image.open(png).convert('RGB')
                if im.width > MOSAIC_W:
                    im = im.resize((MOSAIC_W, round(im.height * MOSAIC_W / im.width)), Image.LANCZOS)
                im.save(dst, quality=78, optimize=True)
            found.setdefault(sample, []).append({'name': base, 'path': rel})
    return found


# ───────────────────────────── page ──────────────────────────────

def build_html(lenses, meta, mosaics, generated):
    payload = json.dumps({'lenses': [{k: v for k, v in L.items()} for L in lenses],
                          'meta': meta, 'mosaics': mosaics, 'generated': generated},
                         default=str)
    # strip absolute disk paths from the payload: the page needs none of them
    for L in lenses:
        for b in L['bands']:
            b.pop('sci_path', None)
            b.pop('mask_path', None)
    payload = json.dumps({'lenses': lenses, 'meta': meta, 'mosaics': mosaics,
                          'generated': generated}, default=str)
    template_path = os.path.join(ws_path, 'scripts', 'dashboard_template.html')
    with open(template_path) as f:
        html = f.read()
    return html.replace('/*__DATA__*/null', payload)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--sample', default='all', help='one sample name or "all" (default)')
    p.add_argument('--out-dir', default=os.path.join(ws_path, 'diagnostics', 'dashboard'))
    p.add_argument('--no-thumbs', action='store_true', help='rewrite index.html only')
    p.add_argument('--force', action='store_true', help='re-render thumbnails that already exist')
    a = p.parse_args()

    samples = mast_target_names.samples() if a.sample == 'all' else [a.sample]
    lenses, meta = collect(samples)
    os.makedirs(a.out_dir, exist_ok=True)

    n_thumb = 0
    if not a.no_thumbs:
        for L in lenses:
            for b in L['bands']:
                dst = os.path.join(a.out_dir, b['thumb'])
                if a.force or not os.path.exists(dst) or os.path.getmtime(dst) < os.path.getmtime(b['sci_path']):
                    render_thumb(b['sci_path'], b['mask_path'], dst)
                    n_thumb += 1
        print(f'thumbnails rendered: {n_thumb}')
    bands_on_disk = {}
    for L in lenses:
        bands_on_disk.setdefault(L['sample'], set()).update(base_filt(b['band']) for b in L['bands'])
    mosaics = downsize_mosaics(samples, a.out_dir, bands_on_disk)

    from datetime import date
    html = build_html(lenses, meta, mosaics, date.today().isoformat())
    with open(os.path.join(a.out_dir, 'index.html'), 'w') as f:
        f.write(html)
    n_bands = sum(len(L['bands']) for L in lenses)
    print(f'{len(lenses)} lenses, {n_bands} band products, {len(samples)} samples -> {a.out_dir}/index.html')


if __name__ == '__main__':
    main()
