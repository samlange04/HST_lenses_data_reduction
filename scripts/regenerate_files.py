#!/usr/bin/env python
"""
Regenerate everything git does NOT carry, for one sample or all of them.

A fresh clone has the scripts, the tracking JSONs in info/, and the finished science
tree (data/cutouts/: stamps, noise maps, hand-drawn masks, positions, PSFs) plus the
QC mosaics in data/mosaics/. Everything else under data/ is gitignored because it is
either large or a plain re-render of a tracked file. This wrapper rebuilds those
ignored products from the tracked ones, in dependency order, by calling the existing
runners -- it adds no reduction logic of its own, so the traps the runners exist to
avoid (WFPC2 align step, split visits, one-instrument-per-band) stay handled there.

Stages, in the order they run:

  cutout-pngs     data/cutouts/**/cutout[_cr].png        LIGHT: redrawn from the tracked
                  sci/noise stamps (make_cutouts.plot_cutouts). No mosaic needed.
  positions-pngs  data/cutouts/**/*_positions.png        LIGHT: make_positions.py
                  --overlays-only, redrawn from the tracked positions JSONs. No GUI.
  drizzle         data/calibrated/ data/reference_files/  HEAVY: MAST download + CRDS +
                  data/drizzle_files/ data/drizzled/      AstroDrizzle via the per-instrument
                  data/run_logs/                           runners (run_acs_all.sh,
                  run_wfc3_all.sh, run_wfpc2_wf3.sh, run_gallery_uvis_all.sh). Only the
                  instruments the sample has products for in info/lens_instrument.json are
                  run (--all-instruments overrides). Tens of GB and hours per sample.
                  NOTE run_wfpc2_wf3.sh also re-cuts the F606W stamps (drizzle -> align ->
                  cutout is one unit there), so tracked F606W stamps may show as modified.
  cutouts         data/cutouts_<S>arcsec/ data/mosaics_<S>arcsec/   HEAVY, only with --size:
                  info/lens_cutout_qc_<S>arcsec.json      stamps at a non-default size
                  (run_cutouts_all.sh <sample> <S> + make_mosaics.py --size <S>). Needs
                  the drizzle stage's mosaics. The 12" tree is tracked and is NOT re-cut.
  psf             data/psf/<sample>/                      HEAVY, OPT-IN (--stages psf):
                  run_psf_all.sh. Needs the mosaics, and it also rewrites the tracked
                  cutout_[cr_]psf*.fits kernels, so it is never part of a default run.

Default stage set: cutout-pngs, positions-pngs, drizzle (plus cutouts when --size is
given). --light runs only the two LIGHT stages (minutes, no network). --stages picks any
subset explicitly. --dry-run prints every command without running anything.

Usage:
    uv run python scripts/regenerate_files.py --sample slacs_gold --light
    uv run python scripts/regenerate_files.py --sample gallery
    uv run python scripts/regenerate_files.py --sample all --size 20
    uv run python scripts/regenerate_files.py --sample bells --stages drizzle psf
    uv run python scripts/regenerate_files.py --sample all --dry-run

Every stage is idempotent (the drizzle scripts skip existing products; the renders just
overwrite), so an interrupted run is resumed by running the same command again. A stage
that fails is reported and the run continues with the remaining stages and samples; the
exit status is non-zero if anything failed.
"""

import argparse
import glob
import os
import subprocess
import sys

import matplotlib
matplotlib.use('Agg')
import numpy as np
from astropy.io import fits

ws_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
scripts_dir = os.path.join(ws_path, 'scripts')
sys.path.insert(0, scripts_dir)
import cutout_paths
import info_json
import mast_target_names
from make_cutouts import plot_cutouts

LIGHT_STAGES = ['cutout-pngs', 'positions-pngs']
HEAVY_STAGES = ['drizzle', 'cutouts', 'psf']
ALL_STAGES = LIGHT_STAGES + HEAVY_STAGES
DEFAULT_STAGES = ['cutout-pngs', 'positions-pngs', 'drizzle']   # + cutouts when --size

# instrument string recorded in info/lens_instrument.json -> the runner that makes it
RUNNERS = [
    ('ACS/WFC',   'run_acs_all.sh'),
    ('WFC3/IR',   'run_wfc3_all.sh'),
    ('WFPC2/WF3', 'run_wfpc2_wf3.sh'),
    ('WFC3/UVIS', 'run_gallery_uvis_all.sh'),
]


def run(cmd, dry_run, cwd=ws_path):
    """Run a command with its output streamed to the console. Returns True on success."""
    print(f"\n$ {' '.join(cmd)}", flush=True)
    if dry_run:
        return True
    return subprocess.run(cmd, cwd=cwd).returncode == 0


def sample_instruments(sample):
    """Instruments with at least one product recorded for the sample, from
    info/lens_instrument.json. Empty when the sample has never been reduced."""
    recs = info_json.load(os.path.join(ws_path, 'info', 'lens_instrument.json')).get(sample, {})
    return {v for lens in recs.values() for v in (lens or {}).values() if v}


# ── LIGHT stages ─────────────────────────────────────────────────────────────────────

def stage_cutout_pngs(sample, dry_run):
    """Redraw the 3-panel cutout[_cr][_<variant>].png beside every tracked sci/noise stamp
    in the default-size tree. Same renderer and the same title fields as make_cutouts.py
    (drizzle pass, size and recentring offset come from info/lens_cutout_qc.json; the size
    falls back to the stamp's CUTSIZE card)."""
    root = os.path.join(cutout_paths.cutouts_root(ws_path), sample)
    if not os.path.isdir(root):
        print(f"  no cutouts under {root}; nothing to render")
        return True
    qc_path = cutout_paths.qc_json_path(ws_path)
    n = bad = 0
    for cutout_dir in sorted(glob.glob(os.path.join(root, '*', '*'))):
        if not os.path.isdir(cutout_dir):
            continue
        lens = os.path.basename(os.path.dirname(cutout_dir))
        filt = os.path.basename(cutout_dir)
        for prefix in cutout_paths.PREFIXES:
            for variant, sci_path in sorted(cutout_paths.list_stamps(cutout_dir, prefix, 'sci').items()):
                noise_path = os.path.join(cutout_dir, cutout_paths.stamp_name(prefix, 'noise', variant))
                if not os.path.exists(noise_path):
                    print(f"  {lens} {filt}: {os.path.basename(sci_path)} has no noise stamp; skipped")
                    bad += 1
                    continue
                png = os.path.join(cutout_dir, cutout_paths.stamp_png_name(prefix, variant))
                if dry_run:
                    print(f"  would write {os.path.relpath(png, ws_path)}")
                    n += 1
                    continue
                sci_hdr = fits.getheader(sci_path)
                qc = info_json.entry(qc_path, sample, lens, filt, default={}) or {}
                drizzle_pass = qc.get('drizzle_pass') or ('cr' if prefix == 'cutout_cr' else 'nocrrej')
                size = qc.get('size_arcsec') or sci_hdr.get('CUTSIZE', cutout_paths.DEFAULT_SIZE)
                offset = qc.get('offset_arcsec')
                title = f"{lens}  {filt}  [{drizzle_pass}]  ({size:g}\" cutout"
                title += f", recentred {offset:.2f}\" from catalogue)" if offset is not None else ")"
                plot_cutouts(fits.getdata(sci_path).astype(float),
                             fits.getdata(noise_path).astype(float), png, title)
                n += 1
    print(f"  {n} cutout PNG(s) {'planned' if dry_run else 'written'}, {bad} skipped")
    return bad == 0


def stage_positions_pngs(sample, dry_run):
    root = os.path.join(cutout_paths.cutouts_root(ws_path), sample)
    if not os.path.isdir(root):
        print(f"  no cutouts under {root}; nothing to render")
        return True
    return run([sys.executable, os.path.join(scripts_dir, 'make_positions.py'),
                '--sample', sample, '--overlays-only'], dry_run)


# ── HEAVY stages ─────────────────────────────────────────────────────────────────────

def stage_drizzle(sample, dry_run, all_instruments=False):
    have = sample_instruments(sample)
    if all_instruments or not have:
        if not have and not all_instruments:
            print(f"  {sample} has no products recorded in info/lens_instrument.json; "
                  f"trying every instrument runner (each lens costs one MAST query per band)")
        wanted = [r for _, r in RUNNERS]
    else:
        wanted = [r for inst, r in RUNNERS if inst in have]
        skipped = [inst for inst, _ in RUNNERS if inst not in have]
        if skipped:
            print(f"  skipping runners for instruments {sample} has no products in: "
                  f"{', '.join(skipped)} (--all-instruments to try them anyway)")
    ok = True
    for runner in wanted:
        ok &= run(['bash', os.path.join(scripts_dir, runner), sample], dry_run)
    return ok


def stage_cutouts(sample, dry_run, sizes):
    if not sizes:
        print("  no --size given; the default 12\" stamps are tracked and are not re-cut")
        return True
    ok = True
    for size in sizes:
        if float(size) == cutout_paths.DEFAULT_SIZE:
            print(f"  --size {size:g} is the tracked default tree; not re-cut")
            continue
        ok &= run(['bash', os.path.join(scripts_dir, 'run_cutouts_all.sh'), sample, f'{size:g}'], dry_run)
        ok &= run([sys.executable, os.path.join(scripts_dir, 'make_mosaics.py'),
                   '--sample', sample, '--size', f'{size:g}'], dry_run)
    return ok


def stage_psf(sample, dry_run):
    return run(['bash', os.path.join(scripts_dir, 'run_psf_all.sh'), sample], dry_run)


# ── driver ────────────────────────────────────────────────────────────────────────────

def main():
    samples = mast_target_names.samples()
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--sample', required=True, choices=list(samples) + ['all'],
                   help='a sample from info/lens_samples.json, or "all" for every sample')
    p.add_argument('--stages', nargs='+', choices=ALL_STAGES, default=None, metavar='STAGE',
                   help=f'stages to run, in pipeline order whatever order given here '
                        f'(choices: {", ".join(ALL_STAGES)}; default: '
                        f'{", ".join(DEFAULT_STAGES)}, plus cutouts when --size is given)')
    p.add_argument('--light', action='store_true', default=False,
                   help=f'only the stages that need nothing but the tracked files: '
                        f'{", ".join(LIGHT_STAGES)}')
    p.add_argument('--size', type=float, nargs='+', default=[], metavar='ARCSEC',
                   help='also cut and mosaic size-variant stamps (data/cutouts_<S>arcsec/) '
                        'at these sizes; needs the drizzle stage\'s mosaics')
    p.add_argument('--all-instruments', action='store_true', default=False,
                   help='drizzle: run every instrument runner, not only the ones the '
                        'sample has products for in info/lens_instrument.json')
    p.add_argument('--dry-run', action='store_true', default=False,
                   help='print every command and output path; run nothing')
    a = p.parse_args()

    if a.light and a.stages:
        p.error('--light and --stages are alternatives')
    if a.light:
        stages = list(LIGHT_STAGES)
    elif a.stages:
        stages = [s for s in ALL_STAGES if s in a.stages]
    else:
        stages = list(DEFAULT_STAGES) + (['cutouts'] if a.size else [])
    if 'cutouts' in stages and not a.size and a.stages:
        p.error('the cutouts stage needs --size <arcsec> (the 12" tree is tracked)')

    todo = list(samples) if a.sample == 'all' else [a.sample]
    print(f"regenerate_files: samples {', '.join(todo)}; stages {', '.join(stages)}"
          f"{'  [DRY RUN]' if a.dry_run else ''}")
    if any(s in HEAVY_STAGES for s in stages) and not a.dry_run:
        print("  heavy stages selected: MAST downloads and AstroDrizzle take hours and tens "
              "of GB per sample. Interrupt freely -- re-running resumes.")

    failed = []
    for sample in todo:
        for stage in stages:
            print(f"\n===== {sample}: {stage} =====", flush=True)
            if stage == 'cutout-pngs':
                ok = stage_cutout_pngs(sample, a.dry_run)
            elif stage == 'positions-pngs':
                ok = stage_positions_pngs(sample, a.dry_run)
            elif stage == 'drizzle':
                ok = stage_drizzle(sample, a.dry_run, a.all_instruments)
            elif stage == 'cutouts':
                ok = stage_cutouts(sample, a.dry_run, a.size)
            elif stage == 'psf':
                ok = stage_psf(sample, a.dry_run)
            if not ok:
                failed.append(f'{sample}:{stage}')
                print(f"===== {sample}: {stage} FAILED (continuing) =====", flush=True)

    print("\n===== regenerate_files done =====")
    if failed:
        print(f"FAILED stages ({len(failed)}): {', '.join(failed)} -- see data/run_logs/ for "
              f"the per-lens logs the runners wrote")
        sys.exit(1)
    print("all stages succeeded")


if __name__ == '__main__':
    main()
