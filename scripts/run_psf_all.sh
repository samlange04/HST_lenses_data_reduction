#!/usr/bin/env bash
# Build a PSF for every drizzled product. make_psf.py defaults to --pass auto (detect
# stars in the LACosmic CR mosaic where one exists, else the no-CR mosaic) and --method
# auto (empirical ePSF, model fallback when a field is too star-poor or its empirical
# wings are poor). For a model-tier build, make_psf.py auto-chains make_psf_inject.py
# (promote=True) so the canonical psf_kernel.fits / cutout_[cr_]psf.fits it leaves behind
# is already the drizzle-broadened injected kernel, not the sharper analytic model (which
# is kept alongside as *_analytic for comparison) -- see make_psf_inject.py's docstring.
#
# Usage: run_psf_all.sh [SAMPLE] [--all|--models-only]
#   (default SAMPLE: mast_target_names.DEFAULT_SAMPLE;  default mode: --all)
#
#   --all          (default) (re)build every product -- empirical AND model.
#   --models-only  skip any product whose EXISTING lens_psf.json record is `empirical`;
#                  (re)build only the model tier (STDPSF / focus-diverse / MAST PSF DB) plus
#                  any product not yet built. Use this when you change only model-PSF code
#                  (e.g. the detector->North-up rotation) and don't want to churn the good
#                  empirical builds, which that code never touches (they are cut from the
#                  North-up mosaic itself; the pedestal subtraction is a ~1e-4 no-op on the
#                  sharp ACS bands). The skip set is "empirical AND already recorded", so a
#                  lens with no entry still builds.
#
# Like run_cutouts_all.sh -- and unlike the drizzle runners -- this globs data/drizzled/
# rather than the roster: a PSF can only be built from a mosaic that exists, so the
# products on disk are the right thing to iterate. Written for macOS bash 3.2.
SD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; WS="$(dirname "$SD")"
LOG="$WS/data/run_logs"; mkdir -p "$LOG"

MODE="all"; SAMPLE_ARG=""
for arg in "$@"; do
  case "$arg" in
    --models-only) MODE="models" ;;
    --all)         MODE="all" ;;
    -*)            echo "unknown flag: $arg (use --all or --models-only)" >&2; exit 2 ;;
    *)             SAMPLE_ARG="$arg" ;;
  esac
done
SAMPLE="$(uv run --project "$WS" python "$SD/mast_target_names.py" ${SAMPLE_ARG:+"$SAMPLE_ARG"} --print-sample)" || exit 1

# --models-only skip set: "<lens>/<key>" for every product recorded as empirical in
# info/lens_psf.json (an empirical record implies a built product). Newline-delimited so the
# per-product test below is a shell `case` match -- no python/grep spawned inside the loop.
EMP_SKIP=""
if [ "$MODE" = models ]; then
  EMP_SKIP="$(uv run --project "$WS" python -c "import json,os; p=os.path.join('$WS','info','lens_psf.json'); d=json.load(open(p)) if os.path.exists(p) else {}; s=d.get('$SAMPLE', {}); print(chr(10).join(l+'/'+k for l in s for k,v in (s[l] or {}).items() if v and str(v.get('method','')).startswith('empirical')))")" || exit 1
fi

# Globs are "$filt"* , not "$filt", so per-visit dirs (f606W_v1, f606W_v2) are included;
# --filt is taken from the directory basename. Any *_sci.fits means the mosaic is buildable
# (ACS/WFPC2 hold *_cr_*; F160W holds *nocrrej* -- make_psf --pass auto picks the right one).
# f438W is gallery's blue UVIS science band; F225W/F275W are deliberately omitted -- they
# are confirmed unusable for lens science sample-wide (arc undetected), so no PSF work is
# done on them (see CLAUDE.md, gallery_uv_bands_unusable). The extra filters are no-ops for
# samples that lack them (the glob matches nothing).
ok=0; fail=0; skip=0
for filt in f606W f814W f555W f160W f438W; do
  for d in "$WS"/data/drizzled/"$SAMPLE"/*/"$filt"*; do
    [ -d "$d" ] || continue
    ls "$d"/*_sci.fits >/dev/null 2>&1 || continue
    lens=$(basename "$(dirname "$d")")
    key=$(basename "$d")
    if [ "$MODE" = models ]; then
      case "$(printf '\n%s\n' "$EMP_SKIP")" in
        *"$(printf '\n%s\n' "$lens/$key")"*)
          skip=$((skip+1)); echo "skip empirical (exists): $lens $key"; continue ;;
      esac
    fi
    if uv run --project "$WS" python "$SD/make_psf.py" --lens "$lens" --filt "$key" \
         --sample "$SAMPLE" > "$LOG/${lens}_${key}_psf.log" 2>&1; then
      ok=$((ok+1))
    else
      fail=$((fail+1)); echo "FAILED: $lens $key (see $LOG/${lens}_${key}_psf.log)"
    fi
  done
done
echo "=== psf $SAMPLE ($MODE) done: $ok ok, $fail failed, $skip skipped  $(date +%H:%M:%S) ==="
