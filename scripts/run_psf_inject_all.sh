#!/usr/bin/env bash
# Build a drizzle-broadened, injection-drizzled model PSF (make_psf_inject.py, the Anderson
# 2016 route) for the MODEL-tier products -- the only ones it improves. An empirical ePSF is
# already the true drizzled PSF (cut from the drizzled mosaic), so injection is a no-op there;
# by default this runner skips empirical products and only (re)builds the model tier.
#
# NOTE: make_psf.py now auto-chains into make_psf_inject.py (promote=True) right after
# building a model-tier product, so run_psf_all.sh already leaves every model-tier lens
# promoted -- this runner is rarely needed for routine use. It's still useful to: re-run
# injection/promotion alone after an injection-code change (without re-running the analytic
# model build), or run --all for the empirical-vs-injected validation comparison.
#
# Usage: run_psf_inject_all.sh [SAMPLE] [--all|--models-only]
#   (default SAMPLE: mast_target_names.DEFAULT_SAMPLE;  default mode: --models-only)
#
#   --models-only  (default) run only products whose lens_psf.json method starts with `model`
#                  or `inject` (STDPSF / focus-diverse / MAST PSF DB, promoted or not) -- the
#                  tier that lacks drizzle broadening on its analytic build. This PROMOTES:
#                  the injected kernel becomes the canonical psf_kernel.fits /
#                  cutout_[cr_]psf.fits, and the current canonical (if not already an
#                  injected build) is moved aside to *_analytic.
#   --all          also run empirical products, so psf_kernel_injected.fits can be compared
#                  against the empirical truth (validation: injected-model should approach the
#                  empirical FWHM, both being drizzled -- see the F160W check in CLAUDE.md).
#                  Promotion never applies to an empirical primary -- these write the
#                  parallel *_injected-suffixed comparison files instead, as before.
#
# Mirrors run_psf_all.sh: globs data/drizzled/ (a PSF needs a mosaic), reads the target set
# from info/lens_psf.json, writes per-product logs, and reports ok/failed/skipped. Requires
# the persisted data/drizzle_files/ inputs (present after a normal drizzle run). bash 3.2.
SD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; WS="$(dirname "$SD")"
LOG="$WS/data/run_logs"; mkdir -p "$LOG"

MODE="models"; SAMPLE_ARG=""
for arg in "$@"; do
  case "$arg" in
    --models-only) MODE="models" ;;
    --all)         MODE="all" ;;
    -*)            echo "unknown flag: $arg (use --all or --models-only)" >&2; exit 2 ;;
    *)             SAMPLE_ARG="$arg" ;;
  esac
done
SAMPLE="$(uv run --project "$WS" python "$SD/mast_target_names.py" ${SAMPLE_ARG:+"$SAMPLE_ARG"} --print-sample)" || exit 1

# In --models-only mode, the RUN set is every product recorded as model-tier in
# lens_psf.json -- 'model...' (not yet promoted) or 'inject...' (already promoted by a
# previous run; still worth re-running after an injection-code change). In --all mode we
# iterate the mosaics directly (below) and let make_psf_inject decide.
MODEL_RUN=""
if [ "$MODE" = models ]; then
  MODEL_RUN="$(uv run --project "$WS" python -c "import json,os; p=os.path.join('$WS','info','lens_psf.json'); d=json.load(open(p)) if os.path.exists(p) else {}; s=d.get('$SAMPLE', {}); print(chr(10).join(l+'/'+k for l in s for k,v in (s[l] or {}).items() if v and (str(v.get('method','')).startswith('model') or str(v.get('method','')).startswith('inject'))))")" || exit 1
fi

ok=0; fail=0; skip=0
# f438W is gallery's blue UVIS science band (mirrors run_psf_all.sh's loop); F225W/F275W are
# omitted -- unusable for lens science sample-wide, so no PSF/injection work is done on them.
# Extra filters are no-ops for samples that lack them (the glob matches nothing).
for filt in f606W f814W f555W f160W f438W; do
  for dd in "$WS"/data/drizzled/"$SAMPLE"/*/"$filt"*; do
    [ -d "$dd" ] || continue
    ls "$dd"/*_sci.fits >/dev/null 2>&1 || continue
    lens=$(basename "$(dirname "$dd")")
    key=$(basename "$dd")
    if [ "$MODE" = models ]; then
      case "$(printf '\n%s\n' "$MODEL_RUN")" in
        *"$(printf '\n%s\n' "$lens/$key")"*) : ;;   # in the model set -> run
        *) skip=$((skip+1)); echo "skip non-model: $lens $key"; continue ;;
      esac
    fi
    if uv run --project "$WS" python "$SD/make_psf_inject.py" --lens "$lens" --filt "$key" \
         --sample "$SAMPLE" > "$LOG/${lens}_${key}_psf_inject.log" 2>&1; then
      ok=$((ok+1)); echo "ok: $lens $key"
    else
      fail=$((fail+1)); echo "FAILED: $lens $key (see $LOG/${lens}_${key}_psf_inject.log)"
    fi
  done
done
echo "=== psf_inject $SAMPLE ($MODE) done: $ok ok, $fail failed, $skip skipped  $(date +%H:%M:%S) ==="
