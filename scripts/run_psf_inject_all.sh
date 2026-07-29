#!/usr/bin/env bash
# Build a drizzle-broadened, injection-drizzled model PSF (make_psf_inject.py, the Anderson
# 2016 route) for the MODEL-tier products -- the only ones it improves. An empirical ePSF is
# already the true drizzled PSF (cut from the drizzled mosaic), so injection is a no-op there;
# by default this runner skips empirical products and only (re)builds the model tier.
#
# Usage: run_psf_inject_all.sh [SAMPLE] [--all|--models-only]
#   (default SAMPLE: mast_target_names.DEFAULT_SAMPLE;  default mode: --models-only)
#
#   --models-only  (default) run only products whose lens_psf.json method starts with `model`
#                  (STDPSF / focus-diverse / MAST PSF DB) -- the tier that lacks drizzle
#                  broadening. This is the set injection is for.
#   --all          also run empirical products, so psf_kernel_injected.fits can be compared
#                  against the empirical truth (validation: injected-model should approach the
#                  empirical FWHM, both being drizzled -- see the F160W check in CLAUDE.md).
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
SAMPLE="$(conda run -n stenv python "$SD/mast_target_names.py" ${SAMPLE_ARG:+"$SAMPLE_ARG"} --print-sample)" || exit 1

# In --models-only mode, the RUN set is every product recorded as a model in lens_psf.json.
# In --all mode we iterate the mosaics directly (below) and let make_psf_inject decide.
MODEL_RUN=""
if [ "$MODE" = models ]; then
  MODEL_RUN="$(conda run -n stenv python -c "import json,os; p=os.path.join('$WS','info','lens_psf.json'); d=json.load(open(p)) if os.path.exists(p) else {}; s=d.get('$SAMPLE', {}); print(chr(10).join(l+'/'+k for l in s for k,v in (s[l] or {}).items() if v and str(v.get('method','')).startswith('model')))")" || exit 1
fi

ok=0; fail=0; skip=0
for filt in f606W f814W f555W f160W; do
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
    if conda run -n stenv python "$SD/make_psf_inject.py" --lens "$lens" --filt "$key" \
         --sample "$SAMPLE" > "$LOG/${lens}_${key}_psf_inject.log" 2>&1; then
      ok=$((ok+1)); echo "ok: $lens $key"
    else
      fail=$((fail+1)); echo "FAILED: $lens $key (see $LOG/${lens}_${key}_psf_inject.log)"
    fi
  done
done
echo "=== psf_inject $SAMPLE ($MODE) done: $ok ok, $fail failed, $skip skipped  $(date +%H:%M:%S) ==="
