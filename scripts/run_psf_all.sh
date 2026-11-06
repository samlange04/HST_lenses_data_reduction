#!/usr/bin/env bash
# Build a PSF for every drizzled product. make_psf.py defaults to --pass auto (the kernel
# is named for the CR pass; stars are detected in the no-CR mosaic where one exists, since
# LACosmic eats star cores) and --method auto, which since 2026-11-03 ships the HYBRID
# kernel wherever a model exists: the model injected into the band's frames and drizzled
# (the grid mean) as the core, the field-star ePSF as the wing; the pure model where the
# field has no usable stars (WFPC2), the ePSF where the band has no model (SL2S F475X /
# F600LP). See AGENTS.md *PSF generation* / *The 2026-11-03 PSF audit*. Run
# make_psf_err_injected.py afterwards for the model-tier / hybrid error maps (and
# psf_star_check.py --error-maps before it, for the measured model term).
#
# Usage: run_psf_all.sh [SAMPLE] [--all|--models-only|--empirical-only]
#   (default SAMPLE: mast_target_names.DEFAULT_SAMPLE;  default mode: --all)
#
#   --all          (default) (re)build every product -- empirical AND model.
#   --models-only  skip any product whose EXISTING lens_psf.json record is `empirical`
#                  (since 2026-11-03 that is only a band with no model: SL2S F475X/F600LP,
#                  or a pre-2026-11-03 record); (re)build the rest plus any product not yet
#                  built. The skip set is "empirical AND already recorded", so a lens with
#                  no entry still builds.
#   --empirical-only  build every product with make_psf.py --method empirical: no model
#                  fallback, so a too-star-poor field FAILS and is listed for review. Used
#                  for the SL2S samples (2026-10-04), whose F475X/F600LP bands have no
#                  STDPSF grid yet; what to do with the failures is decided afterwards.
#
# Like run_stamps_all.sh -- and unlike the drizzle runners -- this globs data/drizzled/
# rather than the roster: a PSF can only be built from a mosaic that exists, so the
# products on disk are the right thing to iterate. Written for macOS bash 3.2.
SD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; WS="$(dirname "$SD")"
LOG="$WS/data/run_logs"; mkdir -p "$LOG"

MODE="all"; SAMPLE_ARG=""
for arg in "$@"; do
  case "$arg" in
    --models-only) MODE="models" ;;
    --all)         MODE="all" ;;
    --empirical-only) MODE="empirical" ;;
    -*)            echo "unknown flag: $arg (use --all, --models-only or --empirical-only)" >&2; exit 2 ;;
    *)             SAMPLE_ARG="$arg" ;;
  esac
done
SAMPLE="$(uv run --project "$WS" python "$SD/mast_target_names.py" ${SAMPLE_ARG:+"$SAMPLE_ARG"} --print-sample)" || exit 1

# --models-only skip set: "<lens>/<key>" for every product recorded as empirical in
# info/provenance/lens_psf.json (an empirical record implies a built product). Newline-delimited so the
# per-product test below is a shell `case` match -- no python/grep spawned inside the loop.
EMP_SKIP=""
if [ "$MODE" = models ]; then
  EMP_SKIP="$(uv run --project "$WS" python -c "import json,os; p=os.path.join('$WS','info', 'provenance', 'lens_psf.json'); d=json.load(open(p)) if os.path.exists(p) else {}; s=d.get('$SAMPLE', {}); print(chr(10).join(l+'/'+k for l in s for k,v in (s[l] or {}).items() if v and str(v.get('method','')).startswith('empirical')))")" || exit 1
fi

# Globs are "$filt"* , not "$filt", so a split lens's secondary-visit dir (f606W_visit2) is
# included alongside its primary/bare f606W;
# --filt is taken from the directory basename. Any *_sci.fits means the mosaic is buildable
# (ACS/WFPC2 hold *_cr_*; F160W holds *nocrrej* -- make_psf --pass auto picks the right one).
# f438W is gallery's blue UVIS science band; F225W/F275W are deliberately omitted -- they
# are confirmed unusable for lens science sample-wide (arc undetected), so no PSF work is
# done on them (see AGENTS.md, gallery_uv_bands_unusable). f475X/f600LP are the SL2S
# WFC3/UVIS bands. The extra filters are no-ops for samples that lack them (the glob
# matches nothing).
METHOD_ARG=""; [ "$MODE" = empirical ] && METHOD_ARG="--method empirical"
ok=0; fail=0; skip=0
for filt in f606W f814W f555W f160W f438W f450W f475X f600LP; do
  for d in "$WS"/data/drizzled/"$SAMPLE"/*/*/*/"$filt"*; do
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
         --sample "$SAMPLE" $METHOD_ARG > "$LOG/${lens}_${key}_psf.log" 2>&1; then
      ok=$((ok+1))
    else
      fail=$((fail+1)); echo "FAILED: $lens $key (see $LOG/${lens}_${key}_psf.log)"
    fi
  done
done
echo "=== psf $SAMPLE ($MODE) done: $ok ok, $fail failed, $skip skipped  $(date +%H:%M:%S) ==="
