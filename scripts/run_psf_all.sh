#!/usr/bin/env bash
# Build a PSF for every drizzled product. make_psf.py defaults to --pass auto (detect
# stars in the LACosmic CR mosaic where one exists, else the no-CR mosaic) and --method
# auto (empirical ePSF, STDPSF model fallback when a field is too star-poor).
#
# Usage: run_psf_all.sh [SAMPLE]    (default: mast_target_names.DEFAULT_SAMPLE)
#
# Like run_cutouts_all.sh -- and unlike the drizzle runners -- this globs data/drizzled/
# rather than the roster: a PSF can only be built from a mosaic that exists, so the
# products on disk are the right thing to iterate. Written for macOS bash 3.2.
SD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; WS="$(dirname "$SD")"
LOG="$WS/data/run_logs"; mkdir -p "$LOG"
SAMPLE="$(conda run -n stenv python "$SD/mast_target_names.py" ${1:+"$1"} --print-sample)" || exit 1
#
# Globs are "$filt"* , not "$filt", so per-visit dirs (f606W_v1, f606W_v2) are included;
# --filt is taken from the directory basename. Any *_sci.fits means the mosaic is buildable
# (ACS/WFPC2 hold *_cr_*; F160W holds *nocrrej* -- make_psf --pass auto picks the right one).
ok=0; fail=0
for filt in f606W f814W f555W f160W; do
  for d in "$WS"/data/drizzled/"$SAMPLE"/*/"$filt"*; do
    [ -d "$d" ] || continue
    ls "$d"/*_sci.fits >/dev/null 2>&1 || continue
    lens=$(basename "$(dirname "$d")")
    key=$(basename "$d")
    if conda run -n stenv python "$SD/make_psf.py" --lens "$lens" --filt "$key" \
         --sample "$SAMPLE" > "$LOG/${lens}_${key}_psf.log" 2>&1; then
      ok=$((ok+1))
    else
      fail=$((fail+1)); echo "FAILED: $lens $key (see $LOG/${lens}_${key}_psf.log)"
    fi
  done
done
echo "=== psf $SAMPLE done: $ok ok, $fail failed  $(date +%H:%M:%S) ==="
