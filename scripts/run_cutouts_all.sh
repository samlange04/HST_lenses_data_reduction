#!/usr/bin/env bash
# Regenerate every cutout. make_cutouts.py defaults to --pass auto, which cuts from
# the LACosmic CR pass where one exists (optical: clean, ~99% core flux retained) and
# falls back to no-CR for F160W, which has no CR pass.
#
# Usage: run_cutouts_all.sh [SAMPLE]    (default: mast_target_names.DEFAULT_SAMPLE)
#
# Unlike the drizzle runners, this one globs data/drizzled/ rather than the sample list,
# and that is correct: a stamp can only be cut from a mosaic that exists, so the products
# on disk -- not the roster of lenses -- are the right thing to iterate.
SD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; WS="$(dirname "$SD")"
LOG="$WS/data/run_logs"; mkdir -p "$LOG"
SAMPLE="$(conda run -n stenv python "$SD/mast_target_names.py" ${1:+"$1"} --print-sample)" || exit 1
#
# The globs are "$filt"* , not "$filt", so the per-visit product directories
# (f606W_v1, f606W_v2 on J0822+2652 and J0728+3835) are included. With a bare
# "$filt" they matched nothing for those lenses and their stamps were silently
# never regenerated -- no error, just two lenses quietly left on stale cutouts.
# The cutout --filt is taken from the directory name for the same reason.
ok=0; fail=0
for filt in f606W f814W f555W f160W; do
  for d in "$WS"/data/drizzled/"$SAMPLE"/*/"$filt"*; do
    [ -d "$d" ] || continue
    # Any final drizzle product (CR or no-CR) means the mosaic exists and is cuttable.
    # Must NOT key on *nocrrej* specifically: ACS/WFPC2 now default to a CR-only pass, so
    # their dirs hold only *_cr_*; F160W (no CR pass) holds only *nocrrej*. make_cutouts
    # --pass auto picks the right one.
    ls "$d"/*_sci.fits >/dev/null 2>&1 || continue
    lens=$(basename "$(dirname "$d")")
    key=$(basename "$d")
    if conda run -n stenv python "$SD/make_cutouts.py" --lens "$lens" --filt "$key" \
         --sample "$SAMPLE" > "$LOG/${lens}_${key}_cut.log" 2>&1; then
      ok=$((ok+1))
    else
      fail=$((fail+1)); echo "FAILED: $lens $key (see $LOG/${lens}_${key}_cut.log)"
    fi
  done
done
echo "=== cutouts $SAMPLE done: $ok ok, $fail failed  $(date +%H:%M:%S) ==="
