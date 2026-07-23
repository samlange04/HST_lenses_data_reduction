#!/usr/bin/env bash
# Regenerate every WFC3/IR F160W product with --align mast (delivered FIT_REL WCS,
# no updatewcs, no TweakReg). J0728+3835 is already done but is re-run for uniformity.
SD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; WS="$(dirname "$SD")"
LOG="$WS/data/run_logs"; mkdir -p "$LOG"
for d in "$WS"/data/calibrated/slacs/*/f160W; do
  [ -d "$d" ] || continue
  ls "$d"/*flt.fits >/dev/null 2>&1 || continue
  lens=$(basename "$(dirname "$d")")
  rm -rf "$WS/data/drizzled/slacs/$lens/f160W" "$WS/data/drizzle_files/slacs/$lens/f160W"
  printf '%-12s f160W ' "$lens"
  if conda run -n stenv python "$SD/drizzle_wfc3_ir.py" --lens "$lens" --filt f160W \
       --sample slacs --align mast > "$LOG/${lens}_f160W_wfc3.log" 2>&1; then echo OK; else echo FAILED; fi
done
echo "=== all WFC3 done $(date +%H:%M:%S) ==="
