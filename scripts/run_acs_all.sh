#!/usr/bin/env bash
# Regenerate every ACS product with --align mast (GSC242 WCS, no updatewcs/TweakReg)
# and --cr (LACosmic CR masking). Both defaults; passed explicitly for the record.
SD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; WS="$(dirname "$SD")"
LOG="$WS/data/run_logs"; mkdir -p "$LOG"
for filt in f814W f555W; do
  for d in "$WS"/data/calibrated/slacs/*/"$filt"; do
    [ -d "$d" ] || continue
    ls "$d"/*flc.fits >/dev/null 2>&1 || continue
    lens=$(basename "$(dirname "$d")")
    rm -rf "$WS/data/drizzled/slacs/$lens/$filt" "$WS/data/drizzle_files/slacs/$lens/$filt"
    printf '%-12s %-6s ' "$lens" "$filt"
    if conda run -n stenv python "$SD/drizzle_acs_wfc.py" --lens "$lens" --filt "$filt" \
         --sample slacs --cr --align mast --cr-method lacosmic \
         > "$LOG/${lens}_${filt}_acs.log" 2>&1; then echo OK; else echo FAILED; fi
  done
done
echo "=== all ACS done $(date +%H:%M:%S) ==="
