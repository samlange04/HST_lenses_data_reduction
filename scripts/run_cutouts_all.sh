#!/usr/bin/env bash
# Regenerate every cutout. make_cutouts.py defaults to --pass auto, which cuts from
# the LACosmic CR pass where one exists (optical: clean, ~99% core flux retained) and
# falls back to no-CR for F160W, which has no CR pass.
SD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; WS="$(dirname "$SD")"
LOG="$WS/data/run_logs"; mkdir -p "$LOG"
ok=0; fail=0
for filt in f606W f814W f555W f160W; do
  for d in "$WS"/data/drizzled/slacs/*/"$filt"; do
    [ -d "$d" ] || continue
    ls "$d"/*nocrrej*sci.fits >/dev/null 2>&1 || continue
    lens=$(basename "$(dirname "$d")")
    if conda run -n stenv python "$SD/make_cutouts.py" --lens "$lens" --filt "$filt" \
         --sample slacs > "$LOG/${lens}_${filt}_cut.log" 2>&1; then
      ok=$((ok+1))
    else
      fail=$((fail+1)); echo "FAILED: $lens $filt (see $LOG/${lens}_${filt}_cut.log)"
    fi
  done
done
echo "=== cutouts done: $ok ok, $fail failed  $(date +%H:%M:%S) ==="
