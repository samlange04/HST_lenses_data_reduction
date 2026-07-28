#!/usr/bin/env bash
# Drizzle every ACS/WFC product in a sample, with --align mast (GSC242/GAIAeDR3 WCS, no
# updatewcs/TweakReg) and --cr (LACosmic CR masking). Both are the script defaults;
# passed explicitly for the record.
#
# Usage: run_acs_all.sh [SAMPLE]        (default: mast_target_names.DEFAULT_SAMPLE)
#
# The lens list comes from info/lens_samples.json via scripts/mast_target_names.py --
# NOT from globbing data/calibrated/. That matters: globbing the download directory only
# ever re-runs lenses already on disk, so after a wipe it silently does nothing, and a
# lens newly added to the sample is never picked up. Every lens is tried on every run and
# drizzle_acs_wfc.py exits 0 having recorded `null` where MAST has no data, so "no data"
# is reported separately and does not count as a failure.
#
# Written for macOS's bash 3.2: no `mapfile`, no associative arrays, no `[[ -v ]]`.
SD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; WS="$(dirname "$SD")"
LOG="$WS/data/run_logs"; mkdir -p "$LOG"

SAMPLE="$(conda run -n stenv python "$SD/mast_target_names.py" ${1:+"$1"} --print-sample)" || exit 1
LENSES=()
while IFS= read -r _l; do
  [ -n "$_l" ] && LENSES+=("$_l")
done < <(conda run -n stenv python "$SD/mast_target_names.py" "$SAMPLE")
[ "${#LENSES[@]}" -gt 0 ] || { echo "No lenses in sample '$SAMPLE'" >&2; exit 1; }
echo "=== ACS: ${#LENSES[@]} lenses in sample '$SAMPLE' ==="

ok=0; nodata=0; blocked=0; fail=0; FAILED=()
for filt in f814W f555W; do
  for lens in "${LENSES[@]}"; do
    log="$LOG/${lens}_${filt}_acs.log"
    rm -rf "$WS/data/drizzled/$SAMPLE/$lens/$filt" "$WS/data/drizzle_files/$SAMPLE/$lens/$filt"
    printf '%-12s %-6s ' "$lens" "$filt"
    if conda run -n stenv python "$SD/drizzle_acs_wfc.py" --lens "$lens" --filt "$filt" \
         --sample "$SAMPLE" --cr --align mast --cr-method lacosmic > "$log" 2>&1; then
      # An exit-0 run that wrote nothing is either "MAST has no data" or "total exptime
      # below BLOCK_EXPTIME" -- both are ordinary outcomes, not a product and not a failure.
      if grep -q '^=== NO DATA:' "$log"; then
        echo "no data"; nodata=$((nodata + 1))
      elif grep -q '^=== BLOCKED (exptime):' "$log"; then
        echo "blocked (exptime)"; blocked=$((blocked + 1))
      else
        if grep -q '^  EXPTIME WARNING:' "$log"; then
          echo "OK (low exptime)"
        else
          echo OK
        fi
        ok=$((ok + 1))
      fi
    else
      echo "FAILED (see $log)"; fail=$((fail + 1)); FAILED+=("$lens $filt")
    fi
  done
done

echo "=== ACS $SAMPLE done $(date +%H:%M:%S): $ok ok, $nodata no data, $blocked blocked (exptime), $fail failed ==="
if [ "$fail" -gt 0 ]; then printf '  FAILED: %s\n' "${FAILED[@]}"; exit 1; fi
