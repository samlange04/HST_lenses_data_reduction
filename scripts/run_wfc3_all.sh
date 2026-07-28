#!/usr/bin/env bash
# Drizzle every WFC3/IR F160W product in a sample, with --align mast (the delivered
# FIT_REL WCS; no updatewcs, no TweakReg) and --dq-refine at its default.
#
# Usage: run_wfc3_all.sh [SAMPLE]       (default: mast_target_names.DEFAULT_SAMPLE)
#
# No --cr: the IR PSF is ~1 px FWHM at 0.1283"/px, so both LACosmic and driz_cr eat
# point sources, and the FLTs are already up-the-ramp CR-rejected. See CLAUDE.md.
#
# The lens list comes from info/lens_samples.json via scripts/mast_target_names.py --
# NOT from globbing data/calibrated/, which only ever re-runs what is already on disk
# and so does nothing at all after a wipe. Every lens is tried on every run;
# drizzle_wfc3_ir.py exits 0 having recorded `null` where MAST has no data, so "no data"
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
echo "=== WFC3/IR: ${#LENSES[@]} lenses in sample '$SAMPLE' ==="

ok=0; nodata=0; blocked=0; fail=0; FAILED=()
for lens in "${LENSES[@]}"; do
  log="$LOG/${lens}_f160W_wfc3.log"
  rm -rf "$WS/data/drizzled/$SAMPLE/$lens/f160W" "$WS/data/drizzle_files/$SAMPLE/$lens/f160W"
  printf '%-12s f160W ' "$lens"
  if conda run -n stenv python "$SD/drizzle_wfc3_ir.py" --lens "$lens" --filt f160W \
       --sample "$SAMPLE" --align mast > "$log" 2>&1; then
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
    echo "FAILED (see $log)"; fail=$((fail + 1)); FAILED+=("$lens")
  fi
done

echo "=== WFC3/IR $SAMPLE done $(date +%H:%M:%S): $ok ok, $nodata no data, $blocked blocked (exptime), $fail failed ==="
if [ "$fail" -gt 0 ]; then printf '  FAILED: %s\n' "${FAILED[@]}"; exit 1; fi
