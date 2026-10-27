#!/usr/bin/env bash
# Drizzle every WFC3/UVIS product in a sample, the gallery and SL2S filters, with the script
# defaults (--align mast, --cr LACosmic, native 0.0396"/px pixfrac 1.0 since 2026-09-29). Passed
# explicitly for the record.
#
# Usage: run_gallery_uvis_all.sh [SAMPLE]   (default: gallery)
#
# Filters come from the sample's `uvis_filters` (gallery: F606W F814W F438W; SL2S: F475X
# F600LP). Gallery coverage is sparse and per-lens: F606W is on all 15 lenses, F814W/F438W
# on 6 each (its F225W/F275W products were deleted 2026-09-29). Every
# lens is tried in every filter; drizzle_wfc3_uvis.py exits 0 having recorded `null`
# where MAST has no data, so the many empty lens+filter combos report as "no data", not
# failures. F606W is drizzled first (the primary band and the stamp --center-band).
# The SL2S samples (sl2s_dinos_a/x, sl2s_wf, sl2s_snap) add F475X (GO-17130 Dinos II, or
# the GO-11588 snapshot) and F600LP (GO-11588 snapshot); each sample just reports "no
# data" for the other set's filters. Which program's F475X is used is decided per lens by
# the `proposals` quirk in info/samples/lens_samples.json, not here. The filter list itself comes
# from the sample's `uvis_filters` (mast_target_names.uvis_filters; default f606W f814W
# f438W): SL2S runs only f475X f600LP, because its F606W/F814W belong to WFPC2/ACS and one
# filter key holds one instrument per lens. Looping F814W over SL2S on 2026-10-04 drizzled
# J0959+0206's ACS FLCs as UVIS and nulled its WFPC2 F606W records.
#
# Unlike run_acs_all.sh this does NOT rm the output dir first: it relies on the drizzle
# script's idempotent skip (final product exists -> skip) so a multi-GB gallery campaign
# is resumable after an interruption. To force a re-run, delete the lens's dir under
# data/drizzled/ (and data/drizzle_files/) as AGENTS.md describes.
#
# The lens list comes from info/samples/lens_samples.json via scripts/mast_target_names.py --
# NOT from globbing data/calibrated/, which does nothing after a wipe.
#
# Written for macOS's bash 3.2: no `mapfile`, no associative arrays, no `[[ -v ]]`.
SD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; WS="$(dirname "$SD")"
LOG="$WS/data/run_logs"; mkdir -p "$LOG"

SAMPLE="$(uv run --project "$WS" python "$SD/mast_target_names.py" "${1:-gallery}" --print-sample)" || exit 1
LENSES=()
while IFS= read -r _l; do
  [ -n "$_l" ] && LENSES+=("$_l")
done < <(uv run --project "$WS" python "$SD/mast_target_names.py" "$SAMPLE")
[ "${#LENSES[@]}" -gt 0 ] || { echo "No lenses in sample '$SAMPLE'" >&2; exit 1; }
echo "=== WFC3/UVIS: ${#LENSES[@]} lenses in sample '$SAMPLE' ==="

ok=0; nodata=0; blocked=0; fail=0; FAILED=()
UVIS_FILTERS="$(uv run --project "$WS" python "$SD/mast_target_names.py" "$SAMPLE" --uvis-filters)" || exit 1
for filt in $UVIS_FILTERS; do
  for lens in "${LENSES[@]}"; do
    log="$LOG/${lens}_${filt}_uvis.log"
    printf '%-12s %-6s ' "$lens" "$filt"
    if uv run --project "$WS" python "$SD/drizzle_wfc3_uvis.py" --lens "$lens" --filt "$filt" \
         --sample "$SAMPLE" --cr --align mast --cr-method lacosmic > "$log" 2>&1; then
      # An exit-0 run that wrote nothing is either "MAST has no data" or "total exptime
      # below BLOCK_EXPTIME" -- both are ordinary outcomes, not a product and not a failure.
      if grep -q '^=== NO DATA:' "$log"; then
        echo "no data"; nodata=$((nodata + 1))
      elif grep -q '^=== BLOCKED (exptime):' "$log"; then
        echo "blocked (exptime)"; blocked=$((blocked + 1))
      elif grep -q 'drizzled products already exist, skipping' "$log"; then
        echo "OK (cached)"; ok=$((ok + 1))
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

echo "=== UVIS $SAMPLE done $(date +%H:%M:%S): $ok ok, $nodata no data, $blocked blocked (exptime), $fail failed ==="
if [ "$fail" -gt 0 ]; then printf '  FAILED: %s\n' "${FAILED[@]}"; exit 1; fi
