#!/usr/bin/env bash
# Drizzle WF3 F606W for every SLACS lens with WFPC2 F606W data, tie each product to
# the GAIA-accurate ACS F814W astrometry, and cut the stamps.
#
# The full documented run order is drizzle -> align_wfpc2_to_acs.py -> make_cutouts.py.
# The align step is NOT optional and was missing from this script until 2026-07-26:
# WFPC2 F606W carries only a GSC240 solution, ~0.3-0.9" off in absolute astrometry, so
# skipping it yields stamps that look perfect in isolation (right centre, calibrated
# noise, no warning) and are misregistered against F814W/F160W. A re-drizzle also
# discards any previous tie, because align_wfpc2_to_acs.py edits CRVAL1/2 in the
# drizzled product -- so the tie must be re-applied after every drizzle, not once.
#
# All 22 lenses are included. J0728+3835 was previously excluded for having only 2
# exposures on one x-phase; with the -COPY visit it has 6. No exclusion list is needed
# because the drizzle script measures each lens's dither coverage itself and skips any
# lens that cannot reach 0.05"/px.
#
# Written for macOS's bash 3.2: no associative arrays, no `[[ -v ]]`.
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WS="$(dirname "$SCRIPT_DIR")"
LOGDIR="$WS/data/run_logs"; mkdir -p "$LOGDIR"
ALIGN_JSON="$WS/info/wfpc2_alignment.json"

LENSES=(J0008-0004 J0029-0055 J0157-0056 J0252+0039 J0330-0020 J0728+3835
        J0822+2652 J0841+3824 J0903+4116 J0936+0913 J0946+1006 J1020+1122
        J1023+4230 J1029+0420 J1032+5322 J1142+1001 J1213+6708 J1218+0830
        J1430+4105 J1432+6317 J1525+3327 J2341+0000)

# ── Split-visit lenses ─────────────────────────────────────────────────────────
# J0728+3835 and J0822+2652 each have two visits at a ~14-16 deg roll difference, i.e.
# two guide-star solutions. They are drizzled as SEPARATE per-visit datasets rather
# than cross-registered, so each visit is single-guide-star and can use MAST alignment.
# Each field is "PA_V3:out_suffix"; the drizzle keeps frames within 1 deg of that PA.
#   J0822+2652 -> f606W_v1 (2x1100s) + f606W_v2 (4x1100s)
#   J0728+3835 -> f606W_v2 only; its 2-frame visit has just 1 x-dither-phase and
#                 cannot reach 0.05", so the drizzle script skips it by design.
# J1142+1001 is deliberately NOT here: its two visits share a roll (PA 119.00 vs
# 118.87), so there is no offset to separate and it stays combined.
split_visits_for() {
  case "$1" in
    J0728+3835) echo "91.374092:_v2" ;;
    J0822+2652) echo "101.8479:_v1 87.918716:_v2" ;;
    *)          echo "" ;;
  esac
}

# Per-lens alignment mode from the core-registration audit. Every lens currently reads
# "mast"; the file is still consulted per lens rather than hardcoded so that a future
# lens whose audit picks "tweakreg" is honoured without editing this script. Read once
# into a lookup table -- one python start-up, not one per lens.
ALIGN_TBL="$(conda run -n stenv python -c "
import json
for k, v in json.load(open('$ALIGN_JSON')).items():
    print(k, v)")" || { echo "cannot read $ALIGN_JSON"; exit 1; }

align_for() {
  local a
  a="$(printf '%s\n' "$ALIGN_TBL" | awk -v l="$1" '$1==l {print $2; exit}')"
  # A lens absent from the audit falls back to mast, which is what all 22 audited
  # lenses chose; tweakreg is never a safe default (it erases the dither).
  [ -n "$a" ] && echo "$a" || echo mast
}

# drizzle -> align -> cutout for one product (one lens, one visit).
# $1 lens, $2 align mode, $3 extra drizzle args, $4 output dir suffix
run_product() {
  local lens="$1" align="$2" extra="$3" suffix="$4"
  local key="f606W${suffix}"
  local log="$LOGDIR/${lens}_${key}_wf3.log"

  printf '  %-10s align=%-8s ' "$key" "$align"
  if ! conda run -n stenv python "$SCRIPT_DIR/drizzle_wfpc2_wf3.py" --lens "$lens" \
         --filt f606W --sample slacs --align "$align" $extra > "$log" 2>&1; then
    # A lens skipped for insufficient dither phase is an intended outcome, not a
    # failure -- the script exits non-zero without writing, so separate the two.
    if grep -q 'cannot drizzle to' "$log"; then
      echo "SKIPPED (dither phase)"; return 0
    fi
    echo "DRIZZLE FAILED (see $log)"; return 1
  fi

  # Absolute-astrometry tie to ACS F814W. Idempotent (re-measures the residual and
  # applies ~0), so it is safe on a product the drizzle skipped as already-existing.
  if ! conda run -n stenv python "$SCRIPT_DIR/align_wfpc2_to_acs.py" --lens "$lens" \
         --f606-dir "$key" >> "$log" 2>&1; then
    echo "ALIGN FAILED (see $log)"; return 1
  fi

  if ! conda run -n stenv python "$SCRIPT_DIR/make_cutouts.py" --lens "$lens" \
         --filt "$key" --sample slacs >> "$log" 2>&1; then
    echo "CUTOUT FAILED (see $log)"; return 1
  fi
  echo "OK"
}

FAILED=()
for lens in "${LENSES[@]}"; do
  echo "=== $lens $(date +%H:%M:%S) ==="
  align="$(align_for "$lens")"
  visits="$(split_visits_for "$lens")"
  if [ -n "$visits" ]; then
    for visit in $visits; do
      pa="${visit%%:*}"; suffix="${visit##*:}"
      run_product "$lens" "$align" "--pa $pa --out-suffix $suffix" "$suffix" \
        || FAILED+=("$lens$suffix")
    done
  else
    run_product "$lens" "$align" "" "" || FAILED+=("$lens")
  fi
done

# ── Retry pass ────────────────────────────────────────────────────────────────
# Absorbed from the retired run_all_lenses.sh. Most failures here are transient MAST
# or CRDS timeouts, which a second attempt clears; a lens that fails twice needs its
# log read. Re-running is safe because every stage is idempotent (the drizzle skips
# existing products, the astrometric tie re-measures and applies ~0).
if [ ${#FAILED[@]} -gt 0 ]; then
  echo ""
  echo "=== Retrying ${#FAILED[@]} failed: ${FAILED[*]} ==="
  RETRY=("${FAILED[@]}"); FAILED=()
  for item in "${RETRY[@]}"; do
    lens="${item%%_v*}"
    align="$(align_for "$lens")"
    echo "=== $lens $(date +%H:%M:%S) ==="
    case "$item" in
      *_v*) suffix="_v${item##*_v}"
            pa="$(split_visits_for "$lens" | tr ' ' '\n' | awk -F: -v s="$suffix" '$2==s {print $1}')"
            run_product "$lens" "$align" "--pa $pa --out-suffix $suffix" "$suffix" \
              || FAILED+=("$item") ;;
      *)    run_product "$lens" "$align" "" "" || FAILED+=("$item") ;;
    esac
  done
fi

echo "=== all done $(date +%H:%M:%S) ==="
if [ ${#FAILED[@]} -gt 0 ]; then
  echo "STILL FAILED after retry (${#FAILED[@]}): ${FAILED[*]}"
  exit 1
fi
