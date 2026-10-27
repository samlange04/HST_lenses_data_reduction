#!/usr/bin/env bash
# Regenerate every stamp. make_stamps.py defaults to --pass auto, which cuts from
# the LACosmic CR pass where one exists (optical: clean, ~99% core flux retained) and
# falls back to no-CR for F160W, which has no CR pass.
#
# Usage: run_stamps_all.sh [SAMPLE] [SIZE_ARCSEC] [--bcfill]
#          (defaults: mast_target_names.DEFAULT_SAMPLE, 12)
#
# --bcfill cuts from the bad-column-filled drizzles (data/drizzled_bcfill/) into the SAME
# data/stamps/ tree under a TAGGED name (stamp12_bcfill_{sci,noise}.fits; see AGENTS.md
# *Bad-column fill* and stamp_paths.py). Only ACS/WFPC2 bands have a bcfill drizzle; a band
# without one globs to nothing, which is harmless. A band keeps ONE stamp: make_stamps.py
# refuses to put a second reduction's stamp beside an existing one unless --force, which
# this runner does not pass -- so on main (standard stamps) a --bcfill run refuses every
# band, and on the bcfill branch a plain run does. That refusal is the point.
#
# SIZE_ARCSEC, if given, is passed to make_stamps.py --size. Omitted, each lens is cut at its
# own modelling size (its stamp_size in info/samples/lens_samples.json, else 12"). A size
# writes stamp<size>_* files beside the modelling ones (see scripts/stamp_paths.py), so
# running at a second size adds a product set rather than replacing the standard one; the
# per-lens logs are tagged with it too, for the same reason.
#
# Unlike the drizzle runners, this one globs data/drizzled/ rather than the sample list,
# and that is correct: a stamp can only be cut from a mosaic that exists, so the products
# on disk -- not the roster of lenses -- are the right thing to iterate.
SD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; WS="$(dirname "$SD")"
LOG="$WS/data/run_logs"; mkdir -p "$LOG"
# Pull the optional --bcfill flag out of the args, leaving SAMPLE=$1 and SIZE=$2. bcfill reads
# data/drizzled_bcfill/ and writes tagged stamps into data/stamps/ (make_stamps.py --bcfill).
BCFILL=""; ARGS=()
for _arg in "$@"; do
  if [ "$_arg" = "--bcfill" ]; then BCFILL="--bcfill"; else ARGS+=("$_arg"); fi
done
set -- ${ARGS[@]+"${ARGS[@]}"}
SAMPLE="$(uv run --project "$WS" python "$SD/mast_target_names.py" ${1:+"$1"} --print-sample)" || exit 1
SIZE="${2:-}"
SIZE_ARG=(); [ -n "$SIZE" ] && SIZE_ARG=(--size "$SIZE")
DZ="drizzled"; [ -n "$BCFILL" ] && DZ="drizzled_bcfill"
# Log-name tag, empty at the default size so existing log paths are unchanged.
TAG=""; [ -z "$SIZE" ] || TAG="_${SIZE}arcsec"
TAG="${TAG}${BCFILL:+_bcfill}"
#
# The globs are "$filt"* , not "$filt", so a split lens's secondary-visit product
# directory (f606W_v2 on J0822+2652; its primary visit and J0728+3835's sole visit are
# now bare f606W, matched by either form) is still included. With a bare "$filt" the
# v2 dir matched nothing and its stamps were silently never regenerated -- no error,
# just left on stale stamps. The stamp --filt is taken from the directory name for
# the same reason.
ok=0; fail=0
# SLACS bands (f606W f814W f555W f160W), the gallery WFC3/UVIS blue band (f438W) and the
# SL2S WFC3/UVIS bands (f475X f600LP). A filter not present in the sample globs to nothing,
# so listing all known bands here is harmless and keeps this one runner correct for every
# sample. make_stamps reads the centre band's *mosaic*, not its stamp, so the order only
# needs that drizzle to exist -- which the drizzle phase guarantees.
for filt in f606W f814W f555W f160W f438W f475X f600LP; do
  for d in "$WS"/data/"$DZ"/"$SAMPLE"/*/"$filt"*; do
    [ -d "$d" ] || continue
    # Any final drizzle product (CR or no-CR) means the mosaic exists and is cuttable.
    # Must NOT key on *nocrrej* specifically: ACS/WFPC2 now default to a CR-only pass, so
    # their dirs hold only *_cr_*; F160W (no CR pass) holds only *nocrrej*. make_stamps
    # --pass auto picks the right one.
    ls "$d"/*_sci.fits >/dev/null 2>&1 || continue
    lens=$(basename "$(dirname "$d")")
    key=$(basename "$d")
    if uv run --project "$WS" python "$SD/make_stamps.py" --lens "$lens" --filt "$key" \
         --sample "$SAMPLE" ${SIZE_ARG[@]+"${SIZE_ARG[@]}"} $BCFILL > "$LOG/${lens}_${key}${TAG}_cut.log" 2>&1; then
      ok=$((ok+1))
    else
      fail=$((fail+1)); echo "FAILED: $lens $key (see $LOG/${lens}_${key}${TAG}_cut.log)"
    fi
  done
done
echo "=== stamps $SAMPLE ${SIZE:-modelling-size} done: $ok ok, $fail failed  $(date +%H:%M:%S) ==="
