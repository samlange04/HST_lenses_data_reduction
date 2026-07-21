#!/usr/bin/env bash
# Drizzle WF3 F606W + cutouts for every SLACS lens with WFPC2 F606W data.
# All 22 are included: J0728+3835 was previously excluded for having only 2
# exposures on one x-phase, but once the -COPY visit is included it has 6 exposures
# with 5 phases per axis. The drizzle script itself measures the dither coverage and
# skips any lens that cannot reach 0.05"/px, so no exclusion list is needed here.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WS="$(dirname "$SCRIPT_DIR")"
LOGDIR="$WS/data/run_logs"; mkdir -p "$LOGDIR"
LENSES=(J0008-0004 J0029-0055 J0157-0056 J0252+0039 J0330-0020 J0728+3835
        J0822+2652 J0841+3824 J0903+4116 J0936+0913 J0946+1006 J1020+1122
        J1023+4230 J1029+0420 J1032+5322 J1142+1001 J1213+6708 J1218+0830
        J1430+4105 J1432+6317 J1525+3327 J2341+0000)
for lens in "${LENSES[@]}"; do
  echo "=== $lens $(date +%H:%M:%S) ==="
  if conda run -n stenv python "$SCRIPT_DIR/drizzle_wfpc2_wf3.py" --lens "$lens" \
       --filt f606W --sample slacs > "$LOGDIR/${lens}_f606W_wf3.log" 2>&1; then
    echo "    drizzle OK"
    if conda run -n stenv python "$SCRIPT_DIR/make_cutouts.py" --lens "$lens" \
         --filt f606W --sample slacs >> "$LOGDIR/${lens}_f606W_wf3.log" 2>&1; then
      echo "    cutout OK"
    else
      echo "    CUTOUT FAILED (see $LOGDIR/${lens}_f606W_wf3.log)"
    fi
  else
    echo "    DRIZZLE FAILED (see $LOGDIR/${lens}_f606W_wf3.log)"
  fi
done
echo "=== all done $(date +%H:%M:%S) ==="
