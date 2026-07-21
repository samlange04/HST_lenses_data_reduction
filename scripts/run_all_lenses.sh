#!/usr/bin/env bash

LENSES=(
    J0008-0004 J0029-0055 J0157-0056 J0216-0813 J0252+0039
    J0330-0020 J0728+3835 J0737+3216 J0822+2652 J0841+3824
    J0903+4116 J0912+0029 J0936+0913 J0946+1006 J0956+5100
    J0959+0410 J1020+1122 J1023+4230 J1029+0420 J1032+5322
    J1142+1001 J1143-0144 J1205+4910 J1213+6708 J1218+0830
    J1250+0523 J1402+6321 J1420+6019 J1430+4105 J1432+6317
    J1451-0239 J1525+3327 J1627-0053 J1630+4520 J2238-0754
    J2300+0022 J2303+1422 J2341+0000
)

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
LOG_DIR="$SCRIPT_DIR/../data/run_logs"
mkdir -p "$LOG_DIR"

run_lens() {
    local lens="$1" log="$LOG_DIR/${1}.log"
    if conda run -n stenv python "$SCRIPT_DIR/drizzle_wfpc2_wf3.py" --lens "$lens" \
            > "$log" 2>&1; then
        echo "  OK"
        return 0
    else
        echo "  FAILED — see $log"
        return 1
    fi
}

# ── First pass ─────────────────────────────────────────────────────────────────
FAILED=()
TOTAL=${#LENSES[@]}
COUNT=0
for lens in "${LENSES[@]}"; do
    COUNT=$((COUNT + 1))
    echo "[$COUNT/$TOTAL] $lens"
    run_lens "$lens" || FAILED+=("$lens")
done

# ── Re-run failures ────────────────────────────────────────────────────────────
if [ ${#FAILED[@]} -gt 0 ]; then
    echo ""
    echo "=== Re-running ${#FAILED[@]} failed lenses ==="
    STILL_FAILED=()
    for lens in "${FAILED[@]}"; do
        echo "  $lens"
        run_lens "$lens" || STILL_FAILED+=("$lens")
    done

    if [ ${#STILL_FAILED[@]} -gt 0 ]; then
        echo ""
        echo "Still failed after retry:"
        for lens in "${STILL_FAILED[@]}"; do
            echo "  $lens"
        done
    fi
fi

echo ""
echo "Done ($TOTAL lenses)"
