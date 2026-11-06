#!/usr/bin/env bash
# Run psf_star_check.py on every band of a sample that has a PSF record (the model against
# the real PSF stars, one twin per star). Usage: run_psf_star_check_all.sh [SAMPLE] [--blur S]
# Logs to data/run_logs/<lens>_<filt>_starcheck.log. Bands with no stars are quick no-ops.
set -u
SD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WS="$(dirname "$SD")"
SAMPLE="${1:-slacs}"
shift || true
cd "$WS"
mkdir -p data/run_logs
ok=0; failed=0
while read -r lens filt; do
  log="data/run_logs/${lens}_${filt}_starcheck.log"
  if uv run python scripts/psf_star_check.py --sample "$SAMPLE" --lens "$lens" --filt "$filt" "$@" > "$log" 2>&1; then
    ok=$((ok + 1)); tail -1 "$log"
  else
    failed=$((failed + 1)); echo "FAILED $lens $filt (see $log)"
  fi
done < <(uv run python - "$SAMPLE" <<'EOF'
import json, sys
d = json.load(open('info/provenance/lens_psf.json')).get(sys.argv[1], {})
for lens, filts in sorted(d.items()):
    for filt, rec in sorted(filts.items()):
        if rec:
            print(lens, filt)
EOF
)
echo "=== star check ($SAMPLE): $ok ok, $failed failed ==="
