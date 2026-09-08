#!/usr/bin/env bash
# Mark lensed-image positions across every sample, in one invocation.
#
# Usage: run_positions_all.sh [SAMPLE] [extra make_positions.py flags...]
#   run_positions_all.sh                     # every sample, best band per lens
#   run_positions_all.sh slacs_gold          # one sample
#   run_positions_all.sh slacs_gold --force  # re-mark, one sample
#   run_positions_all.sh --filt f606W        # every sample, forced band (see note)
#
# make_positions.py already loops every lens WITHIN a sample and marks positions once per
# lens on its best band, broadcasting the result to that lens's other bands (positions are
# band-independent arcsec coords -- see the script header / AGENTS.md). So this wrapper's
# only job is to sweep the samples; with no SAMPLE arg it does all three science samples in
# order. It is INTERACTIVE (a blocking Tk window per lens) and deliberately does NOT
# redirect to log files -- you need to see the GUI and the per-click console output. It is
# resumable: an already-marked lens is skipped, so re-running picks up where you left off.
#
# The first arg is treated as SAMPLE only if it does not start with '-'; anything else
# (and every arg after the sample) is forwarded verbatim to make_positions.py, so --filt,
# --force, --pass, --variant, --size, --search-box-size, --display etc. all pass through.
SD="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; WS="$(dirname "$SD")"

# All science samples with cutouts. Positions are a lens-modelling product, so every sample
# is in scope (unlike the drizzle runners, which are instrument-specific).
ALL_SAMPLES="slacs_gold slacs_other gallery"

if [ $# -gt 0 ] && [ "${1#-}" = "$1" ]; then
  SAMPLES="$1"; shift            # explicit sample; remaining args pass through
else
  SAMPLES="$ALL_SAMPLES"         # no sample given: sweep them all
fi

for sample in $SAMPLES; do
  echo "=== positions: $sample  $(date +%H:%M:%S) ==="
  uv run --project "$WS" python "$SD/make_positions.py" --sample "$sample" "$@" || {
    echo "make_positions.py exited non-zero for $sample" >&2; exit 1; }
done
echo "=== positions done: $SAMPLES ==="
