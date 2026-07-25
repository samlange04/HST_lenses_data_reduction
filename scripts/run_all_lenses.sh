#!/usr/bin/env bash
# ── RETIRED 2026-07-26 — use run_wfpc2_wf3.sh ─────────────────────────────────
# This was a second, independently-maintained WFPC2 F606W driver, and it drifted
# badly from the pipeline it was supposed to run. By the time it was retired it:
#
#   * passed no --align, so it used the drizzle script's then-default 'tweakreg' --
#     the mode the per-lens core-registration audit rejected for all 22 lenses,
#     which scatters frames ~0.7" and splits the deflector core into ~4 knots;
#   * had no split-visit handling, so J0728+3835 and J0822+2652 would have been
#     drizzled as single combined datasets across a ~15 deg roll difference, and
#     would have rewritten the per-visit tracking-JSON keys back to a bogus
#     combined 'f606W' entry;
#   * skipped align_wfpc2_to_acs.py entirely, leaving every product ~0.3-0.9" off
#     the other bands in absolute astrometry -- invisible in a single-band look;
#   * iterated all 38 SLACS lenses, including the 16 with no WFPC2 data at all.
#
# Its one feature the other runner lacked -- a retry pass over failures -- has been
# absorbed into run_wfpc2_wf3.sh. Keeping two drivers is what allowed the drift, so
# there is now one. This stub refuses rather than being deleted so that an old
# invocation fails loudly instead of hitting "command not found" and being retyped.
echo "run_all_lenses.sh is retired: it ran TweakReg alignment (rejected by the" >&2
echo "per-lens audit), skipped the astrometric tie to ACS, and mishandled the" >&2
echo "split-visit lenses. Use scripts/run_wfpc2_wf3.sh, which does drizzle ->" >&2
echo "align_wfpc2_to_acs -> make_cutouts per lens and retries failures." >&2
exit 1
