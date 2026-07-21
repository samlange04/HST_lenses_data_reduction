#!/usr/bin/env bash
#
# Rebuild stenv as a NATIVE arm64 (Apple Silicon) conda env.
#
# WHY: the current `stenv` is an x86_64 env running under Rosetta 2 on this
# arm64 M3 Max. Rosetta-translated processes wedge into an unkillable macOS
# U-state (uninterruptible sleep) on large write() syscalls — this is what
# hangs drizzle_acs_wfc.py / drizzle_wfc3_ir.py (big 68 MB mosaics) while the
# small WFPC2/PC outputs (~2 MB) squeak through. A native arm64 build removes
# the translation layer and the whole class of hangs.
#
# ROOT GOTCHA: `conda config --show subdir` is globally pinned to osx-64, which
# forces x86_64 packages even on arm64 hardware. Every step below therefore
# sets CONDA_SUBDIR=osx-arm64 explicitly and pins it into the new env.
#
# Run this AFTER rebooting to clear the hung x86 process. Safe to run while the
# old `stenv` still exists — this creates a separate env and does not touch it.
set -euo pipefail

NEW_ENV="stenv_arm64"
YAML="$HOME/Downloads/stenv-macOS-ARM64-py3.12-2026.05.08.yaml"   # already downloaded; py3.12 matches current env

echo ">>> 0. Confirm host is arm64 (must print: arm64)"
uname -m

echo ">>> 1. Create the env forcing arm64 package resolution"
CONDA_SUBDIR=osx-arm64 conda env create -n "$NEW_ENV" -f "$YAML"

echo ">>> 2. Pin the env to arm64 so future conda/pip installs stay native"
conda activate "$NEW_ENV"
conda config --env --set subdir osx-arm64

echo ">>> 3. VERIFY the interpreter is truly native (must print: arm64  /  ...-arm64-...)"
python -c "import platform; print(platform.machine()); print(platform.platform())"
file "$(python -c 'import sys; print(sys.executable)')"   # expect: Mach-O 64-bit executable arm64

echo ">>> 4. Smoke-test the STScI stack imports"
python -c "import drizzlepac, stwcs, astropy, numpy; print('drizzlepac', drizzlepac.__version__); print('numpy', numpy.__version__)"

cat <<'NEXT'

>>> DONE. Next steps (manual):
  A. Re-run the ACS lens that hung, in the NEW env, and confirm it does NOT wedge:
       conda run -n stenv_arm64 python scripts/drizzle_acs_wfc.py --lens J0008-0004 --filt f814W --sample slacs
     Watch that the process stays R/S (never U) through the 68 MB drizzle write.
  B. Once confident, retire the x86 env and adopt the arm64 one. Either:
       - update CLAUDE.md + scripts to use `stenv_arm64`, OR
       - `conda env remove -n stenv` then rename:
           conda create -n stenv --clone stenv_arm64 && conda env remove -n stenv_arm64
  C. Consider unpinning the global x86 default so this can't recur:
       conda config --set subdir osx-arm64      # global default -> native
NEXT