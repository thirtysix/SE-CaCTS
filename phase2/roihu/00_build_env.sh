#!/usr/bin/env bash
# Roihu — build the Phase-2 python environment ON THE LOGIN NODE (only it has internet for pip).
# One-off, ~2-3 min. array.slurm / reduce.slurm expect the venv at $VENV.
#
#   ssh roihu 'bash -lc "$(cat)"' < 00_build_env.sh
#
# WHY A VENV, not tykky (2026-07-21): Roihu has NO `python-data` module and NO tykky (checked), unlike
# Puhti/Mahti. But cnrose needs only pyBigWig + numpy + scipy (+ pandas for the reduce), Roihu CPU login
# and compute nodes are x86_64, and pyBigWig ships a manylinux x86_64 wheel — so a plain venv installs in
# minutes with no compilation. Built in /projappl (meant for libraries; better for many small files than
# scratch). cnrose has no >=3.10 syntax, so the system python 3.9 is fine.
set -euo pipefail
# Machine-local config: walk up to the repo root and source .env via secacts_env.sh.
# Depth-independent so this works from a staged copy on the cluster too; if the file is
# absent (partial stage), a real SECACTS_CSC_PROJECT in the environment still satisfies the
# `:?` checks below.
_d="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [ "$_d" != "/" ] && [ ! -f "$_d/secacts_env.sh" ]; do _d="$(dirname "$_d")"; done
[ -f "$_d/secacts_env.sh" ] && . "$_d/secacts_env.sh"
PROJ="${PROJ:-${SECACTS_CSC_PROJECT:?set SECACTS_CSC_PROJECT in .env, or pass PROJ=}}"
VENV="${VENV:-/projappl/$PROJ/secacts_venv}"

if [[ -x "$VENV/bin/python" ]] && "$VENV/bin/python" -c "import pyBigWig,numpy,scipy,pandas" 2>/dev/null; then
  echo "[roihu] venv already good at $VENV"
  "$VENV/bin/python" -c "import pyBigWig; print('  pyBigWig',pyBigWig.__version__)"
  exit 0
fi

echo "[roihu] building venv at $VENV (login node, ~2-3 min)"
python3 -m venv "$VENV"
"$VENV/bin/pip" install --quiet --upgrade pip
"$VENV/bin/pip" install --quiet pyBigWig numpy scipy pandas
"$VENV/bin/python" -c "import pyBigWig,numpy,scipy,pandas; \
  print('[roihu] installed pyBigWig',pyBigWig.__version__,'numpy',numpy.__version__,'scipy',scipy.__version__)"
echo "[roihu] done -> array.slurm/reduce.slurm use \$VENV=$VENV"
