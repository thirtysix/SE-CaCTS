#!/usr/bin/env bash
# Roihu — build the pilot's tykky container ON THE LOGIN NODE (compute nodes have no internet).
# Run this once. If `module spider deeptools` shows a maintained module, skip this and load that instead.
#
#   ssh roihu 'bash -s' < 00_build_env.sh          # or copy up and run there
set -euo pipefail
# Machine-local config: walk up to the repo root and source .env via secacts_env.sh.
# Depth-independent so this works from a staged copy on the cluster too; if the file is
# absent (partial stage), a real SECACTS_CSC_PROJECT in the environment still satisfies the
# `:?` checks below.
_d="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [ "$_d" != "/" ] && [ ! -f "$_d/secacts_env.sh" ]; do _d="$(dirname "$_d")"; done
[ -f "$_d/secacts_env.sh" ] && . "$_d/secacts_env.sh"
PROJ="${PROJ:-${SECACTS_CSC_PROJECT:?set SECACTS_CSC_PROJECT in .env, or pass PROJ=}}"
WORK="/scratch/${PROJ}/se-cacts/pilot"
PREFIX="/projappl/${PROJ}/secacts_tykky"

source /appl/profile/zz-csc-env.sh
# Prefer a maintained module if one exists (cheaper than a container):
if module -t spider 2>/dev/null | grep -qi '^deeptools'; then
  echo "[roihu] deeptools IS available as a module — you can skip tykky and 'module load deeptools' in pilot.slurm"
fi
module load tykky
mkdir -p "$WORK"
cp "$(dirname "$0")/env.yaml" "$WORK/env.yaml" 2>/dev/null || true
echo "[roihu] building container at $PREFIX (login node, ~10-20 min)..."
conda-containerize new --prefix "$PREFIX" "$WORK/env.yaml"
echo "[roihu] done. Use it via:  export PATH=$PREFIX/bin:\$PATH"
