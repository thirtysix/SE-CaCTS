#!/usr/bin/env bash
# Phase 2 — 44: rehearse the ARRAY BODY locally, against the live ChIP-Atlas endpoint. No Roihu, no BU.
#
#   bash 44_local_rehearse.sh [n_samples] [samples_per_task]
#
# Runs phase2/array.slurm itself (not a copy) with SLURM_ARRAY_TASK_ID stubbed, so the manifest chunking,
# resume, failure handling and output shape are exercised exactly as they will be on Roihu. Catches script
# bugs that would otherwise surface one task at a time, hours into a real array.
set -uo pipefail
N="${1:-4}"; SPT="${2:-2}"
HERE="$(cd "$(dirname "$0")" && pwd)"; SEC="$(cd "$HERE/../.." && pwd)"
# shellcheck source=/dev/null
[ -f "$SEC/secacts_env.sh" ] && . "$SEC/secacts_env.sh"
# Scratch for the rehearsal. Override with SCR=<dir>; defaults under the system temp dir.
SCR="${SCR:-${TMPDIR:-/tmp}/se-cacts-rehearse}"
W="$SCR/rehearse_array"; rm -rf "$W"; mkdir -p "$W"/{data,scripts,logs}
cp "$SEC/phase2/scripts/40_pull_one.sh" "$W/scripts/"
cp "$SEC/phase2/data/grid.20.bed" "$W/data/"
ls "$SEC/pilot/data/bigwigs"/*.bw | xargs -n1 basename | sed 's/\.bw$//' | head -"$N" > "$W/data/pull_srx.txt"
# Per-SRX bed20 peak dir. No default location exists — it is fetched by 22_fetch_peaks.sh —
# so require it explicitly rather than guessing at a path that may be someone else's.
BED20="${BED20:-$SEC/phase2/data/bed20}"
[[ -d "$BED20" ]] || { echo "[44] FATAL: no bed20 peak dir at '$BED20'; set BED20=<dir of SRX.20.bed> (see 22_fetch_peaks.sh)" >&2; exit 1; }
ln -sfn "$BED20" "$W/data/bed20"
TASKS=$(( (N + SPT - 1) / SPT ))
echo "[44] $N samples, $SPT/task -> $TASKS tasks, live endpoint"
for (( t=0; t<TASKS; t++ )); do
  echo "--- task $t ---"
  SLURM_ARRAY_TASK_ID=$t PY="$HOME/miniconda3/envs/atac_hdac/bin/python" \
  WORK="$W" CNROSE="$SEC/cnrose" OUT="$W/out" PEAK_DIR="$W/data/bed20" \
  GRIDS="$W/data/grid.20.bed" MANIFEST="$W/data/pull_srx.txt" \
  SAMPLES_PER_TASK="$SPT" BINS=100 LIMIT_RATE=20m \
  bash "$SEC/phase2/array.slurm" 2>&1 | grep -E '^\[array\]|^\[40\]|DONE'
done
echo "[44] --- results ---"
echo "  done=$(find "$W/out" -name '*.done' | wc -l)  failed=$(find "$W/out" -name '*.failed' | wc -l)"
echo "  shards: $(ls "$W/out" | tr '\n' ' ')"
echo "  per-sample files: $(find "$W/out" -name '*.qc.json' | wc -l) qc, $(find "$W/out" -name '*.f32' | wc -l) f32, $(find "$W/out" -name '*.f16.gz' | wc -l) archives"
du -sh "$W/out"
echo "[44] re-running (everything should SKIP):"
SLURM_ARRAY_TASK_ID=0 PY="$HOME/miniconda3/envs/atac_hdac/bin/python" WORK="$W" CNROSE="$SEC/cnrose" \
  OUT="$W/out" PEAK_DIR="$W/data/bed20" GRIDS="$W/data/grid.20.bed" MANIFEST="$W/data/pull_srx.txt" \
  SAMPLES_PER_TASK="$SPT" bash "$SEC/phase2/array.slurm" 2>&1 | grep -E 'DONE'
