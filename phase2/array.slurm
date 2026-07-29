#!/bin/bash
# NOTE: --account is NOT hardcoded. secacts_env.sh exports SBATCH_ACCOUNT from
#       SECACTS_CSC_PROJECT (.env), which sbatch honours. Override per-run with
#       `sbatch --account=...`.
#SBATCH --partition=small
#SBATCH --job-name=secacts-pull
#SBATCH --time=04:00:00
#SBATCH --mem=2500M
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.err
#
# Phase 2 — the heavy per-SRX pull (PULL_DESIGN.md §5). Each array task takes a CHUNK of the manifest and,
# per sample: download the bigWig -> one cnrose pass (SE calls + every grid + genome-wide archive + QC) ->
# provenance -> delete the bigWig (§8.5 retention set). Compute nodes have direct internet (verified), so
# this is a plain array — no login-node feeder.
#
#   sbatch --array=0-182%16 array.slurm        # 183 tasks x 16 samples, 16 concurrent (courtesy throttle)
#
# WHY 1 CPU / 2500M: cnrose is single-threaded. Process RSS is a stable ~1.4 GiB and barely grows with file
# size (measured 1.37-1.44 GiB across chunk sizes on the 1.22 GB worst-case manifest file). SLURM's cgroup
# MaxRSS additionally counts the bigWig's page cache (reclaimable) -> ~1.86 GiB on that file via the real
# download->cnrose path, which COMPLETED cleanly at --mem=2G. 2G is the billing break-even (0.75/core ==
# 0.375x2GiB); we take 2500M (+~25% BU, ~+9 BU total on the run) purely as headroom for the unknown
# large-file tail of the full 2,917 set, so an unattended run never OOMs. Do NOT raise further without
# re-measuring: above 2 GiB every extra GiB is billed. Chunk size does NOT reduce this (page cache dominates).
#
# WHY --time=04:00:00: ~60 s/sample (~30 s download + 15-42 s compute) x SAMPLES_PER_TASK=16 = ~16 min.
# Four hours is ~15x headroom for a slow endpoint; `small` allows 3 days if a chunk needs it.
#
# RESUMABLE: re-submitting the same array skips every sample with a `.done` marker. Nothing depends on a
# live laptop or a valid SSH certificate — submitted jobs survive cert expiry (README: re-sign every 24 h).

# before set -u (the CSC env script is not -u clean); absent off-Roihu so the array body stays
# locally testable — see scripts/44_local_rehearse.sh
source /appl/profile/zz-csc-env.sh 2>/dev/null || true
set -uo pipefail
# Machine-local config: walk up to the repo root and source .env via secacts_env.sh.
# Depth-independent so this works from a staged copy on the cluster too; if the file is
# absent (partial stage), a real SECACTS_CSC_PROJECT in the environment still satisfies the
# `:?` checks below.
_d="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [ "$_d" != "/" ] && [ ! -f "$_d/secacts_env.sh" ]; do _d="$(dirname "$_d")"; done
[ -f "$_d/secacts_env.sh" ] && . "$_d/secacts_env.sh"

PROJ="${PROJ:-${SECACTS_CSC_PROJECT:?set SECACTS_CSC_PROJECT in .env, or pass PROJ=}}"
WORK="${WORK:-/scratch/$PROJ/se-cacts/phase2}"
CNROSE="${CNROSE:-$WORK/cnrose}"
VENV="${VENV:-/projappl/$PROJ/secacts_venv}"     # built on the login node by roihu/00_build_env.sh

MANIFEST="${MANIFEST:-$WORK/data/pull_srx.txt}"
OUT="${OUT:-$WORK/out}"
PEAK_DIR="${PEAK_DIR:-$WORK/data/bed20}"
GRIDS="${GRIDS:-$WORK/data/grid.20.bed $WORK/data/grid.10.bed}"
SAMPLES_PER_TASK="${SAMPLES_PER_TASK:-16}"
BINS="${BINS:-100}"
LIMIT_RATE="${LIMIT_RATE:-}"

mkdir -p "$OUT" logs

# Toolchain: cnrose needs only pyBigWig + numpy + scipy (no R, no deepTools). Roihu has NO python-data
# module and NO tykky, so we use a plain venv built on the login node (roihu/00_build_env.sh); Roihu CPU
# nodes are x86_64, so pyBigWig installs from a manylinux wheel with no compilation.
if [[ -n "${PY:-}" ]]; then
  :                                   # caller supplied an interpreter (local rehearsal)
elif [[ -x "$VENV/bin/python" ]]; then
  PY="$VENV/bin/python"
else
  echo "[array] FATAL: no python with pyBigWig at $VENV (build it: roihu/00_build_env.sh)" >&2; exit 1
fi
echo "[array] task ${SLURM_ARRAY_TASK_ID:-?} on $(hostname) | python=$($PY -c 'import sys;print(sys.version.split()[0])')"

# Fail fast on CONFIGURATION errors. Without this, a missing peak dir or manifest makes every task march
# through its chunk writing .failed markers — cheap, but it looks like a data problem and buries the real
# cause. Per-SAMPLE problems still write .failed and continue; only whole-run misconfiguration aborts.
[[ -s "$MANIFEST" ]] || { echo "[array] FATAL: manifest missing/empty: $MANIFEST" >&2; exit 1; }
if [[ ! -d "$PEAK_DIR" ]] || [[ -z "$(ls -A "$PEAK_DIR" 2>/dev/null | head -1)" ]]; then
  echo "[array] FATAL: peak dir missing or empty: $PEAK_DIR (fetch bed20 on Roihu: 21_grid_from_persrx.sh)" >&2
  exit 1
fi
for g in $GRIDS; do
  [[ -s "$g" ]] || { echo "[array] FATAL: grid missing: $g" >&2; exit 1; }
done

mapfile -t ALL < <(grep -v '^[[:space:]]*$' "$MANIFEST")
N=${#ALL[@]}
START=$(( SLURM_ARRAY_TASK_ID * SAMPLES_PER_TASK ))
END=$(( START + SAMPLES_PER_TASK ))
(( END > N )) && END=$N
echo "[array] manifest=$N samples; this task handles [$START,$END)"

ok=0; skip=0; fail=0
for (( i=START; i<END; i++ )); do
  srx="${ALL[$i]}"
  shard="${srx: -2}"
  if [[ -f "$OUT/$shard/$srx.done" ]]; then skip=$((skip+1)); continue; fi
  OUT="$OUT" GRIDS="$GRIDS" PEAK_DIR="$PEAK_DIR" CNROSE="$CNROSE" PY="$PY" \
    BINS="$BINS" LIMIT_RATE="$LIMIT_RATE" \
    bash "$WORK/scripts/40_pull_one.sh" "$srx"
  if [[ -f "$OUT/$shard/$srx.done" ]]; then ok=$((ok+1)); else fail=$((fail+1)); fi
done

# Report explicitly — never let a partial run look complete (PULL_DESIGN.md §5)
echo "[array] task ${SLURM_ARRAY_TASK_ID:-?} DONE: ok=$ok skipped=$skip failed=$fail (of $((END-START)))"
