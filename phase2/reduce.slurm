#!/bin/bash
# NOTE: --account is NOT hardcoded. secacts_env.sh exports SBATCH_ACCOUNT from
#       SECACTS_CSC_PROJECT (.env), which sbatch honours. Override per-run with
#       `sbatch --account=...`.
#SBATCH --partition=small
#SBATCH --job-name=secacts-reduce
#SBATCH --time=02:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#
# Phase 2 — the barrier reduce. Runs ONCE after the array (PULL_DESIGN.md §5.1 shard-then-reduce), reading
# every per-sample output and writing the shared products a single time. Submit chained to the array:
#
#   aid=$(sbatch --parsable --array=0-182%16 array.slurm)
#   sbatch --dependency=afterok:$aid reduce.slurm
#
# WHY 48G: the fine matrix is n_grid x n_samples float32 — grid.20 at 2,917 samples is 504,855 x 2,917 x 4 B
# = 5.9 GB resident, and `--norm s3norm` holds a second copy plus the strata masks. 48G leaves room for the
# grid.10 reduce too. This is the one step that is NOT embarrassingly parallel.
#
# FIRST: it reports how many samples failed. A reduce over a partially-complete pull is silently wrong, so
# it refuses to run if the failure count exceeds MAX_FAILED unless FORCE=1.

source /appl/profile/zz-csc-env.sh
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
VENV="${VENV:-/projappl/$PROJ/secacts_venv}"
OUT="${OUT:-$WORK/out}"
MANIFEST="${MANIFEST:-$WORK/data/pull_srx.txt}"
GRID="${GRID:-$WORK/data/grid.20.bed}"
GRID_LABEL="${GRID_LABEL:-grid.20}"
PEAK_DIR="${PEAK_DIR:-$WORK/data/bed20}"
MIN_PEAKS="${MIN_PEAKS:-2000}"
MAX_FAILED="${MAX_FAILED:-50}"
RESULTS="${RESULTS:-$WORK/results}"

if [[ -x "$VENV/bin/python" ]]; then PY="$VENV/bin/python"
else echo "[reduce] FATAL: no python with pyBigWig at $VENV (build it: roihu/00_build_env.sh)" >&2; exit 1; fi

mkdir -p "$RESULTS"
n_expect=$(grep -cv '^[[:space:]]*$' "$MANIFEST")
n_done=$(find "$OUT" -name '*.done' | wc -l)
n_failed=$(find "$OUT" -name '*.failed' | wc -l)
echo "[reduce] manifest=$n_expect  done=$n_done  failed=$n_failed"
find "$OUT" -name '*.failed' -printf '%f\n' | sed 's/\.failed$//' | sort > "$RESULTS/failed_srx.txt" 2>/dev/null || true

if (( n_failed > MAX_FAILED )) && [[ "${FORCE:-0}" != "1" ]]; then
  echo "[reduce] REFUSING: $n_failed failures > MAX_FAILED=$MAX_FAILED. Investigate $RESULTS/failed_srx.txt," >&2
  echo "[reduce] re-run the array to retry them (it is resumable), or set FORCE=1 to reduce anyway." >&2
  exit 1
fi
if (( n_done < n_expect )) && [[ "${FORCE:-0}" != "1" ]]; then
  echo "[reduce] REFUSING: only $n_done/$n_expect samples complete. Re-run the array (resumable) or FORCE=1." >&2
  exit 1
fi

# the pull shards outputs into out/<last2>/ — give aggregate.py one flat view without copying
FLAT="$WORK/.flat"
rm -rf "$FLAT"; mkdir -p "$FLAT"
find "$OUT" -mindepth 2 -maxdepth 2 -type f \
     \( -name '*.se.bed' -o -name "*.$GRID_LABEL.f32" -o -name '*.signal.tsv' \) \
     -exec ln -sf {} "$FLAT/" \;
echo "[reduce] linked $(ls "$FLAT" | wc -l) per-sample files"

run () {   # run <out-prefix> <extra args...>
  echo "[reduce] -> $1"
  "$PY" "$WORK/aggregate.py" --se-dir "$FLAT" --signal-dir "$FLAT" --grid "$GRID" \
        --grid-label "$GRID_LABEL" --peak-dir "$PEAK_DIR" --out "$1" "${@:2}"
}

# Four catalogs: {agnostic, CN-corrected} x {unnormalized, S3norm}. S3norm REQUIRES the QC gate — ungated it
# is worse than quantile (DESIGN.md §normalization / PULL_DESIGN.md §8.5).
run "$RESULTS/atlas"
run "$RESULTS/atlas.cn"     --se-suffix .cn
run "$RESULTS/atlas.s3"     --norm s3norm --min-peaks "$MIN_PEAKS"
run "$RESULTS/atlas.s3.cn"  --se-suffix .cn --norm s3norm --min-peaks "$MIN_PEAKS"

echo "[reduce] done -> $RESULTS/atlas{,.cn,.s3,.s3.cn}.*"
echo "[reduce] NOTE: $n_failed samples failed and are NOT in the atlas; see $RESULTS/failed_srx.txt"
