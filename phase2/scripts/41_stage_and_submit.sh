#!/usr/bin/env bash
# Phase 2 — 41: stage the pull to Roihu scratch and submit it. Run FROM THE LAPTOP.
#
#   bash 41_stage_and_submit.sh              # stage + dry-run (prints the sbatch lines, submits nothing)
#   bash 41_stage_and_submit.sh --go         # stage + actually submit array + chained reduce
#
# Requires a working `ssh roihu` (certificate re-signed within 24 h — see pilot/roihu/README.md).
# The per-SRX bed20 peak BEDs are fetched ON ROIHU by 21_grid_from_persrx.sh, not shipped from here.
set -euo pipefail
# Machine-local config: walk up to the repo root and source .env via secacts_env.sh.
# Depth-independent so this works from a staged copy on the cluster too; if the file is
# absent (partial stage), a real SECACTS_CSC_PROJECT in the environment still satisfies the
# `:?` checks below.
_d="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [ "$_d" != "/" ] && [ ! -f "$_d/secacts_env.sh" ]; do _d="$(dirname "$_d")"; done
[ -f "$_d/secacts_env.sh" ] && . "$_d/secacts_env.sh"

HERE="$(cd "$(dirname "$0")" && pwd)"
SECACTS="$(cd "$HERE/../.." && pwd)"
PROJ="${PROJ:-${SECACTS_CSC_PROJECT:?set SECACTS_CSC_PROJECT in .env, or pass PROJ=}}"
REMOTE="${REMOTE:-roihu}"
WORK="/scratch/$PROJ/se-cacts/phase2"
CONC="${CONC:-16}"                 # array concurrency = courtesy throttle AND the transient-disk cap (§6)
SPT="${SPT:-16}"                   # samples per array task
GO=0; [[ "${1:-}" == "--go" ]] && GO=1

MANIFEST="$SECACTS/phase2/data/pull_srx.txt"
N=$(grep -cv '^[[:space:]]*$' "$MANIFEST")
TASKS=$(( (N + SPT - 1) / SPT ))

cat <<EOF
[41] plan
  samples            : $N
  samples/task       : $SPT   -> array 0-$((TASKS-1))
  concurrency        : %$CONC
  transient disk     : ~$(( CONC * 322 / 1024 )) GiB steady-state (concurrent bigWigs)
  est. compute       : ~$(( N * 60 / 3600 )) core-hours   (~60 s/sample: ~30 s download + 15-42 s cnrose)
  est. wall-clock    : ~$(( N * 60 / CONC / 60 )) min at %$CONC (network-bound; see PULL_DESIGN.md §6)
  durable output     : ~60 GB  (PULL_DESIGN.md §8.5.4)
EOF

echo "[41] staging code + inputs to $REMOTE:$WORK"
ssh "$REMOTE" "mkdir -p $WORK/{data,scripts,logs,results}"
rsync -az --delete "$SECACTS/cnrose/" "$REMOTE:$WORK/cnrose/"
# record which cnrose commit produced this run (the staged tree has no .git) -> provenance.json
GITSHA=$(git -C "$SECACTS" rev-parse --short HEAD 2>/dev/null || echo unknown)
ssh "$REMOTE" "echo $GITSHA > $WORK/cnrose.gitsha"
echo "[41] cnrose git = $GITSHA"
rsync -az "$SECACTS/phase2/aggregate.py" "$SECACTS/phase2/s3norm.py" \
          "$SECACTS/phase2/specificity.py" "$SECACTS/phase2/array.slurm" \
          "$SECACTS/phase2/reduce.slurm" "$REMOTE:$WORK/"
rsync -az "$SECACTS/phase2/scripts/40_pull_one.sh" "$REMOTE:$WORK/scripts/"
rsync -az "$SECACTS/phase2/data/grid.20.bed" "$SECACTS/phase2/data/pull_srx.txt" \
          "$SECACTS/phase2/data/pull_set.tsv" "$REMOTE:$WORK/data/"
# grid.10 is gitignored locally (regenerable); ship it if present, else build it on Roihu first
[[ -f "$SECACTS/phase2/data/grid.10.bed" ]] \
  && rsync -az "$SECACTS/phase2/data/grid.10.bed" "$REMOTE:$WORK/data/" \
  || echo "[41] WARN grid.10.bed absent locally — build it on Roihu (21_grid_from_persrx.sh) before submitting"

echo "[41] verifying remote inputs"
ssh "$REMOTE" "cd $WORK && ls -1 data/ && echo '--- bed20 peaks:' && ls data/bed20 2>/dev/null | wc -l"

SUB1="sbatch --parsable --array=0-$((TASKS-1))%$CONC array.slurm"
SUB2="sbatch --dependency=afterok:\$aid reduce.slurm"
if (( GO )); then
  echo "[41] submitting"
  ssh "$REMOTE" "cd $WORK && aid=\$($SUB1) && echo array=\$aid && $SUB2 && squeue -u \$USER"
else
  echo "[41] DRY RUN — would submit (re-run with --go):"
  echo "     cd $WORK && aid=\$($SUB1) && $SUB2"
fi
