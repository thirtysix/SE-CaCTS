#!/bin/bash
#SBATCH --partition=small
#SBATCH --job-name=secacts-score
#SBATCH --time=03:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x-%j-%a.out
#SBATCH --error=logs/%x-%j-%a.err
# One scoring arm of the v2 atlas (see phase2/scripts/70_score_v2.sh for the arm definitions).
source /appl/profile/zz-csc-env.sh 2>/dev/null || true
set -euo pipefail
SB="${SB:?}"; W="${W:?}"; ARM="${ARM:?}"
export SECACTS_DATAROOT="$SB/dataroot" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=/projappl/${PROJ:?}/secacts_venv/bin/python
REPO="$SB/dataroot/002.AI_projects/SE-CaCTS"
OUT="$SB/out"; mkdir -p "$OUT"
DISC="ACH-000232,ACH-000338,ACH-000856,ACH-000957,ACH-000960,ACH-000997,ACH-001106,ACH-000768"
V2=(--signal "$W/results_v2/atlas.s3.se_signal.tsv" --catalog "$W/results_v2/atlas.s3.union_catalog.bed"
    --pull-set "$REPO/phase2/data/pull_set.v2.tsv" --norm none
    --fdr-method permutation --n-perm 1000 --dump-specific 0.10)
case "$ARM" in
  main)    args=("${V2[@]}" --out "$OUT/atlas.s3.perm" --cn-diagnostic --save-matrices) ;;
  nocn)    args=("${V2[@]}" --out "$OUT/atlas.s3.perm.nocn" --no-cn) ;;
  shuffle) args=("${V2[@]}" --out "$OUT/atlas.s3.perm.shuffle" --shuffle-labels 1
                 --levels OncotreeLineage,OncotreePrimaryDisease,OncotreeSubtype) ;;
  nodisc)  args=("${V2[@]}" --out "$OUT/atlas.s3.perm.nodisc" --exclude-keys "$DISC"
                 --levels OncotreeLineage,OncotreePrimaryDisease) ;;
  cmpfirst) args=("${V2[@]}" --pull-set "$REPO/phase2/data/pull_set.v2.cmpfirst.tsv" --out "$OUT/atlas.s3.perm.cmpfirst"
                 --cn-diagnostic --levels OncotreeLineage,OncotreePrimaryDisease) ;;
  anshuf)  args=(--signal "$W/results_v2/atlas.s3.se_signal.tsv" --catalog "$W/results_v2/atlas.s3.union_catalog.bed"
                 --pull-set "$REPO/phase2/data/pull_set.v2.tsv" --norm none --levels OncotreeLineage
                 --dump-specific 0.10 --shuffle-labels 1 --out "$OUT/atlas.s3.analytic.shuffle") ;;
  anreal)  args=(--signal "$W/results_v2/atlas.s3.se_signal.tsv" --catalog "$W/results_v2/atlas.s3.union_catalog.bed"
                 --pull-set "$REPO/phase2/data/pull_set.v2.tsv" --norm none --levels OncotreeLineage
                 --dump-specific 0.10 --out "$OUT/atlas.s3.analytic") ;;
  regress) args=(--signal "$W/results/atlas.s3.se_signal.tsv.gz" --catalog "$W/results/atlas.s3.union_catalog.bed.gz"
                 --pull-set "$REPO/phase2/data/pull_set.tsv" --norm none --levels OncotreeLineage
                 --out "$OUT/regress_v1") ;;
esac
echo "[arm] $ARM on $(hostname) $(date -Is)"
"$PY" "$REPO/phase2/score_pilot.py" "${args[@]}" > "$OUT/$ARM.stdout" 2> "$OUT/$ARM.stderr"
echo "[arm] $ARM done $(date -Is)"
