#!/usr/bin/env bash
# After v3_chain.sh's arms finish on Roihu: fetch the v3 atlas and scores, build the tables the figures read,
# then the poster panels into poster/figures_v3 (v2's stay in poster/figures).
#   bash phase2/scripts/72_finish_v3.sh            all steps
#   STEPS="fetch derive" bash ...                  a subset (fetch derive figures)
# One single-threaded process at a time: the laptop throttles on 2-3 (GOTCHA 103).
set -euo pipefail
SECACTS="$(cd "$(dirname "$0")/../.." && pwd)"
. "$SECACTS/secacts_env.sh"
PROJ="${SECACTS_CSC_PROJECT:?}"; R="/scratch/$PROJ/se-cacts"
PY="${PY:-$HOME/miniconda3/envs/atac_hdac/bin/python}"
RES="$SECACTS/phase2/results_v3"; SC="$SECACTS/phase2/scores_v3"
STEPS="${STEPS:-fetch derive figures}"
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
mkdir -p "$RES" "$SC"

if [[ " $STEPS " == *" fetch "* ]]; then
  ssh roihu "cd $R/phase2/results_v3 && for f in *.tsv *.bed; do [ -f \$f.gz ] || pigz -p 4 -k \$f; done"
  rsync -a --include='*.gz' --include='failed_srx.txt' --exclude='*' "roihu:$R/phase2/results_v3/" "$RES/"
  rsync -a --exclude='*.stdout' "roihu:$R/score_v3/out/" "$SC/"
  rsync -a "roihu:$R/score_v3/dataroot/002.AI_projects/SE-CaCTS/phase2/data/pull_set.v3.tsv" "$SECACTS/phase2/data/"
  rsync -a "roihu:$R/phase2/v3check/stageB_candidate_check.tsv" "$SECACTS/phase2/analysis/out/"
  echo "[72] fetched -> $RES, $SC"
fi
if [[ " $STEPS " == *" derive "* ]]; then
  "$PY" "$SECACTS/phase2/analysis/cn_ablation_calls.py" --corrected "$SC/atlas.s3.perm" \
        --uncorrected "$SC/atlas.s3.perm.nocn" --catalog "$RES/atlas.s3.union_catalog.bed.gz" \
        --out "$SC/atlas.s3.perm.cn_ablation_calls.tsv"
  "$PY" "$SECACTS/phase2/analysis/concordance_bridge2.py" --scores "$SC/atlas.s3.perm" \
        --signal "$RES/atlas.s3.se_signal.tsv.gz" --catalog "$RES/atlas.s3.union_catalog.bed.gz" \
        --pull-set "$SECACTS/phase2/data/pull_set.v3.tsv" --levels OncotreeLineage,OncotreePrimaryDisease \
        --out "$SC/atlas.s3.perm.concordance2"
  echo "[72] derived tables -> $SC"
fi
if [[ " $STEPS " == *" figures "* ]]; then
  POSTER_RES="$RES" POSTER_SC="$SC" POSTER_OUT="$SECACTS/poster/figures_v3" \
    POSTER_PS="$SECACTS/phase2/data/pull_set.v3.tsv" POSTER_PREV_SC="$SECACTS/phase2/scores_v2" \
    POSTER_LABEL=v3 POSTER_PREV_LABEL=v2 \
    "$PY" "$SECACTS/phase2/figures/poster_figures.py"
  echo "[72] figures -> poster/figures_v3"
fi
