#!/usr/bin/env bash
# Phase 4 on the v2 (CN-expanded) atlas: the canonical permutation run plus the three runs that test it.
#
#   main     CN-corrected, permutation null (B=1000)      -> the calls
#   nocn     same, correction OFF                          -> the CN ablation's uncorrected arm
#   shuffle  labels shuffled across lines (seed 1)         -> calibration: a working FDR calls ~nothing
#   nodisc   the 7 cross-source-discordant CN lines + MDA-MB-231 excluded -> sensitivity
#
# ARMS="main nocn" bash phase2/scripts/70_score_v2.sh   runs a subset. Each arm is single-process with
# 2 BLAS threads, so three at once stay inside the laptop's 8-thread default (global CLAUDE.md).
set -euo pipefail
SECACTS="$(cd "$(dirname "$0")/../.." && pwd)"
. "$SECACTS/secacts_env.sh"
PY="${PY:-$HOME/miniconda3/envs/atac_hdac/bin/python}"
RES="${RES:-$SECACTS/phase2/results_v2}"
OUT="${OUT:-$SECACTS/phase2/scores_v2}"
ARMS="${ARMS:-main nocn shuffle}"
mkdir -p "$OUT"

DISCORDANT="ACH-000232,ACH-000338,ACH-000856,ACH-000957,ACH-000960,ACH-000997,ACH-001106,ACH-000768"
COMMON=(--signal "$RES/atlas.s3.se_signal.tsv.gz" --catalog "$RES/atlas.s3.union_catalog.bed.gz"
        --pull-set "$SECACTS/phase2/data/pull_set.v2.tsv" --norm none
        --fdr-method permutation --n-perm 1000 --dump-specific 0.10)

run() {   # run <arm> <prefix> <extra args...>
  echo "[70] $1 -> $2"
  OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2 \
    "$PY" "$SECACTS/phase2/score_pilot.py" "${COMMON[@]}" --out "$2" "${@:3}" \
    > "$2.run.out" 2> "$2.run.log"
  echo "[70] $1 done ($(tail -1 "$2.run.log"))"
}

pids=()
for arm in $ARMS; do
  case "$arm" in
    main)    run main    "$OUT/atlas.s3.perm"         --cn-diagnostic --save-matrices & ;;
    nocn)    run nocn    "$OUT/atlas.s3.perm.nocn"    --no-cn & ;;
    shuffle) run shuffle "$OUT/atlas.s3.perm.shuffle" --shuffle-labels 1 \
                         --levels OncotreeLineage,OncotreePrimaryDisease,OncotreeSubtype & ;;
    nodisc)  run nodisc  "$OUT/atlas.s3.perm.nodisc"  --exclude-keys "$DISCORDANT" \
                         --levels OncotreeLineage,OncotreePrimaryDisease & ;;
    *) echo "[70] unknown arm $arm" >&2; exit 1 ;;
  esac
  pids+=($!)
done
for p in "${pids[@]}"; do wait "$p"; done
echo "[70] all arms finished -> $OUT"
