#!/bin/bash
#SBATCH --partition=small
#SBATCH --job-name=recall-cn
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/%x-%A-%a.out
#SBATCH --error=logs/%x-%A-%a.err
# Calling-time CN re-call from the retained <SRX>.enhancers.tsv (phase2/recall_cn.py; v3.1 item 13, GOTCHAS 59).
# Submit from $SB (score_v3/) as an array split by line:
#   sbatch -A <project> --array=0-15 --export=ALL,PROJ=<project>,NT=16 recall_cn.slurm
export CSC_ENV_INIT_NON_INTERACTIVE=yes; source /etc/profile.d/zz-csc-env.sh 2>/dev/null || true
set -euo pipefail
: "${PROJ:?}"; NT="${NT:-16}"; T="${SLURM_ARRAY_TASK_ID:-0}"
R=/scratch/$PROJ/se-cacts; W=$R/phase2; SB=$R/score_v3; REPO=$SB/dataroot/002.AI_projects/SE-CaCTS
export SECACTS_DATAROOT="$SB/dataroot" OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
PY=/projappl/$PROJ/secacts_venv/bin/python
mkdir -p "$SB/recall"
"$PY" "$REPO/phase2/recall_cn.py" --manifest "${MANIFEST:-$REPO/phase2/data/recall_manifest.v31.tsv}" \
  --out-dir "$W/out" --task "$T" --ntasks "$NT" --summary "$SB/recall/summary${MODE:+.$MODE}.$T.tsv" --mode "${MODE:-cn}" \
  --inferred-bins "$W/../cncal/bins" --blacklist "$W/../cncal/hg38-blacklist.v2.bed.gz" ${EXTRA:-}
