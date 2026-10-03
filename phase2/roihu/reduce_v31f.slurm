#!/bin/bash
#SBATCH --partition=small
#SBATCH --job-name=reduce-v31f
#SBATCH --time=03:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
# v3.1 fused build (FINDINGS §54-55, user 2026-10-03). Submit from $SB (score_v3/):
#   sbatch -A <project> --export=ALL,PROJ=<project> reduce_v31f.slurm
# 1. fused re-call of every v3.1 sample (baseline + treated), forced, with the CN source the pull set scores it on
#    (the first fused run read recall_manifest.v31.tsv, whose source differs for 11 samples)
# 2. baseline atlas: union of agnostic + fused calls, merge on any overlap; presence (agnostic, corrected, fused)
#    and fused label codes by overlap; S3norm pinned to the v1-v3.1 reference              -> $W/results_v31f
# 3. all experiments (variants 4-5) quantified on the SAME catalogue                       -> $W/results_v31fall
export CSC_ENV_INIT_NON_INTERACTIVE=yes; source /etc/profile.d/zz-csc-env.sh 2>/dev/null || true
set -euo pipefail
: "${PROJ:?}"
R=/scratch/$PROJ/se-cacts; W=$R/phase2; SB=$R/score_v3; REPO=$SB/dataroot/002.AI_projects/SE-CaCTS
export SECACTS_DATAROOT="$SB/dataroot" OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
PY=/projappl/$PROJ/secacts_venv/bin/python
D=$REPO/phase2/data
mkdir -p "$SB/recall"
cut -f1,2,3,7,8,9,10 "$D/pull_set.v31.all.tsv" > "$SB/recall/manifest.fused.tsv"   # srx cell cvcl key cn_provider cn_cvcl source
echo "[v31f] fused re-call of $(($(wc -l < "$SB/recall/manifest.fused.tsv") - 1)) samples $(date -Is)"
"$PY" "$REPO/phase2/recall_cn.py" --mode fused --force --manifest "$SB/recall/manifest.fused.tsv" --out-dir "$W/out" \
  --summary "$SB/recall/summary.fused_v31f.tsv" --inferred-bins "$W/../cncal/bins" \
  --blacklist "$W/../cncal/hg38-blacklist.v2.bed.gz"
common=(--out-dir "$W/out" --grid "$W/data/grid.20.bed" --grid-label grid.20 --peak-dir "$W/data/bed20"
        --min-peaks 2000 --norm s3norm --s3-ref SRX16495452)
mkdir -p "$W/results_v31f" "$W/results_v31fall"
echo "[v31f] baseline reduce $(date -Is)"
"$PY" "$REPO/phase2/aggregate_v31.py" --samples "$D/pull_srx.v31.txt" "${common[@]}" --calls agnostic,fused \
      --merge any --out "$W/results_v31f/atlas.s3"
echo "[v31f] all-experiments reduce on the baseline catalogue $(date -Is)"
"$PY" "$REPO/phase2/aggregate_v31.py" --samples "$D/pull_srx.v31.all.txt" "${common[@]}" --fused \
      --catalog-in "$W/results_v31f/atlas.s3.union_catalog.bed" --out "$W/results_v31fall/atlas.s3"
echo "[v31f] done $(date -Is)"
