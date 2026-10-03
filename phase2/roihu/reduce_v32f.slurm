#!/bin/bash
#SBATCH --partition=small
#SBATCH --job-name=reduce-v32f
#SBATCH --time=03:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
# v3.2 fused build (FINDINGS §57): v3.1f with measured CN for lines that were scored on input inference (CCMA WGS
# segments; subsampled WGS reads over a normal reference). Submit from $SB (score_v3/):
#   sbatch -A <project> --export=ALL,PROJ=<project> reduce_v32f.slurm
# 1. the samples whose line changed CN source (recall_redo.v32.tsv) are re-called in both modes (cn, fused), forced,
#    after their v3.1f calls are archived once in $W/backup_v31f_calls.tgz; every other sample keeps its v3.1f calls
# 2. baseline atlas: union of agnostic + fused calls, merge on any overlap, as v3.1f                -> $W/results_v32f
# 3. all experiments (variants 4-5) quantified on the SAME catalogue                              -> $W/results_v32fall
export CSC_ENV_INIT_NON_INTERACTIVE=yes; source /etc/profile.d/zz-csc-env.sh 2>/dev/null || true
set -euo pipefail
: "${PROJ:?}"
R=/scratch/$PROJ/se-cacts; W=$R/phase2; SB=$R/score_v3; REPO=$SB/dataroot/002.AI_projects/SE-CaCTS
export SECACTS_DATAROOT="$SB/dataroot" OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
PY=/projappl/$PROJ/secacts_venv/bin/python
D=$REPO/phase2/data
mkdir -p "$SB/recall"
redo="$D/recall_redo.v32.tsv"
# archive the v3.1f calls of every sample about to be re-called for the first time (later rounds add lines: each
# round's new files go to their own tarball, listed in backup_v31f_calls.list so nothing is archived twice)
touch "$W/backup_v31f_calls.list"
[ -s "$W/backup_v31f_calls.tgz" ] && [ ! -s "$W/backup_v31f_calls.list" ] && tar -tzf "$W/backup_v31f_calls.tgz" > "$W/backup_v31f_calls.list"
tail -n +2 "$redo" | cut -f1 | sort -u | while read -r s; do
  for x in fu.se.bed cn.se.bed cn.enhancers.tsv; do f="out/${s: -2}/$s.$x"; [ -e "$W/$f" ] && echo "$f"; done
done | grep -v -x -F -f "$W/backup_v31f_calls.list" > "$SB/recall/v32_backup.new" || true
if [ -s "$SB/recall/v32_backup.new" ]; then
  n=$(ls "$W"/backup_v31f_calls*.tgz 2>/dev/null | wc -l)
  tar -C "$W" -czf "$W/backup_v31f_calls.$n.tgz" -T "$SB/recall/v32_backup.new"
  cat "$SB/recall/v32_backup.new" >> "$W/backup_v31f_calls.list"
  echo "[v32f] archived $(wc -l < "$SB/recall/v32_backup.new") v3.1f call files -> backup_v31f_calls.$n.tgz"
fi
cn_args=(--inferred-bins "$W/../cncal/bins" --blacklist "$W/../cncal/hg38-blacklist.v2.bed.gz"
         --wgs-bins "$W/../wgs_cn/ratio" --ccma-segments "$D/cn_ccma_wgs.hg38.tsv.gz")
for mode in cn fused; do
  echo "[v32f] $mode re-call of $(($(wc -l < "$redo") - 1)) samples $(date -Is)"
  "$PY" "$REPO/phase2/recall_cn.py" --mode "$mode" --force --manifest "$redo" --out-dir "$W/out" \
    --summary "$SB/recall/summary.$mode.v32f.tsv" "${cn_args[@]}"
done
common=(--out-dir "$W/out" --grid "$W/data/grid.20.bed" --grid-label grid.20 --peak-dir "$W/data/bed20"
        --min-peaks 2000 --norm s3norm --s3-ref SRX16495452)
mkdir -p "$W/results_v32f" "$W/results_v32fall"
echo "[v32f] baseline reduce $(date -Is)"
"$PY" "$REPO/phase2/aggregate_v31.py" --samples "$D/pull_srx.v32.txt" "${common[@]}" --calls agnostic,fused \
      --merge any --out "$W/results_v32f/atlas.s3"
echo "[v32f] all-experiments reduce on the baseline catalogue $(date -Is)"
"$PY" "$REPO/phase2/aggregate_v31.py" --samples "$D/pull_srx.v32.all.txt" "${common[@]}" --fused \
      --catalog-in "$W/results_v32f/atlas.s3.union_catalog.bed" --out "$W/results_v32fall/atlas.s3"
echo "[v32f] done $(date -Is)"
