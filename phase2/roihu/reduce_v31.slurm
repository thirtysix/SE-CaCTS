#!/bin/bash
#SBATCH --partition=small
#SBATCH --job-name=reduce-v31
#SBATCH --time=03:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
# v3.1 reduce (ROADMAP "v3.1 scope, decided 2026-10-02"). Submit from $SB (score_v3/):
#   sbatch -A <project> --export=ALL,PROJ=<project> reduce_v31.slurm
# 1. calling-time CN re-call for pull-set samples that have none yet (recall_cn.py skips existing outputs), and
#    a forced re-call for lines whose CN source moved to DepMap WGS (25_pull_set_v31.py: recall_redo.v31.tsv)
# 2. baseline atlas: union of agnostic + corrected calls, merge on any overlap, presence by overlap,
#    S3norm pinned to the v1/v2/v3 reference                                          -> $W/results_v31
# 3. all experiments (variants 4-5) quantified on the SAME catalogue                  -> $W/results_v31all
export CSC_ENV_INIT_NON_INTERACTIVE=yes; source /etc/profile.d/zz-csc-env.sh 2>/dev/null || true
set -euo pipefail
: "${PROJ:?}"
R=/scratch/$PROJ/se-cacts; W=$R/phase2; SB=$R/score_v3; REPO=$SB/dataroot/002.AI_projects/SE-CaCTS
export SECACTS_DATAROOT="$SB/dataroot" OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2
PY=/projappl/$PROJ/secacts_venv/bin/python
D=$REPO/phase2/data
echo "[v31] recall for samples without .cn.se.bed $(date -Is)"
"$PY" - "$D/pull_set.v31.all.tsv" "$W/out" "$SB/recall/manifest.missing.tsv" <<'EOF'
import os, sys, pandas as pd
ps, out, dst = sys.argv[1:4]
p = pd.read_csv(ps, sep="\t")
miss = p[[not os.path.exists(os.path.join(out, s[-2:], s + ".cn.se.bed")) for s in p.srx]]
miss[["srx", "key", "cn_provider", "cvcl", "cn_cvcl", "cell", "source"]].to_csv(dst, sep="\t", index=False)
print(f"[v31] {len(miss)} samples need a calling-time CN re-call")
EOF
if [ "$(wc -l < "$SB/recall/manifest.missing.tsv")" -gt 1 ]; then
  "$PY" "$REPO/phase2/recall_cn.py" --manifest "$SB/recall/manifest.missing.tsv" --out-dir "$W/out" \
    --summary "$SB/recall/summary.missing.tsv" --inferred-bins "$W/../cncal/bins" \
    --blacklist "$W/../cncal/hg38-blacklist.v2.bed.gz"
fi
echo "[v31] re-call with the new CN source for lines moved to DepMap WGS $(date -Is)"
"$PY" "$REPO/phase2/recall_cn.py" --manifest "$D/recall_redo.v31.tsv" --out-dir "$W/out" --force \
  --summary "$SB/recall/summary.redo_wgs.tsv" --inferred-bins "$W/../cncal/bins" \
  --blacklist "$W/../cncal/hg38-blacklist.v2.bed.gz"
common=(--out-dir "$W/out" --grid "$W/data/grid.20.bed" --grid-label grid.20 --peak-dir "$W/data/bed20"
        --min-peaks 2000 --norm s3norm --s3-ref SRX16495452)
mkdir -p "$W/results_v31" "$W/results_v31all"
echo "[v31] baseline reduce $(date -Is)"
"$PY" "$W/aggregate_v31.py" --samples "$D/pull_srx.v31.txt" "${common[@]}" --calls agnostic,cn --merge any \
      --out "$W/results_v31/atlas.s3"
echo "[v31] all-experiments reduce on the baseline catalogue $(date -Is)"
"$PY" "$W/aggregate_v31.py" --samples "$D/pull_srx.v31.all.txt" "${common[@]}" \
      --catalog-in "$W/results_v31/atlas.s3.union_catalog.bed" --out "$W/results_v31all/atlas.s3"
echo "[v31] done $(date -Is)"
