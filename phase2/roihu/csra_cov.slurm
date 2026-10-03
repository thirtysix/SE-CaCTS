#!/bin/bash
#SBATCH --partition=small
#SBATCH --job-name=csra-cov
#SBATCH --time=02:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x-%A-%a.out
#SBATCH --error=logs/%x-%A-%a.err
# Read counts per 5 kb from an SRA run submitted as a COORDINATE-SORTED alignment (v3.2, FINDINGS §59): the cSRA
# REFERENCE table lists, per 5 kb chunk of the submitter's reference, the ids of the primary alignments starting in
# it, and those ids are contiguous in sorted order, so a chunk's count is the next chunk's first id minus its own.
# Only the first id is kept (vdb-dump -M): no reads are downloaded or realigned. Lifting to hg38 and the pooled
# reference are local (phase1/scripts/27_wgs_cn.py csra). Manifest: name key role run (phase2/data/csra_cn_manifest.tsv).
#   sbatch -A <project> --array=0-<n-1> --export=ALL,PROJ=<project>,MANIFEST=<tsv> csra_cov.slurm
source /appl/profile/zz-csc-env.sh 2>/dev/null || true
set -euo pipefail
: "${PROJ:?}" "${MANIFEST:?}"
R=/scratch/$PROJ/se-cacts; S=$R/stageB; SB=${SB:-$R/score_v3}
SRA="apptainer exec -B /scratch $S/sif/sra-tools.sif"
OUT=$SB/wgs_cn/csra; mkdir -p "$OUT"
IFS=$'\t' read -r name key role run < <(awk -F'\t' -v i=$((SLURM_ARRAY_TASK_ID + 2)) 'NR==i' "$MANIFEST")
[ -s "$OUT/$name.ref.tsv.gz" ] && { echo "[csra] $name done, skipped"; exit 0; }
echo "[csra] $name ($key, $role) $run on $(hostname) $(date -Is)"
for k in 1 2 3; do
  $SRA vdb-dump -T REFERENCE -C NAME,SEQ_START,SEQ_LEN,PRIMARY_ALIGNMENT_IDS -f tab -M 60 "$run" \
    | awk -F'\t' '{split($4, a, ","); print $1 "\t" $2 "\t" $3 "\t" a[1]}' | gzip > "$OUT/$name.ref.tsv.gz.tmp" && break
  echo "[csra] attempt $k failed"; sleep $((60 * k))
done
mv "$OUT/$name.ref.tsv.gz.tmp" "$OUT/$name.ref.tsv.gz"
echo "[csra] $name: $(zcat "$OUT/$name.ref.tsv.gz" | wc -l) chunks $(date -Is)"
