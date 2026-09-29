#!/bin/bash
#SBATCH --partition=small
#SBATCH --job-name=sB-valtest
#SBATCH --time=02:00:00
#SBATCH --mem=36G
#SBATCH --cpus-per-task=18
#SBATCH --output=logs/%x-%j.out
# Re-run one validated experiment (MANIFEST row ROW, mode=validate) through the current steps into
# validate_rerun/, compare every output with the validated copy in validate/, and log the work dir's disk peak
# per step (sampled every 2 s). For changes that should not move a byte: gzip'd FASTQ, compressed unsorted BAM,
# more align threads.
set -euo pipefail
source ./common.sh
export VALIDATE_DIR="$S/validate_rerun" ALIGN_THREADS="${ALIGN_THREADS:-16}" FETCH_THREADS="${FETCH_THREADS:-6}"
i="${ROW:-0}"; srx=$(row "$i" | cut -f1); d="$S/work/$srx"
rm -rf "$d" "${VALIDATE_DIR:?}/$srx"; mkdir -p "$VALIDATE_DIR"
peak() { local p=0 h; while sleep 2; do h=$(du -sk "$d" 2>/dev/null | cut -f1) || h=0; (( ${h:-0} > p )) && p=$h && echo "$p" > "$S/work/.peak.$srx.$1"; done; }
for s in fetch align post; do
  peak "$s" & sp=$!
  step "${s}_one" "$i"
  kill "$sp"; echo "[valtest] $srx $s peak_kb=$(cat "$S/work/.peak.$srx.$s" 2>/dev/null || echo 0)"
done
rm -f "$S/work/.peak.$srx".*
for f in "$S/validate/$srx"/*; do
  b=$(basename "$f"); g="$VALIDATE_DIR/$srx/$b"
  case "$b" in timing.*|*.stageB_qc.json) continue ;; esac
  case "$b" in *.gz) cmp -s <(zcat "$f") <(zcat "$g") ;; *) cmp -s "$f" "$g" ;; esac && echo "IDENTICAL $b" || echo "DIFFERENT $b"
done
echo "validated: $(cat "$S/validate/$srx/timing.align")  rerun (${ALIGN_THREADS} threads): $(cat "$VALIDATE_DIR/$srx/timing.align")"
