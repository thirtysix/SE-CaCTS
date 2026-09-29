#!/bin/bash
# Start or resume the Stage B worker pool (worker.slurm), on Roihu from this directory:
#   PROJ=<project> ./launch.sh <manifest.tsv> <workers> [queue dir]
# First run writes <manifest>.pool.tsv: the rows not yet published, largest FASTQ first (worker.slurm's
# estimate), which the workers claim by index; later runs reuse it, so claims stay valid. Resuming needs no
# worker running: claims of rows neither published nor failed are dropped so they are taken again
# (RETRY_FAILED=1 retries failed rows too).
set -euo pipefail
: "${PROJ:?set PROJ}"; cd "$(dirname "$0")"
M="$1"; NW="$2"; P="${M%.tsv}.pool.tsv"; Q="${3:-$PWD/queue/$(basename "$M" .tsv)}"
W="/scratch/$PROJ/se-cacts/phase2"
[ -z "$(squeue -h -u "$USER" -n sB-worker)" ] || { echo "[launch] sB-worker jobs are still running"; exit 1; }
if [ ! -f "$P" ]; then
  { head -1 "$M"
    tail -n +2 "$M" | while IFS= read -r l; do s=${l%%$'\t'*}; [ -f "$W/out/${s: -2}/$s.done" ] || printf '%s\n' "$l"; done \
      | awk -F'\t' -v b="$(head -1 "$M" | tr '\t' '\n' | grep -nx bases | cut -d: -f1)" \
          '{printf "%.0f\t%s\n", 2 * $b + 58 * $8 * ($2 == "PAIRED" ? 2 : 1), $0}' | sort -t$'\t' -k1,1gr | cut -f2-
  } > "$P"
fi
mkdir -p "$Q"/{claim,reserve,failed,peak}; rm -f "$Q"/reserve/*; rmdir "$Q/lock" 2>/dev/null || true
[ "${RETRY_FAILED:-0}" = 1 ] && rm -f "$Q"/failed/*
dropped=0
for c in "$Q"/claim/*; do
  [ -e "$c" ] || continue; i=$(basename "$c"); s=$(awk -F'\t' -v n=$(( i + 2 )) 'NR==n{print $1}' "$P")
  [ -f "$W/out/${s: -2}/$s.done" ] || [ -e "$Q/failed/$i" ] || [ -e "$Q/failed/$i.post" ] || { rmdir "$c"; dropped=$((dropped+1)); }
done
echo "[launch] $(( $(wc -l < "$P") - 1 )) rows in $P, $(ls "$Q/claim" | wc -l) claimed ($dropped stale claims dropped)"
sbatch --parsable --account="$PROJ" --array=0-$(( NW - 1 )) \
  --export=ALL,PROJ="$PROJ",MANIFEST="$PWD/$P",QUEUE="$Q",MARGIN_GB="${MARGIN_GB:-80}"${FETCH_X:+,FETCH_X=$FETCH_X}${ALIGN_X:+,ALIGN_X=$ALIGN_X}${POST_X:+,POST_X=$POST_X} worker.slurm
