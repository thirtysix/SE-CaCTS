#!/bin/bash
#SBATCH --partition=small
#SBATCH --job-name=sB-worker
#SBATCH --time=24:00:00
#SBATCH --mem=40G
#SBATCH --cpus-per-task=20
#SBATCH --output=logs/%x-%A_%a.out
# One worker of the Stage B pool (launch.sh). It claims the largest unclaimed row whose disk need fits the free
# budget, fetches and aligns it on 16 cores, then posts it on 2 of the other 4 while it takes the next; post is
# single-threaded and takes ~1.5x a 16-thread alignment, so up to MAX_POSTS (2) run at once.
# Claims are mkdir, atomic on Lustre, so any number of workers share one queue and a relaunch resumes it. The
# budget is the live project quota minus MARGIN_GB, so a full disk makes workers wait instead of failing
# (GOTCHAS 110). 20 cores / 40 GB is the 2 GiB/core billing break-even.
set -euo pipefail
source ./common.sh
Q="${QUEUE:?set QUEUE}"; mkdir -p "$Q"/{claim,reserve,failed,peak}
export ALIGN_THREADS=16 FETCH_THREADS=6   # fasterq-dump output is byte-identical at 6 and 1 threads (dumptest.slurm)
LPID=$(lfs project -d "/scratch/$PROJ" | awk '{print $1}')
MARGIN_KB=$(( ${MARGIN_GB:-80} * 1048576 ))
# KB a row may hold, as multiples of its FASTQ size (2 B a base + 58 a read; `bases` is NCBI runinfo summed per
# SRX): FETCH_X while it downloads and unpacks (.sra, plain FASTQ, then gzip), ALIGN_X from then until the
# sort ends (gzipped reads + compressed BAM), POST_X until post ends. valtest.slurm measured 1.12 / 0.54 / 0.52.
fq_kb() { awk -F'\t' -v x="$1" 'NR==1{for (c = 1; c <= NF; c++) if ($c == "bases") b = c; next}
  {printf "%d\n", (2 * $b + 58 * $8 * ($2 == "PAIRED" ? 2 : 1)) * x / 1024}' "$MANIFEST"; }
mapfile -t NEED < <(fq_kb "${FETCH_X:-1.3}"); mapfile -t HOLD < <(fq_kb "${ALIGN_X:-0.65}"); mapfile -t POSTH < <(fq_kb "${POST_X:-0.6}")
N=${#NEED[@]}

lock()   { local n=0; until mkdir "$Q/lock" 2>/dev/null; do sleep 1; (( ++n % 300 )) || rmdir "$Q/lock" 2>/dev/null || true; done; }
unlock() { rmdir "$Q/lock" 2>/dev/null || true; }
# Free budget in KB: quota - margin - used - what reserved rows may still write beyond what they hold now
avail_kb() {
  local used quota k s h f pend=0
  read -r used quota < <(lfs quota -q -p "$LPID" /scratch | tr '\n' ' ' | awk '{gsub(/\*/, ""); print $2, $3}')
  for f in "$Q"/reserve/*; do
    [ -e "$f" ] || continue; read -r k s < "$f"
    h=$(du -sk "$S/work/$s" 2>/dev/null | cut -f1); h=${h:-0}; (( k > h )) && pend=$(( pend + k - h ))
  done
  echo $(( quota - MARGIN_KB - used - pend ))
}
# Print the first unclaimed row (largest first: launch.sh sorts the pool) that fits; fail when none is left
claim() {
  local a i left
  while :; do
    lock; a=$(avail_kb); left=0
    for (( i = 0; i < N; i++ )); do
      [ -d "$Q/claim/$i" ] && continue; left=1
      if (( NEED[i] <= a )); then
        mkdir "$Q/claim/$i"; echo "${NEED[i]} $(row "$i" | cut -f1)" > "$Q/reserve/$i"; unlock; echo "$i"; return 0
      fi
    done
    unlock; (( left )) || return 1
    echo "[worker $SLURM_ARRAY_TASK_ID] $(date +%T) nothing fits in $(( a / 1048576 )) GB free; waiting" >&2; sleep 60
  done
}

posts=()
while i=$(claim); do
  srx=$(row "$i" | cut -f1)
  echo "[worker $SLURM_ARRAY_TASK_ID] $(date +%T) row $i $srx, $(( NEED[i] / 1048576 )) GB reserved, $(( HOLD[i] / 1048576 )) GB once fetched"
  # the peak this row reaches on disk until its sort ends, to check FETCH_X / ALIGN_X against
  ( p=0; while sleep 20; do h=$(du -sk "$S/work/$srx" 2>/dev/null | cut -f1) || h=0; (( ${h:-0} > p )) && p=$h && echo "$p ${NEED[i]}" > "$Q/peak/$srx"; done ) & sampler=$!
  if step fetch_one "$i" && echo "${HOLD[i]} $srx" > "$Q/reserve/$i" && step align_one "$i"; then
    kill "$sampler" 2>/dev/null || true; echo "${POSTH[i]} $srx" > "$Q/reserve/$i"
    while (( ${#posts[@]} >= ${MAX_POSTS:-2} )); do wait "${posts[0]}" || true; posts=("${posts[@]:1}"); done
    { step post_one "$i" || { touch "$Q/failed/$i.post"; rm -rf "$S/work/$srx"; }; rm -f "$Q/reserve/$i"; } & posts+=($!)
  else
    kill "$sampler" 2>/dev/null || true
    touch "$Q/failed/$i"; rm -rf "$S/work/$srx"; rm -f "$Q/reserve/$i"
  fi
done
for p in "${posts[@]}"; do wait "$p" || true; done
echo "[worker $SLURM_ARRAY_TASK_ID] $(date +%T) queue empty"
