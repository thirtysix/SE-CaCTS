#!/bin/bash
#SBATCH --partition=small
#SBATCH --job-name=sB-align
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/%x-%A_%a.out
# Bowtie2 2.2.2 defaults + --no-unal (pairs as pairs), samtools 0.1.19 sort: the multi-threaded part only.
# 8 cores / 16 GB sits at the 2 GiB/core billing break-even. Everything single-threaded is in post.slurm.
set -euo pipefail
source ./common.sh
# Roihu caps a user at 200 queued / 100 running jobs on `small`, so each array task takes SAMPLES_PER_TASK
# consecutive manifest rows, one after another; a failed row never stops the rest of the task.
process() {
  r=$(row "$1"); srx=$(cut -f1 <<<"$r"); layout=$(cut -f2 <<<"$r")
  d="$S/work/$srx"; cd "$d"
  [ -f .aligned ] && { echo "[align] $srx done"; exit 0; }
  [ -f .fetched ] || { echo "[align] $srx not fetched"; exit 1; }
  z=.gz; [ -s r1.fq.gz ] || z=""                        # ENA gives .fq.gz; the SRA route gives plain .fq
  if [ "$layout" = PAIRED ] && [ -s "r2.fq$z" ]; then in=(-1 "r1.fq$z" -2 "r2.fq$z"); else in=(-U "r1.fq$z"); fi
  t0=$(date +%s)
  "$BT2" -p "$SLURM_CPUS_PER_TASK" -t --no-unal -x "$IDX" -q "${in[@]}" 2> bowtieReport.txt | $SAMT view -S -u - > unsrt.bam
  t1=$(date +%s)
  $SAMT sort -@ "$SLURM_CPUS_PER_TASK" -m 800M unsrt.bam dup
  rm -f unsrt.bam r1.fq.gz r2.fq.gz r1.fq r2.fq
  awk '/reads; of these/{print $1}' bowtieReport.txt > nspots
  echo "align_s=$((t1-t0)) sort_s=$(( $(date +%s) - t1 ))" > timing.align
  cat bowtieReport.txt timing.align; touch .aligned
}
K="${SAMPLES_PER_TASK:-1}"; N=$(( $(wc -l < "$MANIFEST") - 1 ))
for (( i = SLURM_ARRAY_TASK_ID * K; i < (SLURM_ARRAY_TASK_ID + 1) * K && i < N; i++ )); do
  ( process "$i" ) || echo "[task $SLURM_ARRAY_TASK_ID] row $i failed"
done
