#!/bin/bash
#SBATCH --partition=small
#SBATCH --job-name=sB-post
#SBATCH --time=12:00:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/%x-%A_%a.out
# Post step only: the single-threaded tail (rmdup, bigWig, MACS2, cnrose) on 2 cores.
set -euo pipefail
source ./common.sh
# One step per array task (post_one in common.sh). Roihu caps a user at 200 queued / 100 running jobs on `small`,
# so each task takes SAMPLES_PER_TASK consecutive manifest rows; a failed row never stops the rest of the task.
K="${SAMPLES_PER_TASK:-1}"; N=$(( $(wc -l < "$MANIFEST") - 1 ))
for (( i = SLURM_ARRAY_TASK_ID * K; i < (SLURM_ARRAY_TASK_ID + 1) * K && i < N; i++ )); do
  step post_one "$i" || echo "[task $SLURM_ARRAY_TASK_ID] row $i failed"
done
