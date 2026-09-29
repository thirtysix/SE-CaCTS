#!/bin/bash
#SBATCH --partition=small
#SBATCH --job-name=sB-align
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=logs/%x-%A_%a.out
# Align step only (Bowtie2 + sort), the multi-threaded part.
# 8 cores / 16 GB sits at the 2 GiB/core billing break-even. Everything single-threaded is in post.slurm.
set -euo pipefail
source ./common.sh
# One step per array task (align_one in common.sh). Roihu caps a user at 200 queued / 100 running jobs on `small`,
# so each task takes SAMPLES_PER_TASK consecutive manifest rows; a failed row never stops the rest of the task.
K="${SAMPLES_PER_TASK:-1}"; N=$(( $(wc -l < "$MANIFEST") - 1 ))
for (( i = SLURM_ARRAY_TASK_ID * K; i < (SLURM_ARRAY_TASK_ID + 1) * K && i < N; i++ )); do
  step align_one "$i" || echo "[task $SLURM_ARRAY_TASK_ID] row $i failed"
done
