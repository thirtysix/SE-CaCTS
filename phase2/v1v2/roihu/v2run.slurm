#!/bin/bash
#SBATCH --partition=small
#SBATCH --job-name=v1v2-v2run
#SBATCH --time=03:00:00
#SBATCH --mem=28G
#SBATCH --cpus-per-task=5
#SBATCH --output=logs/%x-%A_%a.out
# ChIP-Atlas pipeline v2, unmodified (their pipeline-v2.sh in their container), on one sample.
# Runs merged per experiment as ChIP-Atlas does; intermediates in a scratch TMPDIR (not RAM-backed /tmp).
set -euo pipefail
source ./common.sh
row=$(awk -F'\t' -v i=$((SLURM_ARRAY_TASK_ID + 2)) 'NR==i' "$SAMPLES")
srx=$(cut -f1 <<< "$row"); layout=$(cut -f6 <<< "$row")
d="$V2W/fastq/$srx"; out="$V2W/v2/$srx"
[ -s "$out/$srx.stats.tsv" ] && { echo "[v2run] $srx done"; exit 0; }
export TMPDIR="$V2W/tmp/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"; mkdir -p "$TMPDIR"
if [ "$layout" = PAIRED ]; then
  cat "$d"/*_1.fastq.gz > "$TMPDIR/r1.fq.gz"; cat "$d"/*_2.fastq.gz > "$TMPDIR/r2.fq.gz"
  args=(--fastq-fwd "$TMPDIR/r1.fq.gz" --fastq-rev "$TMPDIR/r2.fq.gz")
else
  cat "$d"/*.fastq.gz > "$TMPDIR/r1.fq.gz"; args=(--fastq-fwd "$TMPDIR/r1.fq.gz")
fi
apptainer exec --bind "$V2W" --env TMPDIR="$TMPDIR" "$SIF" bash "$REPO/scripts/pipeline-v2.sh" \
  --sample-id "$srx" "${args[@]}" --genome-fasta "$REF/hg38.fa" --chrom-sizes "$REF/chrom.sizes" \
  --genome-size hs --outdir "$out" --threads "$SLURM_CPUS_PER_TASK"
rm -rf "$TMPDIR"
cat "$out/$srx.stats.tsv"
