#!/bin/bash
#SBATCH --partition=small
#SBATCH --job-name=v1v2-cmp
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x-%A_%a.out
# Per sample: call SEs on the v2 bigWig with the same cnrose + grid the pull used on v1 (and once more on the
# v1 peaks, to separate the coverage effect from the peak-calling effect), then compare_v1_v2.py sample.
# Submit with --dependency=aftercorr:<v2run array> so each task follows its own v2 run.
set -euo pipefail
source ./common.sh
row=$(awk -F'\t' -v i=$((SLURM_ARRAY_TASK_ID + 2)) 'NR==i' "$SAMPLES")
srx=$(cut -f1 <<< "$row"); layout=$(cut -f6 <<< "$row")
v1d="$W/out/${srx: -2}"; v2d="$V2W/v2/$srx"; c="$V2W/cmp/$srx"; mkdir -p "$c"
[ -s "$v2d/$srx.stats.tsv" ] || { echo "[cmp] $srx: no v2 output"; exit 1; }
ALT="$V2W/ref/altLocations.txt.gz"
if [ ! -s "$ALT" ]; then
  curl -sfL --retry 3 -o "$ALT.$SLURM_ARRAY_TASK_ID" \
    https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/altLocations.txt.gz
  mv "$ALT.$SLURM_ARRAY_TASK_ID" "$ALT"
fi
PY="$VENV/bin/python"; export PYTHONPATH="$W/cnrose"
[ -s "$c/$srx.v2.se.bed" ] || "$PY" -m cnrose.cli call --bw "$v2d/$srx.bw" \
  --peaks "$v2d/$srx.20_peaks.narrowPeak" --grid "$W/data/grid.20.bed" --signal-format f32 --bins 100 --qc \
  --out "$c/$srx.v2"
[ -s "$c/$srx.v2bw_v1pk.se.bed" ] || "$PY" -m cnrose.cli call --bw "$v2d/$srx.bw" \
  --peaks "$W/data/bed20/$srx.20.bed" --out "$c/$srx.v2bw_v1pk"
"$PY" "$V2W/compare_v1_v2.py" sample --srx "$srx" --layout "$layout" --v1-bw "$V2W/v1/$srx.bw" \
  --v1-dir "$v1d" --v1-peaks "$W/data/bed20/$srx.20.bed" --v2-dir "$v2d" --cmp-dir "$c" \
  --grid "$W/data/grid.20.bed" --catalog "$V2W/atlas.s3.union_catalog.bed.gz" --alt "$ALT" \
  --replicates "$V2W/replicates.tsv" --pull-out "$W/out" --out "$V2W/cmp/$srx.metrics.json"
