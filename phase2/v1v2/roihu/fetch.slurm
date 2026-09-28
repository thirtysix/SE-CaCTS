#!/bin/bash
#SBATCH --partition=small
#SBATCH --job-name=v1v2-fetch
#SBATCH --time=03:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x-%A_%a.out
# Per sample: ENA FASTQ (the raw reads v2 starts from) and ChIP-Atlas's v1 bigWig (deleted after our atlas pull).
set -euo pipefail
source ./common.sh
row=$(awk -F'\t' -v i=$((SLURM_ARRAY_TASK_ID + 2)) 'NR==i' "$SAMPLES")
srx=$(cut -f1 <<< "$row"); layout=$(cut -f6 <<< "$row"); ftp=$(cut -f10 <<< "$row")
d="$V2W/fastq/$srx"; mkdir -p "$d"
for u in ${ftp//;/ }; do
  f=$(basename "$u")
  case "$layout:$f" in PAIRED:*_1.fastq.gz|PAIRED:*_2.fastq.gz|SINGLE:*) ;; *) continue ;; esac
  if [ ! -s "$d/$f" ]; then curl -sfL --retry 3 -o "$d/$f.part" "https://$u"; mv "$d/$f.part" "$d/$f"; fi
done
if [ ! -s "$V2W/v1/$srx.bw" ]; then
  curl -sfL --retry 3 -o "$V2W/v1/$srx.bw.part" "https://chip-atlas.dbcls.jp/data/hg38/eachData/bw/$srx.bw"
  mv "$V2W/v1/$srx.bw.part" "$V2W/v1/$srx.bw"
fi
ls -la "$d" "$V2W/v1/$srx.bw"
