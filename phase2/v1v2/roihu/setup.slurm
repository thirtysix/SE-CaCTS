#!/bin/bash
#SBATCH --partition=small
#SBATCH --job-name=v1v2-setup
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/%x-%j.out
# Pull ChIP-Atlas's v2 container and code (pinned), and the reference exactly as their prepare-genomes.sh does
# (UCSC hg38.fa.gz, samtools faidx -> chrom.sizes). The bwa-mem2 index is its own job (index.slurm).
set -euo pipefail
source ./common.sh
cd "$V2W"
[ -s "$SIF" ] || apptainer pull "$SIF" "$V2IMG"
[ -d "$REPO/.git" ] || git clone -q https://github.com/inutano/chip-atlas-pipeline-v2 "$REPO"
git -C "$REPO" rev-parse HEAD > "$V2W/v2_repo_commit.txt"
if [ ! -s "$REF/hg38.fa" ]; then
  curl -sL -o "$REF/hg38.fa.gz" https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
  gunzip -f "$REF/hg38.fa.gz"
fi
[ -s "$REF/hg38.fa.fai" ] || apptainer exec --bind "$V2W" "$SIF" samtools faidx "$REF/hg38.fa"
cut -f1,2 "$REF/hg38.fa.fai" > "$REF/chrom.sizes"
apptainer exec "$SIF" bash -c 'bwa-mem2 version 2>&1 | head -1; macs3 --version; fastp --version 2>&1 | head -1' > "$V2W/v2_tool_versions.txt"
echo "[setup] v2 commit $(cat "$V2W/v2_repo_commit.txt"); contigs $(wc -l < "$REF/chrom.sizes")"
