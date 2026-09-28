#!/bin/bash
#SBATCH --partition=small
#SBATCH --job-name=v1v2-index
#SBATCH --time=06:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x-%j.out
# bwa-mem2 index of UCSC hg38.fa inside the v2 container (what prepare-genomes.sh runs). Single-threaded, memory-bound.
set -euo pipefail
source ./common.sh
[ -s "$REF/hg38.fa.bwt.2bit.64" ] || apptainer exec --bind "$V2W" "$SIF" bwa-mem2 index "$REF/hg38.fa"
ls -la "$REF"
