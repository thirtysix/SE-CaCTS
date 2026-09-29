#!/bin/bash
#SBATCH --partition=small
#SBATCH --job-name=sB-dumptest
#SBATCH --time=02:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=6
#SBATCH --output=logs/%x-%j.out
# Is fasterq-dump output byte-identical at --threads 6 vs 1, and how much faster? Also times prefetch alone.
set -euo pipefail
source ./common.sh
run=${RUN:-SRR22377967}; T="$S/work/_dumptest"; rm -rf "$T"; mkdir -p "$T"; cd "$T"
SRA="$AX $S/sif/sra-tools.sif"
t0=$(date +%s); $SRA prefetch --max-size 500G --check-rs no -O "$T" "$run"; t1=$(date +%s)
echo "prefetch_s=$((t1-t0)) sra_bytes=$(du -sb $T/$run | cut -f1)"
for th in 6 1; do
  mkdir -p "t$th"; a=$(date +%s)
  $SRA fasterq-dump --split-files --threads $th -t "$T/t$th" -O "$T/t$th" "$T/$run/$run.sra"
  echo "threads=$th dump_s=$(( $(date +%s) - a )) bytes=$(du -cb $T/t$th/*.fastq | tail -1 | cut -f1)"
  (cd "t$th" && md5sum *.fastq)
done
for f in t1/*.fastq; do cmp -s "$f" "t6/$(basename $f)" && echo "IDENTICAL $(basename $f)" || echo "DIFFERENT $(basename $f)"; done
cd "$S/work"; rm -rf "$T"
