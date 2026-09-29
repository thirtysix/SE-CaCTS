#!/bin/bash
#SBATCH --partition=small
#SBATCH --job-name=sB-fetch
#SBATCH --time=12:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x-%A_%a.out
# Per experiment: ENA FASTQs, runs concatenated per mate (sraTailor.sh concatenates runs per SRX). One core:
# the step is download-bound, so it never holds the aligner's 8 cores idle.
set -euo pipefail
source ./common.sh
# Roihu caps a user at 200 queued / 100 running jobs on `small`, so each array task takes SAMPLES_PER_TASK
# consecutive manifest rows, one after another; a failed row never stops the rest of the task.
process() {
  r=$(row "$1"); srx=$(cut -f1 <<<"$r"); layout=$(cut -f2 <<<"$r"); ftp=$(cut -f3 <<<"$r")
  runs=$(cut -f9 <<<"$r")
  d="$S/work/$srx"; mkdir -p "$d/runs"
  [ -f "$d/.fetched" ] && { echo "[fetch] $srx done"; exit 0; }
  # SRA first: fasterq-dump reproduces the read names ChIP-Atlas v1 aligned, which Bowtie2 uses to break ties
  # between equally good placements, so the output matches v1 exactly (validation 2026-09-29); ENA FASTQs carry
  # other names and differ slightly in repeats. --check-rs no: do not pre-download the references of
  # reference-compressed runs (84 chromosomes, one timed out); fasterq-dump resolves them as needed.
  fetch_sra() {
    local SRA="$AX $S/sif/sra-tools.sif" run
    rm -rf "$d/runs"; mkdir -p "$d/runs"
    for run in ${runs//;/ }; do
      $SRA prefetch --max-size 500G --check-rs no -O "$d/runs" "$run" || return 1
      $SRA fasterq-dump --split-files --threads 1 -t "$d/runs" -O "$d/runs" "$d/runs/$run/$run.sra" || return 1
      rm -rf "$d/runs/$run"
    done
    if [ "$layout" = PAIRED ]; then
      cat $(ls "$d"/runs/*_1.fastq | sort) > "$d/r1.fq" && cat $(ls "$d"/runs/*_2.fastq | sort) > "$d/r2.fq" || return 1
    else
      cat $(ls "$d"/runs/*.fastq | sort) > "$d/r1.fq" || return 1
    fi
    rm -rf "$d/runs"
    [ -s "$d/r1.fq" ] && { [ "$layout" != PAIRED ] || [ -s "$d/r2.fq" ]; }
  }
  if [ -n "$runs" ]; then
    if fetch_sra; then touch "$d/.fetched"; echo "[fetch] $srx via SRA"; ls -la "$d"; exit 0; fi
    echo "[fetch] $srx: SRA route failed; falling back to ENA"; rm -f "$d"/r1.fq "$d"/r2.fq; rm -rf "$d/runs"; mkdir -p "$d/runs"
    [ -n "$ftp" ] || { echo "[fetch] $srx: no ENA FASTQ either"; exit 1; }
  fi
  for u in ${ftp//;/ }; do
    f=$(basename "$u")
    case "$layout:$f" in PAIRED:*_1.fastq.gz|PAIRED:*_2.fastq.gz|SINGLE:*) ;; *) continue ;; esac
    [ -s "$d/runs/$f" ] || { curl -sfL --retry 5 --retry-delay 15 -o "$d/runs/$f.part" "https://$u" && mv "$d/runs/$f.part" "$d/runs/$f"; }
  done
  if [ "$layout" = PAIRED ]; then
    cat $(ls "$d"/runs/*_1.fastq.gz | sort) > "$d/r1.fq.gz"; cat $(ls "$d"/runs/*_2.fastq.gz | sort) > "$d/r2.fq.gz"
  else
    cat $(ls "$d"/runs/*.fastq.gz | sort) > "$d/r1.fq.gz"
  fi
  rm -rf "$d/runs"
  # never mark a partial fetch complete: an empty mate file would be aligned as single-end (validation, 2026-09-29)
  [ -s "$d/r1.fq.gz" ] && { [ "$layout" != PAIRED ] || [ -s "$d/r2.fq.gz" ]; } || { echo "[fetch] $srx incomplete"; exit 1; }
  touch "$d/.fetched"; ls -la "$d"
}
K="${SAMPLES_PER_TASK:-1}"; N=$(( $(wc -l < "$MANIFEST") - 1 ))
for (( i = SLURM_ARRAY_TASK_ID * K; i < (SLURM_ARRAY_TASK_ID + 1) * K && i < N; i++ )); do
  ( process "$i" ) || echo "[task $SLURM_ARRAY_TASK_ID] row $i failed"
done
