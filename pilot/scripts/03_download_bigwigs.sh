#!/usr/bin/env bash
# Step B pilot — 03: download the selected per-experiment bigWigs from ChIP-Atlas.
#
# The ONLY network step. DRY-RUN by default (prints + HEAD-checks each URL); pass --go to download.
# On CSC/Roihu run this on the LOGIN NODE (the only node with internet), writing into project scratch.
#
# Usage:
#   bash 03_download_bigwigs.sh                 # dry-run, default selection ../data/selection.tsv
#   bash 03_download_bigwigs.sh --go            # download
#   OUT="$SECACTS_SCRATCH/bw" bash 03_download_bigwigs.sh --go            # Roihu scratch
set -euo pipefail

SEL="${SEL:-../data/selection.tsv}"
OUT="${OUT:-../data/bigwigs}"
GENOME="${GENOME:-hg38}"
GO=0
for arg in "$@"; do [[ "$arg" == "--go" ]] && GO=1; [[ "$arg" == *.tsv ]] && SEL="$arg"; done

# ChIP-Atlas per-experiment bigWig endpoints (primary first, legacy fallback).
BASE_PRIMARY="https://chip-atlas.dbcls.jp/data/${GENOME}/eachData/bw"
BASE_LEGACY="http://dbarchive.biosciencedbc.jp/kyushu-u/${GENOME}/eachData/bw"

[[ -f "$SEL" ]] || { echo "selection not found: $SEL (run 02_select_pilot.py first)"; exit 1; }
mkdir -p "$OUT"

# SRX = first column, skip header
mapfile -t SRXS < <(tail -n +2 "$SEL" | cut -f1 | grep -E '^SRX[0-9]+$')
echo "[03] ${#SRXS[@]} experiments; genome=$GENOME; out=$OUT; mode=$([[ $GO -eq 1 ]] && echo DOWNLOAD || echo DRY-RUN)"

fail=0
for SRX in "${SRXS[@]}"; do
  url="${BASE_PRIMARY}/${SRX}.bw"
  dst="${OUT}/${SRX}.bw"
  if [[ $GO -eq 1 ]]; then
    if [[ -s "$dst" ]]; then echo "  have  $SRX"; continue; fi
    if wget -c -q -O "$dst" "$url" && [[ -s "$dst" ]]; then
      echo "  ok    $SRX  ($(du -h "$dst" | cut -f1))"
    else
      echo "  RETRY $SRX via legacy endpoint"
      if wget -c -q -O "$dst" "${BASE_LEGACY}/${SRX}.bw" && [[ -s "$dst" ]]; then
        echo "  ok    $SRX  ($(du -h "$dst" | cut -f1)) [legacy]"
      else
        echo "  FAIL  $SRX"; rm -f "$dst"; fail=$((fail+1))
      fi
    fi
  else
    code=$(curl -s -o /dev/null -I -L -w '%{http_code}' "$url" || echo 000)
    printf '  %-12s HTTP %s  %s\n' "$SRX" "$code" "$url"
  fi
done
[[ $GO -eq 1 && $fail -gt 0 ]] && echo "[03] $fail downloads failed" && exit 2
echo "[03] done."
