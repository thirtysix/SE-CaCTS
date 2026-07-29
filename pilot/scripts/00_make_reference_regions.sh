#!/usr/bin/env bash
# Step B pilot — 00 (OPTIONAL, local): build an H3K27ac peak reference for 04's `bed` mode.
#
# Only needed if you want SE-relevant regions instead of genome-wide bins. Runs LOCALLY (needs the 3.1 GB
# AllCell bed); emits a small bed you can stage to HPC. For the default `bins` mode this is not required.
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
SECACTS="$(cd "$HERE/../.." && pwd)"
. "$SECACTS/secacts_env.sh"
need_dataroot
BED="${BED:-$SECACTS_DATAROOT/chip-atlas/00.data/Histone/His.ALL.50.H3K27ac.AllCell.bed.gz}"
OUT="${OUT:-../data/reference_regions.bed}"
N="${N:-50000}"          # subsample to N regions (0 = keep all merged)
SEED="${SEED:-1}"

mkdir -p "$(dirname "$OUT")"
echo "[00] merging H3K27ac peaks from $(basename "$BED") ..."
TMP="$(mktemp)"
zcat "$BED" | grep -vE '^track' | cut -f1-3 \
  | grep -E '^chr([0-9]+|X|Y)\b' \
  | sort -k1,1 -k2,2n \
  | bedtools merge > "$TMP"
total=$(wc -l < "$TMP")
if [[ "$N" -gt 0 && "$total" -gt "$N" ]]; then
  awk -v seed="$SEED" 'BEGIN{srand(seed)} {print rand()"\t"$0}' "$TMP" \
    | sort -k1,1g | head -n "$N" | cut -f2- | sort -k1,1 -k2,2n > "$OUT"
else
  cp "$TMP" "$OUT"
fi
rm -f "$TMP"
echo "[00] merged=$total  ->  kept=$(wc -l < "$OUT")  ->  $OUT"
