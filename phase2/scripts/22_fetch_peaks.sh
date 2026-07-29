#!/usr/bin/env bash
# Phase 2 — 22: fetch ONLY the per-SRX ChIP-Atlas peak BEDs the pull reads, into $OUT (default data/bed20).
#
#   NN=20 LIMIT=40 bash 22_fetch_peaks.sh          # first 40 SRX (smoke+small); LIMIT=0 = all 2,917
#
# This is the download half of 21_grid_from_persrx.sh WITHOUT the grid-building tail. That distinction is
# not cosmetic: script 21 writes `$OUT/grid.$NN.bed`, so pointing its OUT at data/ silently OVERWRITES the
# authoritative full grid.20.bed (504,855 regions) with a subset grid built from only the SRX fetched this
# run. The pull's fixed-grid design depends on that file being immutable — caught 2026-07-21 when preflight
# reported grid.20.bed had shrunk to 82,460 lines. This script only ever writes peak BEDs; it never touches
# a grid. Run ON ROIHU (login or compute node; both have internet). Idempotent: existing beds are skipped.
set -uo pipefail

NN="${NN:-20}"
SRXLIST="${SRXLIST:-data/pull_srx.txt}"
OUT="${OUT:-data/bed$NN}"
JOBS="${JOBS:-8}"                          # keep modest — courtesy to chip-atlas.dbcls.jp (PULL_DESIGN §6)
LIMIT="${LIMIT:-0}"                        # 0 = all; >0 = only first N (staged validation)
PRIMARY="https://chip-atlas.dbcls.jp/data/hg38/eachData/bed$NN"
LEGACY="http://dbarchive.biosciencedbc.jp/kyushu-u/hg38/eachData/bed$NN"

[[ -f "$SRXLIST" ]] || { echo "[22] missing SRX list: $SRXLIST" >&2; exit 1; }
mkdir -p "$OUT"
mapfile -t SRXS < <(grep -E '^[A-Z]RX[0-9]+$' "$SRXLIST")
[[ "$LIMIT" -gt 0 ]] && SRXS=("${SRXS[@]:0:$LIMIT}")
NSRX=${#SRXS[@]}
echo "[22] threshold bed$NN | $NSRX SRX | jobs=$JOBS -> $OUT"

printf '%s\n' "${SRXS[@]}" | xargs -P "$JOBS" -I{} bash -c '
  s="$1"; [ -n "$s" ] || exit 0
  d="'"$OUT"'/$s.'"$NN"'.bed"; [ -s "$d" ] && exit 0
  curl -fsS --noproxy "*" --max-time 60 -o "$d" "'"$PRIMARY"'/$s.'"$NN"'.bed" 2>/dev/null && [ -s "$d" ] && exit 0
  curl -fsS --noproxy "*" --max-time 60 -o "$d" "'"$LEGACY"'/$s.'"$NN"'.bed"  2>/dev/null && [ -s "$d" ] && exit 0
  rm -f "$d"; exit 0
' _ {}

NFILES=$(find "$OUT" -maxdepth 1 -name "*.$NN.bed" | wc -l)
echo "[22] have $NFILES/$NSRX bed$NN files ($(du -sh "$OUT" 2>/dev/null | cut -f1)); $(( NSRX - NFILES )) missing (404 = no peaks at this threshold)"
