#!/usr/bin/env bash
# Phase 2 — 21: build the GRID from PER-SRX ChIP-Atlas peak BEDs, fetched directly on Roihu,
# limited to our analysis-set experiments, at a chosen (permissive) threshold.
#
# Avoids the 107 GB AllCell download: only our 2,917 SRX, only peak coords -> ~0.1-6 GB by threshold.
# Portable: needs only curl + sort + awk (no bedtools/pigz). Run ON ROIHU (login or compute node;
# both have internet).  Endpoint pattern verified 2026-07-18.
#
#   NN=10 SRXLIST=pull_srx.txt OUT=grid_persrx bash 21_grid_from_persrx.sh
#   LIMIT=40 NN=10 ... bash 21_grid_from_persrx.sh      # validation: first 40 SRX only
#
# NN = ChIP-Atlas filename threshold: 05=Q<1E-05(loosest) 10=1E-10 20=1E-20 50=1E-50(stringent)
set -uo pipefail

NN="${NN:-10}"
SRXLIST="${SRXLIST:-pull_srx.txt}"
OUT="${OUT:-grid_persrx}"
JOBS="${JOBS:-16}"
LIMIT="${LIMIT:-0}"                       # 0 = all; >0 = only first N (validation)
BEDDIR="$OUT/bed$NN"
PRIMARY="https://chip-atlas.dbcls.jp/data/hg38/eachData/bed$NN"
LEGACY="http://dbarchive.biosciencedbc.jp/kyushu-u/hg38/eachData/bed$NN"
mkdir -p "$BEDDIR"

[[ -f "$SRXLIST" ]] || { echo "missing SRX list: $SRXLIST"; exit 1; }
mapfile -t SRXS < <(grep -E '^[A-Z]RX[0-9]+$' "$SRXLIST")
[[ "$LIMIT" -gt 0 ]] && SRXS=("${SRXS[@]:0:$LIMIT}")
NSRX=${#SRXS[@]}
echo "[21] threshold bed$NN | $NSRX SRX | jobs=$JOBS | out=$BEDDIR"

# write the (post-LIMIT) work list to a file and drive xargs from it — robust vs array expansion.
# Shell vars are expanded at PARENT time via quote-breaking; each SRX arrives as $1 (no exported fn).
printf '%s\n' "${SRXS[@]}" > "$OUT/srx.work.$NN"
echo "[21] work list: $(wc -l < "$OUT/srx.work.$NN") SRX -> fetching (primary->legacy, 404 = no peaks)"
xargs -P "$JOBS" -I{} bash -c '
  s="$1"; [ -n "$s" ] || exit 0
  d="'"$BEDDIR"'/$s.'"$NN"'.bed"; [ -s "$d" ] && exit 0
  curl -fsS --noproxy "*" --max-time 60 -o "$d" "'"$PRIMARY"'/$s.'"$NN"'.bed" 2>/dev/null && [ -s "$d" ] && exit 0
  curl -fsS --noproxy "*" --max-time 60 -o "$d" "'"$LEGACY"'/$s.'"$NN"'.bed"  2>/dev/null && [ -s "$d" ] && exit 0
  rm -f "$d"; exit 0
' _ {} < "$OUT/srx.work.$NN"

NFILES=$(find "$BEDDIR" -name "*.$NN.bed" ! -name ".$NN.bed" | wc -l)
NMISS=$(( NSRX - NFILES ))
DLBYTES=$(du -sb "$BEDDIR" 2>/dev/null | cut -f1); DLBYTES=${DLBYTES:-0}
printf '[21] fetched %s beds (%.2f GB); %s SRX had no bed%s (no peaks this strict)\n' \
  "$NFILES" "$(awk -v b=$DLBYTES 'BEGIN{print b/1e9}')" "$NMISS" "$NN"

# merge -> grid: canonical chroms, cols 1-3, sort, streaming awk merge (overlap/bookend)
cat "$BEDDIR"/*.$NN.bed \
 | awk -F'\t' '$1 ~ /^chr([1-9][0-9]?|X|Y)$/ {print $1"\t"$2"\t"$3}' \
 | sort -k1,1 -k2,2n --parallel="$JOBS" -S 2G \
 | awk 'BEGIN{OFS="\t"} NR==1{c=$1;s=$2;e=$3;next}
        {if($1==c&&$2<=e){if($3>e)e=$3}else{print c,s,e;c=$1;s=$2;e=$3}} END{if(NR>0)print c,s,e}' \
 > "$OUT/grid.$NN.bed"

REG=$(wc -l < "$OUT/grid.$NN.bed")
BP=$(awk '{s+=$3-$2} END{print s+0}' "$OUT/grid.$NN.bed")
printf '[21] grid.%s.bed : %s regions | %d bp (%.2f%% genome) | SRX with peaks: %d/%d\n' \
  "$NN" "$REG" "$BP" "$(awk -v b=$BP 'BEGIN{print 100*b/3.1e9}')" "$NFILES" "$NSRX"