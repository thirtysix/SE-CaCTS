#!/usr/bin/env bash
# Phase 2 — 30: rehearse the full per-sample -> reduce pipeline LOCALLY on the 13 pilot bigWigs.
#
# Proves the whole Phase-2 mechanic (cnrose call --grid  ->  aggregate.py) end-to-end before any Roihu
# BU is spent. Uses the atac_hdac env (pyBigWig+numpy+scipy). Idempotent: skips work already done.
#
# CN=1 also produces the CN-CORRECTED catalog (dual catalogs, DESIGN.md §3.2): a fast second cnrose pass
# WITHOUT --grid (reuses the uncorrected fine matrix — CN is separable), keyed by each SRX's DepMap
# ModelID from pull_set.tsv, then a second reduce with --se-suffix .cn.
#
#   BED20_DIR=... bash 30_rehearse_local.sh              # agnostic only
#   CN=1 BED20_DIR=... bash 30_rehearse_local.sh         # + CN-corrected dual catalog
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
SECACTS="$(cd "$HERE/../.." && pwd)"
PY="${PY:-${SECACTS_PY:-$HOME/miniconda3/envs/atac_hdac/bin/python}}"
CNROSE="$SECACTS/cnrose"
GRID="$SECACTS/phase2/data/grid.20.bed"
. "$SECACTS/secacts_env.sh"
need_dataroot
DATAROOT="$SECACTS_DATAROOT"

BW_DIR="${BW_DIR:-$SECACTS/pilot/data/bigwigs}"
BED20_DIR="${BED20_DIR:?set BED20_DIR to the dir of <SRX>.20.bed}"
OUT="${OUT:-$SECACTS/phase2/rehearse/out}"
mkdir -p "$OUT"

# CN backend paths (only used when CN=1)
CN="${CN:-0}"
CNCSV="${CNCSV:-$DATAROOT/DepMap/2026q1/OmicsCNGeneWGS.csv}"
GTF="${GTF:-$DATAROOT/0.human_genome/Homo_sapiens.GRCh38.106.chr.gtf.gz}"
GENE_CACHE="${GENE_CACHE:-$OUT/../gene_coords.GRCh38.106.tsv}"
PULLSET="$SECACTS/phase2/data/pull_set.tsv"

srxs=$(ls "$BW_DIR"/*.bw | xargs -n1 basename | sed 's/\.bw$//')
n=$(echo "$srxs" | wc -w)
echo "[30] $n samples | grid=$(basename "$GRID") | CN=$CN | out=$OUT"

i=0
for srx in $srxs; do
  i=$((i+1))
  peaks="$BED20_DIR/$srx.20.bed"
  [[ -s "$peaks" ]] || { echo "  [$i/$n] MISS peaks $srx -> skip"; continue; }

  # expensive pass: agnostic + grid fine-matrix column (once)
  if [[ ! -s "$OUT/$srx.signal.tsv" || ! -s "$OUT/$srx.se.bed" ]]; then
    echo "  [$i/$n] call+grid  $srx"
    PYTHONPATH="$CNROSE" "$PY" -m cnrose.cli call \
      --bw "$BW_DIR/$srx.bw" --peaks "$peaks" --grid "$GRID" --out "$OUT/$srx"
  fi

  # fast CN pass: corrected catalog only, reuse the fine matrix (no --grid)
  if [[ "$CN" == "1" && ! -s "$OUT/$srx.cn.se.bed" ]]; then
    mid=$(awk -F'\t' -v s="$srx" '$1==s{print $4; exit}' "$PULLSET")
    if [[ -n "$mid" ]]; then
      echo "  [$i/$n] call+CN   $srx  (ModelID=$mid)"
      PYTHONPATH="$CNROSE" "$PY" -m cnrose.cli call \
        --bw "$BW_DIR/$srx.bw" --peaks "$peaks" --out "$OUT/$srx" \
        --cn depmap --cn-key "$mid" --cn-gene-csv "$CNCSV" --cn-gtf "$GTF" --cn-gene-cache "$GENE_CACHE"
    else
      echo "  [$i/$n] no ModelID for $srx in pull_set -> agnostic only"
    fi
  fi
done

echo "[30] reduce (agnostic) -> aggregate.py"
"$PY" "$SECACTS/phase2/aggregate.py" \
  --se-dir "$OUT" --signal-dir "$OUT" --grid "$GRID" --out "$SECACTS/phase2/rehearse/pilot"

if [[ "$CN" == "1" ]]; then
  echo "[30] reduce (CN-corrected) -> aggregate.py --se-suffix .cn"
  "$PY" "$SECACTS/phase2/aggregate.py" \
    --se-dir "$OUT" --signal-dir "$OUT" --grid "$GRID" --se-suffix .cn \
    --out "$SECACTS/phase2/rehearse/pilot.cn"
fi

# S3NORM=1 also emits an S3norm-normalized reduce. S3norm is fitted on the FINE grid and re-summed here
# (the transform is nonlinear), so it must happen inside aggregate.py — score with `score_pilot.py --norm none`.
#
# The MIN_PEAKS QC gate is NOT optional here: s3norm fits an expanding exponent (B>1) on low-enrichment
# samples, amplifying their noise into apparent specificity. Ungated it made the pilot WORSE than quantile
# (Breast ESR1 #5 -> #21); gated it matches quantile on biology (ESR1 #2) with the scale properly equalized.
# 2000 sits in a natural gap in the pilot's peak counts (2046 vs 361).
if [[ "${S3NORM:-0}" == "1" ]]; then
  MIN_PEAKS="${MIN_PEAKS:-2000}"
  echo "[30] reduce (s3norm, QC >=${MIN_PEAKS} peaks) -> aggregate.py --norm s3norm"
  "$PY" "$SECACTS/phase2/aggregate.py" \
    --se-dir "$OUT" --signal-dir "$OUT" --grid "$GRID" --norm s3norm \
    --peak-dir "$BED20_DIR" --min-peaks "$MIN_PEAKS" \
    --out "$SECACTS/phase2/rehearse/pilot.s3"
  if [[ "$CN" == "1" ]]; then
    "$PY" "$SECACTS/phase2/aggregate.py" \
      --se-dir "$OUT" --signal-dir "$OUT" --grid "$GRID" --se-suffix .cn --norm s3norm \
      --peak-dir "$BED20_DIR" --min-peaks "$MIN_PEAKS" \
      --out "$SECACTS/phase2/rehearse/pilot.s3.cn"
  fi
fi
echo "[30] done -> $SECACTS/phase2/rehearse/pilot{,.cn,.s3,.s3.cn}.*"
