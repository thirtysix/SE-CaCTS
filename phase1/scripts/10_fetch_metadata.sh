#!/usr/bin/env bash
# Phase 1 — 10: fetch the authoritative metadata (idempotent; skips files already present).
#   - Cellosaurus flat file (name/synonym -> CVCL accession), ~116 MB
#   - ChIP-Atlas experimentList.tab (SRX + antigen + cell type + per-experiment QC), ~328 MB
# Both land in the shared local caches so they're reusable across projects.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
SECACTS="$(cd "$HERE/../.." && pwd)"
. "$SECACTS/secacts_env.sh"
need_dataroot
D="$SECACTS_DATAROOT"
CELL="$D/cellosaurus/cellosaurus.txt"
EXP="$D/chip-atlas/experimentList.tab"

mkdir -p "$D/cellosaurus" "$D/chip-atlas"
[[ -s "$CELL" ]] && echo "[10] have $(du -h "$CELL"|cut -f1) cellosaurus.txt" || {
  echo "[10] fetching cellosaurus.txt ..."
  wget -q -O "$CELL" "https://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt"; }
[[ -s "$EXP" ]] && echo "[10] have $(du -h "$EXP"|cut -f1) experimentList.tab" || {
  echo "[10] fetching experimentList.tab ..."
  wget -q -O "$EXP" "https://chip-atlas.dbcls.jp/data/metadata/experimentList.tab"; }
echo "[10] done."
