#!/usr/bin/env bash
# Phase 2 — 50: pull the retention set OFF Roihu and free the scratch space. Run FROM THE LAPTOP.
#
#   bash 50_fetch_and_clean.sh <local-dest>              # fetch + verify, leave Roihu untouched
#   bash 50_fetch_and_clean.sh <local-dest> --clean      # fetch + verify + DELETE from Roihu scratch
#   ARCHIVE_ONLY=1 bash 50_fetch_and_clean.sh <dest> --clean   # move only the bulky bin100 archive
#
# WHY THIS IS NOT OPTIONAL. Roihu bills scratch at 6 BU/TiB-hour **from the first byte** — unlike Puhti,
# there is no free tier (docs.csc.fi/computing/hpc-billing). The retention set is ~58 GB:
#
#     compute for the whole 2,917-sample pull   ~36 BU   (one-off)
#     leaving that output on scratch            ~7.6 BU/day
#
# so the storage overtakes the entire cost of the pull after **~4.8 days**, and costs ~229 BU/month
# thereafter. The genome-wide archive is ~44 GB of the ~58 GB, i.e. ~75% of that bleed. Fetch it down
# promptly; scratch is also un-backed-up and subject to cleaning, so it was never a home for it anyway.
set -euo pipefail
# Machine-local config: walk up to the repo root and source .env via secacts_env.sh.
# Depth-independent so this works from a staged copy on the cluster too; if the file is
# absent (partial stage), a real SECACTS_CSC_PROJECT in the environment still satisfies the
# `:?` checks below.
_d="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [ "$_d" != "/" ] && [ ! -f "$_d/secacts_env.sh" ]; do _d="$(dirname "$_d")"; done
[ -f "$_d/secacts_env.sh" ] && . "$_d/secacts_env.sh"

DEST="${1:?usage: 50_fetch_and_clean.sh <local-dest> [--clean]}"
CLEAN=0; [[ "${2:-}" == "--clean" ]] && CLEAN=1
PROJ="${PROJ:-${SECACTS_CSC_PROJECT:?set SECACTS_CSC_PROJECT in .env, or pass PROJ=}}"
REMOTE="${REMOTE:-roihu}"
WORK="/scratch/$PROJ/se-cacts/phase2"
ARCHIVE_ONLY="${ARCHIVE_ONLY:-0}"

mkdir -p "$DEST"
echo "[50] remote usage before:"
ssh "$REMOTE" "du -sh $WORK/out $WORK/results 2>/dev/null; echo '--- BU/day at 6 BU/TiB-h:'; \
  du -sb $WORK/out 2>/dev/null | awk '{printf \"    %.1f GB -> %.1f BU/day\\n\", \$1/1e9, 6*(\$1/1024^4)*24}'"

if [[ "$ARCHIVE_ONLY" == "1" ]]; then
  FILTER=(--include='*/' --include='*.bin*.f16.gz' --include='*.bin*.f16.gz.json' --exclude='*')
  echo "[50] fetching ARCHIVE ONLY -> $DEST"
else
  FILTER=()
  echo "[50] fetching the FULL retention set -> $DEST"
fi

rsync -azP --info=progress2 "${FILTER[@]}" "$REMOTE:$WORK/out/" "$DEST/out/"
rsync -azP "$REMOTE:$WORK/results/" "$DEST/results/" 2>/dev/null || true

# ---- verify before deleting anything: counts must match, and every .done must have its files
echo "[50] verifying"
r_done=$(ssh "$REMOTE" "find $WORK/out -name '*.done' | wc -l")
l_done=$(find "$DEST/out" -name '*.done' 2>/dev/null | wc -l)
if [[ "$ARCHIVE_ONLY" == "1" ]]; then
  r_n=$(ssh "$REMOTE" "find $WORK/out -name '*.bin*.f16.gz' | wc -l")
  l_n=$(find "$DEST/out" -name '*.bin*.f16.gz' | wc -l)
  echo "[50]   archives: remote=$r_n local=$l_n"
  [[ "$r_n" == "$l_n" && "$r_n" -gt 0 ]] || { echo "[50] MISMATCH — not cleaning" >&2; exit 1; }
else
  echo "[50]   .done markers: remote=$r_done local=$l_done"
  [[ "$r_done" == "$l_done" && "$r_done" -gt 0 ]] || { echo "[50] MISMATCH — not cleaning" >&2; exit 1; }
  # spot-check byte-for-byte on a sample of files rather than trusting rsync alone
  ssh "$REMOTE" "cd $WORK/out && find . -name '*.grid.20.f32' | head -20 | xargs md5sum" > /tmp/.r.md5
  (cd "$DEST/out" && md5sum -c /tmp/.r.md5 --quiet) \
    && echo "[50]   spot-check: 20 grid columns match byte-for-byte" \
    || { echo "[50] CHECKSUM MISMATCH — not cleaning" >&2; exit 1; }
fi

if (( CLEAN )); then
  if [[ "$ARCHIVE_ONLY" == "1" ]]; then
    echo "[50] deleting archives from Roihu scratch (keeping the exact tier for downstream work)"
    ssh "$REMOTE" "find $WORK/out -name '*.bin*.f16.gz' -delete -o -name '*.bin*.f16.gz.json' -delete"
  else
    echo "[50] deleting the whole per-sample output tree from Roihu scratch"
    ssh "$REMOTE" "rm -rf $WORK/out"
  fi
  echo "[50] remote usage after:"
  ssh "$REMOTE" "du -sh $WORK 2>/dev/null"
else
  echo "[50] NOT cleaning (pass --clean once you are happy). Every day it stays costs ~7.6 BU."
fi

echo "[50] local: $(du -sh "$DEST" | cut -f1) at $DEST"
echo "[50] REMINDER: this is the only copy now — it is not backed up. Consider Allas or an external disk."
