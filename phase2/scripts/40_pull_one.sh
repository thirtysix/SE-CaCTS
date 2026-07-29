#!/usr/bin/env bash
# Phase 2 — 40: the per-sample unit of the pull. Download one bigWig, produce the full §8.5 retention set,
# delete the bigWig. Self-contained and resumable, so it is also runnable standalone for testing:
#
#   OUT=/scratch/.../out GRIDS="grid.20.bed grid.10.bed" bash 40_pull_one.sh SRX067407
#
# Write-safety (PULL_DESIGN.md §5.1): every output is SRX-keyed, so no two tasks write the same file.
# cnrose writes into a staging dir ON THE SAME FILESYSTEM as the destination, then each file is `mv`d into
# place — rename within one FS is atomic, so a reader never sees a half-written file. A cross-filesystem mv
# would be a copy and would NOT be atomic, which is why staging is on scratch even though the bigWig is on
# node-local disk. `.done` is written LAST and is the only resume key: a task killed mid-write leaves no
# marker and is simply redone.
set -uo pipefail

SRX="${1:?usage: 40_pull_one.sh <SRX>}"
OUT="${OUT:?set OUT to the output root}"
GRIDS="${GRIDS:?set GRIDS to a space-separated list of grid BEDs}"
PEAK_DIR="${PEAK_DIR:?set PEAK_DIR to the dir of <SRX>.20.bed}"
CNROSE="${CNROSE:?set CNROSE to the cnrose package root}"
PY="${PY:-python3}"
BINS="${BINS:-100}"
SIGFMT="${SIGFMT:-f32}"
LIMIT_RATE="${LIMIT_RATE:-}"          # e.g. 5m — courtesy throttle for chip-atlas.dbcls.jp (§6)
KEEP_BW="${KEEP_BW:-0}"               # 1 = don't delete (local testing only)

# Lustre metadata: shard outputs so no directory takes thousands of concurrent creates (§5.1)
SHARD="${SRX: -2}"
DEST="$OUT/$SHARD"
mkdir -p "$DEST"

[[ -f "$DEST/$SRX.done" ]] && { echo "[40] $SRX already done -> skip"; exit 0; }

PRIMARY="https://chip-atlas.dbcls.jp/data/hg38/eachData/bw/$SRX.bw"
LEGACY="http://dbarchive.biosciencedbc.jp/kyushu-u/hg38/eachData/bw/$SRX.bw"

# Transient bigWig on node-local disk when available. Deliberately NOT falling back to $TMPDIR or /tmp:
# on CSC those are often a RAM-backed tmpfs, so a 322-570 MB bigWig would count against the job's 2 GiB
# memory limit and OOM the task. $LOCAL_SCRATCH needs `--gres=nvme:N` on some CSC systems; when it is
# absent we use a scratch scratch-dir instead (slower, but Lustre-backed and off the memory budget).
BWDIR="${LOCAL_SCRATCH:-$OUT/.bw.${SLURM_JOB_ID:-local}.$$}"
mkdir -p "$BWDIR"
TMPBW="$BWDIR/$SRX.bw"
STAGE="$DEST/.staging.$SRX.${SLURM_JOB_ID:-local}"
rm -rf "$STAGE"; mkdir -p "$STAGE"
cleanup() { rm -rf "$STAGE"; [[ "$KEEP_BW" == "1" ]] || rm -f "$TMPBW";
            [[ -n "${LOCAL_SCRATCH:-}" ]] || rmdir "$BWDIR" 2>/dev/null || true; }
trap cleanup EXIT

# ---- download (resumable; primary then legacy endpoint)
URL_USED=""
for url in "$PRIMARY" "$LEGACY"; do
  echo "[40] $SRX fetching $url"
  if wget -q -c ${LIMIT_RATE:+--limit-rate="$LIMIT_RATE"} --tries=3 --timeout=60 -O "$TMPBW" "$url"; then
    [[ -s "$TMPBW" ]] && { URL_USED="$url"; break; }
  fi
  rm -f "$TMPBW"
done
if [[ -z "$URL_USED" ]]; then
  # per-task marker, NOT an append to a shared file (O_APPEND is not atomic across nodes on Lustre, §5.1)
  echo "download failed from both endpoints" > "$DEST/$SRX.failed"
  echo "[40] $SRX DOWNLOAD FAILED -> $DEST/$SRX.failed (continuing)"
  exit 0
fi

PEAKS="$PEAK_DIR/$SRX.20.bed"
if [[ ! -s "$PEAKS" ]]; then
  echo "no peak BED at $PEAKS" > "$DEST/$SRX.failed"
  echo "[40] $SRX NO PEAKS -> $DEST/$SRX.failed (continuing)"
  exit 0
fi

# ---- the one bigWig pass: SE calls + every grid + genome-wide archive + QC (cnrose.io.scan_bigwig)
GRID_ARGS=(); for g in $GRIDS; do GRID_ARGS+=(--grid "$g"); done
if ! PYTHONPATH="$CNROSE" "$PY" -m cnrose.cli call \
      --bw "$TMPBW" --peaks "$PEAKS" "${GRID_ARGS[@]}" \
      --signal-format "$SIGFMT" ${BINS:+--bins "$BINS"} --qc \
      --out "$STAGE/$SRX"; then
  echo "cnrose failed" > "$DEST/$SRX.failed"
  echo "[40] $SRX CNROSE FAILED -> $DEST/$SRX.failed (continuing)"
  exit 0
fi

# ---- provenance: what was fetched, from where, and with which code (PULL_DESIGN.md §8.5.4)
SHA=$(sha256sum "$TMPBW" | cut -d' ' -f1)
BYTES=$(stat -c%s "$TMPBW")
# cnrose git SHA: the staged tree is NOT a git repo (rsync copies files, not .git), so read the SHA the
# staging step recorded ($WORK/cnrose.gitsha); fall back to a local git query if we ARE in a repo.
GITSHA=$(cat "${WORK:-.}/cnrose.gitsha" 2>/dev/null || git -C "$CNROSE" rev-parse --short HEAD 2>/dev/null || echo unknown)
{
  printf '{\n'
  printf '  "srx": "%s",\n' "$SRX"
  printf '  "url": "%s",\n' "$URL_USED"
  printf '  "bytes": %s,\n' "$BYTES"
  printf '  "sha256": "%s",\n' "$SHA"
  printf '  "fetched_utc": "%s",\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf '  "peaks": "%s",\n' "$(basename "$PEAKS")"
  printf '  "grids": [%s],\n' "$(for g in $GRIDS; do printf '"%s",' "$(basename "$g")"; done | sed 's/,$//')"
  printf '  "bins": %s,\n' "${BINS:-null}"
  printf '  "cnrose_git": "%s",\n' "$GITSHA"
  printf '  "host": "%s",\n' "$(hostname)"
  printf '  "slurm_job": "%s"\n' "${SLURM_JOB_ID:-local}"
  printf '}\n'
} > "$STAGE/$SRX.provenance.json"

# ---- publish atomically (same filesystem => rename is atomic), then the .done marker LAST
for f in "$STAGE/$SRX".*; do mv -f "$f" "$DEST/$(basename "$f")"; done
touch "$DEST/$SRX.done"
rm -f "$DEST/$SRX.failed"        # clear any marker from an earlier failed attempt

echo "[40] $SRX OK -> $DEST ($(du -sh "$DEST/$SRX".* 2>/dev/null | awk '{s+=$1} END {print NR" files"}'))"
