#!/usr/bin/env bash
# NOTE: --account is NOT hardcoded. secacts_env.sh exports SBATCH_ACCOUNT from
#       SECACTS_CSC_PROJECT (.env), which sbatch honours. Override per-run with
#       `sbatch --account=...`.
#SBATCH --partition=small
#SBATCH --job-name=secacts-gzip
#SBATCH --time=00:30:00
#SBATCH --mem=1G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#
# Phase 2 — 51: gzip the results/ atlas matrices IN A JOB (never on the login node — CSC etiquette; heavy
# CPU/IO belongs on compute nodes). The atlas TSVs are plain text: presence matrices compress ~141x (0/1),
# signal matrices ~2.3x (float text). ~2.1 GB -> ~0.7 GB. Idempotent: skips files already gzipped.
#
#   sbatch phase2/scripts/51_compress_results.sh            # from $WORK on Roihu
#
# Uses `pigz` (parallel gzip) across --cpus-per-task if available, else plain gzip. 4 cores keeps this a
# ~1-minute job; the BU is negligible either way.
set -uo pipefail
# Machine-local config: walk up to the repo root and source .env via secacts_env.sh.
# Depth-independent so this works from a staged copy on the cluster too; if the file is
# absent (partial stage), a real SECACTS_CSC_PROJECT in the environment still satisfies the
# `:?` checks below.
_d="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [ "$_d" != "/" ] && [ ! -f "$_d/secacts_env.sh" ]; do _d="$(dirname "$_d")"; done
[ -f "$_d/secacts_env.sh" ] && . "$_d/secacts_env.sh"
PROJ="${PROJ:-${SECACTS_CSC_PROJECT:?set SECACTS_CSC_PROJECT in .env, or pass PROJ=}}"
WORK="${WORK:-/scratch/$PROJ/se-cacts/phase2}"
RESULTS="${RESULTS:-$WORK/results}"
cd "$RESULTS"

if command -v pigz >/dev/null 2>&1; then ZIP="pigz -6 -p ${SLURM_CPUS_PER_TASK:-4}"; else ZIP="gzip -6"; fi
echo "[51] compressing $RESULTS with: $ZIP"
before=$(du -sb . | cut -f1)
# gzip every .tsv / .bed that isn't already gzipped; run a few in parallel when using plain gzip
shopt -s nullglob
for f in *.tsv *.bed; do
  [[ -f "$f" && ! -f "$f.gz" ]] && { echo "  $f"; $ZIP -f "$f" & }
  while (( $(jobs -r | wc -l) >= ${SLURM_CPUS_PER_TASK:-4} )); do wait -n; done
done
wait
after=$(du -sb . | cut -f1)
echo "[51] done: $(awk "BEGIN{printf \"%.1f GB -> %.1f GB (%.1fx)\", $before/1e9, $after/1e9, $before/$after}")"
ls -lh *.gz 2>/dev/null | awk '{print "  "$5, $9}'
