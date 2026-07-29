#!/usr/bin/env bash
# Phase 2 — 43: run the pull in STAGES, so errors surface on 2 samples rather than 2,917.
#
#   bash 43_run_stage.sh smoke          # dry run: print what would be submitted
#   bash 43_run_stage.sh smoke  --go    #     2 samples,  %1  — mechanics, download rate, seff
#   bash 43_run_stage.sh small  --go    #    32 samples,  %8  — concurrency, throttling, Lustre
#   bash 43_run_stage.sh medium --go    #   256 samples,  %8  — stragglers, failure rate at scale
#   bash 43_run_stage.sh large  --go    #   512 samples, %16  — full-run concurrency, endpoint courtesy
#   bash 43_run_stage.sh full   --go    # 2,917 samples, %16  — the real thing
#   bash 43_run_stage.sh report         # summarise whatever has completed so far
#
# Stages are PREFIXES of the same manifest and share one output tree, so each stage extends the last:
# `.done` markers make later stages skip everything earlier stages finished. Nothing is recomputed, and
# array indices stay stable because a prefix of the manifest indexes identically to the whole.
set -uo pipefail
# Machine-local config: walk up to the repo root and source .env via secacts_env.sh.
# Depth-independent so this works from a staged copy on the cluster too; if the file is
# absent (partial stage), a real SECACTS_CSC_PROJECT in the environment still satisfies the
# `:?` checks below.
_d="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [ "$_d" != "/" ] && [ ! -f "$_d/secacts_env.sh" ]; do _d="$(dirname "$_d")"; done
[ -f "$_d/secacts_env.sh" ] && . "$_d/secacts_env.sh"

# NB: keep this a plain test, not ${1:?...with {braces}...} — braces inside :? terminate the expansion
# early and silently append the trailing text to $1 (bit us 2026-07-21; bash -n does not catch it).
STAGE="${1:-}"
[[ -n "$STAGE" ]] || { echo "usage: 43_run_stage.sh smoke|small|medium|full|report [--go]" >&2; exit 1; }
GO=0; [[ "${2:-}" == "--go" ]] && GO=1
PROJ="${PROJ:-${SECACTS_CSC_PROJECT:?set SECACTS_CSC_PROJECT in .env, or pass PROJ=}}"
REMOTE="${REMOTE:-roihu}"
WORK="/scratch/$PROJ/se-cacts/phase2"

case "$STAGE" in
  smoke)  N=2;    SPT=2;  CONC=1  ;;
  small)  N=32;   SPT=8;  CONC=4  ;;
  medium) N=256;  SPT=16; CONC=8  ;;
  large)  N=512;  SPT=16; CONC=16 ;;   # last check before full: exercises the FULL-run concurrency (%16)
  full)   N=0;    SPT=16; CONC=16 ;;   # N=0 -> whole manifest
  report) N=-1 ;;
  *) echo "unknown stage '$STAGE'" >&2; exit 1 ;;
esac

remote() { ssh -o BatchMode=yes "$REMOTE" "$@"; }

# ---------------------------------------------------------------- report
if [[ "$STAGE" == "report" ]]; then
  remote "cd $WORK && bash -s" <<'EOS'
set -u
tot=$(grep -cv '^[[:space:]]*$' data/pull_srx.txt)
done=$(find out -name '*.done' 2>/dev/null | wc -l)
fail=$(find out -name '*.failed' 2>/dev/null | wc -l)
echo "progress: $done/$tot done, $fail failed"
[[ $fail -gt 0 ]] && { echo "--- failures ---"; find out -name '*.failed' -printf '%f: ' -exec head -c 120 {} \; -exec echo \; | head -20; }
echo "--- disk ---"; du -sh out 2>/dev/null
du -sb out 2>/dev/null | awk '{printf "  %.1f GB on scratch -> %.1f BU/day (6 BU/TiB-h)\n", $1/1e9, 6*($1/1024^4)*24}'
echo "--- per-sample QC spread (dynamic range; the S3norm gate) ---"
find out -name '*.qc.json' | head -400 | xargs -r grep -h '"dynamic_range_p99_over_median"' 2>/dev/null \
  | sed -E 's/.*:[[:space:]]*([0-9.]+).*/\1/' | sort -n | awk '{a[NR]=$1} END {if(NR) printf "  n=%d  min=%.1f  p25=%.1f  median=%.1f  p75=%.1f  max=%.1f\n", NR, a[1], a[int(NR*.25)+1], a[int(NR/2)+1], a[int(NR*.75)+1], a[NR]}'
echo "--- recent slurm accounting (BU proxy: CPUTime, MaxRSS) ---"
sacct -X -n --format=JobID%18,State%12,Elapsed%10,MaxRSS%10,AllocCPUS%4 -S $(date -d '2 days ago' +%F) 2>/dev/null | tail -15
EOS
  exit 0
fi

# ---------------------------------------------------------------- submit
TOT=$(remote "grep -cv '^[[:space:]]*$' $WORK/data/pull_srx.txt" 2>/dev/null || echo 0)
[[ "$TOT" -gt 0 ]] || { echo "cannot read the manifest on $REMOTE — run 42_preflight.sh" >&2; exit 1; }
(( N == 0 )) && N=$TOT
(( N > TOT )) && N=$TOT
TASKS=$(( (N + SPT - 1) / SPT ))
DONE=$(remote "find $WORK/out -name '*.done' 2>/dev/null | wc -l" || echo 0)

cat <<EOF
[43] stage=$STAGE
  samples this stage : $N of $TOT   (already done: $DONE -> those tasks will skip)
  samples/task       : $SPT   -> --array=0-$((TASKS-1))%$CONC
  concurrency        : %$CONC          transient disk ~$(( CONC * 322 / 1024 )) GiB
  est. wall-clock    : ~$(( N * 60 / CONC / 60 )) min      est. compute ~$(( N * 60 / 3600 )) core-h = ~$(awk "BEGIN{printf \"%.1f\", $N*60/3600*0.75}") BU
  est. output        : ~$(( N * 20 / 1000 )) GB  -> ~$(awk "BEGIN{printf \"%.2f\", 6*($N*20e6/1024^4)*24}") BU/day on scratch until fetched
EOF

SUB="sbatch --parsable --array=0-$((TASKS-1))%$CONC --export=ALL,SAMPLES_PER_TASK=$SPT array.slurm"
if (( ! GO )); then
  echo "[43] DRY RUN — would submit (add --go):"
  echo "     cd $WORK && $SUB"
  exit 0
fi

echo "[43] submitting"
aid=$(remote "cd $WORK && $SUB")
echo "[43] array job $aid submitted"
echo "[43] watch:   ssh $REMOTE 'squeue -j $aid'"
echo "[43] progress: bash 43_run_stage.sh report"
if [[ "$STAGE" != "full" ]]; then
  echo "[43] when it finishes, CHECK BEFORE THE NEXT STAGE:"
  echo "       ssh $REMOTE 'seff \$(sacct -X -n -j $aid --format=JobID%20 | head -1)'   # real BU shape"
  echo "       bash 43_run_stage.sh report                                              # failures, QC spread"
else
  echo "[43] then: sbatch --dependency=afterok:$aid reduce.slurm  (on $REMOTE, in $WORK)"
  echo "[43] then: bash 50_fetch_and_clean.sh <local-dest> --clean   (scratch bills 7.6 BU/day)"
fi
