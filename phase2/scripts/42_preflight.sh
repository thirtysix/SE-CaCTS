#!/usr/bin/env bash
# Phase 2 — 42: preflight the Roihu environment BEFORE submitting anything. Read-only, costs no BU.
#
#   bash 42_preflight.sh
#
# Every check here is something that would otherwise fail a task at a time, hours into a 2,917-task array.
# Exit 0 = safe to submit the smoke stage.
set -uo pipefail
# Machine-local config: walk up to the repo root and source .env via secacts_env.sh.
# Depth-independent so this works from a staged copy on the cluster too; if the file is
# absent (partial stage), a real SECACTS_CSC_PROJECT in the environment still satisfies the
# `:?` checks below.
_d="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
while [ "$_d" != "/" ] && [ ! -f "$_d/secacts_env.sh" ]; do _d="$(dirname "$_d")"; done
[ -f "$_d/secacts_env.sh" ] && . "$_d/secacts_env.sh"

PROJ="${PROJ:-${SECACTS_CSC_PROJECT:?set SECACTS_CSC_PROJECT in .env, or pass PROJ=}}"
REMOTE="${REMOTE:-roihu}"
WORK="/scratch/$PROJ/se-cacts/phase2"
fail=0
ok()   { echo "  [ OK ] $*"; }
bad()  { echo "  [FAIL] $*"; fail=$((fail+1)); }
warn() { echo "  [warn] $*"; }

echo "== 0. connectivity =="
if ! timeout 25 ssh -o BatchMode=yes -o ConnectTimeout=15 "$REMOTE" true 2>/dev/null; then
  bad "cannot ssh $REMOTE — re-sign the certificate (valid 24 h): MyCSC > SSH public keys > Sign,"
  echo "         or: python3 ~/.ssh/csc_cert.py -u \$USER ~/.ssh/id_ed25519.pub  ->  ~/.ssh/id_ed25519-cert.pub"
  exit 1
fi
ok "ssh $REMOTE"

remote() { ssh -o BatchMode=yes "$REMOTE" "$@" 2>/dev/null; }

echo "== 1. project, partitions, quotas =="
[[ -n "$(remote "ls -d /scratch/$PROJ")" ]] && ok "/scratch/$PROJ" || bad "/scratch/$PROJ missing"
[[ -n "$(remote "ls -d /projappl/$PROJ")" ]] && ok "/projappl/$PROJ" || bad "/projappl/$PROJ missing"
parts=$(remote "sinfo -h -o %P" | tr '\n' ' ')
echo "         partitions: $parts"
grep -q small <<<"$parts" && ok "partition 'small' present" || bad "no 'small' partition (edit array.slurm)"
echo "         quota:"; remote "csc-workspaces 2>/dev/null | head -20" | sed 's/^/         /'
echo "         billing units left:"; remote "csc-projects -p $PROJ 2>/dev/null | head -10" | sed 's/^/         /'

echo "== 2. toolchain (cnrose needs pyBigWig + numpy + scipy only) =="
VENV="/projappl/$PROJ/secacts_venv"
tool=$(remote "$VENV/bin/python -c 'import pyBigWig,numpy,scipy,pandas;print(\"venv\")' 2>/dev/null")
if [[ "$tool" == "venv" ]]; then
  ok "venv provides pyBigWig+numpy+scipy+pandas ($VENV)"
else
  bad "no python with pyBigWig at $VENV — run: ssh $REMOTE 'bash -lc \"\$(cat)\"' < phase2/roihu/00_build_env.sh"
fi

echo "== 3. staged inputs =="
for f in aggregate.py s3norm.py specificity.py array.slurm reduce.slurm scripts/40_pull_one.sh; do
  [[ -n "$(remote "ls $WORK/$f 2>/dev/null")" ]] && ok "$f" || bad "$f not staged (run 41_stage_and_submit.sh)"
done
for f in grid.20.bed grid.10.bed pull_srx.txt; do
  n=$(remote "wc -l < $WORK/data/$f 2>/dev/null" || echo 0)
  [[ "${n:-0}" -gt 0 ]] && ok "data/$f ($n lines)" || bad "data/$f missing or empty"
done
nbed=$(remote "ls $WORK/data/bed20 2>/dev/null | wc -l" || echo 0)
nsrx=$(remote "wc -l < $WORK/data/pull_srx.txt 2>/dev/null" || echo 0)
if [[ "${nbed:-0}" -ge "${nsrx:-1}" ]]; then ok "bed20 peaks: $nbed for $nsrx samples"
elif [[ "${nbed:-0}" -gt 0 ]]; then warn "bed20 peaks: only $nbed of $nsrx — the rest will write .failed"
else bad "no bed20 peaks — fetch them on Roihu (21_grid_from_persrx.sh)"; fi

echo "== 4. cnrose imports on a compute-visible python =="
imp=$(remote "PYTHONPATH=$WORK/cnrose $VENV/bin/python -c 'from cnrose.io import scan_bigwig; print(\"import-ok\")' 2>/dev/null")
[[ "$imp" == "import-ok" ]] && ok "cnrose imports on the venv python" || bad "cnrose fails to import on the venv"

echo "== 5. outbound internet from a COMPUTE node (the §2 assumption) =="
# Confirmed 2026-07-21: a test-partition probe returned HTTP 200 from a compute node (~0.001 BU). We record
# it so re-running preflight does not keep submitting probes. Set PROBE=1 to force a fresh one.
if [[ "${PROBE:-0}" != "1" && -n "$(remote "cat $WORK/.compute_internet_ok 2>/dev/null")" ]]; then
  ok "compute-node internet confirmed earlier ($(remote "cat $WORK/.compute_internet_ok")); PROBE=1 to re-check"
else
  echo "         submitting a 2-minute test-partition probe..."
  probe=$(remote "cd $WORK && sbatch --account=$PROJ --partition=test --time=00:02:00 --mem=1G \
    --cpus-per-task=1 --wrap='code=\$(curl -sI -o /dev/null -w %{http_code} \
    https://chip-atlas.dbcls.jp/data/hg38/eachData/bw/SRX067407.bw); \
    [ \"\$code\" = 200 ] && date -u +%Y-%m-%dT%H:%MZ > $WORK/.compute_internet_ok; echo HTTP \$code' \
    --parsable" 2>/dev/null)
  [[ -n "$probe" ]] && ok "probe job $probe submitted — check: ssh $REMOTE 'cat $WORK/slurm-$probe.out'" \
                    || warn "could not submit probe (no 'test' partition?); the smoke stage will confirm"
fi

echo
if (( fail )); then echo "PREFLIGHT: $fail blocking issue(s) — fix before submitting."; exit 1; fi
echo "PREFLIGHT: clear. Next: bash 43_run_stage.sh smoke --go"
