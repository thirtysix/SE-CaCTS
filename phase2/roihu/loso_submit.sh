#!/bin/bash
# Leave-one-study-out lineage rescore (FINDINGS §42, §46): one score_arm.slurm job per exclusion set.
# Sets come from phase2/analysis/loso_sets.py, synced to $SB/<setdir>/sets/. Run on Roihu from $SB:
#   PROJ=<project> bash loso_submit.sh [set-name ...]      (default: every *.srx not yet done)
# V31F=1: the v3.1 release build (fused presence, the SE-of-its-group rule in the null): sets in loso_v31f/sets,
# outputs in out_loso_v31f/ (v3: loso/sets -> out_loso/)
set -euo pipefail
: "${PROJ:?}"
R=/scratch/$PROJ/se-cacts; W=$R/phase2; SB=$R/score_v3; REPO=$SB/dataroot/002.AI_projects/SE-CaCTS; D=$REPO/phase2/data
if [ -n "${V31F:-}" ]; then
  SETS=loso_v31f/sets; OUTS=out_loso_v31f
  EXP="V31=1,FUSED=1,RES=$W/results_v31f,PS=$D/pull_set.v31.tsv"
else
  SETS=loso/sets; OUTS=out_loso
  EXP="RES=$W/results_v3,PS=$D/pull_set.v3.tsv"
fi
cd "$SB"
names=("$@")
[ ${#names[@]} -eq 0 ] && mapfile -t names < <(ls "$SETS"/*.srx | xargs -n1 basename | sed 's/\.srx$//')
for n in "${names[@]}"; do
  out=$SB/$OUTS/$n
  [ -s "$out/atlas.s3.perm.OncotreeLineage.specific.tsv.gz" ] && { echo "[loso] $n done, skipped"; continue; }
  j=$(sbatch --parsable -A "$PROJ" --job-name=loso --cpus-per-task=4 --time=01:00:00 \
      --export=ALL,PROJ="$PROJ",SB="$SB",W="$W",ARM=loso,$EXP,EXCL="$SB/$SETS/$n.srx",OUTDIR="$out",PERMIMPL=count,PERMWORKERS=4 \
      score_arm.slurm)
  echo "[loso] $n $j"
done
