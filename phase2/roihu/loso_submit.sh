#!/bin/bash
# Leave-one-study-out lineage rescore (FINDINGS §42): one score_arm.slurm job per exclusion set.
# Sets come from phase2/analysis/loso_sets.py, synced to $SB/loso/sets/. Run on Roihu from $SB:
#   PROJ=<project> bash loso_submit.sh [set-name ...]      (default: every *.srx not yet done)
set -euo pipefail
: "${PROJ:?}"
R=/scratch/$PROJ/se-cacts; W=$R/phase2; SB=$R/score_v3; REPO=$SB/dataroot/002.AI_projects/SE-CaCTS
cd "$SB"
names=("$@")
[ ${#names[@]} -eq 0 ] && mapfile -t names < <(ls loso/sets/*.srx | xargs -n1 basename | sed 's/\.srx$//')
for n in "${names[@]}"; do
  out=$SB/out_loso/$n
  [ -s "$out/atlas.s3.perm.OncotreeLineage.specific.tsv.gz" ] && { echo "[loso] $n done, skipped"; continue; }
  j=$(sbatch --parsable -A "$PROJ" --job-name=loso --cpus-per-task=4 --time=01:00:00 \
      --export=ALL,PROJ="$PROJ",SB="$SB",W="$W",ARM=loso,RES="$W/results_v3",PS="$REPO/phase2/data/pull_set.v3.tsv",EXCL="$SB/loso/sets/$n.srx",OUTDIR="$out",PERMIMPL=count,PERMWORKERS=4 \
      score_arm.slurm)
  echo "[loso] $n $j"
done
