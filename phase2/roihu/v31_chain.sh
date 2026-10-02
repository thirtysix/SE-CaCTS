#!/bin/bash
# v3.1 build on Roihu: reduce -> smoke (real data, 20 permutations) -> every scoring arm, chained afterok.
# Run from $SB (score_v3/) after syncing the repo copy:   PROJ=<project> bash v31_chain.sh
# Outputs: $SB/out_v31 (baseline), $SB/out_v31_pergroup (item 6), $SB/out_v31all (variants 4-5).
set -euo pipefail
: "${PROJ:?}"
R=/scratch/$PROJ/se-cacts; W=$R/phase2; SB=$R/score_v3; REPO=$SB/dataroot/002.AI_projects/SE-CaCTS; D=$REPO/phase2/data
cd "$SB"
rid=$(sbatch --parsable -A "$PROJ" --export=ALL,PROJ="$PROJ" reduce_v31.slurm)
echo "[v31] reduce $rid"

base="ALL,PROJ=$PROJ,SB=$SB,W=$W,V31=1,RES=$W/results_v31,PS=$D/pull_set.v31.tsv,INFKEYS=$D/v31_inferred_keys.txt,RANDKEYS=$D/v31_randexcl_keys.txt,PERMIMPL=count,PERMWORKERS=4"
all="ALL,PROJ=$PROJ,SB=$SB,W=$W,V31=1,RES=$W/results_v31all,PS=$D/pull_set.v31.all.tsv,INFKEYS=$D/v31_inferred_keys.txt,PERMIMPL=count,PERMWORKERS=4"
G3=OncotreeLineage,OncotreePrimaryDisease,OncotreeSubtype

sid=$(sbatch --parsable -A "$PROJ" --cpus-per-task=4 --time=01:00:00 --dependency=afterok:"$rid" \
      --export="$base",ARM=main,NPERM=20,LEVELS=OncotreeLineage,OUTDIR=$SB/out_v31_smoke score_arm.slurm)
sid2=$(sbatch --parsable -A "$PROJ" --cpus-per-task=4 --time=01:00:00 --dependency=afterok:"$rid" \
      --export="$base",ARM=lall,NPERM=5,OUTDIR=$SB/out_v31_smoke score_arm.slurm)
echo "[v31] smoke $sid (group levels, mean path) $sid2 (cell-line studies, units path)"
sid="$sid:$sid2"

sub () {   # sub <export> <arm> <outdir> [extra export] [time]
  local j
  j=$(sbatch --parsable -A "$PROJ" --cpus-per-task=4 --time="${5:-03:00:00}" --dependency=afterok:"$sid" \
      --export="$1,ARM=$2,OUTDIR=$3${4:+,$4}" score_arm.slurm)
  echo "[v31] $2 -> $(basename "$3") $j"
}
for a in main nocn; do sub "$base" "$a" "$SB/out_v31" "LEVELS=$G3"; done
for a in shuffle noinf nodisc infnocn randexcl sw; do sub "$base" "$a" "$SB/out_v31"; done
for a in csub csubshuf hsub hsubshuf cstudy cstudyshuf hline hlineshuf cline clineshuf \
         lall lallshuf lsub lsubshuf ldis ldisshuf llin llinshuf lall_nocn; do
  sub "$base" "$a" "$SB/out_v31" "" 05:00:00
done
for a in csub csubshuf; do sub "$base" "$a" "$SB/out_v31_pergroup" "FDRSCOPE=pergroup" 05:00:00; done
for a in main nocn; do sub "$all" "$a" "$SB/out_v31all" "LEVELS=$G3"; done
sub "$all" shuffle "$SB/out_v31all"
