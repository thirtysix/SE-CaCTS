#!/bin/bash
# v3.1 build on Roihu: reduce -> smoke (real data, 20 permutations) -> every scoring arm, chained afterok.
# Run from $SB (score_v3/) after syncing the repo copy:   PROJ=<project> bash v31_chain.sh
# Outputs: $SB/out_v31 (baseline), $SB/out_v31_pergroup (item 6), $SB/out_v31all (variants 4-5).
# FUSED=1: the fused build (FINDINGS §54-55) -> results_v31f{,all}, out_v31f{,_pergroup,all}, plus uncorrected twins of
# the per-line arms so every call the dashboard lists can be labelled CN-robust / amplicon-driven / CN-unmasked.
# VER=v32 (with FUSED=1): the v3.2 build (FINDINGS §57) from pull_set.v32*, v32_{inferred,randexcl}_keys, reduce_v32f.
set -euo pipefail
: "${PROJ:?}"
R=/scratch/$PROJ/se-cacts; W=$R/phase2; SB=$R/score_v3; REPO=$SB/dataroot/002.AI_projects/SE-CaCTS; D=$REPO/phase2/data
V=${VER:-v31}; T=$V${FUSED:+f}
cd "$SB"
rid=$(sbatch --parsable -A "$PROJ" --export=ALL,PROJ="$PROJ" "reduce_$T.slurm")
echo "[$T] reduce $rid"

base="ALL,PROJ=$PROJ,SB=$SB,W=$W,V31=1,RES=$W/results_$T,PS=$D/pull_set.$V.tsv,INFKEYS=$D/${V}_inferred_keys.txt,RANDKEYS=$D/${V}_randexcl_keys.txt,PERMIMPL=count,PERMWORKERS=4${FUSED:+,FUSED=1}"
all="ALL,PROJ=$PROJ,SB=$SB,W=$W,V31=1,RES=$W/results_${T}all,PS=$D/pull_set.$V.all.tsv,INFKEYS=$D/${V}_inferred_keys.txt,PERMIMPL=count,PERMWORKERS=4${FUSED:+,FUSED=1}"
G3=OncotreeLineage+OncotreePrimaryDisease+OncotreeSubtype   # "+": --export splits on commas

sid=$(sbatch --parsable -A "$PROJ" --cpus-per-task=4 --time=01:00:00 --dependency=afterok:"$rid" \
      --export="$base",ARM=main,NPERM=20,LEVELS=OncotreeLineage,OUTDIR=$SB/out_${T}_smoke score_arm.slurm)
sid2=$(sbatch --parsable -A "$PROJ" --cpus-per-task=4 --time=01:00:00 --dependency=afterok:"$rid" \
      --export="$base",ARM=lall,NPERM=5,OUTDIR=$SB/out_${T}_smoke score_arm.slurm)
echo "[$T] smoke $sid (group levels, mean path) $sid2 (cell-line studies, units path)"
sid="$sid:$sid2"

sub () {   # sub <export> <arm> <outdir> [extra export] [time]
  local j
  j=$(sbatch --parsable -A "$PROJ" --cpus-per-task=4 --time="${5:-03:00:00}" --dependency=afterok:"$sid" \
      --export="$1,ARM=$2,OUTDIR=$3${4:+,$4}" score_arm.slurm)
  echo "[$T] $2 -> $(basename "$3") $j"
}
for a in main nocn; do sub "$base" "$a" "$SB/out_$T" "LEVELS=$G3"; done
for a in shuffle noinf nodisc infnocn randexcl sw; do sub "$base" "$a" "$SB/out_$T"; done
twins=""; [ -n "${FUSED:-}" ] && twins="lsub_nocn ldis_nocn llin_nocn"
for a in csub csubshuf hsub hsubshuf cstudy cstudyshuf hline hlineshuf cline clineshuf \
         lall lallshuf lsub lsubshuf ldis ldisshuf llin llinshuf lall_nocn $twins; do
  sub "$base" "$a" "$SB/out_$T" "" 05:00:00
done
for a in csub csubshuf; do sub "$base" "$a" "$SB/out_${T}_pergroup" "FDRSCOPE=pergroup" 05:00:00; done
for a in main nocn; do sub "$all" "$a" "$SB/out_${T}all" "LEVELS=$G3"; done
sub "$all" shuffle "$SB/out_${T}all"
