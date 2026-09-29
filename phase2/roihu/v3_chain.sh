#!/bin/bash
#SBATCH --partition=small
#SBATCH --job-name=v3-chain
#SBATCH --time=01:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/%x-%j.out
# The v3 build after the Stage B pool (FINDINGS 2026-09-29 §31). Submit from score_v3/ with
#   sbatch -A <project> --dependency=afterany:<stageB pool job> --export=ALL,PROJ=<project> v3_chain.sh
# 1. signal/identity check of the published Stage B experiments, as every other new experiment had (A3b)
# 2. final pull set: pull_set.v3.all.tsv minus Stage B rows that are unpublished or did not pass
# 3. reduce (catalog rebuilt, S3norm pinned to the v1/v2 reference) -> every scoring arm, chained afterok
set -euo pipefail
: "${PROJ:?}"
R=/scratch/$PROJ/se-cacts; W=$R/phase2; SB=$R/score_v3; REPO=$SB/dataroot/002.AI_projects/SE-CaCTS
PY=/projappl/$PROJ/secacts_venv/bin/python
ARMS="${ARMS:-main nocn shuffle noinf nodisc csub csubshuf cstudy cstudyshuf hsub hsubshuf hline hlineshuf cline clineshuf}"

cd "$W/v3check"
find "$W/out" -name '*.done' -printf '%f\n' | sed 's/\.done$//' > stageB_published.txt
awk -F'\t' 'NR==FNR{d[$1]=1; next} FNR==1 || ($1 in d)' stageB_published.txt stageB_candidates.tsv > stageB_candidates.done.tsv
echo "[chain] Stage B published: $(( $(wc -l < stageB_candidates.done.tsv) - 1 )) of $(( $(wc -l < stageB_candidates.tsv) - 1 ))"
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 "$PY" v3_candidate_check.py --work "$W" --meta-dir "$W/v3check" \
  --new-lines v3_new_lines.tsv --candidates stageB_candidates.done.tsv --out stageB_candidate_check.tsv

"$PY" - "$REPO/phase2/data" "$W/v3check/stageB_candidate_check.tsv" "$W/data/pull_srx.v3.txt" <<'EOF'
import sys, pandas as pd
d, chk, man = sys.argv[1:4]
ps = pd.read_csv(f"{d}/pull_set.v3.all.tsv", sep="\t")
c = pd.read_csv(chk, sep="\t")
ok = set(c.loc[c["status"] == "pass", "srx"])
sb = ps["source"] == "stageB"
print(f"[chain] Stage B: {int(sb.sum())} in pull set, {len(c)} checked, {len(ok)} pass; flags "
      f"{c.loc[c['status'] != 'pass', 'flags'].value_counts().head(6).to_dict()}")
ps = ps[~sb | ps["srx"].isin(ok)]
ps.to_csv(f"{d}/pull_set.v3.tsv", sep="\t", index=False)
ps["srx"].sort_values().to_csv(man, index=False, header=False)
print(f"[chain] final pull set: {len(ps)} experiments on {ps['key'].nunique()} keys")
EOF

cd "$W"
rid=$(CATALOGS="atlas atlas.s3" sbatch --parsable -A "$PROJ" \
      --export=ALL,PROJ="$PROJ",MANIFEST="$W/data/pull_srx.v3.txt",RESULTS="$W/results_v3",S3_REF=SRX16495452 reduce.slurm)
echo "[chain] reduce $rid"
cd "$SB"
for arm in $ARMS; do
  j=$(sbatch --parsable -A "$PROJ" --dependency=afterok:"$rid" \
      --export=ALL,PROJ="$PROJ",SB="$SB",W="$W",ARM="$arm",RES="$W/results_v3",PS="$REPO/phase2/data/pull_set.v3.tsv" \
      score_arm.slurm)
  echo "[chain] arm $arm $j"
done
