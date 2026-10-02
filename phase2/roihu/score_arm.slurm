#!/bin/bash
#SBATCH --partition=small
#SBATCH --job-name=secacts-score
#SBATCH --time=03:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x-%j-%a.out
#SBATCH --error=logs/%x-%j-%a.err
# One scoring arm of the v2 atlas (see phase2/scripts/70_score_v2.sh for the arm definitions).
source /appl/profile/zz-csc-env.sh 2>/dev/null || true
set -euo pipefail
SB="${SB:?}"; W="${W:?}"; ARM="${ARM:?}"
export SECACTS_DATAROOT="$SB/dataroot" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
PY=/projappl/${PROJ:?}/secacts_venv/bin/python
REPO="$SB/dataroot/002.AI_projects/SE-CaCTS"
OUT="${OUTDIR:-$SB/out}"; mkdir -p "$OUT"
# ARM=<arm>_nocn runs <arm> with CN correction off, writing <prefix>.nocn.* (the uncorrected atlas side by side)
BASE="${ARM%_nocn}"; NOCN=$([ "$BASE" != "$ARM" ] && echo 1 || true)
DISC="ACH-000232,ACH-000338,ACH-000856,ACH-000957,ACH-000960,ACH-000997,ACH-001106,ACH-000768"
# RES / PS pick the atlas (default v2). v3 adds input-inferred CN lines, so those two flags come along.
RES="${RES:-$W/results_v2}"; PS="${PS:-$REPO/phase2/data/pull_set.v2.tsv}"
V2=(--signal "$RES/atlas.s3.se_signal.tsv" --catalog "$RES/atlas.s3.union_catalog.bed"
    --pull-set "$PS" --norm none
    --fdr-method permutation --n-perm "${NPERM:-1000}" --keep-frac "${KEEPFRAC:-0.05}" --dump-specific 0.10
    --perm-impl "${PERMIMPL:-tail}" --perm-workers "${PERMWORKERS:-1}" --fdr-scope "${FDRSCOPE:-global}"
    --inferred-bins "$W/../cncal/bins" --blacklist "$W/../cncal/hg38-blacklist.v2.bed.gz")
case "$BASE" in
  main)    args=("${V2[@]}" --out "$OUT/atlas.s3.perm" --cn-diagnostic --save-matrices) ;;
  nocn)    args=("${V2[@]}" --out "$OUT/atlas.s3.perm.nocn" --no-cn) ;;
  shuffle) args=("${V2[@]}" --out "$OUT/atlas.s3.perm.shuffle" --shuffle-labels 1
                 --levels OncotreeLineage,OncotreePrimaryDisease,OncotreeSubtype) ;;
  noinf)   args=("${V2[@]}" --out "$OUT/atlas.s3.perm.noinf" --exclude-keys "$(cat "$REPO/phase2/data/v3_inferred_keys.txt")"
                 --levels OncotreeLineage,OncotreePrimaryDisease) ;;
  noinfall) args=("${V2[@]}" --out "$OUT/atlas.s3.perm.noinfall" --exclude-keys "$(cat "$REPO/phase2/data/v3_inferred_keys.txt")") ;;
  nodisc)  args=("${V2[@]}" --out "$OUT/atlas.s3.perm.nodisc" --exclude-keys "$DISC"
                 --levels OncotreeLineage,OncotreePrimaryDisease) ;;
  cmpfirst) args=("${V2[@]}" --pull-set "$REPO/phase2/data/pull_set.v2.cmpfirst.tsv" --out "$OUT/atlas.s3.perm.cmpfirst"
                 --cn-diagnostic --levels OncotreeLineage,OncotreePrimaryDisease) ;;
  layout)   args=("${V2[@]}" --signal "$W/results_v2/atlas.s3.layoutcorr.se_signal.tsv"
                 --out "$OUT/atlas.s3.perm.layoutcorr" --levels OncotreeLineage,OncotreePrimaryDisease) ;;
  # leave-one-study-out (FINDINGS §42): EXCL = SRX list from phase2/analysis/loso_sets.py, lineage level only
  loso)     args=("${V2[@]}" --out "$OUT/atlas.s3.perm" --exclude-srx "${EXCL:?}" --levels OncotreeLineage) ;;
  csub)     args=("${V2[@]}" --out "$OUT/atlas.s3.cons" --agg q25 --levels OncotreeSubtype) ;;
  csubshuf) args=("${V2[@]}" --out "$OUT/atlas.s3.cons.shuffle" --agg q25 --levels OncotreeSubtype --shuffle-labels 1) ;;
  cline)    args=("${V2[@]}" --out "$OUT/atlas.s3.consx" --agg q25 --levels line --line-members experiments --keep-frac 0.01) ;;
  clineshuf) args=("${V2[@]}" --out "$OUT/atlas.s3.consx.shuffle" --agg q25 --levels line --line-members experiments --keep-frac 0.01 --shuffle-labels 1) ;;
  cstudy)   args=("${V2[@]}" --out "$OUT/atlas.s3.conss" --agg q25 --levels line --line-members studies --keep-frac 0.01) ;;
  cstudyshuf) args=("${V2[@]}" --out "$OUT/atlas.s3.conss.shuffle" --agg q25 --levels line --line-members studies --keep-frac 0.01 --shuffle-labels 1) ;;
  hsub)     args=("${V2[@]}" --out "$OUT/atlas.s3.hyb" --agg q25 --levels OncotreeSubtype --null-strata lineage) ;;
  hsubshuf) args=("${V2[@]}" --out "$OUT/atlas.s3.hyb.shuffle" --agg q25 --levels OncotreeSubtype --null-strata lineage --shuffle-labels 1) ;;
  # per-line lists for the dashboard (2026-09-30): vs all lines, and vs relatives at subtype / disease / lineage.
  # Untestable lines (one study; no label at the level; < 4 lines in the stratum) are left out of BH.
  lall)     args=("${V2[@]}" --out "$OUT/atlas.s3.lines.all" --agg q25 --levels line --line-members studies --keep-frac 0.01 --strata-min 4) ;;
  lallshuf) args=("${V2[@]}" --out "$OUT/atlas.s3.lines.all.shuffle" --agg q25 --levels line --line-members studies --keep-frac 0.01 --strata-min 4 --shuffle-labels 1) ;;
  lsub)     args=("${V2[@]}" --out "$OUT/atlas.s3.lines.sub" --agg q25 --levels line --line-members studies --keep-frac 0.01 --strata-min 4 --null-strata subtype) ;;
  lsubshuf) args=("${V2[@]}" --out "$OUT/atlas.s3.lines.sub.shuffle" --agg q25 --levels line --line-members studies --keep-frac 0.01 --strata-min 4 --null-strata subtype --shuffle-labels 1) ;;
  ldis)     args=("${V2[@]}" --out "$OUT/atlas.s3.lines.dis" --agg q25 --levels line --line-members studies --keep-frac 0.01 --strata-min 4 --null-strata disease) ;;
  ldisshuf) args=("${V2[@]}" --out "$OUT/atlas.s3.lines.dis.shuffle" --agg q25 --levels line --line-members studies --keep-frac 0.01 --strata-min 4 --null-strata disease --shuffle-labels 1) ;;
  llin)     args=("${V2[@]}" --out "$OUT/atlas.s3.lines.lin" --agg q25 --levels line --line-members studies --keep-frac 0.01 --strata-min 4 --null-strata lineage) ;;
  llinshuf) args=("${V2[@]}" --out "$OUT/atlas.s3.lines.lin.shuffle" --agg q25 --levels line --line-members studies --keep-frac 0.01 --strata-min 4 --null-strata lineage --shuffle-labels 1) ;;
  hline)    args=("${V2[@]}" --out "$OUT/atlas.s3.hybl" --agg q25 --levels line --line-members studies --null-strata adaptive --keep-frac 0.01) ;;
  hlineshuf) args=("${V2[@]}" --out "$OUT/atlas.s3.hybl.shuffle" --agg q25 --levels line --line-members studies --null-strata adaptive --keep-frac 0.01 --shuffle-labels 1) ;;
  hsmoke)   args=("${V2[@]}" --out "$OUT/hsmoke" --agg q25 --levels OncotreeSubtype,line --line-members studies --null-strata adaptive --n-perm 3 --keep-frac 0.01) ;;
  smoke)    args=("${V2[@]}" --out "$OUT/smoke" --agg q25 --levels OncotreeSubtype,line --line-members experiments --n-perm 3 --keep-frac 0.01) ;;
  anshuf)  args=(--signal "$RES/atlas.s3.se_signal.tsv" --catalog "$RES/atlas.s3.union_catalog.bed"
                 --pull-set "$PS" --norm none --inferred-bins "$W/../cncal/bins" --blacklist "$W/../cncal/hg38-blacklist.v2.bed.gz" --levels OncotreeLineage
                 --dump-specific 0.10 --shuffle-labels 1 --out "$OUT/atlas.s3.analytic.shuffle") ;;
  anreal)  args=(--signal "$RES/atlas.s3.se_signal.tsv" --catalog "$RES/atlas.s3.union_catalog.bed"
                 --pull-set "$PS" --norm none --inferred-bins "$W/../cncal/bins" --blacklist "$W/../cncal/hg38-blacklist.v2.bed.gz" --levels OncotreeLineage
                 --dump-specific 0.10 --out "$OUT/atlas.s3.analytic") ;;
  regress) args=(--signal "$W/results/atlas.s3.se_signal.tsv.gz" --catalog "$W/results/atlas.s3.union_catalog.bed.gz"
                 --pull-set "$REPO/phase2/data/pull_set.tsv" --norm none --levels OncotreeLineage
                 --out "$OUT/regress_v1") ;;
esac
# LEVELS overrides the arm's levels (argparse keeps the last --levels): e.g. NPERM=10000 KEEPFRAC=0.005
# LEVELS=OncotreeSubtype OUTDIR=... re-runs one level off the permutation floor (FINDINGS §37)
[ -n "${LEVELS:-}" ] && args+=(--levels "$LEVELS")
if [ -n "$NOCN" ]; then
  for k in "${!args[@]}"; do [ "${args[$k]}" = "--out" ] && args[$((k + 1))]="${args[$((k + 1))]}.nocn"; done
  args+=(--no-cn)
fi
echo "[arm] $ARM on $(hostname) $(date -Is)"
"$PY" "$REPO/phase2/score_pilot.py" "${args[@]}" > "$OUT/$ARM.stdout" 2> "$OUT/$ARM.stderr"
echo "[arm] $ARM done $(date -Is)"
