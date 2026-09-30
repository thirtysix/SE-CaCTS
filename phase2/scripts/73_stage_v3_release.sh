#!/usr/bin/env bash
# Stage the v3 dashboard release into a copy of docs/ for review; nothing in docs/ changes.
#   bash phase2/scripts/73_stage_v3_release.sh <staging-dir>
# Needs phase2/results_v3, phase2/scores_v3 (with out_lines/, the l* arms of phase2/roihu/score_arm.slurm) and
# phase2/data/pull_set.v3.tsv. Publishing = copying <staging-dir>/data over docs/data after review.
set -euo pipefail
SECACTS="$(cd "$(dirname "$0")/../.." && pwd)"
D="${1:?staging dir}"
PY="${PY:-$HOME/miniconda3/envs/atac_hdac/bin/python}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1        # 8 staging processes x 1 thread
rm -rf "$D"; cp -r "$SECACTS/docs" "$D"; rm -rf "$D/data/lines"
S="$SECACTS/phase2"
SECACTS_DOCS="$D" "$PY" "$S/scripts/60_stage_dashboard.py" --scores "$S/scores_v3" --results "$S/results_v3" \
  --pull-set "$S/data/pull_set.v3.tsv" --release v3 --release-date 2026-09-30 \
  --release-title "Baseline-only atlas: untreated and control experiments" \
  --release-notes "Perturbed experiments, input controls and other marks removed atlas-wide; own processing of SRA data ChIP-Atlas lacks (Stage B); CN from CCLE SNP6 and inferred from ChIP input admits more lines. Every cell line compared four ways (all lines, same subtype, disease, lineage); the Line view is now the Genomic View (IGV), for all lines." \
  --pull-desc "ChIP-Atlas hg38, plus SRA data it lacks run through its v1 pipeline: | 2,422 untreated or control experiments on human cancer cell lines" \
  --main-label "Default: CN-corrected, all 519 lines" \
  --main-desc "Untreated and control experiments; copy-number-corrected signal; every scored line." \
  --variant "noinf|atlas.s3.perm.noinf|Measured copy number only (439 lines)|The 80 lines whose copy number is inferred from ChIP input are left out: that inference corrects less than measured copy number." \
  --variant "nocn|atlas.s3.perm.nocn|No copy-number correction|The same experiments and lines scored without dividing by copy number, so amplicons can look like specificity (see the CN ablation tab)."
SECACTS_DOCS="$D" SECACTS_RES="$S/results_v3" SECACTS_SC="$S/scores_v3" SECACTS_PS="$S/data/pull_set.v3.tsv" \
  STAGE_WORKERS="${STAGE_WORKERS:-8}" "$PY" "$S/scripts/61_stage_lines.py" --lines all
du -sh "$D/data/lines"
echo "[73] staged -> $D"
