#!/usr/bin/env bash
# Stage the v3 dashboard release into a copy of docs/ for review; nothing in docs/ changes.
#   bash phase2/scripts/73_stage_v3_release.sh <staging-dir>
# v3.0.1 (specific SEs must be SEs of their group, 74_called_filter.py):
#   SCORES=phase2/scores_v3c RELEASE=v3.0.1 bash phase2/scripts/73_stage_v3_release.sh <staging-dir>
# Needs phase2/results_v3, phase2/scores_v3 (with out_lines/, the l* arms of phase2/roihu/score_arm.slurm) and
# phase2/data/pull_set.v3.tsv. Publishing = copying <staging-dir>/data over docs/data after review.
set -euo pipefail
SECACTS="$(cd "$(dirname "$0")/../.." && pwd)"
D="${1:?staging dir}"
PY="${PY:-$HOME/miniconda3/envs/atac_hdac/bin/python}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1        # 8 staging processes x 1 thread
rm -rf "$D"; cp -r "$SECACTS/docs" "$D"; rm -rf "$D/data/lines"
S="$SECACTS/phase2"
SC="$SECACTS/${SCORES:-phase2/scores_v3}"
REL="${RELEASE:-v3}"
if [[ "$REL" == "v3" ]]; then
  TITLE="Baseline-only atlas: untreated and control experiments"
  NOTES="Perturbed experiments, input controls and other marks removed atlas-wide; own processing of SRA data ChIP-Atlas lacks (Stage B); CN from CCLE SNP6 and inferred from ChIP input admits more lines. Every cell line compared four ways (all lines, same subtype, disease, lineage); the Line view is now the Genomic View (IGV), for all lines."
else
  TITLE="Specific super-enhancers must be super-enhancers of their group"
  NOTES="A specific call now needs at least one experiment of the lineage, disease, subtype or line to call a super-enhancer overlapping the locus. Before, a locus counted when the group had the most H3K27ac there even if none of its experiments called it a super-enhancer (18% of lineage calls, 74% of per-line calls). Same experiments, scores and FDRs; lineage calls 18,639 -> 15,324, disease 17,031 -> 13,950, per line vs all 601,878 -> 154,343."
fi
SECACTS_DOCS="$D" "$PY" "$S/scripts/60_stage_dashboard.py" --scores "$SC" --results "$S/results_v3" \
  --pull-set "$S/data/pull_set.v3.tsv" --release "$REL" --release-date 2026-09-30 \
  --release-title "$TITLE" --release-notes "$NOTES" \
  --pull-desc "ChIP-Atlas hg38, plus SRA data it lacks run through its v1 pipeline: | 2,422 untreated or control experiments on human cancer cell lines" \
  --main-label "Default: CN-corrected, all 519 lines" \
  --main-desc "Untreated and control experiments; copy-number-corrected signal; every scored line." \
  --variant "noinf|atlas.s3.perm.noinf|Measured copy number only (439 lines)|The 80 lines whose copy number is inferred from ChIP input are left out: that inference corrects less than measured copy number." \
  --variant "nocn|atlas.s3.perm.nocn|No copy-number correction|The same experiments and lines scored without dividing by copy number, so amplicons can look like specificity (see the CN ablation tab)."
SECACTS_DOCS="$D" SECACTS_RES="$S/results_v3" SECACTS_SC="$SC" SECACTS_PS="$S/data/pull_set.v3.tsv" \
  STAGE_WORKERS="${STAGE_WORKERS:-8}" "$PY" "$S/scripts/61_stage_lines.py" --lines all
du -sh "$D/data/lines"
echo "[73] staged -> $D"
