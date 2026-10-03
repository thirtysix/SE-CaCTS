#!/usr/bin/env bash
# Stage the v3.1 dashboard release (the fused build, FINDINGS §53-55) into a copy of docs/ for review; nothing in
# docs/ changes. Publishing = copying <staging-dir>/data over docs/data (and docs/js, docs/index.html) after review.
#   bash phase2/scripts/76_stage_v31_release.sh <staging-dir>
# Needs phase2/results_v31f, phase2/scores_v31f (main, nocn, noinf, l* arms, out_all/, fused_labels.*, the
# cn_ablation_calls and concordance2 tables) and phase2/data/pull_set.v31.tsv. One staging worker by default: the
# laptop reaches 94 C with two (nextsession 2026-09-30).
set -euo pipefail
SECACTS="$(cd "$(dirname "$0")/../.." && pwd)"
D="${1:?staging dir}"
PY="${PY:-$HOME/miniconda3/envs/atac_hdac/bin/python}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# STEPS="60 75" restages into an existing copy (e.g. after a release-notes edit); the default rebuilds the copy
STEPS="${STEPS:-}"
if [ -z "$STEPS" ]; then
  rm -rf "$D"; cp -r "$SECACTS/docs" "$D"; rm -rf "$D/data/lines" "$D"/data/calls_*.tsv "$D"/data/rank_*.tsv "$D"/data/gene_index*.json
  STEPS="60 61 75"
fi
step () { [[ " $STEPS " == *" $1 "* ]]; }
S="$SECACTS/phase2"; SC="$S/scores_v31f"; RES="$S/results_v31f"; PS="$S/data/pull_set.v31.tsv"
TITLE="Copy-number-aware calling, labelled; subtype calls"
NOTES="Super-enhancers are called with a cutoff set on copy-number-corrected signal, so an amplicon no longer raises the bar for every other locus, and SEs that need their amplification to pass are kept and labelled (gain, amplified) instead of dropped. A specific call needs an SE of its group, and the permutation null applies the same rule. Every call is scored with and without copy-number correction: CN-robust calls pass both, CN-unmasked only with correction; calls that pass only without it are listed beside them as amplicon-driven or gain-dependent. Subtypes now have calls; a group of one cell line lists none (its specificity is that line's own, tested at the cell-line level). Also: experiments relabelled (73 no longer baseline), DepMap WGS first wherever it exists, 7 misidentified lines dropped, one merge rule for the catalogue, chrY not scored, and 11 cell lines whose names carry punctuation get their per-line calls back. Lineage calls 15,324 -> 13,754, disease 13,950 -> 15,699, subtype 10,971 in 58 subtypes."
step 60 && SECACTS_DOCS="$D" "$PY" "$S/scripts/60_stage_dashboard.py" --scores "$SC" --results "$RES" \
  --pull-set "$PS" --release v3.1 --release-date "${RELEASE_DATE:-$(date +%F)}" \
  --release-title "$TITLE" --release-notes "$NOTES" --subtype-calls --min-group-lines 2 --line-prefix atlas.s3.lines.all --labels "$SC/fused_labels.groups.tsv.gz" \
  --n-pull "$(grep -c . "$S/data/pull_srx.v31.txt")" \
  --pull-desc "ChIP-Atlas hg38, plus SRA data it lacks run through its v1 pipeline: | 2,188 untreated or control experiments on human cancer cell lines" \
  --main-label "Default: CN-corrected, all 510 lines" \
  --main-desc "Untreated and control experiments; copy-number-corrected signal; every scored line." \
  --variant "noinf|atlas.s3.perm.noinf|Measured copy number only (438 lines)|The 72 lines whose copy number is inferred from ChIP input are left out (lineage and disease only): that inference corrects less than measured copy number." \
  --variant "nocn|atlas.s3.perm.nocn|No copy-number correction|The same experiments and lines scored without dividing by copy number, so amplicons can look like specificity (see the CN ablation tab)." \
  --variant "all|out_all/atlas.s3.perm|All experiments, incl. treated (3,058)|Drug-treated, knocked-down and otherwise perturbed experiments added back (ChIP-Atlas only), on the same super-enhancer loci; copy-number-corrected."
step 61 && SECACTS_DOCS="$D" SECACTS_RES="$RES" SECACTS_SC="$SC" SECACTS_SC_LINES="$SC" SECACTS_PS="$PS" \
  SECACTS_PRES=atlas.s3.se_presence.fu.tsv.gz SECACTS_LABELS="$SC" SECACTS_FLAGS=0 \
  STAGE_WORKERS="${STAGE_WORKERS:-1}" "$PY" "$S/scripts/61_stage_lines.py" --lines "${LINES:-all}"
du -sh "$D/data/lines"
step 75 && "$PY" "$S/scripts/75_stage_se_genes.py" --docs "$D" --catalog "$RES/atlas.s3.union_catalog.bed.gz" \
  --signal "$RES/atlas.s3.se_signal.tsv.gz" --pull-set "$PS" --pairs "$SC/atlas.s3.perm.concordance2.pairs.tsv.gz"
echo "[76] staged -> $D"
