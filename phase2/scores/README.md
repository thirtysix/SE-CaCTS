# phase2/scores — what's canonical, what's void

Outputs of `phase2/score_pilot.py` and the `phase2/analysis/` scripts. **Read `../../RESULTS.md` for the
claims.** This file is only about which artifact is which.

## Naming convention (the one thing to get right)

| prefix | FDR basis | use it? |
|---|---|---|
| **`atlas.s3.perm.*`** | **label-permutation null** | **CANONICAL.** Counts and calls come from here. |
| `atlas.s3.*` (no `perm`, no `nocn`) | analytic normal-approx null | **rankings valid, COUNTS VOID** (gotcha 71) |
| `atlas.s3.nocn.*` | analytic, CN correction OFF | the uncorrected arm of the CN ablation only |

**Rankings are identical between `atlas.s3.*` and `atlas.s3.perm.*`** — same JSD, only the FDR column
differs (rankings are invariant to the null, gotcha 28). So `atlas.s3.*.top_specific.tsv` is safe to read
for *order* (it feeds the CN ablation as the corrected arm); never read its `n_spec_*` counts or its
`.specific.tsv.gz` dump as calls.

Both the canonical (`perm`) and the analytic runs are **CN-corrected** (scoring-time correction is the
default; only `nocn` disables it).

## Per-file

**Canonical (permutation):**
- `atlas.s3.perm.hierarchy_summary.tsv` — per-group specific-SE counts. **The count table to cite.**
- `atlas.s3.perm.{OncotreeLineage,OncotreePrimaryDisease}.specific.tsv.gz` — every FDR≤0.10 call. These two
  levels are the ones the panel supports (gotcha 72).
- `atlas.s3.perm.OncotreeSubtype.specific.tsv.gz` — 1 call total (the panel does not support subtype; kept
  as evidence of that, not as a result).
- `atlas.s3.perm.{level}.top_specific.tsv` — top-15 ranked SEs per group (rankings, all levels incl. line).
- `atlas.s3.perm.concordance2.{pairs.tsv.gz,summary.tsv}` — Phase-6 cross-layer bridge on the perm calls.

**CN ablation (rank-based, null-invariant):**
- `atlas.s3.cn_ablation.tsv` + `.log` — amplicon-driven calls (the MYCN/MYC story).
- `atlas.s3.nocn.*` — the uncorrected arm the ablation compares against.
- `atlas.s3.{level}.top_specific.tsv` — the corrected arm's rankings (analytic run, but rankings are valid).

**Analytic, retained for the calibration comparison only:**
- `atlas.s3.hierarchy_summary.tsv` — its `n_spec_*` columns are the VOID counts; its `genes`/rankings are
  fine. Kept because `RESULTS.md §5` and the calibration story quote analytic-vs-permutation side by side.
- `atlas.s3.{concordance,concordance2}.summary.tsv` — the analytic concordance summaries (the "2.0× that
  became 4.1× under permutation" comparison).
- `*.run.log` — provenance.

## Regenerable, NOT in git (gitignored)

Pruned 2026-07-23 (109 MB → 3 MB). All rebuildable:
- `*.line.specific.tsv.gz` — line level is rankings-only, so a "call" dump there is meaningless.
- `atlas.s3.Oncotree*.specific.tsv.gz` (analytic, no `perm`) — the VOID analytic count dumps.
- `atlas.s3.concordance2.pairs.tsv.gz` — the analytic stage-2 raw pairs (summary is kept).

Rebuild with:
```bash
python score_pilot.py --signal ../results/atlas.s3.se_signal.tsv.gz \
    --catalog ../results/atlas.s3.union_catalog.bed.gz --norm none \
    --fdr-method permutation --n-perm 1000 --dump-specific 0.10 --out scores/atlas.s3.perm
```
