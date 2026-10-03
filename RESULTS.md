# SE-CaCTS — defensible results

**As of 2026-10-03 (v3.1, https://doi.org/10.5281/zenodo.23124053).** This is the claims document: what the
project can currently assert, at what resolution, and with what caveats. It is deliberately narrower than the raw
outputs. The claims document for the v2 atlas is in the git history (`git show c18acf5:RESULTS.md`).

> **The canonical run is the fused build:** `phase2/scores_v31f/atlas.s3.perm.*` on `phase2/results_v31f/`
> (`score_pilot.py --fdr-method permutation --n-perm 1000`, global BH, FDR ≤ 0.10, the "super-enhancer of its
> group" rule inside the permutation null). Counts below are at FDR ≤ 0.10 unless stated.

---

## 1. The atlas

| | v3.1 | v3.0.1 |
|---|---|---|
| Source | ChIP-Atlas hg38 H3K27ac, untreated or control experiments on human cancer lines, plus SRA data it lacks run through its v1 pipeline | same |
| Experiments | 2,188 (each study's context read; 73 v3 experiments no longer count as baseline) | 2,422 |
| Scored (QC pass) | **1,756 experiments → 510 cell lines** | 1,945 → 519 |
| Catalogue | **33,561 loci** (agnostic + fused calls, merged on any overlap); 33,255 scored | 46,443 |
| Copy number | DepMap WGS 312 · CMP WES 103 · DepMap WES 20 · CCLE SNP6 3 · inferred from ChIP input 72 | 282 · 133 · 21 · 3 · 80 |
| Hierarchy | 25 lineages / 56 primary diseases / 103 subtypes | 25 / 56 / 104 |

What changed in v3.1: DepMap WGS copy number is used wherever it exists, 7 misidentified lines are dropped, chrY is
not scored, and calling is copy-number-aware (§5). Each experiment's ROSE cutoff is set on copy-number-corrected
signal and applied to the uncorrected signal, so an amplicon no longer raises the bar for every other locus, and an
SE that passes only through its extra copies is kept and labelled instead of dropped. Of 1,029,900 experiment × locus
calls, 68.6% are SEs with or without correction (core), 12.1% only with it (unmasked), and 18.9% need extra copies:
gain (CN < 2) 15.6%, amplified (2–3) 2.1%, high-level (≥ 3) 1.2%.

## 2. Calibration

Run the whole procedure on **shuffled labels**, where nothing real exists to find: **0 calls at every level** (lineage,
disease, subtype, subtype consensus, within-lineage subtype, each of the four per-line comparisons, and the
all-experiments variant). The normal-approximation null called 6.05% of tests on shuffled labels, so every count
here comes from the permutation null.

Two parts of the procedure exist because a simpler version failed this test. **Multiple testing is corrected over
all groups together:** correcting each group separately called 672 SEs on shuffled subtypes. **The rule "a specific
SE must be an SE of its group" is applied inside the null** (each shuffled group is tested only where its own
experiments call an SE): applying it afterwards, without conditioning the null, gave 155 false calls on shuffled
labels in a simulation.

## 3. What resolution the panel supports

| level | groups | groups with calls | calls |
|---|---:|---:|---:|
| OncotreeLineage | 25 | 23 | **13,754** |
| OncotreePrimaryDisease | 56 | 35 | **15,699** |
| OncotreeSubtype | 103 | 58 | **10,971** |
| cell line vs all lines | 510 | 215 | 78,745 |
| cell line vs its lineage / disease / subtype | 510 | 207 / 190 / 163 | 30,627 / 15,628 / 8,282 |

**A group of one cell line lists no calls.** Its "group" specificity is that line's own, which the cell-line level
tests against a stricter standard (two independent studies). Counting them would add 285 disease and 171 subtype
calls (raw: 15,984 in 47 diseases, 11,142 in 83 subtypes). 41 of 103 subtypes hold one line and 73 hold four or
fewer. **Cell lines are tested only when they have two independent studies**, which is why 215 of 510 have calls.

**Subtypes have calls now, and they are stable.** In v3.0.1 subtype counts swung with how BH pooled the tests, so
subtype was shown as rankings only. With the rule inside the null, BH no longer pools the many loci a subtype never
calls (on the v3.1 build before fusion the rule cut lineage-level tests to about 15%). On that build, removing the 72
inferred-copy-number lines or 72 lineage-matched random lines left 7,358 and 7,105 subtype calls.

## 4. Known biology recovered without supervision

**Pre-specified panel.** 22 master regulators chosen from the literature before looking (12 lineages). For each, the
best-ranked lineage call within 100 kb of the gene:

| | own lineage | other lineages |
|---|---:|---:|
| passes FDR ≤ 0.10 | **21/22** | 1.1% |

Ranks in the lineage: MECOM #1, SOX17 #12 (Ovary); CDX2 #52, HNF4A #336 (Bowel); ESR1 #32, FOXA1 #50, GATA3 #343
(Breast); SPI1 #26, CEBPA #602 (Myeloid); IKZF1 #1, PAX5 #1 (Lymphoid); SOX10 #176, MITF #124 (Skin); PHOX2B #1,
HAND2 #2 (PNS); NKX2-1 #18, ASCL1 #38 (Lung); HNF1A #79 (Liver); AR #82 (Prostate); PAX2 #177 (Kidney); TP63 #21
(Head and Neck). **Not passing:** PAX8 (#162, FDR 0.13). v3.0.1 found 19/22 (PAX8, HNF4A and GATA3 missing); the
uncorrected statistic of v3.1 finds 20/22 (PAX8, HNF4A).

**AR shows why calling had to change.** Its SE is a CN-robust prostate call, but every prostate experiment that calls
an SE there is in an AR-amplified line (11 of 11). Strict calling-time correction removed all of those calls, so on
the v3.1 build before fusion the locus was not a test for Prostate at all.

**Caveat.** The top SE in a group is often not a recognizable identity gene, gene assignment is proximity only, and
nearby rows can tile one SE domain.

## 5. Copy number

Every call is scored twice over the same tests, with and without copy-number correction. The two statistics are never
merged into one FDR:

| level | corrected calls | CN-robust (both) | CN-unmasked (corrected only) | uncorrected only: amplicon-driven (CN ≥ 2) / gain-dependent |
|---|---:|---:|---:|---:|
| lineage | 13,754 | 11,546 | 2,208 | 127 / 677 |
| disease | 15,699 | 12,477 | 3,222 | 168 / 556 |
| subtype | 10,971 | 7,787 | 3,184 | 156 / 261 |

1. **CN-unmasked calls are copy-neutral** (median CN 1.00 at every level). Amplicons in individual lines inflate the
   permutation null; correcting them lets real specificity through.
2. **Calls that pass only without correction are listed beside the atlas, not counted.** The amplicon-driven ones are
   the recurrent lineage amplicons, found with no gene list: **MYCN** with **DDX1** and **NBAS** (2p24, CN 9–58) in
   neuroblastoma and the PNS, **MYC** and **POU5F1B** (8q24) in embryonal and CNS tumours, **OTX2** in embryonal
   tumours. Most are low-level gain instead (lineage: median CN 1.25, 16% at CN ≥ 2), hence the two labels.
3. **The 72 inferred-copy-number lines matter as panel size, not as copy number.** Without them: 11,213 lineage
   calls; without 72 lineage-matched random measured-CN lines instead: 11,248; keeping them but uncorrected: 13,444.

**The copy-number source is a known sensitivity.** On the v2 panel, correcting 228 lines with CMP WES instead of
DepMap WGS moved 11–15% of calls (exome copy number is compressed, SD of log2 CN 0.41 vs 0.68). v3.1 uses DepMap
WGS for every line that has it.

## 6. Robustness

| check | lineage | disease |
|---|---:|---:|
| v3.0.1 calls still called in v3.1 (same group, overlapping locus) | 78.4% | 80.2% |
| library-layout effect removed: calls kept | 94.1% | 94.1% |
| study-weighted line profiles: calls | 14,163 | 16,575 |
| 8 cross-source-discordant lines dropped: calls | 13,585 | 15,473 |
| all experiments, treated included (3,058): calls | 14,764 | 16,950 |

Disease counts and retention in this table include groups of one line (the main arm on that basis: 15,984).

**Library layout.** 75% of SEs read about 1.2× higher in single-end than in paired-end experiments of the same line
(142 lines with both; sign-flip null 0). Removing that within-line effect keeps 94.1% of lineage calls, but the
losses sit in lineages made mostly of single-end experiments (Spearman +0.78 between retention and paired-end share):
those under 25% paired-end keep 85%, the rest 95%. Lowest: Esophagus/Stomach 0.76, Head and Neck 0.79, Breast 0.79,
Pancreas 0.80. Layout is confounded with read length and sequencing era, so it names the batch, not the cause.

**Leave one study out.** For each lineage, its largest study (most lines) was removed and the lineage rescored,
against two removals of random studies of the same lineage covering as many lines. Other lineages barely move
(median retention 0.996). **Kidney and Bladder rest on one study:** Kidney keeps 11% of its calls without the NCI-60
renal panel (PRJNA601191; 6 of its 14 lines exist only there) against 63–66% for the random removals; Bladder keeps
14% against 94–98%. Weaker dependence: Bowel 0.48 (random 0.93–0.98), Esophagus/Stomach 0.53 (0.77–0.85), Skin 0.74
(0.81–0.93), Liver 0.74 (0.99–1.00). Lymphoid, Myeloid, PNS, Prostate and Soft Tissue keep 86–97% of their calls
and 95–100% of their top 100; Breast and Ovary keep more than their random controls (0.84 vs 0.10–0.23, 0.83 vs
0.60–0.70); Bone and Lung lose more of their tail (0.65, 0.58) but keep 99% of their top 100. Cervix, Testis, Eye,
Pleura and Uterus are too small to read.

## 7. Cross-layer validation

Genes within 100 kb of a group-specific SE are themselves specific to that group, scored by CaCTS on DepMap
expression over the same lines and groups:

| level | per pair | background | enrichment | SEs with ≥ 1 concordant gene | nearest gene | shuffled groups |
|---|---:|---:|---:|---:|---:|---:|
| lineage | **16.4%** | 4.11% | **4.0×** | 29.5% | 22.1% | 4.4% |
| disease | **14.3%** | 2.81% | **5.1×** | 26.2% | 19.7% | 3.1% |
| subtype | **16.6%** | 2.54% | **6.5×** | 29.1% | 22.9% | 2.9% |

The controls behave as a local regulatory link must: the **group shuffle** sits at background, and concordance
**decays with distance** (all levels pooled: 24.5% under 10 kb, 7.3% beyond 100 kb, while the shuffle stays at
3–4%). v3.0.1: 3.7× lineage, 5.1× disease.

## 8. What is NOT claimed

- Any specific-SE **count** from the analytic null (§2).
- Group-level calls for **groups of one cell line**; those are the cell-line level's (§3).
- That **amplicon-driven or gain-dependent** calls are lineage-specific; they pass only without correction (§5).
- That every lineage's calls are **independent of study**: Kidney and Bladder rest on one study (§6).
- That the calls are **free of the layout batch**: about 6% of lineage calls depend on it, concentrated in
  single-end lineages (§6).
- That the call set is **independent of the copy-number source** (§5); the inferred correction of 72 lines is weaker.
- **SE → target-gene assignment.** Proximity only; §7 measures concordance in aggregate.
- The **"67 → 6,790" magnitude** of the v1 copy-number rescue: it is a BH threshold effect (at FDR ≤ 0.25 the same
  two arms give 14,278 vs 17,149). Only its direction holds.

## 9. Reproducing (v3.1)

```bash
conda activate atac_hdac
python3 phase1/scripts/25_pull_set_v31.py                 # pull sets (relabelled, WGS first, identity drops)
# on Roihu, from score_v3/ (see each script's header):
#   fused re-call:  sbatch -A <project> --array=0-15 --export=ALL,PROJ=<project>,NT=16,MODE=fused,\
#                     MANIFEST=<repo>/phase2/data/pull_set.v31.all.tsv recall_cn.slurm
#   reduce + arms:  FUSED=1 PROJ=<project> bash v31_chain.sh      -> results_v31f/, out_v31f/
python3 phase2/analysis/fused_labels.py --scores phase2/scores_v31f --results phase2/results_v31f \
    --pull-set phase2/data/pull_set.v31.tsv                 # CN-robust / unmasked / amplicon-driven labels
bash phase2/scripts/76_stage_v31_release.sh <staging-dir>   # dashboard data (copy <staging-dir>/data over docs/data)
# robustness: layout (§6)
python3 phase2/analysis/srx_layout_fill.py --signal phase2/results_v31f/atlas.s3.se_signal.tsv.gz
python3 phase2/analysis/layout_batch.py --results phase2/results_v31f --pull-set phase2/data/pull_set.v31.tsv \
    --scores phase2/scores_v31f --tag v31f
python3 phase2/analysis/layout_correct.py --signal phase2/results_v31f/atlas.s3.se_signal.tsv.gz \
    --effects phase2/analysis/out/layout_sensitive_ses.v31f.tsv.gz --layout phase2/data/srx_layout.tsv \
    --out phase2/results_v31f/atlas.s3.layoutcorr.se_signal.tsv.gz      # then score_arm.slurm ARM=layout V31=1 FUSED=1
# robustness: leave one study out (§6)
python3 phase2/analysis/loso_sets.py --signal phase2/results_v31f/atlas.s3.se_signal.tsv.gz \
    --pull-set phase2/data/pull_set.v31.tsv --out phase2/analysis/out/loso_v31f   # then V31F=1 bash loso_submit.sh
python3 phase2/analysis/loso_eval.py --rule-in-scores --runs phase2/analysis/out/loso_v31f/runs \
    --sets phase2/analysis/out/loso_v31f --base phase2/scores_v31f/atlas.s3.perm.OncotreeLineage.specific.tsv.gz \
    --out phase2/analysis/out/loso_v31f/loso_retention.tsv
```
