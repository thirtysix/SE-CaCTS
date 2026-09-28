# SE-CaCTS — defensible results

**As of 2026-09-28 (v2 atlas).** This is the claims document: what the project can currently assert, at
what resolution, and with what caveats. It is deliberately narrower than the raw outputs.

> **The canonical scoring run is the PERMUTATION one on the v2 atlas** — `phase2/scores_v2/atlas.s3.perm.*`
> (`score_pilot.py --fdr-method permutation --n-perm 1000`). The 282-line v1 atlas (`phase2/scores/`) is
> kept for comparison; §6 shows 97% of its calls reappear in v2. Analytic-null outputs are retained only
> as calibration evidence: **their counts are not usable** (§2).

---

## 1. The atlas

| | v2 (canonical) | v1 |
|---|---|---|
| Source | ChIP-Atlas hg38 H3K27ac, QC-pass human cancer lines **with measured copy number** | same, DepMap WGS lines only |
| Pull | 3,468 experiments (v1 2,916 + 552 new), ~53 BU in total on CSC Roihu | 2,916 |
| Agnostic atlas | 3,468 samples × 48,180 SE loci | 2,916 × 43,931 |
| S3norm atlas (scored) | **2,563 samples × 47,101 SE loci → 386 cell lines** | 2,136 × 42,943 → 282 |
| Copy number | DepMap WGS 282 · CMP WES pureCN 91 · DepMap WES 13 | DepMap WGS 282 |
| Hierarchy | 24 lineages / 46 primary diseases / 84 subtypes | 24 / 44 / 75 |
| Reconstruction | grid→SE `max|err| = 0` | same |

v2 adds the 122 lines that became copy-number-correctable through whole-exome sources
(`phase1/scripts/15_expansion_set.py`); 104 pass the `--min-peaks 2000` QC gate. It reuses v1's caller
code, quantification grid and S3norm reference, so v1's samples normalize bit-identically and v1 and v2
differ only by the added samples. The fixed grid holds 96.4% of the new samples' peak territory; 124
union loci fall outside it entirely and are not scored (46,977 scored loci).

## 2. Calibration: which FDR can be trusted

Run the whole procedure on **shuffled labels**, where nothing real exists to find:

| null | real labels (lineage) | shuffled labels (lineage) | shuffled, all three levels |
|---|---:|---:|---:|
| normal approximation | 82,887 (7.35%) | **58,241 (5.17%)** | — |
| **label permutation, B = 1,000** | 12,994 (1.15%) | **0** | **0 of 7.25 M tests** |

The normal-approximation null calls almost as much on noise as on data (v1: 6.05% vs 7.4%), so every
count here is from the permutation null. The failure reproduces across panels, so it is a property of
that null, not of one dataset (gotchas 70–71).

## 3. What resolution the panel supports

**This is the most important limitation and it should lead any write-up.**

| level | groups | groups with ≥1 call | calls (v2) | v1 |
|---|---:|---:|---:|---:|
| OncotreeLineage | 24 | **23/24** | **12,994** | 6,790 (23/24) |
| OncotreePrimaryDisease | 46 | **43/46** | **11,652** | 4,343 (41/44) |
| OncotreeSubtype | 84 | 5/84 | 10 | 1 (1/75) |
| line | 386 | — | — | permutation is degenerate |

**Report lineage and primary disease. Treat subtype and cell-line level as RANKINGS ONLY.** 27 of 84
subtypes hold a single cell line and 62 hold four or fewer; because the permutation keeps group sizes, a
random handful of lines scores as extreme as the real grouping. Enlarging the panel by 37% nearly tripled
disease-level calls but moved subtype from 1 call to 10. Subtype resolution needs **more lines per
subtype**, and measured copy-number sources are close to exhausted (`phase1/CN_COVERAGE.md §5b`).

## 4. Known biology recovered without supervision

**Pre-specified panel.** 22 master regulators chosen from the literature before looking (12 lineages).
For each, the lineage-best SE within 100 kb of the gene:

| | own lineage | other lineages |
|---|---:|---:|
| passes FDR ≤ 0.10 | **18/22 (82%)** | 23/506 (4.5%) |

Odds ratio 94.5, Fisher p = 1e-18. Passing, with rank in the lineage: MECOM #1, SOX17 #5 (Ovary);
HNF4A #6, CDX2 #28 (Bowel); ESR1 #34, FOXA1 #112 (Breast); SPI1 #7, CEBPA #32 (Myeloid); IKZF1 #38,
PAX5 #151 at FDR 6e-4 (Lymphoid); SOX10, MITF (Skin); PHOX2B #3, HAND2 #27 (PNS); NKX2-1, ASCL1 (Lung);
HNF1A (Liver); AR (Prostate). **Not passing:** PAX8 (#123, FDR 0.16), GATA3 (#1,115), PAX2, TP63.
GATA3 instead passes in the peripheral nervous system, where it belongs to the neuroblastoma core
regulatory circuit. Figure: `phase2/figures/poster_figures.py identity`.

**Stable from v1 to v2** (primary disease): Ovarian Epithelial Tumor MECOM #1 → #1, SOX17 #5 → #5;
Colorectal HNF4A #6 → #7; AML SPI1 #19 → #20, CEBPA #25 → #10, IRF8 #18 → #12; Breast ESR1 #40 → #41.
FDR at the head of a disease group sits on a BH plateau (0.053 for the top 10 of most groups), so rank
order within a plateau is not meaningful (gotcha 28).

**Rankings that are informative but NOT callable** (line level, v1): MCF7 → ESR1 #1, THP-1 → CEBPA #3,
MOLM-13 → IRF8 #2, SKOV3 → MECOM #1, SW48 → CDX2 #5, P12-ICHIKAWA → LEF1 #6.

**Negative controls** (v1, rank-based): six triple-negative breast lines bury ESR1 at ranks
11,000–21,000; lobular carcinoma (near-always ER+) gives ESR1 #22 and FOXA1 #37. The method was never
told which lines were ER+.

**Caveat.** The top SE in a group is often not a recognizable identity gene, gene assignment is
proximity only (gotcha 22), and nearby rows can tile one SE domain.

## 5. Copy-number correction

Call-based ablation under the permutation null (`--no-cn`, then `cn_ablation_calls.py`):

| level | uncorrected | corrected | removed (of which CN > 1.3) | rescued | stable |
|---|---:|---:|---:|---:|---:|
| OncotreeLineage | 10,142 | **12,994** | 699 (456) | 3,551 | 9,443 |
| OncotreePrimaryDisease | 5,934 | **11,652** | 316 (275) | 6,034 | 5,618 |

1. **Removed calls sit on amplicons.** 1,015 in total, median CN 1.66, 72% at CN > 1.3, led by recurrent
   lineage amplicons found with no gene list: **MYCN** at 71× in neuroblastoma with co-amplified
   **DDX1** and **CYRIA** (2p24), **MYC** and **POU5F1B** (8q24) in embryonal tumours, CNS and pleura,
   **OTX2** in embryonal tumours.
2. **Rescued calls are copy-neutral.** 9,585, median CN 1.00, 2% at CN > 1.3. Amplicon spikes in
   individual lines inflate the permutation null's left tail; removing them tightens the null and lets
   genuine specificity through.

**The size of the rescue depends on the panel and the threshold.** v1 reported lineage calls rising
67 → 6,790. That reproduces exactly with today's code, but it is a BH threshold effect: v1's uncorrected
p-values sat just above the bar, and at FDR ≤ 0.25 the same two arms give 14,278 vs 17,149 (+20%). On v2
correction adds 28% at lineage level and doubles disease-level calls. The direction holds at every
threshold on both panels; the 100-fold figure does not and should not be quoted.

**Per-line rank flips** (line level, v1, rankings only): SK-N-BE(2), KELLY and NB1643 lose MYCN SEs at
177–215× from ranks #1–#14; COLO320 loses POU5F1B (8q24, 120×) from #1. Correction is bidirectional: MCF7
ESR1 #5 → #1 and P12-ICHIKAWA LEF1 #598 → #6 improve, Bowel CDX2 #1 → #14 is demoted. **OVCAR3**:
uncorrected, 13 of 15 top calls are the 19q13 amplicon (CN 5.6–9.4×); corrected, 3 of 15.

**MECOM is real, not an amplicon.** Rank #1 in both arms and both panels; in SKOV3 the locus is
CN-neutral (1.055) while MECOM SEs still take ranks #1–#6.

## 6. Robustness

| check | lineage calls kept | disease calls kept |
|---|---:|---:|
| v1 calls reappearing in v2 (overlapping locus, same group) | 97.4% | 97.8% |
| drop the 7 cross-source-discordant lines + MDA-MB-231 | 98.3% (Jaccard 0.97) | 95.8% (0.95) |
| correct 228 lines with CMP WES instead of DepMap WGS | 88.7% (0.87) | 84.8% (0.83) |

**The copy-number source is the main remaining sensitivity.** On the 228 lines with both sources, scored
against each (`cn_source_paired.py`), CMP WES copy number captures less of the signal's copy-number
dependence (per-line Spearman of SE signal vs CN: +0.089 vs +0.146 raw, −0.104 vs −0.060 corrected;
paired p ≈ 1e-31). At SE resolution the two sources agree only moderately (median r = 0.61 of log2 CN),
and CMP's copy number is compressed (SD of log2 CN 0.41 vs 0.68). The shift is largest in Bowel, Myeloid
and Lymphoid. So ~11–15% of calls depend on which measured source a line was corrected with, and lines
corrected with exome copy number carry a slightly different correction than WGS lines.

## 7. Cross-layer validation (Phase 6)

Genes near group-specific SEs are themselves specific to that group, scored by CaCTS on DepMap
expression over the same lines (359 with expression) and groups:

| set | per-pair | background | enrichment | per-SE-any | nearest-gene | shuffled |
|---|---:|---:|---:|---:|---:|---:|
| **v2, lineage** | **16.0%** | 4.27% | **3.7×** | 29.9% | 22.9% | 4.4% |
| **v2, disease** | **15.8%** | 3.04% | **5.2×** | 29.3% | 23.4% | 3.0% |
| v1, lineage | 18.0% | 4.36% | 4.1× | 34.2% | 27.3% | — |
| v1, disease | 18.1% | 3.07% | 5.9× | 33.2% | 27.4% | — |
| v1 analytic, lineage | 8.7% | 4.36% | 2.0× | 18.3% | 12.6% | — |

Controls behave as a local regulatory link must: the **group shuffle** sits at background; concordance
**decays with distance**, 25.8% at < 10 kb → 8.4% at 100–250 kb, with median rho tracking it
(+0.305 → +0.143); and SE signal vs neighbour expression is higher for concordant than discordant pairs
(+0.314 vs +0.171). On v1 the concordance roughly **doubled** on the permutation-filtered set versus the
analytic one, independent evidence that the permutation FDR removes noise rather than signal.

## 8. EMX2 — an honest partial result (v1)

`USE_6049` (v1 id; chr10:117,543,636–117,545,605) sits on the **EMX2 promoter**. Rank #4 for OVCAR3 in
both CN arms; against independent DepMap RNA, rho = **+0.461** (p = 3e-16) for EMX2 vs +0.046 PAX8 and
+0.016 WT1 as controls. **It does not pass the specificity bar** (HGSOC subtype rank #4, FDR 0.173, and
subtype level is unsupported), and it was called as an SE in only 1 of 2,136 experiments. "A specific
H3K27ac element at the EMX2 promoter" is the accurate description.

## 9. What is NOT claimed

- Any specific-SE **count** from the analytic null (§2).
- **Subtype- or cell-line-level** specificity calls (§3).
- That the call set is **independent of the copy-number source**: ~11–15% of calls move with it (§6).
- The **"67 → 6,790" magnitude** of the v1 rescue (§5); only its direction.
- **SE → target-gene assignment.** Proximity only; §7 measures concordance in aggregate.
- Anything from the CN-corrected **calling-time** atlases; they were never built (gotcha 59).

## 10. Reproducing (v2)

```bash
conda activate atac_hdac
python phase1/scripts/15_expansion_set.py                  # the 122-line / 552-experiment expansion
# pull + reduce on Roihu: array.slurm with MANIFEST=pull_srx.expand.txt, then
#   reduce.slurm with MANIFEST=pull_srx.v2.txt CATALOGS="atlas atlas.s3" S3_REF=SRX16495452
ARMS="main nocn shuffle nodisc" bash phase2/scripts/70_score_v2.sh   # ~35 min per arm (run on HPC)
python phase2/analysis/cn_ablation_calls.py --corrected phase2/scores_v2/atlas.s3.perm \
    --uncorrected phase2/scores_v2/atlas.s3.perm.nocn \
    --catalog phase2/results_v2/atlas.s3.union_catalog.bed.gz \
    --out phase2/scores_v2/atlas.s3.perm.cn_ablation_calls.tsv
python phase2/analysis/concordance_bridge2.py --scores phase2/scores_v2/atlas.s3.perm \
    --signal phase2/results_v2/atlas.s3.se_signal.tsv.gz \
    --catalog phase2/results_v2/atlas.s3.union_catalog.bed.gz \
    --pull-set phase2/data/pull_set.v2.tsv --levels OncotreeLineage,OncotreePrimaryDisease \
    --out phase2/scores_v2/atlas.s3.perm.concordance2
python phase2/analysis/cn_source_paired.py
python phase2/analysis/v2_compare.py --a phase2/scores/atlas.s3.perm \
    --a-catalog phase2/results/atlas.s3.union_catalog.bed.gz --b phase2/scores_v2/atlas.s3.perm \
    --b-catalog phase2/results_v2/atlas.s3.union_catalog.bed.gz --out phase2/scores_v2/compare_v1_v2
python phase2/scripts/60_stage_dashboard.py --scores phase2/scores_v2 --results phase2/results_v2 \
    --pull-set phase2/data/pull_set.v2.tsv --pull-bu 53 \
    --release v2 --release-date 2026-09-28 --release-title "Copy-number-expanded panel: 386 cell lines"
python phase2/figures/poster_figures.py
```
