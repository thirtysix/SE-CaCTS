# SE-CaCTS — defensible results

**As of 2026-07-28.** This is the claims document: what the project can currently assert, at what
resolution, and with what caveats. It is deliberately narrower than the raw outputs.

> **The canonical scoring run is the PERMUTATION one** — `phase2/scores/atlas.s3.perm.*`
> (`score_pilot.py --fdr-method permutation --n-perm 1000`). The analytic-null outputs
> (`atlas.s3.*`) are retained for comparison **but their specific-SE counts are not usable**: the
> analytic FDR fails an outright calibration test (gotcha 71). Rankings are shared between the two.

---

## 1. The atlas

| | |
|---|---|
| Source | ChIP-Atlas hg38 H3K27ac, QC-pass, DepMap-joined |
| Pull | 2,916 of 2,917 experiments (one retired SRX), zero OOM, ~43 BU on CSC Roihu |
| Agnostic atlas | **2,916 samples × 43,931 SE loci** |
| S3norm atlas (canonical) | **2,136 samples × 42,943 SE loci → 282 cell lines** |
| Hierarchy | 24 Oncotree lineages / 44 primary diseases / 75 subtypes |
| Reconstruction | grid→SE `max|err| = 0` on both |

The S3norm atlas applies the `--min-peaks 2000` QC gate, which drops 784 weak samples (27%) and 42
lines. That gate is **not optional** — ungated S3norm performs worse than quantile normalization
(gotcha 25).

## 2. What resolution the panel supports

**This is the most important limitation and it should lead any write-up.**

| level | groups | groups with ≥1 call | total calls | median/group |
|---|---|---|---|---|
| OncotreeLineage | 24 | **23/24** | 6,790 | 81.5 |
| OncotreePrimaryDisease | 44 | **41/44** | 4,343 | 15 |
| OncotreeSubtype | 75 | **0/75** (1 call total) | 1 | 0 |
| line | 282 | — | — | permutation is degenerate |

`median/group` is taken over **all** groups at that level; restricted to groups with ≥1 call it is 94
(lineage) and 18 (disease). All five columns recompute from
`phase2/scores/atlas.s3.perm.hierarchy_summary.tsv` (column `n_spec_fdr10`).

**Report lineage and primary disease. Treat subtype and cell-line level as RANKINGS ONLY, never as
calls.** Subtype fails because 29 of 75 subtypes contain a single cell line and 56 contain ≤4; since the
permutation null preserves group size, a random handful of lines produces JSD as extreme as the real
grouping. Calls by subtype group size: n≤6 → 0, n=7 → 1. For contrast, at lineage level groups with 15+
lines contribute 4,064 calls and groups with 1–4 lines contribute 200. Subtype resolution needs **more
lines per subtype**, not more samples per line — i.e. the scope expansion toward ~722 lines.

## 3. Known biology recovered without supervision

All of the following **pass the permutation FDR ≤ 0.10** at the stated level (FDR 0.054–0.096 — passing,
but not by a wide margin):

| level | group | n | gene | SE rank |
|---|---|---|---|---|
| Lineage | Ovary/Fallopian Tube | 9 | **MECOM** | #1, #2, #8 |
| Lineage | Ovary/Fallopian Tube | 9 | **SOX17** | #6 |
| Lineage | Bowel | 15 | **HNF4A** / **CDX2** | #4 / #14 |
| Lineage | Myeloid | 19 | **SPI1** | #9 |
| PrimaryDisease | Ovarian Epithelial Tumor | 9 | **MECOM** / **SOX17** | #1,#2,#8 / #5 |
| PrimaryDisease | Colorectal Adenocarcinoma | 15 | **HNF4A** / **CDX2** | #6 / #13 |

Nothing supplied these genes to the scorer; they are recovered from H3K27ac signal alone.

**Rankings that are informative but NOT callable** (line/subtype level, where the permutation is
degenerate or unsupported): MCF7 → ESR1 #1, THP-1 → CEBPA #3, MOLM-13 → IRF8 #2, SKOV3 → MECOM #1,
SW48 → CDX2 #5, P12-ICHIKAWA → LEF1 #6.

**Negative controls (rank-based, unaffected by the FDR question).** Six triple-negative breast lines bury
ESR1 at ranks 11,000–21,000; the ER-negative subtype group "Breast Invasive Carcinoma, NOS" puts it at
#19,059; lobular carcinoma — near-always ER+ — gives ESR1 #22 and FOXA1 #37. The method was never told
which lines were ER+.

**Caveat on the top-ranked SE per group.** The single most significant SE in a group is frequently *not*
a recognizable identity gene (e.g. Lymphoid → IGLL5, Myeloid → GOLGA8Q, Lung → CYRIA). BH creates FDR
plateaus (gotcha 28) so ordering within a plateau is unstable, and gene assignment is proximity-only
(gotcha 22). Identity genes appear in the passing set, not reliably at its head.

## 4. Copy-number correction

The CN layer is validated two ways: by the SEs it **removes** (rank-based, null-invariant) and by a
call-based ablation under the honest permutation null (`--no-cn` scored with `--fdr-method permutation`,
compared by `cn_ablation_calls.py`). The call-based one is the stronger statement.

**Call-based ablation (permutation null, the supported levels).** Comparing which SE×group tests pass the
permutation FDR ≤ 0.10 with vs without correction:

| level | uncorrected calls | corrected calls | amplicon-driven (removed) | rescued | stable |
|---|---|---|---|---|---|
| OncotreeLineage | 67 | 6,790 | 22 | 6,745 | 45 |
| OncotreePrimaryDisease | 11 | 4,343 | 9 | 4,341 | 2 |

Two findings, both clean:

1. **The amplicon-driven false calls are few, and every one is a named recurrent lineage amplicon.** 31
   calls total, `cn_mean` median **82.9**, 100 % at CN > 1.3: **MYCN** in Peripheral Nervous System +
   Neuroblastoma (cn 83, the defining neuroblastoma amplicon), **OTX2** in CNS/Brain (cn 7, medulloblastoma),
   **FGFR2** in Esophagus/Stomach (cn 12.5, gastric), **ANO1**/11q13 in Head & Neck SCC (cn 4). Called
   "specific" under the uncorrected null and correctly removed by correction — unprompted.
2. **The dominant effect at these resolutions is RESCUE, not removal.** Correction turns 67 → 6,790 calls
   at lineage. The rescued calls sit at `cn_mean` median **1.025** (neutral) — they are real specificity,
   not amplicon artifacts. Mechanism: amplicon spikes in individual lines inflate the *permutation null's*
   left tail, costing power; removing them tightens the null and lets genuine, CN-neutral specificity
   through. **Even the headline biology was being masked** — MECOM, SOX17 (Ovary), CDX2 (Bowel) and SPI1
   (Myeloid) all pass corrected but FAIL uncorrected, i.e. they are rescued calls.

**Per-line rank story (line level, rankings-only — the dramatic flips).** Line-level permutation is
degenerate, but the rank shifts are real and null-invariant:

| line | locus | mean CN | uncorrected rank | after correction |
|---|---|---|---|---|
| SK-N-BE(2) | **MYCN** | 215× | #1–#8 | all removed |
| KELLY | **MYCN** | 208× | #1–#14 | all removed |
| NB1643 | **MYCN** | 177× | #1–#14 | all removed |
| COLO320 | **POU5F1B** (8q24, beside MYC) | 120× | #1 | removed |

Correction is **bidirectional** — 36 identity genes improved, 27 stable, 27 demoted. It is not a penalty
on strong signal: P12-ICHIKAWA LEF1 #598→#6 and MCF7 ESR1 #5→#1 improve, while Bowel **CDX2 #1→#14** is
demoted, meaning the uncorrected arm *overstates* the headline colorectal SE. Surviving corrected calls
sit at `cn_mean` median 0.99–1.02 (neutral).

> The earlier rank-based figure of "2,323 amplicon-driven calls" (`cn_ablation.py`, top-15 per group,
> all levels incl. line/subtype) is superseded for reporting by the call-based numbers above. It counted
> top-N rank drops, not calls; the call-based version restricts to the levels the panel supports.

**Case study — OVCAR3.** Uncorrected, **13 of 15** top calls are 19q13 amplicon (chr19:53.6–54.2 Mb, CN
5.6–9.4×). Corrected, 3 of 15.

**MECOM is real, not an amplicon.** Its rank is invariant across both arms at every level; in SKOV3 the
locus is CN-neutral (1.055) while MECOM SEs still take ranks #1–#6; and it is a cluster of six SEs. 3q26
gain across the panel is real but modest (median 1.20 over 11 ovarian lines).

## 5. Cross-layer validation (Phase 6)

Genes near group-specific SEs are themselves specific to that same group, scored by CaCTS on DepMap
expression over the **same 282 lines and same groups**:

| set | per-pair | background | enrichment | per-SE-any | nearest-gene |
|---|---|---|---|---|---|
| permutation, lineage | **18.0%** | 4.36% | **4.1×** | 34.2% | 27.3% |
| permutation, disease | **18.1%** | 3.07% | **5.9×** | 33.2% | 27.4% |
| analytic, lineage | 8.7% | 4.36% | 2.0× | 18.3% | 12.6% |
| analytic, disease | 6.1% | 3.07% | 2.0× | 13.4% | 8.7% |

Controls behave as a *local regulatory* link must, not as a lineage confound: the **group shuffle** (same
gene, random other group) sits at background; **distance decays monotonically**, 30.4% concordant at
<10 kb → 9.3% at 100–250 kb with median rho tracking it (+0.283 → +0.138); and the **direct Spearman** of
SE signal vs neighbour expression — independent of both CaCTS runs — is higher for concordant than
discordant pairs (+0.294 vs +0.160).

That the concordance roughly **doubles** on the permutation-filtered set is independent evidence that the
permutation FDR removes noise rather than signal. The concordance measure knows nothing about it.

Uncurated top nominations: CDX2/bowel (3 kb, rho 0.643), IKZF1/lymphoid, SERPINB5/pancreas, PKP1/squamous,
EN1 and RARRES1/breast.

## 6. EMX2 — an honest partial result

`USE_6049` (chr10:117,543,636–117,545,605) sits on the **EMX2 promoter**, overlapping the EMX2 TSS and the
3′ end of the antisense lncRNA EMX2OS. It is rank #4 for OVCAR3 in **both** CN arms, and it validates
against fully independent DepMap RNA:

| gene | rho | p | |
|---|---|---|---|
| **EMX2** | **+0.461** | 3.0e-16 | hypothesis |
| SOX17 | +0.222 | 1.7e-04 | control |
| PAX8 | +0.046 | 0.45 | control |
| WT1 | +0.016 | 0.79 | control |

Across the nine ovarian lines, H3K27ac and EMX2 expression rank in near-lockstep (1694/6.65, 1068/6.17,
842/4.60 … 1/0.00). EMX2 is a Müllerian-duct developmental TF and the fallopian tube is Müllerian, so the
association is coherent.

**But it does NOT pass the specificity bar.** At HGSOC subtype it is rank #4 with **FDR 0.173**, and
subtype level is unsupported anyway (§2). The correlation is a genuine, independent measurement; the
"specific SE" claim is not currently supported. It was called as a super-enhancer in only **1 of 2,136**
experiments, so "specific H3K27ac element at the EMX2 promoter" is the accurate description.

## 7. What is NOT claimed

- **Any specific-SE count from the analytic null.** It calls 6.05% of tests on shuffled labels vs 7.4% on
  real data — a non-functional FDR (gotcha 71).
- **Subtype- or cell-line-level specificity calls** (§2).
- **SE → target-gene assignment.** Proximity only (gotcha 22). §5 measures concordance in aggregate; it
  does not establish that any individual SE regulates its nearest gene.
- **Anything from the CN-corrected calling-time atlases** — they were never built (gotcha 59).

## 8. Reproducing

```bash
conda activate atac_hdac
python phase2/score_pilot.py \
    --signal  phase2/results/atlas.s3.se_signal.tsv.gz \
    --catalog phase2/results/atlas.s3.union_catalog.bed.gz \
    --norm none --fdr-method permutation --n-perm 1000 --dump-specific 0.10 \
    --out phase2/scores/atlas.s3.perm                       # ~18 min, ~4.9 GB

python phase2/analysis/cn_ablation.py       --corrected phase2/scores/atlas.s3 \
                                            --uncorrected phase2/scores/atlas.s3.nocn
python phase2/analysis/concordance_bridge2.py --scores phase2/scores/atlas.s3.perm \
                                            --levels OncotreeLineage,OncotreePrimaryDisease
python phase2/analysis/null_calibration.py  --jsd-pickle <saved JSD>
python phase2/analysis/se_drilldown.py USE_6049 --genes EMX2,PAX8,SOX17,WT1
```
