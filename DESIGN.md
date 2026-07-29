# Design & decisions

Everything decided so far, from the scoping discussion (2026-07-16). Where a choice is still open it is marked
**[open]**. See `PRIOR_ART.md` for what exists and `DATA_SOURCES.md` for the data.

## What CaCTS needs, and what each ingredient becomes for a super-enhancer

CaCTS requires four things: (1) a **compendium** of comparable measurements across many lineages; (2) a
**per-feature quantity** comparable across all of them (for a TF: its expression — a gene has one universal
identity); (3) a **specificity statistic** (JSD to a one-hot ideal); (4) an **empirical null + FDR** and a
**target group**. Ingredient (2) is where SEs are harder than TFs, for three reasons:

1. **No stable cross-sample identity** — an SE is a *called region* with sample-specific boundaries. You must
   first build a common coordinate system (a union/atlas of regions) before "the same SE" means anything.
2. **Signal is a property of position, not an object** — H3K27ac density is continuous; the SE is a
   discretization. Specificity is more honest computed on the *signal over a region* than on presence/absence.
3. **Copy-number confounding** (cancer-specific) — more DNA copies → more H3K27ac → a locus looks like a
   stronger/more-specific SE. Uncorrected, the score crowns amplicons (MYC). This is the crux and the novelty.

## The approach landscape (and where we landed)

| # | approach | what it is | verdict |
|---|---|---|---|
| 1 | Region-anchored H3K27ac specificity | union SE catalog; quantitative H3K27ac per region per lineage; a specificity index | **substrate** — this is the data layer |
| 2 | Binary SE-presence specificity | present/absent per lineage; categorical | discarded (throws away magnitude; what SEdb/dbSUPER do) |
| 3 | Target-gene proxy | SE inherits its target gene's CaCTS score | keep as a **cheap cross-check**, not the method |
| 4 | Regulator-content / motif grammar | SE specificity = aggregate specificity of the TFs it carries | keep as an **added feature / validation** (ties to the CRC) |
| 5 | Differential / one-vs-background framing | compare each lineage's SE signal against the atlas; effect + FDR | **the chosen operation** (see below) |
| 6 | (orthogonal) the statistic | tau / JSD / entropy / Gini | secondary — see "Statistic" |

**Decision — the core operation is #5, in its *individualized* form.** "Build the atlas, then compare each
lineage against it" is a one-vs-rest (individualized) computation that yields a score per **(SE, lineage)** —
"how specifically is this SE active in HGSOC vs everything else." This is exactly the CaCTS shape (CaCTS scores
gene *i* against a one-hot ideal centered on cancer type *j*) and is what we actually want ("give me the
HGSOC-specific SEs"). Note the prior SE tools mostly used *general* (dispersion) metrics — tau (Ryu), entropy
(SEA) — which answer "how cell-restricted is this SE" but not *which* lineage. Choosing the individualized form
is both the correct instinct and the less-trodden path.

## The unit

- **Region-anchored union SE catalog** (build a union/consensus of SE calls across the compendium; re-quantify
  H3K27ac signal over it per lineage). SEdb/dbSUPER store *calls*, not signal tracks, so signal must be
  re-quantified from ChIP-Atlas bigWigs/BAMs. **[open]** SE-domain resolution vs **constituent-enhancer**
  resolution (cSEAdb used constituents — finer, arguably better localizes the specific element).
- Grouping resolution: score at both the **cell-line** level and the aggregated **lineage/subtype** level.
  Average replicates to one vector-entry per group.

## The specificity statistic

> **Superseded in part (2026-07-22).** Option 1 was adopted and is what shipped — but its *null* was not.
> The analytic (normal-approximation) empirical-null FDR failed an outright calibration test (6.05% false
> calls on shuffled labels) and was replaced by a **label-permutation null**, `phase2/permutation.py`.
> The statistic (JSD-to-one-hot) is unchanged. See `RESULTS.md` and gotcha 71.

The particular statistic is **secondary** (owner is agnostic; requirement = defensible as the right choice).
Recommendation, in order:

1. **Lead: CaCTS JSD-to-one-hot + empirical-null FDR**, applied to the CN-corrected, normalized H3K27ac vector.
   *Defense:* methodological consistency — "we applied the identical specificity framework to the super-enhancer
   layer that we used for transcription factors." One coherent methodology across the TF and SE layers of the
   same body of work beats any single metric's benchmark ranking. Reuses the pyCaCTS null (permute lineage
   labels for FDR).
2. **Fallback / sensitivity: regression with a copy-number offset.** Model per-SE H3K27ac counts ~ lineage with
   a CN offset → per-lineage effect + FDR. *Defense:* CN correction lives naturally as an offset rather than a
   pre-correction, and this mirrors the standard CN-covariate treatment in contact/coverage models.
3. If a reviewer wants the field-validated index: **tau** is the benchmarked winner (Kryuchkova-Mostacci 2017),
   but it is *general* — report it as a descriptive companion ("how cell-restricted is this SE overall"), not as
   the lineage-specific score.

**The statistic is not the hard part — normalization is (see below).**

## Copy-number correction (the novel core)

Apply the Su et al. 2025 recipe: per cell line, obtain copy-number (DepMap/GISTIC segments where available;
CNVkit from ChIP/ATAC input DNA otherwise), then **divide H3K27ac signal by the local copy-number ratio** in
gained regions (scale up losses) *before* scoring — or carry CN as an **offset/covariate** in the regression
form. Goal: a locus must be H3K27ac-high *beyond what its copy number explains* to count as specific. This is
the inverse of Zhang 2016 (which used the CN↔H3K27ac coupling to find amplified SEs) and directly defuses the
"you just found amplicons" critique. **[open]** how to handle lines with no DepMap CN — CNVkit-from-input,
impute, or exclude.

### Should we even correct? — the trade-off, and the decision (2026-07-16)

**The key realization: copy number is *both* a confounder *and*, sometimes, the mechanism.** Two things are
tangled inside the CN↔H3K27ac coupling, and the divide-by-CN recipe removes **both** because it can't tell them
apart:

- **An assay artifact** — ChIP/ATAC read density scales with template abundance (20 copies → ~20× the DNA to
  IP and sequence). Removing *this* is uncontroversial measurement hygiene, like input/GC/mappability correction.
- **A biological dosage effect** — focal SE amplification is a bona fide lineage-defining oncogenic mechanism
  (Zhang 2016: KLF5, MYC-LASE/ECSE). Here the amplification *is* why the SE is specifically active; correcting it
  away deletes a true positive.

Framed by the question each readout answers: **raw H3K27ac = aggregate SE output** (total activity by any means,
incl. dosage); **per-copy H3K27ac = the regulatory *state* of the element** (is it *wired* as an identity SE).
CaCTS's spirit is identity/master-regulator wiring, so the per-copy state is the more apt readout — a principled
argument *for* correction, not merely a defensive one.

**Pros (correct):** (1) isolates the CN-independent regulatory axis that is the project's novelty; (2) defuses
the dominant "you just re-found amplicons" critique; (3) avoids double-counting the genomics/GISTIC layer, whose
value-add SE-CaCTS is supposed to exceed; (4) restores cross-lineage comparability of the JSD vector (ploidy
differences between lineages otherwise distort it).

**Cons (don't over-correct):** (1) **over-correction erases real biology** — a lineage whose identity is achieved
*through* SE amplification loses a true, functional, lineage-specific hit; (2) the **signal ∝ CN linearity
assumption is only approximate** — breaks under chromatin/antibody saturation, sub-clonal/heterogeneous (non-
integer effective) amplification, allele-specific acetylation, and ecDNA (hyper-active per copy); (3)
**attribution ambiguity** — a ROSE-stitched SE can straddle a CN-segment boundary, so which single CN ratio
applies? (constituent-level is cleaner but noisier; gene-level DepMap CN is too coarse — want segment-level mapped
to the SE interval); (4) **coverage + variance cost** — no-CN lines must be excluded (shrinks the ~295-type
breadth headline) or imputed (adds error), and dividing by a noisy CNR inflates variance most in the gained
regions we care about.

**Decision — compute and report BOTH, as co-primary outputs (not one as a buried "sensitivity analysis").** The
correction is run as an *explicit, ablatable, reported* dimension, never a silent preprocessing step:

1. **Both scores are first-class deliverables:** the CN-corrected specificity score *and* the uncorrected score,
   shipped together. The **delta between them is itself a headline result** ("N 'specific' SEs are amplicon-driven;
   M survive dosage correction") — and that ablation only exists if we *don't* silently correct.
2. **Prefer the CN-offset/covariate form over hard division** — it propagates uncertainty, degrades gracefully
   with noisy/missing CN, and follows the usual CN-covariate precedent. Hard division is the simpler
   cross-check, not the primary.
3. **Label two complementary SE categories rather than declaring one "right":** *dosage-independent* lineage-
   specific SEs (survive correction → the epigenetic-identity story) and *amplification-associated* lineage-
   specific SEs (specific only because amplified → the Zhang-2016 driver story). The correction's job is to
   **label which is which**, not to invalidate either.
4. **MYC/CCAT1 is the worked example** precisely because it is ambiguous (amplified *and* lineage-associated);
   showing what happens to its rank under correction is the demonstration that the two-category split is real.

**[open]** no-CN lines (CNVkit-from-input / impute / exclude); and CN attribution for a stitched SE spanning
heterogeneous CN segments (max vs mean vs constituent-weighted, or correct at constituent resolution).

## Additional data as model features (see `DATA_SOURCES.md`)

- **DepMap (highest value, joins by cell line):** CN (correction input), expression (target-gene readout — lets
  a specific SE be cross-checked against the *gene*-CaCTS score), CRISPR dependency (does the specific SE's
  target actually matter → function). "Lineage-specific + amplification-independent + drives a dependency gene"
  is a hard-to-dismiss SE nomination.
- **Other ChIP-Atlas marks** to sharpen the active-SE state: H3K4me1 (enhancer), H3K4me3 (separate SE from active
  promoter), H3K27me3 (a truly specific SE is H3K27ac-high in the target and Polycomb-repressed elsewhere).
- **SE-defining coactivators** (BRD4, MED1, EP300, CDK7) — corroborate SE identity, guard against H3K27ac
  artifacts.
- **Lineage master-TF ChIP** (PAX8, SOX17, …) — a lineage-specific SE should be bound by the lineage's masters;
  validation *and* feature, and the bridge to the CRC work.
- **ATAC/DNase** — orthogonal, less antibody-batch-prone activity readout.

**Breadth-vs-depth caveat:** every mark you *require* shrinks the cell-type intersection. Keep **H3K27ac as the
primary axis** (max breadth, ~295 types); treat the other marks/factors as *corroborating features where
available* (down-weight/impute where missing), never hard requirements — otherwise the breadth advantage is lost.

## The main technical risk

**Cross-study H3K27ac normalization** across ~1,789 experiments (many labs, antibodies, depths) is where this
project succeeds or fails — far more than the choice of statistic. ChIP-Atlas helps (one uniform MACS2/bigWig
pipeline), but batch/antibody effects remain. Plan: quantile-normalize the region×lineage signal matrix and/or a
reference-based normalization; validate with a small pilot before committing (see `ROADMAP.md`).

## The unifying story (the biggest novelty)

SE-CaCTS lets *one* specificity framework span three layers of the same study — **TF expression** (CaCTS, done),
**SE H3K27ac** (this project), and **target-gene expression** (gene-CaCTS) — triangulating lineage-defining
regulatory elements and their masters. That integration, more than CN or breadth alone, is the thing none of the
prior tools did.

## Novelty framing for an eventual paper

- **Headline:** breadth (largest cancer H3K27ac compendium with a per-lineage SE specificity score) + integration
  with the master-TF/CRC layer.
- **Trust/correctness:** copy-number correction (the unique, defensible-against-"it's-just-amplicons" feature).
- **Method:** the individualized CaCTS-JSD form on SEs (consistency with the TF layer).

## Methods & tooling decisions (tooling survey, 2026-07-17)

A cited best-practice survey of every compute-intensive step (deep-research run; 20 claims verified 3–0; the
CN-correction specifics were verified separately, see `STEP_A_NOVELTY_CHECK.md` sibling notes / status log).
Two decisions **changed** from the pilot's defaults:

| step | decision | why / citation |
|---|---|---|
| **normalization** ⚠️ **changed** → ✅ **BUILT 2026-07-19** | **S3norm** as the Phase-2 normalizer (`phase2/s3norm.py`, `aggregate.py --norm s3norm`); **quantile only as a baseline**; S3V2's background-variance step still open | Quantile *manufactures false positives* by forcing identical distributions — S3norm benchmark: quantile-specific CTCF peaks had the motif **9.2%** vs **80.5%** for consensus (Xiang et al., *NAR* 2020, e43). S3norm co-normalizes depth **and** signal-to-noise; **validated on H3K27ac**. S3V2-IDEAS (2021, PMC8479670) is built to remove **cross-lab** variation — our exact problem. (Pilot finding "batch is normalizable" stands; only the normalizer is upgraded.) **Build notes:** fitted at **fine-grid** resolution and re-summed, since `Σf(xᵢ) ≠ f(Σxᵢ)`; anchors are **signal quantiles, not peak calls** (per-sample peak counts span 53–23,528 → call-based anchors collapse, and the call threshold is itself SNR-dependent, hence circular). Validated: fitted exponent **B tracks dynamic range, Spearman −0.956**; HT29/HCT116 scale **−1.315 → −0.043** log2. |
| **QC filtering** ⚠️ **now COUPLED to normalization** | A peak-count QC gate is **mandatory** when s3norm is on (`aggregate.py --min-peaks`, pilot ≥2000) | Discovered 2026-07-19, not anticipated by the survey: low-enrichment samples have compressed dynamic range, so s3norm fits an **expanding** exponent (B>1) that amplifies their noise into apparent specificity. **Ungated, s3norm was worse than quantile** — Breast ESR1 fell #5 → #21 and identity-gene recovery in the top 15 went 3 → 0. Gated, biology matches quantile (ESR1 **#2**, SOX17 **#8**) with the scale correctly equalized. So the "use Cistrome QC to filter ChIP-Atlas" row below is not an independent nicety — **adopting S3norm without a QC gate is a regression.** |
| **reuse vs reprocess** ⚠️ **new option** | **Use Cistrome DB QC to filter our ChIP-Atlas set**; keep ChIP-Atlas for breadth | Cistrome DB v3.0 (~45k human samples, uniform CHIPS pipeline, **6-metric QC** + regulatory-potential scores; *NAR* 2024, D61) vs ChIP-Atlas (76,217 experiments, ~90% of SRA ChIP-seq). Breadth from ChIP-Atlas, QC rigor from Cistrome. |
| **SE calling** | **ROSE** for the catalog; **evaluate LILY** (may fold CN correction into calling) | ROSE is the de-facto standard (SEdb 3.0, cSEAdb both use it; `rose2`/`se_rose` local). **LILY** (BoevaLab) detects SEs from H3K27ac **with built-in CN correction** via HMCan. |
| **union SE catalog** | **cSEAdb recipe**: merge per-sample ROSE calls where regions overlap **≥25% of width**; work at **constituent-enhancer** resolution | Explicit, reproducible; cSEAdb (*PLOS Comput Biol* 2024). Resolves the `[open]` SE-domain-vs-constituent question toward constituents. |
| **signal quantification** | keep **deepTools `multiBigwigSummary`** over the region set (from bigWig) | Confirmed fit-for-purpose. Rigorous alternative = recount from BAM (DiffBind/csaw). Key caveat: *background-vs-enriched normalization choice matters more than the counting tool.* |
| **CN correction** | **offset/windowed** correction preferred over naive divide-by-ratio; use **HMCan/LILY** to estimate CN **from the input track** for lines lacking DepMap CN | Verified (2026-07-17): Su 2025 *is* a divide-by-CNR scaling (÷CNR in gains, ×CNR in losses) — confirmed, so our offset choice is a deliberate upgrade, not a misread; HMCan estimates CN from the control/input via **Control-FREEC** and corrects each profile position — confirmed (answers the no-CN-line `[open]`); CN "could drive, if not dominate, differential signals" — confirmed (the ~70% figure is G4-ChIP-specific, ~20% ATAC confirmed). |
| **CN correction (alt normalizer)** | **CHIPIN** is a fallback where matched RNA-seq exists | Verified: CHIPIN anchors normalization on *constant-expression* genes' regulatory regions; on ACC H3K27ac it beat LILY (34.75%→1.8% vs 18%). Note: CHIPIN's win uses **quantile normalization constrained to constant genes** — so it's *unconstrained* quantile that's the problem, not quantile per se. |
| **cell-line harmonization** | **Cellosaurus** CVCL as the join key; the **`cellosaurus` R pkg `mapCells()`** maps free-text / DepMap / ATCC IDs → CVCL | Standard (RRIDs, synonyms). Caveat: RRID↔model not always 1:1 (e.g. SJRH30 & RH30 → CVCL_0041) — join with care. Local `cellosaurus/` is an empty stub → must fetch. |
