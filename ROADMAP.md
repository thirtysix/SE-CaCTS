# Roadmap

A side project — no deadline. Sequenced so the two **cheap de-risking steps come first**, before any heavy data
pull. Read `DESIGN.md` and `PRIOR_ART.md` before starting.

## Do these two first (cheap, high-information)

- [x] **A. Lock down the novelty claim.** ✅ **Done 2026-07-16** (deep-research follow-up run; see
      `STEP_A_NOVELTY_CHECK.md`). SEdb 3.0 (NAR 2026) fully checked — still ROSE-on-input-normalized H3K27ac, no
      specificity statistic, no CN correction. No 2024–26 work computes a copy-number-corrected per-SE
      lineage-specificity score. **"To our knowledge, the first copy-number-aware per-SE lineage-specificity
      score" is defensible.** Residual (non-blocking, re-check before manuscript "first" language): the *live*
      SEdb 3.0 site vs its paper, and the ecDNA/pediatric amplicon-SE subfield.
- [x] **B. Normalization pilot.** ✅ **Done 2026-07-16 → GREEN** (`pilot/PILOT.md`, `pilot/results/`). 13 H3K27ac
      experiments across 4 lineages + a 3-study MCF7 batch probe, quantified over genome bins **and** H3K27ac
      peaks under 3 normalizations. **Quantile normalization** makes the cross-study MCF7 replicates cluster
      (cohesion 0.26–0.29 ≪ 1) and gives the only positive lineage silhouette; **median-of-ratios fails** (one
      outlier dominates). ⇒ cross-study batch is normalizable — the dominant risk is defused. **Decisions
      carried forward:** adopt **quantile** as the Phase-2 normalization; score at **subtype resolution** too
      (breast splits luminal vs basal — real biology, not batch). Full local run in minutes; Roihu scaffold
      (`pilot/roihu/`, `$SECACTS_CSC_PROJECT`) staged for the Phase-1 scale pull.

## Phase 1 — data assembly

- [x] **Manifest + harmonization + DepMap join — done 2026-07-17** (`phase1/`). Fetched `experimentList.tab`
      (authoritative, +QC) and Cellosaurus; `phase1/scripts/11_harmonize_and_join.py` → `phase1_manifest.tsv`
      (11,827 hg38 H3K27ac). Cellosaurus CVCL bridge → DepMap `RRID`/Oncotree + CN/dependency flags + QC filter.
      **Result: 2,917 QC-pass, DepMap-joined, CN-correctable experiments on ~332 lines** (`DATA_SOURCES.md`).
      Synonyms verified (OVCAR-3 → CVCL_0465 → ACH-000001). *(Optional: cross-filter with Cistrome QC.)*
- [x] **Per-line identity, CN provenance & lineage recovery — done 2026-07-22** (`phase1/`, scripts 12–14).
      Front door `CELL_LINE_IDENTITY.md`; CN deep-dive `CN_COVERAGE.md`. Scripted the §4 CN-source audit
      (`12`, DR-routing + input audit), recovered lineage for the non-DepMap lines from Cellosaurus `DI` via an
      empirical NCIt→Oncotree crosswalk (`14`), and flattened to a browsable `data/cell_line_status.tsv` (`13`).
      **Fixed two identity-resolution bugs in the join** (`11`): a `cell_desc="NA"` sentinel sink (4,206
      experiments) and non-human/wrong-subline homonym CVCLs (HAP1, KG-1, PC-3, …) — atlas SRX unaffected.
      **Investigable scope = 631 human cancer lines, all with a lineage** (722 resolved − 91 non-cancer/non-human).
- [x] **Pull the QC-pass bigWigs at scale — DONE 2026-07-21** (`phase2/array.slurm`). Full 2,917-sample pull
      on Roihu: 2,916 done, 1 retired SRX, **zero OOM**, ~43 BU. Streamed per-SRX (download → cnrose → retention
      set → delete). The atlas is built (see Phase 2). Staged rollout smoke→small→medium→large→full caught six
      bugs before the full run; live findings in `PULL_DESIGN.md §8.6-8.7`.

## Phase 2 — reference SE atlas + signal matrix

- [x] **Quantification grid built — done 2026-07-18** (`phase2/`). Chose per-SRX ChIP-Atlas peak BEDs fetched
      directly on Roihu (`21_grid_from_persrx.sh`) over the 107 GB AllCell download. Compared thresholds →
      **`bed20` (Q<1E-20): 505k regions, 11% genome, full 2,917/2,917 coverage** (`phase2/data/grid.20.bed`).
      This fixed grid is the single-pass quantification reference (SEs recovered as sums of grid rows).

- [x] **SE caller BUILT + validated — `cnrose` (2026-07-19).** Instead of shelling out to `ROSE2_callSuper` +
      deepTools, we built **`cnrose/`** — a bigWig-native, CN-aware, ROSE-style caller (pure pyBigWig+numpy+scipy)
      **validated bit-for-bit vs ROSE2** on the 13 pilot bigWigs (stitch 13/13 exact, SE Jaccard **1.0000**). It
      is the idealized "CNV-aware ROSE" LILY only approximates — ROSE algorithm, bigWig-native, external CN,
      optional + separable. **LILY rejected** as-is (BAM/HMCan, non-optional CN, needs matched Input). See
      `cnrose/DESIGN.md`. Per-sample: stitch `bed20` (12.5 kb) → quantify → tangent cutoff.
- [x] **Union catalog + signal matrix BUILT — `phase2/aggregate.py` (2026-07-19).** cSEAdb ≥25% reciprocal-overlap
      merge → union SE catalog; grid→SE Σ(mean·len) reconstruction (no 2nd bigWig pass), exact on the pilot.
      **Dual catalogs** (agnostic + CN-corrected). Region × sample → cell-line → lineage via `pyCaCTS.build_rep_matrix`.
- [x] **Cross-study normalization BUILT — `phase2/s3norm.py` (2026-07-19).** S3norm's two-parameter monotone
      power transform `f(x)=A·(x+p)^B`, matching a reference in BOTH the enriched and background strata (the
      ratio equation eliminates A → one monotone root-find in B). Fitted on the **fine grid** and re-summed
      inside `aggregate.py --norm s3norm`, because the transform is nonlinear (Σf(xᵢ) ≠ f(Σxᵢ)); the SE catalog
      is untouched (SE calling is a within-sample threshold). Quantile retained as `--norm quantile` baseline.
      **Anchors are signal-quantile, not MACS2 calls** — the pilot's per-sample peak counts span 53–23,528, so
      call-based common-peak sets collapse for weak samples, and the call threshold is itself SNR-dependent
      (circular). Validated: reference self-normalizes to identity; fitted **B tracks per-sample dynamic range
      (Spearman −0.956, p=3.4e-07)**; cross-study scale equalized (HT29/HCT116 co-called median log2 ratio
      **−1.315 → −0.043**; column means 2.6× apart → 5%).
- [x] **QC gate — `aggregate.py --min-peaks` (2026-07-19).** NOT optional with s3norm: the two DESIGN decisions
      (S3norm; Cistrome-style QC filtering) are **coupled**. Low-enrichment samples have compressed dynamic
      range, so s3norm fits an *expanding* exponent (B>1: T47D 1.61, MDA-MB-231 1.36) that amplifies their noise
      into apparent specificity. Ungated, s3norm was **worse than quantile** (Breast ESR1 #5 → #21, id@15 3 → 0);
      gated at ≥2000 peaks (a natural gap in the pilot: 2046 vs 361) it matches quantile on biology
      (ESR1 **#2**, SOX17 **#8**, id@15 = 2) with the scale properly equalized.
- [ ] **S3V2 background-variance step** — the remaining piece of the S3V2-IDEAS upgrade (S3norm normalizes the
      two strata *means*; S3V2 additionally normalizes background *variance*). Not yet needed on the pilot.

## Phase 3 — copy-number correction (the novel core)

- [ ] Map copy number onto each region per cell line: DepMap `OmicsCNGeneWGS` where available; for lines lacking
      DepMap CN, estimate CN **from the ChIP input track** via **HMCan/Control-FREEC** (LILY's approach; verified
      2026-07-17) — resolves the no-CN-line `[open]`.
- [ ] Prefer a **CN offset/covariate** over naive divide-by-ratio (Su 2025 is a ÷CNR scaling; the offset
      propagates uncertainty and follows the usual CN-covariate precedent). Keep divide-by-ratio as a cross-check.
- [ ] Sanity check on a known amplicon: does correction drop the MYC/CCAT1 SE's apparent specificity where it is
      amplified? (The whole point.)

## Phase 4 — specificity scoring

- [x] **Proof-of-mechanism on the pilot — `phase2/score_pilot.py` (2026-07-19).** Reuses `../pyCaCTS`
      (`cacts_score_matrix` JSD, `empirical_fdr`, `build_rep_matrix`). Scores at **every Oncotree resolution**
      (Lineage → PrimaryDisease → Subtype → cell line), CN-corrected + uncorrected. **Recovers** SOX17 (Ovary),
      ESR1 (Breast, #1 in MCF7); the **CN ablation** (ΔJSD after correction) labels amplification-associated SEs
      (BCAT1 12p, 19q13, BCAS1 20q13). Specificity cutoff = **empirical-null FDR ≤ 0.10** (per-group counts).
      *Caveat:* absolute counts are panel-sensitive (rankings robust). **Revised 2026-07-19:** this is NOT a
      normalization artifact — S3norm fixes the scale and the counts persist (OVCAR3 24 vs Kuramochi 183 among
      QC-pass samples); sample QC and private-SE-call counts were also tested and falsified (p>0.6).
- [x] **Comparable specificity counts — `phase2/specificity.py` (2026-07-20).** Tested the full 2×2 of null
      calibration × BH scope. **The global BH is the fix; the global null is not.** Sharing ONE
      multiple-testing budget across all SE × group tests (rather than an independent budget per group)
      compresses the count spread — max/min **248 → 3.9** at Disease/Subtype, **10.3 → 4.5** at line level —
      removes pathological zero-call groups (Breast 0 → 12), and *strengthens* both identity genes
      (SOX17 −1.73 → −2.36, ESR1 −1.51 → −2.62). A pooled/global **null**, the intuitive fix, is over-conservative
      for tight groups and drops SOX17 below significance (−0.18). Now the default
      (`score_pilot.py --fdr-scope global`); `--fdr-scope pergroup` reproduces pyCaCTS. Rankings are invariant
      throughout. HGSOC gap 11.2× (quantile, per-group BH) → 7.6× (S3norm) → **3.4×** (+ global BH).
- [x] **RUN AT SCALE — 2026-07-22.** `score_pilot.py` on the real S3norm atlas (**42,943 SEs × 2,136 samples
      → 282 lines**, 24 lineages / 44 diseases / 75 subtypes; `--norm none`; 3:11 wall, 4.8 GB).
      Outputs `phase2/scores/atlas.s3.*`. Recovers known biology unsupervised: Ovary **MECOM #1, SOX17 #6**;
      Bowel **HNF4A #4, CDX2 #14** (pilot had CDX2 #404); Myeloid **SPI1 #9**; at line level MCF7→**ESR1 #1**,
      THP-1→CEBPA #3, MOLM-13→IRF8 #2, SW48→CDX2 #5, P12-ICHIKAWA→LEF1 #6. **Negative control:** 6 TNBC lines
      bury ESR1 at rank 11k–21k and the ER-negative subtype group puts it at #19,059, while lobular
      (near-always ER+) gives ESR1 #22. **Counts are now comparable** — 0 zero-call groups of 425, spread
      3.4–5.6×, i.e. the global-BH fix generalizes from the pilot's 3.9–4.5×.
- [x] **A SILENT TOTAL FDR FAILURE, found and fixed (2026-07-22)** — see gotcha 64. The first at-scale run
      returned FDR **exactly 1.0 for every test at 3 of 4 levels** (401/425 groups, zero calls even at
      FDR≤0.25) and read as "no significant results". Four-link chain: CN correction emits negatives →
      pycacts propagates NaN → `null_params` let one NaN void a whole group → `_bh_log10` swallowed it
      because `min(0.0, nan)` is `0.0`. Fixed at every layer, with 4 new regression guards in
      `tests/test_specificity.py` (12/12) and pyCaCTS equivalence still `max|diff|=0`.
- [x] **LABEL-PERMUTATION FDR — DONE 2026-07-22** (`phase2/permutation.py`, `--fdr-method permutation`).
      The analytic null FAILS an outright calibration test: on SHUFFLED labels it makes **62,321 false
      calls of 1,030,632 (6.05%)** vs the permutation's **0** — i.e. it called almost as many SEs
      "specific" on pure noise as on real data (7.4%). Permuting which line carries which group label
      preserves group SIZE exactly (JSD depends strongly on it), stays per-group (keeps gotcha 27) and
      reuses `_bh_log10` (keeps gotcha 64). Median calls/group (over all groups at the level):
      **81.5 lineage / 15 disease / 0 subtype** — 94 / 18 over groups with ≥1 call — vs thousands
      analytic. All anchors still pass at lineage+disease (Ovary MECOM #1 +
      SOX17 #6, Bowel HNF4A #4 + CDX2 #14, Myeloid SPI1 #9). Independently confirmed by Phase 6: the
      concordance bridge roughly DOUBLES on the permutation set (2.0x -> 4.1x/5.9x enrichment).
      **Use B=1000** — at B=50 the p-floor 4.7e-7 sits above the subtype BH bar 3.1e-8 and returns a
      spurious 0; `permutation_fdr` now warns RESOLUTION-LIMITED.
- [x] **PANEL RESOLUTION LIMIT FOUND — subtype level is NOT supported** (1 call across 75 groups, robust
      to B=1000). 29 of 75 subtypes have ONE line, 56 have <=4. Report lineage + primary disease; treat
      subtype and line as RANKINGS ONLY. Fixing this needs more LINES PER SUBTYPE (scope expansion), not
      more samples per line. Consolidated claims now live in **`RESULTS.md`**.
- [ ] Still TODO at scale: the co-primary corrected/uncorrected atlas; consider
      upstreaming `specificity.py` into pyCaCTS once the same effect is confirmed on the TF layer.
- [ ] At scale: **CaCTS JSD-to-one-hot** on **both** CN-corrected **and** uncorrected normalized vectors —
      **co-primary outputs** (decision 2026-07-16, `DESIGN.md`), not one as a buried sensitivity analysis.
- [ ] CN handling: **regression-with-CN-offset is the primary form** (propagates uncertainty, graceful with
      missing/noisy CN); hard divide-by-CN-ratio is the simpler cross-check. Report **tau** (general) as a
      descriptive companion.

## Phase 5 — validation

- [ ] Positive controls: does a lineage recover its known identity SEs (e.g. HGSOC → PAX8/SOX17 SEs)? Does CN
      correction change the MYC SE's rank as expected?
- [ ] Benchmark against cSEAdb (NCI-60), SEA v4.0 (entropy), Ryu (tau): overlap and where CN correction diverges.
- [x] **Ablation DONE AT SCALE — 2026-07-22** (`score_pilot.py --no-cn` + `phase2/analysis/cn_ablation.py`).
      Both arms evaluate CN and report `cn_mean` per call; `--no-cn` only skips APPLYING it, so they differ in
      exactly one step. **It rediscovers textbook amplicons unprompted:** SK-N-BE(2) / KELLY / NB1643 →
      **MYCN at cn 177–215×**, ranked #1–#14 uncorrected and all removed by correction (the canonical
      MYCN-amplified neuroblastoma lines); COLO320 → **POU5F1B (8q24, beside MYC) at 120×** (the classic
      MYC-ecDNA line). 2,323 amplicon-driven calls → `atlas.s3.cn_ablation.tsv`. Correction is genuinely
      **bidirectional** (36 identity genes improved / 27 stable / 27 demoted), so it is not a penalty on high
      signal: P12-ICHIKAWA LEF1 #598→#6 and MCF7 ESR1 #5→#1, but Bowel **CDX2 #1→#14** — the uncorrected arm
      *overstates* CDX2. Surviving corrected calls sit at `cn_mean` median 0.99–1.02 (neutral).
      **OVCAR3 is the case study:** uncorrected, **13/15** top calls are 19q13 amplicon (chr19:53.6–54.2 Mb,
      cn 5.6–9.4×); corrected, only 3/15, and EMX2 (#4 in BOTH arms) is one of just two survivors.

## Phase 6 — integration (the differentiator)

- [x] **Specificity-concordance bridge, STAGE 1 — 2026-07-22** (`phase2/analysis/concordance_bridge.py`).
      Runs CaCTS on DepMap **expression** over the SAME 282 lines and SAME Oncotree groups as the SE atlas,
      then asks whether genes near group-specific SEs are themselves specific to that same group. This is
      what replaces the proximity-only proxy (gotcha 22) with a controlled association.

        level                    concordant   background   shuffled   enrichment
        OncotreeLineage             26.3%        4.4%        3.9%      6.0x  p=1.6e-45
        OncotreePrimaryDisease      24.6%        3.1%        4.6%      8.0x  p=6.1e-94
        OncotreeSubtype             21.5%        2.6%        3.6%      8.1x  p=7.9e-139

      Controls all behave as a *local regulatory* link must, not as a lineage confound: the **group shuffle**
      (same gene, random other group) collapses to background; **distance decays monotonically** over all six
      bins, 47.2% at <10 kb → 11.3% at >250 kb with median rho tracking it (+0.260 → +0.148); and the
      **direct rho** (SE H3K27ac vs neighbour expression, independent of both CaCTS runs) is >0 for 92.6% of
      pairs and 99.4% of concordant ones. Top nominations are uncurated known biology — CDX2/bowel (3 kb,
      rho 0.643), IKZF1/lymphoid, SERPINB5/pancreas, PKP1/squamous, EN1+RARRES1/breast, EMX2/HGSOC.
      **Scope caveat:** stage 1 measures the top-15 SEs per group (1,208 SEs / 2,130 pairs) against the single
      NEAREST gene — both optimistic samples.
- [ ] **Stage 2:** extend to ALL FDR≤0.10 SEs and to every gene in a window (needs a `--dump-specific` option
      on the scorer). Expect the headline rate to fall; that number is the honest genome-wide one.
- [ ] Then: SE–TF binding / the CRC layer, and DepMap **dependency** → the triangulated
      one-framework-across-three-layers story.

## Phase 7 — write-up

- [x] **Browsable results dashboard — DONE 2026-07-24** (`docs/`). Static, no-backend, modeled
      on `../pyCaCTS/dashboard/`: Overview / SE atlas / CN ablation / Concordance / About. Reads the
      PERMUTATION scores; enforces the resolution rule in the UI (calls at lineage+disease; subtype/line
      rankings-only with FDR shown but never called). Staged by `phase2/scripts/60_stage_dashboard.py`.
      Verified headlessly — all tabs, both themes, the guard, and the biology (Lymphoid→IGLL5/ETS1/IKZF1/
      BCL6, AML→CEBPA/IRF8/GFI1) render correctly.
- [x] **DEPLOYED + PUBLIC — 2026-07-28.** Repo made public; dashboard live at
      <https://thirtysix.github.io/SE-CaCTS/> via Pages (`main` + `/docs`). `phase2/dashboard/` was moved to
      **`docs/`** so Pages serves it at the root URL with no redirect shim. A sixth tab (**SE finder**) and
      the atlas→concordance cross-links landed first; the repo is linked from the Overview and the sidebar.
      The short-lived **analytic-FDR toggle was removed entirely** — UI, the `fdr_analytic` column, and its
      staging — so no void number is reachable from the site (gotchas 74–75).
- [ ] Method paper / tool + the deployed public resource. Framing: breadth + integration
      (headline), CN correction (trust), individualized JSD (method). Cite Ryu / cSEAdb / SEA as the lineage;
      claim novel-in-combination + CN-aware.

## Open questions (carried from the discussion)

- ~~SE-domain vs constituent-enhancer resolution?~~ **Decided (tooling survey): constituent** (cSEAdb recipe).
- ~~Best CN reference; how to handle no-CN lines?~~ **Resolved:** DepMap `OmicsCNGeneWGS` where present; else
  **HMCan/Control-FREEC from the ChIP input track** (LILY's approach). Prefer offset over divide-by-ratio.
- ~~How bad is cross-study H3K27ac batch, really?~~ **Answered (step B): normalizable** — quantile makes
  same-line/different-study replicates re-cluster (cohesion ~0.28). Residual: re-confirm at scale (more
  lineages/replicates) since the pilot was 13 experiments.
- JSD vs tau vs regression — which reads as most defensible once we see the data? (statistic is secondary.)
- ~~Does SEdb 3.0 / any 2025–26 work reduce the novelty?~~ **Answered (step A): no.** Two residuals only —
  live SEdb 3.0 site vs paper, and the ecDNA/pediatric amplicon-SE subfield (both non-blocking).
- **[new, from step A — a design question, not a novelty threat]** How to assign a single copy-number ratio to a
  ROSE-stitched SE that spans heterogeneous CN segments (multi-constituent SE over a CN boundary)? Bears on Phase 3.

## Status log

- **2026-07-16** — Project inception. Deep literature search completed (`PRIOR_ART.md`); design and data plan
  drafted (`DESIGN.md`, `DATA_SOURCES.md`). No code/data yet.
- **2026-07-16 (later)** — **Step A done.** Deep-research follow-up confirmed novelty; SEdb 3.0 fully checked
  (`STEP_A_NOVELTY_CHECK.md`). Also scanned the local machine: the **DepMap join is already 100% local & current**
  (`~/Dropbox/manuscripts/0.datasets_visualizations/DepMap/`: `OmicsCNGeneWGS.csv` = CN correction input,
  `CRISPRGeneEffect.csv`, `protein_coding_expr/`, `Model.csv`); the **pyCaCTS JSD engine** (`../pyCaCTS/pycacts/`)
  is reusable for Phase 4; the local 15 GB `chip-atlas/` cache is **peak-call BEDs, not per-experiment bigWigs**,
  so the H3K27ac signal compendium still needs a targeted per-SRX ChIP-Atlas download.
- **2026-07-16 (later still)** — **Step B done → GREEN** (`pilot/`). Both de-risking steps now clear. The
  manifest is fully recoverable from the local AllCell bed (`His.ALL.50.H3K27ac.AllCell.bed.gz` embeds SRX +
  cell line + GSM + antibody per peak — 11,538 H3K27ac experiments). Quantile normalization defuses cross-study
  batch. **Next: Phase 1 (data assembly) — the full ~1,789-experiment manifest + DepMap join, likely on Roihu
  (`$SECACTS_CSC_PROJECT`; scaffold in `pilot/roihu/`).**
- **2026-07-17** — **Methods & tooling survey** done (`DESIGN.md` → Methods & tooling decisions). Two changes
  from pilot defaults: normalization **quantile → S3norm/S3V2** (unconstrained quantile manufactures false
  positives at scale); add **Cistrome DB** QC as a filter atop ChIP-Atlas breadth. SE calling = **ROSE + cSEAdb
  ≥25% union at constituent resolution**; **evaluate LILY** (calling + CN in one tool). No-CN lines →
  **HMCan/Control-FREEC from the ChIP input** (verified). Join picture measured: **~318 CN-correctable lines**
  floor (`DATA_SOURCES.md`). Cellosaurus must be fetched (local stub empty). Process note: ask before
  >500k-token runs (saved to Claude memory).
- **2026-07-17 (later)** — **Phase 1 metadata slice done** (`phase1/`). Fetched Cellosaurus (117 MB) +
  `experimentList.tab` (329 MB) to local caches (`CATALOG.md` updated); wired the Cellosaurus→DepMap join. Final
  compendium picture: **11,827 H3K27ac hg38 → 2,917 QC-pass, CN-correctable experiments on ~332 lines**,
  lineage-diverse, full HGSOC panel, synonyms verified (OVCAR-3→CVCL_0465→ACH-000001). Remaining Phase 1 = stage
  the QC-pass bigWigs (Roihu) → then SE calling + quantification.
- **2026-07-18** — **Phase 2 designed, de-risked, and grid built** (`phase2/PULL_DESIGN.md`). Four things
  settled, mostly by empirical checks on Roihu: **(1)** the bigWig pull is a **stream/batch** design (download →
  quantify over a fixed grid → delete) that stays under the 250 GiB scratch, so no quota request — SEs are
  recovered as sums of grid rows (single pass, each bigWig fetched once). **(2)** **Roihu compute nodes have
  direct internet** (verified — a departure from Puhti/Mahti), collapsing the pull to a plain SLURM array (no
  login-node feeder). **(3)** The **grid** is built from **per-SRX ChIP-Atlas peak BEDs fetched on Roihu**
  (`21_grid_from_persrx.sh`) — not the 107 GB AllCell download (⅔ discard); compared thresholds and chose
  **`bed20`/Q<1E-20 → 505k regions, 11% genome, full 2,917/2,917 coverage**. **(4)** **SE caller = `ROSE2_callSuper`
  on bigWig-quantified signal**, CN-agnostic; **LILY rejected** (BAM-first via HMCan, CN non-optional, needs
  matched Input). Gotcha: ChIP-Atlas uses two threshold numberings (bulk-filename `NN` vs site `10×NN`) and
  serves downloads uncompressed. **Next: scaffold `array.slurm`, then run the pull → signal matrix.**
- **2026-07-19** — **The entire local pipeline built and validated end-to-end on the 13-sample pilot.** Pivoted
  the SE caller from "shell out to ROSE2" to **`cnrose/`** — our own bigWig-native, CN-aware, ROSE-style caller
  (pure pyBigWig+numpy+scipy, no R/deepTools; the idealized "CNV-aware ROSE" LILY only approximates), **validated
  bit-for-bit vs ROSE2** (SE Jaccard 1.0000) and ~3.6x faster than `multiBigwigSummary`. Built the **CN provider
  layer** (source-agnostic `CNProvider`/`CNTrack` + separable `correct()`; DepMap gene-level backend) — CN
  correction demotes amplicon SEs (MCF7 p=1.5e-97) *and* rescues amplicon-masked true SEs (calling-time is
  amplify-only; scoring-time symmetric). Built **`aggregate.py`** (union catalog + dual matrices, exact) and
  **`score_pilot.py`** — hierarchical CaCTS JSD over the full Oncotree hierarchy, reusing `../pyCaCTS`; recovers
  SOX17/ESR1, the CN ablation labels textbook amplicons (BCAT1/BCAS1). Findings: symmetric divide-ratio boosts
  deletions (-> amplify-only at calling); specificity **counts** are batch/panel-sensitive (need S3norm); the
  SE->MTF link is still proximity-only (Phase-6 association scoring is the next science step). **The method is
  de-risked; scale (Roihu pull + S3norm + balanced panel) is what sharpens it.**
- **2026-07-19 (later)** — **S3norm built + a decision-changing negative result** (`phase2/s3norm.py`).
  Implemented S3norm's two-parameter power transform, fitted on the **fine grid** and re-summed in
  `aggregate.py --norm s3norm` (the transform is nonlinear, so grid resolution is the correct place; the SE
  catalog is untouched). Anchors had to move from MACS2 calls to **signal quantiles** — per-sample peak counts
  span 53–23,528, so call-based anchors collapse for weak samples and the threshold is itself SNR-dependent.
  It demonstrably works: fitted **B tracks dynamic range (Spearman −0.956, p=3.4e-07)** and cross-study scale is
  equalized (HT29/HCT116 co-called median log2 ratio **−1.315 → −0.043**). **But S3norm requires a QC gate** —
  on low-enrichment samples it fits B>1 and amplifies noise into apparent specificity; ungated it was *worse
  than quantile* (Breast ESR1 #5 → #21). Added `aggregate.py --min-peaks` (≥2000; natural gap at 2046 vs 361);
  gated, biology matches quantile (ESR1 #2, SOX17 #8) with scale fixed. **Correction to the 2026-07-19 entry
  above:** the specificity-**count** disparity is *not* a normalization artifact and S3norm does not fix it —
  three causes were tested and falsified (cross-study scale: fixed, counts persisted; sample QC: OVCAR3 24 vs
  Kuramochi 183 among good samples; private-SE-call counts: no correlation, p>0.6). Counts are governed by the
  left-tail mass of each group's JSD distribution under a per-group empirical null — i.e. **panel geometry and
  the statistic, not the normalizer**. "Trust rankings, not counts" therefore still stands, for a new reason,
  and a balanced panel (not a better normalizer) is what would lift it. **Next: scaffold `array.slurm` (gated).**
- **2026-07-20** — **Comparable specificity counts: the BH scope, not the null** (`phase2/specificity.py`).
  Tested the full 2×2 of null calibration × BH scope on the pilot. Sharing ONE multiple-testing budget across
  all SE × group tests compresses the per-group count spread (max/min **248 → 3.9** at Disease/Subtype,
  **10.3 → 4.5** at line), eliminates pathological zero-call groups (Breast 0 → 12), and *strengthens* both
  identity genes (SOX17 −1.73 → −2.36, ESR1 −1.51 → −2.62). The intuitive fix — a pooled **global null** — is
  the wrong lever: over-conservative for genuinely tight groups, it drops SOX17 below significance (−0.18).
  Now the default (`--fdr-scope global`); `--fdr-scope pergroup` reproduces pyCaCTS exactly (asserted,
  max|diff|=0). Rankings are invariant throughout. Cumulative: HGSOC count gap **11.2× → 7.6× → 3.4×**.
  **Next: `phase2/array.slurm` + the per-sample retention policy (what each task must emit before the bigWig
  is deleted).**
- **2026-07-20 (later)** — **Retention policy for the destructive pull** (`phase2/PULL_DESIGN.md §8.5`). The
  pull deletes each bigWig after processing, so every per-sample artifact was audited against one question:
  *could we reconstruct this later without the bigWig?* Three things are genuinely irreversible — off-grid
  signal, the **SE calling-time signal**, and signal-based QC metrics. The second was a surprise: `cnrose`
  quantifies stitched regions straight from the bigWig, and stitched regions span the gaps between constituent
  peaks, so `.enhancers.tsv` `SIGNAL` and the sum of overlapping grid rows agree for only **0.3%** of regions
  (median rel. diff 109%) — the calling signal cannot be re-derived from the fine matrix and must be retained.
  Decisions: quantify on **both** `grid.10` and `grid.20` (verified `grid.20` is 100.00% bp-covered by
  `grid.10`, but rows are not aligned so re-aggregation is inexact — emit both while the bigWig is open, which
  is free against a download-bound pipeline); emit a **fixed off-grid background bin set** (true background,
  FRiP and the S3V2 variance step are otherwise unreachable from a peak-only grid); emit a per-sample
  **`qc.json`** (S3norm requires a QC gate and signal-based QC needs the bigWig); store **float32 binary, not
  TSV** (15 MB → 2.0 MB/sample; 110 GB → 14.9 GB at scale). Total durable **~17 GB**. Residual risk stated
  honestly: genome-wide input-inferred CN (Phase 3) would still force a re-pull unless the background set is
  widened to a coarse 10 kb genome-wide binning (+3.6 GB) — recommended, since ~94% of expansion lines lack
  DepMap CN. **Policy is designed, NOT implemented:** `cnrose` needs multi-grid `--grid`, `--bg-bins`, `.f32`
  output, and a `qc.json` writer before `array.slurm` is worth scaffolding.
- **2026-07-20 (later still)** — **Retention policy revised: a genome-wide binned archive retires the
  irreversible grid choice** (`PULL_DESIGN.md §8.5`, evidence in `phase2/analysis/binned_storage_eval.py`).
  Tested whether storing a reduced/binned copy of each bigWig could remove the peak-threshold decision
  entirely. Result is a clean split: signal is an integral, so reconstruction error enters only at partial end
  bins — contiguous spans come back at **0.04–0.05%** (100 bp bins) but a single ~685 bp grid row is only ~7
  bins wide and lands at **2.8–3.2% median, 14–21% p90**. Decisively, the per-row error is **biased** (signed
  ≈ absolute, +3.1 to +4.5%), so it *accumulates* across an SE's rows rather than cancelling: the atlas SE-sum
  reaches 0.5–1.0% median but 11% max. So bins are an **archive, not a replacement** — primary quantification
  stays exact from the bigWig (cnrose is validated bit-for-bit; 1–4% error would forfeit that for a 2 MB
  saving). Final policy: exact tier (grid.20 + grid.10 float32, `.enhancers.tsv`, `qc.json`, calls,
  provenance ≈ 16 GB) **plus** a genome-wide 100 bp float16 archive (~65 GB) → **~81 GB**, park the archive in
  CSC Allas. The archive also resolves the previous draft's main residual risk: genome-wide input-inferred CN
  (Phase 3) needs broad coverage a peak-only grid could never give, and matters for the ~94% of expansion
  lines lacking DepMap CN. Also corrected: ChIP-Atlas has **four** thresholds (bed05/10/20/50); bed05 was
  omitted from the §8 comparison and was never built.
- **2026-07-20 (evening)** — **Retention policy IMPLEMENTED in cnrose, and it costs nothing.** `--grid` is now
  repeatable, plus `--signal-format f32`, `--bins 100` (genome-wide float16 archive + layout json), `--qc`
  (`qc.json`). The naive implementation re-read the same coverage runs three times and took **67 s/sample**,
  more than the ~30 s download — which would have made the whole policy expensive rather than free. Replaced
  by `cnrose.io.scan_bigwig`, which packs each chromosome's runs once and evaluates every grid and every bin
  edge against it: **17 s/sample for the full retention set, less than the 16 s the old single-grid path
  cost.** Verified bit-exact against the original `quantify`, and the ROSE2 contract re-passes **13/13 at
  Jaccard 1.0000** after the refactor. Two bugs caught by testing rather than reasoning: (1) a bigWig's `mean`
  is over COVERED bases only, so cnrose signal is (integral/covered)×span, not the plain integral — the first
  draft disagreed by up to 99%; (2) the QC dynamic-range gate statistic must be computed over the GRID column,
  since genome bins are ~63% zero → median 0 → `null` for every sample. Fidelity re-measured on **all 13**
  pilot bigWigs (was 2): SE-sum error at 100 bp is 0.29–0.98% median but **4.4–20% worst case**, and the
  signed bias is positive in **13/13** — archive-not-replacement confirmed on the full pilot. Measured
  durable footprint **~60 GB** (was estimated 81 GB). **Next: `array.slurm` (download → cnrose → provenance →
  delete) + `provenance.json`; then scope BU/wall-clock for the go/no-go.**
- **2026-07-20 (night)** — **The pull is BUILT and locally verified; only the go/no-go remains.**
  `phase2/array.slurm` (chunked SLURM array, 16 samples/task) → `scripts/40_pull_one.sh` (per-sample:
  download → one cnrose pass → provenance → delete; resumable via `.done`, atomic publish, per-task
  `.failed` markers) → `reduce.slurm` (barrier reduce producing all four catalogs, refusing to run on an
  incomplete pull unless `FORCE=1`) → `scripts/41_stage_and_submit.sh` (laptop-side, dry-run by default).
  `roihu/env.yaml` is now just pyBigWig+numpy+scipy — no R, no deepTools, nothing to recompile for ARM.
  Verified before spending any Roihu BU: a real pull of SRX067407 from the live ChIP-Atlas endpoint
  reproduces the local SE calls exactly (724 SEs, cutoff 1707.72); resume skips in 2 ms; a bad SRX writes
  `.failed` and exits 0; no staging dirs leak; and the reduce over the sharded `.f32` layout reproduces the
  TSV-path atlas **bit-for-bit** (max diff 0). Budgets are now measured rather than estimated: **~60 s/sample,
  ~50 core-hours, ~3.0 h wall at `%16`, ~60 GB durable, 1 CPU + 2 GB**. Cutting peak RSS 2.74 → 1.4 GB via
  chunked interval reading was a billing decision, not a tidiness one (CSC bills `max(cores, mem/mem_per_core)`).
  **Next: confirm the Roihu BU rate + per-core memory share, run a 2-task test array and check `seff`
  (which also measures the real compute-node download rate), then the go/no-go on the full 2,917.**
- **2026-07-21** — **Read the actual CSC billing rates; storage, not compute, is the constraint.**
  Roihu core partitions bill `max(0.75 BU/coreh × cores, 0.375 BU/GiBh × mem) × hours`. The two terms are
  equal at **exactly 2 GiB/core**, so the pull's `--cpus-per-task=1 --mem=2G` is precisely optimal — 4 GiB
  would double the bill and 8 GiB quadruple it, while asking for less saves nothing. The chunked-interval
  work (2.74 → 1.4 GB peak RSS) is what moved the job into that free headroom. **Whole pull ≈ 36 BU** — cheap.
  But Roihu bills scratch at **6 BU/TiB-hour from the first byte** (no free tier, unlike Puhti's 1 TiB), so
  the ~58 GB retention set costs **7.6 BU/day and overtakes the entire cost of the pull after ~4.8 days**;
  the genome-wide archive is ~75% of that. Added `scripts/50_fetch_and_clean.sh` (fetch → verify by count +
  md5 spot-check → delete from scratch, with `ARCHIVE_ONLY=1` to move just the bulky part) and promoted it
  from cleanup to a pipeline step. Also fixed a latent OOM: the transient bigWig must not land in
  `$TMPDIR`/`/tmp`, which on CSC is often RAM-backed tmpfs and would count against the 2 GiB limit.
  **Next: the go/no-go — a 2-task test array + `seff` to confirm compute-node download rate, then launch.**
- **2026-07-21 (staged rollout tooling + rehearsal)** — Built the staged-launch machinery so the pull rolls
  out smoke(2) → small(32) → medium(256) → full(2917) rather than cold: `scripts/42_preflight.sh`
  (read-only Roihu checks — SSH, partitions, quota, toolchain, staged inputs, and a compute-node internet
  probe), `43_run_stage.sh` (submit a stage / `report` progress+failures+QC spread+BU/day), and
  `44_local_rehearse.sh` (run the real `array.slurm` body locally against the live endpoint, no Roihu/BU).
  Rehearsal earned its keep immediately, catching two bugs: (1) a whole-run config error (missing peak dir)
  made every task write `.failed` and masquerade as a data problem — `array.slurm` now fails fast on missing
  manifest/peak-dir/grid while per-sample issues still continue; (2) the `report` stage parsed the QC
  dynamic-range with `sed 's/[^0-9.]//g'`, which leaked the `99` from the key `p99` and turned 29.41 into
  9929.41 — the QC histogram used to decide panel scope was silently wrong. Both fixed and re-verified.
  Runbook: `PULL_DESIGN.md §8.6`. **BLOCKED only on the user re-signing the 24 h Roihu SSH certificate**
  (expired 2026-07-19); everything else is ready and rehearsed locally.
- **2026-07-21 (session close) — PHASE 2 COMPLETE.** The reference SE atlas is built and validated on Roihu:
  agnostic **2,916 × 43,931** + S3norm **2,136 × 42,943**, `max|err|=0`. Full pull ran zero-OOM (2,916/2,917;
  SRX20868733 retired), ~43 BU total. Storage rationalized 88 → 27 GB (results gzipped 2.1 GB → 734 MB;
  genome archive coarsened 100 bp → 1 kb, verified 0.00% at CN scale, 100 bp chucked). The staged rollout
  (smoke→small→medium→large→full) + local rehearsal caught six bugs pre-full-run; the retention set proved its
  worth twice (QC recompute + archive coarsening, both no-re-pull). **Next: Phase 3 (CN inference) — design
  first, validate inferred CN against DepMap-CN lines. Tail: fetch atlas off scratch; regen CN calling-time
  atlases from retained data; drop SRX20868733 from the manifest.**
- **2026-07-22 — ATLAS LANDED + PHASE 4 AT SCALE + PHASE 6 STARTED.** (1) **Atlas fetched off Roihu scratch**
  to `an external SSD` (26 GB `out/` + 734 MB `results/`), verified: 2,916/2,916 `.done`, all 9
  retained file types 2,916x, 20 grid columns md5-identical to Roihu, all 8 `results/` files md5-identical,
  `gzip -t` clean. `results/` mirrored to `phase2/results/` (Dropbox) with the small artifacts tracked in git
  — so catalogs/presence/s3norm-params have a third copy on GitHub. The two signal matrices (423 + 306 MB)
  exceed GitHub's 100 MB limit and live on Dropbox + SSD only. **Scratch deliberately NOT cleaned** (user
  accepts ~3.5 BU/day; it keeps `out/` from being single-copy on a removable disk). (2) **Phase-4 scoring run
  at scale** — see Phase 4 above; a silent total-FDR failure was found and fixed en route (gotcha 64), which
  is the single most important thing to remember from this session. (3) **CN ablation at scale** — see
  Phase 5; it recovers MYCN/MYC amplicons unprompted, which is the strongest validation the CN layer has had.
  (4) **Phase-6 concordance bridge stage 1** — see Phase 6; 6–8x enrichment with the shuffle and distance
  controls behaving correctly. **Next: bridge stage 2** (all specific SEs, all genes in a window). Tail still
  open: regen the CN calling-time atlases from retained data; drop SRX20868733 from `data/pull_srx.txt`.
