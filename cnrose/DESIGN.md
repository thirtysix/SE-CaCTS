# cnrose — design

**A bigWig-native, external-CN, ROSE-style super-enhancer caller.**

**Status:** design agreed 2026-07-18 (all four forks A–D locked with the user). No code yet. This doc is the
build spec; it supersedes the ROSE2-shell-out sketch in `../phase2/PULL_DESIGN.md` §8 for the SE-calling step.

Read first: `../phase2/PULL_DESIGN.md` (the pull this is the per-sample engine for), `../ROADMAP.md`
(Phases 2–4), `../DATA_SOURCES.md` (CN sources), and the CN-source tiering computed 2026-07-18 (below).

---

## 1. Why this exists (and why not stock ROSE2 / LILY)

The Phase-2 pull needs a per-sample super-enhancer caller that runs over **2,917–~8,800 ChIP-Atlas H3K27ac
bigWigs** on Roihu (ARM) and can optionally correct for copy number **from an external source**. Neither stock
tool fits cleanly:

- **ROSE2** is BAM-first; our compendium is bigWig-only. Its callSuper core *can* eat a pre-quantified table,
  but wiring `ROSE2_callSuper.R` + two `multiBigwigSummary` calls + a `bedtools` stitch per sample is heavy,
  R-on-ARM is fiddly to containerize, and it has no CN hook.
- **LILY** is the right *idea* — ROSE-style calling on CN-normalized signal — but it hard-couples three things
  it only couples because it has **no external CN**: it estimates CN from reads via **HMCan** (BAM-first, needs
  a matched Input BAM, CN non-optional). We *have* external CN (measured for 64% of lines, input-inferrable for
  another 30% — see §4), so we can **decouple** what LILY fuses.

`cnrose` is the idealized form LILY only approximates under its no-external-CN constraint: **ROSE's algorithm,
bigWig-native, with copy number supplied from the best available source and correction that is optional and
separable.** It replaces `ROSE2_callSuper.R` + 2×`multiBigwigSummary` + `bedtools` in one `pyBigWig`+`numpy`
process — leaner, ARM-trivial to deploy, and a reusable tool (Phase-7 artifact, sibling to `../../pyCaCTS`).

**Framing for the paper:** we decouple LILY's CN-aware idea into a scalable, bigWig-only, external-CN form; cite
LILY/ROSE as the lineage. `cnrose --cn none` is validated to reproduce ROSE2 (§7).

## 2. The organizing principle — factor along the axes that change

To not reimplement twice, split stable core from pluggable edges.

**Won't change → port once, validate against ROSE2, freeze:**
- constituent **stitching** (merge within 12.5 kb; optional TSS-break),
- the **tangent-line cutoff** (ROSE `calculate_cutoff`, ~20 lines),
- the **bigWig read primitive** (`pyBigWig`).

**Will change → pluggable from day one:**
- **CN source** (≥6 backends, some not yet built — §4),
- **CN entry point** (calling-time vs scoring-time vs both — §3.2),
- **orchestration** (local loop / SLURM array / in-process pool),
- **region set** (per-SRX stitched for calling; fixed grid for the matrix; genome bins as generic fallback).

## 3. Architecture

```
cnrose/
  io.py         quantify(bigwig, regions) -> np.ndarray        # THE one bigWig primitive (pyBigWig)
  stitch.py     stitch(peaks, window=12500, tss=None)          # ROSE port
  callsuper.py  cutoff(signal) -> thresh ; call(table) -> table# ROSE tangent port  <- validated vs ROSE2
  cn/
    base.py     CNProvider (interface), CNTrack (canonical), correct(signal, cn, model)
    depmap.py   DepMapMerged backend
    progenetix.py  Progenetix backend (CVCL-native Beacon API)
    cmp.py      CellModelPassports backend (absolute CN)
    cosmic.py   CosmicCLP backend
    infer.py    InputInferred (HMCan/Control-FREEC/CNVkit from matched input), CoverageInferred, NoCN
  pipeline.py   call_sample(bw, peaks, grid, cn_provider=None, correct_at="none") -> se_bed, signal_tsv
  cli.py        `cnrose call ...`
```

### 3.1 Abstraction 1 — `CNProvider`: source-agnostic CN

The caller never knows where CN came from. It asks a provider for a positional track keyed by cell line
(CVCL), and gets back a **canonical `CNTrack`**: segments → **linear ratio vs. sample ploidy** (1.0 = neutral).

```python
class CNProvider:
    def track(self, cvcl) -> CNTrack | None: ...   # None => this source has no CN for this line

class CNTrack:
    # segment-level: (chrom, start, stop, ratio_linear); ploidy stored for absolute<->relative conversion
    def region_cn(self, chrom, start, stop, agg="wlen") -> float: ...  # §6 aggregation
```

Each backend’s **only** job: resolve the join and emit a canonical `CNTrack`. New sources are **additive
backends**, never a core rewrite. Format heterogeneity is normalized *inside* the backend:
- relative ratio (DepMap, Progenetix log2→`exp`) → used directly,
- absolute integer CN (CMP, COSMIC) → `÷ ploidy` to relative.

A **chained provider** resolves per line in priority order (§4), returning the first hit and recording
`cn_source` + a confidence tier on the result.

### 3.2 Abstraction 2 — `correct()`: CN as a separable transform

One pure function, two call sites — so corrected-vs-uncorrected exists at **both** the catalog and the score
level with zero duplicated code:

- **calling-time** (LILY-style): correct the stitched-region signal *before* the tangent cutoff → changes which
  regions clear the bar → a **CN-corrected catalog** (fork D).
- **scoring-time**: correct the fine grid matrix / SE matrix downstream → **CN-corrected scores** on a fixed
  catalog (the existing Phase-4 plan).

The **fine grid matrix stays uncorrected and single-pass** regardless (CN is separable — the corrected matrix is
derived later by applying `correct()` per region with the proper CN join). So calling-time correction adds a
*second catalog*, not a second bigWig pass; the pull economics of `PULL_DESIGN.md` are unchanged.

## 4. CN-source tiers (measured 2026-07-18 on the 398 net-new lines)

Routing each net-new CVCL through Cellosaurus `DR` cross-refs and the ChIP-Atlas input audit:

| Tier | Source | Lines | % |
|---|---|---:|---:|
| **1 — measured CN** | DepMap-merged / Progenetix / CMP / COSMIC | 254 | 64% |
| **2 — input-inferred** | HMCan / Control-FREEC / CNVkit from a matched input bigWig | 119 | 30% |
| **3 — orphan** | coverage-only or CN-agnostic (flagged) | 25 | 6% |
| | **CN available (1+2)** | **373** | **94%** |

Per-resource reach among the 398: Progenetix 245, CMP 153, DepMap-xref-but-no-WGS 140, COSMIC 95, GDSC 98.
Matched-input coverage is high everywhere (current set 97%, expansion 92% of experiments / 88% of lines).

**Backend resolution order** (least friction / best quality first):

1. **DepMapMerged** — swap the WGS-only `OmicsCNGeneWGS.csv` we have for the **merged `OmicsCNGene.csv`**
   (+ segment-level `OmicsCNSegmentsWGS.csv` for positional CN). One download, same CVCL/`ACH-` join, hg38;
   recovers the 140 net-new lines DepMap already has via WES/SNP. **Cheapest win — do this first.**
2. **Progenetix** — CVCL-native (no ID hop), hg38, CC-BY, GA4GH Beacon v2 API. Biggest reach (245).
3. **CellModelPassports** — CVCL→`SIDM` hop; uniquely **absolute integer CN** (best for ploidy-aware
   correction). hg38.
4. **CosmicCLP** — registration/licence friction → top-up only (95), largely redundant with CMP/GDSC.
5. **InputInferred** — the 119: derive CN from the matched input bigWig (tool TBD — HMCan / Control-FREEC /
   CNVkit; `DATA_SOURCES.md` leans CNVkit-from-input). **Gate on input read depth** (manifest `reads`) — "has an
   input" ≠ "has a *good* input"; coarse on focal amplicons.
6. **CoverageInferred / NoCN** — the 25 orphans; flagged `cn_source=none`, appear in the agnostic catalog only
   (or with a low-confidence coverage estimate).

Join recipe (Cellosaurus `DR` token → native id): `DepMap; ACH-*` → `OmicsCNGene.ModelID`;
`Cell_Model_Passport; SIDM*` → CMP `model_id`; `Cosmic-CLP; <id>` → COSMIC `ID_SAMPLE` (== `Model.csv.COSMICID`
== GDSC id — **never** use `WTSIMasterCellID`); `Progenetix`/`cancercelllines; CVCL_*` → Beacon
`id=cellosaurus:CVCL_*`.

## 5. Decisions (forks A–D, locked 2026-07-18)

- **A. Home & name** — package **`cnrose`**, built inside SE-CaCTS (private repo; scooping risk), structured
  standalone so it can spin out later like `pyCaCTS`.
- **B. Correction model** — **log2-offset is primary** (propagates CN uncertainty, graceful with missing/noisy
  CN, matches the Phase-4 regression-with-offset choice and the HiC-DC+ precedent). Alternatives kept on record
  (§6.1) per the user's request, in case of trouble.
- **C. Multi-segment SE → one CN** — **length-weighted mean over the SE's overlap with each CN segment**
  (default `agg="wlen"`); `max` available as the conservative amplicon-killer (`agg="max"`). See §6.2.
- **D. Dual catalogs** — an **agnostic** *and* a **CN-corrected** SE catalog are **both first-class**
  deliverables (not one as a buried sensitivity). This is what justifies wiring calling-time correction.

## 6. Correction model detail

Signal per region `s` with region copy-number ratio `r` (relative to ploidy, 1.0 = neutral):

### 6.1 The model (and the alternatives we're keeping documented)

- **Primary — log2-offset:** work in `log2(signal+ε)` and subtract a CN offset `β·log2(r)` (a covariate, with
  `β` estimable rather than fixed at 1). Degrades gracefully when `r` is missing (offset→0) or noisy (shrink
  `β`), and is the same functional form as the Phase-4 scoring regression → one story end to end.
- **Cross-check — naive ÷CN-ratio:** `signal / r` (Su/Chen 2025 form). Simpler; kept as a documented
  cross-check (`model="divide"`). Discards uncertainty and is unstable at low `r` / high noise.
- **Other possibilities on record** (reach for only if the above misbehave):
  - **GC co-correction** — LILY does GC+CN jointly via HMCan; input-inferred CN (tier 2) may need a GC term.
    Measured-CN tiers (1) don't.
  - **Segment vs gene-level CN application** — prefer **segment-level** tracks (`OmicsCNSegmentsWGS`,
    Progenetix segments) over gene-level; gene-level is a fallback.
  - **Hard filter vs soft correction** — excluding high-CN regions entirely (vs correcting them) as an
    ablation, to bound how much the correction model itself drives results.
  - **Rank/quantile CN adjustment** — a non-parametric alternative if the log/linear forms distort the tail.
  - **HMCan-style joint read-level modeling** — rejected as the primary path (BAM/input-coupled), but the
    reference for tier-2 input inference.

### 6.2 Region CN from heterogeneous segments (fork C)

For SE region `[s,e]` overlapping CN segments `i` with per-segment ratio `c_i` over **overlap length**
`L_i` (bp of the SE that falls in segment `i`):

```
r_wlen = Σ_i (c_i · L_i) / Σ_i L_i         # default
r_max  = max_i c_i                          # conservative amplicon-killer (option)
```

The weight is the **overlap length within the SE**, not the full segment length — a large segment the SE barely
touches must not dominate. Worked example: 8 kb in CN-2 + 2 kb in CN-10 → `r_wlen = (2·8000+10·2000)/10000 =
3.6`; `r_max = 10`. (This is the `ROADMAP.md` open question "single CN ratio for an SE spanning heterogeneous
segments" — resolved for the caller; the grid rows keep per-region coords so scoring can revisit per-segment.)

### 6.3 Calling-time correction is ASYMMETRIC (amplify-only) — rehearsal finding 2026-07-19

Empirically (13-pilot dual-catalog rehearsal), symmetric ÷ratio/log2-offset at **calling** time misbehaves in
two ways; the caller therefore clamps **CN ≥ `cn_floor` (default 1.0 = amplify-only)** before correcting:

- **Deletion-boosting (artifact, fixed):** symmetric correction inflates *deleted*-region signal
  (`(sig+1)/cn` grows as `cn→0`), promoting deleted enhancers to false SEs (MCF7: 156 gained SEs, 96% in
  CN<0.9). A deleted enhancer must never become a super-enhancer at calling time. `cn_floor=1.0` leaves
  deletions untouched and removes this entirely.
- **Amplicon-masking rescue (real, kept):** in lines with strong focal amplicons the amplicon's huge H3K27ac
  *inflates the ROSE tangent cutoff*, suppressing detection of genuine SEs (THP-1 called only 312; Kuramochi
  535 — implausibly few). Correcting the amplicons down (demoted SEs at 7–11× CN) collapses the inflated
  cutoff and restores the true SE population (THP-1 → 774, Kuramochi → 1133). This is a *second* benefit of
  CN-corrected calling beyond removing false amplicon-SEs, and a reason the dual catalogs are co-primary.

**Scoring-time** correction stays SYMMETRIC (per-copy normalisation across lineages, `cn_floor≈0.1`) — the
two stages use the same `correct()` with different floors. Exposed as `--cn-floor`.

## 7. Validation contract (built alongside the core, not after)

`cnrose --cn none` **must** reproduce `rose2` on the 13 local pilot bigWigs before we trust any CN mode:
- same stitched regions (`-s 12500 -t 0` for H3K27ac, the standard stock-ROSE settings),
- SE-set **Jaccard ≥ ~0.98**, cutoff within a small tolerance,
- ROSE2 signal-table contract for reference: `REGION_ID CHROM START STOP NUM_LOCI CONSTITUENT_SIZE <signal>
  [<control>]`, then callSuper appends `enhancerRank`, `isSuper`; SE = `isSuper==1`; controlless = rank by the
  signal column, floored at 0.

If agnostic mode doesn't match ROSE2, CN mode isn't trustworthy. This is the guardrail against reimplementation
skepticism.

## 8. Build order (nothing thrown away)

1. **Agnostic core + validation** — `io` + `stitch` + `callsuper`, prove `--cn none` ≡ `rose2` on the 13 pilot
   bigWigs. *(No CN, no scope decision needed to start.)*
2. **CN interface + `correct()`** — wire `CNProvider`/`CNTrack`/`correct`, backend = **DepMapMerged** first
   (do the `OmicsCNGene.csv` swap). Add calling-time + scoring-time call sites → dual catalogs.
3. **Add backends as confirmed** — Progenetix, CMP, COSMIC, then InputInferred. Each slots into the frozen core.

## 9. How it plugs into the Phase-2 pull

- **`per_srx`** step becomes: `cnrose call --bw SRX.bw --peaks SRX.20.bed --grid ../data/grid.20.bed
  [--cn <provider>] [--correct-at calling] --out SRX` → writes `SRX.se.bed` (+ a CN-corrected `SRX.se.cn.bed`
  when a provider is given) and `SRX.signal.tsv` (uncorrected grid column).
- **`aggregate.py`** (Phase-2 barrier reduce) is unchanged: union catalog (cSEAdb ≥25% overlap, constituent
  resolution) + fine grid×sample matrix; run once per catalog (agnostic + corrected).
- **`array.slurm`** wraps `per_srx`; resumable via `SRX.done`; node-local scratch; `PULL_DESIGN.md` §5.1
  write-safety. The heavy 2,917–~8,800-task launch stays behind the separate BU/wall-clock go/no-go.
- Env: `pyBigWig` + `numpy` (+ `pysam`/an input-inference tool for tier 2). No R, no deepTools → trivial tykky
  container on Roihu ARM.

## 10. Open items to confirm before / during build

- Exact **26Q1 `OmicsCNGene.csv`** line count + download (portal/figshare were flaky; confirm from the file /
  `OmicsProfiles.csv`).
- **Progenetix** access: Beacon v2 REST route is reliable; verify `pgxRpi` R-client signatures if scripting.
- **Tier-2 inference tool**: HMCan vs Control-FREEC vs CNVkit-from-input; pick one + a read-depth QC gate.
- **Scope decision still pending**: whether to relax the pull from `has_cn==1` (2,917) to the expansion
  (~8,800). Needs the parallel **Cellosaurus lineage** extension (lineage currently rides on the DepMap join)
  and ~3× the Roihu pull — feeds the go/no-go. `cnrose` is built the same way either way.
