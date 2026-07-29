# Data sources & assembly plan

What to pull, how the pieces join, and the practical constraints. Nothing is downloaded yet.

## Primary compendium — ChIP-Atlas H3K27ac

- **Scope:** ~**1,789** H3K27ac cancer-cell-line experiments, spanning **~295 unique "cell type" values** (from
  the ChIP-Atlas experiment CSV, filter "H3K27ac" + "cancer cell line"). ~6 experiments per cell type on average
  (highly uneven — some types many, some 1).
- **Role:** the reference atlas each lineage is scored against. This breadth is ~10× cSEAdb's 28 cancer types and
  is itself a headline novelty.
- **Form:** ChIP-Atlas provides uniformly processed peak calls (MACS2 `bed05`, q<10⁻⁵) and read-density bigWigs
  per experiment. **SE calls elsewhere (SEdb/dbSUPER) store regions, not signal** — so build a union region
  catalog and **re-quantify H3K27ac signal from the ChIP-Atlas bigWigs/BAMs** over it.
- **Access:** ChIP-Atlas bulk download (per-experiment bigWig/bed by SRX; the experiment metadata table drives
  the manifest). Large — plan staged downloads (and HPC if needed; see the parent project's `csc-puhti-job`
  workflow experience).

## Copy number, expression, dependency — DepMap (joins by cell line)

Highest-value add because one source supplies three things, all keyed to cell line:
- **Copy number** — GISTIC/segment-level CN per line = the **correction input** (divide H3K27ac by local
  copy-number ratio; Su/Chen et al. 2025 recipe). For lines absent from DepMap, fall back to **CNVkit from the
  ChIP/ATAC input-DNA** track. **[open]** impute vs exclude no-CN lines.
- **Expression** (`OmicsExpressionProteinCodingGenesTPMLogp1`) — the SE's **target-gene readout**; lets an
  H3K27ac-specific SE be cross-checked against the *gene*-CaCTS score (double support).
- **CRISPR dependency** (Chronos) — **function**: does the specific SE's target gene actually matter in the
  lineage. "Lineage-specific + amplification-independent + drives a dependency" = a strong nomination.

## Optional corroborating marks/factors — also ChIP-Atlas (same platform, joinable by cell line)

Use as *features where available*, never hard requirements (see the breadth-vs-depth caveat):
- **Histone marks:** H3K4me1 (enhancer), H3K4me3 (separate SE from active promoter), H3K27me3 (repression — a
  truly specific SE is H3K27ac-high in the target *and* Polycomb-marked elsewhere). A lightweight chromatin-state
  layer.
- **SE-defining coactivators:** BRD4, MED1/Mediator, EP300, CDK7 — corroborate SE identity, guard against
  H3K27ac artifacts.
- **Lineage master-TF ChIP** (PAX8, SOX17, WT1, MECOM, …) — a lineage-specific SE should be bound by the
  lineage's masters; validation + feature + the bridge to the CRC layer.
- **ATAC/DNase accessibility** (ChIP-Atlas accessibility arm) — orthogonal, less antibody-batch-prone activity.

## Metadata harmonization (a prerequisite, not optional)

ChIP-Atlas "cell type" is free text (synonyms, sublines, mislabels). Map the strings to a controlled
vocabulary so lineages group correctly *and* join to DepMap CN/expression/dependency:
- **Cellosaurus** — cell-line identity/synonyms/RRIDs. **Bridge:** ChIP-Atlas cell-type string → CVCL_* →
  DepMap `RRID` (DepMap `Model.csv` carries the CVCL directly). The **`cellosaurus` R package `mapCells()`**
  maps free-text names / DepMap IDs / ATCC IDs → CVCL in one call. Caveat: RRID↔model isn't always 1:1
  (e.g. SJRH30 & RH30 → CVCL_0041) — join with care. ⚠️ **local `cellosaurus/` is an empty stub → must fetch.**
- **OncoTree / DepMap lineage** — cancer-type grouping for the lineage-level specificity (in `Model.csv`:
  `OncotreeLineage`, `OncotreePrimaryDisease`, `OncotreeSubtype`, `DepmapModelType`).
- This mapping is where these projects usually bleed time; budget for it explicitly.

### The join picture (Cellosaurus-bridged + QC, 2026-07-17)

Built by `phase1/scripts/11_harmonize_and_join.py` from `experimentList.tab` → Cellosaurus CVCL → DepMap:

| | count |
|---|---|
| H3K27ac experiments (hg38, authoritative `experimentList.tab`) | 11,827 |
| → mapped to a Cellosaurus **CVCL** | 5,373 |
| → joined to a **DepMap model** | 3,846 on **478 lines** |
| → on lines **with copy number** (`OmicsCNGeneWGS`) | 3,128 on **332 lines** |
| → **QC-pass** (≥5M reads & ≥1k peaks) **and** DepMap+CN | **2,917 experiments** |

DepMap side: 2,105 models (1,966 with CVCL RRID), **1,118 with WGS CN**, 1,186 with CRISPR dependency. Binding
constraint is H3K27ac availability, not CN. The analysis-ready core — **~2,917 QC-passing, CN-correctable,
lineage-diverse experiments on ~332 cell lines** — dwarfs cSEAdb's 60 lines / 28 types. Cellosaurus resolves the
synonyms string-match missed (verified: **OVCAR-3 → CVCL_0465 → ACH-000001 → HGSOC, CN✓**; K-562, NB-4, all
MDA-MB-*). Full HGSOC panel present (OVCAR-3/4/5/8, Kuramochi, HEYA8, Caov-4, SK-OV-3…). CN lines are
lineage-spread (Lung 60, Lymphoid 47, Breast 33, Stomach 25, Myeloid 23, CNS 16, …).

## Authoritative metadata & a QC-reuse option

- **`experimentList.tab`** (ChIP-Atlas) is directly downloadable (~328 MB current endpoint / ~189 MB legacy
  `dbarchive`) and is the authoritative Phase-1 manifest — SRX, genome, antigen, cell-type-class/type, **plus
  per-experiment QC (read count, % mapped, % dup, peak count)**. Use its QC to quality-filter the H3K27ac set
  (the "depths" batch axis the pilot flagged). *(For the pilot, the local `AllCell.bed.gz` already yielded the
  manifest without this download.)*
- **Cistrome DB v3.0** (*NAR* 2024) — ~45,000 uniformly processed human ChIP/ATAC samples via the CHIPS
  pipeline **with a standardized 6-metric QC framework** + regulatory-potential scores. **Plan:** breadth from
  ChIP-Atlas, QC rigor from Cistrome — use Cistrome's QC to filter/cross-check, not necessarily to replace the
  ChIP-Atlas signal quantification. (See `DESIGN.md` → Methods & tooling decisions.)

## Constraints & notes

- **Breadth vs depth:** H3K27ac alone → ~295 cell types. Requiring H3K27ac+H3K4me1+ATAC → far fewer. Keep H3K27ac
  primary; treat extra marks as optional features (impute/down-weight where missing).
- **Cross-study normalization** is the dominant technical risk (many labs/antibodies/depths) — see `DESIGN.md`
  and the ROADMAP normalization pilot. ChIP-Atlas's uniform pipeline mitigates but does not remove batch.
- **Genome build:** ChIP-Atlas provides hg38; keep everything hg38.
- **Local reference:** the CaCTS/pyCaCTS implementation is at `../pyCaCTS` — a worked example of the TF layer
  this complements.
