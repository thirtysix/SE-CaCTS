# Cell-line identity & status — what each line is, and what we can do with it

Plain-English front door for the cell lines. It answers two questions:

1. **Who is each line?** — identity resolution: turning messy ChIP-Atlas names into one canonical ID.
2. **What can we do with each?** — the data we hold, the copy-number we can get, and the scoring we can run.

**Full per-line table:** [`data/cell_line_status.tsv`](data/cell_line_status.tsv) (722 rows, opens in any
spreadsheet). This doc is the *summary + methods*; [`CN_COVERAGE.md`](CN_COVERAGE.md) is the copy-number
deep-dive. Numbers come from `scripts/11` (identity join), `12` (CN source), `14` (lineage recovery),
`13` (this table).

---

## The simple version

We start from **722 identity-resolved lines** that have H3K27ac ChIP-seq in ChIP-Atlas. We then **exclude 91**
that don't belong in a *human cancer* atlas, leaving **631 investigable lines — every one a human cancer line
with a tumour lineage.**

**Excluded (91):** 83 non-cancer (immortalised-normal lines: HEK293, HaCaT, BEAS-2B, …), 3 that DepMap itself
labels *Non-Cancerous* (MCF 10A, HK-2, BJ-5ta), and 5 **non-human homonyms** (a name that collided with a
mouse/rat entry and has no DepMap anchor to correct it — see the QC catches below).

The 631 sort into four readiness buckets by how quickly we can get their copy-number:

| Bucket | Lines (of 631) | Plain meaning | Copy-number |
|---|---:|---|---|
| **1 — Ready now** | 324 | already pulled and **in the built atlas**; CN from DepMap WGS | ✓ have it |
| **2 — One lookup away** | 248 | not pulled yet, but CN is **already measured in another database** | ✓ a lookup |
| **3 — Need CN inference** | 53 | not pulled; no measured CN, but has a matched input to **infer** from | ~ must infer |
| **4 — Limited** | 6 | no CN from any source → usable only "CN-agnostic" | ✗ none |

- **99% (625/631) have copy-number available in principle.** Only 6 lines can't be CN-corrected at all.
- **100% have a tumour lineage** — 460 from DepMap, **171 recovered from Cellosaurus** (see *Lineage*). The old
  "no lineage" gap is closed: every line that survives the cancer/human filter now has one.

## Capability at a glance (the 631 included lines)

| Bucket | Call super-enhancers | CN-correct | Score by lineage | In the atlas now? |
|---|---|---|---|---|
| **1** (324) | ✓ | ✓ **now** | ✓ | ✓ yes |
| **2** (248) | after a pull | ✓ after a DB lookup | ✓ | needs pull |
| **3** (53) | after a pull | after CN inference | ✓ | needs pull |
| **4** (6) | after a pull | ✗ never | ✓ | needs pull |

"After a pull" = the 307 bucket-2/3/4 lines are the **expansion set**; their H3K27ac signal isn't downloaded
yet. The atlas built in Phase 2 is **bucket 1 only** (324 lines / 2,916 usable experiments). Expanding to all
631 is a second Roihu pull.

## What we CAN do

- **Call super-enhancers for any of the 631** — SE calling is a within-sample threshold on the H3K27ac signal;
  no CN, no cross-sample data needed. (Done for the 324 in the atlas; the rest just need pulling.)
- **Copy-number-correct 625 lines** (buckets 1–3), once buckets 2–3 are wired up (a database backend for 2, an
  inference step for 3). Correction demotes amplicon-driven false super-enhancers and rescues true ones.
- **Score tumour-type specificity for all 631** — every included line carries a lineage, so all can enter the
  lineage/subtype scoring (today, on the 324 in the atlas; the rest once pulled).
- **Score at single-cell-line resolution for all 631** — identity is known for every line.

## What we CANNOT do

- **6 lines have no copy-number from any source** → CN-agnostic only (can't separate amplicon-driven false SEs
  from real ones). Hard limit unless a CN source appears.
- **The 307 expansion lines aren't pulled yet** — no atlas signal exists for them until a second Roihu pull.
- **Copy-number is *relative*** for most sources (DepMap, Progenetix); only CMP gives *absolute* integer CN, so
  ploidy-aware correction is limited outside CMP-covered lines.
- **Recovered *subtype* is best-effort.** For the 171 Cellosaurus-recovered lines the *lineage* is
  high-confidence (crosswalk lineage-purity median 1.00), but the finer Oncotree *subtype* is the modal
  subtype for that disease — trust lineage over subtype for these.

---

## How we resolve identity (the methodology)

**The problem.** ChIP-Atlas labels each experiment with a free-text cell name typed by many submitters.
"OVCAR-3", "OVCAR3", "NIH:OVCAR-3", "ovcar 3" are one line. We must collapse those to a single canonical
identity per line.

**The pipeline** (`scripts/11_harmonize_and_join.py`):

1. **Collect** every hg38 H3K27ac experiment from ChIP-Atlas `experimentList.tab` — **11,827** — each with a
   free-text `cell` name (+ `cell_class`, `cell_desc`) and per-experiment QC.
2. **Normalize** each name: lowercase, strip everything non-alphanumeric. `OVCAR-3` → `ovcar3`. Deliberately
   aggressive, so punctuation/spelling variants collapse. **Guarded:** a not-available sentinel (`NA`, `none`,
   `unknown`, …) or a 1-char token is never used as a match key (see the sink catch below).
3. **Build the dictionary from Cellosaurus** (the authoritative cell-line registry): parse all **168,970**
   entries; for each, take its primary name **and every synonym**, normalize, map to the accession
   (`CVCL_xxxx`). → **261,595** name keys.
4. **Resolve** each experiment by cascade: normalized `cell`, then `cell + cell_class`, then `cell_desc`; first
   hit wins → a **CVCL = the canonical identity**. Fallback: direct match to a DepMap cell-line name.
5. **Canonicalize to DepMap.** When a line is in DepMap, its DepMap **RRID is the authoritative CVCL** — it
   overrides the Cellosaurus name-match, which can land on a same-named homonym (see below).
6. **Result:** **5,373** experiments resolve to a CVCL; among QC-pass experiments, **722 distinct lines**. Of
   those, **463** are cross-confirmed by DepMap (Cellosaurus CVCL == DepMap RRID); **259** rest on the
   Cellosaurus name match alone.

## Lineage: DepMap first, then recovered from Cellosaurus

DepMap supplies an Oncotree lineage for the confirmed lines. The Cellosaurus-only lines were lineage-less
**only because the join never read Cellosaurus' own disease field** — not because the data was missing.
`scripts/14_cellosaurus_lineage.py` recovers it:

1. **Read Cellosaurus per line:** `CA` (category — cancer vs not), `OX` (species), `DI` (NCIt disease code + term).
2. **Build an NCIt → Oncotree crosswalk *empirically***, from the ~1,860 DepMap models that carry **both** an
   `OncotreeCode` and an NCIt disease — no external Oncotree file needed. Vote lineage first (robust), then the
   modal subtype within it. → `data/ncit_oncotree_crosswalk.tsv` (302 codes, lineage-purity reported).
3. **Assign a lineage** by precedence: `depmap` → `cellosaurus_ncit` (crosswalk hit) → `cellosaurus_term` (a
   small curated fallback for NCIt codes DepMap never used) → none.
4. **Gate inclusion** to *human cancer*: `CA == "Cancer cell line"`, species = human, and a lineage was
   resolved. Everything else is excluded with a reason.

Result: **171 lineages recovered** (155 crosswalk, 16 curated), lifting lineage coverage of the investigable
set to **631 (100%)**. Zero cancer lines were left unmapped.

## Identity confidence & the two QC catches

Building the lineage layer surfaced two real identity bugs — both now **fixed in the join**, not just flagged:

- **The sentinel sink (fixed).** `cell_desc` is the literal string **"NA"** (not-available) for thousands of
  under-annotated experiments, and Cellosaurus lists "NA" as a synonym of a Chinese-hamster line (CVCL_E4I6) —
  so the unguarded cascade sank **4,206 experiments** (hundreds of unrelated names: "Frontal cortex", "CD4+ T
  cells", "Tumor", …) into that one bogus "line". Guarding the resolver against sentinel/1-char keys removes
  the sink entirely; those experiments correctly become unmapped. (This is why CVCL-mapped dropped from 9,579
  to 5,373 — pure junk removal; the DepMap-anchored analysis was unaffected.)
- **Non-human homonyms (fixed for the 8 with a DepMap anchor).** The aggressive normalization matched human
  ChIP-Atlas samples to same-named **mouse/rat/rabbit** Cellosaurus entries — "HAP1", "KG-1", "M-14", "MEC-1/2",
  "OCUM-1", "Rh-4" all have famous human cancer lines *and* non-human homonyms. The **RRID-canonicalization**
  (step 5) fixes all 8 by overriding to the DepMap human accession. The same rule also caught **PC-3**, whose
  "PC3" name had matched an *oral-cavity SCC* homonym (CVCL_C8XA) — now correctly the prostate line
  `CVCL_0035`. The **5** non-human names with *no* DepMap anchor (AT1, AT2, 42D, L3.6, RH-3) can't be corrected
  and are excluded (`non_human_homonym`).
- **Remaining confidence tiers:** **463** lines are Cellosaurus+DepMap double-confirmed (high confidence);
  **259** are Cellosaurus-name-only — correct for the large majority, but the residual risk is the
  first-writer-wins collision that produced the homonyms above.

## The per-line file

[`data/cell_line_status.tsv`](data/cell_line_status.tsv) — one row per line. Columns:

| column | meaning |
|---|---|
| `cvcl` | canonical identity (`CVCL_xxxx`; DepMap RRID when the line is in DepMap) |
| `cell_line` | representative display name from ChIP-Atlas |
| `category` | Cellosaurus `CA` (Cancer cell line / Transformed / …) |
| `lineage` / `subtype` | Oncotree tumour type (recovered where DepMap had none) |
| `lineage_source` | `depmap` / `cellosaurus_ncit` (crosswalk) / `cellosaurus_term` (curated) |
| `n_experiments` | H3K27ac experiments we have for it |
| `identity` | `confirmed (Cellosaurus+DepMap)` or `Cellosaurus-only` |
| `data_status` | `in atlas` (pulled) or `needs pull` (expansion) |
| `cn_status` / `cn_from` | CN bucket + where it comes from |
| `tier` | A / B / C / D (= buckets 1–4) |
| `include` | `1` = in the human-cancer analysis; `0` = excluded |
| `exclude_reason` | `non_cancer` / `non_human_homonym` / `depmap_non_cancer` (blank if included) |

`data/lineage_resolved.tsv` carries the full lineage detail; `data/ncit_oncotree_crosswalk.tsv` is the
auditable NCIt→Oncotree map.

Regenerate the whole chain: `cd phase1/scripts && python3 11_harmonize_and_join.py && python3 12_cn_source_audit.py && python3 14_cellosaurus_lineage.py && python3 13_cell_line_status.py`.
