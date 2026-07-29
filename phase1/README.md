# Phase 1 — data assembly

The full H3K27ac × DepMap manifest. Started 2026-07-17 with the **cheap metadata slice** (done); the heavy
signal work (SE calling + quantification over ~thousands of bigWigs) goes to Roihu later.

## Done (metadata harmonization — local, cheap)

- `scripts/10_fetch_metadata.sh` — fetched **Cellosaurus** (`cellosaurus/cellosaurus.txt`, 117 MB) and the
  authoritative **ChIP-Atlas `experimentList.tab`** (329 MB) into the shared local caches (see `CATALOG.md`).
- `scripts/11_harmonize_and_join.py` — the join: `experimentList.tab` (hg38 H3K27ac + QC) → Cellosaurus CVCL →
  DepMap `RRID`/Oncotree lineage + CN/dependency flags → **`data/phase1_manifest.tsv`** (11,827 rows).

### Coverage (see `../DATA_SOURCES.md` for the table)

11,827 H3K27ac hg38 experiments → 9,579 CVCL-mapped → 3,846 DepMap-joined (478 lines) → 3,128 on 332 CN lines
→ **2,917 QC-pass + DepMap + CN** (the analysis-ready core). Synonyms resolve correctly (OVCAR-3 → CVCL_0465 →
ACH-000001, HGSOC, CN✓). Lineage-diverse; full HGSOC panel present.

**Cell-line identity & status — start here.** [`CELL_LINE_IDENTITY.md`](CELL_LINE_IDENTITY.md) is the plain-English
front door: how each messy ChIP-Atlas name is resolved to one canonical CVCL identity, and what we can/cannot do
with each line (call SEs / CN-correct / score by lineage). Headline: 722 identity-resolved → exclude 91
(non-cancer / non-human homonyms / DepMap-Non-Cancerous) → **631 investigable human cancer lines, every one with
a lineage**. Full browsable per-line table: [`data/cell_line_status.tsv`](data/cell_line_status.tsv)
(`scripts/13`).

**Lineage recovery** (`scripts/14_cellosaurus_lineage.py`): DepMap covers 460 lines; the rest were lineage-less
only because the join never read Cellosaurus' disease field. Recovered **171** lineages from Cellosaurus `DI`
(NCIt) via an empirical NCIt→Oncotree crosswalk (`data/ncit_oncotree_crosswalk.tsv`), closing the old 259-line
gap. This work also caught **two identity bugs now fixed in `scripts/11`**: a `cell_desc="NA"` **sentinel sink**
(4,206 unrelated experiments collapsed into one hamster CVCL) and non-human **homonym** misresolutions (HAP1,
KG-1, PC-3, …) — fixed by guarding sentinel keys and canonicalizing to the DepMap RRID.

**Per-line CN provenance** (deep-dive behind the buckets) is in [`CN_COVERAGE.md`](CN_COVERAGE.md), produced by
**`scripts/12_cn_source_audit.py` → `data/cn_source.tsv`**. Among the 631 included lines: CN available 625/631
(99%); buckets A 324 (in atlas) / B 248 lookup / C 53 infer / D 6 none.

## Next (heavier — sequence per `../ROADMAP.md` Phase 1–2)

- [ ] Stage the QC-pass bigWigs (the analysis set) — Roihu `$SECACTS_CSC_PROJECT` for the scale pull.
- [ ] ROSE per-sample SE calls → cSEAdb-style union catalog (≥25% overlap, constituent resolution); evaluate LILY.
- [ ] Quantify H3K27ac over the union catalog (`multiBigwigSummary`) → region × experiment matrix.
- [x] S3norm normalization (not quantile — see `../DESIGN.md`) — built 2026-07-19 as `../phase2/s3norm.py`,
      applied via `aggregate.py --norm s3norm`. **Requires `--min-peaks` QC gating** (the two are coupled).

## Run

```bash
cd phase1/scripts
bash   10_fetch_metadata.sh      # idempotent; skips if caches present
python 11_harmonize_and_join.py  # -> data/phase1_manifest.tsv + coverage report
```
