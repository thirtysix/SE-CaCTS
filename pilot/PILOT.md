# Step B — the normalization pilot

**Goal (the go/no-go for the whole project):** measure how badly cross-study batch contaminates ChIP-Atlas
H3K27ac, and whether a normalization makes **lineage** (biology) dominate **study/lab** (batch) in the signal.
If lineage structure survives normalization, SE-CaCTS is viable; if study batch dominates and won't normalize
out, the approach needs rethinking before any large data pull. See `../DESIGN.md` ("The main technical risk").

## Design — how we separate batch from biology

A balanced grid so *lineage* and *study* are **not** collinear:

- **~4 cancer lineages × ~3 cell lines each ≈ 12 H3K27ac experiments** (`02_select_pilot.py`).
- Lineages drawn from cancer cell lines with abundant ChIP-Atlas H3K27ac; **Ovary/HGSOC is always included**
  (OVCAR-3 etc.) as the reference lineage.
- **A cross-study replicate probe:** one cell line (e.g. K562 or MCF-7) represented by ≥2 experiments from
  *different studies* — the cleanest direct read on batch. If the same line from two labs clusters together →
  batch is mild; if it splits by lab → batch is severe.

The question each figure answers: **do samples cluster by lineage (good) or by study (bad), and does
normalization move them from study-dominated to lineage-dominated?**

## Data — all local except the ~12 bigWigs

- **Manifest source (local, no download):** `chip-atlas/00.data/Histone/His.ALL.50.H3K27ac.AllCell.bed.gz`
  (3.1 GB). Every peak's name field embeds `ID=SRX…;Title=GSM…;Cell group=…;cell line=…;cell type=…;chip
  antibody=…`, so the full H3K27ac experiment manifest is recoverable locally (`01_build_manifest.py`).
- **The only download:** the ~12 selected per-experiment bigWigs from ChIP-Atlas (`03_download_bigwigs.sh`,
  dry-run by default). Nothing is fetched until the selection is approved.
- **Tooling:** conda env **`atac_hdac`** (has `multiBigwigSummary`, `computeMatrix`, `bedtools`, `pyBigWig`,
  numpy/pandas/scipy/matplotlib). PCA is done with numpy SVD — no sklearn needed.

## Pipeline

| step | script | reads | writes | network |
|---|---|---|---|---|
| 1 | `01_build_manifest.py` | local AllCell bed | `data/manifest.tsv` | none |
| 2 | `02_select_pilot.py` | `data/manifest.tsv` | `data/selection.tsv` | none (`--resolve-study` opt-in) |
| 3 | `03_download_bigwigs.sh` | `data/selection.tsv` | `data/bigwigs/*.bw` | **yes** (dry-run unless `--go`) |
| 4 | `04_quantify.py` | bigWigs + local bed | `results/signal.{npz,tab}` | none |
| 5 | `05_normalize_pca.py` | `results/signal.tab` | `results/*.png`, `results/pca_summary.tsv` | none |

## Run

```bash
conda activate atac_hdac
cd pilot/scripts
python 01_build_manifest.py                       # local scan → data/manifest.tsv
python 02_select_pilot.py                          # local → data/selection.tsv  (review this!)
bash   03_download_bigwigs.sh                      # DRY-RUN: prints + HEAD-checks URLs
bash   03_download_bigwigs.sh --go                 # actually download the ~12 bigWigs
python 04_quantify.py                              # multiBigwigSummary over H3K27ac peak regions
python 05_normalize_pca.py                         # raw vs quantile vs reference norm → PCA + batch/lineage split
```

## Decision criteria (what the pilot must show to proceed)

- **Green (proceed to Phase 1):** after normalization, top-PC structure and clustering track **lineage**, and
  the cross-study replicate of the same cell line co-clusters. Lineage silhouette > study silhouette.
- **Amber (fixable):** raw signal is study-dominated but a normalization (quantile / reference) flips it to
  lineage-dominated — record which normalization and carry it into Phase 2.
- **Red (rethink):** study batch dominates even after normalization; same-line-different-lab samples never
  co-cluster. Revisit approach (restrict to fewer uniform studies? require input-normalized bigWigs? ENCODE-only
  subset?) before committing to the full ~1,789-experiment pull.

## Where it runs

Steps 00–02 run **locally** (they need the 3.1 GB local bed and are trivial). Steps 03–05 can run locally
(minutes) **or** on **CSC Roihu** under `$SECACTS_CSC_PROJECT` — see `roihu/README.md` for the setup + a
`stage_and_submit.sh` driver. The pilot compute is light, so Roihu is mainly a *rehearsal* of the HPC flow
before the heavy Phase-1 pull; the real HPC payoff is Phase 1–4.

## RESULT (2026-07-16, run locally) — GREEN with quantile normalization

Ran end-to-end on the approved 13 experiments (4.0 GB bigWigs), quantified two ways — genome-wide 10 kb bins
(310,066) **and** H3K27ac peaks (50,000 merged from the AllCell catalog) — under three normalizations. Both
substrates agree:

| normalization | sil(lineage) ↑ | MCF7 cross-study cohesion ↓ (bins / peaks) | read |
|---|---|---|---|
| raw (log)     | −0.08 | 0.85 / 0.86 | batch/coverage-dominated |
| **quantile**  | **+0.13 / +0.11** | **0.26 / 0.29** | **biology emerges; replicates cluster** |
| reference (median-of-ratios) | −0.15 / −0.11 | 0.87 / 0.86 | distorted by an outlier (THP-1) |

*(cohesion = mean within-MCF7 PC-distance ÷ median over all pairs; < 1 means the same line from different
studies clusters tighter than average.)*

**Findings:**
1. **Cross-study batch is normalizable, not fatal — the key de-risk.** The three MCF7 experiments from three
   *different* studies (GSE172174 / GSE144404 / +1) collapse onto each other after **quantile** normalization
   (cohesion 0.26–0.29), while raw signal smears them along a coverage/depth PC1. The dominant technical risk
   named in `../DESIGN.md` is manageable.
2. **Among the three tested, quantile wins; median-of-ratios fails** — the latter is dominated by one
   high-coverage outlier (THP-1). **⚠️ superseded for Phase 2 (tooling survey, 2026-07-17):** *unconstrained*
   quantile is known to manufacture false positives (9.2% vs 80.5% motif match, S3norm benchmark) → **Phase 2
   adopts S3norm / S3V2-IDEAS** (purpose-built for epigenomic compendia, validated on H3K27ac), with quantile
   kept only as a baseline to benchmark against. See `../DESIGN.md` → Methods & tooling decisions.
3. **Positive but modest lineage silhouette (~0.11–0.13)** partly reflects *real subtype heterogeneity*, not
   residual batch: "Breast" splits into luminal (MCF7/T47D) vs basal/TNBC (MDA-MB-231); homogeneous colorectal
   (all adenocarcinoma) clusters cleanly. → **score at subtype resolution too**, as the design already plans.
4. **Peaks ≈ bins**, so the batch structure is genuine signal, not a background-bin artifact.

**Verdict: proceed (GREEN) with quantile normalization**, carrying subtype-resolution scoring forward. Outputs:
`results/{bins,peaks}/pca_batch_vs_lineage.png`, `results/*/pca_summary.tsv`.

## Status

- 2026-07-16 — Scaffold + full local run complete. `data/manifest.tsv` (11,538 H3K27ac experiments),
  `data/selection.tsv` (13), `results/` figures. Step B **done → GREEN**. Roihu scaffold (`roihu/`) ready but
  unused (pilot ran locally in minutes); it's staged for the Phase-1 scale pull.
