# SE-CaCTS dashboard

A static, no-backend results explorer for the SE-CaCTS lineage-specific super-enhancer atlas. Loads
precomputed TSV/JSON from `data/`; no build step. Modeled on the pyCaCTS dashboard.

## Tabs
- **Overview** — atlas size, and the panel-resolution caveat up front (calls at lineage/disease only).
- **SE atlas** — the most group-specific super-enhancers for any group, with JSD, permutation FDR, the
  copy number at the locus, the nearest gene, a UCSC out-link, and how many experiments called each SE.
  A **Level** toggle switches lineage / primary disease / subtype / cell line. At subtype and cell-line
  level it shows a **rankings-only** banner and does not report calls (see below). TSV download.
- **CN ablation** — the call-based comparison of scoring with vs without CN correction, under the
  permutation null: how many calls are *rescued* by correction vs the few *amplicon-driven* false calls it
  removes (every one a named recurrent lineage amplicon — MYCN, OTX2, FGFR2, ANO1).
- **Concordance** — the Phase-6 cross-layer check: genes next to a group-specific SE are specific to the
  same group far above background, and the concordance decays with SE→gene distance while a shuffled
  control stays flat.
- **SE finder** — look up a gene symbol and see every lineage / primary disease where a super-enhancer
  near it is group-specific, with its rank and whether the gene is also concordant in the expression layer.
- **About & methods** — the permutation FDR, and an explicit "what is NOT claimed" list.

## The one rule this dashboard enforces
It presents the **permutation** basis (`phase2/scores/atlas.s3.perm.*`) and nothing else. No analytic-null
FDR is staged, displayed, or exported — that null failed calibration outright (6% false calls on shuffled
labels), so its values are not in `data/` at all. Specificity **calls** are shown only at
**OncotreeLineage** and **OncotreePrimaryDisease**; **subtype and cell line are rankings only** (too few
lines per group for the permutation null to support a call). See the repo `RESULTS.md` and gotchas 71–72.

## Run locally
A server is required (`fetch()` of local TSVs is blocked over `file://`):
- `python3 -m http.server 8000 --bind 127.0.0.1` in this directory, then open http://127.0.0.1:8000
- or `npx serve -l 8000`

## Regenerate `data/`
```bash
~/miniconda3/envs/atac_hdac/bin/python phase2/scripts/60_stage_dashboard.py   # run from the repo root
```
It reads the permutation scores + CN ablation + concordance summaries and the union catalog, annotates
every specific SE with its nearest gene and coordinates, and writes the 9 files in `data/`. Rerun it
whenever the permutation scoring is rerun.

## Deploy (GitHub Pages)
This directory **is** the Pages root: repository Settings → Pages → Source = "Deploy from a branch",
branch `main`, folder `/docs`. The site is self-contained and static — no secrets, no build step, and no
external requests except the UCSC links the user clicks.
