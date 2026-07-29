# Phase-2 reference SE atlas

Built on Roihu 2026-07-21 (reduce job 303509, 7:21 wall) from the full 2,917-sample ChIP-Atlas pull,
and fetched off scratch 2026-07-22 with `phase2/scripts/50_fetch_and_clean.sh`.

**This is not regenerable without a ~43 BU re-pull.** Treat it as primary data, not as build output.

## Inventory

| file | size | shape | tracked in git? |
|---|---|---|---|
| `atlas.se_signal.tsv.gz` | 423 MB | 43,931 SE loci × 2,916 samples | **no** — >100 MB, Dropbox + SSD only |
| `atlas.se_presence.tsv.gz` | 1.7 MB | same shape, 0/1 | yes |
| `atlas.union_catalog.bed.gz` | 0.6 MB | 43,931 loci | yes |
| `atlas.s3.se_signal.tsv.gz` | 306 MB | 42,943 SE loci × 2,136 samples | **no** — >100 MB, Dropbox + SSD only |
| `atlas.s3.se_presence.tsv.gz` | 1.5 MB | same shape, 0/1 | yes |
| `atlas.s3.union_catalog.bed.gz` | 0.6 MB | 42,943 loci | yes |
| `atlas.s3.s3norm_params.tsv.gz` | 0.1 MB | per-sample fitted A, B | yes |
| `failed_srx.txt` | 12 B | `SRX20868733` | yes |

`atlas.*` = **agnostic** (no normalization beyond the caller). `atlas.s3.*` = **S3norm**, which applies
the `--min-peaks 2000` QC gate — that gate is why it has 2,136 samples rather than 2,916 (784 weak
samples dropped, 27%, consistent with the pilot's ~31%). The gate is not optional: ungated S3norm is
*worse* than quantile normalization.

Both atlases reconstruct grid→SE with `max|err| = 0`.

The **CN-corrected calling-time atlases (`atlas.cn`, `atlas.s3.cn`) are absent** — the pull ran
cnrose agnostic-only. They are recoverable without a re-pull from the retained `.enhancers.tsv` +
DepMap CN (gotcha 59). Note that *scoring-time* CN correction — the one that matters for specificity —
is applied downstream on the agnostic matrix and does not need them.

`failed_srx.txt` records the single drop from the 2,917-sample manifest: **SRX20868733** (HCT 116) is a
retired experiment that permanently 404s. The reduce was re-run with `FORCE=1` to build from the 2,916
available.

## Copies, and where the bulk lives

| location | contents | notes |
|---|---|---|
| `phase2/results/` (this dir, Dropbox) | all of `results/` (734 MB) | synced; git tracks all but the two signal matrices |
| an external SSD | `results/` **+ `out/`** (26 GB) | the full per-sample retention set — 2,916 × 9 files |
| Roihu `/scratch/$SECACTS_CSC_PROJECT/se-cacts/phase2/` | `results/` (+ `out/` until cleaned) | scratch is un-backed-up and billed 6 BU/TiB-h |

The **per-sample retention set is only on the 2 TB SSD** — it is not in Dropbox (26 GB) and not in git.
It holds `.enhancers.tsv` (the calling-time signal, which *cannot* be reconstructed from the fine
matrix — gotcha 30), the exact `grid.20`/`grid.10` float32 columns, the 1 kb genome-wide archives that
Phase-3 CN inference depends on, plus `qc.json` and provenance per sample.

## Verification performed at fetch (2026-07-22)

- `.done` markers: remote 2,916 = local 2,916
- per-extension file counts in `out/`: 2,916 for each of the 9 retained file types
- 20 `grid.20.f32` columns md5-identical to Roihu (script's own spot-check)
- all 8 files in this directory md5-identical to Roihu
- `gzip -t` passes on all 7 `.gz` here and on 40 sampled `.bin1000.f16.gz` archives
- atlas dimensions re-derived from the files and match the build log
