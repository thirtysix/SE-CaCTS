# Phase-2 reference SE atlas — v2 (copy-number-expanded panel)

Built on Roihu 2026-09-28. v1 (`../results/`) plus **552 new experiments** from the **122 lines** that
became copy-number-correctable through CMP WES pureCN 2025 (106 lines) and DepMap 26Q1 MC_WES (16 lines).
The set is built by `phase1/scripts/15_expansion_set.py`; the manifest is `phase2/data/pull_srx.v2.txt`
(3,468 experiments = v1's 2,916 + 552).

How it was made, and what was held fixed so v1 and v2 differ only by the added samples:

- **Same caller code.** The new samples were pulled by `array.slurm` against the staged `cnrose` that
  built v1 (unchanged), with a 1 kb genome archive instead of 100 bp (the retained v1 archives were
  coarsened to 1 kb, so the retention set is uniform).
- **Same quantification grid.** `grid.20.bed` was not rebuilt. It holds 96.4% (bp-weighted) of the new
  samples' own bed20 peak territory; 124 union loci fall entirely outside it and carry zero signal.
- **Same normalization reference.** `reduce.slurm` with `S3_REF=SRX16495452` (v1's medoid), so v1's
  2,136 samples normalize bit-identically (max |ΔA| = max |ΔB| = 0).

## Inventory

| file | shape | tracked in git? |
|---|---|---|
| `atlas.se_signal.tsv.gz` | 48,180 SE loci × 3,468 samples (agnostic) | **no** — >100 MB |
| `atlas.se_presence.tsv.gz` / `atlas.union_catalog.bed.gz` | same | yes |
| `atlas.s3.se_signal.tsv.gz` | **47,101 SE loci × 2,563 samples** (S3norm, `--min-peaks 2000`) | **no** — >100 MB |
| `atlas.s3.se_presence.tsv.gz` / `atlas.s3.union_catalog.bed.gz` | same | yes |
| `atlas.s3.s3norm_params.tsv.gz` | per-sample A, B | yes |
| `failed_srx.txt` | `SRX20868733` (retired, 404) | yes |

The QC gate kept 427 of the 552 new samples and 104 of the 122 new lines, so the scored panel is
**2,563 samples → 386 lines** (282 DepMap WGS + 91 CMP WES + 13 DepMap WES). Grid→SE reconstruction is
exact (`max|err| = 0`).

The two signal matrices live on Dropbox (this directory) and on Roihu scratch; they are not in git.
