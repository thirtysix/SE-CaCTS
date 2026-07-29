# Roihu stream/batch pull — design

**Status:** design only (2026-07-18). Nothing launched. Reviewed decisions marked **[decide]**.
Covers the Phase-1→Phase-2 boundary: staging the QC-pass bigWigs and turning them into a durable
signal matrix **without ever holding ~0.9 TiB on disk at once**.

Read first: `ROADMAP.md` (Phase 1–2),
`pilot/roihu/README.md` (the scaffold this extends), `pilot/scripts/03_download_bigwigs.sh` +
`04_quantify.py` (per-SRX logic reused verbatim).

---

## 1. The problem in one paragraph

The analysis set is **2,917 QC-pass, CN-correctable H3K27ac experiments on 323 cell lines**
(`phase1/data/phase1_manifest.tsv`, `qc_pass==1 && has_cn==1`). At a measured mean of **322 MB/bigWig**
that is **~0.9 TiB** — but Roihu's default scratch is **250 GiB**. The naive "download everything, then
compute" plan is ~4× over budget. The durable *output* we actually want (a region × sample signal
matrix + per-sample SE calls) is only a **few GiB**. So the bigWigs are **transient**: the design streams
them through a bounded on-disk buffer, extracts the durable signal, and deletes each bigWig before the
buffer fills. Steady-state stays well under 250 GiB and **no scratch-quota increase is needed**.

## 2. Constraints (Roihu-specific, from recon)

- **Compute nodes have DIRECT outbound internet — VERIFIED 2026-07-18.** A `srun` HEAD to the live
  ChIP-Atlas bigWig endpoint returned **HTTP 200 in ~1.0 s, no proxy**, on both the `test` pool (`rc4184`)
  and the production `small`/`longrun` pool (`rc6136`). DNS resolves on compute nodes; no `http_proxy` set
  or needed. **This is a Roihu departure from the usual CSC login-node-only model** (Puhti/Mahti, and the
  `csc-puhti-job` skill, assume compute nodes are offline). Consequence: downloads and compute do **not**
  need decoupling — each array task fetches its own bigWig (§5). *(Re-verify if CSC changes network policy.)*
- **Scratch 250 GiB / projappl 15 GiB** default. The whole point of this design is to live inside that.
- **SSH certificate expires every 24 h** (re-sign to *reconnect*; **submitted jobs survive** expiry — so a
  long feeder loop should run inside `tmux`/`screen` and be restart-safe, not depend on a live laptop).
- **deepTools is not a module → tykky** container (`pilot/roihu/00_build_env.sh`, already written).
- Partitions: `small` (3-day default), **`longrun` (10-day)** for the feeder; array jobs on `small`.

## 3. Key architectural choice — one download per bigWig via a *fixed grid*  **[decide]**

The obstacle to a single streaming pass is a dependency inversion: the **union SE catalog** (Phase 2) is
*derived from the data*, so the regions to quantify aren't known until every sample's SEs are called — which
normally forces a **second** pass over all bigWigs (→ ~1.8 TiB of transfer). We break that with a **fixed
quantification grid** fixed *before* the pull:

> **Grid G** = the merged union of all ChIP-Atlas H3K27ac peaks, built locally **now** from
> `…/chip-atlas/00.data/Histone/His.ALL.50.H3K27ac.AllCell.bed.gz` (already on disk, no download).
> ROSE/LILY stitched SEs are unions of peaks, so **every union-SE region is a sum of grid rows** → the
> region × sample matrix is recovered by *aggregating the fine grid*, with **no second bigWig pass**.

**Single-pass loop, per resident bigWig** (both steps while it's on disk, then delete):
1. **Call SEs** (LILY/ROSE-stitch on the bigWig) → small per-sample SE BED. *(feeds the union catalog)*
2. **Quantify over grid G** (`multiBigwigSummary BED-file`) → one **column** of the fine matrix (mean
   signal per grid row; row lengths are known from G, so SE total = Σ meanᵢ·lenᵢ, matching ROSE's density).

**After the pull (no bigWigs):** merge the 2,917 SE BEDs → union catalog (cSEAdb recipe: overlap ≥25% of
width, constituent resolution); then aggregate the fine matrix rows per union SE → the **region × sample
matrix**. Normalize (S3norm/S3V2) and CN-correct (DepMap `OmicsCNGeneWGS`, local at
`…/DepMap/2026q1/`) downstream — those steps never touch a bigWig.

- **Recommended:** fixed-grid single pass (≈0.9 TiB transfer, ~6–37 GiB durable fine matrix depending on
  grid granularity, no re-download).
- **Fallback (simpler, 2× transfer):** two streaming passes — pass 1 SE-call, barrier to build the union,
  pass 2 quantify over the union. Choose this only if building/aggregating grid G proves fiddly.
- **Not recommended:** remote HTTP range-reads (`pyBigWig` over URLs, no download) — technically possible
  but hammers ChIP-Atlas with millions of small requests from the login node; impolite and fragile.

*Trade-off to weigh:* the fixed grid caps SE quantification to peak-covered genome (fine for H3K27ac SEs).
Genome-wide 1 kb bins are the fully-generic alternative (larger matrix, no peak dependency).

## 4. Budgets — MEASURED on the pilot (2026-07-20), not estimated

Per-sample cost, measured end-to-end through `scripts/40_pull_one.sh` against the live endpoint:

| Resource | Measured | At 2,917 samples |
|---|---|---|
| cnrose compute (1 CPU) | **15 s** (155 MB bw) → **42 s** (569 MB bw); ~28 s at the 322 MB mean | **~23 core-h** |
| download | ~30 s at 10 MB/s (322 MB mean) | ~0.94 TiB transfer |
| **wall per sample** | **~60 s** | **~50 core-h** total |
| peak RSS | **1.4 GB** (largest pilot bigWig) → `--mem=2G`, `--cpus-per-task=1` | — |
| durable output | **20 MB/sample** (9 files, §8.5.4) | **~60 GB** |
| transient bigWig | 322 MB × concurrency | ~5 GiB at `%16` |

Wall-clock ≈ `2917 × 60 s / concurrency`: **~3.0 h at `%16`**, ~1.5 h at `%32`. At `%16` the aggregate
download rate (~160 MB/s) and the compute are roughly balanced, so raising concurrency past ~32 mostly just
leans harder on ChIP-Atlas (§6) for little gain.

Steady-state scratch = transient + durable ≈ **65 GB** — comfortably inside the 250 GiB quota.

### 4.1 Actual BU cost (rates from docs.csc.fi/computing/hpc-billing, checked 2026-07-21)

Roihu core-based partitions (`small`, `longrun`, `interactive`, `test`):

```
Total CPU BU = max( 0.75 BU/coreh × cores , 0.375 BU/GiBh × mem ) × runtime-hours
             + 0.02 BU/GiBh × reservedstorage × runtime-hours
```

**`--cpus-per-task=1 --mem=2G` is exactly optimal.** The two terms of the `max` are equal at precisely
**2 GiB per core** (0.75 × 1 = 0.375 × 2), so 2 GiB is free headroom and anything above it is a pure penalty:

| request | BU/hour | vs 1 core |
|---|---|---|
| 1 core + **2 GiB** | **0.75** | **1.0×** |
| 1 core + 4 GiB | 1.50 | 2.0× |
| 1 core + 8 GiB | 3.00 | 4.0× |

Asking for *less* than 2 GiB saves nothing (the core term dominates). This is what the chunked-interval work
bought: peak RSS 2.74 GB → 1.4 GB moved the job from the 4 GiB tier into the free 2 GiB headroom, halving the
bill. We reserve no `--gres` local disk, so the `reservedstorage` term is 0.

**Compute for the whole pull: 48.6 core-hours × 0.75 = ~36 BU.** Negligible.

### 4.2 ⚠️ Storage, not compute, is the real cost

Roihu bills scratch at **6 BU/TiB-hour from the first byte** — unlike Puhti (1 TiB free) and Mahti, Roihu has
**no free tier**. The retention set is ~58 GB:

| what | size | BU/day | BU/month |
|---|---|---|---|
| full retention set | 58.3 GB | **7.6** | 229 |
| ├ genome-wide archive | 43.8 GB | 5.7 | 172 |
| └ exact tier only | 14.6 GB | 1.9 | 57 |

**Leaving the output on scratch overtakes the entire cost of the pull after ~4.8 days**, and the archive is
~75% of that bleed. So `scripts/50_fetch_and_clean.sh` (fetch → verify → delete from scratch) is part of the
pipeline, not an afterthought — and scratch is un-backed-up and subject to cleaning, so it was never a home
for the archive anyway. `ARCHIVE_ONLY=1` moves just the bulky part and leaves the exact tier in place for
downstream work.

Still worth a **2-task test array + `seff`** before the full launch: it confirms the real download rate from
a compute node (the one number not measurable from the laptop) and that nothing about the memory shape
surprises us.

## 5. Orchestration — a plain SLURM array (compute-node internet ⇒ no feeder)  **[simplified]**

Because compute nodes have direct internet (§2, verified), the pull collapses to **one self-contained
SLURM array job** — no login-node feeder, no producer/consumer, no backpressure, no 0.9 TiB funnelled
through a shared login node:

```
sbatch --array=0-2916%N array.slurm      # N = max concurrent tasks (the throttle; see §6)
# each task i, on a compute node, fully self-contained:
#   SRX = pull_set.tsv[i]
#   skip if SRX.done exists (resumable)
#   wget -c the bigWig (primary→legacy fallback) into node-local $LOCAL_SCRATCH or scratch
#   run SE caller + multiBigwigSummary over grid.bed  →  SRX.se.bed, SRX.signal.tsv
#   touch SRX.done ; delete the bigWig
#   on hard download failure: append SRX to failed.txt and exit 0 (don't fail the array)
```

- **Steady-state disk** = (concurrent tasks) × 322 MB. At `%64` that's ~20 GiB — trivially inside 250 GiB.
  Prefer **node-local scratch** (`$LOCAL_SCRATCH`) for the transient bigWig so it never touches shared scratch.
- **Resumable / restart-safe** by construction: re-submitting skips `*.done`; a preempted task just re-runs.
  Nothing depends on a live laptop or the SSH cert (submitted jobs survive cert expiry).
- **Failure handling:** primary+legacy endpoint fallback already lives in `03_download_bigwigs.sh`; append
  hard failures to `failed.txt` and **continue**, then **report the count explicitly** at the end (no silent
  truncation). Expect a handful of 404s (retired SRX) — they drop from the analysis set with a note.
- **Fallback (only if a chosen partition turns out to be offline):** the earlier feeder / producer-consumer
  design (git history of this file) still applies — login-node download into scratch + no-internet array
  consumers. Not needed given the §2 result, but kept in mind.

## 5.1 Write-safety under parallelism

Up to N tasks run at once; correctness rests on **one writer per file, sharded by SRX** — no task writes a
file another task writes. There is **no shared mutable "grid"**: `grid.bed` is built once (locally) and
shipped as immutable, read-only input; tasks only *read* it.

- **Per-sample outputs are SRX-keyed** — task `i` writes only `SRX.se.bed`, `SRX.signal.tsv`, `SRX.done`
  (unique names → disjoint writes, no locks, no contention).
- **Shared products are built post-barrier by ONE process** — the region×sample matrix (`aggregate.py`) and
  the union SE catalog are assembled *after* the array finishes (`sbatch --dependency=afterok:<arrayid>`),
  reading all `SRX.*` files and writing the single output once. We **shard-then-reduce**, never
  incrementally append to one shared file.
- **Crash/preempt safety** — write each output to `SRX.signal.tsv.tmp.$SLURM_JOB_ID`, then `mv` to the final
  name (rename within one FS is **atomic** → readers never see a half-file); write `SRX.done` **last**.
  Resume keys on `SRX.done`, not mere file existence, so a task killed mid-write leaves no marker and is
  redone cleanly.
- **No shared appends** — download failures write a per-task `SRX.failed` marker, **not** `echo >> failed.txt`
  (`O_APPEND` is not atomic across nodes on Lustre). The post-step globs `*.failed`.
- **Lustre metadata (perf, not correctness)** — shard the ~2,917 output files into prefix subdirs
  (`out/<SRX-prefix>/`) so no single directory takes thousands of concurrent creates.

## 6. The remaining external courtesy — don't hammer ChIP-Atlas

The login-node bottleneck is gone; the real constraint is now **being a good citizen to the ChIP-Atlas host**
(`chip-atlas.dbcls.jp`, 133.39.78.111) — hundreds of parallel compute nodes could otherwise DoS it:
- **Throttle the array concurrency** (`--array=…%N`, start modest e.g. N≈16–32) — this doubles as the disk
  cap. Add `wget --limit-rate` if needed.
- `wget -c` (resumable) so preemptions/retries don't refetch whole files.
- Fetch each bigWig **once** — the fixed-grid single pass (§3) already guarantees this (no second pass).
- **Consider the AWS Open-Data mirror** of ChIP-Atlas (DBCLS hosts bigWigs on S3) — `aws s3 cp`/`s3cmd`
  spreads load off the dbcls host and may be faster; verify the per-SRX S3 key layout before relying on it.
- Still **scope compute-BU cost + wall-clock and get a go/no-go before the full 2,917-task launch** (per
  approval for large agent runs) — the token cost is only code-gen/monitoring, but the BU spend and
  day-scale runtime warrant an explicit green light.

## 7. What can be done locally *now* (no Roihu, no download, cheap)

1. **Build grid G** from the local `His.ALL.50.H3K27ac.AllCell.bed.gz` (zcat → filter hg38 → `bedtools
   merge`) → `phase2/data/grid.bed` + a companion lengths file. This is the fixed reference the whole pull
   hinges on; producing it first de-risks §3.
2. **Emit the pull selection**: `qc_pass==1 && has_cn==1` slice of the manifest → `phase2/data/pull_set.tsv`
   (2,917 SRX + cell line + CVCL + ModelID), the feeder's work-list.
3. **Pin the SE caller**: decide LILY (bigWig-native, ROSE+CN in one) vs ROSE2, add to `env.yaml`, and dry-
   run it on the 13 local pilot bigWigs to confirm the per-sample step + its outputs before scaling.
4. **Rehearse end-to-end on the 13 pilot bigWigs** (already local, 4.2 GB): grid-quantify → aggregate →
   matrix, so the aggregation math (Σ meanᵢ·lenᵢ per SE) is verified on a known-good tiny set.

Items 1–2 are pure local file ops; 3–4 reuse the pilot. None need approval.

## 8. Open decisions to confirm before building

- **[decide]** Fixed-grid single pass (§3 recommended) vs two-pass fallback.
- **[decide]** Grid = merged H3K27ac peaks (compact, biological) vs genome-wide 1 kb bins (generic, larger).
- **[decide] Grid peak-threshold — resolved empirically 2026-07-18 (recommend `bed20`).** Built three grids
  from **per-SRX ChIP-Atlas peak BEDs** (see the better path below) and compared:

  | Grid | Q-value | SRX coverage | Regions | Genome % | Download |
  |---|---|---|---|---|---|
  | `grid.bed` (stringent) | <1E-50 | 2,800/2,917 ⚠️ | 304k | 7.6% | AllCell (local) |
  | **`grid.20.bed`** ★ | <1E-20 | **2,917/2,917** | 505k | 11.0% | 2.38 GB |
  | `grid.10.bed` | <1E-10 | 2,917/2,917 | 774k | 20.0% | 4.34 GB |

  **Recommend `bed20`**: closes the 117-gap, captures moderate SE constituents, stays enhancer-focused
  (11% genome) where `bed10`'s 20% pulls in background that dilutes specificity. Governing rule unchanged:
  **grid threshold = SE-caller constituent threshold** — revisit when the caller is pinned.

- **Better grid-build path than the AllCell download (verified 2026-07-18).** ChIP-Atlas serves **per-SRX
  peak BEDs** at `…/data/hg38/eachData/bed{NN}/<SRX>.{NN}.bed` (NN∈{05,10,20,50}; primary
  `chip-atlas.dbcls.jp`, legacy `dbarchive.biosciencedbc.jp`). Fetching only our 2,917 SRX at a threshold is
  **~2–4 GB, ~7 min, directly on a Roihu node** — vs the 107 GB uncompressed AllCell download (⅔ of which is
  cell types we filter out). `phase2/scripts/21_grid_from_persrx.sh` does this (curl + awk merge, no deps,
  resumable). **This supersedes `20_build_grid.sh`/the AllCell route** for grid construction. NB the two
  ChIP-Atlas threshold numberings: bulk-filename `NN` (05/10/20/50 = Q<1E-0N) vs the *site* picker
  (50/100/200/500 = −10·log₁₀Q); they differ by 10× (filename `50` = site `500` = Q<1E-50).
- ~~Orchestration feeder vs producer/consumer~~ **resolved:** plain SLURM array (§5), since compute nodes
  have internet (§2). Only the array concurrency `%N` remains to tune (§6).
- ~~SE caller LILY vs ROSE2~~ **RESOLVED 2026-07-18 → ROSE (via `ROSE2_callSuper` on bigWig signal).**
  Both tools are BAM-first (not bigWig-native), but ROSE2's ranking core is decoupled: `ROSE2_callSuper.R`
  consumes a pre-quantified signal table (`chrom/start/stop/id/length/signal/opt-WCE`) + `calculate_cutoff()`
  works on any numeric vector → feed it **bigWig**-quantified signal. **LILY is out**: needs BAMs (H3K27ac
  *and* Input) via HMCan, CN correction is **non-optional** and input-BAM-dependent (compendium inputs often
  missing) — fails our bigWig-only, separable-CN, missing-input reality. **Per-sample step:** stitch our
  per-SRX `bed20` peaks (=grid source → coverage guaranteed) within 12.5 kb → `multiBigwigSummary` per
  stitched region → `ROSE2_callSuper` (controlless; WCE optional). CN-agnostic; DepMap WGS CN applied at
  scoring keeps the corrected-vs-uncorrected co-primary intact. Paper framing: we decouple LILY's CN-aware
  idea into a scalable bigWig-only + DepMap-CN form (cite LILY/ROSE as lineage).
- **[decide]** Stream-only vs request a ~2 TiB scratch quota anyway (a quota would simplify to a plain
  download-then-compute, at the cost of a BU request + a fatter Dropbox-free footprint). Stream-only is the
  default here.
- **Open (carried):** how to assign a single CN ratio to an SE spanning heterogeneous CN segments
  (`ROADMAP.md` open questions) — a Phase-3 concern, but the grid rows carry per-region coordinates that
  make a per-segment CN join tractable later.

## 8.5 Retention policy — what each task MUST emit before the bigWig is deleted

The pull is destructive by design (download → quantify → **delete**), and re-running it costs ~0.9 TiB of
transfer plus the courtesy budget of §6. So the governing question for every per-sample artifact is not "is
this useful?" but **"could we reconstruct it later without the bigWig?"** Anything answered *no* is
effectively write-once: if it is not emitted during the pull, it is gone until a full re-pull.

### 8.5.1 What is genuinely irreversible

Three things — stated first as they stand **if we keep only a peak-derived grid**, since that was the
starting design. §8.5.2 then retires the first two by adding a genome-wide archive.

1. **Signal at any locus or resolution outside the emitted grid.** At `bed20` the grid is 11% of the genome;
   the other 89%, and any sub-row detail, is gone. *(Retired by the §8.5.2 archive, to ~1% median.)*
2. **The SE calling-time signal.** `cnrose` quantifies *stitched regions* directly from the bigWig
   (`pipeline.py:67`), and a stitched region spans the **gaps between** constituent peaks, which are not grid
   rows. Measured on SRX067407: the `.enhancers.tsv` `SIGNAL` column and the sum of overlapping grid rows agree
   for only **0.3%** of regions (median relative difference 109%, p90 376%). These are two different
   quantities serving two purposes — calling vs. atlas quantification — so **the calling signal cannot be
   reconstructed from the fine matrix**, and without it you can never re-call SEs at a different cutoff or
   audit why a region was or wasn't a super-enhancer. *(Retired by the §8.5.2 archive: stitched regions are
   contiguous spans, which reconstruct at ~0.05%. We keep `.enhancers.tsv` anyway — 350 KB for exactness.)*
3. **Anything needed to decide sample QC.** We learned the hard way (DESIGN.md §normalization) that S3norm is
   only safe behind a QC gate. Peak counts come from the ChIP-Atlas BEDs, but every *signal*-based QC metric
   (dynamic range, FRiP, background level) requires the bigWig. Emitting them costs ~1 KB/sample; omitting them
   means a future QC criterion cannot be applied retroactively. **Not retired by anything — just emit them.**

### 8.5.2 Grid choice is irreversible — but a binned archive makes it recoverable

`bed20` was chosen (§8) for good reasons, and §8's own rule says "grid threshold = SE-caller constituent
threshold — **revisit when the caller is pinned**". The caller is now pinned (`cnrose`), and revisiting after
the pull would mean re-downloading everything. ChIP-Atlas serves **four** thresholds, not three:

| grid | Q-value | regions | genome % |
|---|---|---|---|
| `bed05` | <1E-05 | *not built* | ~30% (est.) |
| `grid.10.bed` | <1E-10 | 774k | 20.0% |
| **`grid.20.bed`** ★ | <1E-20 | 505k | 11.0% |
| `grid.bed` | <1E-50 | 304k | 7.6% |

(Higher number = **stricter** = fewer regions. `bed05` is the loosest.) Note the regress: each looser grid is
"safer" insurance, and taken to its limit the insurance grid becomes the whole genome — at which point **bin
size, not peak threshold, is the only parameter**. That is the right frame, so it was measured directly.

**Measured** (`analysis/binned_storage_eval.py`, **all 13 pilot bigWigs**, exact float16 bins on the global
lattice; ranges are across samples):

| bin | SE-sum median | SE-sum p90 | SE-sum **worst case** | signed bias | f16+gz/sample | ×2917 |
|---|---|---|---|---|---|---|
| 100 bp | **0.29–0.98%** | 1.1–4.2% | **4.4–20%** | +1.6 to +6.4% | ~15 MB | ~0.04 TB |
| 200 bp | 0.84–2.17% | 2.7–9.6% | 12–65% | +3.6 to +18.6% | ~10 MB | ~0.03 TB |
| 500 bp | 2.50–5.01% | 6.9–19% | 22–124% | +7.5 to +52% | ~5 MB | ~0.015 TB |
| 1 kb | 3.97–9.40% | 10–32% | 30–268% | +12 to +70% | ~3 MB | ~0.01 TB |

Two things decide it. First, **single grid rows reconstruct badly at every bin size** — a ~685 bp row is only
~7 bins wide and H3K27ac varies sharply exactly where peak boundaries fall (100 bp: 2.8–3.2% median but
14–21% p90 on individual rows). Second, and decisively, the per-row error is **biased, not random**: the
signed mean ≈ the absolute mean and is **positive in 13/13 samples**, so errors **accumulate** across the ~12
rows of an SE instead of cancelling. (Mechanism: the grid is the union across all 2,917 samples, so most rows
are not a given sample's own peaks, and neighbouring bins bleed signal in.) Even at 100 bp the worst-case SE
lands 20% off — fine for deciding whether a re-pull is warranted, not fine as the atlas.

**Decision: bins are an ARCHIVE, not a replacement.**

- **Primary quantification stays exact**, from the bigWig, during the pull. `cnrose` is validated bit-for-bit
  vs ROSE2; degrading the atlas to 1–4% error would forfeit that for no good reason, when exact costs 2 MB.
- **Also emit `grid.10` exactly** (+3.1 MB/sample, 9 GB): trivial while the bigWig is open, and it makes one
  specific future — moving to the looser constituent threshold — exact rather than approximate.
- **Also emit a genome-wide 100 bp float16 archive** (~22 MB/sample compressed, ~65 GB): this is what
  actually retires the irreversibility. Any future region set — a `bed05` grid, a different stitch window,
  a locus nobody anticipated — is recoverable to ~1% median without re-downloading. Contiguous spans
  (e.g. stitched regions) come back at 0.05%.

100 bp over 200 bp is chosen deliberately: it is the only tier whose error is small enough to *substitute*
for exact in a pinch rather than merely to *evaluate*, and the ~10 GB difference is noise against 0.94 TB of
transfer and the cost of a re-pull. Storage is the cheapest resource in this pipeline.

**Cost: none.** The naive implementation (quantify each grid, then bin the genome, then compute QC) re-reads
the same coverage runs three times and took 67 s/sample — more than the ~30 s download, which would have
broken the "free" claim above. `cnrose.io.scan_bigwig` unpacks each chromosome's runs **once** and evaluates
every grid and every bin edge against it, so the whole retention set costs **17 s/sample — less than the
16 s the single-grid quantification cost before**. Verified bit-exact against the original `quantify`
(`cnrose/tests/test_retention.py`), and the ROSE2 contract still passes 13/13 at Jaccard 1.0000.

### 8.5.3 Emit a fixed off-grid background sample

Superseded in part by §8.5.2: a genome-wide archive contains true off-peak background by construction, so a
separate background bin set is no longer needed for FRiP or the S3V2 variance step. **If the genome-wide
archive is dropped for storage reasons, reinstate this**: ~50k fixed 1 kb bins disjoint from the `bed10`
union, shared across samples, 200 KB/sample. Every current "background" estimate (S3norm's `--lo-band`) is
computed over *weak peaks*, not true background — a documented deviation in `s3norm.py` that only off-grid
data can fix.

### 8.5.4 The policy

| Tier | Artifact | Size/sample | Total (2,917) | Why |
|---|---|---|---|---|
| **KEEP — exact, primary** | `<SRX>.grid20.f32` | 2.0 MB | 5.9 GB | atlas quantification (bit-for-bit) |
| | `<SRX>.grid10.f32` | 3.1 MB | 9.0 GB | makes the looser-threshold future exact (§8.5.2) |
| | `<SRX>.enhancers.tsv` | ~350 KB | 1.0 GB | calling-time signal (§8.5.1); ~0.05% from the archive, but exact here |
| | `<SRX>.qc.json` | ~1 KB | 3 MB | re-QC without re-download (§8.5.1) |
| **KEEP — approximate archive** | `<SRX>.bin100.f16.gz` (+`.json`) | ~15 MB | **~44 GB** | retires the irreversible grid choice (§8.5.2) |
| **KEEP — cheap audit trail** | `<SRX>.se.bed`, `<SRX>.cn.se.bed` | ~60 KB | 175 MB | pins the exact call set + cutoff |
| | `<SRX>.provenance.json` | ~1 KB | 3 MB | source URL, bytes, checksum, ChIP-Atlas version, cnrose git SHA |
| **DISCARD — regenerable** | the bigWig | ~322 MB | (0.9 TiB) | the entire point |
| | union catalog, `se_signal`, `se_presence` | — | — | `aggregate.py` re-derives in minutes |
| | S3norm-normalized matrices | — | — | deterministic from the fine matrix |

**Total durable ≈ 60 GB** (measured on a real sample, not estimated), of which ~44 GB is the genome-wide
archive. That fits the 250 GiB scratch comfortably, but it is too large for Dropbox — park the archive in
**CSC Allas** (object storage) or on a local disk, and keep only the ~16 GB exact tier alongside the code.
Dropping the archive returns you to ~16 GB at the cost of reinstating §8.5.3 and re-accepting the
irreversible grid choice.

`<SRX>.qc.json` carries: total signal; percentiles p1/5/25/50/75/90/95/99/99.9; max; zero fraction;
**dynamic range p99/median**; FRiP proxy; the fitted ROSE **tangent cutoff**; and peak/region/SE counts.

⚠️ **The dynamic-range gate statistic must be computed over the GRID column, not over genome-wide bins.**
It is the number S3norm's fitted exponent tracks (Spearman −0.956) and therefore the number the §8.5 QC gate
turns on. Genome bins are ~63% zero, so their median is 0 and the ratio is undefined — an early
implementation silently emitted `null` for every sample. `qc_stats` now reports grid-based and bin-based
statistics in separate blocks so they cannot be confused; the bin-based variant uses a *nonzero* median.

### 8.5.5 Format: float32 binary, not TSV

The pilot's `<SRX>.signal.tsv` is **15 MB**; the same 504,855 values as float32 are **2.0 MB** — 7.5×. Across
2,917 samples at both grids that is 110 GB of text vs 14.9 GB of binary. Text also loses exactness (the pilot
writes `%.6g`). Emit `.f32` raw little-endian arrays in grid-file row order, with row order pinned by the
committed grid BED and a length check on read; keep TSV only for the small per-sample tables a human reads.

### 8.5.6 What would still force a re-pull

With the genome-wide archive in place, only two things:

1. **A different genome build.** Nothing short of the original reads helps here.
2. **Any analysis needing better than ~1% accuracy at a region set we did not quantify exactly.** In
   practice this means: if a future decision hinges on exact signal at `bed05` rows, the archive tells you
   *whether it matters* (~1% median, 3–4% p90) but not the exact value.

Notably **resolved** by the archive: genome-wide input-inferred CN (HMCan/Control-FREEC, Phase 3) needs broad
coverage, which a peak-only grid could never provide — this was the main residual risk under the previous
draft of this section, and it matters because ~94% of the expansion lines lack DepMap CN. CN segments are
megabase-scale, so 100 bp bins are far finer than required.

## 8.6 Staged rollout — the runbook

A 2,917-task array is not something to launch cold. Stages are **prefixes of the same manifest sharing one
output tree**, so each stage extends the last: `.done` markers make later stages skip finished samples, array
indices stay stable (a prefix indexes identically to the whole), and nothing is ever recomputed.

| stage | samples | `%conc` | ~wall | ~BU | what it is for |
|---|---|---|---|---|---|
| `smoke` | 2 | 1 | ~2 min | 0.03 | mechanics; the real compute-node **download rate**; `seff` memory shape |
| `small` | 32 | 4 | ~8 min | 0.4 | concurrency, endpoint throttling, Lustre behaviour |
| `medium` | 256 | 8 | ~32 min | 3.2 | straggler and failure *rate* at scale; QC spread across lineages |
| `full` | 2,917 | 16 | ~3 h | 36 | the real thing |

```bash
bash scripts/42_preflight.sh                 # read-only; costs nothing; refuses to pass on any blocker
bash scripts/43_run_stage.sh smoke --go
bash scripts/43_run_stage.sh report          # progress, failures, QC spread, disk + BU/day
ssh roihu 'seff <jobid>'                     # confirm 1 core / <=2 GiB — the billing shape (§4.1)
bash scripts/43_run_stage.sh small  --go
bash scripts/43_run_stage.sh medium --go
bash scripts/43_run_stage.sh full   --go
ssh roihu 'cd $WORK && sbatch --dependency=afterok:<aid> reduce.slurm'
bash scripts/50_fetch_and_clean.sh <local-dest> --clean     # scratch bills 7.6 BU/day (§4.2)
```

**Check between stages, not just at the end** — each stage answers a question the previous one could not:

- **after `smoke`** — `seff` must show `MaxRSS` ≲ 2 GiB and 1 CPU. Anything higher silently multiplies the
  bill across all 2,917 tasks (§4.1). Also read the observed download rate out of the task log: it is the one
  number that could not be measured from the laptop, and it sets the real wall-clock.
- **after `small`** — any `.failed`? A few 404s are expected (retired SRX) and are fine; a *pattern* is not.
  Watch whether the endpoint slows under 4 concurrent fetches — if so, lower `%conc` or set `LIMIT_RATE`
  rather than pushing harder (§6).
- **after `medium`** — extrapolate the failure rate to 2,917 and sanity-check the QC spread. `report` prints
  the `dynamic_range_p99_over_median` distribution; if a large fraction sits near the low end, the
  `--min-peaks` gate will drop more samples than the pilot implied and the panel needs re-scoping *before*
  the full run, not after.

`scripts/44_local_rehearse.sh` runs `array.slurm` itself locally against the live endpoint with
`SLURM_ARRAY_TASK_ID` stubbed, so manifest chunking, resume and failure handling are exercised with no Roihu
access and no BU at all. Run it after any edit to the array body.

## 8.7 Live findings from the staged rollout (2026-07-21)

Ran smoke(2) → small(32) → medium(256) on Roihu. Everything the staged rollout was meant to surface, it did:

- **Environment.** Roihu has NO `python-data` module and NO tykky (both were on Puhti/Mahti); CPU nodes are
  **x86_64** (only GPU/GH200 are ARM). So the toolchain is a plain venv (`roihu/00_build_env.sh`), pyBigWig
  from a wheel, no compilation. `module` needs a *login* shell (`ssh roihu 'bash -lc "…"'`).
- **Memory / billing shape CONFIRMED.** MaxRSS 1.47 GB against the 2 GB request — the chunked-interval work
  put the job squarely in the 2 GiB/core free tier (§4.1). No BU surprise.
- **Download rate.** ~16 s/bigWig single-stream from a compute node (faster than the 30 s estimate). Under
  concurrency the per-sample time rises and **stragglers dominate**: at %4, task elapsed ranged 3:14–7:19 for
  the same 8-sample batch — an unlucky draw of large bigWigs. Wall-clock for the full run is set by the
  slowest task, not the mean; each task's worst case (~16×55 s ≈ 15 min) is far under the 4 h limit.
- **Failure rate: 0** across 256. The manifest is pre-filtered to grid-source samples, so 404s are rare.
- **QC bugs found and fixed (diagnostics only; atlas data unaffected).** `dynamic_range_p99_over_median` was
  `None` for samples covering few union-grid rows (median 0) → now uses the nonzero median; `frip_proxy`
  exceeded 1.0 (grid-vs-genome convention mismatch) → now archive-based, in [0,1]. Both were **recomputed for
  the already-pulled samples from the retained archive with no re-download** (`scripts/45_recompute_qc.py`) —
  the retention set doing exactly its job.
- **The `--min-peaks 2000` gate.** On the first 32 (a weak DRX-heavy prefix) it dropped 9/32; the medium
  stage gives the representative rate. Do not judge the gate from the manifest head — it is not randomized.

## 9. How this plugs into the existing scaffold

- `pilot/roihu/00_build_env.sh` / `env.yaml` — reuse; **add the SE caller** (LILY/ROSE2) + `bedtools`.
- `pilot/scripts/03_download_bigwigs.sh` — the per-SRX fetch (primary+legacy, resumable) is the feeder's
  inner loop; wrap it in the round/backpressure controller.
- `pilot/scripts/04_quantify.py` — `--mode bed --ref-bed grid.bed` is exactly step-2 of the loop.
- **BUILT 2026-07-20** (all under `phase2/`): `array.slurm` (chunked SLURM array) → `scripts/40_pull_one.sh`
  (the per-sample unit: download → one cnrose pass → provenance → delete, resumable via `.done`) →
  `reduce.slurm` (barrier reduce, refuses to run on an incomplete pull unless `FORCE=1`).
  `scripts/41_stage_and_submit.sh` stages and submits from the laptop (dry-run by default).
  `roihu/{env.yaml,00_build_env.sh}` is the fallback container — the dependency list is now just
  pyBigWig+numpy+scipy (no R, no deepTools), so `module load python-data` may suffice. No feeder needed (§5).
- Verified locally before any Roihu spend: a full pull of SRX067407 from the live endpoint reproduces the
  local SE calls exactly (724 SEs, cutoff 1707.72); resume skips in 2 ms; a bad SRX writes `.failed` and
  exits 0 without failing the array; no staging dirs leak; and the reduce over the sharded `.f32` layout
  reproduces the TSV-path atlas **bit-for-bit** (max diff 0).

**Next action (proposed):** do §7.1–7.2 locally (build grid + pull_set), then review §8 decisions, then
scaffold `array.slurm`. No heavy Roihu run until the §6 go/no-go.
