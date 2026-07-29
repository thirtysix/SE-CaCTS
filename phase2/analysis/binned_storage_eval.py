#!/usr/bin/env python3
"""Can a genome-wide BINNED archive replace exact grid quantification? (Decision evidence, 2026-07-20.)

Motivating question: the pull deletes each bigWig, so the peak-threshold grid choice is irreversible. If we
instead stored a reduced/binned copy of the whole genome, no future region-set question would need a
re-download and the choice would stop mattering. This measures whether that works.

Signal is an integral (mean x length, cnrose/io.py:56) and additive, so reconstruction error enters only at
the two partial end bins => it scales with bin size / region size. Bins here are EXACT means built from
per-bp values (bw.values) on the global bin lattice, stored as float16 — no zoom-level approximation, so
this isolates binning + float16 loss alone.

RESULT (ALL 13 pilot bigWigs; ranges across samples): binning is excellent at SE scale and poor at grid-row
scale, and crucially the per-row error is BIASED (signed mean ~= absolute mean, POSITIVE in 13/13 samples),
so it ACCUMULATES across the rows of an SE rather than cancelling. The grid is the union across all 2,917
samples, so most rows are not a given sample's own peaks and neighbouring bins bleed signal in —
reconstruction systematically overestimates.

  bin    SE-sum median   SE-sum p90    SE-sum WORST   signed bias
  100    0.29-0.98%      1.1-4.2%      4.4-20%        +1.6 to +6.4%
  200    0.84-2.17%      2.7-9.6%      12-65%         +3.6 to +18.6%
  500    2.50-5.01%      6.9-19%       22-124%        +7.5 to +52%
  1000   3.97-9.40%      10-32%        30-268%        +12 to +70%

=> Bins CANNOT replace exact grid quantification (cnrose is validated bit-for-bit; 1-4% error forfeits that).
=> Bins ARE good insurance: any future region set is recoverable to ~1% median at 100 bp, which is ample to
   decide whether a threshold change justifies a re-pull, and exact for contiguous spans like stitched
   regions. See PULL_DESIGN.md §8.5.

  ~/miniconda3/envs/atac_hdac/bin/python phase2/analysis/binned_storage_eval.py [SRX ...]   # ~5 min

---
Original test description:
The decision-relevant test: can genome-wide bins reproduce the ATLAS quantity?

aggregate.py computes an SE's signal as the sum over the grid rows OVERLAPPING it — NOT the integral over
the contiguous span, because the gaps between constituent peaks are excluded. So per-row reconstruction
errors either cancel across the ~tens of rows in an SE (if unbiased) or accumulate (if biased). Individual
rows reconstruct poorly (p90 14-21% at 100 bp), so this is what decides whether binned storage can replace
storing the grid column outright.

Reports, per bin size: error of the summed-over-rows SE signal, and the SIGNED mean per-row error (the bias
term — if this is non-zero, errors accumulate rather than cancel).
"""
import os
import sys

import numpy as np
import pyBigWig

# The repo itself — derived from this file's location, so no local path is baked in.
BASE = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
BW_DIR = os.path.join(BASE, "pilot/data/bigwigs")
GRID = os.path.join(BASE, "phase2/data/grid.20.bed")
BINS = [100, 200, 500, 1000]


def exact_bins(vals, b):
    n = (len(vals) // b) * b
    return vals[:n].reshape(-1, b).mean(axis=1)


def recon(bm, w0, b, s, e):
    total = 0.0
    for i in range((s - w0) // b, (e - 1 - w0) // b + 1):
        if 0 <= i < len(bm):
            bs, be = w0 + i * b, w0 + (i + 1) * b
            total += float(bm[i]) * (min(e, be) - max(s, bs))
    return total


def main():
    srxs = sys.argv[1:] or ["SRX067407", "SRX8588595"]
    grid = []
    with open(GRID) as fh:
        for line in fh:
            f = line.split("\t")
            grid.append((f[0], int(f[1]), int(f[2])))
    rng = np.random.default_rng(1)
    # SE-like groups: runs of grid rows spanning 15-45 kb (the rows themselves, gaps EXCLUDED)
    groups, i = [], 0
    while len(groups) < 250 and i < len(grid) - 2:
        c, s0, _ = grid[i]
        j = i
        while j < len(grid) - 1 and grid[j + 1][0] == c and grid[j + 1][2] - s0 < 45_000:
            j += 1
        if grid[j][2] - s0 > 15_000 and j > i + 3:
            groups.append(grid[i:j + 1])
        i = j + 1 + int(rng.integers(1, 250))
    nrows = np.mean([len(g) for g in groups])
    print(f"{len(groups)} SE-like groups, mean {nrows:.1f} grid rows each, "
          f"mean span {np.mean([g[-1][2]-g[0][1] for g in groups])/1000:.1f} kb "
          f"(rows cover {np.mean([sum(e-s for _,s,e in g)/(g[-1][2]-g[0][1]) for g in groups])*100:.0f}% of span)\n")

    for srx in srxs:
        bw = pyBigWig.open(os.path.join(BW_DIR, f"{srx}.bw"))
        print(f"===== {srx} =====")
        print(f"  {'bin':>6} {'SE-sum med':>11} {'SE-sum p90':>11} {'SE-sum max':>11} "
              f"{'per-row SIGNED mean':>20} {'per-row |err| med':>18}")
        for b in BINS:
            se_err, row_signed, row_abs = [], [], []
            for g in groups:
                c = g[0][0]
                w0 = (g[0][1] // b) * b
                w1 = ((g[-1][2] - 1) // b + 1) * b
                if w1 > bw.chroms().get(c, 0):
                    continue
                vals = np.nan_to_num(bw.values(c, w0, w1, numpy=True)).astype(np.float64)
                bm = exact_bins(vals, b).astype(np.float16).astype(np.float64)
                truth_sum = rec_sum = 0.0
                for _, s, e in g:
                    t = float(vals[s - w0:e - w0].sum())
                    r = recon(bm, w0, b, s, e)
                    truth_sum += t
                    rec_sum += r
                    if t > 0:
                        row_signed.append((r - t) / t)
                        row_abs.append(abs(r - t) / t)
                if truth_sum > 0:
                    se_err.append(abs(rec_sum - truth_sum) / truth_sum)
            se_err = np.array(se_err)
            print(f"  {b:>6} {np.median(se_err)*100:>10.2f}% {np.percentile(se_err,90)*100:>10.2f}% "
                  f"{se_err.max()*100:>10.2f}% {np.mean(row_signed)*100:>19.2f}% "
                  f"{np.median(row_abs)*100:>17.2f}%")
        bw.close()
        print()


if __name__ == "__main__":
    main()
