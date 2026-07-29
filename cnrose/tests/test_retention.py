#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests for the PULL_DESIGN.md §8.5 retention outputs — the artifacts that make deleting a bigWig safe.

The load-bearing assertion is that `scan_bigwig` reproduces `quantify` EXACTLY. `scan_bigwig` is the
single-pass rewrite that produces the calling signal, every grid column and the genome-wide archive from one
read of the runs; `quantify` is the original per-region bw.stats() loop that cnrose was validated bit-for-bit
against ROSE2 with. If they ever diverge, that validation no longer transfers.

Subtlety worth preserving: a bigWig's `mean` is over COVERED bases only, so the signal is
(integral / covered_bp) x span — NOT the plain integral. Computing the plain integral silently changes every
partially-covered region. The archive uses the opposite convention (uncovered = 0, mean over the full bin),
which is correct for it and is asserted separately below.

  ~/miniconda3/envs/atac_hdac/bin/python cnrose/tests/test_retention.py
"""
from __future__ import annotations

import json
import os
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

from cnrose.io import (quantify, scan_bigwig, read_bed3, write_f32, read_f32,      # noqa: E402
                       write_binned_archive, read_binned_archive, qc_stats, bin_chrom)

SEC = os.path.dirname(HERE)
BW = os.path.join(SEC, "..", "pilot/data/bigwigs/SRX067407.bw")
GRID = os.path.join(SEC, "..", "phase2/data/grid.20.bed")


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}{'  — ' + detail if detail else ''}")
    return bool(ok)


def main():
    if not os.path.exists(BW):
        print(f"SKIP: {BW} not present (re-fetch with pilot/scripts/03_download_bigwigs.sh)")
        return 0
    ok = []

    # a subset of the real grid keeps the test quick while staying on real coordinates
    grid = read_bed3(GRID, canonical_only=True)
    sub = [g for g in grid if g[0] in ("chr20", "chr21")]
    sigs, archive, layout = scan_bigwig(BW, [sub], binsize=100, canonical_only=True)
    ref = quantify(BW, sub)
    rel = np.abs(sigs[0] - ref) / np.maximum(np.abs(ref), 1e-9)
    ok.append(check("scan_bigwig reproduces quantify EXACTLY", float(rel.max()) == 0.0,
                    f"n={len(sub)} regions, max rel diff={rel.max():.2e}"))

    # multiple region sets in one pass must each match their own independent quantification
    sub2 = [(c, s, min(e, s + 200)) for c, s, e in sub[:5000]]
    m, _, _ = scan_bigwig(BW, [sub, sub2], binsize=None, canonical_only=True)
    r2 = quantify(BW, sub2)
    ok.append(check("multi-grid pass: 2nd grid also exact",
                    float(np.abs(m[1] - r2).max()) == 0.0, f"n={len(sub2)}"))
    ok.append(check("multi-grid pass: 1st grid unaffected by the 2nd",
                    float(np.abs(m[0] - sigs[0]).max()) == 0.0))

    ok.append(check("agg='mean' consistent with agg='sum'",
                    np.allclose(scan_bigwig(BW, [sub2], agg="mean", canonical_only=True)[0][0],
                                r2 / np.array([e - s for _, s, e in sub2]), atol=1e-9)))

    with tempfile.TemporaryDirectory() as td:
        # f32 round-trip is exact and length-checked
        p = os.path.join(td, "c.f32")
        write_f32(p, sigs[0])
        back = read_f32(p, n_expected=len(sub))
        ok.append(check("f32 round-trip exact at float32 precision",
                        np.array_equal(back, sigs[0].astype("<f4"))))
        try:
            read_f32(p, n_expected=len(sub) + 1)
            ok.append(check("f32 length mismatch raises", False))
        except ValueError:
            ok.append(check("f32 length mismatch raises", True))

        # archive round-trip + layout
        ap = os.path.join(td, "a.f16.gz")
        meta = write_binned_archive(archive, layout, ap, binsize=100)
        v, m2 = read_binned_archive(ap)
        ok.append(check("archive round-trip", v.size == meta["total_bins"] and m2["binsize"] == 100,
                        f"{v.size} bins"))
        ok.append(check("archive layout offsets are contiguous and complete",
                        all(a["offset"] + a["nbins"] == b["offset"]
                            for a, b in zip(m2["layout"], m2["layout"][1:]))
                        and m2["layout"][-1]["offset"] + m2["layout"][-1]["nbins"] == v.size))
        ok.append(check("archive is float16 (2 bytes/bin)", v.dtype == np.float16))

    # archive convention: bin mean over the FULL span (uncovered = 0), != quantify's covered-only mean
    import pyBigWig
    bw = pyBigWig.open(BW)
    m20 = bin_chrom(bw, "chr20", 100)
    vals = np.nan_to_num(bw.values("chr20", 5_000_000, 5_100_000, numpy=True)).astype(np.float64)
    bw.close()
    truth = vals.reshape(-1, 100).mean(axis=1)
    ok.append(check("archive bins are exact vs per-bp truth (uncovered=0)",
                    float(np.abs(truth - m20[50_000:51_000]).max()) < 1e-6,
                    f"max diff={np.abs(truth - m20[50_000:51_000]).max():.2e}"))

    # QC: the gate statistic must come from the GRID column, not genome bins (bins are >50% zero)
    q = qc_stats(BW, archive, 100, grid_regions=sub, grid_signal=sigs[0], layout=layout)
    dr = q.get("dynamic_range_p99_over_median")
    ok.append(check("QC gate stat (dynamic range) is defined and sane",
                    dr is not None and 1 < dr < 5000, f"p99/nonzero-median={dr:.1f}"))
    # FRiP must be a true fraction in [0,1] (the old grid/genome-total mix could exceed 1)
    fr = q.get("frip_proxy")
    ok.append(check("FRiP proxy is a fraction in [0,1]", fr is not None and 0 <= fr <= 1, f"FRiP={fr:.3f}"))
    # dynamic range must stay defined even when many grid rows are zero (union-grid at scale)
    g2 = sigs[0].copy(); g2[: len(g2) // 2] = 0.0            # force half the rows to zero
    q2 = qc_stats(BW, archive, 100, grid_regions=sub, grid_signal=g2, layout=layout)
    ok.append(check("dynamic range defined despite 50% zero rows (nonzero-median)",
                    q2.get("dynamic_range_p99_over_median") is not None,
                    f"raw={q2.get('dynamic_range_p99_over_median_raw')}"))
    ok.append(check("QC keeps genome-bin stats separate from grid stats",
                    "genome_bins" in q and "grid" in q
                    and q["genome_bins"]["zero_fraction"] > 0.5,
                    f"genome bins {q['genome_bins']['zero_fraction']*100:.0f}% zero — "
                    f"why the median there is unusable"))
    ok.append(check("QC json is serialisable", isinstance(json.dumps(q), str)))

    print(f"\n{sum(ok)}/{len(ok)} passed")
    return 0 if all(ok) else 1


if __name__ == "__main__":
    sys.exit(main())
