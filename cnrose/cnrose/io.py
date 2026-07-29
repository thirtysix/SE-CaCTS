"""bigWig quantification + BED I/O — the one signal primitive the whole caller shares.

ROSE ranks stitched regions by total signal ("total rpm") over the region span; we reproduce that as the
integral = mean-coverage x length via pyBigWig (agg='sum'). agg='mean' returns mean coverage instead.
The absolute scale is per-sample and does not affect the tangent cutoff (scale-invariant selection); the
cross-sample normalisation (S3norm) lives downstream on the fine matrix, not here (DESIGN.md §1, §7).
"""
from __future__ import annotations

import gzip
import json
import os

import numpy as np
import pyBigWig

# canonical human chromosomes (chr1..22, X, Y) — matches the grid build in phase2/scripts/21_*.sh
CANONICAL = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
CANON_ORDER = [f"chr{i}" for i in list(range(1, 23)) + ["X", "Y"]]


def read_bed3(path, canonical_only=False):
    """Read a BED into a list of (chrom, start, end); ignores extra columns and track/comment lines."""
    out = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("track", "#", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            c = f[0]
            if canonical_only and c not in CANONICAL:
                continue
            out.append((c, int(f[1]), int(f[2])))
    return out


def quantify(bw_path, regions, agg="sum"):
    """Signal per region from a bigWig.

    Parameters
    ----------
    bw_path : str
    regions : iterable of (chrom, start, end)
    agg : 'sum' (integral = mean x length, ROSE 'total rpm') | 'mean' (mean coverage)

    Returns np.ndarray[float], one value per region (0.0 for missing chrom / empty span).
    """
    regions = list(regions)
    out = np.zeros(len(regions), dtype=float)
    bw = pyBigWig.open(bw_path)
    try:
        chroms = bw.chroms()
        for i, (c, s, e) in enumerate(regions):
            clen = chroms.get(c)
            if clen is None:
                continue
            s2 = max(0, int(s))
            e2 = min(int(e), clen)
            if e2 <= s2:
                continue
            m = bw.stats(c, s2, e2, type="mean", exact=True)[0]
            if m is None:
                continue
            out[i] = float(m) * (e2 - s2) if agg == "sum" else float(m)
    finally:
        bw.close()
    return out


# --------------------------------------------------------------------------- durable outputs
# The Roihu pull deletes each bigWig after processing (PULL_DESIGN.md §8.5), so everything below exists to
# make that deletion safe: exact grid columns as compact binary, a genome-wide binned archive that lets any
# future region set be re-quantified without a re-download, and the QC numbers that only the bigWig can give.


def write_f32(path, values):
    """Exact grid column as raw little-endian float32, in grid-file row order (7.5x smaller than TSV)."""
    np.asarray(values, dtype="<f4").tofile(path)


def read_f32(path, n_expected=None):
    a = np.fromfile(path, dtype="<f4")
    if n_expected is not None and a.size != n_expected:
        raise ValueError(f"{path}: {a.size} values != expected {n_expected} (grid mismatch)")
    return a


def _pack(starts, ends, vals):
    """Pack coverage runs so any query resolves in O(log n).

    Two running integrals are needed, because a bigWig's `mean` is over COVERED bases only:
      F(x) = integral of value over [0, x)      (signal)
      G(x) = covered bases in [0, x)            (denominator)
    """
    w = ends - starts
    fa = np.cumsum(vals * w)
    ga = np.cumsum(w)
    return (starts, ends, vals,
            np.concatenate(([0.0], fa[:-1])), fa,
            np.concatenate(([0.0], ga[:-1])), ga)


def _eval_at(x, packed):
    """(F(x), G(x)) vectorised — signal integral and covered-base count up to each query point."""
    starts, ends, vals, fb, fa, gb, ga = packed
    k = np.searchsorted(starts, x, side="right") - 1
    F = np.zeros(x.size)
    G = np.zeros(x.size)
    ok = k >= 0
    kk = np.clip(k, 0, len(starts) - 1)
    inside = ok & (x < ends[kk])
    F[ok] = fa[kk[ok]]
    G[ok] = ga[kk[ok]]
    d = x[inside] - starts[kk[inside]]
    F[inside] = fb[kk[inside]] + vals[kk[inside]] * d
    G[inside] = gb[kk[inside]] + d
    return F, G


def _eval_chrom(bw, chrom, clen, pts, chunk):
    """(F, G) at every point in `pts` for one chromosome, reading the runs in windows.

    Reading a whole chromosome's runs at once peaks at ~2.7 GB on a dense bigWig (pyBigWig returns a Python
    list of tuples), which on a 2 GB/core node would bill each 1-core task as ~4 cores. Walking the
    chromosome in `chunk`-sized windows and carrying the running integrals across them holds the same exact
    result at a fraction of the memory — F and G are cumulative, so a per-window offset composes them.
    """
    order = np.argsort(pts, kind="stable")
    sp = pts[order]
    F = np.empty(sp.size)
    G = np.empty(sp.size)
    run_f = run_g = 0.0
    w0 = 0
    while w0 < clen:
        w1 = min(w0 + chunk, clen)
        iv = bw.intervals(chrom, int(w0), int(w1))
        lo = np.searchsorted(sp, w0, side="left")
        hi = np.searchsorted(sp, w1, side="left")
        if iv:
            a = np.asarray(iv, dtype=np.float64)
            st = np.maximum(a[:, 0], w0)            # pyBigWig returns ORIGINAL coords -> clip to the window
            en = np.minimum(a[:, 1], w1)
            keep = en > st
            st, en, vv = st[keep], en[keep], a[:, 2][keep]
            if st.size:
                packed = _pack(st, en, vv)
                if hi > lo:
                    f, g = _eval_at(sp[lo:hi], packed)
                    F[lo:hi] = run_f + f
                    G[lo:hi] = run_g + g
                run_f += float((vv * (en - st)).sum())
                run_g += float((en - st).sum())
            elif hi > lo:
                F[lo:hi] = run_f
                G[lo:hi] = run_g
        elif hi > lo:
            F[lo:hi] = run_f
            G[lo:hi] = run_g
        w0 = w1
    tail = np.searchsorted(sp, clen, side="left")   # points at/after the chromosome end
    if tail < sp.size:
        F[tail:] = run_f
        G[tail:] = run_g
    out_f = np.empty_like(F)
    out_g = np.empty_like(G)
    out_f[order] = F
    out_g[order] = G
    return out_f, out_g


def scan_bigwig(bw_path, region_sets=(), binsize=None, agg="sum", canonical_only=True,
                chunk=5_000_000):
    """ONE pass over the bigWig producing every durable artifact (PULL_DESIGN.md §8.5).

    `bw.intervals()` is the dominant cost, so each chromosome is read once and reused for every grid AND the
    genome-wide archive — quantifying two grids plus the archive separately would read the same runs three
    times (67 s/sample vs 17 s). Runs are read in `chunk`-sized windows to bound peak memory.

    Returns (signals, archive, layout): `signals` is one float64 array per region set, matching `quantify`
    EXACTLY (mean over covered bases x span, the ROSE convention cnrose was validated against); `archive`
    is the concatenated float32 bin means (uncovered = 0) or None.
    """
    region_sets = [list(r) for r in region_sets]
    signals = [np.zeros(len(r), dtype=float) for r in region_sets]
    per_chrom = []
    for regs in region_sets:
        d = {}
        for i, (c, s, e) in enumerate(regs):
            d.setdefault(c, []).append((i, int(s), int(e)))
        per_chrom.append(d)

    bw = pyBigWig.open(bw_path)
    parts, layout, off = [], [], 0
    try:
        chroms = bw.chroms()
        names = [c for c in CANON_ORDER if c in chroms] if canonical_only else list(chroms)
        for c in names:
            clen = chroms[c]
            # every query point for this chromosome, in one batch: region starts/ends + bin edges
            blocks, meta = [], []
            for si, d in enumerate(per_chrom):
                items = d.get(c)
                if not items:
                    continue
                idx = np.fromiter((i for i, _, _ in items), np.int64, len(items))
                st = np.clip(np.fromiter((s for _, s, _ in items), np.float64, len(items)), 0, clen)
                en = np.clip(np.fromiter((e for _, _, e in items), np.float64, len(items)), 0, clen)
                meta.append((si, idx, len(blocks), st, en))
                blocks.append(st)
                blocks.append(en)
            nb = int(np.ceil(clen / binsize)) if binsize else 0
            edges = (np.minimum(np.arange(nb + 1, dtype=np.float64) * binsize, clen)
                     if binsize else np.empty(0))
            if binsize:
                blocks.append(edges)
            if not blocks:
                continue
            sizes = [b.size for b in blocks]
            F, G = _eval_chrom(bw, c, clen, np.concatenate(blocks), chunk)
            cuts = np.cumsum([0] + sizes)
            for si, idx, bi, st, en in meta:
                Fs, Gs = F[cuts[bi]:cuts[bi + 1]], G[cuts[bi]:cuts[bi + 1]]
                Fe, Ge = F[cuts[bi + 1]:cuts[bi + 2]], G[cuts[bi + 1]:cuts[bi + 2]]
                integral, covered, span = Fe - Fs, Ge - Gs, en - st
                mean_cov = np.where(covered > 0, integral / np.maximum(covered, 1e-12), 0.0)
                signals[si][idx] = mean_cov * span if agg == "sum" else mean_cov
            if binsize:
                Fe = F[cuts[len(blocks) - 1]:]
                wdt = np.diff(edges)
                # archive convention: uncovered = 0, i.e. mean over the FULL bin span
                m = np.where(wdt > 0, np.diff(Fe) / np.maximum(wdt, 1e-9), 0.0).astype(np.float32)
                parts.append(m)
                layout.append({"chrom": c, "offset": off, "nbins": int(m.size), "length": int(clen)})
                off += int(m.size)
    finally:
        bw.close()
    archive = np.concatenate(parts) if parts else None
    return signals, archive, layout


def bin_chrom(bw, chrom, binsize, chunk=5_000_000):
    """EXACT per-bin mean coverage for one chromosome (uncovered = 0), no per-bp arrays."""
    clen = bw.chroms().get(chrom)
    if not clen:
        return None
    nb = int(np.ceil(clen / binsize))
    edges = np.minimum(np.arange(nb + 1, dtype=np.float64) * binsize, clen)
    F, _ = _eval_chrom(bw, chrom, clen, edges, chunk)
    w = np.diff(edges)
    return np.where(w > 0, np.diff(F) / np.maximum(w, 1e-9), 0.0).astype(np.float32)


def write_binned_archive(values, layout, out_path, binsize=100):
    """Write the genome-wide binned archive: float16, gzipped, chroms concatenated in CANON_ORDER.

    This is the artifact that retires the irreversible grid-threshold choice (PULL_DESIGN.md §8.5.2). It is
    INSURANCE, not a replacement for exact quantification: per-row reconstruction error is biased and
    accumulates (see phase2/analysis/binned_storage_eval.py). Sidecar .json records the layout.
    """
    allv = np.asarray(values, dtype=np.float16)
    with gzip.open(out_path, "wb", compresslevel=6) as fh:
        fh.write(allv.tobytes())
    meta = {"binsize": binsize, "dtype": "float16", "order": "C", "total_bins": int(allv.size),
            "layout": layout}
    with open(out_path + ".json", "w") as fh:
        json.dump(meta, fh, indent=1)
    return meta


def read_binned_archive(path):
    """Return (values float16 array, meta dict) for an archive written by write_binned_archive."""
    with open(path + ".json") as fh:
        meta = json.load(fh)
    with gzip.open(path, "rb") as fh:
        v = np.frombuffer(fh.read(), dtype=np.float16)
    if v.size != meta["total_bins"]:
        raise ValueError(f"{path}: {v.size} bins != meta {meta['total_bins']}")
    return v, meta


def coarsen_archive(values, layout, factor):
    """Re-bin a genome-wide archive to `factor`x coarser bins, WITHIN each chromosome (never across a
    boundary). A coarse bin's value is the MEAN of its constituent fine bins; since the fine bins are
    equal-width means with uncovered=0, that mean is the coarse bin's mean over its span, preserving the
    archive's convention. The last (partial) coarse bin on each chromosome averages only its ACTUAL fine
    bins (not zero-padded), so telomere ends are not biased low.

    This derives e.g. a 1 kb archive (factor 10 from 100 bp) from the RETAINED archive with no re-pull —
    ~1000x finer than megabase-scale CN needs, at 1/`factor` the size (PULL_DESIGN.md §8.5.6). Returns
    (new_values float16, new_layout); new binsize = old binsize * factor.
    """
    v = np.asarray(values, dtype=np.float32)
    parts, new_layout, off = [], [], 0
    for L in layout:
        seg = v[L["offset"]: L["offset"] + L["nbins"]]
        n = seg.size
        full = (n // factor) * factor
        chunks = []
        if full:
            chunks.append(seg[:full].reshape(-1, factor).mean(axis=1))
        if full < n:
            chunks.append(np.array([seg[full:].mean()], dtype=np.float32))   # partial tail bin
        m = (np.concatenate(chunks) if chunks else np.zeros(0)).astype(np.float16)
        parts.append(m)
        new_layout.append({"chrom": L["chrom"], "offset": off, "nbins": int(m.size), "length": L["length"]})
        off += int(m.size)
    return (np.concatenate(parts) if parts else np.zeros(0, np.float16)), new_layout


def _frip_from_archive(binned, binsize, layout, grid_regions):
    """Fraction of total binned signal that falls in grid regions — numerator and denominator both from the
    archive (uncovered=0, full-bin mean), so the result is a true fraction in [0,1]. Avoids the convention
    mismatch of dividing a ROSE covered-mean grid score by a full-bin genome total (which can exceed 1)."""
    total = float(np.asarray(binned, dtype=np.float64).sum())
    if total <= 0:
        return None
    off = {L["chrom"]: (L["offset"], L["nbins"]) for L in layout}
    by_chrom = {}
    for c, s, e in grid_regions:
        by_chrom.setdefault(c, []).append((int(s), int(e)))
    mask = np.zeros(len(binned), dtype=bool)
    for c, ivs in by_chrom.items():
        if c not in off:
            continue
        base, nb = off[c]
        for s, e in ivs:
            i0 = max(0, s // binsize)
            i1 = min(nb, (e - 1) // binsize + 1)
            if i1 > i0:
                mask[base + i0: base + i1] = True
    return float(np.asarray(binned, dtype=np.float64)[mask].sum() / total)


def qc_stats(bw_path, binned, binsize, grid_regions=None, grid_signal=None, layout=None):
    """Per-sample QC that ONLY the bigWig can provide — emit before deleting it (PULL_DESIGN.md §8.5.1).

    The **operative** QC gate is the MACS2 peak count (`n_peaks`, added by the caller) via
    `aggregate.py --min-peaks`; it is unambiguous at any grid scale. `dynamic_range_p99_over_median` is a
    secondary quality DIAGNOSTIC: on the pilot the *plain* p99/median tracked the S3norm exponent (Spearman
    -0.956), i.e. a low value flags a compressed, low-enrichment sample. But the plain median breaks at scale:
    the production grid is the UNION of all 2,917 samples' peaks, so any one sample has signal in only a
    fraction of its rows and the median hits 0 (undefined ratio) — seen on the first live samples. So the
    reported value uses the NONZERO median (always defined; == the pilot value when few rows are zero, so the
    -0.956 was really on nonzero≈plain there); `dynamic_range_p99_over_median_raw` keeps the plain value for
    continuity. Treat it as an ordinal quality hint, not the gate.

    `frip_proxy` is fraction-of-signal-in-grid, computed from the ARCHIVE for BOTH numerator and denominator
    so it stays in [0,1]. It must NOT mix the grid column (ROSE covered-mean x span) with the genome total
    (full-bin mean): those are different conventions and their ratio can exceed 1 (seen live: FRiP 1.14).
    Genome-bin stats live in a separate `genome_bins` block so they are never confused with the grid.
    """
    bw = pyBigWig.open(bw_path)
    try:
        hdr = bw.header()
    finally:
        bw.close()
    out = {"bigwig_bytes": os.path.getsize(bw_path),
           "bigwig_header": {k: hdr[k] for k in ("nBasesCovered", "sumData", "maxVal", "minVal") if k in hdr}}
    out.update(qc_from_retained(binned, binsize, grid_regions=grid_regions, grid_signal=grid_signal,
                                layout=layout))
    return out


def qc_from_retained(binned, binsize, grid_regions=None, grid_signal=None, layout=None):
    """The bigWig-free core of qc_stats — computes every statistic from RETAINED artifacts alone (the
    archive, the grid f32 column, the grid BED, the archive layout). This is what lets QC be recomputed for
    already-pulled samples without re-downloading (PULL_DESIGN.md §8.5): `scripts/45_recompute_qc.py` calls
    it. Does NOT include the bigwig_header / bigwig_bytes fields (those need the bigWig); merge those from the
    existing qc.json if you want to preserve them."""
    pcts = [1, 5, 25, 50, 75, 90, 95, 99, 99.9]
    out = {}
    if binned is not None:
        v = np.asarray(binned, dtype=np.float64)
        q = np.percentile(v, pcts)
        nz = v[v > 0]
        med_nz = float(np.median(nz)) if nz.size else 0.0
        out["genome_bins"] = {
            "bin_size": binsize,
            "n_bins": int(v.size),
            "total_signal": float(v.sum() * binsize),
            "mean": float(v.mean()),
            "max": float(v.max()),
            "zero_fraction": float((v == 0).mean()),
            "percentiles": {str(p): float(x) for p, x in zip(pcts, q)},
            "dynamic_range_p99_over_nonzero_median": (float(q[pcts.index(99)] / med_nz)
                                                      if med_nz > 0 else None),
        }
        out["zero_fraction"] = out["genome_bins"]["zero_fraction"]
        out["total_signal"] = out["genome_bins"]["total_signal"]

    if grid_signal is not None and len(grid_signal):
        g = np.asarray(grid_signal, dtype=np.float64)
        gq = np.percentile(g, pcts)
        gmed = float(gq[pcts.index(50)])
        gnz = g[g > 0]
        gmed_nz = float(np.median(gnz)) if gnz.size else 0.0
        p99 = float(gq[pcts.index(99)])
        out["grid"] = {
            "n_rows": int(g.size),
            "nonzero_fraction": float((g > 0).mean()),
            "bp": int(sum(e - s for _, s, e in grid_regions)) if grid_regions is not None else None,
            "signal_sum": float(g.sum()),
            "percentiles": {str(p): float(x) for p, x in zip(pcts, gq)},
        }
        # THE gate statistic: p99 over the NONZERO median (well-defined on the union grid; == the pilot-
        # validated p99/median when few rows are zero). Spearman -0.956 vs the S3norm exponent.
        out["dynamic_range_p99_over_median"] = (p99 / gmed_nz) if gmed_nz > 0 else None
        out["dynamic_range_p99_over_median_raw"] = (p99 / gmed) if gmed > 0 else None
        if binned is not None and layout and grid_regions is not None:
            out["frip_proxy"] = _frip_from_archive(binned, binsize, layout, grid_regions)
    return out
