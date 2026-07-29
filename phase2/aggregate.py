#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase-2 barrier reduce: per-sample cnrose outputs -> union SE catalog + region x sample matrices.

Runs ONCE after the array pull (DESIGN.md §9, PULL_DESIGN.md §5.1 shard-then-reduce). Reads every
<SRX>.se.bed (per-sample super-enhancer calls) and <SRX>.signal.tsv (grid-row signal, one column per
sample over the FIXED grid), plus the grid BED, and writes:

  <out>.union_catalog.bed     union SE loci (cSEAdb-style: >=25% RECIPROCAL-overlap single-linkage merge)
  <out>.se_signal.tsv         union SE x sample signal = sum of the SE's grid-row signals (Σ meanᵢ·lenᵢ)
  <out>.se_presence.tsv       union SE x sample 0/1: did that sample call an SE overlapping this locus

Because the grid is fixed and every per-sample SE is a union of grid rows, the region x sample matrix is
recovered by summing fine-matrix rows — no second bigWig pass (PULL_DESIGN.md §3). Separable-CN: the fine
matrix stays UNCORRECTED here; CN correction is applied downstream (DESIGN.md §3.2), and the dual
(agnostic vs CN-corrected) catalogs are produced by running this once per per-sample catalog.

`--norm s3norm` applies cross-study normalization to the FINE matrix before the SE sum (writes
<out>.s3norm_params.tsv). It belongs here, not downstream, because the S3norm transform is nonlinear:
Σf(xᵢ) ≠ f(Σxᵢ), so it must act at grid-row resolution and be re-summed. See s3norm.py.

Run with atac_hdac (numpy):  ~/miniconda3/envs/atac_hdac/bin/python phase2/aggregate.py --help
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np


def read_grid(path):
    """Return (chroms_list, starts, ends) as parallel arrays in file order (the fine-matrix row order)."""
    chroms, starts, ends = [], [], []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("track", "#", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            chroms.append(f[0]); starts.append(int(f[1])); ends.append(int(f[2]))
    return chroms, np.asarray(starts), np.asarray(ends)


def grid_index(chroms, starts, ends):
    """Per-chrom sorted view for fast interval->row lookup. Grid rows are non-overlapping & sorted."""
    idx = {}
    order = np.arange(len(chroms))
    chroms = np.asarray(chroms)
    for c in np.unique(chroms):
        m = chroms == c
        rows = order[m]
        # already sorted within chrom (grid build sorts by start); keep the global row ids
        idx[c] = (starts[m], ends[m], rows)
    return idx


def rows_for_region(gidx, chrom, s, e):
    """Global grid-row ids overlapping [s,e) on chrom (half-open BED coords)."""
    if chrom not in gidx:
        return np.empty(0, dtype=int)
    gs, ge, rows = gidx[chrom]
    lo = np.searchsorted(ge, s, side="right")   # first row whose end > s
    hi = np.searchsorted(gs, e, side="left")    # first row whose start >= e
    return rows[lo:hi]


def read_se_bed(path):
    out = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("track", "#", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            out.append((f[0], int(f[1]), int(f[2])))
    return out


def read_signal_col(path, n_expected):
    """Grid column from either <SRX>.signal.tsv (text) or <SRX>.<label>.f32 (raw little-endian float32)."""
    if path.endswith(".f32"):
        a = np.fromfile(path, dtype="<f4")
        if a.size != n_expected:
            raise ValueError(f"{path}: {a.size} rows != grid {n_expected} (grid mismatch)")
        return a.astype(np.float32)
    return _read_signal_tsv(path, n_expected)


def _read_signal_tsv(path, n_expected):
    """Read the SIGNAL column of a <SRX>.signal.tsv (header CHROM START STOP SIGNAL). Assert row count."""
    vals = np.empty(n_expected, dtype=np.float32)
    i = 0
    with open(path) as fh:
        header = fh.readline()
        for line in fh:
            vals[i] = float(line.rstrip("\n").rsplit("\t", 1)[1])
            i += 1
    if i != n_expected:
        raise ValueError(f"{path}: {i} rows != grid {n_expected} (grid mismatch)")
    return vals


def build_union(se_by_sample, recip=0.25):
    """cSEAdb-style union: single-linkage merge of per-sample SEs with >=recip RECIPROCAL overlap.

    Returns list of (chrom, start, end, sample_set, n_members), sorted by (chrom, start).
    """
    by_chrom = {}
    for si, ses in enumerate(se_by_sample):
        for c, s, e in ses:
            by_chrom.setdefault(c, []).append((s, e, si))

    unions = []
    for chrom, ivs in by_chrom.items():
        ivs.sort()
        n = len(ivs)
        parent = list(range(n))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for i in range(n):
            si, ei, _ = ivs[i]
            li = ei - si
            for j in range(i + 1, n):
                sj, ej, _ = ivs[j]
                if sj >= ei:            # sorted by start -> no further overlaps with i
                    break
                ov = min(ei, ej) - sj   # sj >= si so max(si,sj)=sj
                if ov <= 0:
                    continue
                if ov / li >= recip and ov / (ej - sj) >= recip:
                    ra, rb = find(i), find(j)
                    if ra != rb:
                        parent[rb] = ra

        comp = {}
        for i in range(n):
            comp.setdefault(find(i), []).append(i)
        for members in comp.values():
            s = min(ivs[k][0] for k in members)
            e = max(ivs[k][1] for k in members)
            samps = {ivs[k][2] for k in members}
            unions.append((chrom, s, e, samps, len(members)))

    unions.sort(key=lambda u: (u[0], u[1]))
    return unions


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--se-dir", required=True, help="dir of <SRX>{suffix}.se.bed")
    ap.add_argument("--signal-dir", required=True, help="dir of <SRX>.signal.tsv")
    ap.add_argument("--grid", required=True, help="the fixed grid BED (grid.20.bed)")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--recip", type=float, default=0.25, help="reciprocal overlap threshold (default 0.25)")
    ap.add_argument("--se-suffix", default="", help="SE-catalog suffix: '' (agnostic) or '.cn' (corrected). "
                    "Signal always comes from the UNCORRECTED <SRX>.signal.tsv (CN is separable).")
    ap.add_argument("--norm", default="none", choices=["none", "s3norm"],
                    help="cross-study normalization applied to the FINE matrix before the SE sum. "
                         "s3norm is nonlinear, so it must be applied at grid resolution and re-summed "
                         "(Sum f(x) != f(Sum x)) — see s3norm.py. The SE catalog is unaffected.")
    ap.add_argument("--s3-ref", default="medoid", help="s3norm reference: 'medoid', 'mean', or a sample id")
    ap.add_argument("--grid-label", default=None,
                    help="read <SRX>.<label>.f32 exact columns instead of <SRX>.signal.tsv (e.g. 'grid.20'). "
                         "This is what the pull writes — see PULL_DESIGN.md §8.5.5.")
    ap.add_argument("--peak-dir", help="dir of <SRX>.20.bed, for the --min-peaks QC gate")
    ap.add_argument("--min-peaks", type=int, default=0,
                    help="drop samples with fewer than N MACS2 peaks (requires --peak-dir). STRONGLY "
                         "RECOMMENDED with --norm s3norm: low-enrichment samples have compressed dynamic "
                         "range, so s3norm fits an expanding exponent (B>1) that amplifies their noise into "
                         "apparent specificity. Pilot: a natural gap at 2046 vs 361 peaks; without the gate "
                         "s3norm pushed Breast ESR1 from rank #5 to #21, with it ESR1 is #2.")
    a = ap.parse_args()

    chroms, starts, ends = read_grid(a.grid)
    n_grid = len(chroms)
    gidx = grid_index(chroms, starts, ends)
    print(f"[aggregate] grid: {n_grid} rows", file=sys.stderr)

    # samples = SRX with BOTH a se.bed (of the requested catalog) and a signal.tsv.
    # SE regions come from <SRX>{suffix}.se.bed; signal always from the UNCORRECTED <SRX>.signal.tsv.
    se_tail = f"{a.se_suffix}.se.bed"
    se_files = {}
    for p in glob.glob(os.path.join(a.se_dir, f"*{se_tail}")):
        sample = os.path.basename(p)[:-len(se_tail)]
        if a.se_suffix == "" and sample.endswith(".cn"):     # don't let the agnostic glob grab .cn files
            continue
        se_files[sample] = p
    # accept either the historical <SRX>.signal.tsv or the compact exact <SRX>.<gridlabel>.f32 the pull
    # writes (PULL_DESIGN.md §8.5.5). f32 wins when both are present — same values, 7.5x smaller.
    sig_files = {os.path.basename(p)[:-len(".signal.tsv")]: p
                 for p in glob.glob(os.path.join(a.signal_dir, "*.signal.tsv"))}
    if a.grid_label:
        tail = f".{a.grid_label}.f32"
        sig_files.update({os.path.basename(p)[:-len(tail)]: p
                          for p in glob.glob(os.path.join(a.signal_dir, f"*{tail}"))})
    samples = sorted(set(se_files) & set(sig_files))
    if not samples:
        sys.exit("[aggregate] no samples with both .se.bed and .signal.tsv")

    if a.min_peaks:
        if not a.peak_dir:
            sys.exit("[aggregate] --min-peaks requires --peak-dir")
        kept, dropped = [], []
        for srx in samples:
            p = os.path.join(a.peak_dir, f"{srx}.20.bed")
            n = sum(1 for _ in open(p)) if os.path.exists(p) else 0
            (kept if n >= a.min_peaks else dropped).append((srx, n))
        for srx, n in dropped:
            print(f"[aggregate] QC drop {srx}: {n} peaks < {a.min_peaks}", file=sys.stderr)
        samples = [s for s, _ in kept]
        if not samples:
            sys.exit(f"[aggregate] --min-peaks {a.min_peaks} dropped every sample")
    print(f"[aggregate] {len(samples)} samples", file=sys.stderr)

    # fine matrix (n_grid x n_samples), float32
    M = np.empty((n_grid, len(samples)), dtype=np.float32)
    for j, srx in enumerate(samples):
        M[:, j] = read_signal_col(sig_files[srx], n_grid)

    # cross-study normalization at GRID resolution (before the SE sum — the transform is nonlinear)
    if a.norm == "s3norm":
        from s3norm import s3norm_matrix
        M, s3params = s3norm_matrix(M, ref=a.s3_ref, srx=samples)
        s3params.to_csv(a.out + ".s3norm_params.tsv", sep="\t", index=False)

    # union catalog
    se_by_sample = [read_se_bed(se_files[srx]) for srx in samples]
    unions = build_union(se_by_sample, recip=a.recip)
    print(f"[aggregate] union catalog: {len(unions)} SE loci "
          f"(from {sum(len(x) for x in se_by_sample)} per-sample SE calls)", file=sys.stderr)

    # SE x sample signal (Σ grid-row signal) + presence
    se_signal = np.zeros((len(unions), len(samples)), dtype=np.float64)
    grid_rows_per_se = []
    for i, (c, s, e, samps, _n) in enumerate(unions):
        rows = rows_for_region(gidx, c, s, e)
        grid_rows_per_se.append(rows)
        if rows.size:
            se_signal[i, :] = M[rows, :].sum(axis=0)

    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    catalog_bed = a.out + ".union_catalog.bed"
    with open(catalog_bed, "w") as fh:
        for i, (c, s, e, samps, nmem) in enumerate(unions):
            fh.write(f"{c}\t{s}\t{e}\tUSE_{i}\t{len(samps)}\t{nmem}\t{grid_rows_per_se[i].size}\n")

    sig_tsv = a.out + ".se_signal.tsv"
    with open(sig_tsv, "w") as fh:
        fh.write("se_id\t" + "\t".join(samples) + "\n")
        for i in range(len(unions)):
            fh.write(f"USE_{i}\t" + "\t".join(f"{v:.6g}" for v in se_signal[i, :]) + "\n")

    pres_tsv = a.out + ".se_presence.tsv"
    with open(pres_tsv, "w") as fh:
        fh.write("se_id\t" + "\t".join(samples) + "\n")
        for i, (c, s, e, samps, _n) in enumerate(unions):
            fh.write(f"USE_{i}\t" + "\t".join("1" if j in samps else "0" for j in range(len(samples))) + "\n")

    # ---- self-validation of the aggregation plumbing ----
    rng = np.random.default_rng(0)
    checks = rng.choice(len(unions), size=min(200, len(unions)), replace=False)
    max_err = 0.0
    for i in checks:
        rows = grid_rows_per_se[i]
        direct = M[rows, :].sum(axis=0).astype(np.float64) if rows.size else np.zeros(len(samples))
        max_err = max(max_err, float(np.abs(direct - se_signal[i, :]).max()))
    empty = sum(1 for r in grid_rows_per_se if r.size == 0)
    print(f"[aggregate] wrote {catalog_bed}, {sig_tsv}, {pres_tsv}", file=sys.stderr)
    print(f"[aggregate] VALIDATION: recompute max|err|={max_err:.3g} over {len(checks)} SEs; "
          f"{empty} SE loci map to 0 grid rows", file=sys.stderr)


if __name__ == "__main__":
    main()
