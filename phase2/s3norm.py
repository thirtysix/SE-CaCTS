#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""S3norm — cross-study normalization of the Phase-2 fine grid matrix (DESIGN.md §"normalization").

Replaces quantile normalization, which *manufactures false positives* by forcing every sample onto an
identical distribution (S3norm benchmark: quantile-specific CTCF peaks matched the motif 9.2% of the time
vs 80.5% for consensus peaks — Xiang et al., NAR 2020, e43). S3norm instead fits a **two-parameter
monotone power transform** per sample

    f(x) = A * (x + pseudo)^B

choosing (A, B) so that the sample matches the reference in BOTH strata simultaneously — the enriched
("peak") stratum and the background stratum. Two anchors, two unknowns. Because the transform has only
2 degrees of freedom (vs quantile norm's one-per-row), it co-normalizes sequencing depth AND
signal-to-noise while leaving the *shape* of the signal distribution — hence the biology — intact.

Solving. With P = common-enriched rows and B_ = common-background rows, the A cancels in the ratio:

    mean_P((x+p)^B) / mean_B_((x+p)^B)  =  mean_P(ref) / mean_B_(ref)

a single monotone equation in B (Brent root-find), after which A follows in closed form. Values are
rescaled by their P-mean before exponentiation so the powers stay near 1.0 and cannot overflow.

Anchor strata (`--anchor`, default `quantile`). The obvious choice — MACS2 peak calls, as in the paper —
is NOT usable on this compendium: at ChIP-Atlas Q<1e-20 the pilot's per-sample peak counts span 53 to
23,528, so the weakest samples' common-peak sets collapse to a handful of rows and a 2-parameter fit on
them is noise. Worse, the peak-call threshold is itself SNR-dependent, so anchoring on it is circular —
exactly the cross-lab variation we are trying to remove. `quantile` mode therefore defines the strata by
each sample's OWN signal quantiles (top `--hi-frac`; background band `--lo-band`), which is SNR-relative
by construction and yields 1.1k-3.6k common-enriched rows per pair on the pilot. `macs2` mode is kept for
faithfulness/comparison and warns when a stratum is too small to fit.

Where this belongs in the pipeline. The transform is NONLINEAR, so Σf(xᵢ) ≠ f(Σxᵢ): it must be applied at
GRID-ROW resolution and the SE-level signal re-summed afterwards, not applied to SE sums. `aggregate.py
--norm s3norm` does exactly that. The SE *catalog* is deliberately left alone — SE calling is a
within-sample threshold (the ROSE tangent) and is not a cross-sample comparison.

  python s3norm.py --signal-dir <dir of SRX.signal.tsv> --out <prefix>      # normalize + diagnostics
  python s3norm.py --npy fine.npy --srx fine_srx.json --out <prefix>        # from a cached matrix

atac_hdac env (numpy + scipy + pandas).
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import brentq

B_BRACKET = (0.02, 25.0)      # plausible range for the power exponent
MIN_ANCHOR = 200              # refuse to fit a 2-parameter transform on fewer rows than this


# --------------------------------------------------------------------------- strata


def strata_quantile(M, hi_frac=0.01, lo_band=(0.20, 0.50)):
    """Boolean (n x k) masks (enriched, background) from each sample's OWN signal quantiles."""
    n, k = M.shape
    hi = np.zeros((n, k), dtype=bool)
    lo = np.zeros((n, k), dtype=bool)
    for j in range(k):
        col = M[:, j]
        hi[:, j] = col >= np.quantile(col, 1.0 - hi_frac)
        a, b = np.quantile(col, lo_band)
        lo[:, j] = (col >= a) & (col <= b)
    return hi, lo


def strata_macs2(M, peak_mask):
    """Strata from real MACS2 peak calls: enriched = called peak, background = everything else."""
    return peak_mask, ~peak_mask


# --------------------------------------------------------------------------- fit


def _fit_one(x, ref, hi_t, lo_t, hi_r, lo_r, pseudo):
    """Fit (A, B) mapping sample x onto ref. Returns (A, B, n_hi, n_lo)."""
    P = hi_t & hi_r          # common-enriched
    Bg = lo_t & lo_r         # common-background
    n_hi, n_lo = int(P.sum()), int(Bg.sum())
    if n_hi < MIN_ANCHOR or n_lo < MIN_ANCHOR:
        return np.nan, np.nan, n_hi, n_lo

    xp = x.astype(np.float64) + pseudo
    rp = ref.astype(np.float64) + pseudo
    # rescale by the enriched-stratum mean so (x/m)^B stays O(1) regardless of B
    m = xp[P].mean()
    u_hi, u_lo = xp[P] / m, xp[Bg] / m
    target = np.log(rp[P].mean() / rp[Bg].mean())

    def g(b):
        return np.log((u_hi ** b).mean() / (u_lo ** b).mean()) - target

    lo_b, hi_b = B_BRACKET
    g_lo, g_hi = g(lo_b), g(hi_b)
    if g_lo > 0 or g_hi < 0:                       # target ratio outside what any B can produce
        b = lo_b if g_lo > 0 else hi_b
    else:
        b = brentq(g, lo_b, hi_b, xtol=1e-8)
    A = rp[P].mean() / ((u_hi ** b).mean() * (m ** b))
    return float(A), float(b), n_hi, n_lo


def pick_reference(M, mode="medoid", srx=None):
    """Reference column index. 'medoid' = the real sample closest to the average log-quantile profile."""
    if mode == "mean":
        return -1                                   # sentinel: synthetic pointwise-mean reference
    if srx is not None and mode in srx:
        return srx.index(mode)
    if mode != "medoid":
        raise SystemExit(f"--ref: expected 'medoid', 'mean', or a sample id; got {mode!r}")
    qs = np.linspace(0.01, 0.99, 99)
    prof = np.column_stack([np.quantile(np.log1p(M[:, j]), qs) for j in range(M.shape[1])])
    d = np.linalg.norm(prof - prof.mean(axis=1, keepdims=True), axis=0)
    return int(np.argmin(d))


def s3norm_matrix(M, ref="medoid", srx=None, hi_frac=0.01, lo_band=(0.20, 0.50), pseudo=1.0,
                  peak_mask=None, anchor="quantile", verbose=True):
    """Normalize an (n_rows x n_samples) signal matrix. Returns (normalized float32 matrix, params DataFrame)."""
    M = np.asarray(M)
    n, k = M.shape
    srx = list(srx) if srx is not None else [f"s{j}" for j in range(k)]

    ri = pick_reference(M, ref, srx)
    ref_col = M.mean(axis=1) if ri < 0 else M[:, ri]
    ref_name = "<mean>" if ri < 0 else srx[ri]

    if anchor == "macs2":
        if peak_mask is None:
            raise SystemExit("--anchor macs2 requires per-sample peak calls (--peak-dir)")
        hi, lo = strata_macs2(M, peak_mask)
        hi_r = ref_col >= np.quantile(ref_col, 1 - hi_frac) if ri < 0 else hi[:, ri]
        lo_r = ~hi_r if ri < 0 else lo[:, ri]
    else:
        hi, lo = strata_quantile(M, hi_frac, lo_band)
        if ri < 0:
            hi_r = ref_col >= np.quantile(ref_col, 1 - hi_frac)
            a, b = np.quantile(ref_col, lo_band)
            lo_r = (ref_col >= a) & (ref_col <= b)
        else:
            hi_r, lo_r = hi[:, ri], lo[:, ri]

    out = np.empty_like(M, dtype=np.float32)
    rows = []
    for j in range(k):
        A, B, n_hi, n_lo = _fit_one(M[:, j], ref_col, hi[:, j], lo[:, j], hi_r, lo_r, pseudo)
        if np.isnan(A):
            A, B = 1.0, 1.0                          # identity; flagged in the params table
            print(f"[s3norm] WARN {srx[j]}: anchors too small (enriched={n_hi}, bg={n_lo} < {MIN_ANCHOR}) "
                  f"-> left UNNORMALIZED", file=sys.stderr)
        out[:, j] = A * (M[:, j].astype(np.float64) + pseudo) ** B
        rows.append(dict(sample=srx[j], A=A, B=B, n_enriched=n_hi, n_background=n_lo,
                         is_ref=(j == ri), fitted=not np.isnan(B) and (n_hi >= MIN_ANCHOR)))
    params = pd.DataFrame(rows)

    if verbose:
        print(f"[s3norm] {n} rows x {k} samples | ref={ref_name} | anchor={anchor} "
              f"hi_frac={hi_frac} lo_band={lo_band} pseudo={pseudo}", file=sys.stderr)
        for _, r in params.iterrows():
            tag = "  <- REF" if r["is_ref"] else ""
            print(f"[s3norm]   {r['sample']:<14} A={r['A']:<10.4g} B={r['B']:<7.4f} "
                  f"anchors={int(r['n_enriched'])}/{int(r['n_background'])}{tag}", file=sys.stderr)
    return out, params


# --------------------------------------------------------------------------- io / cli


def read_signal_dir(sig_dir):
    """Read <SRX>.signal.tsv files (CHROM START STOP SIGNAL) into (matrix, srx list)."""
    files = sorted(glob.glob(os.path.join(sig_dir, "*.signal.tsv")))
    if not files:
        raise SystemExit(f"no *.signal.tsv in {sig_dir}")
    srx = [os.path.basename(f)[: -len(".signal.tsv")] for f in files]
    cols = [np.loadtxt(f, delimiter="\t", skiprows=1, usecols=3, dtype=np.float32) for f in files]
    ns = {c.size for c in cols}
    if len(ns) != 1:
        raise SystemExit(f"signal files disagree on row count: {ns}")
    return np.column_stack(cols), srx


def read_peak_masks(peak_dir, grid_bed, srx):
    """Boolean (n_grid x k) mask: does this sample's MACS2 peak set cover each grid row (by midpoint)?"""
    import bisect
    grid = {}
    order = []
    with open(grid_bed) as fh:
        for i, line in enumerate(fh):
            c, s, e = line.split("\t")[:3]
            order.append((c, (int(s) + int(e)) // 2))
    by_chrom = {}
    for i, (c, mid) in enumerate(order):
        by_chrom.setdefault(c, []).append((mid, i))
    mask = np.zeros((len(order), len(srx)), dtype=bool)
    for j, s in enumerate(srx):
        path = os.path.join(peak_dir, f"{s}.20.bed")
        if not os.path.exists(path):
            print(f"[s3norm] WARN no peaks for {s}", file=sys.stderr)
            continue
        iv = {}
        with open(path) as fh:
            for line in fh:
                f = line.split("\t")
                iv.setdefault(f[0], []).append((int(f[1]), int(f[2])))
        for c, lst in iv.items():
            lst.sort()
            starts = [a for a, _ in lst]
            for mid, i in by_chrom.get(c, ()):
                p = bisect.bisect_right(starts, mid) - 1
                if p >= 0 and lst[p][1] > mid:
                    mask[i, j] = True
        grid[s] = int(mask[:, j].sum())
    return mask


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--signal-dir", help="dir of <SRX>.signal.tsv (the Phase-2 fine matrix columns)")
    src.add_argument("--npy", help="cached (n_grid x k) float32 matrix")
    ap.add_argument("--srx", help="json list of sample ids (required with --npy)")
    ap.add_argument("--out", required=True, help="output prefix (writes .npy + .params.tsv)")
    ap.add_argument("--ref", default="medoid", help="'medoid' (default), 'mean', or an explicit sample id")
    ap.add_argument("--anchor", default="quantile", choices=["quantile", "macs2"])
    ap.add_argument("--hi-frac", type=float, default=0.01, help="top fraction = enriched stratum")
    ap.add_argument("--lo-band", type=float, nargs=2, default=[0.20, 0.50], help="background percentile band")
    ap.add_argument("--pseudo", type=float, default=1.0, help="pseudocount added before the power transform")
    ap.add_argument("--peak-dir", help="dir of <SRX>.20.bed (required for --anchor macs2)")
    ap.add_argument("--grid", help="grid BED (required for --anchor macs2)")
    a = ap.parse_args()

    if a.signal_dir:
        M, srx = read_signal_dir(a.signal_dir)
    else:
        if not a.srx:
            raise SystemExit("--npy requires --srx")
        M, srx = np.load(a.npy), json.load(open(a.srx))

    peak_mask = read_peak_masks(a.peak_dir, a.grid, srx) if a.anchor == "macs2" else None
    out, params = s3norm_matrix(M, ref=a.ref, srx=srx, hi_frac=a.hi_frac, lo_band=tuple(a.lo_band),
                                pseudo=a.pseudo, peak_mask=peak_mask, anchor=a.anchor)
    np.save(a.out + ".npy", out)
    params.to_csv(a.out + ".params.tsv", sep="\t", index=False)
    json.dump(srx, open(a.out + ".srx.json", "w"))
    print(f"[s3norm] wrote {a.out}.npy + .params.tsv + .srx.json", file=sys.stderr)


if __name__ == "__main__":
    main()
