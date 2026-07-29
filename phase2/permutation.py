#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Label-permutation FDR for SE x group specificity — the fix for gotcha 70.

`specificity.py` builds its null by fitting a NORMAL to each group's own JSD column. JSD-to-one-hot is
bounded at ln(2) and heavily left-skewed, so that normal understates the left tail by ~32x at z=-3 and
~770x at z=-4; the resulting p-values are far too small and ~7% of the catalog gets called "specific"
(99% of SEs specific to >=1 subtype, median SE specific to 6 subtypes of 75).

This replaces the assumed null with a MEASURED one. Shuffle which cell line carries which group label,
rebuild the per-group means, recompute JSD, and repeat. The resulting scores are what this panel produces
when there is no true relationship between lines and labels — so they inherit the real distribution shape,
the ln(2) bound, the group-size structure and the SE-to-SE signal heterogeneity, none of which the normal
captures.

Two properties make the permutation the right null here rather than merely a stricter one:

  * GROUP SIZES ARE PRESERVED EXACTLY. Shuffling the label vector keeps each label's count, so permuted
    "Breast" still averages 32 lines. This matters because JSD depends strongly on group size — a pooled
    or size-agnostic null would be conservative for big groups and permissive for small ones.
  * IT IS PER-GROUP BY CONSTRUCTION, so the per-group null of gotcha 27 is retained. Only the SHAPE
    assumption is dropped. The BH scope stays global (gotcha 27 established that is what makes counts
    comparable), and `_bh_log10` from `specificity.py` is reused so the step-up and the NaN hardening of
    gotcha 64 are shared.

**DEGENERATE AT `line` LEVEL.** When every group is a single cell line, `build_rep_matrix` returns the line
matrix itself and permuting labels only RENAMES columns — the null equals the observed distribution and
every p-value is uniform by construction. `permutation_fdr` refuses that level rather than returning
meaningless numbers.

Resolution: with B permutations over n SEs the smallest attainable p is 1/(B*n+1) — for B=50, n=42,943
that is 4.7e-7, comfortably below the BH bar of alpha*k/N at the few-hundred-calls-per-group scale we
expect. p-values above `keep_frac` are censored to 1.0, which is conservative (inflating large p can only
remove calls, never add them) and keeps memory bounded.
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
import pandas as pd

_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
sys.path.insert(0, _ROOT)
from secacts_env import DATAROOT                      # noqa: E402  (local paths live in .env)
sys.path.insert(0, os.path.join(DATAROOT, "002.AI_projects", "pyCaCTS"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pycacts.score import cacts_score_matrix          # noqa: E402
from pycacts.grouping import build_rep_matrix         # noqa: E402
from specificity import _bh_log10                     # noqa: E402


def permuted_nulls(lines, model, level, n_perm=50, seed=0, keep_frac=0.05, min_group_n=1, verbose=True):
    """{group: sorted np.float32 array of null JSD values} plus the censoring cutoff per group.

    lines : SE x ModelID matrix (the same one scored for the observed JSD)
    """
    if level == "line":
        raise ValueError("permutation is degenerate at 'line' level — permuting labels only renames "
                         "single-line groups, so the null equals the observed distribution")
    rng = np.random.default_rng(seed)
    present = [c for c in lines.columns if c in model.index]
    base = model.loc[present, level].values.copy()
    m2 = model.copy()
    n_se = lines.shape[0]
    # Only the LEFT TAIL is ever consulted (BH thresholds here live near 1e-4 and below), so keep the K
    # smallest draws and prune as we go. Holding every draw would be n_perm x n_se x n_groups floats —
    # 12.9 GB at n_perm=1000 on the subtype level — for information that is never read.
    K = max(int(math.ceil(keep_frac * n_perm * n_se)), 1000)
    kept, buf = {}, {}
    flush_every = max(1, min(25, n_perm))

    def _flush():
        for g, chunks in buf.items():
            parts = chunks if g not in kept else [kept[g]] + chunks
            merged = np.concatenate(parts)
            merged.sort(kind="stable")
            kept[g] = merged[:K].copy()
        buf.clear()

    for b in range(n_perm):
        perm = base.copy()
        rng.shuffle(perm)
        m2.loc[present, level] = perm
        rep, _ = build_rep_matrix(lines, m2, level, min_group_n=min_group_n)
        rep.columns = [str(c) for c in rep.columns]
        j = cacts_score_matrix(rep)
        for g in j.columns:
            buf.setdefault(g, []).append(j[g].values.astype(np.float32))
        if (b + 1) % flush_every == 0:
            _flush()
        if verbose and (b + 1) % 100 == 0:
            print(f"[perm]   {b + 1}/{n_perm} permutations", file=sys.stderr, flush=True)
    _flush()
    m_total = n_perm * n_se
    cutoff = {g: (float(v[-1]) if v.size >= K else float("inf")) for g, v in kept.items()}
    return kept, cutoff, m_total


def permutation_fdr(jsd, lines, model, level, n_perm=50, seed=0, keep_frac=0.05, scope="global",
                    verbose=True):
    """log10(FDR) for an SE x group JSD DataFrame, with the null MEASURED by label permutation.

    Lower JSD = more specific, so the p-value is the left-tail fraction of permuted scores at or below
    the observed one, with the standard +1/+1 correction (a permutation p is never 0).
    """
    nulls, cutoff, m_total = permuted_nulls(lines, model, level, n_perm=n_perm, seed=seed,
                                            keep_frac=keep_frac, verbose=verbose)
    lnp = {}
    for g in jsd.columns:
        x = jsd[g].values.astype(float)
        v = nulls.get(g)
        if v is None:
            lnp[g] = np.zeros_like(x)
            continue
        cnt = np.searchsorted(v, x, side="right").astype(float)
        p = (cnt + 1.0) / (m_total + 1.0)
        p = np.where(x > cutoff[g], 1.0, p)                            # censored region -> conservative
        lnp[g] = np.log(np.clip(p, 1e-300, 1.0))
    groups = list(jsd.columns)
    if scope == "global":
        flat = _bh_log10(np.concatenate([lnp[g] for g in groups])).reshape(len(groups), -1)
        out = {g: flat[i] for i, g in enumerate(groups)}
    else:
        out = {g: _bh_log10(lnp[g]) for g in groups}
    if verbose:
        n_tests = jsd.size
        print(f"[perm] {m_total:,} null draws/group; smallest attainable p = {1.0 / (m_total + 1):.2e}; "
              f"BH bar at k=1 is {0.10 / n_tests:.2e} over {n_tests:,} tests"
              f"{'  <-- RESOLUTION-LIMITED' if 1.0 / (m_total + 1) > 0.10 / n_tests else ''}",
              file=sys.stderr, flush=True)
    return pd.DataFrame(out, index=jsd.index)[groups]


def calibration_check(lines, model, level, n_perm=20, seed=1, fdr=0.10, verbose=True):
    """THE test that matters: run the whole procedure on data where the null is TRUE (labels shuffled).
    A calibrated method should call ~nothing. Returns (n_calls, n_tests)."""
    rng = np.random.default_rng(seed)
    present = [c for c in lines.columns if c in model.index]
    perm = model.loc[present, level].values.copy()
    rng.shuffle(perm)
    m2 = model.copy()
    m2.loc[present, level] = perm
    rep, _ = build_rep_matrix(lines, m2, level, min_group_n=1)
    rep.columns = [str(c) for c in rep.columns]
    jsd_fake = cacts_score_matrix(rep)
    F = permutation_fdr(jsd_fake, lines, m2, level, n_perm=n_perm, seed=seed + 1, verbose=verbose)
    calls = int((np.power(10.0, F.values) <= fdr).sum())
    return calls, int(F.size)
