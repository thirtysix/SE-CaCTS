#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Specificity thresholding: empirical-null FDR over an SE x group JSD matrix.

`pyCaCTS.empirical_fdr` answers "is this SE specific to group g?" one group at a time — it builds a null from
that group's own score column AND runs Benjamini-Hochberg within that column. Those are two separate choices,
and on the pilot the second one is what makes per-group specific-SE COUNTS incomparable: each group spends its
own independent multiple-testing budget, so a group whose scores happen to be tightly distributed converts a
modest left tail into many calls, while a group with a degenerate null can return literally zero.

Measured on the pilot (13 samples, QC-gated S3norm, CN-corrected), max/min specific-SE count across groups:

    null      BH scope   Lineage   Disease/Subtype   line     SOX17     ESR1
    pergroup  pergroup     207.0             248.0   10.3    -1.73    -1.51     <- pyCaCTS default
    pergroup  global        13.8               3.9    4.5    -2.36    -2.62     <- this module's default
    global    pergroup     195.0               3.8   36.0    -0.18 x  -2.79
    global    global         3.1               2.4   12.6    -1.10    -3.64

(SOX17/ESR1 = log10 FDR of the identity-gene SE for Ovary and Breast; <= -1 passes at FDR 0.10.)

**Global BH wins on both axes at once** — it compresses the count spread AND strengthens both known master-TF
SEs, and it removes the pathological zero-call groups. A *global null* (pooling all groups to estimate mu/sig)
is the intuitive fix but is the wrong lever: it is over-conservative for groups whose score distribution is
genuinely tight, and it drops SOX17 below significance entirely. So: **keep the per-group null, share the BH.**

Rankings are invariant to all of this — only which SEs pass threshold changes.

`fdr_matrix(..., null="pergroup", scope="pergroup")` reproduces `pyCaCTS.empirical_fdr` column-by-column
(asserted in tests/test_specificity.py). This module is the natural thing to upstream into pyCaCTS once the
same effect is confirmed on the TF layer.
"""
from __future__ import annotations

import math
import warnings

import numpy as np
import pandas as pd

_erf = np.vectorize(math.erf)
_LN2PI = math.log(2 * math.pi)


def _ln_left_tail(z):
    """ln of the standard-normal left tail Phi(z), stable for very negative z (matches pyCaCTS)."""
    tiny = np.finfo(float).tiny
    central = np.log(np.clip(0.5 * (1.0 + _erf(z / math.sqrt(2))), tiny, None))
    asymp = -0.5 * z * z - np.log(-np.minimum(z, -tiny)) - 0.5 * _LN2PI
    return np.where(z > -6, central, asymp)


def null_params(x):
    """Robust null (median, one-sided MAD from the NON-specific right side) — lower score = more specific.

    NaN-robust: the null is estimated from the FINITE scores only. This matters — `np.median` of a column
    holding a single NaN returns NaN, after which `x >= mu` is all-False, `dev` is empty and `sig` collapses
    to the 1e-9 fallback, so every z in that column becomes NaN. One NaN would otherwise void a whole group.
    """
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return float("nan"), 1e-9
    mu = float(np.median(x))
    dev = np.abs(x[x >= mu] - mu)
    sig = 1.4826 * float(np.median(dev)) if dev.size else 0.0
    return mu, (sig or 1e-9)


def _bh_log10(lnp):
    """BH step-up on natural-log p-values -> log10(FDR), clipped at 0.

    Non-finite p-values are EXCLUDED from the procedure (they are not tests) and come back as NaN rather
    than 0. Both halves of that matter: `n` must count only real tests or the correction is inflated, and
    the old code returned 0 -> FDR=1.0 for them because Python's `min(0.0, nan)` is 0.0, so a NaN column
    silently reported "nothing is significant" instead of failing. Never let a NaN read as a result.
    """
    lnp = np.asarray(lnp, dtype=float)
    out = np.full(lnp.size, np.nan)
    ok = np.flatnonzero(np.isfinite(lnp))
    n = ok.size
    if n == 0:
        return out
    order = ok[np.argsort(lnp[ok], kind="stable")]
    prev = 0.0
    for k in range(n, 0, -1):
        i = order[k - 1]
        prev = min(prev, lnp[i] + math.log(n) - math.log(k))
        out[i] = min(0.0, prev)
    return out / math.log(10)


def fdr_matrix(jsd, null="pergroup", scope="global"):
    """log10(empirical-null FDR) for an SE x group JSD DataFrame (lower JSD = more specific).

    null  : 'pergroup' (default, recommended) estimates mu/sig from each group's own column;
            'global' pools every group's scores into one null — over-conservative for tight groups.
    scope : 'global' (default, recommended) runs BH once across all SE x group tests, so groups share one
            multiple-testing budget and their counts are comparable;
            'pergroup' runs BH within each column — this is pyCaCTS's behaviour.
    """
    if null not in ("pergroup", "global") or scope not in ("pergroup", "global"):
        raise ValueError(f"bad null={null!r} / scope={scope!r}")
    groups = list(jsd.columns)
    n_bad = int((~np.isfinite(jsd.values)).sum())
    if n_bad:
        # Loud on purpose. A non-finite JSD means the upstream score matrix had non-finite or NEGATIVE
        # input (pycacts propagates NaN by design, score.py), which usually means a correction step pushed
        # signal below zero. Those cells are dropped from the BH here, but the caller should fix the source.
        warnings.warn(f"fdr_matrix: {n_bad:,} of {jsd.size:,} JSD values are non-finite "
                      f"({100.0 * n_bad / jsd.size:.3f}%); excluded from BH and returned as NaN. "
                      f"Check for negative values in the representative matrix.", RuntimeWarning, stacklevel=2)
    gmu, gsig = null_params(jsd.values.ravel())

    lnp = {}
    for g in groups:
        x = jsd[g].values.astype(float)
        mu, sig = null_params(x) if null == "pergroup" else (gmu, gsig)
        lnp[g] = _ln_left_tail((x - mu) / sig)

    if scope == "global":
        flat = _bh_log10(np.concatenate([lnp[g] for g in groups])).reshape(len(groups), -1)
        out = {g: flat[i] for i, g in enumerate(groups)}
    else:
        out = {g: _bh_log10(lnp[g]) for g in groups}
    return pd.DataFrame(out, index=jsd.index)[groups]


def specific_counts(jsd, threshold=0.10, **kw):
    """Per-group count of SEs at FDR <= threshold."""
    return (fdr_matrix(jsd, **kw) <= math.log10(threshold)).sum(axis=0)
