"""Super-enhancer cutoff — a faithful port of ROSE2_callSuper.R's calculate_cutoff.

ROSE2 sorts the signal ascending, floors negatives to 0, and slides a diagonal of slope
(max-min)/n, finding the tangent point x that minimises the number of points below the line
(numPts_below_line). The cutoff is the sorted signal at floor(x); super-enhancers are the regions with
signal STRICTLY GREATER than that cutoff. See ../DESIGN.md §7.

R indexes 1..n and `myVector[x]` truncates x toward zero, so the objective is a step function of
floor(x); we reproduce that exactly and minimise it with a bounded Brent search (R's `optimize`).
Validated by feeding byte-identical signal tables to this and ROSE2_callSuper.R on the 13 pilot samples.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize_scalar


def _num_pts_below_line(x, v, slope):
    """ROSE numPts_below_line: yPt = v[floor(x)] (1-indexed); count v_i <= i*slope + (yPt - slope*x)."""
    n = v.shape[0]
    xi = int(np.floor(x))
    xi = 1 if xi < 1 else (n if xi > n else xi)
    y_pt = v[xi - 1]
    b = y_pt - slope * x
    x_pts = np.arange(1, n + 1)
    return int(np.count_nonzero(v <= x_pts * slope + b))


def calculate_cutoff(signal):
    """Return ROSE's absolute super-enhancer cutoff for a 1-D signal vector."""
    v = np.sort(np.asarray(signal, dtype=float))
    v[v < 0] = 0.0
    n = v.shape[0]
    if n == 0:
        return 0.0
    slope = (v[-1] - v[0]) / n
    if slope <= 0:                     # flat vector -> no separation; ROSE would put cutoff at the value
        return float(v[-1])
    res = minimize_scalar(
        _num_pts_below_line, bounds=(1, n), args=(v, slope),
        method="bounded", options={"xatol": 1e-4},
    )
    x_pt = int(np.floor(res.x))
    x_pt = 1 if x_pt < 1 else (n if x_pt > n else x_pt)
    return float(v[x_pt - 1])


def call_super(signal):
    """Return (cutoff, is_super) where is_super[i] = 1 iff signal[i] > cutoff (ROSE `>` semantics)."""
    sig = np.asarray(signal, dtype=float)
    cutoff = calculate_cutoff(sig)
    is_super = (sig > cutoff).astype(int)
    return cutoff, is_super
