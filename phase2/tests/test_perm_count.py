#!/usr/bin/env python3
"""permutation_fdr_count must reproduce permutation_fdr: same calls, and the same p-values wherever the slow
path does not censor. Synthetic panel with planted group-specific loci, small enough to run in seconds.
    python3 phase2/tests/test_perm_count.py
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, _ROOT)
from secacts_env import pycacts_path                                 # noqa: E402
sys.path.insert(0, pycacts_path())
from pycacts.grouping import build_rep_matrix                        # noqa: E402
from pycacts.score import cacts_score_matrix                         # noqa: E402
from permutation import permutation_fdr, permutation_fdr_count      # noqa: E402


def panel(n_se=3000, n_line=90, n_grp=18, seed=3):
    rng = np.random.default_rng(seed)
    lines = [f"ACH-{i:06d}" for i in range(n_line)]
    sizes = rng.integers(1, 9, n_grp); sizes = np.maximum(1, (sizes / sizes.sum() * n_line).astype(int))
    lab = np.concatenate([[f"g{k}"] * s for k, s in enumerate(sizes)])[:n_line]
    lab = np.concatenate([lab, [f"g{n_grp - 1}"] * (n_line - len(lab))])
    lab = lab.astype(object); lab[rng.choice(n_line, 4, replace=False)] = np.nan       # unlabelled lines move too
    X = rng.gamma(2.0, 40.0, (n_se, n_line))
    for k in range(n_grp):                                           # planted specific loci
        idx = np.where(lab == f"g{k}")[0]
        X[np.ix_(rng.choice(n_se, 25, replace=False), idx)] *= 6.0
    lines_df = pd.DataFrame(X, index=[f"USE_{i}" for i in range(n_se)], columns=lines)
    model = pd.DataFrame({"OncotreeSubtype": lab}, index=lines)
    return lines_df, model


def main():
    lines, model = panel()
    lev = "OncotreeSubtype"
    rep, _ = build_rep_matrix(lines, model, lev, min_group_n=1)
    rep.columns = [str(c) for c in rep.columns]
    jsd = cacts_score_matrix(rep)
    kw = dict(n_perm=200, seed=5, keep_frac=0.05, verbose=False)
    slow = permutation_fdr(jsd, lines, model, lev, **kw)
    for workers in (1, 3):
        fast = permutation_fdr_count(jsd, lines, model, lev, n_workers=workers, **kw)
        cs, cf = slow.values <= np.log10(0.10), fast.values <= np.log10(0.10)
        assert (cs == cf).all(), f"calls differ with {workers} workers: slow {cs.sum()}, fast {cf.sum()}"
        both = (slow.values < np.log10(0.10)) | (fast.values < np.log10(0.10))
        assert np.allclose(slow.values[both], fast.values[both], rtol=0, atol=1e-12), "FDR differs where it matters"
        print(f"ok: {workers} worker(s), {int(cs.sum())} calls identical, {int(both.sum())} FDRs below 0.10 identical")


if __name__ == "__main__":
    main()


def units_main():
    """permutation_fdr_units_count vs permutation_fdr_units: q25 aggregation over units, strata, excluded groups."""
    from permutation import permutation_fdr_units, permutation_fdr_units_count, rep_agg
    rng = np.random.default_rng(11)
    n_se, n_grp = 2500, 40
    sizes = rng.integers(1, 6, n_grp)                                  # units (studies) per group (line)
    labels = np.concatenate([[f"L{k}"] * s for k, s in enumerate(sizes)])
    groups = [f"L{k}" for k in range(n_grp)]
    X = rng.gamma(2.0, 40.0, (n_se, len(labels)))
    for k in range(n_grp):
        X[np.ix_(rng.choice(n_se, 20, replace=False), np.where(labels == f"L{k}")[0])] *= 6.0
    lineage = np.array([f"lin{k % 6}" for k in range(n_grp)])
    strata = np.array([lineage[int(l[1:])] for l in labels])
    codes = np.array([groups.index(l) for l in labels])
    jsd = cacts_score_matrix(pd.DataFrame(rep_agg(X, codes, n_grp, "q25"), columns=groups))
    exclude = {g for g, s in zip(groups, sizes) if s < 2} | {"L3"}
    kw = dict(n_perm=150, seed=7, keep_frac=0.05, verbose=False, strata=strata, exclude=exclude)
    slow = permutation_fdr_units(jsd, X, labels, "q25", **kw)
    for workers in (1, 3):
        fast = permutation_fdr_units_count(jsd, X, labels, "q25", n_workers=workers, **kw)
        cs, cf = slow.values <= np.log10(0.10), fast.values <= np.log10(0.10)
        assert (cs == cf).all(), f"units calls differ with {workers} workers: slow {cs.sum()}, fast {cf.sum()}"
        both = (slow.values < np.log10(0.10)) | (fast.values < np.log10(0.10))
        assert np.allclose(slow.values[both], fast.values[both], rtol=0, atol=1e-12), "units FDR differs"
        assert (fast[sorted(exclude)].values == 0).all(), "excluded groups must have FDR 1"
        print(f"ok (units, strata, {len(exclude)} excluded): {workers} worker(s), {int(cs.sum())} calls identical")


if __name__ == "__main__":
    units_main()
