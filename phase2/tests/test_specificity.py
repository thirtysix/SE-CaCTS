#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests for phase2/specificity.py.

The load-bearing one is BACKWARD COMPATIBILITY: fdr_matrix(null='pergroup', scope='pergroup') must reproduce
`pyCaCTS.empirical_fdr` column-by-column. If that drifts, the new default silently changes more than intended.

  ~/miniconda3/envs/atac_hdac/bin/python phase2/tests/test_specificity.py
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, _ROOT)
from secacts_env import pycacts_path                                 # noqa: E402  (paths live in .env)
sys.path.insert(0, pycacts_path())

from specificity import fdr_matrix, null_params, specific_counts     # noqa: E402
from pycacts.stats import empirical_fdr_log10                        # noqa: E402


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}{'  — ' + detail if detail else ''}")
    return bool(ok)


def main():
    rng = np.random.default_rng(0)
    n, k = 3000, 6
    # groups with deliberately different score dispersions — the situation that breaks per-group BH
    jsd = pd.DataFrame(
        {f"g{j}": np.clip(rng.normal(0.53, 0.06 + 0.03 * j, n), 0.01, 0.99) for j in range(k)},
        index=[f"USE_{i}" for i in range(n)])
    for j in range(k):                                    # give each group a real left tail
        jsd.iloc[rng.choice(n, 40, replace=False), j] *= 0.25

    ok = []
    mine = fdr_matrix(jsd, null="pergroup", scope="pergroup")
    theirs = pd.DataFrame({g: empirical_fdr_log10(jsd[g]) for g in jsd.columns})
    ok.append(check("matches pyCaCTS.empirical_fdr (pergroup/pergroup)",
                    np.allclose(mine.values, theirs.values, atol=1e-9),
                    f"max|diff|={np.abs(mine.values - theirs.values).max():.2e}"))

    mu, sig = null_params(jsd["g0"].values)
    ok.append(check("null uses the right (non-specific) side", sig > 0 and 0.4 < mu < 0.7,
                    f"mu={mu:.3f} sig={sig:.3f}"))

    # rankings must be invariant to null/scope — only the threshold moves. BH's step-up creates plateaus of
    # equal FDR, so assert MONOTONICITY along each group's JSD ranking rather than exact index order.
    inv = []
    for null in ("pergroup", "global"):
        for scope in ("pergroup", "global"):
            f = fdr_matrix(jsd, null=null, scope=scope)
            inv.append(all(np.all(np.diff(f.loc[jsd[g].sort_values().index, g].values) >= -1e-12)
                           for g in jsd.columns))
    ok.append(check("FDR is monotone along the JSD ranking, for every null/scope", all(inv)))

    cnt_pg = specific_counts(jsd, null="pergroup", scope="pergroup")
    cnt_gl = specific_counts(jsd, null="pergroup", scope="global")
    # the pathology this module exists to fix: per-group BH starves groups into ZERO calls, because each
    # group spends an independent testing budget. Global BH shares one budget and recovers them.
    z_pg, z_gl = int((cnt_pg == 0).sum()), int((cnt_gl == 0).sum())
    ok.append(check("global BH yields fewer zero-call groups", z_gl <= z_pg,
                    f"zero-call groups {z_pg} -> {z_gl}  (per-group {cnt_pg.tolist()} -> global {cnt_gl.tolist()})"))

    cv = lambda c: float(c.std() / max(c.mean(), 1e-9))                # noqa: E731
    ok.append(check("global BH does not worsen count dispersion", cv(cnt_gl) <= cv(cnt_pg) + 1e-9,
                    f"CV {cv(cnt_pg):.2f} -> {cv(cnt_gl):.2f}"))

    ok.append(check("global BH is not simply more permissive overall",
                    cnt_gl.sum() <= cnt_pg.sum() * 1.5,
                    f"total calls {int(cnt_pg.sum())} -> {int(cnt_gl.sum())}"))

    f = fdr_matrix(jsd)
    ok.append(check("default is pergroup null + global BH",
                    np.allclose(f.values, fdr_matrix(jsd, null="pergroup", scope="global").values)))
    ok.append(check("log10 FDR never exceeds 0", float(f.values.max()) <= 0.0))

    # ---- NaN handling. This is a REGRESSION GUARD for a silent, total failure found 2026-07-22 on the real
    # atlas: scoring-time CN correction pushed a few near-zero cells negative, pycacts turned those into NaN
    # JSD, `np.median` made the whole group's null NaN, and `min(0.0, nan) == 0.0` in Python meant BH
    # reported FDR = 1.0 for EVERY test at 3 of 4 hierarchy levels. It looked like "no significant results".
    jn = jsd.copy()
    jn.iloc[0, 0] = np.nan                                        # ONE NaN, in one group
    mu_n, sig_n = null_params(jn["g0"].values)
    ok.append(check("null_params ignores NaN instead of returning NaN",
                    np.isfinite(mu_n) and np.isfinite(sig_n) and sig_n > 0,
                    f"mu={mu_n:.3f} sig={sig_n:.3f}"))

    fn = fdr_matrix(jn)
    finite_col = np.isfinite(fn["g0"].values)
    ok.append(check("one NaN does not void its whole group column",
                    finite_col.sum() == n - 1 and not np.isfinite(fn.iloc[0, 0]),
                    f"finite entries in g0: {int(finite_col.sum())}/{n}"))
    # the calls in the untouched groups must be unchanged, and g0 must keep essentially all of its own
    ok.append(check("a NaN cell does not suppress calls elsewhere",
                    int((fn.iloc[:, 1:] <= math.log10(0.10)).values.sum())
                    == int((f.iloc[:, 1:] <= math.log10(0.10)).values.sum())
                    and int((fn["g0"] <= math.log10(0.10)).sum()) >= int((f["g0"] <= math.log10(0.10)).sum()) - 1,
                    f"g0 calls {int((f['g0'] <= math.log10(0.10)).sum())} -> "
                    f"{int((fn['g0'] <= math.log10(0.10)).sum())}"))

    # the exact shape of the old bug: an ALL-NaN column must read as NaN, never as FDR = 1.0
    jz = jsd.copy()
    jz["g0"] = np.nan
    fz = fdr_matrix(jz)
    ok.append(check("an all-NaN group returns NaN, not FDR=1.0",
                    bool(np.all(~np.isfinite(fz["g0"].values)))
                    and int((fz.iloc[:, 1:] <= math.log10(0.10)).values.sum()) > 0,
                    f"surviving calls in the other {k - 1} groups: "
                    f"{int((fz.iloc[:, 1:] <= math.log10(0.10)).values.sum())}"))

    print(f"\n{sum(ok)}/{len(ok)} passed")
    return 0 if all(ok) else 1


if __name__ == "__main__":
    sys.exit(main())
