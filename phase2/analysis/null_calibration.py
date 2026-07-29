#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Is the empirical-null FDR actually calibrated? (Answer, 2026-07-22: NO — it massively over-calls.)

`specificity.py` builds a per-group null by fitting a NORMAL to the group's own JSD column (median +
one-sided MAD), converting to left-tail p-values, and running BH. That is only an FDR if the normal
describes the distribution of JSD scores arising WITHOUT real specificity. It does not:

  * JSD-to-one-hot is bounded above by ln(2) = 0.6931 and the observed columns are heavily left-skewed,
    so the fitted normal is a poor description of the left tail.
  * Measured on the real atlas, the observed left tail runs ~32x heavier than the fitted normal at z=-3
    and ~770x at z=-4. The p-values are therefore far too small.
  * BH is correct GIVEN the p-values, so a handful of astronomically small ones drag the step-up
    threshold up to p ~ 8.5e-3 (z ~ -2.4) — a very weak per-test bar — and ~3,000 SEs per group pass.

Consequences, and what survives:
  BROKEN   — any COUNT of "specific SEs". 99% of the catalog is called specific to >=1 subtype and the
             median SE is called specific to 6 subtypes of 75. The word stops meaning anything.
  INTACT   — the RANKINGS (gotcha 28: rankings are invariant to null/BH choice), so MECOM #1, MCF7
             ESR1 #1, EMX2 #4 and the CN ablation (a comparison between two arms) are unaffected.

Two independent lines converge on a few hundred SEs per group rather than a few thousand: FDR<=1e-12
yields a median of ~256/group, and Phase-6 concordance holds above 20% only to about SE rank 200
(29.1% at rank 1-15 -> 6.2% beyond rank 1000).

The principled fix is a LABEL-PERMUTATION null (already on the roadmap): permute cell lines across group
labels, recompute JSD, and take the null from the permuted scores — which inherits the true shape instead
of assuming normality. Degenerate at `line` level, where every group is a single line.

  ~/miniconda3/envs/atac_hdac/bin/python phase2/analysis/null_calibration.py
"""
from __future__ import annotations

import argparse
import math
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))
SECACTS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(SECACTS, "phase2"))

from specificity import null_params, fdr_matrix       # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jsd-pickle", required=True, help="a saved SE x group JSD DataFrame")
    ap.add_argument("--fdr", type=float, default=0.10)
    a = ap.parse_args()
    jsd = pd.read_pickle(a.jsd_pickle)
    F = fdr_matrix(jsd)
    print(f"JSD matrix {jsd.shape[0]:,} SEs x {jsd.shape[1]} groups; ln(2) ceiling = {math.log(2):.4f}\n")

    print("=== 1. tail heaviness vs the fitted normal (pooled over groups) ===")
    print(f"  {'z cut':>7} {'observed':>12} {'normal pred':>13} {'obs/pred':>10}")
    Z = np.concatenate([(jsd[g].values - null_params(jsd[g].values)[0]) / null_params(jsd[g].values)[1]
                        for g in jsd.columns])
    for c in [-1.5, -2, -2.5, -3, -4, -5]:
        obs = int((Z < c).sum())
        pred = norm.cdf(c) * Z.size
        print(f"  {c:>7.1f} {obs:>12,} {pred:>13,.0f} {obs / max(pred, 1e-9):>9.1f}x")

    print(f"\n=== 2. the per-test bar that FDR<={a.fdr} actually applies ===")
    called = np.power(10.0, F.values) <= a.fdr
    zc = Z.reshape(jsd.shape[1], -1).T[called] if called.any() else np.array([])
    if zc.size:
        print(f"  calls={called.sum():,}/{called.size:,}   least-extreme call z={zc.max():+.2f}  "
              f"p={norm.cdf(zc.max()):.2e}")

    print("\n=== 3. calls vs threshold ===")
    for thr in [a.fdr, 0.01, 1e-3, 1e-6, 1e-12]:
        m = np.power(10.0, F.values) <= thr
        print(f"  FDR<={thr:<8g} total={int(m.sum()):>9,}   median/group={int(np.median(m.sum(axis=0))):>6,}"
              f"   ({100 * np.median(m.sum(axis=0)) / jsd.shape[0]:.1f}% of catalog)")


if __name__ == "__main__":
    main()
