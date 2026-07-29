#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Property tests for phase2/s3norm.py — synthetic distortions with KNOWN ground truth.

Each column is built from a shared latent signal (real compendia are biologically correlated, which is what
populates the common-enriched anchor set) and then distorted in a way whose inverse we know analytically:

  depth scaling  x -> c*x      must recover B = 1,     A = 1/c
  SNR compression x -> x^g     must recover B = 1/g    (the exact inverse exponent)
  undistorted                  must recover A = B = 1  (identity)

Plus the invariants: the reference self-normalizes to identity, the transform is rank-preserving (it is
monotone by construction — this is what keeps S3norm from manufacturing the false positives that quantile
normalization does), and post-normalization every column sits at the reference's scale.

  ~/miniconda3/envs/atac_hdac/bin/python phase2/tests/test_s3norm.py
"""
from __future__ import annotations

import os
import sys

import numpy as np
from scipy.stats import spearmanr

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from s3norm import MIN_ANCHOR, s3norm_matrix          # noqa: E402

N = 50_000
GAMMA = 0.7                                            # compression exponent -> expect B = 1/GAMMA
DEPTH = 4.0                                            # depth factor         -> expect A = 1/DEPTH


def build():
    """4 correlated columns: [reference, depth-scaled, SNR-compressed, undistorted]."""
    rng = np.random.default_rng(0)
    latent = np.abs(rng.lognormal(3, 1.5, size=N))
    cols = [latent * np.exp(rng.normal(0, 0.25, N)) for _ in range(4)]
    X = np.column_stack(cols).astype(np.float32)
    X[:, 1] = X[:, 1] * DEPTH
    X[:, 2] = X[:, 2] ** GAMMA
    return X, ["ref", "depth", "compressed", "clean"]


def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}{'  — ' + detail if detail else ''}")
    return bool(ok)


def main():
    X, names = build()
    out, p = s3norm_matrix(X, ref="ref", srx=names, verbose=False)
    p = p.set_index("sample")
    print(p[["A", "B", "n_enriched", "n_background"]].to_string(), "\n")

    ok = []
    ok.append(check("reference self-normalizes to identity",
                    abs(p.loc["ref", "A"] - 1) < 1e-9 and abs(p.loc["ref", "B"] - 1) < 1e-9))
    ok.append(check("depth scaling -> B = 1", abs(p.loc["depth", "B"] - 1) < 0.05,
                    f"B={p.loc['depth','B']:.4f}"))
    ok.append(check(f"depth scaling -> A = 1/{DEPTH:g}", abs(p.loc["depth", "A"] - 1 / DEPTH) < 0.05,
                    f"A={p.loc['depth','A']:.4f} vs {1/DEPTH:.4f}"))
    ok.append(check(f"SNR compression -> B = 1/{GAMMA:g}", abs(p.loc["compressed", "B"] - 1 / GAMMA) < 0.15,
                    f"B={p.loc['compressed','B']:.4f} vs {1/GAMMA:.4f}"))
    ok.append(check("undistorted column -> identity",
                    abs(p.loc["clean", "A"] - 1) < 0.05 and abs(p.loc["clean", "B"] - 1) < 0.05,
                    f"A={p.loc['clean','A']:.4f} B={p.loc['clean','B']:.4f}"))

    ranks = [spearmanr(X[:, j], out[:, j]).statistic for j in range(X.shape[1])]
    ok.append(check("transform is rank-preserving (monotone)", min(ranks) > 0.9999,
                    f"min rho={min(ranks):.6f}"))

    ratios = [float(np.median(np.log2((out[:, j] + 1) / (out[:, 0] + 1)))) for j in range(X.shape[1])]
    ok.append(check("all columns on the reference scale after norm", max(abs(r) for r in ratios) < 0.10,
                    f"max |median log2 ratio|={max(abs(r) for r in ratios):.4f}"))

    # anchor guard: independent (uncorrelated) columns have no common-enriched rows -> refuse, don't guess
    rng = np.random.default_rng(1)
    Y = np.abs(rng.lognormal(3, 1.5, size=(N, 2))).astype(np.float32)
    _, q = s3norm_matrix(Y, ref="a", srx=["a", "b"], verbose=False)
    ok.append(check(f"guard: refuses to fit below MIN_ANCHOR={MIN_ANCHOR}",
                    not bool(q.set_index("sample").loc["b", "fitted"]),
                    "uncorrelated columns left unnormalized"))

    print(f"\n{sum(ok)}/{len(ok)} passed")
    return 0 if all(ok) else 1


if __name__ == "__main__":
    sys.exit(main())
