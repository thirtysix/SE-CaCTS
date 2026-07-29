#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Drill down on ONE SE locus from the atlas: where it is, which lines carry it, and — the part that
matters — whether an INDEPENDENT assay agrees.

The specificity score says an SE discriminates a group. It does not say the SE is real, or that the gene
the SE was assigned to is the gene it acts on (assignment is proximity-only; see gotcha 22). This script
supplies the orthogonal check: correlate the SE's H3K27ac signal across cell lines against DepMap RNA
expression of the nearby gene. A promoter/enhancer element that drives its neighbour should track it, and
crucially should track it BETTER than it tracks other lineage-identity genes — that contrast is the test,
not the raw correlation.

  ~/miniconda3/envs/atac_hdac/bin/python phase2/analysis/se_drilldown.py USE_6049 --genes EMX2,PAX8,SOX17

DepMap join gotcha (16): in both OmicsCNGeneWGS.csv and the expression matrix, ModelID is a COLUMN, not
the index, and the canonical row per model is IsDefaultEntryForModel == "Yes" (the string, not a bool).
"""
from __future__ import annotations

import sys

import argparse
import gzip
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, _ROOT)
from secacts_env import DATAROOT                      # noqa: E402  (local paths live in .env)
SECACTS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # .../SE-CaCTS


def read_se_row(signal_gz, se_id):
    opener = gzip.open if str(signal_gz).endswith(".gz") else open
    with opener(signal_gz, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")[1:]
        for line in fh:
            if line.startswith(se_id + "\t"):
                return pd.Series(np.array(line.rstrip("\n").split("\t")[1:], dtype=float), index=header)
    raise SystemExit(f"{se_id} not found in {signal_gz}")


def load_expression(genes):
    f = f"{DATAROOT}/DepMap/2026q1/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
    cols = pd.read_csv(f, nrows=0).columns
    want = {g: c for c in cols for g in genes if c.split(" (")[0] == g}
    if not want:
        return pd.DataFrame()
    E = pd.read_csv(f, usecols=["ModelID", "IsDefaultEntryForModel"] + list(want.values()), low_memory=False)
    E = E[E["IsDefaultEntryForModel"] == "Yes"].set_index("ModelID")     # gotcha 16
    return E[list(want.values())].rename(columns={v: k for k, v in want.items()}).apply(pd.to_numeric,
                                                                                        errors="coerce")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("se_id")
    ap.add_argument("--signal", default=f"{SECACTS}/phase2/results/atlas.s3.se_signal.tsv.gz")
    ap.add_argument("--catalog", default=f"{SECACTS}/phase2/results/atlas.s3.union_catalog.bed.gz")
    ap.add_argument("--pull-set", default=f"{SECACTS}/phase2/data/pull_set.tsv")
    ap.add_argument("--genes", default="", help="comma-separated genes to correlate; the FIRST is the "
                                                "hypothesis, the rest are controls")
    ap.add_argument("--lineage", default=None, help="print the per-line table for this OncotreeLineage")
    ap.add_argument("--top", type=int, default=10)
    a = ap.parse_args()

    S = read_se_row(a.signal, a.se_id)
    ps = pd.read_csv(a.pull_set, sep="\t")
    sm = dict(zip(ps["srx"], ps["model_id"]))
    S = S.groupby(pd.Series([sm.get(s) for s in S.index], index=S.index)).mean()
    md = pd.read_csv(f"{DATAROOT}/DepMap/2026q1/Model.csv", index_col="ModelID")

    opener = gzip.open if a.catalog.endswith(".gz") else open
    with opener(a.catalog, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[3] == a.se_id:
                print(f"{a.se_id}  {f[0]}:{int(f[1]):,}-{int(f[2]):,}  ({int(f[2]) - int(f[1]):,} bp)  "
                      f"called as an SE in {f[4]} of {len(S)} lines' experiments")
                break

    print(f"\nTop {a.top} lines by mean H3K27ac signal (of {len(S)}):")
    for m, v in S.nlargest(a.top).items():
        if m in md.index:
            print(f"  {md.loc[m, 'StrippedCellLineName']:<14} {v:9.1f}  "
                  f"{str(md.loc[m, 'OncotreeLineage'])[:22]:<22} {str(md.loc[m, 'OncotreeSubtype'])[:34]}")
    print(f"  ...median across all lines = {S.median():.1f}")

    genes = [g for g in a.genes.split(",") if g]
    if genes:
        E = load_expression(genes)
        both = S.index.intersection(E.index)
        print(f"\nORTHOGONAL CHECK — SE signal vs DepMap RNA, {len(both)} lines with both:")
        x = S.loc[both].values
        for i, g in enumerate(genes):
            if g not in E.columns:
                print(f"  {g:<8} not in the expression matrix")
                continue
            y = E.loc[both, g].values
            m = np.isfinite(x) & np.isfinite(y)
            r, p = spearmanr(x[m], y[m])
            tag = "  <- hypothesis" if i == 0 else "  (control)"
            print(f"  {g:<8} rho={r:+.3f}  p={p:.2e}  n={int(m.sum())}{tag}")

    if a.lineage:
        sel = [m for m in S.index if m in md.index and md.loc[m, "OncotreeLineage"] == a.lineage]
        E = load_expression(genes) if genes else pd.DataFrame()
        rows = []
        for m in sel:
            r = dict(line=md.loc[m, "StrippedCellLineName"],
                     subtype=str(md.loc[m, "OncotreeSubtype"])[:32], SE_signal=round(float(S[m]), 0))
            if genes and genes[0] in getattr(E, "columns", []) and m in E.index:
                r[f"{genes[0]}_expr"] = round(float(E.loc[m, genes[0]]), 2)
            rows.append(r)
        print(f"\n{a.lineage} lines:")
        print(pd.DataFrame(rows).sort_values("SE_signal", ascending=False).to_string(index=False))


if __name__ == "__main__":
    main()
