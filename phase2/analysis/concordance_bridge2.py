#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 6 STAGE 2 — the concordance bridge without the optimistic sampling.

Stage 1 (`concordance_bridge.py`) measured concordance on the reported top-15 SEs per group against the
single NEAREST gene. Both choices flatter the result: the top-15 are the strongest calls, and the nearest
gene is the most favourable assignment. This runs the same question over **every SE x group test at
FDR <= 0.10** and **every gene within a window**, which is the honest genome-wide number.

Three rates are reported because they answer different questions:
  PER-PAIR   — of all (specific SE, nearby gene) pairs, how many are gene-specific to the same group?
               Lowest, and the right number for "is a given SE-gene link real".
  PER-SE ANY — of specific SEs, how many have AT LEAST ONE concordant gene in the window? The number that
               matters for "does this SE call point at group-relevant biology".
  NEAREST    — restricted to the nearest gene, i.e. stage 1's definition, for continuity.

Controls, as in stage 1: per-group background rate + hypergeometric p; a GROUP SHUFFLE; distance decay; and
a direct Spearman of SE signal vs gene expression across lines, computed vectorized (rank-transform once,
then Pearson per pair) so millions of pairs stay tractable.

  ~/miniconda3/envs/atac_hdac/bin/python phase2/analysis/concordance_bridge2.py --window-kb 100
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import hypergeom, rankdata

_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, _ROOT)
from secacts_env import DATAROOT, cache_path                      # noqa: E402  (local paths live in .env)
HERE = os.path.dirname(os.path.abspath(__file__))
SECACTS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(DATAROOT, "002.AI_projects", "pyCaCTS"))
sys.path.insert(0, os.path.join(SECACTS, "phase2"))
sys.path.insert(0, os.path.join(SECACTS, "cnrose"))

from pycacts.score import cacts_score_matrix          # noqa: E402
from pycacts.grouping import build_rep_matrix         # noqa: E402
from specificity import fdr_matrix                    # noqa: E402
from cnrose.cn.depmap import load_gene_coords         # noqa: E402

EXPR = f"{DATAROOT}/DepMap/2026q1/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
MODEL = f"{DATAROOT}/DepMap/2026q1/Model.csv"
GTF = f"{DATAROOT}/0.human_genome/Homo_sapiens.GRCh38.106.chr.gtf.gz"
GENE_CACHE = cache_path("gene_coords.GRCh38.106.tsv")


def load_expression(model_ids):
    """genes x lines, log2(TPM+1). Gotcha 16/65: ModelID is a COLUMN and the canonical row is
    IsDefaultEntryForModel == 'Yes' (string). index_col=0 gives a silent 0-row join."""
    E = pd.read_csv(EXPR, low_memory=False)
    E = E[E["IsDefaultEntryForModel"] == "Yes"].set_index("ModelID")
    E = E.drop(columns=[c for c in E.columns if "(" not in c])
    E.columns = [c.split(" (")[0] for c in E.columns]
    E = E.loc[[m for m in model_ids if m in E.index]]
    E = E.loc[:, ~E.columns.duplicated()]
    return E.T.apply(pd.to_numeric, errors="coerce").dropna(how="any")


def rank_rows(A):
    """Row-wise rank transform; Pearson on these == Spearman."""
    R = np.apply_along_axis(rankdata, 1, A)
    R = R - R.mean(axis=1, keepdims=True)
    n = np.sqrt((R ** 2).sum(axis=1))
    n[n == 0] = 1.0
    return R / n[:, None]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scores", default=f"{SECACTS}/phase2/scores/atlas.s3")
    ap.add_argument("--signal", default=f"{SECACTS}/phase2/results/atlas.s3.se_signal.tsv.gz")
    ap.add_argument("--catalog", default=f"{SECACTS}/phase2/results/atlas.s3.union_catalog.bed.gz")
    ap.add_argument("--pull-set", default=f"{SECACTS}/phase2/data/pull_set.tsv")
    ap.add_argument("--levels", default="OncotreeLineage,OncotreePrimaryDisease,OncotreeSubtype")
    ap.add_argument("--fdr", type=float, default=0.10)
    ap.add_argument("--window-kb", type=int, default=100)
    ap.add_argument("--out", default=f"{SECACTS}/phase2/scores/atlas.s3.concordance2")
    a = ap.parse_args()
    W = a.window_kb * 1000

    ps = pd.read_csv(a.pull_set, sep="\t")
    sm = dict(zip(ps["srx"], ps["model_id"]))
    md = pd.read_csv(MODEL, index_col="ModelID")

    # ---- SE signal, collapsed to cell lines (gotcha 68: the panels must match)
    print("[bridge2] loading SE signal", flush=True)
    S = pd.read_csv(a.signal, sep="\t", index_col=0)
    col_model = pd.Series([sm.get(c) for c in S.columns], index=S.columns).dropna()
    SL = S[col_model.index].T.groupby(col_model).mean()                  # lines x SE
    del S
    model_ids = list(SL.index)
    E = load_expression(model_ids)                                       # genes x lines
    lines = [m for m in model_ids if m in E.columns]
    SL, E = SL.loc[lines], E[lines]
    print(f"[bridge2] {SL.shape[1]:,} SEs x {len(lines)} lines ; {E.shape[0]:,} genes", flush=True)

    # ---- SE -> genes within the window
    coords = {}
    with gzip.open(a.catalog, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            coords[f[3]] = (f[0], int(f[1]), int(f[2]))
    gc = load_gene_coords(GTF, cache_path=GENE_CACHE)
    gc = {g: v for g, v in gc.items() if g in E.index}
    bychrom = {}
    for g, (c, s, e) in gc.items():
        bychrom.setdefault(c, []).append((s, e, g))
    for c in bychrom:
        bychrom[c].sort()
    starts = {c: np.array([x[0] for x in v]) for c, v in bychrom.items()}

    def genes_near(se):
        c, s, e = coords[se]
        v = bychrom.get(c)
        if not v:
            return []
        lo = np.searchsorted(starts[c], s - W - 3_000_000, "left")
        hi = np.searchsorted(starts[c], e + W, "right")
        out = []
        for gs, ge, g in v[lo:hi]:
            if ge >= s - W and gs <= e + W:
                mid = (s + e) // 2
                out.append((g, 0 if (ge >= s and gs <= e) else min(abs(gs - mid), abs(ge - mid))))
        return out

    # ---- rank transforms once, for the vectorized Spearman
    print("[bridge2] rank-transforming for vectorized Spearman", flush=True)
    se_index = {v: i for i, v in enumerate(SL.columns)}
    gene_index = {v: i for i, v in enumerate(E.index)}
    SR = rank_rows(SL.values.T.astype(float))                            # SE x lines
    ER = rank_rows(E.values.astype(float))                               # genes x lines

    rows_all, summary = [], []
    for lev in a.levels.split(","):
        f = f"{a.scores}.{lev}.specific.tsv.gz"
        if not os.path.exists(f):
            print(f"  [skip] {lev}: no {os.path.basename(f)} — run score_pilot.py --dump-specific 0.10")
            continue
        SP = pd.read_csv(f, sep="\t")
        SP = SP[SP["fdr"] <= a.fdr]

        rep, _ = build_rep_matrix(E, md, lev, min_group_n=1)
        rep.columns = [str(c) for c in rep.columns]
        gj = cacts_score_matrix(rep)
        gf = np.power(10.0, fdr_matrix(gj, null="pergroup", scope="global"))
        gspec = (gf <= a.fdr)
        bg = float(gspec.values.mean())
        gpos = {g: i for i, g in enumerate(gspec.columns)}
        GS = gspec.values
        gidx_of = {g: i for i, g in enumerate(gspec.index)}

        SP = SP[SP["group"].isin(gpos)]
        print(f"\n===== {lev}: {len(SP):,} specific SE x group tests; "
              f"{gj.shape[0]:,} genes x {gj.shape[1]} groups; background {100 * bg:.2f}% =====", flush=True)

        near_cache = {}
        rng = np.random.default_rng(0)
        allg = list(gspec.columns)
        recs = []
        for se, grp in zip(SP["se"].values, SP["group"].values):
            if se not in near_cache:
                near_cache[se] = genes_near(se)
            gl = near_cache[se]
            if not gl:
                continue
            gi = gpos[grp]
            shuf = gpos[allg[int(rng.integers(len(allg)))]]
            nearest = min(gl, key=lambda t: t[1])[0]
            for g, d in gl:
                r = gidx_of.get(g)
                if r is None:
                    continue
                recs.append((se, grp, g, d, bool(GS[r, gi]), bool(GS[r, shuf]), g == nearest))
        if not recs:
            continue
        P = pd.DataFrame(recs, columns=["se", "group", "gene", "dist_bp", "concordant", "shuffled", "is_nearest"])
        # vectorized Spearman for every pair
        si = np.array([se_index[s] for s in P["se"]])
        gi2 = np.array([gene_index[g] for g in P["gene"]])
        P["rho"] = np.einsum("ij,ij->i", SR[si], ER[gi2])
        P["level"] = lev
        rows_all.append(P)

        n_pair, n_conc = len(P), int(P["concordant"].sum())
        per_se = P.groupby(["se", "group"])["concordant"].any()
        nearP = P[P["is_nearest"]]
        p = hypergeom.sf(n_conc - 1, int(GS.size), int(GS.sum()), n_pair)
        print(f"  PER-PAIR    {n_conc:,}/{n_pair:,} = {100 * n_conc / n_pair:5.1f}%   "
              f"(background {100 * bg:.2f}%, enrichment {n_conc / n_pair / max(bg, 1e-9):.1f}x, p={p:.2e})")
        print(f"  PER-SE ANY  {int(per_se.sum()):,}/{len(per_se):,} = {100 * per_se.mean():5.1f}%   "
              f"(>=1 concordant gene within {a.window_kb} kb)")
        print(f"  NEAREST     {int(nearP['concordant'].sum()):,}/{len(nearP):,} = "
              f"{100 * nearP['concordant'].mean():5.1f}%   (stage-1 definition, all specific SEs)")
        print(f"  SHUFFLED    {int(P['shuffled'].sum()):,}/{n_pair:,} = {100 * P['shuffled'].mean():5.1f}%")
        summary.append(dict(level=lev, n_pairs=n_pair, per_pair=n_conc / n_pair, per_se_any=float(per_se.mean()),
                            nearest=float(nearP["concordant"].mean()), shuffled=float(P["shuffled"].mean()),
                            background=bg, enrichment=n_conc / n_pair / max(bg, 1e-9), p=p))

    if rows_all:
        A = pd.concat(rows_all, ignore_index=True)
        A["dist_kb"] = A["dist_bp"] // 1000
        print("\n===== DISTANCE DEPENDENCE, all levels pooled =====")
        bins = [0, 10, 25, 50, 100, 250, 10 ** 9]
        lab = ["<10kb", "10-25", "25-50", "50-100", "100-250", ">250kb"]
        A["bin"] = pd.cut(A["dist_kb"], bins=bins, labels=lab, right=False)
        d = A.groupby("bin", observed=True).agg(n=("concordant", "size"), conc=("concordant", "mean"),
                                                shuf=("shuffled", "mean"), rho=("rho", "median"))
        for b, r in d.iterrows():
            print(f"  {str(b):<9} n={int(r['n']):<8} concordant={100 * r['conc']:5.1f}%   "
                  f"shuffled={100 * r['shuf']:4.1f}%   median rho={r['rho']:+.3f}")
        print(f"\n  overall median rho={A['rho'].median():+.3f} ; concordant pairs "
              f"{A[A.concordant]['rho'].median():+.3f} vs discordant {A[~A.concordant]['rho'].median():+.3f}")
        A.drop(columns=["bin"]).to_csv(f"{a.out}.pairs.tsv.gz", sep="\t", index=False, compression="gzip")
        pd.DataFrame(summary).to_csv(f"{a.out}.summary.tsv", sep="\t", index=False)
        print(f"\n[bridge2] wrote {a.out}.pairs.tsv.gz ({len(A):,} rows) + .summary.tsv")


if __name__ == "__main__":
    main()
