#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Scoping: what Hi-C SE -> gene links (hi-c_analyses P1) would do to the master-regulator figure (77_stage_master_tfs.py).

The figure links an SE to a gene when it lies within 100 kb of the gene body. P1 links each SE to the gene, among TSSs
within 1 Mb, with the highest Hi-C contact averaged over the 14 deep ENCODE lines (balanced / that line's expected at
100 kb, 10 kb bins; a TSS within one bin of the SE ranks first). Candidate SEs per gene under each rule:

  100kb          the figure's rule (poster_figures.identity_table)
  +nearest       100kb plus every SE whose nearest TSS within 1 Mb is the gene (the distance-only widening: the control)
  +hic           100kb plus every SE whose top contact gene is the gene (the proposal)
  hic            SEs whose top contact gene is the gene, alone

For each rule and gene: the SE best in the gene's own lineage (lowest FDR, ties by JSD, as the figure), whether it is
called there, and how often the same SE is called in the other 11 lineages of the figure (the figure's off-target
rate). Widening the window can only add own-lineage calls, so the selection-matched null is the "any" rate: for each
other lineage, is ANY candidate SE called there (best-in-that-lineage, the same selection the own column gets)?

  ~/miniconda3/envs/atac_hdac/bin/python phase2/analysis/master_tfs_hic.py [--hic ../hi-c_analyses]
"""
from __future__ import annotations

import argparse
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

SECACTS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FDR = 0.10


def contact_top(H):
    """per SE: its top gene by mean deep-line contact (P1 hic_all, --deep-only; 78_hic_se_links.py) and its nearest-TSS gene."""
    spec = importlib.util.spec_from_file_location("l78", os.path.join(SECACTS, "phase2", "scripts", "78_hic_se_links.py"))
    l78 = importlib.util.module_from_spec(spec); spec.loader.exec_module(l78)
    pr, g = l78.contact_scores(H)
    top = pr.sort_values(["se", "hic", "dist"], ascending=[True, False, True]).groupby("se").head(1)
    near = pr.sort_values(["se", "dist"]).groupby("se").head(1)
    print(f"[hic] {len(g)} deep lines; {pr['se'].nunique():,} SEs; top contact gene = nearest TSS for "
          f"{100 * (top.set_index('se')['gene'] == near.set_index('se')['gene']).mean():.1f}% of SEs")
    return pr, top.set_index("se")["gene"], near.set_index("se")["gene"]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--hic", default=os.path.join(SECACTS, "..", "hi-c_analyses"))
    ap.add_argument("--scores", default=os.path.join(SECACTS, "phase2", "scores_v31f"))
    ap.add_argument("--results", default=os.path.join(SECACTS, "phase2", "results_v31f"))
    a = ap.parse_args()

    spec = importlib.util.spec_from_file_location("pf", os.path.join(SECACTS, "phase2", "figures", "poster_figures.py"))
    pf = importlib.util.module_from_spec(spec); spec.loader.exec_module(pf)
    sys.path.insert(0, os.path.join(SECACTS, "cnrose")); sys.path.insert(0, SECACTS)
    from secacts_env import cache_path
    from cnrose.cn.depmap import load_gene_coords
    gc = load_gene_coords(None, cache_path=cache_path("gene_coords.GRCh38.106.tsv"))
    lev = "OncotreeLineage"
    J = pd.read_csv(os.path.join(a.scores, f"atlas.s3.perm.{lev}.jsd.tsv.gz"), sep="\t", index_col=0)
    F = pd.read_csv(os.path.join(a.scores, f"atlas.s3.perm.{lev}.fdr.tsv.gz"), sep="\t", index_col=0)
    cat = pd.read_csv(os.path.join(a.results, "atlas.s3.union_catalog.bed.gz"), sep="\t", header=None,
                      usecols=[0, 1, 2, 3], names=["chrom", "start", "end", "se"]).set_index("se")
    pr, top, near = contact_top(a.hic)
    cols = [grp for grp, _ in pf.IDENTITY]

    def cands(g):
        c, s, e = gc[g]
        w = set(cat[(cat.chrom == c) & (cat.end >= s - 100_000) & (cat.start <= e + 100_000)].index)
        h, n = set(top.index[top == g]), set(near.index[near == g])
        sets = {"100kb": w, "+nearest": w | n, "+hic": w | h, "hic": h}
        return {k: sorted(x for x in v if x in F.index) for k, v in sets.items()}

    dist = pr.set_index(["se", "gene"])["dist"]
    rows = []
    for grp, genes in pf.IDENTITY:
        for g in genes:
            for rule, S in cands(g).items():
                if not S:
                    rows.append(dict(gene=g, lineage=grp, rule=rule, n=0, own=False)); continue
                best = min(S, key=lambda x: (F.loc[x, grp], J.loc[x, grp]))
                oth = [c for c in cols if c != grp]
                rows.append(dict(gene=g, lineage=grp, rule=rule, n=len(S), se=best, own=F.loc[best, grp] <= FDR,
                                 fdr=round(float(F.loc[best, grp]), 4), rank=int((J[grp] < J.loc[best, grp]).sum()) + 1,
                                 kb=round(dist.get((best, g), np.nan) / 1000) if (best, g) in dist.index else None,
                                 off=int((F.loc[best, oth] <= FDR).sum()), any_off=int((F.loc[S, oth] <= FDR).any().sum()),
                                 any_own=bool((F.loc[S, grp] <= FDR).any())))
    R = pd.DataFrame(rows)
    pd.set_option("display.width", 220); pd.set_option("display.max_rows", 200)
    print("\n== per rule (22 genes x 11 other lineages = 242 off-target cells)")
    for rule, d in R.groupby("rule", sort=False):
        print(f"{rule:9s} candidates median {int(d.n.median()):3d} (total {d.n.sum():4d}) | own called {int(d.own.sum()):2d}/22 | "
              f"figure off-target {int(d.off.sum()):3d}/242 ({100 * d.off.sum() / 242:4.1f}%) | "
              f"any-SE null: own {int(d.any_own.sum())}/22 vs other {int(d.any_off.sum())}/242 ({100 * d.any_off.sum() / 242:4.1f}%)")
    print("\n== genes whose plotted SE changes from the 100 kb rule")
    W = R.pivot_table(index=["gene", "lineage"], columns="rule", values=["se", "rank", "fdr", "kb", "off", "n"], aggfunc="first",
                      sort=False)
    for (g, grp), r in W.iterrows():
        if any(r[("se", k)] != r[("se", "100kb")] for k in ("+nearest", "+hic")):
            print(f"{g:7s} {grp[:22]:22s} " + " | ".join(
                f"{k}: n={int(r[('n', k)]) if pd.notna(r[('n', k)]) else 0} {r[('se', k)]} r{r[('rank', k)]} fdr{r[('fdr', k)]} "
                f"{r[('kb', k)]}kb off{r[('off', k)]}" for k in ("100kb", "+nearest", "+hic")))
    out = os.path.join(SECACTS, "phase2", "analysis", "out", "master_tfs_hic.tsv")
    R.to_csv(out, sep="\t", index=False); print(f"\n[out] {out}")


if __name__ == "__main__":
    main()
