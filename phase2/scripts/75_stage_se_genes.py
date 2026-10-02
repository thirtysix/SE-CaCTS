#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Stage every gene within 100 kb of each super-enhancer, for the dashboard (run after 60/61/73 staging).

The dashboard used to name one gene per SE (the nearest). The nearest gene is not always the target, nor the
only one, so this writes the whole neighbourhood and lets every view list and search it:

  <docs>/data/se_genes.json        {se_id: [[gene, kb], ...]} for every catalogue locus, nearest first. A gene
                                   is listed when its body lies within --window-kb of the SE's edges (the rule
                                   and the gene universe of the concordance analysis, concordance_bridge2.py);
                                   kb = gap between SE and gene body (0 = overlaps, else rounded, at least 1; outside the window rounded up).
                                   With no gene in the window: the single nearest gene at any distance, marked
                                   [gene, kb, 1] (outside the window).
  <docs>/data/expr_specific.json   {lineage|disease: {group: [genes]}}: genes group-specific in DepMap
                                   expression (CaCTS JSD, FDR <= 0.10, exactly as concordance_bridge2.py). Checked
                                   against the concordance pairs table: every flag must agree.
  <docs>/data/gene_index{,.v}.json rebuilt from the staged calls_*.tsv: every listed gene of every call, so the
                                   finder matches any gene near an SE, not only the nearest (d = kb, n = nearest).

  ~/miniconda3/envs/atac_hdac/bin/python phase2/scripts/75_stage_se_genes.py --docs docs
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import sys

import numpy as np
import pandas as pd

SECACTS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, SECACTS)
from secacts_env import DATAROOT                                          # noqa: E402
sys.path.insert(0, os.path.join(SECACTS, "phase2", "analysis"))
sys.path.insert(0, os.path.join(SECACTS, "phase2"))
from concordance_bridge2 import EXPR, GENE_CACHE, GTF, MODEL, load_expression  # noqa: E402
from cnrose.cn.depmap import load_gene_coords                              # noqa: E402
from pycacts.grouping import build_rep_matrix                              # noqa: E402
from pycacts.score import cacts_score_matrix                               # noqa: E402
from specificity import fdr_matrix                                         # noqa: E402

LEVELS = {"lineage": "OncotreeLineage", "disease": "OncotreePrimaryDisease"}


def window_genes(catalog, genes, W):
    """{se: [(gene, gap_bp, outside)]} sorted by gap; gap 0 = overlap. Empty window -> the nearest gene anywhere."""
    bychrom = {}
    for g, (c, s, e) in genes.items():
        bychrom.setdefault(c, []).append((s, e, g))
    for c in bychrom:
        bychrom[c].sort()
    starts = {c: np.array([x[0] for x in v]) for c, v in bychrom.items()}
    longest = max(e - s for s, e, _ in (x for v in bychrom.values() for x in v))
    out = {}
    for se, (c, s, e) in catalog.items():
        v = bychrom.get(c)
        if not v:
            out[se] = []
            continue
        lo = np.searchsorted(starts[c], s - W - longest, "left")
        hi = np.searchsorted(starts[c], e + W, "right")
        hits = []
        for gs, ge, g in v[lo:hi]:
            if ge >= s - W and gs <= e + W:
                hits.append((g, max(0, gs - e, s - ge), 0))
        if not hits:                                      # gene desert: the nearest gene, wherever it is
            gap = [(g, max(0, gs - e, s - ge), 1) for gs, ge, g in v]
            hits = [min(gap, key=lambda t: t[1])]
        out[se] = sorted(hits, key=lambda t: (t[1], t[0]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--docs", default=os.path.join(SECACTS, "docs"))
    ap.add_argument("--catalog", default=f"{SECACTS}/phase2/results_v3/atlas.s3.union_catalog.bed.gz")
    ap.add_argument("--signal", default=f"{SECACTS}/phase2/results_v3/atlas.s3.se_signal.tsv.gz")
    ap.add_argument("--pull-set", default=f"{SECACTS}/phase2/data/pull_set.v3.tsv")
    ap.add_argument("--pairs", default=f"{SECACTS}/phase2/scores_v3c/atlas.s3.perm.concordance2.pairs.tsv.gz")
    ap.add_argument("--window-kb", type=int, default=100)
    ap.add_argument("--fdr", type=float, default=0.10)
    a = ap.parse_args()
    data = os.path.join(a.docs, "data")

    # ---- the lines scored (the concordance panel): QC-pass experiments -> DepMap models
    with gzip.open(a.signal, "rt") as fh:
        srx = set(fh.readline().rstrip("\n").split("\t")[1:])
    ps = pd.read_csv(a.pull_set, sep="\t")
    model_ids = sorted(set(ps.loc[ps["srx"].isin(srx), "model_id"].dropna()))
    E = load_expression(model_ids)                                         # genes x lines
    md = pd.read_csv(MODEL, index_col="ModelID")
    print(f"[75] expression: {E.shape[0]:,} genes x {E.shape[1]} lines", flush=True)

    # ---- SE -> genes, same universe as the concordance analysis (expressed protein-coding genes with coords)
    gc = {g: v for g, v in load_gene_coords(GTF, cache_path=GENE_CACHE).items() if g in E.index}
    catalog = {}
    with gzip.open(a.catalog, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            catalog[f[3]] = (f[0], int(f[1]), int(f[2]))
    W = a.window_kb * 1000
    near = window_genes(catalog, gc, W)
    kb = lambda d: 0 if d == 0 else max(1, int(d / 1000 + 0.5))   # outside the window: rounded UP, so never "100 kb"
    se_genes = {se: [[g, kb(d)] if not o else [g, -(-d // 1000), 1] for g, d, o in v] for se, v in near.items() if v}
    n_in = [sum(1 for x in v if len(x) == 2) for v in se_genes.values()]
    print(f"[75] {len(se_genes):,} of {len(catalog):,} loci; genes within {a.window_kb} kb: median "
          f"{int(np.median(n_in))}, max {max(n_in)}; {sum(x == 0 for x in n_in):,} loci with none (nearest kept)")

    # ---- genes group-specific in expression (concordance_bridge2's gspec), and a parity check on its pairs
    expr = {}
    pairs = pd.read_csv(a.pairs, sep="\t")
    for lv, col in LEVELS.items():
        rep, _ = build_rep_matrix(E, md, col, min_group_n=1)
        rep.columns = [str(c) for c in rep.columns]
        gf = np.power(10.0, fdr_matrix(cacts_score_matrix(rep), null="pergroup", scope="global"))
        gspec = gf <= a.fdr
        expr[lv] = {g: sorted(gspec.index[gspec[g].values]) for g in gspec.columns}
        P = pairs[pairs["level"] == col]
        mine = [g in expr[lv].get(grp, ()) for grp, g in zip(P["group"], P["gene"])]
        bad = int((np.array(mine) != P["concordant"].values).sum())
        print(f"[75] {lv}: {sum(map(len, expr[lv].values())):,} group-specific genes; parity vs pairs "
              f"{len(P) - bad:,}/{len(P):,}", flush=True)
        if bad:
            sys.exit(f"[75] FAIL: {bad} concordance flags disagree with {a.pairs}")

    def dump(name, obj):
        p = os.path.join(data, name)
        with open(p, "w") as fh:
            json.dump(obj, fh, separators=(",", ":"), allow_nan=False)
        print(f"[75] wrote {p} ({os.path.getsize(p) / 1e6:.1f} MB)")

    dump("se_genes.json", se_genes)
    dump("expr_specific.json", expr)

    # ---- the finder index, per analysis variant, from the staged call tables
    for suffix in ["", ".noinf", ".nocn"]:
        idx, missing = {}, 0
        for lv in LEVELS:
            f = os.path.join(data, f"calls_{lv}{suffix}.tsv")
            if not os.path.exists(f):
                continue
            C = pd.read_csv(f, sep="\t")
            for r in C.itertuples():
                gl = se_genes.get(r.se)
                if not gl:
                    missing += 1
                    continue
                spec = set(expr[lv].get(r.group, ()))
                for i, x in enumerate(gl):
                    h = {"lv": lv, "g": r.group, "r": int(r.rank), "fdr": r.fdr, "cn": r.cn_mean,
                         "c": int(x[0] in spec), "d": x[1], "n": int(i == 0)}
                    if len(x) > 2:
                        h["o"] = 1                                 # nearest gene, outside the window
                    idx.setdefault(x[0], []).append(h)
        if not idx:
            continue
        for g in idx:
            idx[g].sort(key=lambda x: (x["r"], x["d"]))
        print(f"[75] gene_index{suffix}: {len(idx):,} genes; calls without a catalogue locus: {missing}")
        dump(f"gene_index{suffix}.json", idx)


if __name__ == "__main__":
    main()
