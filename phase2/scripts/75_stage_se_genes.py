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
  <docs>/data/se_lnc.json          {se_id: [[lncRNA, kb], ...]}: HGNC-named lncRNA genes (GENCODE v36, gene_type
                                   lncRNA with an hgnc_id; CCAT1 is not annotated there) within --window-kb, nearest
                                   first. Listed after the protein-coding genes; never "nearest"; no expression mark
                                   (DepMap expression is protein-coding only). v3.1, ROADMAP item 15.
  <docs>/data/gene_index{,.v}.json rebuilt from the staged calls_*.tsv: every listed gene of every call, so the
                                   finder matches any gene near an SE, not only the nearest (d = kb, n = nearest or tied with it;
                                   h = 1: the SE's Hi-C contact gene, added when it lies beyond the window).
  <docs>/data/se_hic.json          {se_id: [gene, kb, contact, n_genes, in_window]}: the gene with the highest measured Hi-C
                                   contact with the SE among the >= --hic-min-genes genes with a TSS within 1 Mb
                                   (78_hic_se_links.py; 14 ENCODE cancer-line maps). Left out where a TSS lies within one
                                   10 kb bin of the SE: contact cannot be measured there and P1 ranks that gene first anyway.
  <docs>/data/gene_lines/*.json    the finder's cell-line index, by the gene's first character: {gene: [[line, rank, pass,
                                   kb, flags], ...]} from the per-line tables (FDR <= 0.10 vs all lines / lineage / disease /
                                   subtype, pass = A/L/D/S as the atlas); flags 1 Hi-C contact gene, 2 lncRNA, 4 nearest gene
                                   outside the window, 8 nearest (or tied); lines.json = [[key, name, lineage, atlas group], ...].

  ~/miniconda3/envs/atac_hdac/bin/python phase2/scripts/75_stage_se_genes.py --docs docs
"""
from __future__ import annotations

import argparse
import glob
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
SUBTYPE = ("subtype", "OncotreeSubtype")        # a call level from v3.1: added when the concordance pairs carry it


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


def lnc_genes(gtf):
    """{name: (chrom, start, end)} for HGNC-named lncRNA genes in a GENCODE GTF (chrom already 'chr'-prefixed)."""
    import re
    out = {}
    opener = gzip.open if gtf.endswith(".gz") else open
    with opener(gtf, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 9)
            if f[2] != "gene" or 'gene_type "lncRNA"' not in f[8] or "hgnc_id" not in f[8]:
                continue
            out[re.search(r'gene_name "([^"]+)"', f[8]).group(1)] = (f[0], int(f[3]) - 1, int(f[4]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--docs", default=os.path.join(SECACTS, "docs"))
    ap.add_argument("--catalog", default=f"{SECACTS}/phase2/results_v3/atlas.s3.union_catalog.bed.gz")
    ap.add_argument("--signal", default=f"{SECACTS}/phase2/results_v3/atlas.s3.se_signal.tsv.gz")
    ap.add_argument("--pull-set", default=f"{SECACTS}/phase2/data/pull_set.v3.tsv")
    ap.add_argument("--pairs", default=f"{SECACTS}/phase2/scores_v3c/atlas.s3.perm.concordance2.pairs.tsv.gz")
    ap.add_argument("--window-kb", type=int, default=100)
    ap.add_argument("--lnc-gtf", default=os.path.join(DATAROOT, "0.human_genome/gencode/v36/gencode.v36.annotation.gtf.gz"),
                    help="GENCODE GTF for HGNC-named lncRNAs ('' = none)")
    ap.add_argument("--fdr", type=float, default=0.10)
    ap.add_argument("--hic-links", default=f"{SECACTS}/phase2/data/hic_se_links.tsv.gz", help="78_hic_se_links.py; 'none' = no Hi-C genes")
    ap.add_argument("--hic-min-genes", type=int, default=2)
    ap.add_argument("--line-scores", default=None, help="dir with atlas.s3.lines.{all,lin,dis,sub}.line.specific.tsv.gz "
                                                         "(the cell-line finder index); omitted = none")
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
    se_lnc = {}
    if a.lnc_gtf:
        # a name that is a protein-coding gene here (DepMap expression) and an lncRNA in GENCODE v36 (TMEM78, CCDC39,
        # SLFN12L) is listed once, as the protein-coding gene
        lg = lnc_genes(a.lnc_gtf)
        dup = sorted(set(lg) & set(E.index))
        ln = window_genes(catalog, {g: v for g, v in lg.items() if g not in E.index}, W)
        print(f"[75] lncRNAs: {len(dup)} names also protein-coding here, left out as lncRNA: {dup[:8]}")
        se_lnc = {se: [[g, kb(d)] for g, d, o in v if not o] for se, v in ln.items()}
        se_lnc = {se: v for se, v in se_lnc.items() if v}
        print(f"[75] lncRNAs: {len(se_lnc):,} loci have >= 1 HGNC-named lncRNA within {a.window_kb} kb")
    n_in = [sum(1 for x in v if len(x) == 2) for v in se_genes.values()]
    print(f"[75] {len(se_genes):,} of {len(catalog):,} loci; genes within {a.window_kb} kb: median "
          f"{int(np.median(n_in))}, max {max(n_in)}; {sum(x == 0 for x in n_in):,} loci with none (nearest kept)")

    # ---- genes group-specific in expression (concordance_bridge2's gspec), and a parity check on its pairs
    expr = {}
    pairs = pd.read_csv(a.pairs, sep="\t")
    if (pairs["level"] == SUBTYPE[1]).any():
        LEVELS[SUBTYPE[0]] = SUBTYPE[1]
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
    if se_lnc:
        dump("se_lnc.json", se_lnc)
    dump("expr_specific.json", expr)
    se_hic = {}
    if a.hic_links != "none":
        T = pd.read_csv(a.hic_links, sep="\t")
        T = T[(T["n_genes"] >= a.hic_min_genes) & np.isfinite(T["top_contact"]) & T["se"].isin(catalog)]
        for r in T.itertuples():
            w = any(x[0] == r.top_gene and len(x) == 2 for x in se_genes.get(r.se, ()))
            se_hic[r.se] = [r.top_gene, kb(r.top_kb * 1000), round(float(r.top_contact), 3), int(r.n_genes), int(w)]
        print(f"[75] Hi-C contact genes: {len(se_hic):,} loci ({sum(1 for v in se_hic.values() if not v[4]):,} beyond "
              f"{a.window_kb} kb)")
        dump("se_hic.json", se_hic)

    # ---- the finder index, per analysis variant, from the staged call tables
    suffixes = sorted({os.path.basename(f)[len("calls_lineage"):-len(".tsv")]
                       for f in glob.glob(os.path.join(data, "calls_lineage*.tsv"))} - {".cndep"}) or [""]
    for suffix in suffixes:
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
                hg = se_hic.get(r.se)
                for x in gl:
                    h = {"lv": lv, "g": r.group, "r": int(r.rank), "fdr": r.fdr, "cn": r.cn_mean,
                         "c": int(x[0] in spec), "d": x[1], "n": int(x[1] == gl[0][1])}   # ties count as nearest
                    if len(x) > 2:
                        h["o"] = 1                                 # nearest gene, outside the window
                    if hg and hg[0] == x[0]:
                        h["h"] = 1                                 # also its Hi-C contact gene
                    idx.setdefault(x[0], []).append(h)
                if hg and not any(x[0] == hg[0] for x in gl):      # the Hi-C contact gene beyond the window
                    idx.setdefault(hg[0], []).append({"lv": lv, "g": r.group, "r": int(r.rank), "fdr": r.fdr,
                                                      "cn": r.cn_mean, "c": int(hg[0] in spec), "d": hg[1], "n": 0, "h": 1})
                for x in se_lnc.get(r.se, ()):                     # lncRNAs: listed, never nearest, no ⇌
                    idx.setdefault(x[0], []).append({"lv": lv, "g": r.group, "r": int(r.rank), "fdr": r.fdr,
                                                     "cn": r.cn_mean, "c": 0, "d": x[1], "n": 0, "t": "lnc"})
        if not idx:
            continue
        for g in idx:
            idx[g].sort(key=lambda x: (x["r"], x["d"]))
        print(f"[75] gene_index{suffix}: {len(idx):,} genes; calls without a catalogue locus: {missing}")
        dump(f"gene_index{suffix}.json", idx)

    # ---- the cell-line finder index: every per-line call (any of the four comparisons), sharded by first character
    if a.line_scores:
        L = {}
        for f in glob.glob(os.path.join(data, "lines", "*.json")):
            j = json.load(open(f))
            if "rows" in j:               # the scoring names DepMap lines by stripped name, others by their name (61 score_name)
                L[j["group"] if j["key"].startswith("ACH-") else j["name"]] = [j["key"], j["name"], j["lineage"], j["group"]]
        calls = {}
        for c, letter in (("all", "A"), ("lin", "L"), ("dis", "D"), ("sub", "S")):
            d = pd.read_csv(os.path.join(a.line_scores, f"atlas.s3.lines.{c}.line.specific.tsv.gz"), sep="\t")
            d = d[d["fdr"] <= a.fdr]
            unknown = sorted(set(d["group"]) - set(L))
            if unknown:
                sys.exit(f"[75] FAIL: per-line groups with no staged line: {unknown[:5]}")
            for r in d.itertuples():
                k = (r.group, r.se)
                if k in calls:
                    calls[k][1] += letter
                else:
                    calls[k] = [int(r.rank), letter]
        keys = sorted(L, key=lambda g: L[g][1].lower())
        li = {g: i for i, g in enumerate(keys)}
        shards = {}
        for (g, se), (rank, ps) in calls.items():
            ps = "".join(x for x in "ALDS" if x in ps)
            gl, hg = se_genes.get(se, []), se_hic.get(se)
            ent = []
            for x in gl:
                ent.append((x[0], x[1], (4 if len(x) > 2 else 0) | (8 if x[1] == gl[0][1] else 0) | (1 if hg and hg[0] == x[0] else 0)))
            if hg and not any(x[0] == hg[0] for x in gl):
                ent.append((hg[0], hg[1], 1))
            ent += [(x[0], x[1], 2) for x in se_lnc.get(se, ())]
            for gene, d, fl in ent:
                ch = gene[0].upper() if gene[0].isalpha() else "_"
                shards.setdefault(ch, {}).setdefault(gene, []).append([li[g], rank, ps, d, fl])
        os.makedirs(os.path.join(data, "gene_lines"), exist_ok=True)
        for f in glob.glob(os.path.join(data, "gene_lines", "*.json")):
            os.remove(f)
        dump("gene_lines/lines.json", [L[g] for g in keys])
        n = 0
        for ch, idx in shards.items():
            for v in idx.values():
                v.sort(key=lambda x: (x[1], x[3]))
                n += len(v)
            dump(f"gene_lines/{ch}.json", idx)
        print(f"[75] cell-line index: {len(calls):,} line x SE calls in {len(L)} lines; {n:,} gene entries in {len(shards)} shards")


if __name__ == "__main__":
    main()
