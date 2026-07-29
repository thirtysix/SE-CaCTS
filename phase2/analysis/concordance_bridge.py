#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 6 (stage 1) — the SPECIFICITY-CONCORDANCE BRIDGE between the SE layer and the gene layer.

Everything so far links an SE to a gene by PROXIMITY only (gotcha 22): `score_pilot.py` reports the
nearest gene and whether a hand-curated identity gene happens to sit within 100 kb. That is a validation
proxy, not an association score, and it cannot say whether a given call means anything.

This asks a falsifiable question instead. Run CaCTS on DepMap EXPRESSION over the SAME cell lines and the
SAME Oncotree groups as the SE atlas, giving gene x group specificity alongside SE x group specificity.
Then: are genes sitting near group-specific SEs themselves SPECIFIC TO THAT SAME GROUP, more often than
chance?

The controls are the whole point:
  * BACKGROUND       — per group, the rate at which ALL tested genes are specific to that group. The
                       enrichment is the observed rate over this, with a hypergeometric p-value.
  * GROUP SHUFFLE    — score each SE's neighbour gene against a DIFFERENT group. Concordance must collapse.
  * DISTANCE         — concordance should decay as the SE gets further from the gene. If it does not, the
                       "association" is a lineage-wide confound rather than a local regulatory link.
  * DIRECT rho       — SE H3K27ac signal vs that gene's expression across all lines, which is independent
                       of both specificity computations.

  ~/miniconda3/envs/atac_hdac/bin/python phase2/analysis/concordance_bridge.py
  ... --levels OncotreeLineage,OncotreeSubtype --fdr 0.10
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import hypergeom, spearmanr

_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, _ROOT)
from secacts_env import DATAROOT                      # noqa: E402  (local paths live in .env)
HERE = os.path.dirname(os.path.abspath(__file__))
SECACTS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(DATAROOT, "002.AI_projects", "pyCaCTS"))
sys.path.insert(0, os.path.join(SECACTS, "phase2"))

from pycacts.score import cacts_score_matrix          # noqa: E402
from pycacts.grouping import build_rep_matrix         # noqa: E402
from specificity import fdr_matrix                    # noqa: E402

EXPR = f"{DATAROOT}/DepMap/2026q1/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
MODEL = f"{DATAROOT}/DepMap/2026q1/Model.csv"


def load_expression(model_ids):
    """genes x lines, log2(TPM+1), canonical profile only. Gotcha 16: ModelID is a COLUMN and the
    canonical row is IsDefaultEntryForModel == 'Yes' (string) — this applies to expression too."""
    E = pd.read_csv(EXPR, low_memory=False)
    E = E[E["IsDefaultEntryForModel"] == "Yes"].set_index("ModelID")
    drop = [c for c in E.columns if "(" not in c]
    E = E.drop(columns=drop)
    E.columns = [c.split(" (")[0] for c in E.columns]
    E = E.loc[[m for m in model_ids if m in E.index]]
    E = E.loc[:, ~E.columns.duplicated()]
    return E.T.apply(pd.to_numeric, errors="coerce").dropna(how="any")     # genes x lines


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scores", default=f"{SECACTS}/phase2/scores/atlas.s3")
    ap.add_argument("--signal", default=f"{SECACTS}/phase2/results/atlas.s3.se_signal.tsv.gz")
    ap.add_argument("--pull-set", default=f"{SECACTS}/phase2/data/pull_set.tsv")
    ap.add_argument("--levels", default="OncotreeLineage,OncotreePrimaryDisease,OncotreeSubtype")
    ap.add_argument("--fdr", type=float, default=0.10)
    ap.add_argument("--out", default=f"{SECACTS}/phase2/scores/atlas.s3.concordance")
    a = ap.parse_args()

    ps = pd.read_csv(a.pull_set, sep="\t")
    sm0 = dict(zip(ps["srx"], ps["model_id"]))
    # Apples-to-apples: the gene layer must be scored over the SAME cell lines as the SE layer, or the two
    # CaCTS runs are grouping different panels. The S3norm atlas is QC-gated (--min-peaks 2000), so it holds
    # 282 lines, not the 324 in pull_set — taking the lines from pull_set silently added 41 and changed the
    # subtype group count (78 vs 75). Derive them from the atlas columns instead.
    with gzip.open(a.signal, "rt") as fh:
        atlas_srx = fh.readline().rstrip("\n").split("\t")[1:]
    model_ids = sorted({sm0[s] for s in atlas_srx if isinstance(sm0.get(s), str)})
    md = pd.read_csv(MODEL, index_col="ModelID")
    E = load_expression(model_ids)
    print(f"[bridge] expression: {E.shape[0]:,} genes x {E.shape[1]} lines "
          f"(matched to the {len(model_ids)} lines in the SE atlas)", flush=True)

    # SE signal per line, for the SEs we will actually test (one pass over the big matrix)
    se_needed = set()
    tops = {}
    for lev in a.levels.split(","):
        f = f"{a.scores}.{lev}.top_specific.tsv"
        if os.path.exists(f):
            t = pd.read_csv(f, sep="\t")
            tops[lev] = t
            se_needed |= set(t["se"])
    print(f"[bridge] pulling signal for {len(se_needed):,} SEs", flush=True)
    with gzip.open(a.signal, "rt") as fh:
        cols = fh.readline().rstrip("\n").split("\t")[1:]
        rows = {}
        for line in fh:
            sid = line.split("\t", 1)[0]
            if sid in se_needed:
                rows[sid] = np.array(line.rstrip("\n").split("\t")[1:], dtype=float)
    sm = dict(zip(ps["srx"], ps["model_id"]))
    SS = pd.DataFrame(rows, index=cols)
    SS = SS.groupby(pd.Series([sm.get(s) for s in SS.index], index=SS.index)).mean()     # lines x SE

    summary, pairs = [], []
    for lev, T in tops.items():
        rep, _ = build_rep_matrix(E, md, lev, min_group_n=1)
        rep.columns = [str(c) for c in rep.columns]
        gj = cacts_score_matrix(rep)
        gf = np.power(10.0, fdr_matrix(gj, null="pergroup", scope="global"))
        gspec = gf <= a.fdr                                   # gene x group boolean
        base = gspec.mean(axis=0)                             # per-group background rate
        print(f"\n===== {lev}: {gj.shape[0]:,} genes x {gj.shape[1]} groups; "
              f"background specific-gene rate {base.mean() * 100:.2f}% =====", flush=True)

        groups = [g for g in T["group"].unique() if g in gspec.columns]
        hit = tot = 0
        shuf_hit = shuf_tot = 0
        rng = np.random.default_rng(0)
        for g in groups:
            sub = T[(T["group"] == g) & (T["fdr"] <= a.fdr)]
            other = [x for x in gspec.columns if x != g]
            for _, r in sub.iterrows():
                gene = r["nearest_gene"]
                if not isinstance(gene, str) or gene not in gspec.index:
                    continue
                tot += 1
                is_spec = bool(gspec.loc[gene, g])
                hit += is_spec
                # group shuffle: same gene, a RANDOM OTHER group
                shuf_tot += 1
                shuf_hit += bool(gspec.loc[gene, rng.choice(other)])
                rho = np.nan
                if r["se"] in SS.columns and gene in E.index:
                    common = SS.index.intersection(E.columns)
                    x, y = SS.loc[common, r["se"]].values, E.loc[gene, common].values
                    m = np.isfinite(x) & np.isfinite(y)
                    if m.sum() > 10:
                        rho = spearmanr(x[m], y[m]).statistic
                pairs.append(dict(level=lev, group=g, se=r["se"], gene=gene, dist_kb=r["dist_kb"],
                                  se_rank=r["rank"], se_fdr=r["fdr"],
                                  gene_fdr=float(gf.loc[gene, g]), gene_specific=is_spec, rho=rho))
        if not tot:
            continue
        bg = float(base.mean())
        p = hypergeom.sf(hit - 1, int(gspec.size), int(gspec.values.sum()), tot)
        print(f"  SE-proximal genes specific to the SAME group : {hit}/{tot} = {100 * hit / tot:.1f}%")
        print(f"  background rate (all genes, all groups)      : {100 * bg:.1f}%")
        print(f"  GROUP-SHUFFLED control                       : {shuf_hit}/{shuf_tot} = "
              f"{100 * shuf_hit / max(shuf_tot, 1):.1f}%")
        print(f"  enrichment {hit / tot / max(bg, 1e-9):.1f}x   hypergeometric p = {p:.2e}")
        summary.append(dict(level=lev, n_pairs=tot, concordant=hit, rate=hit / tot, background=bg,
                            shuffled=shuf_hit / max(shuf_tot, 1), enrichment=hit / tot / max(bg, 1e-9), p=p))

    P = pd.DataFrame(pairs)
    if not P.empty:
        P.to_csv(f"{a.out}.pairs.tsv", sep="\t", index=False)
        pd.DataFrame(summary).to_csv(f"{a.out}.summary.tsv", sep="\t", index=False)
        print("\n===== DISTANCE DEPENDENCE (a local regulatory link should decay; a confound will not) =====")
        bins = [0, 10, 25, 50, 100, 250, 10 ** 9]
        lab = ["<10kb", "10-25", "25-50", "50-100", "100-250", ">250kb"]
        P["bin"] = pd.cut(P["dist_kb"], bins=bins, labels=lab, right=False)
        d = P.groupby("bin", observed=True).agg(n=("gene_specific", "size"),
                                                concordant=("gene_specific", "mean"),
                                                median_rho=("rho", "median"))
        for b, r in d.iterrows():
            print(f"  {str(b):<9} n={int(r['n']):<5} concordant={100 * r['concordant']:5.1f}%   "
                  f"median rho={r['median_rho']:+.3f}")
        print(f"\n  overall median rho (SE signal vs neighbour expression) = {P['rho'].median():+.3f}")
        print(f"  rho>0 in {100 * (P['rho'] > 0).mean():.1f}% of pairs   "
              f"(concordant pairs: {100 * (P[P.gene_specific]['rho'] > 0).mean():.1f}%)")
        print(f"\n[bridge] wrote {a.out}.pairs.tsv ({len(P):,} rows) + .summary.tsv")


if __name__ == "__main__":
    main()
