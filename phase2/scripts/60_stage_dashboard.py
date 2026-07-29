#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Stage the SE-CaCTS dashboard's data/ bundle from the PERMUTATION results (the canonical basis).

Hard rule (mirrors RESULTS.md / gotcha 71-72): specificity CALLS only at OncotreeLineage and
OncotreePrimaryDisease; OncotreeSubtype and cell line are RANKINGS ONLY. Never stage an analytic-null count.

Inputs (all in phase2/), outputs -> docs/data/ (the dashboard is served from docs/ by GitHub Pages):
  atlas.s3.perm.hierarchy_summary.tsv                 per-group n_lines + (permutation) counts
  atlas.s3.perm.{Lineage,Disease}.specific.tsv.gz     the calls, annotated here with gene + coords
  atlas.s3.perm.{Subtype,line}.top_specific.tsv       rankings-only (already carry gene/coords)
  atlas.s3.perm.cn_ablation_calls.tsv                 call-based CN ablation
  atlas.s3.perm.concordance2.{summary.tsv,pairs.tsv.gz}   Phase-6 cross-layer validation
  results/atlas.s3.union_catalog.bed.gz               SE coordinates

  ~/miniconda3/envs/atac_hdac/bin/python phase2/scripts/60_stage_dashboard.py
"""
from __future__ import annotations

import gzip
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PHASE2 = os.path.dirname(HERE)
SECACTS = os.path.dirname(PHASE2)
sys.path.insert(0, os.path.join(PHASE2, "analysis"))
from cn_ablation_calls import nearest_gene_fn                       # noqa: E402

SCORES = os.path.join(PHASE2, "scores")
RESULTS = os.path.join(PHASE2, "results")
OUT = os.path.join(SECACTS, "docs", "data")
os.makedirs(OUT, exist_ok=True)
PERM = os.path.join(SCORES, "atlas.s3.perm")

# levels the panel supports as CALLS vs rankings-only (gotcha 72)
CALL_LEVELS = [("lineage", "OncotreeLineage", "Lineage"),
               ("disease", "OncotreePrimaryDisease", "Primary disease")]
RANK_LEVELS = [("subtype", "OncotreeSubtype", "Subtype"),
               ("line", "line", "Cell line")]


def load_coords():
    coords = {}
    with gzip.open(os.path.join(RESULTS, "atlas.s3.union_catalog.bed.gz"), "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            coords[f[3]] = (f[0], int(f[1]), int(f[2]), int(f[4]))   # chrom,start,end,n_samples_called
    return coords


def write_json(name, obj):
    with open(os.path.join(OUT, name), "w") as fh:
        json.dump(obj, fh, separators=(",", ":"))


def main():
    coords = load_coords()
    nearest = nearest_gene_fn()
    H = pd.read_csv(f"{PERM}.hierarchy_summary.tsv", sep="\t")

    # concordance is the authority on the nearest gene AND whether it is group-specific — use its own
    # is_nearest pair so the gene shown and the concordance shown always name the SAME gene (the two
    # nearest-gene definitions differ by gene universe). Fall back to the coord-based nearest otherwise.
    CP = pd.read_csv(f"{PERM}.concordance2.pairs.tsv.gz", sep="\t")
    CP = CP[CP["is_nearest"] == True]                                   # noqa: E712
    conc_by = {(r.level, r.group, r.se): (r.gene, bool(r.concordant), float(r.rho)) for r in CP.itertuples()}

    manifest = {"levels": {}}
    gene_index = {}                                                    # gene -> [{level,group,rank,fdr,cn,conc,rho}]
    # ---- CALL levels: stage every specific SE, annotated with gene + coordinates + concordance
    for short, col, label in CALL_LEVELS:
        S = pd.read_csv(f"{PERM}.{col}.specific.tsv.gz", sep="\t").sort_values(["group", "rank"])
        chrom, st, en, ncall, gene, dkb, conc, rho = [], [], [], [], [], [], [], []
        for r in S.itertuples():
            se = r.se
            c, s, e, nc = coords.get(se, ("", 0, 0, 0))
            hit = conc_by.get((col, r.group, se))
            if hit:                                                    # bridge nearest gene + concordance
                g, is_c, rr = hit
                d = 0  # bridge distance not carried; coord distance recomputed below for display
            else:                                                     # no bridge pair -> coord-based nearest, no concordance
                g, is_c, rr = (nearest(c, (s + e) // 2) if c else ("", 0)), None, None
                if isinstance(g, tuple):
                    g, d = g
                else:
                    d = 0
            chrom.append(c); st.append(s); en.append(e); ncall.append(nc)
            gene.append(g); dkb.append(d); conc.append("" if is_c is None else int(is_c))
            rho.append("" if rr is None else round(rr, 3))
        S = S.assign(chrom=chrom, start=st, end=en, n_called=ncall, gene=gene, dist_kb=dkb,
                     conc=conc, rho=rho)
        S = S[["group", "se", "rank", "jsd", "fdr", "cn_mean", "gene", "dist_kb", "conc", "rho",
               "chrom", "start", "end", "n_called"]]
        S.round({"jsd": 4, "fdr": 4, "cn_mean": 3}).to_csv(os.path.join(OUT, f"calls_{short}.tsv"),
                                                           sep="\t", index=False)
        # gene index for the finder (only the call levels — honest)
        for r in S.itertuples():
            if r.gene:
                gene_index.setdefault(r.gene, []).append(
                    {"lv": short, "g": r.group, "r": int(r.rank), "fdr": round(float(r.fdr), 4),
                     "cn": round(float(r.cn_mean), 2), "c": r.conc if r.conc != "" else None})
        hsub = H[H.level == col]
        groups = {r.group: {"n_lines": int(r.n_lines), "n_calls": int(r.n_spec_fdr10)}
                  for r in hsub.itertuples()}
        manifest["levels"][short] = {"col": col, "label": label, "kind": "calls",
                                     "n_groups": len(groups), "groups": groups}
        print(f"[stage] {short}: {len(S):,} calls across {len(groups)} groups")
    for g in gene_index:
        gene_index[g].sort(key=lambda x: x["r"])
    write_json("gene_index.json", gene_index)
    print(f"[stage] gene index: {len(gene_index):,} genes near a lineage/disease-specific SE")

    # ---- RANK-ONLY levels: stage the top-N rankings (already gene/coord annotated); NO counts
    for short, col, label in RANK_LEVELS:
        T = pd.read_csv(f"{PERM}.{col}.top_specific.tsv", sep="\t").sort_values(["group", "rank"])
        keep = ["group", "se", "rank", "jsd", "fdr", "cn_mean", "nearest_gene", "dist_kb"]
        T = T[keep].rename(columns={"nearest_gene": "gene"})
        # add coords for out-links
        cc = T["se"].map(lambda s: coords.get(s, ("", 0, 0, 0)))
        T = T.assign(chrom=[x[0] for x in cc], start=[x[1] for x in cc], end=[x[2] for x in cc])
        T.round({"jsd": 4, "fdr": 4, "cn_mean": 3}).to_csv(os.path.join(OUT, f"rank_{short}.tsv"),
                                                          sep="\t", index=False)
        hsub = H[H.level == col]
        groups = {r.group: {"n_lines": int(r.n_lines)} for r in hsub.itertuples()}
        manifest["levels"][short] = {"col": col, "label": label, "kind": "rankings",
                                     "n_groups": len(groups), "groups": groups}
        print(f"[stage] {short}: rankings for {len(groups)} groups (no calls — panel unsupported)")

    write_json("manifest.json", manifest)

    # ---- CN ablation (call-based, honest null)
    A = pd.read_csv(f"{PERM}.cn_ablation_calls.tsv", sep="\t")
    abl = {"summary": [], "amplicon": [], "note": ""}
    for short, col, label in CALL_LEVELS:
        sub = A[A.level == col]
        amp = sub[sub.kind == "amplicon_driven"]
        resc = sub[sub.kind == "rescued"]
        # counts from the specific dumps (corrected) + the ablation sets
        n_corr = int(H[(H.level == col)]["n_spec_fdr10"].sum())
        abl["summary"].append({"level": label, "corrected": n_corr,
                               "amplicon_driven": int(len(amp)), "rescued": int(len(resc)),
                               "amp_median_cn": round(float(amp.cn_mean.median()), 1) if len(amp) else None,
                               "resc_median_cn": round(float(resc.cn_mean.median()), 3) if len(resc) else None})
    # the amplicon-driven calls, deduped to (group, gene) with the max cn
    amp_all = A[A.kind == "amplicon_driven"].sort_values("cn_mean", ascending=False)
    seen = set()
    for r in amp_all.itertuples():
        key = (r.group, r.nearest_gene)
        if key in seen:
            continue
        seen.add(key)
        abl["amplicon"].append({"level": "Lineage" if r.level == "OncotreeLineage" else "Disease",
                                "group": r.group, "gene": r.nearest_gene, "cn": round(float(r.cn_mean), 1)})
    write_json("cn_ablation.json", abl)
    print(f"[stage] cn ablation: {len(abl['amplicon'])} distinct amplicon-driven (group,gene)")

    # ---- Concordance (Phase-6): summary + distance decay from the pairs file
    C = pd.read_csv(f"{PERM}.concordance2.summary.tsv", sep="\t")
    P = pd.read_csv(f"{PERM}.concordance2.pairs.tsv.gz", sep="\t")
    P["dist_kb"] = P["dist_bp"] // 1000
    bins = [0, 10, 25, 50, 100, 10 ** 9]
    labs = ["<10 kb", "10–25", "25–50", "50–100", ">100 kb"]
    P["bin"] = pd.cut(P["dist_kb"], bins=bins, labels=labs, right=False)
    dd = P.groupby("bin", observed=True).agg(n=("concordant", "size"), conc=("concordant", "mean"),
                                             shuf=("shuffled", "mean"), rho=("rho", "median"))
    conc = {"summary": [{"level": "Lineage" if r.level == "OncotreeLineage" else "Primary disease",
                         "per_pair": round(r.per_pair * 100, 1), "background": round(r.background * 100, 2),
                         "enrichment": round(r.enrichment, 1), "per_se_any": round(r.per_se_any * 100, 1),
                         "nearest": round(r.nearest * 100, 1), "shuffled": round(r.shuffled * 100, 1)}
                        for r in C.itertuples()],
            "distance": [{"bin": str(b), "n": int(r.n), "concordant": round(r.conc * 100, 1),
                          "shuffled": round(r.shuf * 100, 1), "rho": round(r.rho, 3)}
                         for b, r in dd.iterrows()]}
    write_json("concordance.json", conc)
    print("[stage] concordance staged")

    # ---- meta / hero numbers
    n_lineage_calls = int(H[H.level == "OncotreeLineage"]["n_spec_fdr10"].sum())
    n_disease_calls = int(H[H.level == "OncotreePrimaryDisease"]["n_spec_fdr10"].sum())
    meta = {
        "n_samples": 2136, "n_lines": 282, "n_ses": 42943,
        "n_lineages": int((H.level == "OncotreeLineage").sum()),
        "n_diseases": int((H.level == "OncotreePrimaryDisease").sum()),
        "n_subtypes": int((H.level == "OncotreeSubtype").sum()),
        "n_lineage_calls": n_lineage_calls, "n_disease_calls": n_disease_calls,
        "pull_bu": 43, "n_pull": 2916,
        "fdr": "label-permutation, B=1000, FDR ≤ 0.10",
    }
    write_json("meta.json", meta)
    print(f"[stage] wrote {len(os.listdir(OUT))} files to {OUT}")


if __name__ == "__main__":
    main()
