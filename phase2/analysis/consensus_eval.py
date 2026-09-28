#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Known-biology test for the SMALL-GROUP levels (subtype, cell line), written BEFORE the consensus runs
were read (2026-09-28), so the gene lists are not tuned to the results.

For each (group, gene) below: does any SE within 100 kb of the gene pass FDR <= 0.10 in its OWN group, and
how often does the same gene pass in the OTHER groups of that level (background)? A score that finds
identity biology should pass on-target far more often than off-target.

  ~/miniconda3/envs/atac_hdac/bin/python phase2/analysis/consensus_eval.py \\
      --scores phase2/scores_v2/atlas.s3.cons --level OncotreeSubtype
  ~/miniconda3/envs/atac_hdac/bin/python phase2/analysis/consensus_eval.py \\
      --scores phase2/scores_v2/atlas.s3.consx --level line
"""
from __future__ import annotations

import argparse
import gzip
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

HERE = os.path.dirname(os.path.abspath(__file__))
SECACTS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, SECACTS)
sys.path.insert(0, os.path.join(SECACTS, "cnrose"))
from secacts_env import DATAROOT, cache_path                                  # noqa: E402
from cnrose.cn.depmap import load_gene_coords                                  # noqa: E402

SUBTYPE_ID = {
    "High-Grade Serous Ovarian Cancer": ["PAX8", "WT1", "SOX17", "MECOM"],
    "Small Cell Lung Cancer": ["ASCL1", "NEUROD1", "INSM1"],
    "Lung Adenocarcinoma": ["NKX2-1"],
    "Colon Adenocarcinoma": ["CDX2", "CDX1", "HNF4A"],
    "Diffuse Large B-Cell Lymphoma, NOS": ["PAX5", "BCL6", "POU2AF1"],
    "Acute Myeloid Leukemia": ["SPI1", "CEBPA", "IRF8"],
    "Neuroblastoma": ["PHOX2B", "HAND2", "GATA3", "ISL1"],
    "T-Cell Acute Lymphoblastic Leukemia": ["TCF7", "LEF1", "BCL11B", "GATA3"],
    "B-Cell Acute Lymphoblastic Leukemia": ["PAX5", "EBF1"],
    "Ewing Sarcoma": ["NKX2-2", "NR0B1"],
    "Plasma Cell Myeloma": ["IRF4", "PRDM1", "XBP1"],
    "Melanoma": ["SOX10", "MITF"],
    "Breast Invasive Lobular Carcinoma": ["ESR1", "FOXA1"],
    "Hepatocellular Carcinoma": ["HNF1A", "HNF4A"],
    "Prostate Adenocarcinoma": ["AR", "NKX3-1", "HOXB13", "FOXA1"],
    "Mantle Cell Lymphoma": ["CCND1", "SOX11"],
    "Esophageal Squamous Cell Carcinoma": ["TP63", "SOX2"],
    "Lung Squamous Cell Carcinoma": ["TP63", "SOX2"],
    "Renal Clear Cell Carcinoma": ["PAX8", "PAX2"],
    "Chronic Myeloid Leukemia, BCR-ABL1+": ["GATA1", "TAL1"],
    "Osteosarcoma": ["RUNX2"],
    "Bladder Urothelial Carcinoma": ["GATA3", "PPARG", "FOXA1"],
}
LINE_ID = {                                  # keyed by DepMap StrippedCellLineName (the line-level group)
    "MCF7": ["ESR1", "GATA3", "FOXA1"], "T47D": ["ESR1", "GATA3", "FOXA1"],
    "THP1": ["CEBPA", "SPI1"], "MOLM13": ["IRF8", "SPI1"], "HL60": ["CEBPA", "SPI1"],
    "SKOV3": ["MECOM", "PAX8"], "SW48": ["CDX2"], "LOVO": ["CDX2"],
    "P12ICHIKAWA": ["LEF1", "TCF7"], "JURKAT": ["TCF7", "LEF1", "TAL1"],
    "K562": ["GATA1", "TAL1"], "HEPG2": ["HNF1A", "HNF4A"],
    "22RV1": ["AR"], "VCAP": ["AR"], "KELLY": ["PHOX2B", "HAND2"], "SKNBE2": ["PHOX2B", "HAND2"],
    "A375": ["SOX10", "MITF"], "NCIH2171": ["ASCL1"], "MM1S": ["IRF4", "PRDM1"], "U266B1": ["IRF4", "PRDM1"],
    "NALM6": ["PAX5", "EBF1"], "KARPAS422": ["BCL6", "PAX5"], "DOHH2": ["BCL6", "PAX5"],
}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scores", required=True, help="scores prefix (…/atlas.s3.cons or …/atlas.s3.consx)")
    ap.add_argument("--level", required=True, choices=["OncotreeSubtype", "line"])
    ap.add_argument("--catalog", default=os.path.join(SECACTS, "phase2/results_v2/atlas.s3.union_catalog.bed.gz"))
    ap.add_argument("--window", type=int, default=100_000)
    ap.add_argument("--min-members", type=int, default=2,
                    help="drop groups with fewer members (their consensus calls are not group properties)")
    a = ap.parse_args()

    ids = SUBTYPE_ID if a.level == "OncotreeSubtype" else LINE_ID
    H = pd.read_csv(f"{a.scores}.hierarchy_summary.tsv", sep="\t")
    Hl = H[H.level == a.level]
    groups = list(Hl[Hl.n_lines >= a.min_members].group)
    if a.level == "line":                    # the line-level output names lines by StrippedCellLineName
        md = pd.read_csv(os.path.join(DATAROOT, "DepMap/2026q1/Model.csv"))
        name = dict(zip(md.StrippedCellLineName, md.StrippedCellLineName))
        groups = [name.get(g, g) for g in groups]
    D = pd.read_csv(f"{a.scores}.{a.level}.specific.tsv.gz", sep="\t")
    D = D[(D.fdr <= 0.10) & D.group.isin(groups)]
    called = D.groupby("group")["se"].apply(set).to_dict()

    cat = {}
    with gzip.open(a.catalog, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            cat[f[3]] = (f[0], int(f[1]), int(f[2]))
    gc = load_gene_coords(None, cache_path=cache_path("gene_coords.GRCh38.106.tsv"))
    near = {}
    for gene in {g for v in ids.values() for g in v}:
        c, s, e = gc[gene]
        near[gene] = {se for se, (cc, ss, ee) in cat.items() if cc == c and ee >= s - a.window and ss <= e + a.window}

    rows, off_hit, off_n = [], 0, 0
    for grp, genes in ids.items():
        if grp not in groups:
            print(f"  [skip] {grp}: not a scored group", file=sys.stderr)
            continue
        n = int(H[(H.level == a.level) & (H.group == grp)].n_lines.iloc[0])
        for gene in genes:
            own = bool(near[gene] & called.get(grp, set()))
            others = [g for g in groups if g != grp]
            hits = sum(bool(near[gene] & called.get(g, set())) for g in others)
            off_hit += hits; off_n += len(others)
            rows.append(dict(group=grp, n=n, gene=gene, passes_own=own, passes_other=f"{hits}/{len(others)}"))
    R = pd.DataFrame(rows)
    on = int(R.passes_own.sum())
    print(R.to_string(index=False))
    odds, p = fisher_exact([[on, len(R) - on], [off_hit, off_n - off_hit]])
    print(f"\n{a.level}: own group {on}/{len(R)} = {100 * on / len(R):.0f}%  vs  other groups "
          f"{off_hit}/{off_n} = {100 * off_hit / max(off_n, 1):.1f}%   (OR {odds:.1f}, Fisher p = {p:.1e})")
    print(f"groups with >= 1 call: {sum(1 for g in groups if called.get(g))}/{len(groups)}; "
          f"total calls {len(D):,}")


if __name__ == "__main__":
    main()
