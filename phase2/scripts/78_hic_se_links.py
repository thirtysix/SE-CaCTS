#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Hi-C super-enhancer -> gene links from the 14 deep ENCODE cancer-line Hi-C maps (the sibling hi-c_analyses project,
pilot P1), reduced to one row per SE for the master-regulator figure (77_stage_master_tfs.py).

Candidates are the genes whose TSS lies within 1 Mb of the SE (P1 pairs.parquet; SE-CaCTS v3.1 catalogue). Contact =
balanced Hi-C at 10 kb divided by that line's genome-wide expected contact at 100 kb (keeps the distance decay, scale-
free across lines), averaged over the lines; a TSS within one 10 kb bin of the SE cannot be measured and ranks first
(inf). The top gene is the highest contact, ties by distance. P1 found this beats the nearest TSS only for SEs > 50 kb
from any TSS, and lineage-matched lines no better than others, so it is shared 3D contact, not a lineage-specific loop.

  phase2/data/hic_se_links.tsv.gz      se, n_genes, top_gene, top_contact, top_kb, next_gene, next_contact,
                                       nearest_gene, nearest_kb (kb = SE edge to TSS)
  phase2/data/hic_se_links.lines.tsv   the Hi-C maps used (ENCODE experiment and file accessions)

  ~/miniconda3/envs/atac_hdac/bin/python phase2/scripts/78_hic_se_links.py --hic <hi-c_analyses dir>
"""
from __future__ import annotations

import argparse
import os
import warnings

import numpy as np
import pandas as pd

SECACTS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RES, REF_BP = 10_000, 100_000


def contact_scores(H):
    """P1 pairs (SE x TSS within 1 Mb) with `hic` = mean deep-line contact (inf within one bin, -inf unmeasured), and the lines."""
    pr = pd.read_parquet(f"{H}/results/p1/pairs.parquet")
    g = pd.read_csv(f"{H}/manifests/p3_groups.tsv", sep="\t")
    r = pd.read_csv(f"{H}/manifests/p3_retained.tsv", sep="\t")
    g = g[(g["include"] == 1) & g["acc"].isin(r.loc[r["map_res_kb"] == 5.0, "acc"])]      # P1 evaluate.py --deep-only
    C = []
    for acc in g["acc"]:
        p = pd.read_parquet(f"{H}/results/p1/contacts/{acc}.pairs.parquet")
        e = pd.read_parquet(f"{H}/results/p1/contacts/{acc}.expected.parquet")
        e = e[e["offset"] == REF_BP // RES]
        assert (p["se"].to_numpy() == pr["se"].to_numpy()).all() and (p["gene"].to_numpy() == pr["gene"].to_numpy()).all()
        C.append((p["bal"] / float((e["expected"] * e["n_valid"]).sum() / e["n_valid"].sum())).to_numpy())
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        s = np.nanmean(np.vstack(C), axis=0)
    b0, b1, bt = pr["se_start"] // RES, (pr["se_end"] - 1) // RES, pr["tss"] // RES
    s = np.where((bt >= b0 - 1) & (bt <= b1 + 1), np.inf, s)
    pr["hic"] = np.where(np.isnan(s), -np.inf, s)
    return pr, g


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--hic", required=True, help="the hi-c_analyses project directory")
    ap.add_argument("--out", default=os.path.join(SECACTS, "phase2", "data", "hic_se_links.tsv.gz"))
    a = ap.parse_args()
    pr, lines = contact_scores(a.hic)
    o = pr.sort_values(["se", "hic", "dist"], ascending=[True, False, True])
    o["k"] = o.groupby("se").cumcount()
    top, nxt = o[o["k"] == 0].set_index("se"), o[o["k"] == 1].set_index("se")
    near = pr.sort_values(["se", "dist"]).groupby("se").head(1).set_index("se")
    T = pd.DataFrame({"n_genes": pr.groupby("se").size(), "top_gene": top["gene"], "top_contact": top["hic"],
                      "top_kb": (top["dist"] / 1000).round(1), "next_gene": nxt["gene"], "next_contact": nxt["hic"],
                      "nearest_gene": near["gene"], "nearest_kb": (near["dist"] / 1000).round(1)})
    T[["top_contact", "next_contact"]] = T[["top_contact", "next_contact"]].round(4).replace(-np.inf, np.nan)
    T.index.name = "se"
    T.to_csv(a.out, sep="\t", compression="gzip")
    sel = pd.read_csv(f"{a.hic}/manifests/p3_encode_selection.tsv", sep="\t")[["experiment", "file"]]
    lines[["line", "lineage", "acc"]].merge(sel, left_on="acc", right_on="file").drop(columns="acc").to_csv(
        a.out.replace(".tsv.gz", ".lines.tsv"), sep="\t", index=False)
    print(f"[78] {len(T):,} SEs from {len(lines)} Hi-C maps; top contact gene = nearest TSS for "
          f"{100 * (T.top_gene == T.nearest_gene).mean():.1f}%; single-candidate SEs {int((T.n_genes == 1).sum()):,} -> {a.out}")


if __name__ == "__main__":
    main()
