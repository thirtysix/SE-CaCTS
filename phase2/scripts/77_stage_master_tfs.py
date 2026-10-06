#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Known lineage master regulators against the atlas (the poster's Fig. 5, redrawn per release), for the dashboard
and the README.

The 22 transcription factors chosen from the literature before scoring (poster_figures.IDENTITY). For each, the
super-enhancer best in the gene's own lineage (lowest FDR, ties by JSD; poster_figures.identity_table says why not the
lowest JSD) among those within 100 kb of the gene and those linked to it by Hi-C (78_hic_se_links.py: the gene is the
SE's top contact among >= --hic-min-genes genes with a TSS within 1 Mb; with one candidate the "link" is forced, e.g.
GATA3's gene desert), with its JSD rank, permutation FDR and CN status in every lineage. The same at primary-disease
and subtype level: columns are the groups of >= 2 lines in the figure's 12 lineages (a group of one line has no calls),
and each gene's SE is the best over its lineage's groups. Each level answers its own question: a subtype is scored
against every other subtype, so a program its sibling subtypes share (PAX8, SOX17 in ovary) is not subtype-specific.

  <docs>/data/master_tfs.json        {release, n_ses, hic: {maps, min_genes} or null, levels: [{key, label,
                                      cols: [[group, lineage, n_lines]], own: [genes called in >= 1 own group, of],
                                      other: [other-lineage cells called, of],
                                      rows: [{gene, lineage, se, locus, own: [col], cells: [[rank, fdr, cn], ...], hic?}]}]}
                                      cn = "r" CN-robust / "u" CN-unmasked, on called cells (FDR <= 0.10) only;
                                      hic = {kb, n, c, next, nc} when the SE is linked by Hi-C, not within 100 kb
  <fig-out>/master-regulators{,-disease,-subtype}.{png,svg,pdf}  the static figures; without --fig-out the JSON only

  ~/miniconda3/envs/atac_hdac/bin/python phase2/scripts/77_stage_master_tfs.py --docs docs \
      --scores phase2/scores_v31f --results phase2/results_v31f --fig-out assets/figures
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys

SECACTS = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
FDR = 0.10
# key (the card's), atlas level, Oncotree column, label, figure name, figure geometry; lit = the literature table's level
# (disease / subtype master TFs, each gene owning its own group) instead of the lineage TFs over the lineage's groups
LEVELS = [dict(key="lineage", lv="lineage", col="OncotreeLineage", label="Lineage", fig="master-regulators", geom={}),
          dict(key="disease", lv="disease", col="OncotreePrimaryDisease", label="Primary disease",
               fig="master-regulators-disease", geom=dict(width=300)),
          dict(key="subtype", lv="subtype", col="OncotreeSubtype", label="Subtype", fig="master-regulators-subtype",
               geom=dict(width=470, cell_fs=9, rot=55)),
          dict(key="disease_tfs", lv="disease", col="OncotreePrimaryDisease", label="Disease TFs", lit="disease",
               fig="master-regulators-disease-tfs", geom=dict(width=430, cell_fs=9, rot=55)),
          dict(key="subtype_tfs", lv="subtype", col="OncotreeSubtype", label="Subtype TFs", lit="subtype",
               fig="master-regulators-subtype-tfs", geom=dict(width=470, cell_fs=9, rot=55))]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--docs", required=True)
    ap.add_argument("--scores", required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--labels", help="fused_labels.groups table (default <scores>/fused_labels.groups.tsv.gz if present)")
    ap.add_argument("--hic-links", default=os.path.join(SECACTS, "phase2", "data", "hic_se_links.tsv.gz"),
                    help="78_hic_se_links.py output; 'none' = the 100 kb rule only")
    ap.add_argument("--hic-min-genes", type=int, default=2)
    ap.add_argument("--lit", default=os.path.join(SECACTS, "phase2", "data", "master_tfs_literature.tsv"),
                    help="disease / subtype master TFs from the literature (level, group, gene, pmid, first_author_year, evidence)")
    ap.add_argument("--fig-out")
    a = ap.parse_args()
    labels = a.labels or os.path.join(a.scores, "fused_labels.groups.tsv.gz")
    labels = labels if os.path.exists(labels) else None

    spec = importlib.util.spec_from_file_location("poster_figures", os.path.join(SECACTS, "phase2", "figures", "poster_figures.py"))
    pf = importlib.util.module_from_spec(spec); spec.loader.exec_module(pf)
    import pandas as pd
    lit = pd.read_csv(a.lit, sep="\t") if os.path.exists(a.lit) else None
    genes = {g for _, gs in pf.IDENTITY for g in gs} | (set(lit["gene"]) if lit is not None else set())
    links, hic = None, None
    if a.hic_links != "none":
        T = pd.read_csv(a.hic_links, sep="\t")
        T = T[(T["n_genes"] >= a.hic_min_genes) & T["top_gene"].isin(genes)]
        links = {}
        for r in T.itertuples():
            links.setdefault(r.top_gene, {})[r.se] = dict(kb=r.top_kb, n=int(r.n_genes), c=r.top_contact, next=r.next_gene,
                                                          nc=r.next_contact)
        nmaps = len(pd.read_csv(a.hic_links.replace(".tsv.gz", ".lines.tsv"), sep="\t"))
        hic = dict(maps=nmaps, min_genes=a.hic_min_genes)
    n_ses = len(pd.read_csv(os.path.join(a.scores, "atlas.s3.perm.OncotreeLineage.fdr.tsv.gz"), sep="\t", usecols=[0]))
    meta = json.load(open(os.path.join(a.docs, "data", "meta.json")))
    man = json.load(open(os.path.join(a.docs, "data", "manifest.json")))
    H = pd.DataFrame(man["hierarchy"], columns=["line", "lineage", "disease", "subtype"])
    code = {"robust": "r", "unmasked": "u"}
    if a.fig_out:
        pf.OUT = a.fig_out
    levels = []
    lin_order = [ln for ln, _ in pf.IDENTITY] + sorted(set(H["lineage"]) - {ln for ln, _ in pf.IDENTITY})
    for L in LEVELS:
        short, lev, label, fname, fig = L["lv"], L["col"], L["label"], L["fig"], L["geom"]
        if "lit" in L and lit is None:
            continue
        groups, identity, blocks, ref = None, None, None, {}
        if short != "lineage":                  # groups that can carry calls (>= 2 lines)
            G = man["levels"][short]["groups"]
            lin_of = H.groupby(short)["lineage"].agg(lambda x: x.mode()[0])
        if "lit" in L:                          # each listed group owns its own column; bracketed by lineage
            T = lit[lit["level"] == L["lit"]]
            bad = sorted(set(T["group"]) - {g for g in G if G[g]["n_lines"] >= 2})
            if bad:
                sys.exit(f"[77] FAIL: literature groups not among the {short} groups of >= 2 lines: {bad}")
            gl = sorted(dict.fromkeys(T["group"]), key=lambda g: (lin_order.index(lin_of[g]), -G[g]["n_lines"], g))
            identity = [(g, list(dict.fromkeys(T.loc[T["group"] == g, "gene"]))) for g in gl]
            groups, blocks = [(g, g) for g in gl], {g: lin_of[g] for g in gl}
            ref = {(r.group, r.gene): [r.first_author_year, str(r.pmid)] for r in T.itertuples()}
        elif short != "lineage":                # the lineage TFs over their lineage's groups, in the figure's lineages
            groups = [(g, ln) for ln, _ in pf.IDENTITY for g in
                      sorted((g for g in G if lin_of.get(g) == ln and G[g]["n_lines"] >= 2), key=lambda g: (-G[g]["n_lines"], g))]
        cols, recs = pf.identity_table(a.scores, a.results, select="fdr", labels=labels, links=links, level=lev, groups=groups,
                                       identity=identity)
        lin = blocks or (dict(groups) if groups else {c: c for c in cols})
        own, other, rows = [0, 0], [0, 0], []
        for r in recs:
            cells = [[r["rank"][j], round(float(r["fdr"][j]), 4), code.get(r["cn"][j]) if r["fdr"][j] <= FDR else None]
                     for j in range(len(cols))]
            oj = [cols.index(o) for o in r["own"]]
            own[0] += any(cells[j][1] <= FDR for j in oj); own[1] += 1
            for j, c in enumerate(cols):
                if j not in oj:
                    other[0] += cells[j][1] <= FDR; other[1] += 1
            rows.append(dict(gene=r["gene"], lineage=r["lineage"], se=r["se"], own=oj,
                             locus=f"{r['chrom']}:{r['start'] + 1}-{r['end']}", cells=cells))
            if r["link"]:
                rows[-1]["hic"] = {k: (round(float(v), 4) if isinstance(v, float) else v) for k, v in r["link"].items()}
            if ref:
                rows[-1]["ref"] = ref[(r["lineage"], r["gene"])]
        n_lines = H.groupby(short).size() if short != "lineage" else H.groupby("lineage").size()
        levels.append(dict(key=L["key"], lv=short, label=label, kind="tfs" if "lit" in L else "lineage_tfs",
                           cols=[[c, lin[c], int(n_lines.get(c, 0))] for c in cols], own=own, other=other, rows=rows))
        print(f"[77] {L['key']}: {len(cols)} groups; {own[0]}/{own[1]} genes called in their own lineage's groups; other-lineage "
              f"cells {other[0]}/{other[1]} ({100 * other[0] / max(other[1], 1):.1f}%); missed "
              f"{[r['gene'] for r in rows if not any(r['cells'][j][1] <= FDR for j in r['own'])]}; Hi-C rows "
              f"{[(r['gene'], r['se']) for r in rows if 'hic' in r]}")
        if a.fig_out:
            pf.fig_identity(table=(cols, recs), name=fname, cn_mark=labels is not None, contrast_text=True,
                            groups=None if blocks else groups, blocks=blocks,
                            fdr_label=f"FDR ({'primary disease' if short == 'disease' else short})", **fig)
    out = dict(release=meta.get("release", {}).get("version"), n_ses=n_ses, hic=hic, levels=levels)
    with open(os.path.join(a.docs, "data", "master_tfs.json"), "w") as fh:
        json.dump(out, fh, separators=(",", ":"))


if __name__ == "__main__":
    main()
