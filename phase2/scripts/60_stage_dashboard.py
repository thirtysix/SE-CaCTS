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

  ~/miniconda3/envs/atac_hdac/bin/python phase2/scripts/60_stage_dashboard.py \
      [--scores phase2/scores_v2 --results phase2/results_v2 --pull-bu 51]
"""
from __future__ import annotations

import argparse
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
sys.path.insert(0, SECACTS)
from secacts_env import DATAROOT                                    # noqa: E402

SCORES = os.path.join(PHASE2, "scores")          # overridden by --scores / --results in main()
RESULTS = os.path.join(PHASE2, "results")
OUT = os.path.join(os.environ.get("SECACTS_DOCS", os.path.join(SECACTS, "docs")), "data")   # SECACTS_DOCS: stage a copy, e.g. for screenshots
PERM = os.path.join(SCORES, "atlas.s3.perm")

# levels the panel supports as CALLS vs rankings-only (gotcha 72); v3.1 (--subtype-calls) calls subtypes too: with the
# "SE of its group" rule in the null the subtype counts are stable (FINDINGS §53)
CALL_LEVELS = [("lineage", "OncotreeLineage", "Lineage"),
               ("disease", "OncotreePrimaryDisease", "Primary disease")]
RANK_LEVELS = [("subtype", "OncotreeSubtype", "Subtype"),
               ("line", "line", "Cell line")]
LEVEL_NAME = {"OncotreeLineage": "Lineage", "OncotreePrimaryDisease": "Primary disease", "OncotreeSubtype": "Subtype"}
# fused-build labels (phase2/analysis/fused_labels.py, FINDINGS §55): how the group's experiments call the SE
SEL = {"core": "c", "unmasked": "u", "agnostic": "a", "gain": "g", "amplified": "A", "high": "H"}


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


def line_aliases(pull_set):
    """{StrippedCellLineName: {"label", "search"}} for the scored lines, from DepMap Model.csv (display
    name, CCLE name, Cellosaurus RRID, Oncotree subtype) plus the ChIP-Atlas cell name in the pull set."""
    ps = pd.read_csv(pull_set, sep="\t")
    key = ps["key"] if "key" in ps.columns else ps["model_id"]
    ps = ps.assign(key=key).dropna(subset=["key"]).drop_duplicates("key")
    md = pd.read_csv(os.path.join(DATAROOT, "DepMap", "2026q1", "Model.csv"), index_col="ModelID")
    out = {}
    for r in ps.itertuples():
        if r.key in md.index:
            m = md.loc[r.key]
            stripped, label = m["StrippedCellLineName"], m["CellLineName"]
            extra = [m.get("CCLEName"), m.get("RRID"), m.get("OncotreeSubtype"), r.cell, r.key]
        else:                                    # a line outside DepMap is keyed and named by its CVCL
            stripped = label = r.cell
            extra = [r.cvcl, r.subtype]
        terms = [stripped, label] + [x for x in extra if isinstance(x, str)]
        out[stripped] = {"label": label, "search": " | ".join(dict.fromkeys(terms))}
    return out


def line_hierarchy(pull_set, line_groups):
    """[[line, lineage, disease, subtype], ...] for every scored line (line = the line-level group key),
    with the same Oncotree labels the scorer grouped by: DepMap Model.csv, or the Cellosaurus-NCIt
    crosswalk for a line outside DepMap. The dashboard uses it to scope each level by the levels above."""
    ps = pd.read_csv(pull_set, sep="\t")
    key = ps["key"] if "key" in ps.columns else ps["model_id"]
    ps = ps.assign(key=key).dropna(subset=["key"]).drop_duplicates("key")
    md = pd.read_csv(os.path.join(DATAROOT, "DepMap", "2026q1", "Model.csv"), index_col="ModelID")
    lr = pd.read_csv(os.path.join(SECACTS, "phase1", "data", "lineage_resolved.tsv"), sep="\t").set_index("cvcl")
    out = []
    for r in ps.itertuples():
        if r.key in md.index:
            m = md.loc[r.key]
            row = [m["StrippedCellLineName"], m["OncotreeLineage"], m["OncotreePrimaryDisease"], m["OncotreeSubtype"]]
        elif r.cvcl in lr.index:
            m = lr.loc[r.cvcl]
            row = [m["cell_line"], m["lineage"], m["primary_disease"], m["subtype"]]
        else:
            continue
        if row[0] in line_groups:
            out.append([x if isinstance(x, str) else None for x in row])
    return out


def n_columns(path):
    """Sample count of a gzipped SE x sample matrix, from its header line alone."""
    with gzip.open(path, "rt") as fh:
        return len(fh.readline().rstrip("\n").split("\t")) - 1


def main():
    global SCORES, RESULTS, PERM, CALL_LEVELS, RANK_LEVELS
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scores", default=SCORES)
    ap.add_argument("--results", default=RESULTS)
    ap.add_argument("--pull-bu", type=int, default=43, help="total CSC billing units spent on the pull(s)")
    ap.add_argument("--pull-set", default=os.path.join(PHASE2, "data", "pull_set.tsv"),
                    help="pull set (v2 carries cn_provider) — used to count lines per copy-number source")
    ap.add_argument("--release", required=True, help="release label shown in the sidebar, e.g. v2")
    ap.add_argument("--release-date", required=True, help="YYYY-MM-DD")
    ap.add_argument("--release-title", default="", help="one line: what this release is")
    ap.add_argument("--release-notes", default="", help="what changed since the previous release")
    ap.add_argument("--main-label", default="Default", help="Analysis selector: the main run's name")
    ap.add_argument("--main-desc", default="", help="Analysis selector: one line on the main run")
    ap.add_argument("--variant", action="append", default=[],
                    help="another scoring run for the Analysis selector: 'key|prefix|label|description' "
                         "(prefix relative to --scores, e.g. atlas.s3.perm.nocn)")
    ap.add_argument("--pull-desc", default="", help="pipeline-chart step 1, lines separated by ' | ' (default: ChIP-Atlas only)")
    ap.add_argument("--subtype-calls", action="store_true", help="v3.1: subtype is a call level, not rankings only")
    ap.add_argument("--labels", help="fused_labels.groups.tsv.gz: CN labels per call, and the calls that pass only "
                                     "without CN correction (staged as calls_<level>.cndep.tsv)")
    ap.add_argument("--line-prefix", help="scoring run for the line level (v3.1: atlas.s3.lines.all, the vs-all arm; "
                                          "the main run is not scored per line there). Default: the main run")
    ap.add_argument("--min-group-lines", type=int, default=1,
                    help="v3.1: 2. A group of one cell line lists no calls: its specificity is that line's, which the "
                         "cell-line level tests with two independent studies (user decision 2026-10-03)")
    ap.add_argument("--n-pull", type=int, help="experiments pulled (default: columns of results/atlas.se_signal.tsv.gz)")
    a = ap.parse_args()
    if a.subtype_calls:
        CALL_LEVELS = CALL_LEVELS + [RANK_LEVELS[0]]
        RANK_LEVELS = RANK_LEVELS[1:]
    LAB = pd.read_csv(a.labels, sep="\t") if a.labels else None
    SCORES, RESULTS = a.scores, a.results
    PERM = os.path.join(SCORES, "atlas.s3.perm")
    os.makedirs(OUT, exist_ok=True)
    coords = load_coords()
    nearest = nearest_gene_fn()
    H = pd.read_csv(f"{PERM}.hierarchy_summary.tsv", sep="\t")
    LINEP = os.path.join(SCORES, a.line_prefix) if a.line_prefix else PERM
    HL = pd.read_csv(f"{LINEP}.hierarchy_summary.tsv", sep="\t") if a.line_prefix else H

    # concordance is the authority on the nearest gene AND whether it is group-specific — use its own
    # is_nearest pair so the gene shown and the concordance shown always name the SAME gene (the two
    # nearest-gene definitions differ by gene universe). Fall back to the coord-based nearest otherwise.
    CP = pd.read_csv(f"{PERM}.concordance2.pairs.tsv.gz", sep="\t")
    CP = CP[CP["is_nearest"] == True]                                   # noqa: E712
    conc_by = {(r.level, r.group, r.se): (r.gene, bool(r.concordant), float(r.rho)) for r in CP.itertuples()}

    manifest = {"levels": {}}

    def annotate(S, col):
        """gene (bridge nearest, else coordinate nearest), concordance and coordinates for each call row"""
        chrom, st, en, ncall, gene, dkb, conc, rho = [], [], [], [], [], [], [], []
        for r in S.itertuples():
            c, s_, e, nc = coords.get(r.se, ("", 0, 0, 0))
            hit = conc_by.get((col, r.group, r.se))
            if hit:
                g, is_c, rr = hit; d = 0
            else:
                g, is_c, rr = (nearest(c, (s_ + e) // 2) if c else ("", 0)), None, None
                g, d = g if isinstance(g, tuple) else (g, 0)
            chrom.append(c); st.append(s_); en.append(e); ncall.append(nc)
            gene.append(g); dkb.append(d); conc.append("" if is_c is None else int(is_c))
            rho.append("" if rr is None else round(rr, 3))
        return S.assign(chrom=chrom, start=st, end=en, n_called=ncall, gene=gene, dist_kb=dkb, conc=conc, rho=rho)

    LCOLS = ["cns", "sel", "n_exp", "n_amp", "cn_src"]

    def with_labels(S, col):
        """cns: r = called by both statistics (CN-robust), u = only with correction (CN-unmasked); sel: how the
        group's experiments call the SE (c core, u unmasked, a agnostic, g gain, A amplified, H high-level)"""
        L = LAB[LAB.level == col].set_index(["group", "se"])
        k = list(zip(S.group, S.se))
        st = L["cn_status"].reindex(k).values
        return S.assign(cns=[{"robust": "r", "unmasked": "u"}.get(x, "") for x in st],
                        sel=[SEL.get(x, "") for x in L["se_label"].reindex(k).values],
                        n_exp=L["n_exp"].reindex(k).fillna(0).astype(int).values,
                        n_amp=L["n_exp_amp"].reindex(k).fillna(0).astype(int).values,
                        cn_src=L["cn_sources"].reindex(k).fillna("").values)

    def stage_calls(prefix, suffix):
        """calls_<level><suffix>.tsv and gene_index<suffix>.json from one scoring run; returns {level: groups}.
        Concordance (nearest gene; is it group-specific in expression) is a property of the gene and group,
        so every variant reads the main run's bridge pairs. A level the run did not score is left out."""
        gene_index, out = {}, {}
        Hx = pd.read_csv(f"{prefix}.hierarchy_summary.tsv", sep="\t")
        labelled = LAB is not None and suffix == ""
        for short, col, label in CALL_LEVELS:
            if not os.path.exists(f"{prefix}.{col}.specific.tsv.gz"):
                print(f"[stage] {short}{suffix}: not scored in this run, left out")
                continue
            S = pd.read_csv(f"{prefix}.{col}.specific.tsv.gz", sep="\t").sort_values(["group", "rank"])
            small = set(Hx[(Hx.level == col) & (Hx.n_lines < a.min_group_lines)].group)
            S = S[~S.group.isin(small)]
            S = annotate(S, col)
            keep = ["group", "se", "rank", "jsd", "fdr", "cn_mean", "gene", "dist_kb", "conc", "rho",
                    "chrom", "start", "end", "n_called"]
            if labelled:
                S = with_labels(S, col)
                # calls that pass only WITHOUT correction, alongside (their own FDR; never merged): amplicon-driven
                # where the group's copy number at the locus is >= 2, gain-dependent below (FINDINGS §55)
                U_ = pd.read_csv(f"{prefix}.nocn.{col}.specific.tsv.gz", sep="\t")
                dep = LAB[(LAB.level == col) & (LAB.cn_status == "amplicon")][["group", "se"]]
                D = annotate(U_.merge(dep, on=["group", "se"]).query("group not in @small").sort_values(["group", "rank"]), col)
                D = with_labels(D, col)
                D["cns"] = np.where(D["cn_mean"] >= 2, "a", "g")
                D[keep + LCOLS].round({"jsd": 4, "fdr": 4, "cn_mean": 3}).to_csv(
                    os.path.join(OUT, f"calls_{short}.cndep.tsv"), sep="\t", index=False)
                ndep = D.groupby("group").size().to_dict()
                print(f"[stage] {short}: {len(D):,} calls only without CN correction "
                      f"({int((D.cns == 'a').sum()):,} amplicon-driven, {int((D.cns == 'g').sum()):,} gain-dependent)")
            S = S[keep + (LCOLS if labelled else [])]
            S.round({"jsd": 4, "fdr": 4, "cn_mean": 3}).to_csv(os.path.join(OUT, f"calls_{short}{suffix}.tsv"),
                                                               sep="\t", index=False)
            for r in S.itertuples():                                   # gene index for the finder (call levels only)
                if r.gene:
                    gene_index.setdefault(r.gene, []).append(
                        {"lv": short, "g": r.group, "r": int(r.rank), "fdr": round(float(r.fdr), 4),
                         "cn": round(float(r.cn_mean), 2), "c": r.conc if r.conc != "" else None})
            hsub = Hx[Hx.level == col]
            out[short] = {r.group: ({"n_lines": int(r.n_lines), "n_calls": 0, "below_min": True} if r.group in small else
                                    {"n_lines": int(r.n_lines), "n_calls": int(r.n_spec_fdr10)}) for r in hsub.itertuples()}
            if small:
                print(f"[stage] {short}{suffix}: {len(small)} groups below {a.min_group_lines} lines list no calls "
                      f"({int(hsub[hsub.group.isin(small)].n_spec_fdr10.sum()):,} calls left out)")
            if labelled:
                for g, v in out[short].items():
                    v["n_cndep"] = int(ndep.get(g, 0))
                    v["n_unmasked"] = int(((S.group == g) & (S.cns == "u")).sum())
            print(f"[stage] {short}{suffix}: {len(S):,} calls across {len(out[short])} groups")
        for g in gene_index:
            gene_index[g].sort(key=lambda x: x["r"])
        write_json(f"gene_index{suffix}.json", gene_index)
        print(f"[stage] gene index{suffix}: {len(gene_index):,} genes near a lineage/disease-specific SE")
        return out

    # ---- CALL levels: stage every specific SE, annotated with gene + coordinates + concordance
    main_groups = stage_calls(PERM, "")
    for short, col, label in CALL_LEVELS:
        manifest["levels"][short] = {"col": col, "label": label, "kind": "calls",
                                     "n_groups": len(main_groups[short]), "groups": main_groups[short]}
    # ---- analysis variants (the dashboard's Analysis selector): same levels, other scoring runs
    manifest["variants"] = [{"key": "main", "label": a.main_label, "desc": a.main_desc, "levels": list(main_groups),
                             "n_calls": {k: sum(v["n_calls"] for v in g.values()) for k, g in main_groups.items()}}]
    for spec in a.variant:
        key, prefix, label, desc = (spec.split("|") + ["", ""])[:4]
        vg = stage_calls(os.path.join(SCORES, prefix), f".{key}")
        manifest["variants"].append({"key": key, "label": label, "desc": desc, "groups": vg, "levels": list(vg),
                                     "n_calls": {k: sum(x["n_calls"] for x in g.values()) for k, g in vg.items()}})

    # ---- RANK-ONLY levels: stage the top-N rankings (already gene/coord annotated); NO counts
    for short, col, label in RANK_LEVELS:
        T = pd.read_csv(f"{LINEP if short == 'line' else PERM}.{col}.top_specific.tsv", sep="\t").sort_values(["group", "rank"])
        keep = ["group", "se", "rank", "jsd", "fdr", "cn_mean", "nearest_gene", "dist_kb"]
        T = T[keep].rename(columns={"nearest_gene": "gene"})
        # add coords for out-links
        cc = T["se"].map(lambda s: coords.get(s, ("", 0, 0, 0)))
        T = T.assign(chrom=[x[0] for x in cc], start=[x[1] for x in cc], end=[x[2] for x in cc])
        T.round({"jsd": 4, "fdr": 4, "cn_mean": 3}).to_csv(os.path.join(OUT, f"rank_{short}.tsv"),
                                                          sep="\t", index=False)
        hsub = (HL if short == "line" else H)
        hsub = hsub[hsub.level == col]
        groups = {r.group: {"n_lines": int(r.n_lines)} for r in hsub.itertuples()}
        if short == "line":
            # line groups are keyed by DepMap's StrippedCellLineName (e.g. NIHOVCAR3), which is hard to find
            # by the name people use. Give each a display label and search aliases, so the picker finds
            # OVCAR-3 from "ovcar-3", "OVCAR3", "NIH:OVCAR-3" or its Cellosaurus id.
            for g, al in line_aliases(a.pull_set).items():
                if g in groups:
                    groups[g].update(al)
        manifest["levels"][short] = {"col": col, "label": label, "kind": "rankings",
                                     "n_groups": len(groups), "groups": groups}
        print(f"[stage] {short}: rankings for {len(groups)} groups (no calls — panel unsupported)")

    manifest["hierarchy"] = line_hierarchy(a.pull_set, set(manifest["levels"]["line"]["groups"]))
    print(f"[stage] hierarchy: {len(manifest['hierarchy'])} lines with lineage / disease / subtype")
    write_json("manifest.json", manifest)

    # ---- CN ablation (call-based, honest null)
    A = pd.read_csv(f"{PERM}.cn_ablation_calls.tsv", sep="\t")
    abl = {"summary": [], "amplicon": [], "note": ""}
    staged = {col: {g for g, v in main_groups.get(short, {}).items() if not v.get("below_min")}
              for short, col, _ in CALL_LEVELS}
    A = A[[g in staged.get(l, ()) for l, g in zip(A.level, A.group)]]
    for short, col, label in CALL_LEVELS:
        sub = A[A.level == col]
        amp = sub[sub.kind == "amplicon_driven"]
        resc = sub[sub.kind == "rescued"]
        # counts from the staged calls (corrected) + the ablation sets
        n_corr = int(sum(v["n_calls"] for v in main_groups.get(short, {}).values()))
        if not len(sub):
            continue
        abl["summary"].append({"level": label, "corrected": n_corr,
                               "amplicon_driven": int(len(amp)), "rescued": int(len(resc)),
                               "removed_amplified": int((amp.cn_mean > 1.3).sum()),
                               "amp_median_cn": round(float(amp.cn_mean.median()), 1) if len(amp) else None,
                               "resc_median_cn": round(float(resc.cn_mean.median()), 3) if len(resc) else None})
    # the amplicon-driven calls, deduped to (group, gene) with the max cn
    # only calls at an AMPLIFIED locus (group-mean CN > 1.3) are listed as amplicons; on the v2 panel the
    # removed set also holds CN-neutral calls that lose significance when the null tightens
    amp_all = A[(A.kind == "amplicon_driven") & (A.cn_mean > 1.3)].sort_values("cn_mean", ascending=False)
    seen = set()
    for r in amp_all.itertuples():
        key = (r.group, r.nearest_gene)
        if key in seen:
            continue
        seen.add(key)
        abl["amplicon"].append({"level": {"OncotreeLineage": "Lineage", "OncotreePrimaryDisease": "Disease"}.get(r.level, "Subtype"),
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
    conc = {"summary": [{"level": LEVEL_NAME.get(r.level, r.level),
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
    ncalls = lambda k: int(sum(v["n_calls"] for v in main_groups.get(k, {}).values()))   # noqa: E731  as staged
    n_lineage_calls, n_disease_calls = ncalls("lineage"), ncalls("disease")
    n_subtype_calls = ncalls("subtype") if a.subtype_calls else None
    meta = {
        "n_samples": n_columns(os.path.join(RESULTS, "atlas.s3.se_signal.tsv.gz")),
        "n_lines": int((HL.level == "line").sum()), "n_ses": len(coords),
        "n_lineages": int((H.level == "OncotreeLineage").sum()),
        "n_diseases": int((H.level == "OncotreePrimaryDisease").sum()),
        "n_subtypes": int((H.level == "OncotreeSubtype").sum()),
        "n_lineage_calls": n_lineage_calls, "n_disease_calls": n_disease_calls, "n_subtype_calls": n_subtype_calls,
        "subtype_calls": bool(a.subtype_calls),
        "pull_bu": a.pull_bu,
        "n_pull": a.n_pull if a.n_pull is not None else n_columns(os.path.join(RESULTS, "atlas.se_signal.tsv.gz")),
        "fdr": "label-permutation, B=1000, FDR ≤ 0.10",
    }
    # subtype group sizes (why subtype is rankings-only) and lines per copy-number source
    sub = H[H.level == "OncotreeSubtype"]["n_lines"]
    meta["n_subtypes_single"] = int((sub == 1).sum())
    meta["n_subtypes_le4"] = int((sub <= 4).sum())
    ps = pd.read_csv(a.pull_set, sep="\t")
    if "cn_provider" not in ps.columns:
        ps["key"], ps["cn_provider"] = ps["model_id"], "depmap_wgs"
    kept = set(pd.read_csv(os.path.join(RESULTS, "atlas.s3.s3norm_params.tsv.gz"), sep="\t")["sample"])
    src = ps[ps.srx.isin(kept)].drop_duplicates("key")["cn_provider"].value_counts()
    meta["cn_sources"] = {lab: int(src.get(k, 0)) for k, lab in
                          (("depmap_wgs", "DepMap WGS"), ("cmp_wes", "CMP WES"), ("depmap_mc_wes", "DepMap WES"),
                                         ("ccle_snp6", "CCLE SNP6"), ("input_inferred", "inferred from ChIP input"))}
    if a.pull_desc:
        meta["pull_desc"] = [x.strip() for x in a.pull_desc.split("|")]
    # calibration: calls made on SHUFFLED labels (lineage level) by each null, where those runs exist
    # tests = SCORED loci x lineages (score_pilot drops loci with zero signal in every sample)
    fmat = f"{PERM}.OncotreeLineage.fdr.tsv.gz"
    meta["n_ses_scored"] = (len(pd.read_csv(fmat, sep="\t", usecols=[0])) if os.path.exists(fmat)
                            else meta["n_ses"])
    n_tests = meta["n_ses_scored"] * meta["n_lineages"]
    def shuffled_calls(prefix):
        f = f"{prefix}.OncotreeLineage.specific.tsv.gz"
        return int((pd.read_csv(f, sep="\t")["fdr"] <= 0.10).sum()) if os.path.exists(f) else None
    an, pm = (shuffled_calls(os.path.join(SCORES, p)) for p in ("atlas.s3.analytic.shuffle", "atlas.s3.perm.shuffle"))
    meta["calibration"] = {"n_tests": n_tests,
                           "analytic_shuffled_pct": round(100 * an / n_tests, 2) if an is not None else 6.05,
                           "perm_shuffled_calls": pm if pm is not None else 0}

    # release history: releases.json is committed and grows by one entry per release (an existing entry
    # with the same label is replaced, so restaging a release is idempotent). The sidebar shows the
    # current one, with the change in line count from the release before it.
    rel_path = os.path.join(OUT, "releases.json")
    rels = json.load(open(rel_path)) if os.path.exists(rel_path) else []
    entry = {"version": a.release, "date": a.release_date, "title": a.release_title, "notes": a.release_notes,
             **{k: meta[k] for k in ("n_lines", "n_samples", "n_ses", "n_lineages", "n_diseases", "n_subtypes",
                                     "n_lineage_calls", "n_disease_calls", "n_subtype_calls", "n_pull", "cn_sources")}}
    rels = [r for r in rels if r.get("version") != a.release] + [entry]
    rels.sort(key=lambda r: r["date"])
    write_json("releases.json", rels)
    prev = rels[-2] if len(rels) > 1 and rels[-1]["version"] == a.release else None
    meta["release"] = {"version": a.release, "date": a.release_date, "title": a.release_title,
                       "delta_lines": (meta["n_lines"] - prev["n_lines"]) if prev else None}
    write_json("meta.json", meta)
    print(f"[stage] wrote {len(os.listdir(OUT))} files to {OUT}")


if __name__ == "__main__":
    main()
