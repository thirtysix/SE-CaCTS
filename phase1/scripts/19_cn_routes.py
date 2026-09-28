#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 1 — 19: which copy-number route could reach the lines held back for lack of CN (FINDINGS §21).

Lines = registry rows with decision "not pulled" (ChIP-Atlas H3K27ac on one of the 631 human cancer lines,
no CN track) plus the identity-gap lines. For each line, every route this project has not used yet:

  ccle2019       CCLE 2019 SNP6 segments (cBioPortal ccle_broad_2019, hg19, continuous log2), via DEPMAPID
  relative       a Cellosaurus parent/child (HI) or same-individual (OI) line with a CN source
  true identity  Cellosaurus says the line is misidentified/contaminated; the real line has a CN source
  depmap-no-cn   the line has a DepMap model but no CN in 26Q1 (a newer/older DepMap release might)
CN sources counted for relatives: DepMap WGS, DepMap MC_WES, CMP WES pureCN 2025, CCLE 2019 SNP6.

Output: phase2/analysis/out/cn_routes.tsv (one row per line, best route first) + a summary on stderr.
"""
import os
import re
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
from secacts_env import DATAROOT, cache_path                     # noqa: E402

P1 = os.path.join(ROOT, "phase1", "data")
OUT = os.path.join(ROOT, "phase2", "analysis", "out")
DEPMAP = os.path.join(DATAROOT, "DepMap", "2026q1")
CCLE = os.path.join(DATAROOT, "cancer_cellline_encyclopedia", "ccle_broad_2019")


def cellosaurus():
    rec, out = {}, {}
    with open(os.path.join(DATAROOT, "cellosaurus", "cellosaurus.txt"), encoding="utf-8", errors="replace") as fh:
        for line in fh:
            tag, val = line[:2], line[5:].rstrip("\n")
            if tag == "ID":
                rec = {"id": val, "hi": [], "oi": [], "prob": "", "prob_cvcl": [], "human": False}
            elif tag == "AC":
                rec["ac"] = val.strip()
            elif tag == "HI":
                rec["hi"].append(val.split("!")[0].strip())
            elif tag == "OI":
                rec["oi"].append(val.split("!")[0].strip())
            elif tag == "OX" and "NCBI_TaxID=9606" in val:
                rec["human"] = True
            elif tag == "CC" and val.startswith("Problematic cell line:"):
                rec["prob"] = val
                rec["prob_cvcl"] = re.findall(r"Cellosaurus=(CVCL_\w+)", val)
            elif tag == "//":
                if rec.get("human") and rec.get("ac"):
                    out[rec["ac"]] = rec
                rec = {}
    return out


def main():
    reg = pd.read_csv(os.path.join(P1, "h3k27ac_registry.tsv"), sep="\t", low_memory=False)
    held = reg[reg.decision == "not pulled"]
    idg = pd.read_csv(os.path.join(OUT, "identity_gap_lines.tsv"), sep="\t")
    lines = held.groupby(["cvcl", "cell_line", "lineage"]).srx.nunique().rename("n_exp").reset_index()
    print(f"[19] held back for lack of CN: {len(lines)} lines, {int(lines.n_exp.sum())} experiments "
          f"({len(set(lines.cvcl) & set(idg.cvcl))} already routed as identity-gap)", file=sys.stderr)

    m = pd.read_csv(os.path.join(DEPMAP, "Model.csv"), usecols=["ModelID", "RRID", "CellLineName"])
    rrid_of = dict(zip(m.ModelID, m.RRID))
    wgs = {rrid_of.get(x) for x in pd.read_csv(os.path.join(DEPMAP, "OmicsCNGeneWGS.csv"), usecols=["ModelID"]).ModelID}
    mc = pd.read_csv(os.path.join(DEPMAP, "ModelCondition.csv"), usecols=["ModelConditionID", "ModelID"])
    mcw_mc = set(pd.read_csv(os.path.join(DEPMAP, "OmicsCNGeneMC_WES.csv"), usecols=["ModelConditionID"]).ModelConditionID)
    mcw = {rrid_of.get(x) for x in mc[mc.ModelConditionID.isin(mcw_mc)].ModelID}
    cmp_ids = set(open(cache_path("cmp_wes_models.txt")).read().split())
    ml = pd.read_csv(os.path.join(DATAROOT, "CellModelPassports", "model_list_20240110.csv"), usecols=["model_id", "RRID"])
    cmp = set(ml[ml.model_id.isin(cmp_ids)].RRID.dropna())
    cs = pd.read_csv(os.path.join(CCLE, "data_clinical_sample.txt"), sep="\t", comment="#", usecols=["SAMPLE_ID", "DEPMAPID"])
    seg_ids = set(pd.read_csv(os.path.join(CCLE, "data_cna_hg19.seg"), sep="\t", usecols=["ID"]).ID)
    ccle = {rrid_of.get(d) for s, d in zip(cs.SAMPLE_ID, cs.DEPMAPID) if s in seg_ids}
    src = {"DepMap WGS": wgs, "DepMap MC_WES": mcw, "CMP WES": cmp, "CCLE 2019 SNP6": ccle}
    has = lambda c: [k for k, v in src.items() if c in v]
    dm_rrid = set(m.RRID.dropna())

    cel = cellosaurus()
    kids, same = {}, {}
    for ac, r in cel.items():
        for p in r["hi"]:
            kids.setdefault(p, []).append(ac)
        for o in r["oi"]:
            same.setdefault(o, set()).add(ac)
            same.setdefault(ac, set()).add(o)
    name = lambda c: cel[c]["id"] if c in cel else c

    rows = []
    for r in lines.itertuples():
        c = r.cvcl
        routes = []
        if c in ccle:
            routes.append(("ccle2019", "CCLE 2019 SNP6 (own)"))
        rec = cel.get(c, {})
        for t in rec.get("prob_cvcl", []):
            if has(t):
                routes.append(("true identity", f"{name(t)} {t}: {'/'.join(has(t))}"))
        rel = [(p, "parent") for p in rec.get("hi", [])] + [(k, "child") for k in kids.get(c, [])] + \
              [(o, "same individual") for o in same.get(c, set())]
        for t, how in rel:
            if has(t):
                routes.append(("relative", f"{how} {name(t)} {t}: {'/'.join(has(t))}"))
        if c in dm_rrid and not routes:
            routes.append(("depmap-no-cn", "DepMap model without 26Q1 CN"))
        rows.append(dict(cvcl=c, cell_line=r.cell_line, lineage=r.lineage, n_exp=r.n_exp,
                         identity_gap=c in set(idg.cvcl), problematic=rec.get("prob", "")[:160],
                         best_route=routes[0][0] if routes else "none",
                         routes=" | ".join(f"{a}: {b}" for a, b in routes)))
    t = pd.DataFrame(rows).sort_values("n_exp", ascending=False)
    t.to_csv(os.path.join(OUT, "cn_routes.tsv"), sep="\t", index=False)
    s = t.groupby("best_route").agg(lines=("cvcl", "size"), experiments=("n_exp", "sum"))
    print(s.to_string(), file=sys.stderr)
    print(t.head(40)[["cell_line", "lineage", "n_exp", "best_route", "routes"]].to_string(index=False), file=sys.stderr)


if __name__ == "__main__":
    main()
