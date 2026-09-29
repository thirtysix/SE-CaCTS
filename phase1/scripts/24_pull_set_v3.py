#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 1 — 24: the v3 pull set, baseline-only (ROADMAP A3, FINDINGS 2026-09-29 §29-30).

Rows: the v2 pull set, plus experiments that passed the v3 signal checks (v3, v3b, v3c, v3d
candidate_check.tsv, status pass), plus Stage B experiments published so far (--stageb-done). An experiment
stays only if its treatment label is untreated, control or unclear (atlas_treatment_labels.tsv, the
atlas_labels/labels_1*.tsv chunks, stageB_treatment_manual.tsv); perturbed, input and other-mark rows go.

CN provider per line, first match wins: a v3b identity join (v3b_candidates.tsv `key`: the same individual as
a scored line, e.g. BE(2)-C -> SK-N-BE(2)) takes that line's key and CN; the line's v2 provider; its phase-1
admission (cn_admissions.tsv:
relatives, CCLE SNP6, input-inferred); else a track that actually builds, in the v2 order DepMap WGS ->
CMP WES -> DepMap MC_WES -> CCLE SNP6 (script 15). Lines with no CN are left out and listed.

Outputs: phase2/data/pull_set.v3.tsv (v2 columns + cn_cvcl, source), pull_srx.v3.txt,
v3_inferred_keys.txt (keys scored on input-inferred CN, for the sensitivity arm).
Reproduce:  python3 phase1/scripts/24_pull_set_v3.py --stageb-done <file of published Stage B SRX>
"""
import argparse
import glob
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "cnrose"))
from secacts_env import DATAROOT, cache_path                      # noqa: E402
from cnrose.cn.cmp import CellModelPassportsWesCN                  # noqa: E402
from cnrose.cn.depmap import DepMapMcWesCN, load_gene_coords       # noqa: E402

P1 = os.path.join(ROOT, "phase1", "data")
P2 = os.path.join(ROOT, "phase2", "data")
OUT = os.path.join(ROOT, "phase2", "analysis", "out")
DEPMAP = os.path.join(DATAROOT, "DepMap", "2026q1")
CMP = os.path.join(DATAROOT, "CellModelPassports")
GTF = os.path.join(DATAROOT, "0.human_genome", "Homo_sapiens.GRCh38.106.chr.gtf.gz")
KEEP = {"untreated", "control", "unclear"}
# LNCaP (parental) has no DepMap model; LNCaP clone FGC is the same individual (as for C4-2 in script 23)
JOINS = {"CVCL_0395": "CVCL_1379"}


def labels():
    lab = pd.read_csv(os.path.join(OUT, "atlas_treatment_labels.tsv"), sep="\t", usecols=["srx", "class"])
    extra = [pd.read_csv(p, sep="\t", usecols=["srx", "class"])
             for p in sorted(glob.glob(os.path.join(OUT, "atlas_labels", "labels_1[0-9].tsv")))]
    sb = pd.read_csv(os.path.join(OUT, "stageB_treatment_manual.tsv"), sep="\t", usecols=["srx", "class"])
    bad = pd.read_csv(os.path.join(OUT, "stageB_not_h3k27ac.txt"), header=None, names=["srx"])
    bad["class"] = "not_h3k27ac"
    allab = pd.concat([lab, *extra, sb, bad]).drop_duplicates("srx", keep="last")
    return allab.set_index("srx")["class"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stageb-done", help="file with one published Stage B SRX per line")
    a = ap.parse_args()
    cls = labels()
    keep = lambda s: cls.get(s) in KEEP                              # noqa: E731

    v2 = pd.read_csv(os.path.join(P2, "pull_set.v2.tsv"), sep="\t")
    v2["cn_cvcl"], v2["source"] = pd.NA, "v2"
    unl = ~v2["srx"].isin(cls.index)
    print(f"[24] v2 pull set {len(v2)}: {int(unl.sum())} unlabelled (dropped), "
          f"{int(v2['srx'].map(keep).sum())} baseline", file=sys.stderr)
    v2 = v2[v2["srx"].map(keep)]

    cand = []
    for name in ("v3", "v3b", "v3c", "v3d"):
        c = pd.read_csv(os.path.join(OUT, f"{name}_candidate_check.tsv"), sep="\t")
        c = c[c["status"] == "pass"][["srx", "cvcl", "cell_line"]].assign(source=name)
        cand.append(c)
    if a.stageb_done:
        done = set(open(a.stageb_done).read().split())
        sb = pd.read_csv(os.path.join(ROOT, "phase2", "stageB", "manifest.baseline.tsv"), sep="\t")
        cand.append(sb[sb["srx"].isin(done)][["srx", "cvcl", "cell_line"]].assign(source="stageB"))
    cand = pd.concat(cand).drop_duplicates("srx")
    cand = cand[~cand["srx"].isin(v2["srx"])]
    n0 = len(cand)
    cand = cand[cand["srx"].map(keep)]
    print(f"[24] candidates passing checks, not in v2: {n0}; baseline {len(cand)}", file=sys.stderr)

    # line metadata: v2 provider first, then phase-1 admissions, then resolve
    v2line = pd.read_csv(os.path.join(P2, "pull_set.v2.tsv"), sep="\t").drop_duplicates("cvcl").set_index("cvcl")
    adm = pd.read_csv(os.path.join(OUT, "cn_admissions.tsv"), sep="\t").drop_duplicates("cvcl").set_index("cvcl")
    model = pd.read_csv(os.path.join(DEPMAP, "Model.csv"), usecols=["ModelID", "RRID", "OncotreeLineage",
                                                                      "OncotreeSubtype"])
    mid_of = model.dropna(subset=["RRID"]).drop_duplicates("RRID").set_index("RRID")["ModelID"]
    lin = pd.read_csv(os.path.join(P1, "lineage_resolved.tsv"), sep="\t").set_index("cvcl")
    v2key = v2line.reset_index().drop_duplicates("key").set_index("key")
    v3bkey = pd.read_csv(os.path.join(OUT, "v3b_candidates.tsv"), sep="\t").dropna(subset=["key"]).set_index("srx")["key"]
    # a join the v3b signal check accepted holds for the line, so Stage B experiments on it join too
    v3bc = pd.read_csv(os.path.join(OUT, "v3b_candidate_check.tsv"), sep="\t")
    ok_cvcl = set(v3bc.loc[v3bc["status"] == "pass", "cvcl"])
    join_of = {c: k for c, k in pd.read_csv(os.path.join(OUT, "v3b_candidates.tsv"), sep="\t").dropna(subset=["key"])
               [["cvcl", "key"]].itertuples(index=False) if c in ok_cvcl}
    need = sorted((set(cand["cvcl"]) - set(v2line.index) - set(adm.index)) | set(JOINS.values()))
    mids = {c: mid_of.get(c) for c in need}
    wgs_ids = set(pd.read_csv(os.path.join(DEPMAP, "OmicsCNGeneWGS.csv"), usecols=[0]).iloc[:, 0])
    cmp = CellModelPassportsWesCN(os.path.join(CMP, "WES_pureCN_CNV_genes_latest.csv.gz"),
                                  os.path.join(CMP, "model_list_20240110.csv"), cache_dir=cache_path("cmp_wes"))
    cmp.preload(need)
    gc = load_gene_coords(GTF, cache_path=cache_path("gene_coords.GRCh38.106.tsv"))
    mcw = DepMapMcWesCN(os.path.join(DEPMAP, "OmicsCNGeneMC_WES.csv"), os.path.join(DEPMAP, "ModelCondition.csv"), gc)
    mcw.preload([m for m in mids.values() if isinstance(m, str)])
    ccle_keys = set(pd.read_csv(os.path.join(P2, "cn_ccle_snp6.hg38.tsv.gz"), sep="\t", usecols=["key"])["key"])

    def resolve(c):
        m = mids.get(c)
        m = m if isinstance(m, str) else None
        key = m or c
        if m and m in wgs_ids:
            return key, m, "depmap_wgs"
        if cmp.track(c) is not None:
            return key, m, "cmp_wes"
        if m and mcw.track(m) is not None:
            return key, m, "depmap_mc_wes"
        if key in ccle_keys:
            return key, m, "ccle_snp6"
        return None

    rows, nocn = [], {}
    for r in cand.itertuples(index=False):
        c = r.cvcl
        jk = v3bkey.get(r.srx) if r.source == "v3b" else (join_of.get(c) if r.source == "stageB" else None)
        if isinstance(jk, str) and jk in v2key.index:
            v = v2key.loc[jk]
            key, mid, prov, cn_cvcl = jk, v["model_id"], v["cn_provider"], v["cvcl"]
        elif c in JOINS and resolve(JOINS[c]) is not None:
            key, mid, prov = resolve(JOINS[c])
            cn_cvcl = JOINS[c]
        elif c in v2line.index:
            v = v2line.loc[c]
            key, mid, prov, cn_cvcl = v["key"], v["model_id"], v["cn_provider"], pd.NA
        elif c in adm.index:
            d = adm.loc[c]
            key, prov = d["key"], d["cn_provider"]
            mid = key if str(key).startswith("ACH-") else pd.NA
            cn_cvcl = d["cn_cvcl"] if isinstance(d["cn_cvcl"], str) and d["cn_cvcl"] else pd.NA
        else:
            got = resolve(c)
            if got is None:
                nocn[c] = r.cell_line
                continue
            key, mid, prov = got
            cn_cvcl = pd.NA
        lm = model.set_index("ModelID").reindex([mid]).iloc[0] if isinstance(mid, str) else None
        rows.append(dict(srx=r.srx, cell=r.cell_line, cvcl=c, model_id=mid,
                         lineage=(lm["OncotreeLineage"] if lm is not None else lin["lineage"].get(c)),
                         subtype=(lm["OncotreeSubtype"] if lm is not None else lin["subtype"].get(c)),
                         key=key, cn_provider=prov, cn_cvcl=cn_cvcl, source=r.source))
    new = pd.DataFrame(rows)
    v3 = pd.concat([v2, new], ignore_index=True)
    v3.to_csv(os.path.join(P2, "pull_set.v3.tsv"), sep="\t", index=False)
    v3["srx"].sort_values().to_csv(os.path.join(P2, "pull_srx.v3.txt"), index=False, header=False)
    inf = sorted(set(v3.loc[v3["cn_provider"] == "input_inferred", "key"]))
    open(os.path.join(P2, "v3_inferred_keys.txt"), "w").write(",".join(inf) + "\n")

    print(f"[24] pull_set.v3: {len(v3)} experiments on {v3['key'].nunique()} keys "
          f"(v2 {len(v2)}, new {len(new)}; by source {new['source'].value_counts().to_dict() if len(new) else {}})")
    print(f"[24] CN providers (lines): {v3.drop_duplicates('key')['cn_provider'].value_counts().to_dict()}")
    print(f"[24] new keys vs v2: {len(set(v3['key']) - set(v2line['key']))}; input-inferred keys {len(inf)}")
    print(f"[24] left out, no CN: {len(nocn)} lines: {', '.join(sorted(nocn.values()))}")


if __name__ == "__main__":
    main()
