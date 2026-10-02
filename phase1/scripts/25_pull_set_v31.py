#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 1 — 25: the v3.1 pull sets (ROADMAP "v3.1 scope, decided 2026-10-02"; FINDINGS §36, §47).

Baseline (`pull_set.v31.tsv`), from v3's pull set:
  1. `phase2/data/v31_label_overrides.tsv`: drop (50), admit (2), own_line (NCI/ADR-RES, its own key and CN, never
     pooled with OVCAR-8), keep (LN-229 tier-2 fallback);
  2. the v3.1 relabel of v3 "controls" in studies with treated arms (`out/v31_relabel/final_labels.tsv`): rows no
     longer untreated/control leave; tier-2 rows stay only where their line has no tier-1 experiment passing QC;
  3. EA.hy926 (HUVEC x A549 fusion hybrid, no lineage) leaves;
  4. identity problems leave (user, 2026-10-02): HSMM (ENCODE normal myoblasts mapped to the HS-MM sarcoma line),
     RMS (CVCL_W527, a catch-all for RD, Rh-4 and others), HNE-1 and HONE-1 (HeLa-contaminated), RERF-LC-OK
     (a Marcus derivative), PC-14 (a PC-9 derivative: the same individual twice), SK-HEP-1 (endothelial);
  5. CN source: a line with a DepMap WGS track (OmicsCNGeneWGS, default entry) is scored on it (the WGS-first
     rule; 33 lines were on CMP WES / DepMap WES / input inference, 3 of them through an ID-join bug). Their
     calling-time re-calls are listed in `recall_redo.v31.tsv` for reduce_v31.slurm to redo.
All experiments (`pull_set.v31.all.tsv`, variants 4-5): the baseline plus every pulled ChIP-Atlas experiment
labelled perturbed on a v3.1 line, if it passed its own checks (in the v2 pull set, or `status == pass` in its
candidate check). Inputs, other marks, out-of-scope rows and derived sublines (relabel note `derivative:`) stay out.

    python3 phase1/scripts/25_pull_set_v31.py
"""
import glob, os, sys
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
from secacts_env import DATAROOT                                  # noqa: E402
P2 = os.path.join(ROOT, "phase2", "data")
OUT = os.path.join(ROOT, "phase2", "analysis", "out")
KEEP = {"untreated", "control"}
DROP_LINES = {"EA.hy926", "EA-hy926", "EA.HY926"}
DROP_CVCL = {"CVCL_W472": "HSMM is normal myoblasts", "CVCL_W527": "RMS catch-all", "CVCL_0308": "HNE-1 HeLa",
             "CVCL_8706": "HONE-1 HeLa", "CVCL_3154": "RERF-LC-OK Marcus derivative",
             "CVCL_1640": "PC-14 = PC-9 derivative", "CVCL_0525": "SK-HEP-1 endothelial"}


def main():
    ps = pd.read_csv(os.path.join(P2, "pull_set.v3.tsv"), sep="\t")
    qc = set(open(os.path.join(P2, "qc_srx.v3.txt")).read().split())
    log = []
    n0 = len(ps)

    # 1. overrides
    ov = pd.read_csv(os.path.join(P2, "v31_label_overrides.tsv"), sep="\t")
    drop = set(ov.loc[ov.action == "drop", "srx"])
    log.append(("override_drop", int(ps.srx.isin(drop).sum())))
    ps = ps[~ps.srx.isin(drop)]
    lines = ps.drop_duplicates("cell").set_index("cell")
    adm = []
    for r in ov[ov.action == "admit"].itertuples():
        if r.cell in lines.index:
            row = lines.loc[r.cell].to_dict(); row.update(srx=r.srx, cell=r.cell, source="v31_admit"); adm.append(row)
        else:
            log.append((f"admit_skipped_no_line:{r.cell}", 1))
    own = ov[ov.action == "own_line"]
    if len(own):
        log.append(("own_line_left_out_no_CN", len(own)))   # NCI/ADR-RES: CMP lists it with no data; not in DepMap or CCLE
    ps = pd.concat([ps, pd.DataFrame(adm)], ignore_index=True)
    log.append(("override_admit", len(adm)))

    # 2. relabel
    F = pd.read_csv(os.path.join(OUT, "v31_relabel", "final_labels.tsv"), sep="\t", dtype=str)
    out_rel = set(F.loc[~F.final_class.isin(KEEP), "srx"])
    log.append(("relabel_not_baseline", int(ps.srx.isin(out_rel).sum())))
    ps = ps[~ps.srx.isin(out_rel)]
    tier2 = set(F.loc[F.final_class.isin(KEEP) & (F.final_tier == "2"), "srx"])
    keep_t2 = set(ov.loc[ov.action == "keep", "srx"])
    ps["tier"] = ps.srx.map(lambda s: 2 if (s in tier2 and s not in keep_t2) else 1)
    t1_qc = ps[(ps.tier == 1) & ps.srx.isin(qc)].groupby("key").size()
    fallback = ps.key.map(t1_qc).fillna(0).eq(0)
    drop_t2 = (ps.tier == 2) & ~fallback
    log.append(("tier2_dropped_line_has_tier1", int(drop_t2.sum())))
    log.append(("tier2_kept_fallback", int(((ps.tier == 2) & fallback).sum())))
    ps = ps[~drop_t2]

    # 3. EA.hy926
    ea = ps.cell.isin(DROP_LINES)
    log.append(("EA.hy926", int(ea.sum())))
    ps = ps[~ea]

    # 4. identity problems
    idp = ps.cvcl.isin(DROP_CVCL)
    log.append(("identity_drop " + ", ".join(sorted({DROP_CVCL[c] for c in ps.loc[idp, "cvcl"]})), int(idp.sum())))
    ps = ps[~idp]

    # 5. WGS first
    D = os.path.join(DATAROOT, "DepMap", "2026q1")
    wc = pd.read_csv(os.path.join(D, "OmicsCNGeneWGS.csv"), usecols=["ModelID", "IsDefaultEntryForModel"])
    wgs = set(wc.loc[wc.IsDefaultEntryForModel == "Yes", "ModelID"])
    up = ps.key.isin(wgs) & (ps.cn_provider != "depmap_wgs")
    log.append(("cn_upgraded_to_depmap_wgs (lines) " + str(ps.loc[up].drop_duplicates("key").cn_provider.value_counts().to_dict()),
                int(ps.loc[up, "key"].nunique())))
    up_keys = set(ps.loc[up, "key"])
    ps.loc[up, "cn_provider"] = "depmap_wgs"
    ps.loc[up, "cn_cvcl"] = pd.NA

    ps.drop(columns="tier").to_csv(os.path.join(P2, "pull_set.v31.tsv"), sep="\t", index=False)
    ps.srx.sort_values().to_csv(os.path.join(P2, "pull_srx.v31.txt"), index=False, header=False)
    inf = sorted(set(ps.loc[ps.cn_provider == "input_inferred", "key"]))
    open(os.path.join(P2, "v31_inferred_keys.txt"), "w").write(",".join(inf) + "\n")

    # all experiments (variants 4-5)
    lab = pd.read_csv(os.path.join(OUT, "atlas_treatment_labels.tsv"), sep="\t", usecols=["srx", "cell_line", "class"])
    extra = [pd.read_csv(p, sep="\t", usecols=["srx", "class"])
             for p in sorted(glob.glob(os.path.join(OUT, "atlas_labels", "labels_1[0-9]*.tsv")))]
    L = pd.concat([lab, *extra]).drop_duplicates("srx", keep="last")
    L = L.merge(lab[["srx", "cell_line"]], on="srx", how="left", suffixes=("", "_l"))
    pert = set(L.loc[L["class"] == "perturbed", "srx"]) | set(F.loc[F.final_class == "perturbed", "srx"])
    deriv = set(F.loc[F.final_note.fillna("").str.startswith("derivative"), "srx"])
    checked = set(pd.read_csv(os.path.join(P2, "pull_set.v2.tsv"), sep="\t").srx)
    for n in ("v3", "v3b", "v3c", "v3d"):
        c = pd.read_csv(os.path.join(OUT, f"{n}_candidate_check.tsv"), sep="\t")
        checked |= set(c.loc[c.status == "pass", "srx"])
    cell_of = dict(zip(lab.srx, lab.cell_line))
    lines = ps.drop_duplicates("cell").set_index("cell")
    rows = []
    for s in sorted(pert - set(ps.srx) - deriv - drop):
        c = cell_of.get(s)
        if c in lines.index and s in checked:
            row = lines.loc[c].to_dict(); row.update(srx=s, cell=c, source="treated"); rows.append(row)
    al = pd.concat([ps.drop(columns="tier"), pd.DataFrame(rows)], ignore_index=True)
    al.to_csv(os.path.join(P2, "pull_set.v31.all.tsv"), sep="\t", index=False)
    redo = al[al.key.isin(up_keys)][["srx", "key", "cn_provider", "cvcl", "cn_cvcl", "cell", "source"]]
    redo.to_csv(os.path.join(P2, "recall_redo.v31.tsv"), sep="\t", index=False)
    log.append(("recall_redo_samples", len(redo)))
    al.srx.sort_values().to_csv(os.path.join(P2, "pull_srx.v31.all.txt"), index=False, header=False)

    for k, v in log:
        print(f"[25] {k}: {v}")
    print(f"[25] baseline v3 {n0} -> v3.1 {len(ps)} experiments on {ps.key.nunique()} keys "
          f"({int(ps.srx.isin(qc).sum())} were QC-pass in v3); input-inferred keys {len(inf)}")
    print(f"[25] all experiments: {len(al)} ({len(rows)} treated added on {pd.DataFrame(rows).key.nunique() if rows else 0} "
          f"lines; derived sublines left out {len(deriv)})")


if __name__ == "__main__":
    main()
