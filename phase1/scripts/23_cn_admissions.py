#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 1 — 23: admit the lines held back for lack of CN, each on its best CN source (FINDINGS §21-24).

Priority per line (measured before inferred):
  relative         a Cellosaurus parent/child/same-individual line with DepMap WGS / MC_WES / CMP WES CN.
                   If that relative is already scored, the line JOINS it (same individual, same CN; the
                   signal check must confirm the relative is the best-matching line); otherwise the line is
                   scored under the relative's DepMap ModelID.
  ccle_snp6        the line's own CCLE 2019 SNP6 segments (validated r 0.870 vs WGS, §23)
  input_inferred   CN from its most spatially coherent ChIP input, admitted if lag-1 autocorrelation >= 0.5
                   (§24), applied with the global slope 0.81
Plus the phase-1 homonym fix: ChIP-Atlas "C4-2" is LNCaP C4-2 (CVCL_4782), not CVCL_WI17 -> relative of
LNCaP clone FGC.

Outputs (phase2/analysis/out/): cn_admissions.tsv (one row per line), v3d_candidates.tsv (per experiment, for
v3_candidate_check.py); phase2/data/pull_srx.v3d.txt (QC-pass experiments not yet pulled).
"""
import os
import re
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
from secacts_env import DATAROOT                                  # noqa: E402

P1 = os.path.join(ROOT, "phase1", "data")
P2 = os.path.join(ROOT, "phase2", "data")
OUT = os.path.join(ROOT, "phase2", "analysis", "out")
HOMONYM_FIX = {"CVCL_WI17": ("CVCL_4782", "LNCaP C4-2", "CVCL_1379")}   # resolved -> (true, name, CN relative)
LAG1_MIN = 0.5


def main():
    reg = pd.read_csv(os.path.join(P1, "h3k27ac_registry.tsv"), sep="\t", low_memory=False)
    held = reg[(reg.decision == "not pulled") & (reg.qc_pass_meta == 1)]
    routes = pd.read_csv(os.path.join(OUT, "cn_routes.tsv"), sep="\t")
    gate = pd.read_csv(os.path.join(OUT, "cn_input_gate.inputs.tsv"), sep="\t")
    gate = gate[gate.role == "target"].sort_values("lag1_autocorr", ascending=False).drop_duplicates("key")
    ps = pd.read_csv(os.path.join(P2, "pull_set.v2.tsv"), sep="\t")
    atlas_keys = set(ps.key)
    lin = pd.read_csv(os.path.join(P1, "lineage_resolved.tsv"), sep="\t").drop_duplicates("cvcl").set_index("cvcl")
    m = pd.read_csv(os.path.join(DATAROOT, "DepMap", "2026q1", "Model.csv"), usecols=["ModelID", "RRID"])
    model_of = m.dropna().drop_duplicates("RRID").set_index("RRID").ModelID

    rows = []
    for r in routes.itertuples():
        cv = r.cvcl
        if cv in HOMONYM_FIX:
            true, name, relcv = HOMONYM_FIX[cv]
            mid = model_of.get(relcv)
            rows.append(dict(cvcl=cv, true_cvcl=true, cell_line=name, lineage="Prostate", cn_provider="depmap_wgs",
                             key=mid, relation=f"homonym fix; child of LNCaP clone FGC ({mid})",
                             joins=mid in atlas_keys))
            continue
        routes_s = str(r.routes) if isinstance(r.routes, str) else ""
        first = next((x.strip() for x in routes_s.split(" | ") if x.strip().startswith("relative:")), "")
        rel = re.match(r"relative: (parent|child|same individual) (.+?) (CVCL_\w+): (.+)$", first)
        if r.best_route == "relative" and rel:
            relcv, srcs = rel.group(3), rel.group(4)
            mid = model_of.get(relcv)
            prov = "depmap_wgs" if "DepMap WGS" in srcs else "depmap_mc_wes" if "DepMap MC_WES" in srcs \
                else "cmp_wes" if "CMP WES" in srcs else "ccle_snp6"
            key = mid if mid else relcv
            rows.append(dict(cvcl=cv, cell_line=r.cell_line, lineage=r.lineage, cn_provider=prov, key=key,
                             cn_cvcl=relcv, relation=f"{rel.group(1)} {rel.group(2)}", joins=key in atlas_keys))
            continue
        if r.best_route == "ccle2019":
            mid = model_of.get(cv)
            rows.append(dict(cvcl=cv, cell_line=r.cell_line, lineage=r.lineage, cn_provider="ccle_snp6",
                             key=mid or cv, relation="own CCLE 2019 SNP6", joins=False))
            continue
        gl = gate[gate.key == cv]
        if len(gl) and gl.lag1_autocorr.iloc[0] >= LAG1_MIN:
            rows.append(dict(cvcl=cv, cell_line=r.cell_line, lineage=r.lineage, cn_provider="input_inferred",
                             key=model_of.get(cv) or cv, cn_input=gl.input_srx.iloc[0],
                             relation=f"input {gl.input_srx.iloc[0]}, lag-1 {gl.lag1_autocorr.iloc[0]:.2f}",
                             joins=False))
    A = pd.DataFrame(rows)
    A["n_exp"] = A.cvcl.map(held.groupby("cvcl").srx.nunique()).fillna(0).astype(int)
    A.to_csv(os.path.join(OUT, "cn_admissions.tsv"), sep="\t", index=False)
    print(A.groupby("cn_provider").agg(lines=("cvcl", "size"), exps=("n_exp", "sum"),
                                        joins=("joins", "sum")).to_string(), file=sys.stderr)

    ex = held.merge(A, on="cvcl", suffixes=("", "_a"))
    c = pd.DataFrame({"srx": ex.srx, "cvcl": ex.cvcl, "cell_line": ex.cell_line_a, "lineage": ex.lineage_a,
                      "route": "cn-" + ex.cn_provider,
                      "line_status": ex.joins.map({True: "scored_v2", False: "new"}),
                      "key": ex.key.where(ex.joins, "")}).drop_duplicates("srx")
    c.to_csv(os.path.join(OUT, "v3d_candidates.tsv"), sep="\t", index=False)
    open(os.path.join(P2, "pull_srx.v3d.txt"), "w").write("\n".join(sorted(c.srx)) + "\n")
    print(f"[23] {len(A)} lines admitted, {c.srx.nunique()} QC-pass experiments to pull "
          f"({int(c.line_status.eq('scored_v2').sum())} join a scored relative)", file=sys.stderr)


if __name__ == "__main__":
    main()
