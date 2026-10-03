#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 1 — 28: the v3.2 pull sets (FINDINGS §57): v3.1 plus measured CN for lines scored on input inference, and the
experiments the 2026-09 ChIP-Atlas metadata adds.

From `pull_set.v31.tsv` / `pull_set.v31.all.tsv`:
  1. CN source: lines in the CCMA segment file (phase1/scripts/26) -> `ccma_wgs`; lines accepted by the WGS evaluation
     (`phase2/data/wgs_cn_gate.tsv`: key, accept; phase1/scripts/27) -> `wgs_reads`. Only lines now on
     `input_inferred` move: a measured source never gives way to these. Their samples are listed in
     `recall_redo.v32.tsv` for the fused re-call.
  2. New experiments (`phase1/data/v32_new_h3k27ac.tsv`, QC-pass, in scope) enter only with a baseline label
     (`phase2/analysis/out/v32_relabel/final_labels.tsv`, class untreated / control) and only on lines already in
     v3.1; until those labels exist, none enter and the log says so. They need pulling: `pull_srx.v32.new.txt`.

    python3 phase1/scripts/28_pull_set_v32.py
"""
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
P1, P2 = os.path.join(ROOT, "phase1", "data"), os.path.join(ROOT, "phase2", "data")
OUT = os.path.join(ROOT, "phase2", "analysis", "out")
KEEP = {"untreated", "control"}


def main():
    ps = pd.read_csv(os.path.join(P2, "pull_set.v31.tsv"), sep="\t")
    al = pd.read_csv(os.path.join(P2, "pull_set.v31.all.tsv"), sep="\t")
    log = []

    # 1. CN source
    ccma = set(pd.read_csv(os.path.join(P2, "cn_ccma_wgs.hg38.tsv.gz"), sep="\t", usecols=["key"])["key"])
    gate_f = os.path.join(P2, "wgs_cn_gate.tsv")
    wgs = set()
    if os.path.exists(gate_f):
        g = pd.read_csv(gate_f, sep="\t")
        wgs = set(g.loc[g["accept"].astype(str).str.lower().isin(["yes", "true", "1"]), "key"])
    else:
        log.append(("no wgs_cn_gate.tsv: no wgs_reads lines", 0))
    moved = {}
    for df in (ps, al):
        inf = df.cn_provider.eq("input_inferred")
        for keys, prov in ((ccma, "ccma_wgs"), (wgs - ccma, "wgs_reads")):
            m = inf & df.key.isin(keys)
            df.loc[m, "cn_provider"] = prov
            df.loc[m, "cn_cvcl"] = pd.NA
            moved.update({k: prov for k in df.loc[m, "key"]})
    log.append((f"cn_moved_from_input_inferred (lines) {pd.Series(moved).value_counts().to_dict()}", len(moved)))

    # 2. new experiments from the 2026-09 metadata
    new = pd.read_csv(os.path.join(P1, "v32_new_h3k27ac.tsv"), sep="\t")
    lab_f = os.path.join(OUT, "v32_relabel", "final_labels.tsv")
    add = []
    if os.path.exists(lab_f):
        F = pd.read_csv(lab_f, sep="\t", dtype=str)
        ok = set(F.loc[F.final_class.isin(KEEP), "srx"])
        lines = ps.drop_duplicates("cvcl").set_index("cvcl")
        for r in new.itertuples():
            if r.srx in ok and r.cvcl in lines.index and r.srx not in set(ps.srx):
                row = lines.loc[r.cvcl].to_dict()
                row.update(srx=r.srx, source="v32_metadata")
                add.append(row)
        log.append(("new_experiments_labelled_baseline", len(ok & set(new.srx))))
    else:
        log.append(("no v32 labels yet: new experiments held back", len(new)))
    log.append(("new_experiments_added (on v3.1 lines)", len(add)))
    if add:
        ps = pd.concat([ps, pd.DataFrame(add)], ignore_index=True)
        al = pd.concat([al, pd.DataFrame(add)], ignore_index=True)

    ps.to_csv(os.path.join(P2, "pull_set.v32.tsv"), sep="\t", index=False)
    ps.srx.sort_values().to_csv(os.path.join(P2, "pull_srx.v32.txt"), index=False, header=False)
    al.to_csv(os.path.join(P2, "pull_set.v32.all.tsv"), sep="\t", index=False)
    al.srx.sort_values().to_csv(os.path.join(P2, "pull_srx.v32.all.txt"), index=False, header=False)
    pd.Series(sorted(a["srx"] for a in add), dtype=str).to_csv(os.path.join(P2, "pull_srx.v32.new.txt"),
                                                               index=False, header=False)
    inf = sorted(set(ps.loc[ps.cn_provider == "input_inferred", "key"]))
    open(os.path.join(P2, "v32_inferred_keys.txt"), "w").write(",".join(inf) + "\n")
    redo = al[al.key.isin(moved)][["srx", "key", "cn_provider", "cvcl", "cn_cvcl", "cell", "source"]]
    redo.to_csv(os.path.join(P2, "recall_redo.v32.tsv"), sep="\t", index=False)
    log.append(("recall_redo_samples", len(redo)))
    for k, v in log:
        print(f"[28] {k}: {v}")
    print(f"[28] v3.2 baseline {len(ps)} experiments on {ps.key.nunique()} keys; all {len(al)}; "
          f"input-inferred keys {len(inf)} (v3.1: 72)")


if __name__ == "__main__":
    main()
