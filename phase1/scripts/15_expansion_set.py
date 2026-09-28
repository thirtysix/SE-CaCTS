#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Phase 1 — 15: the CN-expansion pull set (ROADMAP Phase 3b).

Selects the human cancer lines that are NOT yet in the atlas but now carry measured copy number from a
second source, in the provider priority order the scorer uses (DepMap WGS -> CMP WES pureCN 2025 ->
DepMap MC_WES; CN_COVERAGE.md §5/§5b), and every QC-passing H3K27ac experiment on them.

A line counts only if a CN TRACK is actually built for it — a cross-reference is not the data
(FINDINGS 2026-08-18 §4).

Outputs:
  phase1/data/expansion_lines.tsv   one row per new line: identity, CN provider, Oncotree labels, n exp
  phase2/data/pull_srx.expand.txt   the new SRX only (the array manifest for this pull)
  phase2/data/pull_set.v2.tsv       old pull set (minus the retired SRX) + new, same columns plus
                                    `key` (ModelID, or CVCL when the line has no DepMap model) and
                                    `cn_provider`
  phase2/data/pull_srx.v2.txt       the full v2 manifest the reduce counts against

Reproduce:  python3 phase1/scripts/15_expansion_set.py
"""
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "cnrose"))
from secacts_env import DATAROOT, cache_path                     # noqa: E402
from cnrose.cn.cmp import CellModelPassportsWesCN                 # noqa: E402
from cnrose.cn.depmap import DepMapMcWesCN, load_gene_coords      # noqa: E402

P1 = os.path.join(ROOT, "phase1", "data")
P2 = os.path.join(ROOT, "phase2", "data")
RETIRED = {"SRX20868733"}                     # permanent 404 (phase2/results/failed_srx.txt)

DEPMAP = os.path.join(DATAROOT, "DepMap", "2026q1")
CMP = os.path.join(DATAROOT, "CellModelPassports")
GTF = os.path.join(DATAROOT, "0.human_genome", "Homo_sapiens.GRCh38.106.chr.gtf.gz")


def main():
    status = pd.read_csv(os.path.join(P1, "cell_line_cn_status_2026-08-15.tsv"), sep="\t")
    lin = pd.read_csv(os.path.join(P1, "lineage_resolved.tsv"), sep="\t").set_index("cvcl")
    cnsrc = pd.read_csv(os.path.join(P1, "cn_source.tsv"), sep="\t").set_index("cvcl")
    man = pd.read_csv(os.path.join(P1, "phase1_manifest.tsv"), sep="\t")

    cand = status[status["data_status"] == "needs pull"].copy()
    cand["model_id"] = cand["cvcl"].map(cnsrc["model_id"])
    print(f"[15] candidate lines not in the atlas: {len(cand)}", file=sys.stderr)

    # provider 2: CMP WES pureCN 2025 — a track must actually build
    cmp = CellModelPassportsWesCN(os.path.join(CMP, "WES_pureCN_CNV_genes_latest.csv.gz"),
                                  os.path.join(CMP, "model_list_20240110.csv"), cache_dir=cache_path("cmp_wes"))
    cmp.preload(cand["cvcl"])
    has_cmp = {c for c in cand["cvcl"] if cmp.track(c) is not None}

    # provider 3: DepMap MC_WES, only for lines CMP did not reach (order of acquisition matters)
    rest = cand[~cand["cvcl"].isin(has_cmp) & cand["model_id"].notna()]
    gc = load_gene_coords(GTF, cache_path=cache_path("gene_coords.GRCh38.106.tsv"))
    mcw = DepMapMcWesCN(os.path.join(DEPMAP, "OmicsCNGeneMC_WES.csv"),
                        os.path.join(DEPMAP, "ModelCondition.csv"), gc)
    mcw.preload(list(rest["model_id"]))
    has_mcw = {c for c, m in zip(rest["cvcl"], rest["model_id"]) if mcw.track(m) is not None}

    cand["cn_provider"] = cand["cvcl"].map(
        lambda c: "cmp_wes" if c in has_cmp else ("depmap_mc_wes" if c in has_mcw else ""))
    new = cand[cand["cn_provider"] != ""].copy()
    new["key"] = new["model_id"].where(new["model_id"].notna(), new["cvcl"])
    for col in ("lineage", "primary_disease", "subtype", "lineage_source"):
        new[col] = new["cvcl"].map(lin[col])

    exp = man[(man["qc_pass"] == 1) & man["cvcl"].isin(new["cvcl"])].copy()
    new["n_qc_exp"] = new["cvcl"].map(exp.groupby("cvcl").size()).fillna(0).astype(int)
    new = new[new["n_qc_exp"] > 0]
    cols = ["cvcl", "cell_line", "key", "model_id", "sidm", "cn_provider", "lineage", "primary_disease",
            "subtype", "lineage_source", "n_qc_exp"]
    new[cols].sort_values(["lineage", "cell_line"]).to_csv(
        os.path.join(P1, "expansion_lines.tsv"), sep="\t", index=False)

    # pull sets: old (minus retired) + new, with the scoring key and CN provider carried per SRX
    old = pd.read_csv(os.path.join(P2, "pull_set.tsv"), sep="\t")
    old = old[~old["srx"].isin(RETIRED)].copy()
    old["key"], old["cn_provider"] = old["model_id"], "depmap_wgs"
    exp = exp.merge(new[["cvcl", "key", "cn_provider"]], on="cvcl")
    exp = exp[~exp["srx"].isin(old["srx"]) & ~exp["srx"].isin(RETIRED)]
    exp = exp.rename(columns={"cell": "cell"})[["srx", "cell", "cvcl", "model_id", "lineage", "subtype",
                                                "key", "cn_provider"]]
    exp["lineage"] = exp["cvcl"].map(lin["lineage"])
    exp["subtype"] = exp["cvcl"].map(lin["subtype"])
    v2 = pd.concat([old, exp], ignore_index=True)
    v2.to_csv(os.path.join(P2, "pull_set.v2.tsv"), sep="\t", index=False)
    exp["srx"].sort_values().to_csv(os.path.join(P2, "pull_srx.expand.txt"), index=False, header=False)
    v2["srx"].sort_values().to_csv(os.path.join(P2, "pull_srx.v2.txt"), index=False, header=False)

    print(f"[15] new CN-correctable lines with >=1 QC-pass experiment: {len(new)} "
          f"(cmp_wes {int((new.cn_provider == 'cmp_wes').sum())}, "
          f"depmap_mc_wes {int((new.cn_provider == 'depmap_mc_wes').sum())}; "
          f"{int(new['model_id'].isna().sum())} without a DepMap ModelID)", file=sys.stderr)
    print(f"[15] new experiments: {len(exp)}  |  v2 manifest: {len(v2)} "
          f"(old {len(old)} + new {len(exp)})", file=sys.stderr)
    print(new.groupby("lineage").agg(lines=("cvcl", "size"), exps=("n_qc_exp", "sum"))
          .sort_values("lines", ascending=False).to_string(), file=sys.stderr)


if __name__ == "__main__":
    main()
