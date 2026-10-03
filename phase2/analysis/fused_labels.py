#!/usr/bin/env python3
"""Labels for the calls of the fused build (FINDINGS §54-55; user decision 2026-10-03).

Two specificity statistics, never merged into one FDR: each comparison is scored on CN-corrected and on uncorrected
signal over the SAME tests (FUSED=1: eligibility from the fused calls in both arms). Each (group, locus) call gets
  cn_status   robust (called by both statistics), unmasked (corrected only), amplicon (uncorrected only).
              The headline set is robust + unmasked, i.e. the corrected arm's calls; amplicon-driven calls are listed
              alongside, with their own FDR.
and, from the per-experiment fused labels (aggregate_v31.py `se_label.fu.tsv`), how the group's own experiments call it:
  se_label    the most CN-robust label among the group's experiments that call an SE there: core / unmasked /
              agnostic / gain (region CN < 2) / amplified (2-3) / high (>= 3)
  n_exp, n_exp_amp   the group's experiments calling it, and how many of those only through amplification (gain+)
  cn_sources  CN sources of the lines whose experiments call it (input_inferred is the weakest evidence)

    ~/miniconda3/envs/atac_hdac/bin/python phase2/analysis/fused_labels.py --scores phase2/scores_v31f \\
        --results phase2/results_v31f --pull-set phase2/data/pull_set.v31.tsv
writes <scores>/fused_labels.<comparison>.tsv.gz and prints counts per comparison.
"""
from __future__ import annotations

import argparse, os, sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SECACTS = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, SECACTS)
sys.path.insert(0, os.path.join(SECACTS, "phase2"))
from secacts_env import DATAROOT                                    # noqa: E402
from aggregate_v31 import LABELS                                    # noqa: E402

CODE = {v: k for k, v in LABELS.items()}
GROUP_LEVELS = ["OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype"]
COMPARISONS = [("groups", "atlas.s3.perm", GROUP_LEVELS)] + \
              [(f"lines.{c}", f"atlas.s3.lines.{c}", ["line"]) for c in ("all", "sub", "dis", "lin")]


def read_dump(path):
    if not os.path.exists(path):
        return None
    return pd.read_csv(path, sep="\t", usecols=["group", "se", "rank", "fdr", "cn_mean"])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scores", required=True, help="scores dir of the fused build (corrected + .nocn arms)")
    ap.add_argument("--results", required=True, help="results dir with atlas.s3.se_label.fu.tsv(.gz)")
    ap.add_argument("--pull-set", required=True)
    ap.add_argument("--model", default=os.path.join(DATAROOT, "DepMap/2026q1/Model.csv"))
    ap.add_argument("--lines-meta", default=os.path.join(SECACTS, "phase1/data/lineage_resolved.tsv"))
    ap.add_argument("--only", help="comma list of comparisons (groups, lines.all, lines.sub, lines.dis, lines.lin)")
    a = ap.parse_args()

    lab_path = os.path.join(a.results, "atlas.s3.se_label.fu.tsv")
    lab_path = lab_path if os.path.exists(lab_path) else lab_path + ".gz"
    hdr = pd.read_csv(lab_path, sep="\t", nrows=0).columns
    L = pd.read_csv(lab_path, sep="\t", index_col=0, dtype={c: np.uint8 for c in hdr[1:]})
    ps = pd.read_csv(a.pull_set, sep="\t").drop_duplicates("srx").set_index("srx")
    samples = [s for s in L.columns if s in ps.index]
    L = L[samples]
    key_of = ps.loc[samples, "key"]
    src_of_key = ps.dropna(subset=["key"]).drop_duplicates("key").set_index("key")["cn_provider"]

    model = pd.read_csv(a.model, index_col="ModelID")
    extra = sorted(set(key_of.dropna()) - set(model.index))
    if extra:                                                       # as score_pilot.py: Cellosaurus-NCIt crosswalk
        lm = pd.read_csv(a.lines_meta, sep="\t").set_index("cvcl").loc[extra]
        model = pd.concat([model, pd.DataFrame(
            {"StrippedCellLineName": lm["cell_line"], "OncotreeLineage": lm["lineage"],
             "OncotreePrimaryDisease": lm["primary_disease"], "OncotreeSubtype": lm["subtype"]},
            index=pd.Index(extra, name="ModelID"))])
    meta = model.reindex(key_of.values)
    meta.index = key_of.index                                       # srx -> labels
    meta["line"] = meta["StrippedCellLineName"]
    meta["src"] = key_of.map(src_of_key).values
    print(f"[fused] {L.shape[0]} loci x {len(samples)} samples, {key_of.nunique()} lines; label codes "
          f"{ {CODE[c]: int((L.values == c).sum()) for c in CODE} }", file=sys.stderr, flush=True)

    Lv = L.values
    loc_of = {se: i for i, se in enumerate(L.index)}
    want = set(a.only.split(",")) if a.only else None
    for comp, prefix, levels in COMPARISONS:
        if want and comp not in want:
            continue
        out = []
        for lev in levels:
            C = read_dump(os.path.join(a.scores, f"{prefix}.{lev}.specific.tsv.gz"))
            U = read_dump(os.path.join(a.scores, f"{prefix}.nocn.{lev}.specific.tsv.gz"))
            if C is None or U is None:
                print(f"[fused] {comp} {lev}: missing {'corrected' if C is None else 'uncorrected'} dump, skipped",
                      file=sys.stderr)
                continue
            M = C.merge(U, on=["group", "se"], how="outer", suffixes=("_c", "_u"))
            M["cn_status"] = np.where(M["fdr_c"].notna() & M["fdr_u"].notna(), "robust",
                                      np.where(M["fdr_c"].notna(), "unmasked", "amplicon"))
            M["cn_mean"] = M["cn_mean_c"].fillna(M["cn_mean_u"])
            M.insert(0, "level", lev)
            col = meta[lev]
            se_lab, n_exp, n_amp, srcs = [], [], [], []
            for g, idx in M.groupby("group").groups.items():
                cols = np.flatnonzero((col == g).values)
                rows = np.array([loc_of.get(s, -1) for s in M.loc[idx, "se"]])
                if len(cols) == 0 or (rows < 0).any():
                    raise SystemExit(f"[fused] {comp} {lev} group {g!r}: {len(cols)} samples, "
                                     f"{int((rows < 0).sum())} loci not in the label matrix")
                sub = Lv[np.ix_(rows, cols)]
                called = sub > 0
                lo = np.where(called, sub, 255).min(axis=1)
                sv = meta["src"].values[cols]
                src = [",".join(sorted(set(sv[called[i]]))) for i in range(len(rows))]
                se_lab += list(zip(idx, [CODE.get(int(x), "none") for x in lo]))
                n_exp += list(zip(idx, called.sum(axis=1)))
                n_amp += list(zip(idx, (sub >= LABELS["gain"]).sum(axis=1)))
                srcs += list(zip(idx, src))
            for name, pairs in (("se_label", se_lab), ("n_exp", n_exp), ("n_exp_amp", n_amp), ("cn_sources", srcs)):
                M[name] = pd.Series(dict(pairs))
            out.append(M[["level", "group", "se", "cn_status", "rank_c", "fdr_c", "rank_u", "fdr_u", "cn_mean",
                          "se_label", "n_exp", "n_exp_amp", "cn_sources"]])
        if not out:
            continue
        R = pd.concat(out, ignore_index=True)
        path = os.path.join(a.scores, f"fused_labels.{comp}.tsv.gz")
        R.to_csv(path, sep="\t", index=False, float_format="%.6g")
        print(f"\n== {comp} -> {path}")
        for lev, d in R.groupby("level", sort=False):
            st = d["cn_status"].value_counts()
            print(f"  {lev:<24} corrected {int(st.get('robust', 0) + st.get('unmasked', 0)):>7,} "
                  f"(robust {int(st.get('robust', 0)):,}, unmasked {int(st.get('unmasked', 0)):,}); "
                  f"amplicon-driven {int(st.get('amplicon', 0)):,}; groups {d['group'].nunique()}")
            ct = pd.crosstab(d["cn_status"], d["se_label"]).reindex(columns=[c for c in LABELS if c in
                                                                             set(d["se_label"])])
            print("    " + ct.to_string().replace("\n", "\n    "))
            amp = d[d["cn_status"] == "amplicon"]["cn_mean"]
            if len(amp):
                print(f"    amplicon-driven group CN: median {amp.median():.2f}, >= 2 {100 * (amp >= 2).mean():.0f}%")


if __name__ == "__main__":
    main()
