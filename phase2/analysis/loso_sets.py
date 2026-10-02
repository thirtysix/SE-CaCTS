#!/usr/bin/env python3
"""Leave-one-study-out exclusion sets for the lineage rescore (FINDINGS §42).

Per lineage with >= 2 studies among its QC-pass experiments:
  top  = every experiment of that lineage from its largest study (most lines covered; ties -> most experiments)
  randN = matched controls: other studies of the same lineage drawn at random until they cover as many lines
         as the top study (the draw that lands closest wins), so the loss of power is comparable
Only the lineage's own experiments are removed; a study's lines in other lineages stay. Each set is scored at
the lineage level with the v3 settings (score_arm.slurm ARM=loso), and retention of the lineage's v3 calls is
compared between the two (loso_eval.py).

    python3 phase2/analysis/loso_sets.py --out phase2/analysis/out/loso
"""
import argparse, gzip, importlib.util, os, re
import numpy as np
import pandas as pd

SECACTS = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
spec = importlib.util.spec_from_file_location("cf", os.path.join(SECACTS, "phase2/scripts/74_called_filter.py"))
cf = importlib.util.module_from_spec(spec); spec.loader.exec_module(cf)
from secacts_env import DATAROOT                                         # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--signal", default=os.path.join(SECACTS, "phase2/results_v3/atlas.s3.se_signal.tsv.gz"))
    ap.add_argument("--pull-set", default=os.path.join(SECACTS, "phase2/data/pull_set.v3.tsv"))
    ap.add_argument("--srx-study", default=os.path.join(SECACTS, "phase2/data/srx_study.tsv"))
    ap.add_argument("--model", default=os.path.join(DATAROOT, "DepMap/2026q1/Model.csv"))
    ap.add_argument("--lines-meta", default=os.path.join(SECACTS, "phase1/data/lineage_resolved.tsv"))
    ap.add_argument("--seed", type=int, default=20261002)
    ap.add_argument("--n-rand", type=int, default=2, help="matched random controls per lineage")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    with gzip.open(a.signal, "rt") as fh:
        qc = fh.readline().rstrip("\n").split("\t")[1:]
    ps = pd.read_csv(a.pull_set, sep="\t")
    ps = ps[ps["srx"].isin(qc)].copy()
    ps.loc[ps["srx"].eq("SRX10809652"), "key"] = "ACH-000768"           # as score_pilot.py EXTRA_MODEL
    L = cf.labels(ps, a.model, a.lines_meta)
    ps["lin"] = ps["key"].map(L["OncotreeLineage"])
    st = pd.read_csv(a.srx_study, sep="\t").set_index("srx")["study"]
    ps["study"] = ps["srx"].map(st)
    ps["study"] = ps["study"].fillna(ps["srx"])                         # unresolved = own study (as §42)
    ps = ps.dropna(subset=["lin"])

    rng = np.random.default_rng(a.seed)
    rows = []
    for lin, d in ps.groupby("lin"):
        cover = d.groupby("study").agg(lines=("key", "nunique"), exps=("srx", "size"))
        if len(cover) < 2:
            continue
        cover = cover.sort_values(["lines", "exps"], ascending=False)
        top = cover.index[0]
        target = int(cover.iloc[0]["lines"])
        others = list(cover.index[1:])
        arms = [("top", [top])]
        for r in range(1, a.n_rand + 1):
            picked, best, best_gap = [], None, None
            for s in rng.permutation(others):
                picked.append(s)
                n = d.loc[d["study"].isin(picked), "key"].nunique()
                gap = abs(n - target)
                if best_gap is None or gap < best_gap:
                    best, best_gap = list(picked), gap
                if n >= target:
                    break
            arms.append((f"rand{r}", best))
        slug = re.sub(r"[^A-Za-z0-9]+", "_", lin).strip("_")
        for arm, studies in arms:
            x = d[d["study"].isin(studies)]
            left = d[~d["study"].isin(studies)]
            x["srx"].to_csv(os.path.join(a.out, f"{slug}.{arm}.srx"), index=False, header=False)
            rows.append({"lineage": lin, "slug": slug, "arm": arm, "studies": ",".join(studies),
                         "n_studies": len(studies), "exps_removed": len(x), "lines_touched": x["key"].nunique(),
                         "lines_lost": len(set(d["key"]) - set(left["key"])),
                         "lines_total": d["key"].nunique(), "exps_total": len(d), "studies_total": len(cover)})
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(a.out, "loso_sets.tsv"), sep="\t", index=False)
    print(out.to_string(index=False))


if __name__ == "__main__":
    main()
