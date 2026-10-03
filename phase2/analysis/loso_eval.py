#!/usr/bin/env python3
"""Retention of each lineage's v3.0.1 calls when its largest study is left out, against matched random removals.

Each run (phase2/analysis/loso_sets.py -> score_arm.slurm ARM=loso) is filtered with the v3.0.1 rule using only
the experiments that remain, so an SE called only by the removed study also stops counting as an SE of the
lineage. Retention = share of the lineage's v3.0.1 calls (same locus) still called; the same for the top 100.
Other lineages are reported too: their calls should barely move, since only the target lineage lost experiments.

    python3 phase2/analysis/loso_eval.py --runs phase2/analysis/out/loso/runs --sets phase2/analysis/out/loso

v3.1 (loso_submit.sh V31F=1) applies the rule inside the permutation null, with the remaining experiments, so the runs
are read as written:
    python3 phase2/analysis/loso_eval.py --rule-in-scores --runs phase2/analysis/out/loso_v31f/runs \
        --sets phase2/analysis/out/loso_v31f --base phase2/scores_v31f/atlas.s3.perm.OncotreeLineage.specific.tsv.gz \
        --out phase2/analysis/out/loso_v31f/loso_retention.tsv
"""
import argparse, importlib.util, os, sys
import numpy as np
import pandas as pd
import scipy.sparse as sp

SECACTS = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
spec = importlib.util.spec_from_file_location("cf", os.path.join(SECACTS, "phase2/scripts/74_called_filter.py"))
cf = importlib.util.module_from_spec(spec); spec.loader.exec_module(cf)
from secacts_env import DATAROOT                                         # noqa: E402
LIN = "OncotreeLineage"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", required=True, help="dir of <slug>.<arm>/ run outputs")
    ap.add_argument("--sets", required=True, help="dir with loso_sets.tsv and the .srx lists")
    ap.add_argument("--base", default=os.path.join(SECACTS, "phase2/scores_v3c/atlas.s3.perm.OncotreeLineage.specific.tsv.gz"))
    ap.add_argument("--results", default=os.path.join(SECACTS, "phase2/results_v3"))
    ap.add_argument("--pull-set", default=os.path.join(SECACTS, "phase2/data/pull_set.v3.tsv"))
    ap.add_argument("--out", default=os.path.join(SECACTS, "phase2/analysis/out/loso/loso_retention.tsv"))
    ap.add_argument("--rule-in-scores", action="store_true",
                    help="the runs already apply the SE-of-its-group rule (v3.1): no post-hoc filter")
    a = ap.parse_args()
    if a.rule_in_scores:
        return evaluate(a, lambda r, f: pd.read_csv(f, sep="\t"))

    ps = pd.read_csv(a.pull_set, sep="\t")
    lab = cf.labels(ps, os.path.join(DATAROOT, "DepMap/2026q1/Model.csv"),
                    os.path.join(SECACTS, "phase1/data/lineage_resolved.tsv"))
    cat = pd.read_csv(os.path.join(a.results, "atlas.s3.union_catalog.bed.gz"), sep="\t", header=None,
                      usecols=[0, 1, 2, 3], names=["chrom", "start", "end", "se"])
    se_idx = {s: i for i, s in enumerate(cat.se)}
    pres = pd.read_csv(os.path.join(a.results, "atlas.s3.se_presence.tsv.gz"), sep="\t", index_col=0)
    pres = pres.reindex(cat.se).fillna(0)
    P = sp.csr_matrix((pres.values > 0).astype(np.float32))
    A = cf.overlap_graph(cat)
    srx_key = dict(zip(ps["srx"], ps["key"])); srx_key.update(cf.EXTRA_MODEL)
    srx_lin = [lab.at[srx_key[s], LIN] if srx_key.get(s) in lab.index else None for s in pres.columns]

    def filtered(r, f):
        drop = set(open(os.path.join(a.sets, f"{r.slug}.{r.arm}.srx")).read().split())
        g = [None if s in drop else x for s, x in zip(pres.columns, srx_lin)]
        groups = sorted({x for x in g if isinstance(x, str)})
        cov = cf.coverage(P, A, g, groups)
        return cf.filter_calls(pd.read_csv(f, sep="\t"), cov, se_idx, {x: j for j, x in enumerate(groups)})[0]
    evaluate(a, filtered)


def evaluate(a, read_run):
    base = pd.read_csv(a.base, sep="\t")
    sets = pd.read_csv(os.path.join(a.sets, "loso_sets.tsv"), sep="\t")
    rows = []
    for r in sets.itertuples():
        f = os.path.join(a.runs, f"{r.slug}.{r.arm}", f"atlas.s3.perm.{LIN}.specific.tsv.gz")
        if not os.path.exists(f):
            print(f"[eval] missing {r.slug}.{r.arm}", file=sys.stderr)
            continue
        run = read_run(r, f)
        for lin, b in base.groupby("group"):
            l = run[run.group == lin]
            bs, ls = set(b.se), set(l.se)
            top = set(b.nsmallest(100, "rank").se)
            rows.append(dict(lineage=r.lineage, arm=r.arm, group=lin, target=lin == r.lineage,
                             lines_lost=r.lines_lost, exps_removed=r.exps_removed,
                             base_calls=len(bs), run_calls=len(ls), retained=len(bs & ls),
                             retention=round(len(bs & ls) / len(bs), 3) if bs else np.nan,
                             top100_retention=round(len(top & ls) / len(top), 3) if top else np.nan,
                             new=len(ls - bs)))
    R = pd.DataFrame(rows)
    R.to_csv(a.out, sep="\t", index=False)
    T = R[R.target].pivot_table(index="lineage", columns="arm",
                                values=["base_calls", "run_calls", "retention", "top100_retention", "lines_lost"],
                                aggfunc="first")
    with pd.option_context("display.width", 250, "display.max_rows", 200):
        print(T.to_string())
        o = R[~R.target]
        print(f"\nnon-target lineages: median retention {o.retention.median():.3f} "
              f"(p05 {o.retention.quantile(0.05):.3f}) over {len(o)} lineage x run pairs")


if __name__ == "__main__":
    main()
