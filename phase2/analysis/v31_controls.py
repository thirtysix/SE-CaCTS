#!/usr/bin/env python3
"""v3.1 item 7: a matched control for the 80 input-inferred-CN lines (FINDINGS §37).

Draws as many measured-CN lines as there are input-inferred lines, lineage by lineage (a lineage short of
measured lines is topped up from the rest), so `score_arm.slurm ARM=randexcl` removes a panel-matched set of
lines with GOOD copy number. If dropping the inferred lines changes calls more than dropping these does, the
change is about their copy number; if not, it is about losing that many lines.

    python3 phase2/analysis/v31_controls.py --pull-set phase2/data/pull_set.v31.tsv --out phase2/data/v31_randexcl_keys.txt
"""
import argparse
import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pull-set", required=True)
    ap.add_argument("--qc", default="phase2/data/qc_srx.v3.txt")
    ap.add_argument("--seed", type=int, default=20261002)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    ps = pd.read_csv(a.pull_set, sep="\t")
    qc = set(open(a.qc).read().split())
    lines = ps[ps.srx.isin(qc)].drop_duplicates("key")[["key", "lineage", "cn_provider"]]
    inf = lines[lines.cn_provider == "input_inferred"]
    meas = lines[lines.cn_provider != "input_inferred"]
    rng = np.random.default_rng(a.seed)
    pick = []
    for lin, n in inf.lineage.value_counts().items():
        pool = meas[(meas.lineage == lin) & ~meas.key.isin(pick)].key.tolist()
        pick += list(rng.choice(pool, size=min(n, len(pool)), replace=False))
    short = len(inf) - len(pick)
    if short > 0:
        pool = meas[~meas.key.isin(pick)].key.tolist()
        pick += list(rng.choice(pool, size=short, replace=False))
    open(a.out, "w").write(",".join(sorted(pick)) + "\n")
    got = meas.set_index("key").loc[pick].lineage.value_counts()
    print(f"{len(pick)} measured-CN lines for {len(inf)} input-inferred ({short} topped up across lineages)")
    print(pd.DataFrame({"inferred": inf.lineage.value_counts(), "control": got}).fillna(0).astype(int).to_string())


if __name__ == "__main__":
    main()
