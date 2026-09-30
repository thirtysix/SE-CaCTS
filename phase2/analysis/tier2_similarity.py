#!/usr/bin/env python3
"""Do borderline "tier-2" control arms (dox-induced reporters, dCas9-only, cutting controls, defined media) look like
their line's clean baseline, or like its perturbed experiments?

For every experiment e of a line L, similarity = mean Pearson correlation between e's line-specific profile and
the profiles of L's CLEAN experiments from OTHER studies (a study shares batch, lab and design, so same-study pairs
would flatter anything). Profile = log2(signal + 1) minus the atlas-wide mean at each locus, over the most variable
loci, so the shared super-enhancer landscape does not dominate. Calibration: the same score for L's perturbed
experiments (v2 matrix only; v3 is baseline-only), and whether e's nearest other-study neighbour is its own line.

  python3 phase2/analysis/tier2_similarity.py --matrix phase2/results_v2/atlas.s3.se_signal.tsv.gz --out <prefix>
"""
import argparse
import csv
import glob
import os
from collections import defaultdict

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SE = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = os.path.join(HERE, "out")
TIER2 = {"SRX1048466", "SRX12578555", "SRX12578556", "SRX12877691", "SRX12877692", "SRX4645309", "SRX4645311",
         "SRX4645313", "SRX4645315", "SRX18955857", "SRX18955858", "SRX3590862", "SRX3590863", "SRX11243163",
         "SRX11243164", "SRX3746148", "SRX6853170", "SRX3141136", "SRX3141137", "SRX16888109", "SRX16888110"}
KEEP = {"untreated", "control"}


def labels():
    """srx -> (class, study, cell line). Sonnet's class, with Opus's §36 corrections and every row that raised a
    policy question removed from 'clean' (they are neither clean nor clearly perturbed)."""
    lab = {x["srx"]: x["class"].strip() for p in glob.glob(os.path.join(OUT, "atlas_labels/labels_??.tsv"))
           for x in csv.DictReader(open(p), delimiter="\t")}
    meta = {x["srx"]: (x["group"], x["cell_line"]) for p in glob.glob(os.path.join(OUT, "atlas_labels/input_*.tsv"))
            for x in csv.DictReader(open(p), delimiter="\t")}
    for x in csv.DictReader(open(os.path.join(OUT, "stageB_treatment_manual.tsv")), delimiter="\t"):
        lab.setdefault(x["srx"], x["class"]); meta.setdefault(x["srx"], (x["bioproject"], x["cell_line"]))
    adj = {x["srx"]: x for x in csv.DictReader(open(os.path.join(OUT, "atlas_labels/adjudication/v31_actions.tsv")), delimiter="\t")}
    out = {}
    for s, c in lab.items():
        if s in TIER2:
            k = "tier2"
        elif s in adj and adj[s]["action"] == "remove":
            k = "perturbed"
        elif s in adj and adj[s]["rule_note"].strip():
            k = "policy"
        elif c in KEEP:
            k = "clean"
        elif c == "perturbed":
            k = "perturbed"
        else:
            k = "other"
        out[s] = (k, *meta.get(s, (s, "")))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--top", type=int, default=10000, help="most variable loci used")
    a = ap.parse_args()

    ps = pd.concat([pd.read_csv(p, sep="\t", usecols=["srx", "key", "cell"]) for p in
                    (os.path.join(SE, "phase2/data/pull_set.v2.tsv"), os.path.join(SE, "phase2/data/pull_set.v3.tsv"))])
    line = dict(zip(ps.srx, ps.key)); cellname = dict(zip(ps.srx, ps.cell))
    lab = labels()
    X = pd.read_csv(a.matrix, sep="\t", index_col=0, dtype={0: str}).astype(np.float32)
    X = X[[c for c in X.columns if c in line and c in lab]]
    L = np.log2(X.values + 1.0)
    L -= L.mean(axis=1, keepdims=True)                                   # remove the shared landscape
    var = L.var(axis=1)
    L = L[np.argsort(var)[::-1][:a.top]]
    L -= L.mean(axis=0, keepdims=True)
    L /= np.linalg.norm(L, axis=0, keepdims=True)
    R = L.T @ L                                                          # Pearson between samples
    cols = list(X.columns)
    kind = [lab[c][0] for c in cols]; study = [lab[c][1] for c in cols]; ln = [line[c] for c in cols]
    by_line = defaultdict(list)
    for i, c in enumerate(cols):
        by_line[ln[i]].append(i)
    rows = []
    for i, c in enumerate(cols):
        clean_other = [j for j in by_line[ln[i]] if j != i and kind[j] == "clean" and study[j] != study[i]]
        other_study = [j for j in range(len(cols)) if study[j] != study[i]]
        nn = max(other_study, key=lambda j: R[i, j]) if other_study else None
        rows.append({"srx": c, "line": ln[i], "cell": cellname.get(c, ""), "study": study[i], "kind": kind[i],
                     "n_clean_other_study": len(clean_other),
                     "sim_clean": float(np.mean(R[i, clean_other])) if clean_other else np.nan,
                     "nn_same_line": (ln[nn] == ln[i]) if nn is not None else None})
    D = pd.DataFrame(rows)
    # within-line z against the line's own clean experiments (each clean scored the same way, leave-one-out)
    D["z"] = np.nan
    for l, g in D.groupby("line"):
        ref = g[(g.kind == "clean") & g.sim_clean.notna()].sim_clean
        if len(ref) >= 3 and ref.std() > 0:
            D.loc[g.index, "z"] = (g.sim_clean - ref.mean()) / ref.std()
    D.to_csv(a.out + ".per_experiment.tsv", sep="\t", index=False)
    print(f"{a.matrix}: {len(cols)} experiments, {a.top} loci")
    print(f"{'kind':10s} {'n':>5s} {'scored':>6s} {'sim_clean median':>17s} {'z median':>9s} {'z < -2':>7s} {'nn same line':>13s}")
    for k in ("clean", "tier2", "policy", "perturbed"):
        g = D[D.kind == k]; s = g[g.sim_clean.notna()]
        if not len(g):
            continue
        z = s.z.dropna()
        print(f"{k:10s} {len(g):5d} {len(s):6d} {s.sim_clean.median():17.3f} {z.median() if len(z) else float('nan'):9.2f} "
              f"{(z < -2).sum():4d}/{len(z):<3d} {g.nn_same_line.mean():12.1%}")
    # lines holding both clean and perturbed experiments: paired per-line medians
    pl = []
    for l, g in D.groupby("line"):
        c, p = g[(g.kind == "clean")].sim_clean.dropna(), g[(g.kind == "perturbed")].sim_clean.dropna()
        if len(c) >= 2 and len(p) >= 1:
            pl.append(c.median() - p.median())
    if pl:
        pl = np.array(pl)
        print(f"\nlines with clean and perturbed experiments: {len(pl)}; clean minus perturbed similarity, median "
              f"{np.median(pl):.3f}, positive in {(pl > 0).mean():.0%}")
    t = D[D.kind == "tier2"].sort_values("cell")
    print("\ntier-2 experiments:")
    for r in t.itertuples():
        print(f"  {r.srx:12s} {r.cell[:12]:12s} sim_clean {r.sim_clean:6.3f}  z {r.z:6.2f}  nn_same_line {r.nn_same_line}  "
              f"(clean refs from other studies: {r.n_clean_other_study})")


if __name__ == "__main__":
    main()
