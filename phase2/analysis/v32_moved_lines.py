#!/usr/bin/env python3
"""v3.2: what measured CN changed for the lines that left input inference (FINDINGS §57-58).

For each moved line and each per-line comparison (vs all / lineage / disease / subtype): calls in v3.1f and v3.2f,
matched by locus overlap (the two catalogues number their loci differently), and the CN at the calls: how many
v3.1f calls sit where the MEASURED CN is >= 2 (amplicon calls the inferred track missed) and whether they survive.

    python3 phase2/analysis/v32_moved_lines.py
"""
import os
import sys

import numpy as np
import pandas as pd

S = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, S)
from secacts_env import DATAROOT                                   # noqa: E402
ARMS = {"all": "atlas.s3.lines.all", "lin": "atlas.s3.lines.lin", "dis": "atlas.s3.lines.dis", "sub": "atlas.s3.lines.sub"}


def cat(res):
    c = pd.read_csv(os.path.join(S, res, "atlas.s3.union_catalog.bed.gz"), sep="\t", header=None, usecols=[0, 1, 2, 3],
                    names=["chrom", "start", "end", "se"])
    return c.set_index("se")


def calls(sc, arm, keys):
    d = pd.read_csv(os.path.join(S, sc, f"{arm}.line.specific.tsv.gz"), sep="\t")
    return d[(d.fdr <= 0.10) & d.group.isin(keys)]


def overlap_kept(a, b):
    """a, b: call tables with chrom/start/end; share of a's loci overlapping any b locus (same group)."""
    hit = np.zeros(len(a), dtype=bool)
    for (g, ch), x in a.groupby(["group", "chrom"]):
        y = b[(b.group == g) & (b.chrom == ch)]
        if y.empty:
            continue
        ys, ye = y.start.values, y.end.values
        for i, (s, e) in zip(x.index, zip(x.start, x.end)):
            hit[a.index.get_loc(i)] = bool(((ys < e) & (ye > s)).any())
    return hit


def main():
    p32 = pd.read_csv(os.path.join(S, "phase2/data/pull_set.v32.tsv"), sep="\t")
    p31 = pd.read_csv(os.path.join(S, "phase2/data/pull_set.v31.tsv"), sep="\t")
    src31 = p31.drop_duplicates("key").set_index("key").cn_provider
    moved = p32.drop_duplicates("key").set_index("key")
    moved = moved[(moved.cn_provider != "input_inferred") & (moved.index.map(src31) == "input_inferred")]
    # per-line tables name a line as score_pilot does: DepMap's StrippedCellLineName for a DepMap key (MUTZ-3 ->
    # MUTZ3), the cell name for a line outside DepMap (keyed by CVCL)
    stripped = pd.read_csv(os.path.join(DATAROOT, "DepMap/2026q1/Model.csv"), index_col="ModelID").StrippedCellLineName
    moved["cell"] = [stripped.get(k, c) for k, c in zip(moved.index, moved.cell)]
    moved = moved.reset_index().set_index("cell")
    moved.index.name = "cell"
    keys = list(moved.index)
    c31, c32 = cat("phase2/results_v31f"), cat("phase2/results_v32f")
    rows = []
    for k, arm in ARMS.items():
        a = calls("phase2/scores_v31f", arm, keys).join(c31, on="se").reset_index(drop=True)
        b = calls("phase2/scores_v32f", arm, keys).join(c32, on="se").reset_index(drop=True)
        ka = overlap_kept(a, b)
        kb = overlap_kept(b, a)
        # measured CN at the v3.1f calls: cn_mean in v3.2f's own table is the measured CN of the line at ITS loci,
        # so read the v3.1f loci's measured CN off the matching v3.2f call when there is one; else from v3.2f cn_mean
        for g in keys:
            ag, bg = a.group == g, b.group == g
            rows.append({"line": g, "key": moved.at[g, "key"], "source": moved.at[g, "cn_provider"], "comparison": k,
                         "v31f": int(ag.sum()), "v32f": int(bg.sum()), "kept": int(ka[ag].sum()),
                         "new": int((~kb[bg]).sum()),
                         "v31f_cn_mean_ge2": int((a.cn_mean[ag] >= 2).sum()),
                         "v32f_cn_mean_ge2": int((b.cn_mean[bg] >= 2).sum())})
    T = pd.DataFrame(rows)
    out = os.path.join(S, "phase2/analysis/out/v32_moved_lines.tsv")
    T.to_csv(out, sep="\t", index=False)
    with pd.option_context("display.width", 200, "display.max_rows", 100):
        print(T[T.comparison == "all"].drop(columns=["key", "comparison"]).to_string(index=False))
        s = T.groupby("comparison")[["v31f", "v32f", "kept", "new", "v31f_cn_mean_ge2", "v32f_cn_mean_ge2"]].sum()
        print(s.to_string())


if __name__ == "__main__":
    main()
