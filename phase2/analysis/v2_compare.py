#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compare two permutation-scored atlases — the 282-line panel (v1) against the CN-expanded panel (v2).

The two runs have DIFFERENT union catalogs (adding samples re-clusters the SE loci), so SE ids are not
comparable; calls are matched by genomic overlap instead. Reports, per hierarchy level:

  PANEL      groups, group-size profile, calls, groups with >=1 call, median calls/group
  SUBTYPE    calls by group size — the question the expansion was run to answer
  STABILITY  of A's calls in groups present in both runs, the share B still calls (any overlapping SE,
             same group), and the share of B's calls in those groups that are new
  ANCHORS    best rank of each curated identity gene per group (hierarchy_summary `genes`)

  ~/miniconda3/envs/atac_hdac/bin/python phase2/analysis/v2_compare.py \\
      --a phase2/scores/atlas.s3.perm    --a-catalog phase2/results/atlas.s3.union_catalog.bed.gz \\
      --b phase2/scores_v2/atlas.s3.perm --b-catalog phase2/results_v2/atlas.s3.union_catalog.bed.gz \\
      --out phase2/scores_v2/compare_v1_v2
"""
from __future__ import annotations

import argparse
import gzip
import os

import numpy as np
import pandas as pd

LEVELS = ["OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype"]
SIZE_BINS = [0, 1, 4, 6, 9, 14, 10 ** 6]
SIZE_LABELS = ["1", "2-4", "5-6", "7-9", "10-14", "15+"]


def load_catalog(path):
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            rows.append((f[3], f[0], int(f[1]), int(f[2])))
    return pd.DataFrame(rows, columns=["se", "chrom", "start", "end"]).set_index("se")


def overlap_map(cat_a, cat_b):
    """{a_se: [b_se, ...]} for every B locus overlapping each A locus (any overlap)."""
    out = {}
    for chrom, A in cat_a.groupby("chrom"):
        B = cat_b[cat_b["chrom"] == chrom].sort_values("start")
        if B.empty:
            continue
        bs, be, bid = B["start"].values, B["end"].values, B.index.values
        emax = np.maximum.accumulate(be)                     # running max end, for the left bound
        for se, s, e in zip(A.index.values, A["start"].values, A["end"].values):
            hi = np.searchsorted(bs, e, "left")              # B starts before A ends
            lo = np.searchsorted(emax, s, "right")           # first B whose running max end passes A start
            hits = [bid[i] for i in range(lo, hi) if be[i] > s]
            if hits:
                out[se] = hits
    return out


def parse_genes(txt):
    if not isinstance(txt, str) or not txt:
        return {}
    return {g: int(r) for g, r in (x.split(":") for x in txt.split(";"))}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--a", required=True, help="scores prefix of run A (e.g. .../atlas.s3.perm)")
    ap.add_argument("--a-catalog", required=True)
    ap.add_argument("--b", required=True, help="scores prefix of run B")
    ap.add_argument("--b-catalog", required=True)
    ap.add_argument("--fdr", type=float, default=0.10)
    ap.add_argument("--label-a", default="v1")
    ap.add_argument("--label-b", default="v2")
    ap.add_argument("--out", required=True, help="output prefix")
    a = ap.parse_args()
    la, lb = a.label_a, a.label_b

    HA = pd.read_csv(f"{a.a}.hierarchy_summary.tsv", sep="\t")
    HB = pd.read_csv(f"{a.b}.hierarchy_summary.tsv", sep="\t")
    cat_a, cat_b = load_catalog(a.a_catalog), load_catalog(a.b_catalog)
    omap = overlap_map(cat_a, cat_b)
    print(f"[compare] catalogs: {la} {len(cat_a):,} SEs, {lb} {len(cat_b):,} SEs; "
          f"{len(omap):,} {la} SEs ({100 * len(omap) / len(cat_a):.1f}%) overlap a {lb} SE")

    panel, bysize, stab, anchors = [], [], [], []
    for lev in LEVELS:
        for lab, H in ((la, HA), (lb, HB)):
            h = H[H.level == lev]
            n = h["n_lines"]
            panel.append(dict(level=lev, run=lab, groups=len(h), groups_n1=int((n == 1).sum()),
                              groups_le4=int((n <= 4).sum()), groups_ge7=int((n >= 7).sum()),
                              calls=int(h["n_spec_fdr10"].sum()),
                              groups_with_call=int((h["n_spec_fdr10"] > 0).sum()),
                              median_calls_per_group=float(h["n_spec_fdr10"].median()),
                              lines_total=int(n.sum())))
            cut = pd.cut(n, SIZE_BINS, labels=SIZE_LABELS)
            g = h.assign(size=cut).groupby("size", observed=False).agg(
                groups=("group", "size"), calls=("n_spec_fdr10", "sum"),
                groups_with_call=("n_spec_fdr10", lambda x: int((x > 0).sum())))
            for sz, r in g.iterrows():
                bysize.append(dict(level=lev, run=lab, group_size=sz, groups=int(r.groups),
                                   calls=int(r.calls), groups_with_call=int(r.groups_with_call)))

        # stability — only where both runs wrote a call dump
        fa, fb = f"{a.a}.{lev}.specific.tsv.gz", f"{a.b}.{lev}.specific.tsv.gz"
        if os.path.exists(fa) and os.path.exists(fb):
            SA = pd.read_csv(fa, sep="\t"); SA = SA[SA.fdr <= a.fdr]
            SB = pd.read_csv(fb, sep="\t"); SB = SB[SB.fdr <= a.fdr]
            shared = sorted(set(HA[HA.level == lev].group) & set(HB[HB.level == lev].group))
            b_calls = set(zip(SB.group, SB.se))
            b_matched = set()
            for grp in shared:
                sa = SA[SA.group == grp]
                kept = 0
                for se in sa.se:
                    hit = [b for b in omap.get(se, []) if (grp, b) in b_calls]
                    if hit:
                        kept += 1
                        b_matched.update((grp, b) for b in hit)
                nb = int((SB.group == grp).sum())
                nb_new = nb - sum(1 for (g2, _) in b_matched if g2 == grp)
                na_lines = int(HA[(HA.level == lev) & (HA.group == grp)].n_lines.iloc[0])
                nb_lines = int(HB[(HB.level == lev) & (HB.group == grp)].n_lines.iloc[0])
                stab.append(dict(level=lev, group=grp, lines_a=na_lines, lines_b=nb_lines,
                                 calls_a=len(sa), calls_b=nb, a_retained=kept,
                                 retained_frac=round(kept / len(sa), 4) if len(sa) else np.nan,
                                 b_new=nb_new))

        for grp in sorted(set(HA[HA.level == lev].group) | set(HB[HB.level == lev].group)):
            ra = HA[(HA.level == lev) & (HA.group == grp)]
            rb = HB[(HB.level == lev) & (HB.group == grp)]
            ga = parse_genes(ra.genes.iloc[0]) if len(ra) else {}
            gb = parse_genes(rb.genes.iloc[0]) if len(rb) else {}
            for gene in sorted(set(ga) | set(gb)):
                anchors.append(dict(level=lev, group=grp, gene=gene, rank_a=ga.get(gene), rank_b=gb.get(gene)))

    P, Z, S, N = (pd.DataFrame(x) for x in (panel, bysize, stab, anchors))
    P.to_csv(f"{a.out}.panel.tsv", sep="\t", index=False)
    Z.to_csv(f"{a.out}.by_group_size.tsv", sep="\t", index=False)
    S.to_csv(f"{a.out}.stability.tsv", sep="\t", index=False)
    N.to_csv(f"{a.out}.anchors.tsv", sep="\t", index=False)

    pd.set_option("display.width", 200)
    print("\n== PANEL ==\n" + P.to_string(index=False))
    print("\n== CALLS BY GROUP SIZE ==\n" + Z.pivot_table(index=["level", "group_size"], columns="run",
                                                            values=["groups", "calls", "groups_with_call"],
                                                            observed=False).to_string())
    if len(S):
        for lev, s in S.groupby("level"):
            tot_a, kept = s.calls_a.sum(), s.a_retained.sum()
            print(f"\n== STABILITY {lev}: {len(s)} shared groups; {la} calls retained in {lb}: "
                  f"{kept:,}/{tot_a:,} = {100 * kept / max(tot_a, 1):.1f}%; {lb} calls there: "
                  f"{s.calls_b.sum():,} ({s.b_new.sum():,} with no {la} counterpart) ==")
            print(s.sort_values("lines_b", ascending=False).head(30).to_string(index=False))
    print(f"\n[compare] wrote {a.out}.{{panel,by_group_size,stability,anchors}}.tsv")


if __name__ == "__main__":
    main()
