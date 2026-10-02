#!/usr/bin/env python3
"""Compare union-catalogue merge rules before v3.1 (ROADMAP v3.1 item 10, FINDINGS §28/§39).

v3 merges per-sample SEs by >= 25% RECIPROCAL overlap (cSEAdb), which keeps a small SE apart from the large SE
around it: 37% of loci overlap another and 20% are < 2 kb. Candidates, on agnostic and on agnostic + corrected
(calling-time CN) calls:
  recip25   v3's rule
  contain50 recip25, plus merge when the overlap covers >= 50% of the SMALLER interval (nest-aware)
  any       single-linkage on any overlap (no overlapping loci left; may chain neighbours into long loci)
For each: loci, size quantiles, loci < 2 kb, loci > 500 kb, members of the largest locus, overlapping pairs left.

    python3 phase2/catalog_eval.py --srx-list pull_srx.txt --out-dir $W/out --out catalog_eval.tsv
"""
import argparse, os, sys
import numpy as np
import pandas as pd


def read_bed(path):
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return [(f[0], int(f[1]), int(f[2])) for f in (l.split("\t") for l in fh) if len(f) >= 3]


def merge(ivs, rule):
    """ivs: list of (chrom, start, end). Returns merged loci [(chrom, start, end, n_members)]."""
    out = []
    by = {}
    for c, s, e in ivs:
        by.setdefault(c, []).append((s, e))
    for c, v in by.items():
        v.sort()
        n = len(v)
        parent = list(range(n))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        for i in range(n):
            si, ei = v[i]
            for j in range(i + 1, n):
                sj, ej = v[j]
                if sj >= ei:
                    break
                ov = min(ei, ej) - sj
                if ov <= 0:
                    continue
                li, lj = ei - si, ej - sj
                ok = (rule == "any" or (ov / li >= 0.25 and ov / lj >= 0.25)
                      or (rule == "contain50" and ov / min(li, lj) >= 0.5))
                if ok:
                    a, b = find(i), find(j)
                    if a != b:
                        parent[b] = a
        comp = {}
        for i in range(n):
            comp.setdefault(find(i), []).append(i)
        for m in comp.values():
            out.append((c, min(v[k][0] for k in m), max(v[k][1] for k in m), len(m)))
    return out


def overlapping_pairs(loci):
    n = 0
    by = {}
    for c, s, e, _ in loci:
        by.setdefault(c, []).append((s, e))
    for v in by.values():
        v.sort()
        for i in range(len(v)):
            j = i + 1
            while j < len(v) and v[j][0] < v[i][1]:
                n += 1
                j += 1
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--srx-list", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    srx = [s.strip() for s in open(a.srx_list) if s.strip()]
    ag, cn, miss = [], [], 0
    for s in srx:
        base = os.path.join(a.out_dir, s[-2:], s)
        x = read_bed(base + ".se.bed")
        miss += not x
        ag += x
        cn += read_bed(base + ".cn.se.bed")
    print(f"[cat] {len(srx)} samples ({miss} without SEs): {len(ag)} agnostic, {len(cn)} corrected calls",
          file=sys.stderr, flush=True)
    rows = []
    for src, ivs in (("agnostic", ag), ("agnostic+cn", ag + cn)):
        for rule in ("recip25", "contain50", "any"):
            L = merge(ivs, rule)
            w = np.array([e - s for _, s, e, _ in L])
            rows.append(dict(calls=src, rule=rule, loci=len(L), median_kb=np.median(w) / 1e3,
                             p90_kb=np.quantile(w, 0.9) / 1e3, p99_kb=np.quantile(w, 0.99) / 1e3, max_kb=w.max() / 1e3,
                             lt2kb=int((w < 2000).sum()), gt500kb=int((w > 500_000).sum()),
                             max_members=max(m for *_, m in L), overlapping_pairs=overlapping_pairs(L)))
            print(rows[-1], file=sys.stderr, flush=True)
    pd.DataFrame(rows).round(2).to_csv(a.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
