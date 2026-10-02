#!/usr/bin/env python3
"""v3.1 reduce: union catalogue from agnostic + calling-time-CN calls, a chosen merge rule, presence by overlap.

What changes against aggregate.py (which v1-v3 used and which stays as it was, so those releases rebuild exactly):
  * the catalogue can be built from several per-sample call sets (`--calls agnostic,cn`: `<SRX>.se.bed` and the
    calling-time-corrected `<SRX>.cn.se.bed` from recall_cn.py; GOTCHAS 59, v3.1 item 13);
  * the merge rule is selectable (`--merge`): recip25 (v3: >= 25% reciprocal overlap), contain50 (also merge
    when the overlap covers >= 50% of the smaller locus) or any (single linkage on any overlap) -- item 10;
  * presence is written per call set and by OVERLAP (a sample "calls" a locus when one of its SEs overlaps it),
    which is the rule v3.0.1 applied (74_called_filter.py) and the one scoring will require (item 9);
  * `--catalog-in` quantifies samples on a FIXED catalogue (variants 4-5 reuse the v3.1 baseline loci, so the
    dashboard can switch analyses per locus);
  * samples come from a list and are read from the sharded pull output (`out/<last 2 chars>/`) directly.
Signal, the grid sum and S3norm (pinned reference) are unchanged from aggregate.py.

    python3 phase2/aggregate_v31.py --samples pull_srx.v31.txt --out-dir $W/out --grid $W/data/grid.20.bed \\
        --grid-label grid.20 --peak-dir $W/data/bed20 --min-peaks 2000 --norm s3norm --s3-ref SRX16495452 \\
        --calls agnostic,cn --merge any --out $W/results_v31/atlas.s3
"""
from __future__ import annotations

import argparse, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from aggregate import read_grid, grid_index, rows_for_region, read_se_bed, read_signal_col   # noqa: E402

SUFFIX = {"agnostic": ".se.bed", "cn": ".cn.se.bed"}


def merge(ivs, rule):
    """ivs: [(chrom, start, end)]. Single-linkage merge under `rule`. Returns sorted [(chrom, start, end, n)]."""
    by = {}
    for c, s, e in ivs:
        by.setdefault(c, []).append((s, e))
    out = []
    for c, v in by.items():
        v.sort()
        n = len(v)
        parent = list(range(n))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        if rule == "any":                                   # sweep: a run of overlapping intervals is one locus
            cur_s, cur_e, k = v[0][0], v[0][1], 1
            for s, e in v[1:]:
                if s < cur_e:
                    cur_e, k = max(cur_e, e), k + 1
                else:
                    out.append((c, cur_s, cur_e, k)); cur_s, cur_e, k = s, e, 1
            out.append((c, cur_s, cur_e, k))
            continue
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
                if (ov / li >= 0.25 and ov / lj >= 0.25) or (rule == "contain50" and ov / min(li, lj) >= 0.5):
                    a, b = find(i), find(j)
                    if a != b:
                        parent[b] = a
        comp = {}
        for i in range(n):
            comp.setdefault(find(i), []).append(i)
        out += [(c, min(v[k][0] for k in m), max(v[k][1] for k in m), len(m)) for m in comp.values()]
    out.sort(key=lambda u: (u[0], u[1]))
    return out


class LocusIndex:
    """Overlap queries against a catalogue whose loci may overlap one another."""

    def __init__(self, loci):
        self.by = {}
        for i, (c, s, e) in enumerate(loci):
            self.by.setdefault(c, []).append((s, e, i))
        for c, v in self.by.items():
            v.sort()
            st = np.array([x[0] for x in v]); en = np.array([x[1] for x in v]); ii = np.array([x[2] for x in v])
            self.by[c] = (st, en, ii, int((en - st).max()))

    def hits(self, c, s, e):
        if c not in self.by:
            return ()
        st, en, ii, ml = self.by[c]
        lo, hi = np.searchsorted(st, s - ml, side="left"), np.searchsorted(st, e, side="left")
        m = en[lo:hi] > s
        return ii[lo:hi][m]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", required=True, help="SRX list")
    ap.add_argument("--out-dir", required=True, help="the pull's sharded out/")
    ap.add_argument("--grid", required=True)
    ap.add_argument("--grid-label", default="grid.20")
    ap.add_argument("--peak-dir", required=True)
    ap.add_argument("--min-peaks", type=int, default=2000)
    ap.add_argument("--norm", default="s3norm", choices=["none", "s3norm"])
    ap.add_argument("--s3-ref", default="SRX16495452")
    ap.add_argument("--calls", default="agnostic,cn", help="call sets that build the union catalogue")
    ap.add_argument("--merge", default="any", choices=["recip25", "contain50", "any"])
    ap.add_argument("--catalog-in", help="fixed catalogue BED (chrom start end id ...): no union is built")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    chroms, starts, ends = read_grid(a.grid)
    gidx = grid_index(chroms, starts, ends)
    want = [s.strip() for s in open(a.samples) if s.strip()]
    path = lambda s, tail: os.path.join(a.out_dir, s[-2:], s + tail)          # noqa: E731
    samples, dropped = [], []
    for s in want:
        if not os.path.exists(path(s, f".{a.grid_label}.f32")) or not os.path.exists(path(s, ".se.bed")):
            dropped.append((s, "missing")); continue
        p = os.path.join(a.peak_dir, f"{s}.20.bed")
        n = sum(1 for _ in open(p)) if os.path.exists(p) else 0
        if n < a.min_peaks:
            dropped.append((s, f"{n} peaks")); continue
        samples.append(s)
    print(f"[agg31] {len(samples)} of {len(want)} samples ({len(dropped)} dropped: "
          f"{sum(d[1] == 'missing' for d in dropped)} missing, the rest < {a.min_peaks} peaks)", file=sys.stderr,
          flush=True)
    calls = {k: [read_se_bed(path(s, SUFFIX[k])) if os.path.exists(path(s, SUFFIX[k])) else None for s in samples]
             for k in SUFFIX}
    no_cn = [s for s, x in zip(samples, calls["cn"]) if x is None]
    if no_cn:
        print(f"[agg31] {len(no_cn)} samples without .cn.se.bed (their corrected presence is their agnostic)",
              file=sys.stderr, flush=True)
    calls["cn"] = [x if x is not None else y for x, y in zip(calls["cn"], calls["agnostic"])]

    if a.catalog_in:
        loci = []
        with open(a.catalog_in) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                loci.append((f[0], int(f[1]), int(f[2]), int(f[5]) if len(f) > 5 else 0))
        print(f"[agg31] fixed catalogue {a.catalog_in}: {len(loci)} loci", file=sys.stderr, flush=True)
    else:
        ivs = [iv for k in a.calls.split(",") for per in calls[k] for iv in per]
        loci = merge(ivs, a.merge)
        w = np.array([e - s for _, s, e, _ in loci])
        print(f"[agg31] union ({a.calls}, merge={a.merge}): {len(loci)} loci from {len(ivs)} calls; "
              f"median {np.median(w) / 1e3:.1f} kb, p99 {np.quantile(w, .99) / 1e3:.0f} kb, max {w.max() / 1e3:.0f} kb, "
              f"< 2 kb {int((w < 2000).sum())}", file=sys.stderr, flush=True)

    # presence by overlap, per call set
    L = LocusIndex([(c, s, e) for c, s, e, _ in loci])
    pres = {}
    for k in SUFFIX:
        P = np.zeros((len(loci), len(samples)), dtype=np.uint8)
        for j, per in enumerate(calls[k]):
            for c, s, e in per:
                P[L.hits(c, s, e), j] = 1
        pres[k] = P

    M = np.empty((len(chroms), len(samples)), dtype=np.float32)
    for j, s in enumerate(samples):
        M[:, j] = read_signal_col(path(s, f".{a.grid_label}.f32"), len(chroms))
    if a.norm == "s3norm":
        from s3norm import s3norm_matrix
        M, s3params = s3norm_matrix(M, ref=a.s3_ref, srx=samples)
        s3params.to_csv(a.out + ".s3norm_params.tsv", sep="\t", index=False)
    sig = np.zeros((len(loci), len(samples)), dtype=np.float64)
    nrows = np.zeros(len(loci), dtype=int)
    for i, (c, s, e, _) in enumerate(loci):
        r = rows_for_region(gidx, c, s, e)
        nrows[i] = r.size
        if r.size:
            sig[i] = M[r].sum(axis=0)

    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    with open(a.out + ".union_catalog.bed", "w") as fh:
        for i, (c, s, e, n) in enumerate(loci):
            fh.write(f"{c}\t{s}\t{e}\tUSE_{i}\t{int(pres['agnostic'][i].sum())}\t{n}\t{nrows[i]}\t"
                     f"{int(pres['cn'][i].sum())}\n")
    hdr = "se_id\t" + "\t".join(samples) + "\n"
    with open(a.out + ".se_signal.tsv", "w") as fh:
        fh.write(hdr)
        for i in range(len(loci)):
            fh.write(f"USE_{i}\t" + "\t".join(f"{v:.6g}" for v in sig[i]) + "\n")
    for k, tail in (("agnostic", ".se_presence.tsv"), ("cn", ".se_presence.cn.tsv")):
        with open(a.out + tail, "w") as fh:
            fh.write(hdr)
            for i in range(len(loci)):
                fh.write(f"USE_{i}\t" + "\t".join(map(str, pres[k][i])) + "\n")
    with open(a.out + ".samples_dropped.tsv", "w") as fh:
        fh.write("srx\treason\n" + "".join(f"{s}\t{r}\n" for s, r in dropped))
    print(f"[agg31] wrote {a.out}.{{union_catalog.bed,se_signal.tsv,se_presence.tsv,se_presence.cn.tsv}}; "
          f"{int((nrows == 0).sum())} loci map to 0 grid rows", file=sys.stderr)


if __name__ == "__main__":
    main()
